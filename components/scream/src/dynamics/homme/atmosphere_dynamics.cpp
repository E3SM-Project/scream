#include "atmosphere_dynamics.hpp"

// HOMMEXX Includes
#include "Context.hpp"
#include "Elements.hpp"
#include "ElementsState.hpp"
#include "HommexxEnums.hpp"
#include "SimulationParams.hpp"
#include "ElementsForcing.hpp"
#include "EulerStepFunctor.hpp"
#include "Diagnostics.hpp"
#include "DirkFunctor.hpp"
#include "ForcingFunctor.hpp"
#include "CaarFunctor.hpp"
#include "VerticalRemapManager.hpp"
#include "HyperviscosityFunctor.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"
#include "mpi/ConnectivityHelpers.hpp"
#include "utilities/MathUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "ExecSpaceDefs.hpp"
#include "Types.hpp"

// Scream includes
#include "dynamics/homme/physics_dynamics_remapper.hpp"
#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/util//scream_column_ops.hpp"
#include "share/field/field_property_checks/field_lower_bound_check.hpp"
#include "share/field/field_property_checks/field_nan_check.hpp"
#include "share/field/field_utils.hpp"

// Ekat includes
#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos//ekat_subview_utils.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream
{

HommeDynamics::HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // This class needs Homme's context, so register as a user
  HommeContextUser::singleton().add_user();
}

void HommeDynamics::set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  // Init prim structures
  // TODO: they should not be inited yet; should we error out if they are?
  //       I'm gonna say 'no', for now, cause it might be a pb with unit tests.
  if (!is_data_structures_inited_f90()) {
    prim_init_data_structures_f90 ();
  }

  // Note: time levels are just an expedient used by Homme to
  //  store temporaries in the RK timestepping schemes.
  //  It is best to have this extra array dimension (rather than,
  //  say, having NTL separate arrays) because of memory locality.
  //  At the end of Homme's timestep, only one of those slices
  //  will be meaningful. The phys-dyn remapper will use Homme's
  //  TimeLevel structure to know exactly where to copy data from/to
  //  during the remap.

  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using FL  = FieldLayout;

  constexpr int NGP = HOMMEXX_NP;
  constexpr int NTL = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int QTL = HOMMEXX_Q_NUM_TIME_LEVELS;
  constexpr int N   = HOMMEXX_PACK_SIZE;

  // Some units
  const auto nondim = Units::nondimensional();
  const auto rho = kg/(m*m*m);
  auto Q = kg/kg;
  Q.set_string("kg/kg");

  const auto dgn = "Dynamics";
  m_dyn_grid = grids_manager->get_grid(dgn);
  m_ref_grid = grids_manager->get_reference_grid();
  const auto& rgn = m_ref_grid->name();

  const int ncols = m_ref_grid->get_num_local_dofs();
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int nlev_mid = m_dyn_grid->get_num_vertical_levels();
  const int nlev_int = nlev_mid+1;

  // Sanity check for the grid. This should *always* pass, since Homme builds the grids
  EKAT_REQUIRE_MSG(get_num_local_elems_f90()==nelem,
      "Error! The number of elements computed from the Dynamics grid num_dof()\n"
      "       does not match the number of elements internal in Homme.\n");

  /*
     Explanation of what Homme needs, including what it needs for BFB restarts

     From the pt of view of the atmosphere, Homme has v,w,dp,T,phi,Q on physics grid
     as both inputs and outputs.
     Homme stores *prognostic states* on the dyn grid (v, w, dp, vtheta_dp, Qdp, phi),
     which have to be saved to file for BFB restart.
     Additionally, Homme computes forcing as:
       - FT (temperature forcing): FT_dyn=PD( (T_phys_new - T_phys_old)/dt )
       - FM (momentum forcing): FM_dyn=PD( (v_phys_new-v_phys_old)/dt) )
       - FQ (tracers forcing):
          - ftype==FORCING_DEBUG(0): FQ_dyn =  PD(Q_phys_new) (used as a hard adjustment
          - ftype==FORCING_2(2): FQ_dyn = PD(Q_phys_new - DP(Q_dyn_old)) (an increment for Q)
     Here, PD(x) is quantity x remapped from physics to dynamics grid (viceversa for DP(x)).
     So for BFB restart, we need so store T_phys_old, [u,v,w]_phys_old, and, if
     forcing type is FORCING_2, also Q_dyn_old.
     However, Homme computes Q=Qdp/dp at the end of the time step, so for BFB restarts we need to
     save Qdp only, and recompute Q_dyn_old=Qdp/dp.
   */
  // Notes:
  //  - physics will update T_mid, Q, horiz_winds, but not w_i, phi_i, and pseudo_density.
  //    Dyn will have to back out tendencies, as well as convert T->VTheta_dp.
  //    The simplest way is to store T, uv, and Q at the end of the previous dyn step.
  //    To make SCREAM more similar to what EAM does, we store Q_prev on dyn grid,
  //    while we keep uv_prev and T_prev on the phys grid (better for storage as well).
  //  - physics *does* update phi_i after each parametrization, but dyn discards it,
  //    in favor of reconstructing it from the FT, FQ, and the other thermodyn vars.

  // Note: the field T_mid will contain Temperature (at midpoints) at entry/exit.
  //       However, Homme uses VTheta*dp as state variable. Instead of keeping both
  //       T and VTheta*dp on dyn grid, we remap T_mid from phys directly
  //       into the n0 time-slice of Homme's vtheta_dp, and then do the conversion
  //       T_mid->VTheta_dp in place.

  // Note: qv is needed to transform T<->Theta
  add_field<Updated> ("horiz_winds",   FL({COL,CMP, LEV},{ncols,2,nlev_mid}),m/s,   rgn,N);
  add_field<Updated> ("T_mid",         FL({COL,     LEV},{ncols,  nlev_mid}),K,     rgn,N);
  add_field<Updated> ("w_int",         FL({COL,    ILEV},{ncols,  nlev_int}),m/s,   rgn,N);
  add_field<Updated> ("phi_int",       FL({COL,    ILEV},{ncols,  nlev_int}),Pa/rho,rgn,N);
  add_field<Updated> ("pseudo_density",FL({COL,     LEV},{ncols,  nlev_mid}),Pa,    rgn,N);
  add_field<Updated> ("ps",            FL({COL         },{ncols           }),Pa,    rgn);
  add_field<Required>("qv",            FL({COL,     LEV},{ncols,  nlev_mid}),Q,     rgn,"tracers",N);
  add_field<Computed>("p_int",         FL({COL,    ILEV},{ncols,  nlev_int}),Pa,    rgn,N);
  add_field<Computed>("p_mid",         FL({COL,     LEV},{ncols,  nlev_mid}),Pa,    rgn,N);
  add_group<Updated>("tracers",rgn,N, Bundling::Required);

  // Dynamics grid states
  create_helper_field("v_dyn",        {EL,TL,CMP,GP,GP,LEV}, {nelem,NTL,2,NP,NP,nlev_mid}, dgn);
  create_helper_field("vtheta_dp_dyn",{EL,TL,    GP,GP,LEV}, {nelem,NTL,  NP,NP,nlev_mid}, dgn);
  create_helper_field("dp3d_dyn",     {EL,TL,    GP,GP,LEV}, {nelem,NTL,  NP,NP,nlev_mid}, dgn);
  create_helper_field("w_int_dyn",    {EL,TL,    GP,GP,ILEV},{nelem,NTL,  NP,NP,nlev_int}, dgn);
  create_helper_field("phi_int_dyn",  {EL,TL,    GP,GP,ILEV},{nelem,NTL,  NP,NP,nlev_int}, dgn);
  create_helper_field("ps_dyn",       {EL,TL,    GP,GP},     {nelem,NTL,  NP,NP         }, dgn);
  create_helper_field("Qdp_dyn",      {EL,TL,CMP,GP,GP,LEV}, {nelem,QTL,HOMMEXX_QSIZE_D,NP,NP,nlev_mid},dgn);

  // For BFB restart, we need to read in the state on the dyn grid. The state above has NTL time slices,
  // but only one is really needed for restart. Therefore, we create "dynamic" subfields for
  // the state fields. This allows to save/read only the single slice needed
  // NOTE: the fcn init_time_level_c in Homme should really init also the qdp time levels,
  //       but it doesn't. So let's update them here. Notice that the qdp time levels update
  //       is "idempotent", meaning that it's safe to call it many times. In fact, all it
  //       does it to compute np1 and n0 from nsteps, so if nsteps is unchanged, calling
  //       the fcn twice won't "update" n0 and np1.
  auto& tl = Homme::Context::singleton().get<Homme::TimeLevel>();
  const auto& params = Homme::Context::singleton().get<Homme::SimulationParams>();
  tl.update_tracers_levels(params.qsplit);

  add_internal_field (m_helper_fields.at("v_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("vtheta_dp_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("dp3d_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("w_int_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("phi_int_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("ps_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("Qdp_dyn").subfield(1,tl.n0_qdp,true));

  // Dynamics backs out tendencies from the states, and passes those to Homme.
  // After Homme completes, we remap the updated state to the ref grid.
  // Thus, is more convenient to use two different remappers: the pd remapper
  // will remap into Homme's forcing views, while the dp remapper will remap
  // from Homme's states.
  m_p2d_remapper = grids_manager->create_remapper(m_ref_grid,m_dyn_grid);
  m_d2p_remapper = grids_manager->create_remapper(m_dyn_grid,m_ref_grid);

  // Create separate remapper for Initial Conditions
  m_ic_remapper = grids_manager->create_remapper(m_ref_grid,m_dyn_grid);
}

int HommeDynamics::requested_buffer_size_in_bytes() const
{
  using namespace Homme;

  auto& c = Context::singleton();
  const auto num_elems   = c.get<Elements>().num_elems();
  auto& params  = c.get<SimulationParams>();

  auto& caar = c.create_if_not_there<CaarFunctor>(num_elems,params);
  auto& esf  = c.create_if_not_there<EulerStepFunctor>(num_elems);
  auto& hvf  = c.create_if_not_there<HyperviscosityFunctor>(num_elems, params);
  auto& ff   = c.create_if_not_there<ForcingFunctor>(num_elems, num_elems, params.qsize);
  auto& diag = c.create_if_not_there<Diagnostics> (num_elems,params.theta_hydrostatic_mode);
  auto& vrm  = c.create_if_not_there<VerticalRemapManager>(num_elems);

  const bool need_dirk = (params.time_step_type==TimeStepType::IMEX_KG243 ||
                          params.time_step_type==TimeStepType::IMEX_KG254 ||
                          params.time_step_type==TimeStepType::IMEX_KG255 ||
                          params.time_step_type==TimeStepType::IMEX_KG355);
  if (need_dirk) {
    // Create dirk functor only if needed
    c.create_if_not_there<DirkFunctor>(num_elems);
  }

  // Request buffer sizes in FunctorsBuffersManager and then
  // return the total bytes using the calculated buffer size.
  auto& fbm  = c.create_if_not_there<FunctorsBuffersManager>();
  fbm.request_size(caar.requested_buffer_size());
  fbm.request_size(esf.requested_buffer_size());
  fbm.request_size(hvf.requested_buffer_size());
  fbm.request_size(diag.requested_buffer_size());
  fbm.request_size(ff.requested_buffer_size());
  fbm.request_size(vrm.requested_buffer_size());
  if (need_dirk) {
    const auto& dirk = c.get<DirkFunctor>();
    fbm.request_size(dirk.requested_buffer_size());
  }

  return fbm.allocated_size()*sizeof(Real);
}

void HommeDynamics::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Buffers size not sufficient.\n");

  using namespace Homme;
  auto& c = Context::singleton();
  auto& fbm  = c.get<FunctorsBuffersManager>();

  // Reset Homme buffer to use AD buffer memory.
  // Internally, homme will actually initialize its own buffers.
  EKAT_REQUIRE(buffer_manager.allocated_bytes()%sizeof(Real)==0); // Sanity check
  fbm.allocate(buffer_manager.get_memory(), buffer_manager.allocated_bytes()/sizeof(Real));
}

void HommeDynamics::initialize_impl (const RunType run_type)
{
  const auto& dgn = m_dyn_grid->name();
  const auto& rgn = m_ref_grid->name();

  // Use common/shorter names for tracers.
  alias_group_in  ("tracers",rgn,"Q");
  alias_group_out ("tracers",rgn,"Q");

  // ------ Sanity checks ------- //

  // Nobody should claim to be a provider for dp, w_i.
  // WARNING! If the assumption on 'pseudo_density' ceasaes to be true, you have to revisit
  //          how you restart homme. In particular, p_mid is restarted from pseudo_density,
  //          as it is read from restart file. If other procs update it, the restarted value
  //          might no longer match the end-of-homme-step value, which is what you need
  //          to compute p_mid. Hence, if this assumption goes away, you need to restart
  //          p_mid by first remapping the reestarted dp3d_dyn back to ref grid, and using
  //          that value to compute p_mid.
  const auto& rho_track = get_field_out("pseudo_density").get_header().get_tracking();
  const auto& w_i_track = get_field_out("w_int").get_header().get_tracking();
  EKAT_REQUIRE_MSG (
      rho_track.get_providers().size()==1,
      "Error! Someone other than Dynamics is trying to update the pseudo_density.\n");
  EKAT_REQUIRE_MSG (
      w_i_track.get_providers().size()==1,
      "Error! Someone other than Dynamics is trying to update the vertical velocity.\n");

  // The groups 'tracers' and 'tracers_mass_dyn' should contain the same fields
  auto Q = get_group_out("Q",rgn);
  EKAT_REQUIRE_MSG(not Q.m_info->empty(),
    "Error! There should be at least one tracer (qv) in the tracers group.\n");

  // Create remaining internal fields
  const auto& c = Homme::Context::singleton();
  const auto& params = c.get<Homme::SimulationParams>();
  constexpr int NGP  = HOMMEXX_NP;
  const int qsize = params.qsize;
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int ncols = m_ref_grid->get_num_local_dofs();
  const int nlevs = m_dyn_grid->get_num_vertical_levels();

  using namespace ShortFieldTagsNames;
  create_helper_field("FQ_dyn",{EL,CMP,GP,GP,LEV},{nelem,qsize,NGP,NGP,nlevs},dgn);
  create_helper_field("FT_dyn",{EL,    GP,GP,LEV},{nelem,      NGP,NGP,nlevs},dgn);
  create_helper_field("FM_dyn",{EL,CMP,GP,GP,LEV},{nelem,    3,NGP,NGP,nlevs},dgn);
  create_helper_field("Q_dyn" ,{EL,CMP,GP,GP,LEV},{nelem,qsize,NGP,NGP,nlevs},dgn);
  // Tendencies for temperature and momentum computed on ref grid
  create_helper_field("FT_ref",{COL,LEV},    {ncols,  nlevs},rgn);
  create_helper_field("FM_ref",{COL,CMP,LEV},{ncols,3,nlevs},rgn);

  // Setup the p2d and d2p remappers
  m_p2d_remapper->registration_begins();
  m_d2p_remapper->registration_begins();

  // ftype==FORCING_DEBUG:
  //  1) remap Q_rgn->Q_dyn
  //  2) compute FQ_dyn=(Q_dyn-Q_dyn_old)/dt
  // ftype!=FORCING_DEBUG:
  //  1) compute FQ_rgn=Q_rgn-Q_rgn_old
  //  2) remap FQ_rgn->FQ_dyn
  //  3) Q_dyn=Q_dyn_old+FQ_dyn
  if (params.ftype==Homme::ForcingAlg::FORCING_2) {
    // Need a tmp for dQ_rgn
    using namespace ShortFieldTagsNames;
    create_helper_field("FQ_ref",{COL,CMP,LEV},{ncols,qsize,nlevs},rgn);
    m_p2d_remapper->register_field(m_helper_fields.at("FQ_ref"),m_helper_fields.at("FQ_dyn"));
  } else {
    // Can remap Q directly into FQ, tendency computed in pre_process step
    m_p2d_remapper->register_field(*get_group_out("Q",rgn).m_bundle,m_helper_fields.at("FQ_dyn"));
  }
  m_p2d_remapper->register_field(m_helper_fields.at("FT_ref"),m_helper_fields.at("FT_dyn"));
  m_p2d_remapper->register_field(m_helper_fields.at("FM_ref"),m_helper_fields.at("FM_dyn"));

  // NOTE: for states, if/when we can remap subfields, we can remap the corresponding internal fields,
  //       which are subviews of the corresponding helper field at time slice np1
  m_d2p_remapper->register_field(m_helper_fields.at("vtheta_dp_dyn"),get_field_out("T_mid"));
  m_d2p_remapper->register_field(m_helper_fields.at("v_dyn"),get_field_out("horiz_winds"));
  m_d2p_remapper->register_field(m_helper_fields.at("dp3d_dyn"), get_field_out("pseudo_density"));
  m_d2p_remapper->register_field(m_helper_fields.at("phi_int_dyn"), get_field_out("phi_int"));
  m_d2p_remapper->register_field(m_helper_fields.at("w_int_dyn"), get_field_out("w_int"));
  m_d2p_remapper->register_field(m_helper_fields.at("ps_dyn"), get_field_out("ps"));
  m_d2p_remapper->register_field(m_helper_fields.at("Q_dyn"),*get_group_out("Q",rgn).m_bundle);

  m_p2d_remapper->registration_ends();
  m_d2p_remapper->registration_ends();

  // Sets the scream views into the hommexx internal data structures
  init_homme_views ();

  // Import I.C. from the ref grid to the dyn grid.
  if (run_type==RunType::Initial) {
    initialize_homme_state ();
  } else {
    restart_homme_state ();
  }

  // For BFB restarts, set nstep counter in Homme's TimeLevel to match
  // what's in the timestamp (which, for restarted runs, is
  set_homme_param("num_steps",timestamp().get_num_steps());
  Homme::Context::singleton().get<Homme::TimeLevel>().nstep = timestamp().get_num_steps();

  // Complete homme model initialization
  prim_init_model_f90 ();
}

void HommeDynamics::run_impl (const int dt)
{
  try {
    Kokkos::fence();
    homme_pre_process (dt);

    // Note: Homme's step lasts homme_dt*max(dt_remap_factor,dt_tracers_factor), and it must divide dt.
    // We neeed to compute dt/homme_dt, and subcycle homme that many times
    const int nsplit = get_homme_nsplit_f90(dt);
    for (int subiter=0; subiter<nsplit; ++subiter) {
      Kokkos::fence();
      prim_run_f90 ();
    }

    // Post process Homme's output, to produce what the rest of Atm expects
    Kokkos::fence();
    homme_post_process ();
  } catch (std::exception& e) {
    EKAT_ERROR_MSG(e.what());
  } catch (...) {
    EKAT_ERROR_MSG("Something went wrong, but we don't know what.\n");
  }
}

void HommeDynamics::finalize_impl (/* what inputs? */)
{
  prim_finalize_f90();

  // This class is done needing Homme's context, so remove myself as customer
  Homme::Context::singleton().finalize_singleton();
}

void HommeDynamics::set_computed_group_impl (const FieldGroup<Real>& group)
{
  if (group.m_info->m_group_name=="tracers") {
    // Set runtime number of tracers in Homme
    const auto& c = Homme::Context::singleton();
    auto& params = c.get<Homme::SimulationParams>();
    auto& tracers = c.get<Homme::Tracers>();
    const int qsize = group.m_info->size();
    params.qsize = qsize;           // Set in the CXX data structure
    set_homme_param("qsize",qsize); // Set in the F90 module
    tracers.init(tracers.num_elems(),qsize);
  }
}

void HommeDynamics::homme_pre_process (const int dt) {
  // T and uv tendencies are backed out on the ref grid.
  // Homme takes care of turning the FT tendency into a tendency for VTheta_dp.
  using KT = KokkosTypes<DefaultDevice>;

  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using Pack = ekat::Pack<Real,N>;

  const int ncols = m_ref_grid->get_num_local_dofs();
  const int nlevs = m_ref_grid->get_num_vertical_levels();
  const int npacks = ekat::PackInfo<N>::num_packs(nlevs);

  const auto& rgn = m_ref_grid->name();

  // At the beginning of the step, FT and FM store T_prev and V_prev,
  // the temperature and (3d) velocity at the end of the previous
  // Homme step
  auto T  = get_field_in("T_mid").get_view<const Pack**>();
  auto v  = get_field_in("horiz_winds").get_view<const Pack***>();
  auto w  = get_field_in("w_int").get_view<const Pack**>();
  auto FT = m_helper_fields.at("FT_ref").get_view<Pack**>();
  auto FM = m_helper_fields.at("FM_ref").get_view<Pack***>();

  // If there are other atm procs updating the vertical velocity,
  // then we need to compute forcing for w as well
  const bool has_w_forcing = get_field_out("w_int").get_header().get_tracking().get_providers().size()>1;
  Kokkos::parallel_for(KT::RangePolicy(0,ncols*npacks),
                       KOKKOS_LAMBDA(const int& idx) {
    const int icol = idx / npacks;
    const int ilev = idx % npacks;

    // Temperature forcing
    // Note: Homme takes care of converting ft into a forcing for vtheta
    const auto& t_new =  T(icol,ilev);
    const auto& t_old = FT(icol,ilev);
          auto& ft    = FT(icol,ilev);
    ft = (t_new - t_old) / dt;

    // Horizontal velocity forcing
    const auto& u_new =  v(icol,0,ilev);
    const auto& u_old = FM(icol,0,ilev);
          auto& fu    = FM(icol,0,ilev);
    fu = (u_new - u_old) / dt;

    const auto& v_new =  v(icol,1,ilev);
    const auto& v_old = FM(icol,1,ilev);
          auto& fv    = FM(icol,1,ilev);
    fv = (v_new - v_old) / dt;

    if (has_w_forcing) {
      // Vertical velocity forcing
      // Recall: fm(2) stores forcing for w_i at [0,num_int_levels-1],
      //         since w_i at surf is determined with no penetration bc.
      const auto& w_new =  w(icol,ilev);
      const auto& w_old = FM(icol,2,ilev);
            auto& fw    = FM(icol,2,ilev);
      fw = (w_new - w_old) / dt;
    }
  });

  using namespace Homme;
  const auto& c = Context::singleton();
  const auto& params = c.get<SimulationParams>();
  const auto ftype = params.ftype;
  if (ftype==Homme::ForcingAlg::FORCING_2) {
    const auto Q  = get_group_in("Q",rgn).m_bundle->get_view<const Pack***>();
    const auto FQ = m_helper_fields.at("FQ_ref").get_view<Pack***>();
    const int qsize = Q.extent(1);
    Kokkos::parallel_for(Kokkos::RangePolicy<>(0,ncols*qsize*npacks),
                         KOKKOS_LAMBDA (const int idx) {
      const int icol =  idx / (qsize*npacks);
      const int iq   = (idx / (npacks)) % qsize;
      const int ilev = idx % npacks;
      // Recall: at this point FQ stores Q_ref at the end of previous homme step
      FQ(icol,iq,ilev) = Q(icol,iq,ilev) - FQ(icol,iq,ilev);
    });
  }

  // Remap FT, FM, and Q (or FQ, depending on ftype)
  m_p2d_remapper->remap(true);

  auto& tl = c.get<TimeLevel>();
  auto& ff = c.get<ForcingFunctor>();

  const auto& tracers = c.get<Tracers>();
  const auto& state = c.get<ElementsState>();
  auto Q    = tracers.Q;
  auto FQ   = tracers.fq;
  auto Qdp  = tracers.qdp;
  auto dp3d = state.m_dp3d;

  // Note: np1_qdp and n0_qdp are 'deduced' from tl.nstep, so the
  //       following call may not even change them (i.e., they are
  //       not updated regardless). So if they were already up-to-date,
  //       the following call will do nothing.
  tl.update_tracers_levels(params.qsplit);

  const auto n0 = tl.n0;  // The time level where pd coupling remapped into
  constexpr int NVL = HOMMEXX_NUM_LEV;
  const int qsize = params.qsize;
  const auto n0_qdp = tl.n0_qdp;
  const auto Q_dyn = m_helper_fields.at("Q_dyn").get_view<Homme::Scalar*****>();

  // At this point, FQ contains Qnew (coming from physics).
  // Depending on ftype, we are going to do different things:
  //  ftype=0: FQ = dp*(Qnew-Qold) / dt
  //  ftype=2: qdp = dp*Qnew, with Qnew stored in FQ
  switch(ftype) {
    case ForcingAlg::FORCING_0:
      // Back out tracers tendency for Qdp
      Kokkos::parallel_for(Kokkos::RangePolicy<>(0,Q.size()),KOKKOS_LAMBDA(const int idx) {
        const int ie = idx / (qsize*NP*NP*NVL);
        const int iq = (idx / (NP*NP*NVL)) % qsize;
        const int ip = (idx / (NP*NVL)) % NP;
        const int jp = (idx / NVL) % NP;
        const int k  =  idx % NVL;

        // fq is currently storing q_new
        const auto& q_prev = Q(ie,iq,ip,jp,k);
        const auto& dp     = dp3d(ie,n0,ip,jp,k);
              auto& fq     = FQ(ie,iq,ip,jp,k);

        fq -= q_prev;
        fq /= dt;
        fq *= dp;
      });
      Kokkos::fence();
      break;
    case ForcingAlg::FORCING_2:
      Kokkos::parallel_for(Kokkos::RangePolicy<>(0,Q.size()),KOKKOS_LAMBDA(const int idx) {
        const int ie = idx / (qsize*NP*NP*NVL);
        const int iq = (idx / (NP*NP*NVL)) % qsize;
        const int ip = (idx / (NP*NVL)) % NP;
        const int jp = (idx / NVL) % NP;
        const int k  =  idx % NVL;

        // So far, FQ contains the remap of Q_ref_new - Q_ref_old. Add Q_dyn_old to get new state.
        FQ(ie,iq,ip,jp,k) += Q_dyn(ie,iq,ip,jp,k);
      });

      // Hard adjustment of qdp
      ff.tracers_forcing(dt,n0,n0_qdp,true,params.moisture);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unexpected/unsupported forcing algorithm.\n"
                      "  ftype: " + std::to_string(Homme::etoi(ftype)) + "\n");
  }
}

void HommeDynamics::homme_post_process () {
  const auto& rgn = m_ref_grid->name();
  const auto& c = Homme::Context::singleton();
  const auto& tl = c.get<Homme::TimeLevel>();

  // Remap outputs to ref grid
  m_d2p_remapper->remap(true);

  using KT = KokkosTypes<DefaultDevice>;
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using Pack = ekat::Pack<Real,N>;
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;

  // The internal fields are dynamic subfields of the homme states.
  // We need to update the slice they are subviewing
  get_internal_field("v_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("vtheta_dp_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("dp3d_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("w_int_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("phi_int_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("ps_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("Qdp_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.np1_qdp);

  // Convert VTheta_dp->T, store T,uv, and possibly w in FT, FM,
  // compute p_int on ref grid.
  const auto dp_view = get_field_out("pseudo_density").get_view<Pack**>();
  const auto p_mid_view = get_field_out("p_mid").get_view<Pack**>();
  const auto p_int_view = get_field_out("p_int").get_view<Pack**>();
  const auto Q_view   = get_group_out("Q",rgn).m_bundle->get_view<Pack***>();

  const auto T_view  = get_field_out("T_mid").get_view<Pack**>();
  const auto v_view  = get_field_out("horiz_winds").get_view<Pack***>();
  const auto w_view  = get_field_out("w_int").get_view<Pack**>();
  const auto T_prev_view = m_helper_fields.at("FT_ref").get_view<Pack**>();
  const auto V_prev_view = m_helper_fields.at("FM_ref").get_view<Pack***>();

  const auto ncols = m_ref_grid->get_num_local_dofs();
  const auto nlevs = m_ref_grid->get_num_vertical_levels();
  const auto npacks= ekat::PackInfo<N>::num_packs(nlevs);

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(ncols,npacks);

  // If there are other atm procs updating the vertical velocity,
  // then we need to store w_int_old (to compute forcing for w next iteration)
  const bool has_w_forcing = get_field_out("w_int").get_header().get_tracking().get_providers().size()>1;

  // Establish the boundary condition for the TOA
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();

    // Compute p_int and p_mid
    auto dp = ekat::subview(dp_view,icol);
    auto p_mid = ekat::subview(p_mid_view,icol);
    auto p_int = ekat::subview(p_int_view,icol);

    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    team.team_barrier();

    // Convert VTheta_dp->VTheta->Theta->T
    auto T   = ekat::subview(T_view,icol);
    auto v   = ekat::subview(v_view,icol);
    auto w   = ekat::subview(w_view,icol);
    auto qv  = ekat::subview(Q_view,icol,0);

    auto T_prev = ekat::subview(T_prev_view,icol);
    auto V_prev = ekat::subview(V_prev_view,icol);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,npacks),
                         [&](const int ilev) {
      // VTheta_dp->VTheta->Theta->T
      auto& T_val = T(ilev);
      T_val /= dp(ilev);
      T_val = PF::calculate_temperature_from_virtual_temperature(T_val,qv(ilev));
      T_val = PF::calculate_T_from_theta(T_val,p_mid(ilev));

      // Store T, v (and possibly w) at end of the dyn timestep (to back out tendencies later)
      T_prev(ilev) = T_val;
      V_prev(0,ilev) = v(0,ilev);
      V_prev(1,ilev) = v(1,ilev);
      if (has_w_forcing) {
        V_prev(2,ilev) = w(ilev);
      }
    });
  });

  auto s = Homme::Context::singleton().get<Homme::ElementsState>();
  auto tr = Homme::Context::singleton().get<Homme::Tracers>();
  auto vdyn = m_helper_fields.at("v_dyn");
  
  // If ftype==FORCING_2, also set FQ_ref=Q_ref. Next step's Q_dyn will be set
  // as Q_dyn = Q_dyn_old + PD_remap(Q_ref-Q_ref_old)
  const auto ftype = c.get<Homme::SimulationParams>().ftype;
  if (ftype==Homme::ForcingAlg::FORCING_2) {
    m_helper_fields.at("FQ_ref").deep_copy(*get_group_out("Q",rgn).m_bundle);
  }
}

void HommeDynamics::
create_helper_field (const std::string& name,
                     const std::vector<FieldTag>& tags,
                     const std::vector<int>& dims,
                     const std::string& grid) {
  using namespace ekat::units;
  FieldIdentifier id(name,FieldLayout{tags,dims},Units::nondimensional(),grid);

  const auto lt = get_layout_type(id.get_layout().tags());

  // Only request packed field for 3d quantities
  int pack_size = 1;
  if (lt==LayoutType::Scalar3D || lt==LayoutType::Vector3D || lt==LayoutType::Tensor3D) {
    pack_size = HOMMEXX_PACK_SIZE;
  }

  field_type f(id);
  f.get_header().get_alloc_properties().request_allocation<Real>(pack_size);
  f.allocate_view();
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[name] = f;
}

void HommeDynamics::init_homme_views () {

  const auto& c = Homme::Context::singleton();
  auto& state  = c.get<Homme::ElementsState>();
  auto& tracers = c.get<Homme::Tracers>();
  auto& forcing = c.get<Homme::ElementsForcing>();

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int QTL  = HOMMEXX_Q_NUM_TIME_LEVELS;
  constexpr int QSZ  = HOMMEXX_QSIZE_D;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int NVLI = HOMMEXX_NUM_LEV_P;

  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int qsize = tracers.num_tracers();

  // Print homme's parameters, so user can see whether something wasn't set right.
  if (get_comm().am_i_root()) c.get<Homme::SimulationParams>().print();

  // ------------ Set views in Homme ------------- //
  // Velocity
  auto v_in = m_helper_fields.at("v_dyn").get_view<Homme::Scalar*[NTL][2][NP][NP][NVL]>();
  using v_type = std::remove_reference<decltype(state.m_v)>::type;
  state.m_v = v_type (v_in.data(),nelem);

  // Virtual potential temperature
  auto vtheta_in = m_helper_fields.at("vtheta_dp_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
  using vtheta_type = std::remove_reference<decltype(state.m_vtheta_dp)>::type;
  state.m_vtheta_dp = vtheta_type(vtheta_in.data(),nelem);

  // Geopotential
  auto phi_in = m_helper_fields.at("phi_int_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
  using phi_type = std::remove_reference<decltype(state.m_phinh_i)>::type;
  state.m_phinh_i = phi_type(phi_in.data(),nelem);

  // Vertical velocity
  auto w_in = m_helper_fields.at("w_int_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
  using w_type = std::remove_reference<decltype(state.m_w_i)>::type;
  state.m_w_i = w_type(w_in.data(),nelem);

  // Pseudo-density
  auto dp3d_in = m_helper_fields.at("dp3d_dyn").template get_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
  using dp3d_type = std::remove_reference<decltype(state.m_dp3d)>::type;
  state.m_dp3d = dp3d_type(dp3d_in.data(),nelem);

  // Surface pressure
  auto ps_in = m_helper_fields.at("ps_dyn").template get_view<Real*[NTL][NP][NP]>();
  using ps_type = std::remove_reference<decltype(state.m_ps_v)>::type;
  state.m_ps_v = ps_type(ps_in.data(),nelem);

  // Tracers mixing ratio
  auto q_in = m_helper_fields.at("Q_dyn").template get_view<Homme::Scalar**[NP][NP][NVL]>();
  using q_type = std::remove_reference<decltype(tracers.Q)>::type;
  tracers.Q = q_type(q_in.data(),nelem,qsize);

  // Tracers mass
  auto qdp_in = m_helper_fields.at("Qdp_dyn").template get_view<Homme::Scalar*[QTL][QSZ][NP][NP][NVL]>();
  using qdp_type = std::remove_reference<decltype(tracers.qdp)>::type;
  tracers.qdp = qdp_type(qdp_in.data(),nelem);

  // Tracers forcing
  auto fq_in = m_helper_fields.at("FQ_dyn").template get_view<Homme::Scalar**[NP][NP][NVL]>();
  using fq_type = std::remove_reference<decltype(tracers.fq)>::type;
  tracers.fq = fq_type(fq_in.data(),nelem,qsize);

  // Temperature forcing
  auto ft_in = m_helper_fields.at("FT_dyn").template get_view<Homme::Scalar*[NP][NP][NVL]>();
  using ft_type = std::remove_reference<decltype(forcing.m_ft)>::type;
  forcing.m_ft = ft_type(ft_in.data(),nelem);

  // Momentum forcing
  auto fm_in = m_helper_fields.at("FM_dyn").template get_view<Homme::Scalar**[NP][NP][NVL]>();
  using fm_type = std::remove_reference<decltype(forcing.m_fm)>::type;
  forcing.m_fm = fm_type(fm_in.data(),nelem);
}

void HommeDynamics::restart_homme_state () {
  // Safety checks: internal fields *should* have been restarted (and therefore have a valid timestamp)
  for (auto& f : get_internal_fields()) {
    auto ts = f.get_header().get_tracking().get_time_stamp();
    EKAT_REQUIRE_MSG(ts.is_valid(),
        "Error! Found HommeDynamics internal field not restarted.\n"
        "  - field name: " + f.get_header().get_identifier().name() + "\n");
  }

  const auto& dgn = m_dyn_grid->name();
  const auto& rgn = m_ref_grid->name();

  // All internal fields should have been read from restart file.
  // We need to remap Q_dyn, v_dyn, w_dyn, T_dyn back to ref grid,
  // to handle the backing out of the tendencies
  // Note: Q_ref only in case of ftype!=FORCING_DEBUG.
  // TODO: p2d remapper does not support subfields, so we need to create temps
  //       to remap v_dyn and w_dyn separately, then copy into v3d_prev.
  const auto& c = Homme::Context::singleton();
  auto& params = c.get<Homme::SimulationParams>();

  constexpr int NGP = HOMMEXX_NP;
  const int nlevs = m_ref_grid->get_num_vertical_levels();
  const int ncols = m_ref_grid->get_num_local_dofs();
  const int nelem = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);
  const int qsize = params.qsize;

  // NOTE: when restarting stuff like T_prev, and other "previous steps" quantities that HommeDynamics
  //       uses for tendencies calculation, we need to compute them in the *exact same way* as they
  //       were computed during the original simulation (in homme_post_process).
  //       E.g., we read vtheta_dp(dyn) from restart file, and need to recompute T_prev. Inside
  //       homme_post_process, we use qv(ref), but that's the qv obtained by remapping qv(dyn)
  //       to ref grid *right after homme ran*. Here, we cannot use qv(ref) as read from restart
  //       file, since that's qv(ref) *at the end of the timestep* in the original simulation.
  //       Therefore, we need to remap the end-of-homme-step qv from dyn to ref grid, and use that one.
  //       Another field we need is dp3d(ref), but Homme *CHECKS* that no other atm proc updates
  //       dp3d(ref), so the value read from restart file *coincides* with the value at the end
  //       of the last Homme run. So we can safely recompute pressure using dp(ref) as read from restart.
  update_pressure();

  // Copy all restarted dyn states on all timelevels.
  copy_dyn_states_to_all_timelevels ();

  // Restart end-of-homme-step Q as Qdp/dp. That's what Homme does at the end of the timestep,
  // and by writing/loading only Qdp in the restart file, we save space.
  auto Qdp_dyn_view = get_internal_field("Qdp_dyn",dgn).get_view<Real*****>();
  auto Q_dyn_view = m_helper_fields.at("Q_dyn").get_view<Real*****>();
  auto dp_dyn_view = get_internal_field("dp3d_dyn",dgn).get_view<Real****>();
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*qsize*NGP*NGP*nlevs),
                       KOKKOS_LAMBDA (const int idx) {
    const int ie =  idx / (qsize*NGP*NGP*nlevs);
    const int iq = (idx / (NGP*NGP*nlevs)) % qsize;
    const int ip = (idx / (NGP*nlevs)) % NGP;
    const int jp = (idx / nlevs) % NGP;
    const int k  =  idx % nlevs;
    Q_dyn_view(ie,iq,ip,jp,k) = Qdp_dyn_view(ie,iq,ip,jp,k) / dp_dyn_view(ie,ip,jp,k);
  });

  // TODO: if/when PD remapper supports remapping directly to/from subfields,
  //       you could remap v_dyn and w_int into the proper slice of FM_ref,
  //       and do away with the uv_prev and w_prev helper fields.
  // NOTE: actually, that might be true for uv_prev, but w_int is defined on
  //       interfaces, so you can't remap directly into FM, which is defined
  //       on midpoints.
  using namespace ShortFieldTagsNames;
  create_helper_field("uv_prev",{COL,CMP,LEV},{ncols,2,nlevs},rgn);
  create_helper_field("w_prev", {COL,    LEV},{ncols,  nlevs+1},rgn);

  m_ic_remapper->registration_begins();
  m_ic_remapper->register_field(m_helper_fields.at("FT_ref"),m_helper_fields.at("vtheta_dp_dyn"));
  m_ic_remapper->register_field(m_helper_fields.at("uv_prev"),m_helper_fields.at("v_dyn"));
  m_ic_remapper->register_field(m_helper_fields.at("w_prev"),m_helper_fields.at("w_int_dyn"));
  auto qv_prev_ref = std::make_shared<Field<Real>>();
  auto Q_dyn = m_helper_fields.at("Q_dyn");
  if (params.ftype==Homme::ForcingAlg::FORCING_2) {
    // Recall, we store Q_old in FQ_ref, and do FQ_ref = Q_new - FQ_ref during pre-process
    // Q_old is the tracers at the end of last step, which we can recompute by remapping
    // Q_dyn, which was part of the restart
    auto dQ = m_helper_fields.at("FQ_ref");  
    m_ic_remapper->register_field(dQ,Q_dyn);

    // Grab qv_ref_old from dQ
    *qv_prev_ref = dQ.get_component(0);
  } else {
    // NOTE: we need the end-of-homme-step qv on the ref grid, to do the Theta->T conversion
    //       to compute T_prev *in the same way as homme_post_process* did in the original run.
    //       If PD remapper supported subfields, we would not need qv_prev_dyn, and could use
    //       the proper subfield of Q_dyn instead.
    create_helper_field("qv_prev_ref",{COL,LEV},{ncols,nlevs},rgn);
    create_helper_field("qv_prev_dyn",{EL,GP,GP,LEV},{nelem,NGP,NGP,nlevs},dgn);
    m_ic_remapper->register_field(m_helper_fields.at("qv_prev_ref"),m_helper_fields.at("qv_prev_dyn"));

    // Copy qv from the qsize-sized Q_dyn array, to the individual-tracer field qv_prev_dyn.
    auto qv_prev_dyn = m_helper_fields.at("qv_prev_dyn");
    qv_prev_dyn.deep_copy(Q_dyn.get_component(0));

    *qv_prev_ref = m_helper_fields.at("qv_prev_ref");
  }
  m_ic_remapper->registration_ends();
  m_ic_remapper->remap(/*forward = */false);
  m_ic_remapper = nullptr; // Can clean up the IC remapper now.

  // Copy uv_prev and w_prev into FM_ref. Also, FT_ref contains vtheta_dp,
  // so convert it to actual temperature
  using KT = KokkosTypes<DefaultDevice>;
  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto policy = ESU::get_default_team_policy(ncols,nlevs);
  const bool has_w_forcing = get_field_out("w_int").get_header().get_tracking().get_providers().size()>1;

  auto uv_view     = m_helper_fields.at("uv_prev").get_view<Real***>();
  auto w_view      = m_helper_fields.at("w_prev").get_view<Real**>();
  auto V_prev_view = m_helper_fields.at("FM_ref").get_view<Real***>();
  auto T_prev_view = m_helper_fields.at("FT_ref").get_view<Real**>();
  auto dp_view     = get_field_in("pseudo_density",rgn).get_view<const Real**>();
  auto p_mid_view  = get_field_out("p_mid").get_view<Real**>();
  auto qv_view     = qv_prev_ref->get_view<Real**>();

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team){
    const int icol = team.league_rank();

    auto p_mid = ekat::subview(p_mid_view,icol);
    auto qv    = ekat::subview(qv_view,icol);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,nlevs),
                         [&](const int& ilev) {
      // Init v_prev from uv and, possibly, w
      V_prev_view(icol,0,ilev) = uv_view(icol,0,ilev);
      V_prev_view(icol,1,ilev) = uv_view(icol,1,ilev);
      if (has_w_forcing) {
        V_prev_view(icol,2,ilev) = w_view (icol,  ilev);
      } else {
        V_prev_view(icol,2,ilev) = 0;
      }

      // T_prev as of now contains vtheta_dp. Convert to temperature
      auto& T_prev = T_prev_view(icol,ilev);
      T_prev /= dp_view(icol,ilev);
      T_prev = PF::calculate_temperature_from_virtual_temperature(T_prev,qv(ilev));
      T_prev = PF::calculate_T_from_theta(T_prev,p_mid(ilev));
    });
  });
  Kokkos::fence();

  // We can now erase the uv_prev and w_prev fields.
  // Erase also qv_prev_ref/qv_prev_dyn (if we created them).
  m_helper_fields.erase("uv_prev");
  m_helper_fields.erase("w_prev");
  if (m_helper_fields.find("qv_prev_ref")!=m_helper_fields.end()) {
    m_helper_fields.erase("qv_prev_ref");
    m_helper_fields.erase("qv_prev_dyn");
  }
  qv_prev_ref = nullptr;
}

void HommeDynamics::initialize_homme_state () {
  const auto& rgn = m_ref_grid->name();

  // Import IC from ref grid to dyn grid
  // NOTE: if/when PD remapper supports remapping directly to/from subfields,
  //       you can use get_internal_field (which have a single time slice) rather than
  //       the helper fields (which have NTL time slices).
  m_ic_remapper->registration_begins();
  m_ic_remapper->register_field(get_field_out("w_int",rgn),m_helper_fields.at("w_int_dyn"));
  m_ic_remapper->register_field(get_field_out("phi_int",rgn),m_helper_fields.at("phi_int_dyn"));
  m_ic_remapper->register_field(get_field_out("horiz_winds",rgn),m_helper_fields.at("v_dyn"));
  m_ic_remapper->register_field(get_field_out("pseudo_density",rgn),m_helper_fields.at("dp3d_dyn"));
  m_ic_remapper->register_field(get_field_out("ps",rgn),m_helper_fields.at("ps_dyn"));
  m_ic_remapper->register_field(get_field_out("T_mid",rgn),m_helper_fields.at("vtheta_dp_dyn"));
  m_ic_remapper->register_field(*get_group_out("Q",rgn).m_bundle,m_helper_fields.at("Q_dyn"));
  m_ic_remapper->registration_ends();
  m_ic_remapper->remap(true);

  // Some types
  using KT = KokkosTypes<DefaultDevice>;
  using Pack = ekat::Pack<Real,HOMMEXX_PACK_SIZE>;
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;
  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;

  const auto& c = Homme::Context::singleton();
  auto& params = c.get<Homme::SimulationParams>();

  // Extents
  constexpr int NGP  = HOMMEXX_NP;
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int nlevs = m_dyn_grid->get_num_vertical_levels();
  const int qsize = params.qsize;
  const int npacks_mid = ekat::PackInfo<HOMMEXX_PACK_SIZE>::num_packs(nlevs);
  const int npacks_int = ekat::PackInfo<HOMMEXX_PACK_SIZE>::num_packs(nlevs+1);

  // Homme states
  const auto dp_view  = m_helper_fields.at("dp3d_dyn").get_view<Pack*****>();
  const auto vth_view = m_helper_fields.at("vtheta_dp_dyn").get_view<Pack*****>();
  const auto Q_view   = m_helper_fields.at("Q_dyn").get_view<Pack*****>();

  // State time slices
  const int n0  = c.get<Homme::TimeLevel>().n0;
  const int n0_qdp  = c.get<Homme::TimeLevel>().n0_qdp;

  // ps0 is the TOM boundary condition when we compute p_int(z)=integral_top^z dp(s)ds
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(nelem*NGP*NGP,npacks_mid);

  // Need two temporaries, for pi_mid and pi_int
  ekat::WorkspaceManager<Pack,DefaultDevice> wsm(npacks_int,2,policy);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int ie  =  team.league_rank() / (NGP*NGP);
    const int igp = (team.league_rank() / NGP) % NGP;
    const int jgp =  team.league_rank() % NGP;

    // Compute p_mid
    auto ws = wsm.get_workspace(team);
    auto p_int = ws.take("p_int");
    auto p_mid = ws.take("p_mid");

    auto dp = ekat::subview(dp_view,ie,n0,igp,jgp);

    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    team.team_barrier();
    
    // Convert T->Theta->VTheta->VTheta*dp in place
    auto T      = ekat::subview(vth_view,ie,n0,igp,jgp);
    auto vTh_dp = ekat::subview(vth_view,ie,n0,igp,jgp);
    auto qv     = ekat::subview(Q_view,ie,0,igp,jgp);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,npacks_mid),
                         [&](const int ilev) {
      const auto th = PF::calculate_theta_from_T(T(ilev),p_mid(ilev));
      vTh_dp(ilev) = PF::calculate_virtual_temperature(th,qv(ilev))*dp(ilev);
    });

    // Release the scratch mem
    ws.release(p_int);
    ws.release(p_mid);
  });

  // Update internal fields time stamp
  for (const auto& it : get_internal_fields()) {
    // Unfortunately, get_internal_fields() returns a list of const fields,
    // so grab the name and grid name, then call get_internal_field(name,grid)
    // It's a bit clunky, but not that bad
    const auto& name = it.get_header().get_identifier().name();
    const auto& grid = it.get_header().get_identifier().get_grid_name();
    auto& f = get_internal_field(name,grid);
    f.get_header().get_tracking().update_time_stamp(timestamp());
  }

  // Forcings are computed as some version of "value coming in from AD
  // minus value at the end of last HommeDynamics run".
  // At the first time step, we don't have a value at the end of last
  // HommeDynamics run, so init with the initial conditions.
  // NOTE: for FM, we can't deep copy w_int, since w_int and FM_w
  //       have different number of levels. For u,v, we could, but
  //       we cannot (11/2021) subview 2 slices of FM together, so
  //       we'd need to also subview horiz_winds. Since we
  const auto ncols = m_ref_grid->get_num_local_dofs();
  if (params.ftype==Homme::ForcingAlg::FORCING_2) {
    m_helper_fields.at("FQ_ref").deep_copy(*get_group_out("Q",rgn).m_bundle);
  }
  m_helper_fields.at("FT_ref").deep_copy(get_field_in("T_mid",rgn));
  auto FM_ref = m_helper_fields.at("FM_ref").get_view<Real***>();
  auto horiz_winds = get_field_out("horiz_winds",rgn).get_view<Real***>();
  auto w_int = get_field_out("w_int",rgn).get_view<Real**>();
  const bool has_w_forcing = get_field_out("w_int").get_header().get_tracking().get_providers().size()>1;
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,ncols*nlevs),
                       KOKKOS_LAMBDA (const int idx) {
    const int icol = idx / nlevs;
    const int ilev = idx % nlevs;

    FM_ref(icol,0,ilev) = horiz_winds(icol,0,ilev);
    FM_ref(icol,1,ilev) = horiz_winds(icol,1,ilev);
    if (has_w_forcing) {
      FM_ref(icol,2,ilev) = w_int(icol,ilev);
    } else {
      // Unfortunately, Homme *does* use FM_w even in hydrostatic mode (missing ifdef).
      // Later, it computes diags on w, so if FM_w contains NaN's, repro sum in Homme
      // will crap out. To prevent that, set FM_w=0
      FM_ref(icol,2,ilev) = 0;
    }
  });

  // For initial runs, it's easier to prescribe IC for Q, and compute Qdp = Q*dp
  auto& tracers = c.get<Homme::Tracers>();
  const auto qdp = tracers.qdp;
  const auto q   = tracers.Q;
  const auto dp  = c.get<Homme::ElementsState>().m_dp3d;
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*qsize*NGP*NGP*npacks_mid),
                       KOKKOS_LAMBDA (const int idx) {
    const int ie =  idx / (qsize*NGP*NGP*npacks_mid);
    const int iq = (idx / (NGP*NGP*npacks_mid)) % qsize;
    const int ip = (idx / (NGP*npacks_mid)) % NGP;
    const int jp = (idx / npacks_mid) % NGP;
    const int k  = idx % npacks_mid;

    qdp(ie,n0_qdp,iq,ip,jp,k) = q(ie,iq,ip,jp,k) * dp(ie,n0,ip,jp,k);
  });

  // Initialize p_mid/p_int
  update_pressure ();

  // Copy IC states on all timelevel slices
  copy_dyn_states_to_all_timelevels ();

  // Can clean up the IC remapper now.
  m_ic_remapper = nullptr;
}
// =========================================================================================
void HommeDynamics::
copy_dyn_states_to_all_timelevels () {
  const auto& c = Homme::Context::singleton();

  // State time slices
  const int nm1     = c.get<Homme::TimeLevel>().nm1;
  const int n0      = c.get<Homme::TimeLevel>().n0;
  const int np1     = c.get<Homme::TimeLevel>().np1;
  const int n0_qdp  = c.get<Homme::TimeLevel>().n0_qdp;
  const int np1_qdp = c.get<Homme::TimeLevel>().np1_qdp;

  // States
  const auto dp3d      = m_helper_fields.at("dp3d_dyn");
  const auto ps        = m_helper_fields.at("ps_dyn");
  const auto v         = m_helper_fields.at("v_dyn");
  const auto w_i       = m_helper_fields.at("w_int_dyn");
  const auto phinh_i   = m_helper_fields.at("phi_int_dyn");
  const auto vtheta_dp = m_helper_fields.at("vtheta_dp_dyn");
  const auto qdp       = m_helper_fields.at("Qdp_dyn");

  // Note: it might be somewhat faster to do a single parallel region,
  //       rather than 13 deep copies. However, this is much clearer
  //       to read, and since this method is called only during
  //       initialization, we prefer a more readable method
  dp3d.subfield(1,nm1).deep_copy(dp3d.subfield(1,n0));
  dp3d.subfield(1,np1).deep_copy(dp3d.subfield(1,n0));
  ps.subfield(1,nm1).deep_copy(ps.subfield(1,n0));
  ps.subfield(1,np1).deep_copy(ps.subfield(1,n0));
  v.subfield(1,nm1).deep_copy(v.subfield(1,n0));
  v.subfield(1,np1).deep_copy(v.subfield(1,n0));
  w_i.subfield(1,nm1).deep_copy(w_i.subfield(1,n0));
  w_i.subfield(1,np1).deep_copy(w_i.subfield(1,n0));
  phinh_i.subfield(1,nm1).deep_copy(phinh_i.subfield(1,n0));
  phinh_i.subfield(1,np1).deep_copy(phinh_i.subfield(1,n0));
  vtheta_dp.subfield(1,nm1).deep_copy(vtheta_dp.subfield(1,n0));
  vtheta_dp.subfield(1,np1).deep_copy(vtheta_dp.subfield(1,n0));
  qdp.subfield(1,np1_qdp).deep_copy(qdp.subfield(1,n0_qdp));
}
// =========================================================================================
void HommeDynamics::update_pressure() {
  using KT = KokkosTypes<DefaultDevice>;
  using Pack = ekat::Pack<Real,HOMMEXX_PACK_SIZE>;
  using ColOps = ColumnOps<DefaultDevice,Real>;

  const auto ncols = m_ref_grid->get_num_local_dofs();
  const auto nlevs = m_ref_grid->get_num_vertical_levels();
  const auto npacks= ekat::PackInfo<HOMMEXX_PACK_SIZE>::num_packs(nlevs);

  const auto& c = Homme::Context::singleton();
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

  const auto dp_view    = get_field_out("pseudo_density").get_view<Pack**>();
  const auto p_int_view = get_field_out("p_int").get_view<Pack**>();
  const auto p_mid_view = get_field_out("p_mid").get_view<Pack**>();

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(ncols,npacks);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();

    auto dp = ekat::subview(dp_view,icol);
    auto p_mid = ekat::subview(p_mid_view,icol);
    auto p_int = ekat::subview(p_int_view,icol);

    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
  });
}
// =========================================================================================
void HommeDynamics::
check_computed_fields_impl () {
    // Note: We are seeing near epsilon negative values in a handful of places,
    // The strategy is to 
    // 1. First check that no values are sufficiently negative as to require an error.
    //    i.e. below some tolerance.
    // 2. Clip all negative values to zero.
    // TODO: Construct a more robust check that compares the value of Q against an
    // average value or maximum value over each column.  That way we can use a relative
    // error as our threshold, rather than an arbitrary tolerance.
   
    // Grab the pointer to the tracer group. 
    const auto& rgn = m_ref_grid->name();
          auto& Q   = *get_group_out("Q",rgn).m_bundle;
    // Create a local copy of a lower bound check to ensure we are not encountering truly
    // bad negative tracer values.
    Real tol = -1e-20;
    auto lower_bound_check = std::make_shared<FieldLowerBoundCheck<Real>>(tol);
    lower_bound_check->check(Q);
    // Now repair negative tracers using a lower bounds check at 0.0
    auto lower_bound_repair = std::make_shared<FieldLowerBoundCheck<Real>>(0.0);
    lower_bound_repair->repair(Q);
  
}
// =========================================================================================

} // namespace scream
