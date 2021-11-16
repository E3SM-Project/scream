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
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/util/ekat_units.hpp"
#include "utilities/SubviewUtils.hpp"
#include "VerticalRemapManager.hpp"
#include "HyperviscosityFunctor.hpp"
#include "ExecSpaceDefs.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_workspace.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"
#include "mpi/ConnectivityHelpers.hpp"
#include "Types.hpp"
#include "utilities/MathUtils.hpp"

// Scream includes
#include "dynamics/homme/physics_dynamics_remapper.hpp"
#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/util//scream_column_ops.hpp"
#include "share/field/field_property_checks/field_lower_bound_check.hpp"

// Ekat includes
#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos//ekat_subview_utils.hpp"

// #include <iostream>
// #include <iomanip>

namespace scream
{

HommeDynamics::HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
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
  // using FID = FieldIdentifier;
  using FL  = FieldLayout;
  // using FR  = FieldRequest;

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int QTL  = HOMMEXX_Q_NUM_TIME_LEVELS;

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
          - ftype==FORCING_DEBUG(0): FQ_dyn =( PD(Q_phys_new) - Q_dyn_old )/dt (a tendency value)
          - ftype==FORCING_2(2): FQ_dyn = Q_dyn_old + ( PD(Q_phys_new - Q_phys_old) ) (a new state value)
     Here, PD(x) is quantity x remapped from physics to dynamics grid.
     So for BFB restart, we need so store T_phys_old, v_phys_old, and Q_dyn_old.
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
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
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
  // NOTE: we start with the slice index set to n0 (rather than np1),
  //       so that if this is a restarted run, the IO class will read into the
  //       "previous" timestep slice (i.e., n0).

  const auto& tl = Homme::Context::singleton().get<Homme::TimeLevel>();
  add_internal_field (m_helper_fields.at("v_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("vtheta_dp_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("dp3d_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("w_int_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("phi_int_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("ps_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("Qdp_dyn").subfield(1,tl.n0_qdp,true));
  // SubviewInfo state_sv_info(1,tl.n0,NTL,true);
  // FR v_dyn         (FID("v_dyn",        FL({EL,CMP,GP,GP,LEV}, {nelem,2,NP,NP,nlev_mid}),m/s,     dgn),"state_v_dyn"        ,state_sv_info);
  // FR vtheta_dp_dyn (FID("vtheta_dp_dyn",FL({EL,    GP,GP,LEV}, {nelem,  NP,NP,nlev_mid}),K*Pa,    dgn),"state_vtheta_dp_dyn",state_sv_info);
  // FR dp3d_dyn      (FID("dp3d_dyn",     FL({EL,    GP,GP,LEV}, {nelem,  NP,NP,nlev_mid}),Pa,      dgn),"state_dp3d_dyn"     ,state_sv_info);
  // FR w_int_dyn     (FID("w_int_dyn",    FL({EL,    GP,GP,ILEV},{nelem,  NP,NP,nlev_int}),m/s,     dgn),"state_w_int_dyn"    ,state_sv_info);
  // FR phi_int_dyn   (FID("phi_int_dyn",  FL({EL,    GP,GP,ILEV},{nelem,  NP,NP,nlev_int}),Pa/rho,  dgn),"state_phi_int_dyn"  ,state_sv_info);
  // FR ps_dyn        (FID("ps_dyn",       FL({EL,    GP,GP},     {nelem,  NP,NP         }),Pa,      dgn),"state_ps_dyn"       ,state_sv_info);
  // add_field<Internal>(v_dyn);
  // add_field<Internal>(vtheta_dp_dyn);
  // add_field<Internal>(dp3d_dyn);
  // add_field<Internal>(w_int_dyn);
  // add_field<Internal>(phi_int_dyn);
  // add_field<Internal>(ps_dyn);

  // Import the whole tracers group from the ref grid to the dyn grid.
  // NOTE: this group is in the FM to allow for BFB restarts
  // add_group<Internal>("tracers_mass_dyn",m_dyn_grid->name(),N, Bundling::Required,
  //                     DerivationType::Import, "tracers", rgn);

  // Dynamics backs out tendencies from the states, and passes those to Homme.
  // After Homme completes, we remap the updates state to the ref grid.
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

  // if (run_type==RunType::Restarted) {
  //   // The PD remapper, when remapping p->d for dyn states, remaps the phys field into the n0
  //   // time slice of the dyn states.
  //   auto vd = get_internal_field("v_dyn");
  //   vd.sync_to_host();
  //   auto vdvh = vd.get_view<Real******,Host>();
  //   auto& tl = Homme::Context::singleton().get<Homme::TimeLevel>();
  //   auto un0 = ekat::subview(vdvh,0,tl.n0,0,0,0);
  //   auto unp1 = ekat::subview(vdvh,0,tl.np1,0,0,0);
  //   std::cout << "un0:";
  //   for (int i=0; i<m_dyn_grid->get_num_vertical_levels(); ++i) {
  //     std::cout << " " << un0[i];
  //   }
  //   std::cout << "\nunp1:";
  //   for (int i=0; i<m_dyn_grid->get_num_vertical_levels(); ++i) {
  //     std::cout << " " << un0[i];
  //   }
  //   std::cout << "\n";
  // }

  // Use common/shorter names for tracers.
  alias_group_in  ("tracers",rgn,"Q");
  alias_group_out ("tracers",rgn,"Q");

  // ------ Sanity checks ------- //

  // Nobody should claim to be a provider for dp, w_i.
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
  // Used to back out tendencies FT and FM
  create_helper_field("T_prev",{COL,LEV},    {ncols,  nlevs},rgn);
  create_helper_field("v3d_prev",{COL,CMP,LEV},{ncols,3,nlevs},rgn);

  // Setup the p2d and d2p remappers
  m_p2d_remapper->registration_begins();
  m_d2p_remapper->registration_begins();

  // ftype==FORCING_DEBUG:
  //  1) remap Q_rgn->Q_dyn
  //  2) compute FQ_dyn=(Q_dyn-Q_dyn_old)/dt
  // ftype!=FORCING_DEBUG:
  //  1) compute dQ_rgn=Q_rgn-Q_rgn_old
  //  2) remap dQ_rgn->dQ_dyn
  //  3) Q_dyn=Q_dyn_old+dQ_dyn
  // Note: calling 'register_field_from_src' also *creates* the tgt field
  const auto ftype = params.ftype;
  if (ftype==Homme::ForcingAlg::FORCING_2) {
    // Need a tmp for dQ_rgn
    using namespace ShortFieldTagsNames;
    create_helper_field("dQ_ref",{COL,CMP,LEV},{ncols,qsize,nlevs},rgn);
    m_p2d_remapper->register_field(m_helper_fields.at("dQ_ref"),m_helper_fields.at("FQ_dyn"));
  } else {
    // Can remap Q directly into FQ, tendency computed in pre_process step
    m_p2d_remapper->register_field(*get_group_out("Q",rgn).m_bundle,m_helper_fields.at("FQ_dyn"));
  }
  m_p2d_remapper->register_field(m_helper_fields.at("T_prev"),m_helper_fields.at("FT_dyn"));
  m_p2d_remapper->register_field(m_helper_fields.at("v3d_prev"),m_helper_fields.at("FM_dyn"));

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
  process_initial_conditions (run_type);
  // auto s = Homme::Context::singleton().get<Homme::ElementsState>();
  // auto nm1 = Homme::Context::singleton().get<Homme::TimeLevel>().nm1;
  // auto n0  = Homme::Context::singleton().get<Homme::TimeLevel>().n0;
  // auto np1 = Homme::Context::singleton().get<Homme::TimeLevel>().np1;
  // std::cout << "||v_init(nm1)|| = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,nm1)) << "\n";
  // std::cout << "||v_init(n0) || = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,n0)) << "\n";
  // std::cout << "||v_init(np1)|| = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,np1)) << "\n";

  // Update p_int and p_mid. Other models might import these values from
  // SurfaceCoupling during initialization. They will be overwritten
  // during postprocessing.
  update_pressure();

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
    // Prepare inputs for homme
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
    auto s = Homme::Context::singleton().get<Homme::ElementsState>();
    auto tr = Homme::Context::singleton().get<Homme::Tracers>();
    auto nm1 = Homme::Context::singleton().get<Homme::TimeLevel>().nm1;
    auto n0  = Homme::Context::singleton().get<Homme::TimeLevel>().n0;
    auto np1 = Homme::Context::singleton().get<Homme::TimeLevel>().np1;
    auto n0_qdp  = Homme::Context::singleton().get<Homme::TimeLevel>().n0_qdp;
    auto np1_qdp = Homme::Context::singleton().get<Homme::TimeLevel>().np1_qdp;
    std::cout << "||v(" << nm1 << ")|| = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,nm1)) << "\n";
    std::cout << "||v(" << n0 << ") || = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,n0)) << "\n";
    std::cout << "||v(" << np1 << ")|| = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,np1)) << "\n";
    std::cout << "||tr(" << n0_qdp << ") || = " << Homme::frobenius_norm(Homme::subview(tr.qdp,0,n0_qdp,0)) << "\n";
    std::cout << "||tr(" << np1_qdp << ")|| = " << Homme::frobenius_norm(Homme::subview(tr.qdp,0,np1_qdp,0)) << "\n";
  //   auto t = timestamp();
  //   t += dt;
  // if (t.get_num_steps()==1) {
    // std::cout << "saved Qdp\n";
    // auto Qdp_homme = Homme::Context::singleton().get<Homme::Tracers>().qdp;
    // const int nelem = Qdp_homme.extent(0);
    // const int NGP = 4;
    // const int nlevs = m_dyn_grid->get_num_vertical_levels();
    // const int np1 = Homme::Context::singleton().get<Homme::TimeLevel>().np1_qdp;
    //     std::cout << std::setprecision(17);
    //     for (int ie=0; ie<nelem; ++ie) {
    //       for (int ip=0; ip<NGP; ++ip) {
    //         for (int jp=0; jp<NGP; ++jp) {
    //           std::cout << "q_dyn(" << ie << "," << ip << "," << jp << ",:):";
    //           for (int k=0; k<nlevs; ++k) {
    //             std::cout << " " << Qdp_homme(ie,np1,0,ip,jp,k)[0];
    //           }
    //           std::cout << "\n";
    //         }
    //       }
    //     }
  //   std::cout << "saved u dyn\n";
  // const auto& c = Homme::Context::singleton();
  // const int NGP = 4;
  //   const int  nm1 = c.get<Homme::TimeLevel>().nm1;
  //   const int  n0  = c.get<Homme::TimeLevel>().n0;
  //   const int  np1 = c.get<Homme::TimeLevel>().np1;
  //   auto elems = c.get<Homme::ElementsState>();
  //   const int nelem = elems.num_elems();
  //   const int nlevs = m_dyn_grid->get_num_vertical_levels();
  //   auto v_homme = elems.m_v;
  //       std::cout << std::setprecision(17);
  //       for (int ie=0; ie<nelem; ++ie) {
  //         for (int ip=0; ip<NGP; ++ip) {
  //           for (int jp=0; jp<NGP; ++jp) {
  //             std::cout << "u_dyn(" << ie << "," << ip << "," << jp << ",:):";
  //             for (int k=0; k<nlevs; ++k) {
  //               std::cout << " " << v_homme(ie,n0,0,ip,jp,k)[0];
  //             }
  //             std::cout << "\n";
  //           }
  //         }
  //       }
  //       std::cout << "nm1, n0, np1: " << nm1 << ", " << n0 << ", " << np1 << "\n";


    // const auto& context = Homme::Context::singleton();
    // const auto& conn = context.get<Homme::Connectivity>();
    // auto h_c = conn.get_connections<Homme::HostMemSpace>();
    // auto h_ccs = Homme::subview(h_c,23);


  // Edges
  // constexpr auto SOUTH = Homme::ConnectionName::SOUTH ;
  // constexpr auto NORTH = Homme::ConnectionName::NORTH ;
  // constexpr auto WEST  = Homme::ConnectionName::WEST  ;
  // constexpr auto EAST  = Homme::ConnectionName::EAST  ;
  // constexpr auto SWEST  = Homme::ConnectionName::SWEST  ;
  // constexpr auto SEAST  = Homme::ConnectionName::SEAST  ;
  // constexpr auto NWEST  = Homme::ConnectionName::NWEST  ;
  // constexpr auto NEAST  = Homme::ConnectionName::NEAST  ;
  // Homme::ConnectionHelpers helpers;
    
  //   int ic = -1;
  //   auto k = -1;
  //   auto i = 3;
  //   auto j = 3;
  //   // auto lev = 51;
  //   for (Homme::ConnectionName cn : {SOUTH, NORTH, WEST, EAST, SWEST, SEAST, NWEST, NEAST} ) {
  //     auto idx_c = Homme::etoi(cn);
  //     const auto& this_c = h_ccs(idx_c);
  //     bool has_33 = false;
  //     if (this_c.sharing == Homme::etoi(Homme::ConnectionSharing::MISSING)) {
  //       std::cout << "connection " << idx_c << " is missing.\n";
  //       continue;
  //     }
  //     auto kind = Homme::etoi(helpers.CONNECTION_KIND[idx_c]);
  //     for (int ii=0; ii<helpers.CONNECTION_SIZE[kind]; ++ii) {
  //       const auto& l_pts = helpers.CONNECTION_PTS[Homme::etoi(Homme::Direction::FORWARD)][this_c.local.pos];
  //       const auto& gp = l_pts[ii];
  //       if (gp.ip==i && gp.jp==j) {
  //         has_33 = true;
  //         break;
  //       }
  //     }

  //     if (has_33) {
  //       int kk=-1;
  //       const auto& l_pts = helpers.CONNECTION_PTS[Homme::etoi(Homme::Direction::FORWARD)][this_c.local.pos];
  //       for (int ii=0; ii<4; ++ii) {
  //         if (l_pts[ii].ip==i && l_pts[ii].jp==j) {
  //           kk = ii;
  //           break;
  //         }
  //       }

  //       EKAT_REQUIRE_MSG (kk>=0, "WTF!\n");

  //       const auto& r_pts = helpers.CONNECTION_PTS[this_c.direction][this_c.remote.pos];
  //           std::cout << "connection idx: " << idx_c << "\n";
  //           std::cout << "  remote lid: " << this_c.remote.lid << "\n";
  //           std::cout << "  remote (i,j): " << r_pts[kk].ip << ", " << r_pts[kk].jp << "\n";
  //     }

  //   }
    // for (Homme::ConnectionName cn : {SOUTH, NORTH, WEST, EAST} ) {
    //   const auto& this_c = h_ccs(Homme::etoi(cn));
    //   EKAT_REQUIRE_MSG (this_c.local.lid==23, "WHAT?!?\n");
    //   auto l_pts = helpers.CONNECTION_PTS[Homme::etoi(Homme::Direction::FORWARD)][this_c.local.pos];
    //   for (int ii=0; ii<4; ++ii) {
    //     if (l_pts[ii].ip==i && l_pts[ii].jp==j) {
    //       k = ii;
    //       ic = Homme::etoi(cn);
    //       break;
    //     }
    //   }
    //   if (k>=0) { break;}
    // }
    // EKAT_REQUIRE_MSG (ic>=0 && k>=0, "WTF!\n");
    //   const auto& this_c = h_ccs(ic);
    //   // const auto& r = this_c.remote;
    //   // auto dir = this_c.direction;
    //   // auto r_pts = helpers.CONNECTION_PTS[dir][this_c.local.pos];

    //   // auto tgt = Qdp_homme(23,np1,0,i,j,lev)[0];
    //   // auto remote = Qdp_homme(r.lid,np1,0,r_pts[k].ip,r_pts[k].jp,lev)[0];
    //   // std::cout << "qdp(23," << np1 << ",0,3,3,52): " << tgt << "\n";
    //   // std::cout << "qdp(" << r.lid << "," << np1 << ",0," << r_pts[k].ip << "," << r_pts[k].jp << ",52): " << remote << "\n";

    //   auto gids = m_dyn_grid->get_dofs_gids();
    //   auto lids = m_dyn_grid->get_lid_to_idx_map();
    //   int idx = -1;
    //   for (int ilid=0; ilid<lids.extent_int(0); ++ilid) {
    //     auto native = ekat::subview(lids,ilid);
    //     if (native(0)==this_c.local.lid && native(1)==i && native(2)==j) {
    //       idx = ilid;
    //       break;
    //     }
    //   }
    //   EKAT_REQUIRE_MSG(idx>=0, "You serious?\n");
    //   auto gid = gids(idx);
    //   auto gids_p = m_ref_grid->get_dofs_gids();
    //   int idx_p = -1;
    //   for (int ii=0; ii<gids_p.extent_int(0); ++ii) {
    //     if (gids_p(ii)==gid) {
    //       idx_p = ii;
    //       break;
    //     }
    //   }
    //   EKAT_REQUIRE_MSG(idx_p>=0, "Oh no...\n");
    //   std::cout << "dof gid, lid_p: " << gid << ", " << idx_p << "\n";


    // for (int i=0
  // }
  } catch (std::exception& e) {
    EKAT_ERROR_MSG(e.what());
  } catch (...) {
    EKAT_ERROR_MSG("Something went wrong, but we don't know what.\n");
  }
}

void HommeDynamics::finalize_impl (/* what inputs? */)
{
  Homme::Context::singleton().finalize_singleton();
  prim_finalize_f90();
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

    // // Create a helper field to store the Qdp state, and an internal subfield
    // // to view the time slice to be read/written for restarts
    // using namespace ShortFieldTagsNames;
    // constexpr auto NGP = HOMMEXX_NP;
    // constexpr auto QTL = HOMMEXX_Q_NUM_TIME_LEVELS;
    // const auto nelem = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);
    // const auto nlevs = m_dyn_grid->get_num_vertical_levels();
    // const auto& dgn = m_dyn_grid->name();
    // create_helper_field("Qdp",{EL,TL,CMP,GP,GP,LEV},{nelem,QTL,qsize,NGP,NGP,nlevs},dgn);

    // // Note: only qdp_restart will be in the FM.
    // // NOTE: we start with the slice index set to n0_qdp (rather than np1_qdp),
    // //       so that if this is a restarted run, the IO class will read into the
    // //       "previous" timestep slice (i.e., n0_qdp).
    // const auto n0_qdp = c.get<Homme::TimeLevel>().n0_qdp;
    // Field<Real> qdp_restart = m_helper_fields.at("Qdp").subfield(1,n0_qdp,true);
    // add_internal_field(qdp_restart);
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

  // We can compute FT and FM in place, using T_prev and v3d_prev fields
  auto T  = get_field_in("T_mid").get_view<const Pack**>();
  auto v  = get_field_in("horiz_winds").get_view<const Pack***>();
  auto w  = get_field_in("w_int").get_view<const Pack**>();
  auto FT = m_helper_fields.at("T_prev").get_view<Pack**>();
  auto FM = m_helper_fields.at("v3d_prev").get_view<Pack***>();

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
    const auto dQ = m_helper_fields.at("dQ_ref").get_view<Pack***>();
    const int qsize = Q.extent(1);
    Kokkos::parallel_for(Kokkos::RangePolicy<>(0,ncols*qsize*npacks),
                         KOKKOS_LAMBDA (const int idx) {
      const int icol =  idx / (qsize*npacks);
      const int iq   = (idx / (npacks)) % qsize;
      const int ilev = idx % npacks;
      dQ(icol,iq,ilev) = Q(icol,iq,ilev) - dQ(icol,iq,ilev);
    });
  }

  // Remap FT, FM, and Q (or dQ, depending on ftype)
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
  //  ftype=2: qdp = dp*Qnew
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

  // Remap outputs to ref grid
  m_d2p_remapper->remap(true);

  using KT = KokkosTypes<DefaultDevice>;
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using Pack = ekat::Pack<Real,N>;
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto& c = Homme::Context::singleton();
  const auto& tl = c.get<Homme::TimeLevel>();

  // The internal fields are dynamic subfields of the homme states.
  // We need to update the slice they are subviewing
  std::cout << "setting dynamic slices:\n"
    << "nm1: " << tl.nm1 << "\n"
    << "n0: " << tl.n0 << "\n"
    << "np1: " << tl.np1 << "\n"
    << "n0_qdp: " << tl.n0_qdp << "\n"
    << "np1_qdp: " << tl.np1_qdp << "\n";
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
  const auto T_prev_view = m_helper_fields.at("T_prev").get_view<Pack**>();
  const auto V_prev_view = m_helper_fields.at("v3d_prev").get_view<Pack***>();

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
  std::cout << "has_w_forcing: " << has_w_forcing << "\n";
  std::cout << "||uv|| = " << Homme::frobenius_norm(ekat::scalarize(v_view)) << "\n";
  auto vdyn = get_internal_field("v_dyn").get_view<Real*****>();
  std::cout << "||v_dyn|| = " << Homme::frobenius_norm(vdyn) << "\n";
  // for (int ie=0; ie<vdyn.extent_int(0); ++ie) {
  //   for (int idim=0; idim<vdyn.extent_int(1); ++idim) {
  //     for (int ip=0; ip<vdyn.extent_int(2); ++ip) {
  //       for (int jp=0; jp<vdyn.extent_int(3); ++jp) {
  //         for (int k=0; k<vdyn.extent_int(4); ++k) {
  //           printf("v(%d,%d,%d,%d,%d) = %f\n",ie,idim,ip,jp,k,vdyn(ie,idim,ip,jp,k));
  //         }
  //       }
  //     }
  //   }
  // }
  std::cout << "||v_prev|| = " << Homme::frobenius_norm(ekat::scalarize(V_prev_view)) << "\n";

  // Finally, for BFB restarts we need Qdp.
  // NOTE: if we could save just a subview of homme's Qdp (the right time-level, and possibly
  //       the right qsize), this local copy would not be necessary
  // const auto Qdp_save = get_internal_field("Qdp_dyn").get_view<Homme::Scalar*****>();
  // const auto qdp = c.get<Homme::Tracers>().qdp;
  // const int nelem = Qdp_save.extent(0);
  // const int qsize = Qdp_save.extent(1);
  // const int np    = Qdp_save.extent(2);
  // const int np1_qdp = c.get<Homme::TimeLevel>().np1_qdp;
  // Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*qsize*np*np*npacks),
  //                      KOKKOS_LAMBDA (const int idx) {
  //   const int ie =  idx / (qsize*np*np*npacks);
  //   const int iq = (idx / (np*np*npacks)) % qsize;
  //   const int ip = (idx / (np*npacks)) % np;
  //   const int jp = (idx / npacks) % np;
  //   const int k  =  idx % npacks;

  //   Qdp_save(ie,iq,ip,jp,k) = qdp(ie,np1_qdp,iq,ip,jp,k);
  // });
  
  // If ftype==FORCING_2, also set dQ_ref=Q_ref. Next step's Q_dyn will be set
  // as Q_dyn = Q_dyn_old + PD_remap(Q_ref-Q_ref_old)
  const auto ftype = c.get<Homme::SimulationParams>().ftype;
  if (ftype==Homme::ForcingAlg::FORCING_2) {
    auto dQ_ref = m_helper_fields.at("dQ_ref").get_view<Real***>();
    Kokkos::deep_copy(dQ_ref,get_group_out("Q",rgn).m_bundle->get_view<Real***>());
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

void HommeDynamics::process_initial_conditions (const RunType run_type) {
  const auto& c = Homme::Context::singleton();
  auto& params = c.get<Homme::SimulationParams>();

  const auto& dgn = m_dyn_grid->name();
  const auto& rgn = m_ref_grid->name();
  const auto nlevs = m_ref_grid->get_num_vertical_levels();
  if (run_type==RunType::Restarted) {
    // Safety checks: internal fields *should* have been restarted (and therefore have a valid timestamp)
    for (auto& f : get_internal_fields()) {
      auto ts = f.get_header().get_tracking().get_time_stamp();
      EKAT_REQUIRE_MSG(ts.is_valid(), "WHAT?!?\n");
    }
    // for (auto& g : get_internal_groups()) {
    //   auto ts = g.m_bundle->get_header().get_tracking().get_time_stamp();
    //   EKAT_REQUIRE_MSG(ts.is_valid(), "WHAT?!?\n");
    // }
    auto s = Homme::Context::singleton().get<Homme::ElementsState>();
    auto tr = Homme::Context::singleton().get<Homme::Tracers>();
    auto tl = Homme::Context::singleton().get<Homme::TimeLevel>();

    // All internal fields should have been read from restart file.
    // We need to remap Q_dyn, v_dyn, w_dyn, T_dyn back to ref grid,
    // to handle the backing out of the tendencies
    // Note: Q_ref only in case of ftype!=FORCING_DEBUG.
    // TODO: p2d remapper does not support subfields, so we need to create temps
    //       to remap v_dyn and w_dyn separately, then copy into v3d_prev.
    using namespace ShortFieldTagsNames;
    const auto ncols = m_ref_grid->get_num_local_dofs();

    // Copy Homme states from timelevel n0 to nm1,np1, and init p_mid/p_int
    // Do not compute Theta from T_mid
////////////////////////////////////
    //TODO: IMPORTANT! YOU MUST restart homme's timelevel struct! Need to save tl's nm1, n0, np1 to file,
    //      or else create u_dyn_save as the np1-th timeslice of homme's m_v, so only the "current"
    //      time slice is saved.
    // Also, do not copy states (not necessary).
////////////////////////////////////
    std::cout << "||v_init(nm1)|| = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,tl.nm1)) << "\n";
    std::cout << "||v_init(n0) || = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,tl.n0)) << "\n";
    std::cout << "||v_init(np1)|| = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,tl.np1)) << "\n";
    auto vdyn = get_internal_field("v_dyn").get_view<Real*****>();
    std::cout << "||v_dyn_init|| = " << Homme::frobenius_norm(vdyn) << "\n";
    // for (int ie=0; ie<vdyn.extent_int(0); ++ie) {
    //   for (int idim=0; idim<vdyn.extent_int(1); ++idim) {
    //     for (int ip=0; ip<vdyn.extent_int(2); ++ip) {
    //       for (int jp=0; jp<vdyn.extent_int(3); ++jp) {
    //         for (int k=0; k<vdyn.extent_int(4); ++k) {
    //           printf("v(%d,%d,%d,%d,%d) = %f\n",ie,idim,ip,jp,k,vdyn(ie,idim,ip,jp,k));
    //         }
    //       }
    //     }
    //   }
    // }
    copy_states_and_init_pressure (false);
    std::cout << "||v_init(nm1)|| = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,tl.nm1)) << "\n";
    std::cout << "||v_init(n0) || = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,tl.n0)) << "\n";
    std::cout << "||v_init(np1)|| = " << Homme::frobenius_norm(Homme::subview(s.m_v,0,tl.np1)) << "\n";
    std::cout << "||qv_init(n0) || = " << Homme::frobenius_norm(Homme::subview(tr.qdp,0,tl.n0_qdp,0)) << "\n";
    std::cout << "||qv_init(np1)|| = " << Homme::frobenius_norm(Homme::subview(tr.qdp,0,tl.np1_qdp,0)) << "\n";

    create_helper_field("uv_prev",{COL,CMP,LEV},{ncols,2,nlevs},rgn);
    create_helper_field("w_prev", {COL,    LEV},{ncols,  nlevs+1},rgn);
    m_ic_remapper->registration_begins();
    // NOTE: if/when PD remapper supports remapping directly to/from subfields,
    //       you can use get_internal_field (which have a single time slice) rather than
    //       the helper fields (which have NTL time slices).
    m_ic_remapper->register_field(m_helper_fields.at("T_prev"),m_helper_fields.at("vtheta_dp_dyn"));
    m_ic_remapper->register_field(m_helper_fields.at("uv_prev"),m_helper_fields.at("v_dyn"));
    m_ic_remapper->register_field(m_helper_fields.at("w_prev"),m_helper_fields.at("w_int_dyn"));
    if (params.ftype==Homme::ForcingAlg::FORCING_2) {
      // Recall, we store Q_old in dQ_ref, and do dQ_ref = Q_new - dQ_ref during pre-process
      // Q_old is the tracers at the end of last step, which we can recompute by remapping
      // Q_dyn, which was part of the restart
      auto dQ = m_helper_fields.at("dQ_ref");  
      auto Q_dyn = m_helper_fields.at("Q_dyn");
      m_ic_remapper->register_field(dQ,Q_dyn);
    }
    m_ic_remapper->registration_ends();
    m_ic_remapper->remap(/*forward = */false);
    std::cout << "||uv_prev|| = " << Homme::frobenius_norm(m_helper_fields.at("uv_prev").get_view<Real***>()) << "\n";

    // Can clean up the IC remapper now.
    m_ic_remapper = nullptr;

    // Copy uv_prev and w_prev into v3d_prev. Also, T_prev contains vtheta_dp,
    // so convert it to actual temperature

    using KT = KokkosTypes<DefaultDevice>;
    using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
    using PF = PhysicsFunctions<DefaultDevice>;

    const auto policy = ESU::get_thread_range_parallel_scan_team_policy(ncols,nlevs);

    auto uv_view     = m_helper_fields.at("uv_prev").get_view<Real***>();
    auto w_view      = m_helper_fields.at("w_prev").get_view<Real**>();
    auto V_prev_view = m_helper_fields.at("v3d_prev").get_view<Real***>();
    auto T_prev_view = m_helper_fields.at("T_prev").get_view<Real**>();
    auto dp_view     = get_field_in("pseudo_density",rgn).get_view<const Real**>();
    auto p_mid_view  = get_field_out("p_mid").get_view<Real**>();
    auto Q_view      = get_group_out("Q",rgn).m_bundle->get_view<Real***>();
    const bool has_w_forcing = get_field_out("w_int").get_header().get_tracking().get_providers().size()>1;

    Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team){
      const int icol = team.league_rank();

      auto p_mid = ekat::subview(p_mid_view,icol);
      auto qv    = ekat::subview(Q_view,icol,0);
      // ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
      // team.team_barrier();
      // ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
      // team.team_barrier();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,nlevs),
                           [&](const int& ilev) {
        // Init v3d from uv and w
        V_prev_view(icol,0,ilev) = uv_view(icol,0,ilev);
        V_prev_view(icol,1,ilev) = uv_view(icol,1,ilev);
        if (has_w_forcing) {
          V_prev_view(icol,2,ilev) = w_view (icol,  ilev);
        }

        // T_prev as of now contains vtheta_dp. Convert to temperature
        auto& T_prev = T_prev_view(icol,ilev);
        T_prev /= dp_view(icol,ilev);
        T_prev = PF::calculate_temperature_from_virtual_temperature(T_prev,qv(ilev));
        T_prev = PF::calculate_T_from_theta(T_prev,p_mid(ilev));
      });
    });
    Kokkos::fence();
    std::cout << "has_w_forcing: " << has_w_forcing << "\n";
    std::cout << "||v_prev|| = " << Homme::frobenius_norm(V_prev_view) << "\n";

    // We can now erase the uv_prev and w_prev fields
    m_helper_fields.erase("uv_prev");
    m_helper_fields.erase("w_prev");

    // Compute Q from Qdp (that was loaded from file)
    // // Need to copy the Qdp internal group into Qdp in homme
    // // NOTE: if we found a way to load the restart field directly into a time slice
    // //       of Homme's Qdp, we could make the internal Qdp group alias the Qdp view
    // //       in Homme directly. Instead, we use a (nelem,qsize,np,np,nlev) temporary,
    // //       load the restart field in it, and then copy into Homme's qdp
    auto Qdp = get_internal_field("Qdp_dyn",dgn).get_view<Real*****>();
    auto Q_dyn = m_helper_fields.at("Q_dyn").get_view<Real*****>();
    auto dp_dyn = get_internal_field("dp3d_dyn",dgn).get_view<Real****>();
    const auto& tracers = c.get<Homme::Tracers>();
    const int nelem = tracers.num_elems();
    const int qsize = tracers.num_tracers();
    constexpr int NGP = HOMMEXX_NP;
    // const auto& Qdp_homme = c.get<Homme::Tracers>().qdp;
    // const int  n0 = c.get<Homme::TimeLevel>().n0;
    // const int qn0 = c.get<Homme::TimeLevel>().n0_qdp;
    // std::cout << "initializeing homme tracers: qsize=" << qsize << "\n";

    // constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
    const auto npacks = ekat::PackInfo<HOMMEXX_PACK_SIZE>::num_packs(nlevs);
    Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*qsize*NGP*NGP*npacks),
                         KOKKOS_LAMBDA (const int idx) {
      const int ie =  idx / (qsize*NGP*NGP*npacks);
      const int iq = (idx / (NGP*NGP*npacks)) % qsize;
      const int ip = (idx / (NGP*npacks)) % NGP;
      const int jp = (idx / npacks) % NGP;
      const int k  =  idx % npacks;
      // const int ilev = k / HOMMEXX_PACK_SIZE;
      // const int ivec = k % HOMMEXX_PACK_SIZE;
      // TODO: should we copy only into tl.n0_qdp?
      // Qdp(ie,iq,ip,jp,ilev)[ivec] = Qdp_homme(ie,1,iq,ip,jp,ilev)[ivec] = Qdp_loaded(ie,iq,ip,jp,k);
      // printf("(ie,iq,i,j,k) = (%d,%d,%d,%d,%d)\n",ie,iq,ip,jp,k);
      Q_dyn(ie,iq,ip,jp,k) = Qdp(ie,iq,ip,jp,k) / dp_dyn(ie,ip,jp,k);
    });
    // std::cout << "restarted u dyn\n";
    // auto v_homme = c.get<Homme::ElementsState>().m_v;
    //     std::cout << std::setprecision(17);
    //     for (int ie=0; ie<nelem; ++ie) {
    //       for (int ip=0; ip<NGP; ++ip) {
    //         for (int jp=0; jp<NGP; ++jp) {
    //           std::cout << "u_dyn(" << ie << "," << ip << "," << jp << ",:):";
    //           for (int k=0; k<nlevs; ++k) {
    //             std::cout << " " << v_homme(ie,n0,0,ip,jp,k)[0];
    //           }
    //           std::cout << "\n";
    //         }
    //       }
    //     }
    // std::cout << "restarted Qdp\n";
    //     std::cout << std::setprecision(17);
    //     for (int ie=0; ie<nelem; ++ie) {
    //       for (int ip=0; ip<NGP; ++ip) {
    //         for (int jp=0; jp<NGP; ++jp) {
    //           std::cout << "q_dyn(" << ie << "," << ip << "," << jp << ",:):";
    //           for (int k=0; k<nlevs; ++k) {
    //             std::cout << " " << Qdp_homme(ie,0,0,ip,jp,k)[0];
    //           }
    //           std::cout << "\n";
    //         }
    //       }
    //     }
  //   const auto& context = Homme::Context::singleton();
  //   const auto& conn = context.get<Homme::Connectivity>();
  //   auto h_c = conn.get_connections<Homme::HostMemSpace>();
  //   auto h_ccs = Homme::subview(h_c,23);


  // // Edges
  // constexpr auto SOUTH = Homme::ConnectionName::SOUTH ;
  // constexpr auto NORTH = Homme::ConnectionName::NORTH ;
  // constexpr auto WEST  = Homme::ConnectionName::WEST  ;
  // constexpr auto EAST  = Homme::ConnectionName::EAST  ;
  // Homme::ConnectionHelpers helpers;
    
  //   int ic = -1;
  //   auto k = -1;
  //   auto i = 3;
  //   auto j = 3;
  //   auto lev = 51;
  //   for (Homme::ConnectionName cn : {SOUTH, NORTH, WEST, EAST} ) {
  //     const auto& this_c = h_ccs(Homme::etoi(cn));
  //     EKAT_REQUIRE_MSG (this_c.local.lid==23, "WHAT?!?\n");
  //     auto l_pts = helpers.CONNECTION_PTS[Homme::etoi(Homme::Direction::FORWARD)][this_c.local.pos];
  //     for (int ii=0; ii<4; ++ii) {
  //       if (l_pts[ii].ip==i && l_pts[ii].jp==j) {
  //         k = ii;
  //         ic = Homme::etoi(cn);
  //         break;
  //       }
  //     }
  //     if (k>=0) { break;}
  //   }
  //   const int  n0_qdp = c.get<Homme::TimeLevel>().n0_qdp;
  //   EKAT_REQUIRE_MSG (ic>=0 && k>=0, "WTF!\n");
  //     const auto& this_c = h_ccs(ic);
  //     const auto& r = this_c.remote;
  //     auto dir = this_c.direction;
  //     auto r_pts = helpers.CONNECTION_PTS[dir][this_c.local.pos];

  //     auto tgt = Qdp_homme(23,n0_qdp,0,i,j,lev)[0];
  //     auto remote = Qdp_homme(r.lid,n0_qdp,0,r_pts[k].ip,r_pts[k].jp,lev)[0];
  //     std::cout << "qdp(23," << n0_qdp << ",0,3,3,52): " << tgt << "\n";
  //     std::cout << "qdp(" << r.lid << "," << n0_qdp << ",0," << r_pts[k].ip << "," << r_pts[k].jp << ",52): " << remote << "\n";

    return;
  }

  // Import IC from ref grid to dyn grid internal fields
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

  // Spatial extents of views
  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);

  // Initialize p_mid/p_int, and also copy IC states on all timelevel slices
  copy_states_and_init_pressure (true);

  // Init previous timestep quantity to the corresponding IC value
  const auto ncols = m_ref_grid->get_num_local_dofs();
  if (params.ftype==Homme::ForcingAlg::FORCING_2) {
    m_helper_fields.at("dQ_ref").deep_copy(*get_group_out("Q",rgn).m_bundle);
  }
  m_helper_fields.at("T_prev").deep_copy(get_field_in("T_mid",rgn));
  auto v_prev = m_helper_fields.at("v3d_prev").get_view<Real***>();
  auto horiz_winds = get_field_out("horiz_winds",rgn).get_view<Real***>();
  auto w_int = get_field_out("w_int",rgn).get_view<Real**>();
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,ncols*nlevs),
                       KOKKOS_LAMBDA (const int idx) {
    const int icol = idx / nlevs;
    const int ilev = idx % nlevs;
    v_prev(icol,0,ilev) = horiz_winds(icol,0,ilev);
    v_prev(icol,1,ilev) = horiz_winds(icol,1,ilev);
    v_prev(icol,2,ilev) = w_int(icol,ilev);
  });

  // For initial runs, it's easier to prescribe IC for Q,
  // and compute Qdp = Q*dp
  auto& tracers = c.get<Homme::Tracers>();
  const auto qdp = tracers.qdp;
  const auto q   = tracers.Q;
  const auto dp  = c.get<Homme::ElementsState>().m_dp3d;
  const int n0_qdp = c.get<Homme::TimeLevel>().n0_qdp;
  const int n0  = c.get<Homme::TimeLevel>().n0;
  const int qsize = tracers.num_tracers(); 
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*qsize*NP*NP*NVL),
                       KOKKOS_LAMBDA (const int idx) {
    const int ie =  idx / (qsize*NP*NP*NVL);
    const int iq = (idx / (NP*NP*NVL)) % qsize;
    const int ip = (idx / (NP*NVL)) % NP;
    const int jp = (idx / NVL) % NP;
    const int k  = idx % NVL;

    qdp(ie,n0_qdp,iq,ip,jp,k) = q(ie,iq,ip,jp,k) * dp(ie,n0,ip,jp,k);
  });

  // Can clean up the IC remapper now.
  m_ic_remapper = nullptr;
}
// void HommeDynamics::
// store_prev_fields (const view_Nd_type<3>& uv, const view_Nd_type<2>& w) {
//   using KT = KokkosTypes<DefaultDevice>;
//   using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
//   using KT = KokkosTypes<DefaultDevice>;
//   using ColOps = ColumnOps<DefaultDevice,Real>;

//   const auto& rgn = m_ref_grid->name();
//   const auto nlevs = m_ref_grid->get_num_vertical_levels();
//   const auto ncols = m_ref_grid->get_num_local_dofs();

//   const auto dp_view     = get_field_in("pseudo_density",rgn).get_view<Real**>();
//   const auto p_int_view  = get_field_out("p_int").get_view<Real**>();
//   const auto p_mid_view  = get_field_out("p_mid").get_view<Real**>();
//   const auto T_view      = get_field_out("T_mid").get_view<Real**>();
//   const auto T_prev_view = m_helper_fields.at("T_prev").get_view<Real**>();
//   const auto V_prev_view = m_helper_fields.at("v3d_prev").get_view<Real***>();
//   const auto Q_view      = get_group_out("Q",rgn).m_bundle->get_view<Real***>();

//   const auto policy = ESU::get_thread_range_parallel_scan_team_policy(ncols,nlevs);
//   ekat::WorkspaceManager<Real,DefaultDevice> wsm(nlevs+1,2,policy);
//   Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
//     const int& icol = team.league_rank();

//     // Compute p_int and p_mid
//     auto dp = ekat::subview(dp_view,icol);
//     auto p_mid = ekat::subview(p_mid_view,icol);
//     auto p_int = ekat::subview(p_int_view,icol);

//     // Convert VTheta_dp->VTheta->Theta->T
//     auto T   = ekat::subview(T_view,icol);
//     auto qv  = ekat::subview(Q_view,icol,0);

//     auto T_prev = ekat::subview(T_prev_view,icol);
//     auto V_prev = ekat::subview(V_prev_view,icol);

//     Kokkos::parallel_for(Kokkos::TeamThreadRange(team,nlevs),
//                          [&](const int ilev) {

//       // VTheta_dp->VTheta->Theta->T
//       auto& T_val = T_prev(ilev);
//       T_val /= dp(ilev);
//       T_val = PF::calculate_temperature_from_virtual_temperature(T_val,qv(ilev));
//       T_val = PF::calculate_T_from_theta(T_val,p_mid(ilev));

//       // Store T, v (and possibly w) at end of the dyn timestep (to back out tendencies later)
//       T_prev(ilev) = T_val;
//       V_prev(0,ilev) = v(0,ilev);
//       V_prev(1,ilev) = v(1,ilev);
//       if (has_w_forcing) {
//         V_prev(2,ilev) = w(ilev);
//       }
//     });
//   });
// }
// =========================================================================================
void HommeDynamics::
copy_states_and_init_pressure (const bool compute_theta_from_T) {
  // Some types
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using KT = KokkosTypes<DefaultDevice>;
  using Pack = ekat::Pack<Real,N>;
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;

  // Fields extents
  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int NVLI = HOMMEXX_NUM_LEV_P;
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int nlevs = m_dyn_grid->get_num_vertical_levels();

  const auto& c = Homme::Context::singleton();
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

  // Homme states
  const auto dp3d      = m_helper_fields.at("dp3d_dyn").get_view<Pack*****>();
  const auto ps        = m_helper_fields.at("ps_dyn").get_view<Real****>();
  const auto v         = m_helper_fields.at("v_dyn").get_view<Pack******>();
  const auto w_i       = m_helper_fields.at("w_int_dyn").get_view<Pack*****>();
  const auto phinh_i   = m_helper_fields.at("phi_int_dyn").get_view<Pack*****>();
  const auto vtheta_dp = m_helper_fields.at("vtheta_dp_dyn").get_view<Pack*****>();
  const auto qdp       = m_helper_fields.at("Qdp_dyn").get_view<Pack******>();
  const auto Q         = m_helper_fields.at("Q_dyn").get_view<Pack*****>();

  // State time slices
  const int nm1 = c.get<Homme::TimeLevel>().nm1;
  const int n0  = c.get<Homme::TimeLevel>().n0;
  const int np1 = c.get<Homme::TimeLevel>().np1;
  const int n0_qdp  = c.get<Homme::TimeLevel>().n0_qdp;
  const int np1_qdp = c.get<Homme::TimeLevel>().np1_qdp;

  const int qsize = c.get<Homme::Tracers>().num_tracers();

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(nelem*NP*NP,NVL);

  // Need two temporaries, for pi_mid and pi_int
  ekat::WorkspaceManager<Pack,DefaultDevice> wsm(NVLI,2,policy);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int ie  =  team.league_rank() / (NP*NP);
    const int igp = (team.league_rank() / NP) % NP;
    const int jgp =  team.league_rank() % NP;

    auto ws = wsm.get_workspace(team);

    // Compute p_mid
    auto dp = ekat::subview(dp3d,ie,n0,igp,jgp);

    auto p_int = ws.take("p_int");
    auto p_mid = ws.take("p_mid");
    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    team.team_barrier();
    
    auto T      = ekat::subview(vtheta_dp,ie,n0,igp,jgp);
    auto vTh_dp = ekat::subview(vtheta_dp,ie,n0,igp,jgp);
    auto qv     = ekat::subview(Q,ie,0,igp,jgp);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,NVL),
                         [&](const int ilev) {
      if (compute_theta_from_T) {
        // Convert T->Theta->Theta*dp->VTheta*dp
        const auto vthdp = dp(ilev)*PF::calculate_theta_from_T(T(ilev),p_mid(ilev));
        vTh_dp(ilev) = PF::calculate_virtual_temperature(vthdp,qv(ilev));
      }

      // Copy states from the n0 timelevel to all the other ones
      dp3d(ie,nm1,igp,jgp,ilev) = dp3d(ie,n0,igp,jgp,ilev);
      dp3d(ie,np1,igp,jgp,ilev) = dp3d(ie,n0,igp,jgp,ilev);
      
      v(ie,nm1,0,igp,jgp,ilev) = v(ie,n0,0,igp,jgp,ilev);
      v(ie,nm1,1,igp,jgp,ilev) = v(ie,n0,1,igp,jgp,ilev);
      v(ie,np1,0,igp,jgp,ilev) = v(ie,n0,0,igp,jgp,ilev);
      v(ie,np1,1,igp,jgp,ilev) = v(ie,n0,1,igp,jgp,ilev);

      vtheta_dp(ie,nm1,igp,jgp,ilev) = vtheta_dp(ie,n0,igp,jgp,ilev);
      vtheta_dp(ie,np1,igp,jgp,ilev) = vtheta_dp(ie,n0,igp,jgp,ilev);

      w_i(ie,nm1,igp,jgp,ilev) = w_i(ie,n0,igp,jgp,ilev);
      w_i(ie,np1,igp,jgp,ilev) = w_i(ie,n0,igp,jgp,ilev);

      phinh_i(ie,nm1,igp,jgp,ilev) = phinh_i(ie,n0,igp,jgp,ilev);
      phinh_i(ie,np1,igp,jgp,ilev) = phinh_i(ie,n0,igp,jgp,ilev);

      // Note: this has bad access, especially on GPU, but it's an init routine,
      //       so we accept the performance hit
      for (int iq=0; iq<qsize; ++iq) {
        qdp(ie,np1_qdp,iq,igp,jgp,ilev) = qdp(ie,n0_qdp,iq,igp,jgp,ilev);
      }
    });

    // Copy ps (and last interface of w_i and phinh_i)  from the n0 timelevel to all the other ones
    Kokkos::single(Kokkos::PerTeam(team),[&](){
      ps(ie,nm1,igp,jgp) = ps(ie,n0,igp,jgp);
      ps(ie,np1,igp,jgp) = ps(ie,n0,igp,jgp);
      w_i(ie,nm1,igp,jgp,NVLI-1) = w_i(ie,n0,igp,jgp,NVLI-1);
      w_i(ie,np1,igp,jgp,NVLI-1) = w_i(ie,n0,igp,jgp,NVLI-1);

      phinh_i(ie,nm1,igp,jgp,NVLI-1) = phinh_i(ie,n0,igp,jgp,NVLI-1);
      phinh_i(ie,np1,igp,jgp,NVLI-1) = phinh_i(ie,n0,igp,jgp,NVLI-1);
    });

    // Release the scratch mem
    ws.release(p_int);
    ws.release(p_mid);
  });
}
// =========================================================================================
void HommeDynamics::update_pressure() {
  using KT = KokkosTypes<DefaultDevice>;
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using Pack = ekat::Pack<Real,N>;
  using ColOps = ColumnOps<DefaultDevice,Real>;

  const auto ncols = m_ref_grid->get_num_local_dofs();
  const auto nlevs = m_ref_grid->get_num_vertical_levels();
  const auto npacks= ekat::PackInfo<N>::num_packs(nlevs);

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
