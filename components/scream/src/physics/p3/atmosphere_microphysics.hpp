#ifndef SCREAM_P3_MICROPHYSICS_HPP
#define SCREAM_P3_MICROPHYSICS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "physics/p3/p3_functions.hpp"
#include "share/util/scream_common_physics_functions.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the atmosphere microphysics
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
 *  Note: for now, scream is only going to accommodate P3 as microphysics
*/

  using namespace p3;
  using P3F          = Functions<Real, DefaultDevice>;
  using Spack        = typename P3F::Spack;
  using Smask        = typename P3F::Smask;
  using Pack         = ekat::Pack<Real,Spack::n>;
  using PF           = scream::PhysicsFunctions<DefaultDevice>;
  using KT           = ekat::KokkosTypes<DefaultDevice>;
  using WSM          = ekat::WorkspaceManager<Spack, KT::Device>;

  using view_1d  = typename P3F::view_1d<Real>;
  using view_2d  = typename P3F::view_2d<Spack>;
  using view_2d_const  = typename P3F::view_2d<const Spack>;
  using sview_2d = typename KokkosTypes<DefaultDevice>::template view_2d<Real>;

  using uview_1d  = Unmanaged<view_1d>;
  using uview_2d  = Unmanaged<view_2d>;
  using suview_2d = Unmanaged<sview_2d>;

class P3Microphysics : public AtmosphereProcess
{
public:
  // Constructors
  P3Microphysics (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Microphysics"; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  /*--------------------------------------------------------------------------------------------*/
  // Most individual processes have a pre-processing step that constructs needed variables from
  // the set of fields stored in the field manager.  A structure like this defines those operations,
  // which can then be called during run_imple in the main .cpp code.
  // NOTE: the use of the ekat command "set" to copy local values into the variables of interest.
  // This is an important step to avoid having those variables just share pointers to memory.
  // Structure to handle the local generation of data needed by p3_main in run_impl
  struct p3_preamble {
    p3_preamble() = default;
    // Functor for Kokkos loop to pre-process every run step
    KOKKOS_INLINE_FUNCTION
    void operator()(const int icol) const {
      for (int ipack=0;ipack<m_npack;ipack++) {
        // The ipack slice of input variables used more than once
        const Spack& pmid_pack(pmid(icol,ipack));
        const Spack& T_atm_pack(T_atm(icol,ipack));
        const Spack& cld_frac_t_pack(cld_frac_t(icol,ipack));
        // Ensure that the pseudo_density and pmid used by P3 is dry
        pdel_dry(icol,ipack) = PF::calculate_drymmr_from_wetmmr(pseudo_density(icol,ipack),qv(icol,ipack));
        pmid_dry(icol,ipack) = PF::calculate_drymmr_from_wetmmr(pmid_pack,qv(icol,ipack));
        // Exner
        const auto& exner = PF::exner_function(pmid_dry(icol,ipack));
        inv_exner(icol,ipack) = 1.0/exner;
        // Potential temperature
        th_atm(icol,ipack) = PF::calculate_theta_from_T(T_atm_pack,pmid_dry(icol,ipack));
        // DZ
        dz(icol,ipack) = PF::calculate_dz(pdel_dry(icol,ipack), pmid_dry(icol,ipack), T_atm_pack, qv(icol,ipack));
        // Cloud fraction
        // Set minimum cloud fraction - avoids division by zero
        cld_frac_l(icol,ipack) = ekat::max(cld_frac_t_pack,mincld);
        cld_frac_i(icol,ipack) = ekat::max(cld_frac_t_pack,mincld);
        cld_frac_r(icol,ipack) = ekat::max(cld_frac_t_pack,mincld);

        /*----------------------------------------------------------------------------------------------------------------------
         *Wet to dry mixing ratios:
         *-------------------------
         *Since state constituents from the host model (or AD) are  wet mixing ratios and P3 needs
         *these constituents in dry mixing ratios, we convert the wet mixing ratios to dry mixing ratios.

         *NOTE:Function calculate_drymmr_from_wetmmr takes 2 arguments: ( wet mmr and "wet" water vapor mixing ratio)

         *IMPORTANT:Convert "qv wet mmr" to "qv dry mmr" after converting all other constituents to dry mmr as "qv" (as wet mmr)
         * is an input for converting all other constituent5Bs to have dry mmr.

         *----------------------------------------------------------------------------------------------------------------------
         */

        //Since "qv" has a wet mixing ratio, we can use "qv" to compute dry mixing ratios of the following constituents:
        //Units of all constituents below are [kg/kg(dry-air)] for mass and [#/kg(dry-air)] for number
        qc(icol, ipack)      = PF::calculate_drymmr_from_wetmmr(qc(icol,ipack),qv(icol,ipack)); //Cloud liquid mass
        nc(icol, ipack)      = PF::calculate_drymmr_from_wetmmr(nc(icol,ipack),qv(icol,ipack)); //Cloud liquid numbe
        qr(icol, ipack)      = PF::calculate_drymmr_from_wetmmr(qr(icol,ipack),qv(icol,ipack)); //Rain mass
        nr(icol, ipack)      = PF::calculate_drymmr_from_wetmmr(nr(icol,ipack),qv(icol,ipack)); //Rain number
        qi(icol, ipack)      = PF::calculate_drymmr_from_wetmmr(qi(icol,ipack),qv(icol,ipack)); //Cloud ice mass
        ni(icol, ipack)      = PF::calculate_drymmr_from_wetmmr(ni(icol,ipack),qv(icol,ipack)); //Cloud ice number
        qm(icol, ipack)      = PF::calculate_drymmr_from_wetmmr(qm(icol,ipack),qv(icol,ipack)); //Rimmed ice mass
        bm(icol, ipack)      = PF::calculate_drymmr_from_wetmmr(bm(icol,ipack),qv(icol,ipack)); //Rimmed ice number
        //Water vapor from previous time step
        qv_prev(icol, ipack) = PF::calculate_drymmr_from_wetmmr(qv_prev(icol,ipack),qv(icol,ipack));

        // ^^ Ensure that qv is "wet mmr" till this point ^^
        //NOTE: Convert "qv" to dry mmr in the end after converting all other constituents to dry mmr
        qv(icol, ipack)      = PF::calculate_drymmr_from_wetmmr(qv(icol,ipack),qv(icol,ipack));


        // update rain cloud fraction given neighboring levels using max-overlap approach.
        for (int ivec=0;ivec<Spack::n;ivec++)
        {
          // Hard-coded max-overlap cloud fraction calculation.  Cycle through the layers from top to bottom and determine if the rain fraction needs to
          // be updated to match the cloud fraction in the layer above.  It is necessary to calculate the location of the layer directly above this one,
          // labeled ipack_m1 and ivec_m1 respectively.  Note, the top layer has no layer above it, which is why we have the kstr index in the loop.
          Int lev = ipack*Spack::n + ivec;  // Determine the level at this pack/vec location.
          Int ipack_m1 = (lev - 1) / Spack::n;
          Int ivec_m1  = (lev - 1) % Spack::n;
          if (lev != 0) { /* Not applicable at the very top layer */
            cld_frac_r(icol,ipack)[ivec] = cld_frac_t(icol,ipack_m1)[ivec_m1]>cld_frac_r(icol,ipack)[ivec] ?
                                              cld_frac_t(icol,ipack_m1)[ivec_m1] :
                                              cld_frac_r(icol,ipack)[ivec];
          }
        }
        //
      }
    } // operator
    // Local variables
    int m_ncol, m_npack;
    Real mincld = 0.0001;  // TODO: These should be stored somewhere as more universal constants.  Or maybe in the P3 class hpp
    view_2d_const pmid;
    view_2d_const pseudo_density;
    view_2d       T_atm;
    view_2d_const cld_frac_t;
    view_2d       qv;
    view_2d       qc;
    view_2d       nc;
    view_2d       qr;
    view_2d       nr;
    view_2d       qi;
    view_2d       qm;
    view_2d       ni;
    view_2d       bm;
    view_2d       qv_prev;
    view_2d       inv_exner;
    view_2d       th_atm;
    view_2d       cld_frac_l;
    view_2d       cld_frac_i;
    view_2d       cld_frac_r;
    view_2d       dz;
    view_2d       pdel_dry;
    view_2d       pmid_dry;
    // Assigning local variables
    void set_variables(const int ncol, const int npack,
           const view_2d_const& pmid_, const view_2d_const& pseudo_density_, const view_2d& T_atm_,
           const view_2d_const& cld_frac_t_, const view_2d& qv_, const view_2d& qc_,
           const view_2d& nc_, const view_2d& qr_, const view_2d& nr_, const view_2d& qi_,
           const view_2d& qm_, const view_2d& ni_, const view_2d& bm_, const view_2d& qv_prev_,
           const view_2d& inv_exner_, const view_2d& th_atm_, const view_2d& cld_frac_l_,
           const view_2d& cld_frac_i_, const view_2d& cld_frac_r_, const view_2d& dz_, const view_2d& pdel_dry_, const view_2d& pmid_dry_
           )
    {
      m_ncol = ncol;
      m_npack = npack;
      // IN
      pmid           = pmid_;
      pseudo_density = pseudo_density_;
      T_atm          = T_atm_;
      cld_frac_t     = cld_frac_t_;
      // OUT
      qv             = qv_;
      qc             = qc_;
      nc             = nc_;
      qr             = qr_;
      nr             = nr_;
      qi             = qi_;
      qm             = qm_;
      ni             = ni_;
      bm             = bm_;
      qv_prev        = qv_prev_;
      inv_exner = inv_exner_;
      th_atm = th_atm_;
      cld_frac_l = cld_frac_l_;
      cld_frac_i = cld_frac_i_;
      cld_frac_r = cld_frac_r_;
      dz = dz_;
      pdel_dry = pdel_dry_;
      pmid_dry = pmid_dry_;
    } // set_variables
  }; // p3_preamble
  /* --------------------------------------------------------------------------------------------*/
  // Most individual processes have a post-processing step that derives variables needed by the rest
  // of the model, using outputs from this process.
  // Structure to handle the generation of data needed by the rest of the model based on output from
  // p3_main.
  struct p3_postamble {
    p3_postamble() = default;
    // Functor for Kokkos loop to pre-process every run step
    KOKKOS_INLINE_FUNCTION
    void operator()(const int icol) const {
      for (int ipack=0;ipack<m_npack;ipack++) {
        // Update the atmospheric temperature and the previous temperature.
        T_atm(icol,ipack)  = PF::calculate_T_from_theta(th_atm(icol,ipack),pmid(icol,ipack));
        T_prev(icol,ipack) = T_atm(icol,ipack);

        /*----------------------------------------------------------------------------------------------------------------------
         *DRY-TO-WET MMRs:
         *-----------------
         *Since the host model (or AD) needs wet mixing ratios, we need to convert dry mixing ratios from P3 to
         *wet mixing ratios.

         *NOTE: Function calculate_wetmmr_from_drymmr takes 2 arguments: ( dry mmr and "dry" water vapor mixing ratio)

         *IMPORTANT:Convert "qv dry mmr" to "qv wet mmr" after converting all other constituents to wet mmr as "qv" (as dry mmr)
         * is an input for converting all other constituents to have wet mmr.
         *----------------------------------------------------------------------------------------------------------------------
         */
        //Units of all constituents below are [kg/kg(wet-air)] for mass and [#/kg(wet-air)] for number
        qc(icol,ipack) = PF::calculate_wetmmr_from_drymmr(qc(icol,ipack), qv(icol,ipack));//Cloud liquid mass
        nc(icol,ipack) = PF::calculate_wetmmr_from_drymmr(nc(icol,ipack), qv(icol,ipack));//Cloud liquid number
        qr(icol,ipack) = PF::calculate_wetmmr_from_drymmr(qr(icol,ipack), qv(icol,ipack));//Rain mass
        nr(icol,ipack) = PF::calculate_wetmmr_from_drymmr(nr(icol,ipack), qv(icol,ipack));//Rain number
        qi(icol,ipack) = PF::calculate_wetmmr_from_drymmr(qi(icol,ipack), qv(icol,ipack));//Cloud ice mass
        ni(icol,ipack) = PF::calculate_wetmmr_from_drymmr(ni(icol,ipack), qv(icol,ipack));//Cloud ice number
        qm(icol,ipack) = PF::calculate_wetmmr_from_drymmr(qm(icol,ipack), qv(icol,ipack));//Rimmed ice mass
        bm(icol,ipack) = PF::calculate_wetmmr_from_drymmr(bm(icol,ipack), qv(icol,ipack));//Rimmed ice number

        // ^^ Ensure that qv is "dry mmr" till this point ^^
        //NOTE:Convert "qv" to wet mmr in the end after converting all other constituents to wet mmr
        qv(icol,ipack) = PF::calculate_wetmmr_from_drymmr(qv(icol,ipack), qv(icol,ipack));

        // Update qv_prev with qv(which should now be a wet mmr) so that qv_prev is in wet mmr
        qv_prev(icol,ipack) = qv(icol,ipack);

        // Rescale effective radius' into microns
        diag_eff_radius_qc(icol,ipack) *= 1e6;
        diag_eff_radius_qi(icol,ipack) *= 1e6;
      } // for ipack
    } // operator
    // Local variables
    int m_ncol, m_npack;
    view_2d       T_atm;
    view_2d_const pmid;
    view_2d       th_atm;
    view_2d       T_prev;
    view_2d       qv;
    view_2d       qc;
    view_2d       nc;
    view_2d       qr;
    view_2d       nr;
    view_2d       qi;
    view_2d       qm;
    view_2d       ni;
    view_2d       bm;
    view_2d       qv_prev;
    view_2d       diag_eff_radius_qc;
    view_2d       diag_eff_radius_qi;
    // Assigning local values
    void set_variables(const int ncol, const int npack,
                    const view_2d& th_atm_, const view_2d_const& pmid_, const view_2d& T_atm_, const view_2d& T_prev_,
                    const view_2d& qv_, const view_2d& qc_, const view_2d& nc_, const view_2d& qr_, const view_2d& nr_,
                    const view_2d& qi_, const view_2d& qm_, const view_2d& ni_, const view_2d& bm_,
                    const view_2d& qv_prev_, const view_2d& diag_eff_radius_qc_,
                    const view_2d& diag_eff_radius_qi_)
    {
      m_ncol  = ncol;
      m_npack = npack;
      // IN
      th_atm      = th_atm_;
      pmid        = pmid_;
      qv          = qv_;
      qc          = qc_;
      nc          = nc_;
      qr          = qr_;
      nr          = nr_;
      qi          = qi_;
      qm          = qm_;
      ni          = ni_;
      bm          = bm_;
      // OUT
      T_atm              = T_atm_;
      T_prev             = T_prev_;
      qv_prev            = qv_prev_;
      diag_eff_radius_qc = diag_eff_radius_qc_;
      diag_eff_radius_qi = diag_eff_radius_qi_;
      // TODO: This is a list of variables not yet defined for post-processing, but are
      // defined in the F90 p3 interface code.  So this list will need to be checked as
      // new processes come online to make sure their requirements from p3 are being met.
      // qme, vap_liq_exchange
      // ENERGY Conservation: prec_str, snow_str
      // RAD Vars: icinc, icwnc, icimrst, icwmrst
      // COSP Vars: flxprc, flxsnw, flxprc, flxsnw, cvreffliq, cvreffice, reffrain, reffsnow
    } // set_variables
  }; // p3_postamble
  /* --------------------------------------------------------------------------------------------*/

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    // 1d view scalar, size (ncol)
    static constexpr int num_1d_scalar = 0; //no 2d vars now, but keeping 1d struct for future expansion
    // 2d view packed, size (ncol, nlev_packs)
    static constexpr int num_2d_vector = 11;
    static constexpr int num_2dp1_vector = 2;

    uview_2d inv_exner;
    uview_2d th_atm;
    uview_2d cld_frac_l;
    uview_2d cld_frac_i;
    uview_2d cld_frac_r;
    uview_2d dz;
    uview_2d pdel_dry;
    uview_2d pmid_dry;
    uview_2d qv2qi_depos_tend;
    uview_2d rho_qi;
    uview_2d precip_liq_flux; //nlev+1
    uview_2d precip_ice_flux; //nlev+1
    uview_2d unused;

    suview_2d col_location;

    Spack* wsm_data;
  };

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const int dt);
  void finalize_impl   ();

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  // Keep track of field dimensions and the iteration count
  Int m_num_cols;
  Int m_num_levs;
  Int m_nk_pack;

  // Struct which contains local variables
  Buffer m_buffer;

  // Store the structures for each arguement to p3_main;
  P3F::P3PrognosticState   prog_state;
  P3F::P3DiagnosticInputs  diag_inputs;
  P3F::P3DiagnosticOutputs diag_outputs;
  P3F::P3HistoryOnly       history_only;
  P3F::P3LookupTables      lookup_tables;
  P3F::P3Infrastructure    infrastructure;
  p3_preamble              p3_preproc;
  p3_postamble             p3_postproc;

  // WSM for internal local variables
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr;

  std::shared_ptr<const AbstractGrid>   m_grid;
  // Iteration count is internal to P3 and keeps track of the number of times p3_main has been called.
  // infrastructure.it is passed as an arguement to p3_main and is used for identifying which iteration an error occurs.

}; // class P3Microphysics

} // namespace scream

#endif // SCREAM_P3_MICROPHYSICS_HPP
