#ifndef SCREAM_SHOC_FUNCTIONS_F90_HPP
#define SCREAM_SHOC_FUNCTIONS_F90_HPP

#include "share/scream_types.hpp"
#include "physics/share/physics_test_data.hpp"

#include "shoc_functions.hpp"

#include <vector>
#include <array>
#include <utility>

//
// Bridge functions to call fortran version of shoc functions from C++
//

namespace scream {
namespace shoc {

//Create data structure to hold data for shoc_grid
struct SHOCGridData : public PhysicsTestData {
  // Inputs
  Real *zt_grid, *zi_grid, *pdel;

  // In/out
  Real *dz_zt, *dz_zi, *rho_zt;

  SHOCGridData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&zt_grid, &dz_zt, &pdel, &rho_zt}, {&zi_grid, &dz_zi}) {}

  DATA_COPY_CONS(SHOCGridData, 3);
  ASSIGN(SHOCGridData, 0);
};

//Create data structure to hold data for integ_column_stability
struct SHOCColstabData : public PhysicsTestData {
  // Inputs
  Real *dz_zt, *pres, *brunt;

  // Output
  Real *brunt_int;

  SHOCColstabData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&dz_zt, &pres, &brunt}, {&brunt_int}) {}

  DATA_COPY_CONS(SHOCColstabData, 2);
  ASSIGN(SHOCColstabData, 0);
};//SHOCColstabData

//Create data structure to hold data for compute_shr_prod
struct SHOCTkeshearData : public PhysicsTestData {
  // Inputs
  Real *dz_zi, *u_wind, *v_wind;

  // In/out
  Real *sterm;

  //functions to initialize data
  SHOCTkeshearData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&u_wind, &v_wind}, {&dz_zi, &sterm}) {}

  DATA_COPY_CONS(SHOCTkeshearData, 3);
  ASSIGN(SHOCTkeshearData, 0);
};//SHOCTkeshearData

//Create data structure to hold data for isotropic_ts
struct SHOCIsotropicData : public PhysicsTestData {
  // Inputs
  Real *tke, *a_diss, *brunt, *brunt_int;

  // Output
  Real *isotropy;

  //functions to initialize data
  SHOCIsotropicData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&tke, &a_diss, &brunt, &isotropy}, {&brunt_int}) {}

  DATA_COPY_CONS(SHOCIsotropicData, 2);
  ASSIGN(SHOCIsotropicData, 0);
};//SHOCIsotropicData

//Create data structure to hold data for adv_sgs_tke
struct SHOCAdvsgstkeData : public PhysicsTestData {
  // Inputs
  Real dtime;
  Real *shoc_mix, *wthv_sec, *sterm_zt, *tk;

  // In/out
  Real *tke;

  // Outputs
  Real *a_diss;

  //functions to initialize data
  SHOCAdvsgstkeData(Int shcol_, Int nlev_, Real dtime_) :
    PhysicsTestData(shcol_, nlev_, {&shoc_mix, &wthv_sec, &sterm_zt, &tk, &tke, &a_diss}), dtime(dtime_) {}

  DATA_COPY_CONS(SHOCAdvsgstkeData, 3);
  ASSIGN(SHOCAdvsgstkeData, 1, dtime);
};//SHOCAdvsgstkeData

//Create data structure to hold data for eddy_diffusivities
struct SHOCEddydiffData : public PhysicsTestData {
  // Inputs
  Real *pblh, *obklen, *zt_grid, *shoc_mix, *sterm_zt,
        *isotropy, *tke;

  // Output
  Real *tk, *tkh;

  //functions to initialize data
  SHOCEddydiffData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&zt_grid, &shoc_mix, &isotropy, &tke, &tk, &tkh, &sterm_zt}, {&obklen, &pblh}) {}

  DATA_COPY_CONS(SHOCEddydiffData, 2);
  ASSIGN(SHOCEddydiffData, 0);
};//SHOCEddydiffData


//create data structure for update_host_dse
struct SHOCEnergydseData : public PhysicsTestData {
  // Inputs
  Real *thlm, *shoc_ql, *exner, *zt_grid, *phis;

  // Output
  Real *host_dse;

  //functions to initialize data
  SHOCEnergydseData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&thlm, &shoc_ql, &exner, &zt_grid, &host_dse}, {&phis}) {}

  DATA_COPY_CONS(SHOCEnergydseData, 2);
  ASSIGN(SHOCEnergydseData, 0);
};//SHOCEnergydseData

//create data structure for shoc_energy_integrals
struct SHOCEnergyintData : public PhysicsTestData {
  // Inputs
  Real *host_dse, *pdel, *rtm, *rcm, *u_wind, *v_wind;

  // Output
  Real *se_int, *ke_int, *wv_int, *wl_int;

  //functions to initialize data
  SHOCEnergyintData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&host_dse, &pdel, &rtm, &rcm, &u_wind, &v_wind}, {&se_int, &ke_int, &wv_int, &wl_int}) {}

  DATA_COPY_CONS(SHOCEnergyintData, 2);
  ASSIGN(SHOCEnergyintData, 0);
};//SHOCEnergyintData

//Create data structure for shoc_energy_total_fixer
struct SHOCEnergytotData : public PhysicsTestData {
  // Inputs
  Int nadv;
  Real dtime;
  Real *zt_grid, *zi_grid, *se_b, *ke_b, *wv_b, *wl_b, *se_a;
  Real *ke_a, *wv_a, *wl_a, *wthl_sfc, *wqw_sfc, *rho_zt;

  // Output
  Real *te_a, *te_b;

  //functions to initialize data for shoc_energy_total_fixer
  SHOCEnergytotData(Int shcol_, Int nlev_, Int nlevi_, Real dtime_, Int nadv_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&zt_grid, &rho_zt}, {&zi_grid}, {&se_b, &ke_b, &wv_b, &wl_b, &se_a, &ke_a, &wv_a, &wl_a, &wthl_sfc, &wqw_sfc, &te_a, &te_b}), nadv(nadv_), dtime(dtime_) {}

  DATA_COPY_CONS(SHOCEnergytotData, 5);
  ASSIGN(SHOCEnergytotData, 2, dtime, nadv);
};//SHOCEnergytotData

//create data structure for shoc_energy_threshold_fixer
struct SHOCEnergythreshfixerData : public PhysicsTestData {
  // Inputs
  Real *pint, *tke, *te_a, *te_b;

  // In/out
  Real *se_dis;
  Int *shoctop;

  //functions to initialize data
  SHOCEnergythreshfixerData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&tke}, {&pint}, {&se_dis, &te_a, &te_b}, {&shoctop}) {}

  DATA_COPY_CONS(SHOCEnergythreshfixerData, 3);
  ASSIGN(SHOCEnergythreshfixerData, 0);
};//SHOCEnergythreshfixerData

//create data structure for shoc_energy_dse_fixer
struct SHOCEnergydsefixerData : public PhysicsTestData {
  // Inputs
  Real *se_dis;
  Int *shoctop;

  // In/out
  Real *host_dse;

  //functions to initialize data
  SHOCEnergydsefixerData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&host_dse}, {&se_dis}, {&shoctop}) {}

  DATA_COPY_CONS(SHOCEnergydsefixerData, 2);
  ASSIGN(SHOCEnergydsefixerData, 0);
};//SHOCEnergydsefixerData

//Create data structure to hold data for calc_shoc_vertflux
struct SHOCVertfluxData : public PhysicsTestData {
  // Inputs
  Real *tkh_zi, *dz_zi, *invar;

  // In/out
  Real *vertflux;

  SHOCVertfluxData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&invar}, {&tkh_zi, &dz_zi, &vertflux}) {}

  DATA_COPY_CONS(SHOCVertfluxData, 3);
  ASSIGN(SHOCVertfluxData, 0);
}; //SHOCVertfluxData

//Create data structure to hold data for calc_shoc_varorcovar
struct SHOCVarorcovarData : public PhysicsTestData {
  // Inputs
  Real tunefac;
  Real *tkh_zi, *dz_zi, *isotropy_zi, *invar1, *invar2;

  // In/out
  Real *varorcovar;

  SHOCVarorcovarData(Int shcol_, Int nlev_, Int nlevi_, Real tunefac_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&invar1, &invar2}, {&tkh_zi, &dz_zi, &isotropy_zi, &varorcovar}), tunefac(tunefac_) {}

  DATA_COPY_CONS(SHOCVarorcovarData, 4);
  ASSIGN(SHOCVarorcovarData, 1, tunefac);
};//SHOCVarorcovarData

//Create data structure to hold data for compute_brunt_shoc_length
struct SHOCBruntlengthData : public PhysicsTestData {
  // Inputs
  Real *dz_zt, *thv, *thv_zi;

  // In/out
  Real *brunt;

  SHOCBruntlengthData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&dz_zt, &thv, &brunt}, {&thv_zi}) {}

  DATA_COPY_CONS(SHOCBruntlengthData, 3);
  ASSIGN(SHOCBruntlengthData, 0)
};//SHOCBruntlengthData

//Create data structure to hold data for compute_l_inf_shoc_length
struct SHOCInflengthData : public PhysicsTestData {
  // Inputs
  Real *zt_grid, *dz_zt, *tke;

  // In/out
  Real *l_inf;

  SHOCInflengthData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&zt_grid, &dz_zt, &tke}, {&l_inf}) {}

  DATA_COPY_CONS(SHOCInflengthData, 2);
  ASSIGN(SHOCInflengthData, 0);
};//SHOCInflengthData

//Create data structure to hold data for compute_vel_shoc_length
struct SHOCConvvelData : public PhysicsTestData {
  // Inputs
  Real *pblh, *zt_grid, *dz_zt, *thv, *wthv_sec;

  // In/out
  Real *conv_vel;

  SHOCConvvelData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&zt_grid, &dz_zt, &thv, &wthv_sec}, {&conv_vel, &pblh}) {}

  DATA_COPY_CONS(SHOCConvvelData, 2);
  ASSIGN(SHOCConvvelData, 0);
};//SHOCConvvelData

//Create data structure to hold data for compute_conv_time_shoc_length
struct SHOCConvtimeData : public PhysicsTestData {
  // Inputs
  Real *pblh, *conv_vel;

  // In/out
  Real *tscale;

  SHOCConvtimeData(Int shcol_) :
    PhysicsTestData(shcol_, {&conv_vel, &pblh, &tscale}) {}

  DATA_COPY_CONS(SHOCConvtimeData, 1);
  ASSIGN(SHOCConvtimeData, 0);
};//SHOCConvtimeData

//Create data structure to hold data for compute_shoc_mix_shoc_length
struct SHOCMixlengthData : public PhysicsTestData {
  // Inputs
  Real *tke, *brunt, *tscale, *zt_grid, *l_inf;

  // In/out
  Real *shoc_mix;

  SHOCMixlengthData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&tke, &brunt, &zt_grid, &shoc_mix}, {&l_inf, &tscale}) {}

  DATA_COPY_CONS(SHOCMixlengthData, 2);
  ASSIGN(SHOCMixlengthData, 0);
};//SHOCMixlengthData

//Create data structure to hold data for check_length_scale_shoc_length
struct SHOCMixcheckData : public PhysicsTestData {
  // Inputs
  Real *host_dx, *host_dy;

  // In/out
  Real *shoc_mix;

  SHOCMixcheckData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&shoc_mix}, {&host_dx, &host_dy}) {}

  DATA_COPY_CONS(SHOCMixcheckData, 2);
  ASSIGN(SHOCMixcheckData, 0);
};//SHOCMixcheckData

//Create data structure to hold data for clipping_diag_third_shoc_moments
struct SHOCClipthirdmomsData : public SHOCDataBase {
  // Inputs
  Real *w_sec_zi;

  // In/out
  Real *w3;

  SHOCClipthirdmomsData(Int shcol_, Int nlevi_) :
    SHOCDataBase(shcol_, 0, nlevi_, {}, {&w_sec_zi, &w3}, {}) {}
  SHOCClipthirdmomsData(const SHOCClipthirdmomsData &rhs) :
    SHOCDataBase(rhs, {}, {&w_sec_zi, &w3}, {}) {}
  SHOCClipthirdmomsData &operator=(const SHOCClipthirdmomsData &rhs)
  { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCClipthirdmomsData

struct SHOCSecondMomentSrfData : public PhysicsTestData {
  // Inputs
  Real *wthl, *uw, *vw;

  // out
  Real *ustar2, *wstar;

  SHOCSecondMomentSrfData(Int shcol_) :
    PhysicsTestData(shcol_, {&wthl, &uw, &vw, &ustar2, &wstar}) {}

  DATA_COPY_CONS(SHOCSecondMomentSrfData, 1);
  ASSIGN(SHOCSecondMomentSrfData, 0);
};

//Create data structure to hold data for linear_interp
struct SHOCLinearintData : public PhysicsTestData {
  // Inputs
  Real minthresh;
  Real *x1, *x2, *y1;

  // In/out
  Real *y2;

  SHOCLinearintData(Int shcol_, Int nlev_, Int nlevi_, Real minthresh_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&x1, &y1}, {&x2, &y2}), minthresh(minthresh_) {}

  DATA_COPY_CONS(SHOCLinearintData, 4);
  ASSIGN(SHOCLinearintData, 1, minthresh);
};//SHOCLinearintData

struct SHOCPblintdInitPotData : public SHOCDataBase {
  // inputs
  Real *thl, *ql, *q;

  // outputs
  Real *thv;

  SHOCPblintdInitPotData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&thl, &ql, &q, &thv}, {}, {}) {}
  SHOCPblintdInitPotData(const SHOCPblintdInitPotData &rhs) : SHOCDataBase(rhs, {&thl, &ql, &q, &thv}, {}, {}) {}
  SHOCPblintdInitPotData &operator=(const SHOCPblintdInitPotData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
}; 

//Create data structure to hold data for shoc_assumed_pdf_tilda_to_real
struct SHOCPDFtildaData
{
  // inputs
  Real w_first, sqrtw2;
  
  // outputs
  Real w1;
};

// Create data structure to hold data for shoc_assumed_pdf_vv_parameters
struct SHOCPDFvvparamData
{
  // inputs
  Real w_first, w_sec, w3var;
  
  // outputs
  Real Skew_w, w1_1, w1_2, w2_1, w2_2, a;
};

// Create data structure to hold data for shoc_assumed_pdf_thl_parameters
struct SHOCPDFthlparamData
{
  // inputs
  Real wthlsec, sqrtw2, sqrtthl, thlsec, thl_first, w1_1, w1_2, Skew_w, a;
  bool dothetal_skew;
  
  // outputs
  Real thl1_1, thl1_2, thl2_1, thl2_2, sqrtthl2_1, sqrtthl2_2;
};

// Create data structure to hold data for shoc_assumed_pdf_qw_parameters
struct SHOCPDFqwparamData
{
  // inputs
  Real wqwsec, qwsec, sqrtw2, sqrtqt, qw_first, w1_1, w1_2, Skew_w, a;
  
  // outputs
  Real qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2;
};

// Create data structure to hold data for shoc_assumed_pdf_inplume_correlations
struct SHOCPDFinplumeData
{
  // inputs
  Real sqrtqw2_1,sqrtthl2_1,a,sqrtqw2_2,sqrtthl2_2;
  Real qwthlsec,qw1_1,qw_first,thl1_1,thl_first,qw1_2,thl1_2;
  
  // outputs
  Real r_qwthl_1;
};

//Create data structure to hold data for shoc_assumed_pdf_compute_temperature
struct SHOCPDFcomptempData
{
  // inputs
  Real thl1, basepres, pval;
  
  // outputs
  Real Tl1;
};

//Create data structure to hold data for shoc_assumed_pdf_compute_qs
struct SHOCPDFcompqsData
{
  // inputs
  Real Tl1_1, Tl1_2, pval;
  
  // outputs
  Real qs1, beta1, qs2, beta2;
};

//Create data structure to hold data for shoc_assumed_pdf_compute_s
struct SHOCPDFcompsData
{
  // inputs
  Real qw1, qs1, beta, pval, thl2, qw2, sqrtthl2, sqrtqw2, r_qwthl;
  
  // outputs
  Real s, std_s, qn, C;
};

//Create data structure to hold data for shoc_assumed_pdf_compute_sgs_liquid
struct SHOCPDFcompsgsliqData
{
  // inputs
  Real a, ql1, ql2;
  
  // outputs
  Real shoc_ql;
};

//Create data structure to hold data for shoc_assumed_pdf_compute_cloud_liquid_variance
struct SHOCPDFcompcloudvarData
{
  // inputs
  Real a, s1, ql1, C1, std_s1, s2, ql2, C2, std_s2, shoc_ql;
  
  // outputs
  Real shoc_ql2;
};

//Create data structure to hold data for shoc_assumed_pdf_compute_liquid_water_flux
struct SHOCPDFcompliqfluxData
{
  // inputs
  Real a, w1_1, w_first, ql1, w1_2, ql2;
  
  // outputs
  Real wqls;
};

//Create data structure to hold data for shoc_assumed_pdf_compute_buoyancy_flux
struct SHOCPDFcompbuoyfluxData
{
  // inputs
  Real wthlsec, epsterm, wqwsec, pval, wqls;
  
  // outputs
  Real wthv_sec;
};

struct SHOCSecondMomentUbycondData : public SHOCDataBase {
  // Outputs
  Real *thl, *qw, *wthl, *wqw, *qwthl, *uw, *vw, *wtke;

  SHOCSecondMomentUbycondData(Int shcol_) :
     SHOCDataBase(shcol_, 1, 0, {&thl, &qw, &wthl, &wqw, &qwthl, &uw, &vw, &wtke}, {}, {}) {}
  SHOCSecondMomentUbycondData(const SHOCSecondMomentUbycondData &rhs) :
     SHOCDataBase(rhs, {&thl, &qw, &wthl, &wqw, &qwthl, &uw, &vw, &wtke}, {}, {}) {}
  SHOCSecondMomentUbycondData &operator=(const SHOCSecondMomentUbycondData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};

//
// Glue functions to call fortran from from C++ with the Data struct
//

// This function initialzes the grid used by shoc. Given the
// locations of the cell center (location of thermodynaics quantities), cell
// interfaces, and pressure gradient the functon returns dz_zi, dz_zt,
// and density.
void shoc_grid                      (SHOCGridData &d);
void update_host_dse                (SHOCEnergydseData &d);
void shoc_energy_integrals          (SHOCEnergyintData &d);
void shoc_energy_total_fixer        (SHOCEnergytotData &d);
void shoc_energy_threshold_fixer    (SHOCEnergythreshfixerData &d);
void shoc_energy_dse_fixer          (SHOCEnergydsefixerData &d);
void calc_shoc_vertflux             (SHOCVertfluxData &d);
void calc_shoc_varorcovar           (SHOCVarorcovarData &d);
void integ_column_stability         (SHOCColstabData &d);
void compute_shr_prod               (SHOCTkeshearData &d);
void isotropic_ts                   (SHOCIsotropicData &d);
void adv_sgs_tke                    (SHOCAdvsgstkeData &d);
void eddy_diffusivities             (SHOCEddydiffData &d);
void compute_brunt_shoc_length      (SHOCBruntlengthData &d);
void compute_l_inf_shoc_length      (SHOCInflengthData &d);
void compute_conv_vel_shoc_length   (SHOCConvvelData &d);
void compute_conv_time_shoc_length  (SHOCConvtimeData &d);
void compute_shoc_mix_shoc_length   (SHOCMixlengthData &d);
void check_length_scale_shoc_length (SHOCMixcheckData &d);
void shoc_diag_second_moments_srf   (SHOCSecondMomentSrfData& d);
void linear_interp                  (SHOCLinearintData &d);

//
// _f functions decls
//
extern "C" {

void calc_shoc_varorcovar_f(Int shcol, Int nlev, Int nlevi, Real tunefac,
                            Real *isotropy_zi, Real *tkh_zi, Real *dz_zi,
                            Real *invar1, Real *invar2, Real *varorcovar);
void calc_shoc_vertflux_f(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux);
void shoc_diag_second_moments_srf_f(Int shcol, Real* wthl, Real* uw, Real* vw,
                         Real* ustar2, Real* wstar);
void shoc_pblintd_init_pot_f(Int shcol, Int nlev, Real* thl, Real* ql, Real* q, Real* thv);

void shoc_diag_second_moments_ubycond_f(Int shcol, Real* thl, Real* qw, Real* wthl,
                          Real* wqw, Real* qwthl, Real* uw, Real* vw, Real* wtke);

}

}  // namespace shoc
}  // namespace scream

#endif // SCREAM_SHOC_FUNCTIONS_F90_HPP
