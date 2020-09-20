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

  SHOCGridData(const SHOCGridData &rhs) : PhysicsTestData(rhs, {&zt_grid, &dz_zt, &pdel, &rho_zt, &zi_grid, &dz_zi}) {}

  SHOCGridData &operator=(const SHOCGridData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
};

//Create data structure to hold data for integ_column_stability
struct SHOCColstabData : public PhysicsTestData {
  // Inputs
  Real *dz_zt, *pres, *brunt;

  // Output
  Real *brunt_int;

  SHOCColstabData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&dz_zt, &pres, &brunt}, {&brunt_int}) {}
  SHOCColstabData(const SHOCColstabData &rhs) : PhysicsTestData(rhs, {&dz_zt, &pres, &brunt, &brunt_int}) {}
  SHOCColstabData &operator=(const SHOCColstabData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
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
  SHOCTkeshearData(const SHOCTkeshearData &rhs) : PhysicsTestData(rhs, {&u_wind, &v_wind, &dz_zi, &sterm}) {}
  SHOCTkeshearData &operator=(const SHOCTkeshearData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
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
  SHOCIsotropicData(const SHOCIsotropicData &rhs) : PhysicsTestData(rhs, {&tke, &a_diss, &brunt, &isotropy, &brunt_int}) {}
  SHOCIsotropicData &operator=(const SHOCIsotropicData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
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
  SHOCAdvsgstkeData(const SHOCAdvsgstkeData &rhs) : PhysicsTestData(rhs, {&shoc_mix, &wthv_sec, &sterm_zt, &tk, &tke, &a_diss}), dtime(rhs.dtime) {}
  SHOCAdvsgstkeData &operator=(const SHOCAdvsgstkeData &rhs)
  { PhysicsTestData::operator=(rhs); dtime = rhs.dtime; return *this; }
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
  SHOCEddydiffData(const SHOCEddydiffData &rhs) : PhysicsTestData(rhs, {&zt_grid, &shoc_mix, &isotropy, &tke, &tk, &tkh, &sterm_zt, &obklen, &pblh}) {}
  SHOCEddydiffData &operator=(const SHOCEddydiffData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
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
  SHOCEnergydseData(const SHOCEnergydseData &rhs) : PhysicsTestData(rhs, {&thlm, &shoc_ql, &exner, &zt_grid, &host_dse, &phis}) {}
  SHOCEnergydseData &operator=(const SHOCEnergydseData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
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
  SHOCEnergyintData(const SHOCEnergyintData &rhs) : PhysicsTestData(rhs, {&host_dse, &pdel, &rtm, &rcm, &u_wind, &v_wind, &se_int, &ke_int, &wv_int, &wl_int}) {}
  SHOCEnergyintData &operator=(const SHOCEnergyintData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
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
  SHOCEnergytotData(const SHOCEnergytotData &rhs) : PhysicsTestData(rhs, {&zt_grid, &rho_zt, &zi_grid, &se_b, &ke_b, &wv_b, &wl_b, &se_a, &ke_a, &wv_a, &wl_a, &wthl_sfc, &wqw_sfc, &te_a, &te_b}) {}
  SHOCEnergytotData &operator=(const SHOCEnergytotData &rhs) { PhysicsTestData::operator=(rhs); dtime = rhs.dtime; nadv = rhs.nadv; return *this; }
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
  SHOCEnergythreshfixerData(const SHOCEnergythreshfixerData &rhs) : PhysicsTestData(rhs, {&tke, &pint, &se_dis, &te_a, &te_b}, {&shoctop}) {}
  SHOCEnergythreshfixerData &operator=(const SHOCEnergythreshfixerData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
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
  SHOCEnergydsefixerData(const SHOCEnergydsefixerData &rhs) : PhysicsTestData(rhs, {&host_dse, &se_dis}, {&shoctop}) {}
  SHOCEnergydsefixerData &operator=(const SHOCEnergydsefixerData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
};//SHOCEnergydsefixerData

//Create data structure to hold data for calc_shoc_vertflux
struct SHOCVertfluxData : public PhysicsTestData {
  // Inputs
  Real *tkh_zi, *dz_zi, *invar;

  // In/out
  Real *vertflux;

  SHOCVertfluxData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&invar}, {&tkh_zi, &dz_zi, &vertflux}) {}
  SHOCVertfluxData(const SHOCVertfluxData &rhs) : PhysicsTestData(rhs, {&invar, &tkh_zi, &dz_zi, &vertflux}) {}
  SHOCVertfluxData &operator=(const SHOCVertfluxData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
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
  SHOCVarorcovarData(const SHOCVarorcovarData &rhs) :
    PhysicsTestData(rhs, {&invar1, &invar2, &tkh_zi, &dz_zi, &isotropy_zi, &varorcovar}), tunefac(rhs.tunefac) {}
  SHOCVarorcovarData &operator=(const SHOCVarorcovarData &rhs)
  { PhysicsTestData::operator=(rhs); tunefac = rhs.tunefac; return *this; }
};//SHOCVarorcovarData

//Create data structure to hold data for compute_brunt_shoc_length
struct SHOCBruntlengthData : public PhysicsTestData {
  // Inputs
  Real *dz_zt, *thv, *thv_zi;

  // In/out
  Real *brunt;

  SHOCBruntlengthData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData(shcol_, nlev_, nlevi_, {&dz_zt, &thv, &brunt}, {&thv_zi}) {}
  SHOCBruntlengthData(const SHOCBruntlengthData &rhs) :
    PhysicsTestData(rhs, {&dz_zt, &thv, &brunt, &thv_zi}) {}
  SHOCBruntlengthData &operator=(const SHOCBruntlengthData &rhs)
  { PhysicsTestData::operator=(rhs); return *this; }
};//SHOCBruntlengthData

//Create data structure to hold data for compute_l_inf_shoc_length
struct SHOCInflengthData : public PhysicsTestData {
  // Inputs
  Real *zt_grid, *dz_zt, *tke;

  // In/out
  Real *l_inf;

  SHOCInflengthData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&zt_grid, &dz_zt, &tke}, {&l_inf}) {}
  SHOCInflengthData(const SHOCInflengthData &rhs) :
    PhysicsTestData(rhs, {&zt_grid, &dz_zt, &tke, &l_inf}) {}
  SHOCInflengthData &operator=(const SHOCInflengthData &rhs)
  { PhysicsTestData::operator=(rhs); return *this; }
};//SHOCInflengthData

//Create data structure to hold data for compute_vel_shoc_length
struct SHOCConvvelData : public PhysicsTestData {
  // Inputs
  Real *pblh, *zt_grid, *dz_zt, *thv, *wthv_sec;

  // In/out
  Real *conv_vel;

  SHOCConvvelData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&zt_grid, &dz_zt, &thv, &wthv_sec}, {&conv_vel, &pblh}) {}
  SHOCConvvelData(const SHOCConvvelData &rhs) :
    PhysicsTestData(rhs, {&zt_grid, &dz_zt, &thv, &wthv_sec, &conv_vel, &pblh}) {}
  SHOCConvvelData &operator=(const SHOCConvvelData &rhs)
  { PhysicsTestData::operator=(rhs); return *this; }
};//SHOCConvvelData

//Create data structure to hold data for compute_conv_time_shoc_length
struct SHOCConvtimeData : public PhysicsTestData {
  // Inputs
  Real *pblh, *conv_vel;

  // In/out
  Real *tscale;

  SHOCConvtimeData(Int shcol_) :
    PhysicsTestData(shcol_, {&conv_vel, &pblh, &tscale}) {}
  SHOCConvtimeData(const SHOCConvtimeData &rhs) :
    PhysicsTestData(rhs, {&conv_vel, &pblh, &tscale}) {}
  SHOCConvtimeData &operator=(const SHOCConvtimeData &rhs)
  { PhysicsTestData::operator=(rhs); return *this; }
};//SHOCConvtimeData

//Create data structure to hold data for compute_shoc_mix_shoc_length
struct SHOCMixlengthData : public PhysicsTestData {
  // Inputs
  Real *tke, *brunt, *tscale, *zt_grid, *l_inf;

  // In/out
  Real *shoc_mix;

  SHOCMixlengthData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&tke, &brunt, &zt_grid, &shoc_mix}, {&l_inf, &tscale}) {}
  SHOCMixlengthData(const SHOCMixlengthData &rhs) :
    PhysicsTestData(rhs, {&tke, &brunt, &zt_grid, &shoc_mix, &l_inf, &tscale}) {}
  SHOCMixlengthData &operator=(const SHOCMixlengthData &rhs)
  { PhysicsTestData::operator=(rhs); return *this; }
};//SHOCMixlengthData

//Create data structure to hold data for check_length_scale_shoc_length
struct SHOCMixcheckData : public PhysicsTestData {
  // Inputs
  Real *host_dx, *host_dy;

  // In/out
  Real *shoc_mix;

  SHOCMixcheckData(Int shcol_, Int nlev_) :
    PhysicsTestData(shcol_, nlev_, {&shoc_mix}, {&host_dx, &host_dy}) {}
  SHOCMixcheckData(const SHOCMixcheckData &rhs) :
    PhysicsTestData(rhs, {&shoc_mix, &host_dx, &host_dy}) {}
  SHOCMixcheckData &operator=(const SHOCMixcheckData &rhs)
  { PhysicsTestData::operator=(rhs); return *this; }
};//SHOCMixcheckData

struct SHOCSecondMomentSrfData : public PhysicsTestData {
  // Inputs
  Real *wthl, *uw, *vw;

  // out
  Real *ustar2, *wstar;

  SHOCSecondMomentSrfData(Int shcol_) :
  PhysicsTestData(shcol_, {&wthl, &uw, &vw, &ustar2, &wstar}) {}
  SHOCSecondMomentSrfData(const SHOCSecondMomentSrfData &rhs) : PhysicsTestData(rhs, {&wthl, &uw, &vw, &ustar2, &wstar}) {}
  SHOCSecondMomentSrfData &operator=(const SHOCSecondMomentSrfData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
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
  SHOCLinearintData(const SHOCLinearintData &rhs) :
    PhysicsTestData(rhs, {&x1, &y1, &x2, &y2}), minthresh(rhs.minthresh) {}
  SHOCLinearintData &operator=(const SHOCLinearintData &rhs)
  { PhysicsTestData::operator=(rhs); minthresh = rhs.minthresh; return *this; }
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
void shoc_grid(SHOCGridData &d);
void update_host_dse(SHOCEnergydseData &d);
void shoc_energy_integrals(SHOCEnergyintData &d);
void shoc_energy_total_fixer(SHOCEnergytotData &d);
void shoc_energy_threshold_fixer(SHOCEnergythreshfixerData &d);
void shoc_energy_dse_fixer(SHOCEnergydsefixerData &d);
void calc_shoc_vertflux(SHOCVertfluxData &d);
void calc_shoc_varorcovar(SHOCVarorcovarData &d);
void integ_column_stability(SHOCColstabData &d);
void compute_shr_prod(SHOCTkeshearData &d);
void isotropic_ts(SHOCIsotropicData &d);
void adv_sgs_tke(SHOCAdvsgstkeData &d);
void eddy_diffusivities(SHOCEddydiffData &d);
void compute_brunt_shoc_length(SHOCBruntlengthData &d);
void compute_l_inf_shoc_length(SHOCInflengthData &d);
void compute_conv_vel_shoc_length(SHOCConvvelData &d);
void compute_conv_time_shoc_length(SHOCConvtimeData &d);
void compute_shoc_mix_shoc_length(SHOCMixlengthData &d);
void check_length_scale_shoc_length(SHOCMixcheckData &d);
void shoc_diag_second_moments_srf(SHOCSecondMomentSrfData& d);
void linear_interp(SHOCLinearintData &d);
void shoc_pblintd_init_pot(SHOCPblintdInitPotData &d);
void shoc_assumed_pdf_tilda_to_real(SHOCPDFtildaData &d);
void shoc_assumed_pdf_vv_parameters(SHOCPDFvvparamData &d);
void shoc_assumed_pdf_thl_parameters(SHOCPDFthlparamData &d);
void shoc_assumed_pdf_qw_parameters(SHOCPDFqwparamData &d);
void shoc_assumed_pdf_inplume_correlations(SHOCPDFinplumeData &d);
void shoc_assumed_pdf_compute_temperature(SHOCPDFcomptempData &d);
void shoc_assumed_pdf_compute_qs(SHOCPDFcompqsData &d);
void shoc_assumed_pdf_compute_s(SHOCPDFcompsData &d);
void shoc_assumed_pdf_compute_sgs_liquid(SHOCPDFcompsgsliqData &d);
void shoc_assumed_pdf_compute_cloud_liquid_variance(SHOCPDFcompcloudvarData &d);
void shoc_assumed_pdf_compute_liquid_water_flux(SHOCPDFcompliqfluxData &d);
void shoc_assumed_pdf_compute_buoyancy_flux(SHOCPDFcompbuoyfluxData &d);
void shoc_diag_second_moments_ubycond(SHOCSecondMomentUbycondData& d);

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
