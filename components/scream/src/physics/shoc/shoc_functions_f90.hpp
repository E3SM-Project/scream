#ifndef SCREAM_SHOC_FUNCTIONS_F90_HPP
#define SCREAM_SHOC_FUNCTIONS_F90_HPP

#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

#include "shoc_functions.hpp"

#include <vector>
#include <array>
#include <utility>

//
// Bridge functions to call fortran version of shoc functions from C++
//

namespace scream {
namespace shoc {

// Base class for common SHOC data setup
struct SHOCDataBase
{
  Int shcol, nlev, nlevi;

  SHOCDataBase(Int shcol_, Int nlev_, Int nlevi_,
               const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i);

  SHOCDataBase(const SHOCDataBase &rhs) = delete;
  SHOCDataBase(const SHOCDataBase &rhs, const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i);
  SHOCDataBase &operator=(const SHOCDataBase &rhs);

  template <util::TransposeDirection::Enum D>
  void transpose() {
    std::vector<Real> data(m_data.size());

    // Transpose on the zt grid
    for (size_t i = 0; i < m_ptrs.size(); ++i) {
      util::transpose<D>(*(m_ptrs[i]), data.data() + (m_total*i) , shcol, nlev);
    }

    // Transpose on the zi grid
    for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
      util::transpose<D>(*(m_ptrs_i[i]), data.data() + (m_ptrs.size()*m_total) + (m_totali*i), shcol, nlevi);
    }

    m_data = data;
  }

  void randomize();

  Int total() const { return m_total; }
  Int totali() const { return m_totali; }

 private:
  void init_ptrs();

  // Internals
  Int m_total, m_totali;
  std::vector<Real**> m_ptrs, m_ptrs_i;
  std::vector<Real> m_data;
};

struct SHOCGridData : public SHOCDataBase {
  // Inputs
  Real *zt_grid, *zi_grid, *pdel;

  // In/out
  Real *dz_zt, *dz_zi, *rho_zt;

  SHOCGridData(Int shcol_, Int nlev_, Int nlevi_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&zt_grid, &dz_zt, &pdel, &rho_zt}, {&zi_grid, &dz_zi}) {}

  SHOCGridData(const SHOCGridData &rhs) : SHOCDataBase(rhs, {&zt_grid, &dz_zt, &pdel, &rho_zt}, {&zi_grid, &dz_zi}) {}

  SHOCGridData &operator=(const SHOCGridData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};

struct SHOCColstabData : public SHOCDataBase {
  // Inputs
  Real *dz_zt, *pres, *brunt;

  // Output
  Real *brunt_int;

  SHOCColstabData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 1, {&dz_zt, &pres, &brunt}, {&brunt_int}) {}
  SHOCColstabData(const SHOCColstabData &rhs) : SHOCDataBase(rhs, {&dz_zt, &pres, &brunt}, {&brunt_int}) {}
  SHOCColstabData &operator=(const SHOCColstabData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCColstabData

struct SHOCTkeshearData : public SHOCDataBase {
  // Inputs
  Real *dz_zi, *u_wind, *v_wind;

  // In/out
  Real *sterm;

  //functions to initialize data
  SHOCTkeshearData(Int shcol_, Int nlev_, Int nlevi_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&u_wind, &v_wind}, {&dz_zi, &sterm}) {}
  SHOCTkeshearData(const SHOCTkeshearData &rhs) : SHOCDataBase(rhs, {&u_wind, &v_wind}, {&dz_zi, &sterm}) {}
  SHOCTkeshearData &operator=(const SHOCTkeshearData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCTkeshearData

struct SHOCIsotropicData : public SHOCDataBase {
  // Inputs
  Real *tke, *a_diss, *brunt, *brunt_int;

  // Output
  Real *isotropy;

  //functions to initialize data
  SHOCIsotropicData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 1, {&tke, &a_diss, &brunt, &isotropy}, {&brunt_int}) {}
  SHOCIsotropicData(const SHOCIsotropicData &rhs) : SHOCDataBase(rhs, {&tke, &a_diss, &brunt, &isotropy}, {&brunt_int}) {}
  SHOCIsotropicData &operator=(const SHOCIsotropicData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCIsotropicData

struct SHOCAdvsgstkeData : public SHOCDataBase {
  // Inputs
  Real dtime;
  Real *shoc_mix, *wthv_sec, *sterm_zt, *tk;

  // In/out
  Real *tke;

  // Outputs
  Real *a_diss;

  //functions to initialize data
  SHOCAdvsgstkeData(Int shcol_, Int nlev_, Real dtime_) :
    SHOCDataBase(shcol_, nlev_, 0, {&shoc_mix, &wthv_sec, &sterm_zt, &tk, &tke, &a_diss}, {}), dtime(dtime_) {}
  SHOCAdvsgstkeData(const SHOCAdvsgstkeData &rhs) : SHOCDataBase(rhs, {&shoc_mix, &wthv_sec, &sterm_zt, &tk, &tke, &a_diss}, {}), dtime(rhs.dtime) {}
  SHOCAdvsgstkeData &operator=(const SHOCAdvsgstkeData &rhs)
  { SHOCDataBase::operator=(rhs); dtime = rhs.dtime; return *this; }
};//SHOCAdvsgstkeData

struct SHOCEddydiffData : public SHOCDataBase {
  // Inputs
  Real *pblh, *obklen, *zt_grid, *shoc_mix, *sterm_zt,
        *isotropy, *tke;

  // Output
  Real *tk, *tkh;

  //functions to initialize data
  SHOCEddydiffData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 1, {&zt_grid, &shoc_mix, &isotropy, &tke, &tk, &tkh, &sterm_zt}, {&obklen, &pblh}) {}
  SHOCEddydiffData(const SHOCEddydiffData &rhs) : SHOCDataBase(rhs, {&zt_grid, &shoc_mix, &isotropy, &tke, &tk, &tkh, &sterm_zt}, {&obklen, &pblh}) {}
  SHOCEddydiffData &operator=(const SHOCEddydiffData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCEddydiffData

//Create data structure to hold data for calc_shoc_vertflux
struct SHOCVertfluxData : public SHOCDataBase {
  // Inputs
  Real *tkh_zi, *dz_zi, *invar;

  // In/out
  Real *vertflux;

  SHOCVertfluxData(Int shcol_, Int nlev_, Int nlevi_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&invar}, {&tkh_zi, &dz_zi, &vertflux}) {}
  SHOCVertfluxData(const SHOCVertfluxData &rhs) : SHOCDataBase(rhs, {&invar}, {&tkh_zi, &dz_zi, &vertflux}) {}
  SHOCVertfluxData &operator=(const SHOCVertfluxData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
}; //SHOCVertfluxData

struct SHOCVarorcovarData : public SHOCDataBase {
  // Inputs
  Real tunefac;
  Real *tkh_zi, *dz_zi, *isotropy_zi, *invar1, *invar2;

  // In/out
  Real *varorcovar;

  SHOCVarorcovarData(Int shcol_, Int nlev_, Int nlevi_, Real tunefac_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&invar1, &invar2}, {&tkh_zi, &dz_zi, &isotropy_zi, &varorcovar}), tunefac(tunefac_) {}
  SHOCVarorcovarData(const SHOCVarorcovarData &rhs) :
    SHOCDataBase(rhs, {&invar1, &invar2}, {&tkh_zi, &dz_zi, &isotropy_zi, &varorcovar}), tunefac(rhs.tunefac) {}
  SHOCVarorcovarData &operator=(const SHOCVarorcovarData &rhs)
  { SHOCDataBase::operator=(rhs); tunefac = rhs.tunefac; return *this; }
};//SHOCVarorcovarData

struct MomSrfData {
  static constexpr size_t NUM_ARRAYS = 5;

  Int shcol;

  Real *wthl, *uw, *vw, *ustar2, *wstar;

  MomSrfData(Int shcol_, const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  MomSrfData(const MomSrfData& rhs);

  Int nk() const { return m_nk; }

  private:
  // data
  Int m_nk;
  std::vector<Real> m_data;
};

//
// Glue functions to call fortran from from C++ with the Data struct
//

// This function initialzes the grid used by shoc. Given the
// locations of the cell center (location of thermodynaics quantities), cell
// interfaces, and pressure gradient the functon returns dz_zi, dz_zt,
// and density.
void shoc_grid(SHOCGridData &d);
void calc_shoc_vertflux(SHOCVertfluxData &d);
void calc_shoc_varorcovar(SHOCVarorcovarData &d);
void integ_column_stability(SHOCColstabData &d);
void compute_shr_prod(SHOCTkeshearData &d);
void isotropic_ts(SHOCIsotropicData &d);
void adv_sgs_tke(SHOCAdvsgstkeData &d);
void eddy_diffusivities(SHOCEddydiffData &d);
void shoc_diag_second_moments_srf(MomSrfData& d);

//
// _f functions decls
//
extern "C" {

void calc_shoc_vertflux_f(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux);

void shoc_diag_second_moments_srf_f(Int shcol, Real* wthl, Real* uw, Real* vw, Real* ustar2, Real* wstar);

}

}  // namespace shoc
}  // namespace scream

#endif
