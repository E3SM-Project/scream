#include "p3_f90.hpp"
#include "p3_constants.hpp"
#include "p3_ic_cases.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"

using scream::Real;
using scream::Int;
extern "C" {
  void micro_p3_utils_init_c(Real Cpair, Real Rair, Real RH2O, Real RhoH2O, 
                 Real MWH2O, Real MWdry, Real gravit, Real LatVap, Real LatIce, 
                 Real CpLiq, Real Tmelt, Real Pi, Int iulog, bool masterproc);
  void p3_init_c(const char** lookup_file_dir, int* info);
  void p3_main_c(Real* qc, Real* nc, Real* qr, Real* nr, Real* th_old, Real* th,
                 Real* qv_old, Real* qv, Real dt, Real* qitot, Real* qirim,
                 Real* nitot, Real* birim, Real* ssat, Real* pres,
                 Real* dzq, Real* npccn, Real* naai, Int it, Real* prt_liq, Real* prt_sol, Int its,
                 Int ite, Int kts, Int kte, Real* diag_ze,
                 Real* diag_effc, Real* diag_effi, Real* diag_vmi,
                 Real* diag_di, Real* diag_rhoi,
                 bool log_predictNc,
                 Real* pdel, Real* exner, Real* cmeiout, Real* prain,
                 Real* nevapr, Real* prer_evap,
                 Real* rflx, Real* sflx, // 1 extra column size
                 Real* rcldm, Real* lcldm, Real* icldm, Real* p3_tend_out);
}

namespace scream {
namespace p3 {

FortranData::FortranData (Int ncol_, Int nlev_)
  : ncol(ncol_), nlev(nlev_)
{

  dt = -1; // model time step, s; set to invalid -1
  it = 1;  // seems essentially unused
  // In/out
  qc = Array2("cloud liquid water mixing ratio, kg/kg", ncol, nlev);
  nc = Array2("cloud liquid drop number, #/kg", ncol, nlev);
  qr = Array2("rain water mixing ratio, kg/kg", ncol, nlev);
  nr = Array2("rain drop number, #/kg", ncol, nlev);
  qitot = Array2("total ice mass mixing ratio, kg/kg", ncol, nlev);
  nitot = Array2("total ice number, #/kg", ncol, nlev);
  qirim = Array2("rime ice mass mixing ratio, kg/kg", ncol, nlev);
  birim = Array2("rime ice volume mixing ratio, m3/kg", ncol, nlev);
  ssat = Array2("supersaturation (qv - qs), kg/kg", ncol, nlev);
  qv = Array2("water vapor mixing ratio, kg/kg", ncol, nlev);
  th = Array2("potential temperature, K", ncol, nlev);
  qv_old = Array2("qv at beginning of timestep, kg/kg", ncol, nlev);
  th_old = Array2("theta at beginning of timestep, K", ncol, nlev);
  pres = Array2("pressure, Pa", ncol, nlev);
  dzq = Array2("vertical grid spacing, m", ncol, nlev);
  npccn = Array2("ccn activated number tendency, kg-1 s-1", ncol, nlev);
  naai = Array2("activated nuclei concentration, kg-1", ncol, nlev);
  pdel = Array2("pressure thickness, Pa", ncol, nlev);
  exner = Array2("Exner expression", ncol, nlev);
  // Out
  prt_liq = Array1("precipitation rate, liquid  m/s", ncol);
  prt_sol = Array1("precipitation rate, solid   m/s", ncol);
  diag_ze = Array2("equivalent reflectivity, dBZ", ncol, nlev);
  diag_effc = Array2("effective radius, cloud, m", ncol, nlev);
  diag_effi = Array2("effective radius, ice, m", ncol, nlev);
  diag_vmi = Array2("mass-weighted fall speed of ice, m/s", ncol, nlev);
  diag_di = Array2("mean diameter of ice, m", ncol, nlev);
  diag_rhoi = Array2("bulk density of ice, kg/m", ncol, nlev);
  cmeiout = Array2("qitend due to deposition/sublimation ", ncol, nlev);
  prain = Array2("Total precipitation (rain + snow)", ncol, nlev);
  nevapr = Array2("evaporation of total precipitation (rain + snow)", ncol, nlev);
  prer_evap = Array2("evaporation of rain", ncol, nlev);
  rflx = Array2("grid-box average rain flux (kg m^-2 s^-1), pverp", ncol, nlev+1);
  sflx = Array2("grid-box average ice/snow flux (kg m^-2 s^-1), pverp", ncol, nlev+1);
  rcldm = Array2("Rain cloud fraction", ncol, nlev);
  lcldm = Array2("Liquid cloud fraction", ncol, nlev);
  icldm = Array2("Ice cloud fraction", ncol, nlev);
  p3_tend_out = Array3("Microphysics Tendencies", ncol, nlev, 49);
}

FortranDataIterator::FortranDataIterator (const FortranData::Ptr& d) {
  init(d);
}

void FortranDataIterator::init (const FortranData::Ptr& dp) {
  d_ = dp;
#define fdipb(name)                                                     \
  fields_.push_back({#name,                                             \
        2,                                                              \
        {d_->name.extent_int(0), d_->name.extent_int(1), d_->name.extent_int(2)}, \
        d_->name.data(),                                                \
        d_->name.size()})
  fdipb(qv); fdipb(th); fdipb(qv_old); fdipb(th_old); fdipb(pres);
  fdipb(dzq); fdipb(npccn); fdipb(naai); fdipb(qc); fdipb(nc); fdipb(qr); fdipb(nr);
  fdipb(ssat); fdipb(qitot); fdipb(nitot);
  fdipb(qirim); fdipb(birim); fdipb(prt_liq); fdipb(prt_sol);
  fdipb(diag_ze); fdipb(diag_effc); fdipb(diag_effi);
  fdipb(diag_vmi); fdipb(diag_di); fdipb(diag_rhoi);
  fdipb(pdel); fdipb(exner); fdipb(cmeiout); fdipb(prain);
  fdipb(nevapr); fdipb(prer_evap);
  fdipb(rflx); fdipb(sflx);
  fdipb(rcldm); fdipb(lcldm); fdipb(icldm); fdipb(p3_tend_out);
#undef fdipb
}

const FortranDataIterator::RawArray&
FortranDataIterator::getfield (Int i) const {
  scream_assert(i >= 0 || i < nfield());
  return fields_[i];
}

void micro_p3_utils_init () {
  using c = Constants<Real>;
  micro_p3_utils_init_c(c::Cpair, c::Rair, c::RH2O, c::RhoH2O, 
                 c::MWH2O, c::MWdry, c::gravit, c::LatVap, c::LatIce, 
                 c::CpLiq, c::Tmelt, c::Pi, c::iulog, c::masterproc);
}

void p3_init () {
  micro_p3_utils_init();
  static const char* dir = ".";
  Int info;
  p3_init_c(&dir, &info);
  scream_require_msg(info == 0, "p3_init_c returned info " << info);
}

void p3_main (const FortranData& d) {
  p3_main_c(d.qc.data(), d.nc.data(), d.qr.data(), d.nr.data(), d.th_old.data(),
            d.th.data(), d.qv_old.data(), d.qv.data(), d.dt, d.qitot.data(),
            d.qirim.data(), d.nitot.data(), d.birim.data(), d.ssat.data(),
            d.pres.data(), d.dzq.data(), d.npccn.data(), d.naai.data(), d.it, d.prt_liq.data(),
            d.prt_sol.data(), 1, d.ncol, 1, d.nlev, d.diag_ze.data(),
            d.diag_effc.data(), d.diag_effi.data(), d.diag_vmi.data(),
            d.diag_di.data(), d.diag_rhoi.data(),
            d.log_predictnc,
            d.pdel.data(), d.exner.data(), d.cmeiout.data(), d.prain.data(),
            d.nevapr.data(), d.prer_evap.data(),
            d.rflx.data(), d.sflx.data(),
            d.rcldm.data(), d.lcldm.data(), d.icldm.data(),d.p3_tend_out.data());
}

Int check_against_python (const FortranData& d) {
  Int nerr = 0;
  if (util::is_single_precision<Real>::value) {
    const double tol = 0;
    if (util::reldif<double>(d.birim(0,d.nlev-1), 7.237245824853744e-08) > tol)
      ++nerr;
    if (util::reldif<double>(d.qirim(0,d.nlev-1), 9.047746971191373e-06) > tol)
      ++nerr;
    if (util::reldif<double>(d.nr(0,d.nlev-1), 3.177030468750000e+04) > tol)
      ++nerr;
  }
  return nerr;
}

int test_FortranData () {
  FortranData d(11, 72);
  return 0;
}

int test_p3_init () {
  p3_init();
  return 0;
}

int test_p3_ic () {
  const auto d = ic::Factory::create(ic::Factory::mixed);
  p3_init();
  p3_main(*d);
  return check_against_python(*d);
}

} // namespace p3
} // namespace scream
