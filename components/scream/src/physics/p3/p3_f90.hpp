#ifndef SCREAM_P3_F90_HPP
#define SCREAM_P3_F90_HPP

#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

#include <memory>
#include <vector>

namespace scream {
namespace p3 {

// Data format we can use to communicate with Fortran version.
struct FortranData {
  typedef std::shared_ptr<FortranData> Ptr;

  using KT     = KokkosTypes<HostDevice>;
  using Scalar = Real;

  using Array1 = typename KT::template lview<Scalar*>;
  using Array2 = typename KT::template lview<Scalar**>;
  using Array3 = typename KT::template lview<Scalar***>;

  bool log_predictNc;
  const Int ncol, nlev;

  // In
  Real dt;
  Int it;
  Array2 qv, th, pres, dzq, npccn, naai, qc_relvar, qc, nc, qr, nr,  qitot,
    nitot, qirim, birim, pdel, exner;
  // Out
  Array1 prt_liq, prt_sol;
  Array2 diag_ze, diag_effc, diag_effi, diag_vmi, diag_di, diag_rhoi, cmeiout, prain, nevapr, prer_evap, rflx, sflx, rcldm, lcldm, icldm;
  Array2 pratot, prctot;
  Array3 p3_tend_out;
  Array2 mu_c, lamc;
  Array2 liq_ice_exchange,vap_liq_exchange,vap_ice_exchange,vap_cld_exchange;

  FortranData(Int ncol, Int nlev);
};

// Iterate over a FortranData's arrays. For examples, see Baseline::write, read.
struct FortranDataIterator {
  struct RawArray {
    std::string name;
    Int dim;
    Int extent[3];
    FortranData::Scalar* data;
    FortranData::Array1::size_type size;
  };

  explicit FortranDataIterator(const FortranData::Ptr& d);

  Int nfield () const { return fields_.size(); }
  const RawArray& getfield(Int i) const;

private:
  FortranData::Ptr d_;
  std::vector<RawArray> fields_;

  void init(const FortranData::Ptr& d);
};

void p3_init(bool use_fortran=false);
void p3_main(const FortranData& d);

// We will likely want to remove these checks in the future, as we're not tied
// to the exact implementation or arithmetic in P3. For now, these checks are
// here to establish that the initial regression-testing code gives results that
// match the python f2py tester, without needing a data file.
Int check_against_python(const FortranData& d);

int test_FortranData();
int test_p3_init(bool use_fortran);
int test_p3_ic(bool use_fortran);

}  // namespace p3
}  // namespace scream

#endif
