#ifndef SCREAM_P3_F90_HPP
#define SCREAM_P3_F90_HPP

#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"

#include <memory>
#include <vector>

namespace scream {
namespace p3 {

// Data format we can use to communicate with Fortran version.
struct FortranData {
  typedef std::shared_ptr<FortranData> Ptr;

  typedef Kokkos::HostSpace ExeSpace;
  typedef Kokkos::LayoutLeft Layout;
  typedef Real Scalar;

  using Array1 = Kokkos::View<Scalar*, Layout, ExeSpace>;
  using Array2 = Kokkos::View<Scalar**, Layout, ExeSpace>;
  using Array3 = Kokkos::View<Scalar***, Layout, ExeSpace>;

  static constexpr bool log_predictnc = true;

  const Int ncol, nlev;

  // In
  Real dt;
  Int it;
  Array2 qv, th, qv_old, th_old, pres, dzq, npccn, naai, qc, nc, qr, nr, ssat, qitot, nitot, qirim, birim, pdel, exner, ast;
  // Out
  Array1 prt_liq, prt_sol;
  Array2 diag_ze, diag_effc, diag_effi, diag_vmi, diag_di, diag_rhoi, cmeiout, prain, nevapr, prer_evap, rflx, sflx, rcldm, lcldm, icldm;
  Array3 p3_tend_out;

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

void p3_init();
void p3_main(const FortranData& d);

// We will likely want to remove these checks in the future, as we're not tied
// to the exact implementation or arithmetic in P3. For now, these checks are
// here to establish that the initial regression-testing code gives results that
// match the python f2py tester, without needing a data file.
Int check_against_python(const FortranData& d);

int test_FortranData();
int test_p3_init();
int test_p3_main();
int test_p3_ic();

}  // namespace p3
}  // namespace scream

#endif
