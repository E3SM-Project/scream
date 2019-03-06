#include "share/scream_session.hpp"
#include "share/util/file_utils.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"
#include "share/scream_assert.hpp"

#include "physics/p3/p3_f90.hpp"
#include "physics/p3/p3_ic_cases.hpp"

#include <vector>

namespace {
using namespace scream;
using namespace scream::util;
using namespace scream::p3;

template <typename Scalar>
static Int compare (const std::string& label, const Scalar* a,
                    const Scalar* b, const Int& n, const Real& tol) {
  Int nerr = 0;
  Real den = 0;
  for (Int i = 0; i < n; ++i)
    den = std::max(den, std::abs(a[i]));
  Real worst = 0;
  for (Int i = 0; i < n; ++i) {
    const auto num = std::abs(a[i] - b[i]);
    if (num > tol*den) {
      ++nerr;
      worst = std::max(worst, num);
    }
  }
  if (nerr)
    std::cout << label << " nerr " << nerr << " worst " << (worst/den)
              << " with denominator " << den << "\n";
  return nerr;
}

Int compare (const std::string& label, const double& tol,
             const FortranData::Ptr& ref, const FortranData::Ptr& d) {
  Int nerr = 0;
  FortranDataIterator refi(ref), di(d);
  scream_assert(refi.nfield() == di.nfield());
  for (Int i = 0, n = refi.nfield(); i < n; ++i) {
    const auto& fr = refi.getfield(i);
    const auto& fd = di.getfield(i);
    scream_assert(fr.size == fd.size);
    nerr += compare(label + std::string(" ") + fr.name,
                    fr.data, fd.data, fr.size, tol);
  }
  return nerr;
}

struct Baseline {
  Baseline () {
    params_.push_back({ic::Factory::mixed, 1800});
  }

  Int generate_baseline (const std::string& filename) {
    auto fid = FILEPtr(fopen(filename.c_str(), "w"));
    scream_require_msg( fid, "generate_baseline can't write " << filename);
    Int nerr = 0;
    bool first = true;
    for (auto ps : params_) {
      // Run reference p3 on this set of parameters.
      const auto d = ic::Factory::create(ps.ic);
      d->dt = ps.dt;
      p3_init();
      p3_main(*d);
      if (first) {
        first = false; 
        nerr += check_against_python(*d);
        if (nerr) std::cout << "Spot-check against Python failed.\n";
      }
      // Save the fields to the baseline file.
      write(fid, d);
    }
    return nerr;
  }

  Int run_and_cmp (const std::string& filename, const double& tol) {
    auto fid = FILEPtr(fopen(filename.c_str(), "r"));
    scream_require_msg( fid, "generate_baseline can't read " << filename);
    Int nerr = 0, ne;
    for (auto ps : params_) {
      // Read the reference impl's data from the baseline file.
      const auto d_ref = ic::Factory::create(ps.ic);
      d_ref->dt = ps.dt;
      read(fid, d_ref);
      // Now run a sequence of other impls. This includes the reference
      // implementation b/c it's likely we'll want to change it as we go.
      {
        const auto d = ic::Factory::create(ps.ic);
        p3_init();
        p3_main(*d);
        ne = compare("ref", tol, d_ref, d);
        if (ne) std::cout << "Ref impl failed.\n";
        nerr += ne;
      }
    }
    return nerr;
  }

private:
  struct ParamSet {
    ic::Factory::IC ic;
    Real dt;
  };
  
  std::vector<ParamSet> params_;

  static void write (const FILEPtr& fid, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      util::write(&f.dim, 1, fid);
      util::write(f.extent, f.dim, fid);
      util::write(f.data, f.size, fid);    
    }
  }

  static void read (const FILEPtr& fid, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      int dim, ds[3];
      util::read(&dim, 1, fid);
      scream_require_msg(dim == f.dim,
                      "For field " << f.name << " read expected dim " <<
                      f.dim << " but got " << dim);
      util::read(ds, dim, fid);
      for (int i = 0; i < dim; ++i)
        scream_require_msg(ds[i] == f.extent[i],
                        "For field " << f.name << " read expected dim "
                        << i << " to have extent " << f.extent[i] << " but got "
                        << ds[i]);
      util::read(f.data, f.size, fid);
    }
  }
};

void expect_another_arg (int i, int argc) {
  scream_require_msg(i != argc-1, "Expected another cmd-line arg.");
}

} // namespace anon

int main (int argc, char** argv) {

  int nerr = 0;

  if (argc == 1) {
    std::cout <<
      argv[0] << " [options] baseline-filename\n"
      "Options:\n"
      "  -g        Generate baseline file.\n"
      "  -t <tol>  Tolerance for relative error.\n";
    return -1;
  }

  bool generate = false;
  scream::Real tol = 0;
  for (int i = 1; i < argc-1; ++i) {
    if (util::eq(argv[i], "-g", "--generate")) generate = true;
    if (util::eq(argv[i], "-t", "--tol")) {
      expect_another_arg(i, argc);
      ++i;
      tol = std::atof(argv[i]);
    }
  }

  // Decorate baseline name with precision.
  std::string baseline_fn(argv[argc-1]);
  baseline_fn += std::to_string(sizeof(scream::Real));

  scream::initialize_scream_session(argc, argv); {
    Baseline bln;
    if (generate) {
      std::cout << "Generating to " << baseline_fn << "\n";
      nerr += bln.generate_baseline(baseline_fn);
    } else {
      printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
      nerr += bln.run_and_cmp(baseline_fn, tol);
    }    
  } scream::finalize_scream_session();

  return nerr;
}
