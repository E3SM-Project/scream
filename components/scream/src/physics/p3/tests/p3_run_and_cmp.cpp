#include "share/scream_session.hpp"
#include "share/util/file_utils.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"
#include "share/scream_assert.hpp"

#include "physics/p3/p3_f90.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/p3/p3_ic_cases.hpp"

#include "share/simpleio/simple_io_mod.hpp"

#include <vector>

namespace {
using namespace scream;
using namespace scream::util;
using namespace scream::p3;
using namespace scream::simpleio;

template <typename Scalar>
static Int compare (const std::string& label, const Scalar* a,
                    const Scalar* b, const Int& n, const Real& tol) {
  Int nerr = 0;
  Real den = 0;
  for (Int i = 0; i < n; ++i)
    den = std::max(den, std::abs(a[i]));
  Real worst = 0;
  for (Int i = 0; i < n; ++i) {
    if (std::isnan(a[i]) || std::isinf(a[i]) ||
        std::isnan(b[i]) || std::isinf(b[i])) {
      ++nerr;
      continue;
    }
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
    std::string setname = "";
    Int myit = 0;
    for (const bool log_predictNc : {true, false}) {
      for (const int it : {1, 2}) {
        params_.push_back({ic::Factory::mixed, 1800, it, log_predictNc, setname + std::to_string(myit)});
        myit+=1;
      }
    }
  }

  Int generate_baseline (const std::string& filename, bool use_fortran) {
    Int nerr = 0;
    for (auto ps : params_) {
      // Run reference p3 on this set of parameters.
      const auto d = ic::Factory::create(ps.ic, ic_ncol);
      set_params(ps, *d);
      p3_init(use_fortran);
      p3_main(*d);
      // Add simple io stuff
      std::string mode("baseline_");
      int ncid = open_nc_file(filename,mode + ps.setname,"write");
      reg_fields_nc(ncid,d);
      write_fields_nc(ncid,d);
    }
    
    return nerr;
  }

  Int run_and_cmp (const std::string& filename, const double& tol, bool use_fortran) {
    Int nerr = 0, ne;
    for (auto ps : params_) {
      // Read the reference impl's data from the baseline file.
      const auto d_ref = ic::Factory::create(ps.ic, ic_ncol);
      set_params(ps, *d_ref);
      std::string mode("baseline_");
      auto ncid_in = open_nc_file(filename,mode + ps.setname,"read");
      read_fields_nc(ncid_in,d_ref);
      // Now run a sequence of other impls. This includes the reference
      // implementation b/c it's likely we'll want to change it as we go.
      {
        const auto d = ic::Factory::create(ps.ic, ic_ncol);
        set_params(ps, *d);
        p3_init(use_fortran);
        p3_main(*d);
        ne = compare("ref", tol, d_ref, d);
        if (ne) std::cout << "Ref impl failed.\n";
        nerr += ne;
        // Write output to netCDF to be used with cprnc
        std::string mode("test_");
        auto ncid = open_nc_file(filename,mode + ps.setname,"write");
        reg_fields_nc(ncid,d);
        write_fields_nc(ncid,d);
      }
    }
    return nerr;
  }

private:
  static Int ic_ncol;

  struct ParamSet {
    ic::Factory::IC ic;
    Real dt;
    Int it;
    bool log_predictNc;
    std::string setname;
  };

  static void set_params (const ParamSet& ps, FortranData& d) {
    d.dt = ps.dt;
    d.it = ps.it;
    d.log_predictNc = ps.log_predictNc;
  }

  std::vector<ParamSet> params_;

  /*----------------------------------------------------------------------*/
  // Simple io routines:
  int  open_nc_file(const std::string filename, const std::string runtype, const std::string mode) {
      const int ndims = 4;
      std::string dimnames[] = {"column","level","ilevel","fields"};
      std::string nc_filename = "";
      nc_filename.append(filename.c_str());
      nc_filename.append("_");
      nc_filename.append(runtype);
      nc_filename.append(".nc");
      int dimrng[] = {ic_ncol,72,73,49};
      int ncid = -999;
      if (mode == "write") {
        ncid = init_output1(nc_filename, ndims, dimnames, dimrng);
      } else if (mode == "read") {
        ncid = init_input(nc_filename);
      }
      scream_require_msg( ncid != -999, "Error opening netcdf file " << nc_filename);
      return ncid;
  }
  static void reg_fields_nc (const Int ncid, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    std::string dimnames[] = {"column","level","ilevel","fields"};
    std::string units="unitless";
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      if ( (f.extent[1]==72) && (f.extent[2]>1) ) {
        dimnames[1] = "level";
        dimnames[2] = "fields";
        regfield(ncid,f.name,NC_REAL,3,dimnames,units);
      } else if (f.extent[1]==72) {
        dimnames[1] = "level";
        regfield(ncid,f.name,NC_REAL,2,dimnames,units);
      } else if (f.extent[1]==73) {
        dimnames[1] = "ilevel";
        regfield(ncid,f.name,NC_REAL,2,dimnames,units);
      } else {
        regfield(ncid,f.name,NC_REAL,1,dimnames,units);
      }
    }
    init_output2(ncid);
  }
  static void write_fields_nc (const Int ncid, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      writefield(ncid,f.name,*f.data,-1);
    }
    finalize_io(ncid);
  }
  static void read_fields_nc (const Int ncid, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      readfield(ncid,f.name,*f.data,-1);
    }
    finalize_io(ncid);
  }
  /*----------------------------------------------------------------------*/
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

Int Baseline::ic_ncol = 3;

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
      "  -f        Use fortran impls instead of c++.\n"
      "  -t <tol>  Tolerance for relative error.\n";
    return 1;
  }

  bool generate = false, use_fortran = false;
  scream::Real tol = 0;
  for (int i = 1; i < argc-1; ++i) {
    if (util::eq(argv[i], "-g", "--generate")) generate = true;
    if (util::eq(argv[i], "-f", "--fortran")) use_fortran = true;
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
      nerr += bln.generate_baseline(baseline_fn, use_fortran);
    } else {
      printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
      nerr += bln.run_and_cmp(baseline_fn, tol, use_fortran);
    }
    P3GlobalForFortran::deinit();
  } scream::finalize_scream_session();

  return nerr != 0 ? 1 : 0;
}
