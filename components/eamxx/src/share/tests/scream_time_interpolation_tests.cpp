#include <catch2/catch.hpp>

#include "share/io/scream_output_manager.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

#include "share/util/scream_time_interpolation.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"

#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/util/ekat_test_utils.hpp"
namespace scream {

constexpr Real tol               = 1e-8;
constexpr int num_output_steps   = 5;
constexpr int dt_data            = 3600;
constexpr int freq               = 3;
const     std::string freq_units = "nsteps";

// Subclass of OutputManager that will is needed to track the set of output files
// written in the test setup.
class OutputManager4Test : public scream::OutputManager
{
public:
  OutputManager4Test()
   : OutputManager()
  {
    // Do Nothing
  }

  void update_file_list() {
    if (std::find(m_list_of_files.begin(),m_list_of_files.end(), m_output_file_specs.filename) == m_list_of_files.end()) {
      m_list_of_files.push_back(m_output_file_specs.filename);
    }
  }
  std::vector<std::string> get_files() { return m_list_of_files; }

protected:
  std::vector<std::string> m_list_of_files;
};

// We use a linear function to set test data to make the expected result
// of time interpolation independent of data frequency.
Real test_func(const Real t, const Real y0) 
{
  return 2*t + y0;
}

void update_field (const Field& f, const Field& f_0, const Real t) {
  auto data     = f.get_internal_view_data<Real,Host>();
  auto data_0   = f_0.get_internal_view_data<Real,Host>();
  auto nscalars = f.get_header().get_alloc_properties().get_num_scalars();
  for (int i=0; i<nscalars; ++i) {
    data[i] = test_func(t,data_0[i]);
  }
}

util::TimeStamp get_t0 () {
  return util::TimeStamp({2023,5,4},{0,0,0});
}

std::shared_ptr<const GridsManager>
get_gm (const ekat::Comm& comm)
{
  const int nlcols = 3;
  const int nlevs = 4;
  const int ngcols = nlcols*comm.size();
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",ngcols);
  gm_params.set("number_of_vertical_levels",nlevs);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();
  return gm;
}

std::shared_ptr<FieldManager>
get_fm (const std::shared_ptr<const AbstractGrid>& grid,
        const util::TimeStamp& t0, const int seed)
{
  using FL  = FieldLayout;
  using FID = FieldIdentifier;
  using namespace ShortFieldTagsNames;

  // Random number generation stuff
  // NOTES
  //  - Use integers, so we can check answers without risk of
  //    non bfb diffs due to different order of sums.
  //  - Uniform_int_distribution returns an int, and the randomize
  //    util checks that return type matches the Field data type.
  //    So wrap the int pdf in a lambda, that does the cast.
  std::mt19937_64 engine(seed); 
  auto my_pdf = [&](std::mt19937_64& engine) -> Real {
    std::uniform_int_distribution<int> pdf (0,100);
    Real v = pdf(engine);
    return v;
  };

  const int nlcols = grid->get_num_local_dofs();
  const int nlevs  = grid->get_num_vertical_levels();

  std::vector<FL> layouts =
  {
    FL({COL         }, {nlcols        }),
    FL({COL,     LEV}, {nlcols,  nlevs}),
    FL({COL,CMP,ILEV}, {nlcols,2,nlevs+1})
  };

  auto fm = std::make_shared<FieldManager>(grid);
  fm->registration_begins();
  fm->registration_ends();
  
  const auto units = ekat::units::Units::nondimensional();
  for (const auto& fl : layouts) {
    FID fid("f_"+std::to_string(fl.size()),fl,units,grid->name());
    Field f(fid);
    f.allocate_view();
    randomize (f,engine,my_pdf);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }

  return fm;
}

// Returns fields after initialization
void write_test_files (const std::shared_ptr<const GridsManager>& gm, const std::shared_ptr<FieldManager>& fm, const int seed, const ekat::Comm& comm, std::vector<std::string>& list_of_files)
{

  const int dt = dt_data;
  // Time advance parameters
  auto t0 = get_t0();

  auto grid = gm->get_grid("Point Grid");
  auto fm_0 = get_fm(grid,t0,seed); // should match fm

  // Create some fields
  std::vector<std::string> fnames;
  for (auto it : *fm) {
    fnames.push_back(it.second->name());
  }

  // Create output params
  ekat::ParameterList om_pl;
  om_pl.set("MPI Ranks in Filename",true);
  om_pl.set("filename_prefix",std::string("time_interp_source"));
  om_pl.set("Field Names",fnames);
  om_pl.set("Averaging Type", std::string("INSTANT"));
  om_pl.set("Max Snapshots Per File",num_output_steps);
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",freq_units);
  ctrl_pl.set("Frequency",freq);
  ctrl_pl.set("MPI Ranks in Filename",true);
  ctrl_pl.set("save_grid_data",false);

  // Create Output manager
  OutputManager4Test om;
  om.setup(comm,om_pl,fm,gm,t0,t0,false);

  // Time loop: ensure we always hit at least 2 files
  const int nsteps = num_output_steps*freq*2-1;
  auto t = t0;
  for (int n=0; n<nsteps; ++n) {
    // Update time
    t += dt;

    // Add to all fields entries
    for (const auto& n : fnames) {
      auto f   = fm->get_field(n);
      auto f_0 = fm_0->get_field(n);
      auto tt  = t-t0;
      update_field(f,f_0,tt);
    }

    // Run output manager
    om.run (t);
    om.update_file_list();
  }

  // Close file and cleanup
  om.finalize();

  list_of_files = om.get_files();
}


TEST_CASE ("scream_time_interpolation") {

  using namespace util;

  // Set up basics of the tests
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::eam_init_pio_subsystem(comm);
  auto seed = get_random_test_seed(&comm);

  // Create grid
  auto gm = get_gm(comm);
  auto grid = gm->get_grid("Point Grid");

  // Create fields for testing
  // We create two field managers, one for writing output, the other to track
  // for testing.  Note that get_fm will produce identical fields for the same seed
  // so at this point fm_write and fm_track match.  The function write_test_files
  // will update fm_write so we need the copy for testing.
  auto t0 = get_t0();
  auto fm_write     = get_fm(grid,t0,seed); // TODO: probably don't need this.
  auto fm_t0        = get_fm(grid,t0,seed);
  auto fm_reference = get_fm(grid,t0,seed); 
  auto fm_track     = get_fm(grid,t0,seed);
  std::vector<std::string> files;
  write_test_files(gm,fm_write,seed,comm,files);

  std::vector<std::string> fnames;
  for (auto it : *fm_reference) {
    fnames.push_back(it.second->name());
  }

  // Create Interpolator object
  auto time_interp = TimeInterpolation(grid);

  // Set list of variables and files to feed to time_interp
  const auto num_of_fields = fm_track->size();
  for (auto ii = fm_track->begin(); ii != fm_track->end(); ii++) {
    auto ff = ii->second;
    time_interp.add_field(*ff);
  }
  // List files in backwards order to force sorting.
  std::reverse(files.begin(),files.end());
  time_interp.init(t0,files);

  // Now run for each timestep and interpolate.  The interpolated fields
  // are stored in fm_track, we can compare them to the test function (test_func),
  // which is linear and thus should match to within machine precision.
  const int nsteps = num_output_steps*freq*2-1;
  auto t = t0;
  const int dt = dt_data/freq;
  for (int nn=0; nn<nsteps; nn++) {
    t += dt;
    auto interp_fields = time_interp.perform_time_interpolation(t);
    for (const auto& n : fnames) {
      auto f   = fm_reference->get_field(n);
      auto f_0 = fm_t0->get_field(n);
      auto f_t = interp_fields.at(n);
      auto tt  = t-t0;
      update_field(f,f_0,tt);
      f_0.sync_to_host();
      f_t.sync_to_host();
      bool same_locally = true;
      const auto& l1 = f_0.get_header().get_identifier().get_layout();
      const auto& dims = l1.dims();
      switch (l1.rank()) {
        case 1:
          {
            auto v1 = f_0.template get_view<Real*,Host>();
            auto v2 = f_t.template get_view<Real*,Host>();
            for (int i=0; i<dims[0]; ++i) {
              if (std::abs(v1(i) - v2(i)) < tol) {
                same_locally = false;
                break;
              }
            }
          }
          break;
        case 2:
          {
            auto v1 = f_0.template get_view<Real**,Host>();
            auto v2 = f_t.template get_view<Real**,Host>();
            for (int i=0; same_locally && i<dims[0]; ++i) {
              for (int j=0; j<dims[1]; ++j) {
                if (std::abs(v1(i,j) - v2(i,j)) < tol) {
                  same_locally = false;
                  break;
                }
            }}
          }
          break;
        case 3:
          {
            auto v1 = f_0.template get_view<Real***,Host>();
            auto v2 = f_t.template get_view<Real***,Host>();
            for (int i=0; same_locally && i<dims[0]; ++i) {
              for (int j=0; same_locally && j<dims[1]; ++j) {
                for (int k=0; k<dims[2]; ++k) {
                  if (std::abs(v1(i,j,k) - v2(i,j,k)) < tol) {
                    same_locally = false;
                    break;
                  }
            }}}
          }
          break;
        default:
          EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
      }
      REQUIRE(same_locally);
    }
  }


  // All done with IO
  scorpio::eam_pio_finalize();
}

} // namespace scream
