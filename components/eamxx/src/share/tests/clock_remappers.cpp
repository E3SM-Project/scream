#include <catch2/catch.hpp>

#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/refining_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_timing.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include <ekat/util/ekat_test_utils.hpp>

namespace scream {

template<typename ViewT>
typename ViewT::HostMirror
cmvc (const ViewT& v) {
  auto vh = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(vh,v);
  return vh;
}

Field create_field (const std::string& name, const LayoutType lt, const AbstractGrid& grid)
{
  const auto u = ekat::units::Units::nondimensional();
  const auto CMP = ShortFieldTagsNames::CMP;
  const auto& gn = grid.name();
  const auto  ndims = 2;
  Field f;
  switch (lt) {
    case LayoutType::Scalar2D:
      f = Field(FieldIdentifier(name,grid.get_2d_scalar_layout(),u,gn));  break;
    case LayoutType::Vector2D:
      f = Field(FieldIdentifier(name,grid.get_2d_vector_layout(CMP,ndims),u,gn));  break;
    case LayoutType::Scalar3D:
      f = Field(FieldIdentifier(name,grid.get_3d_scalar_layout(true),u,gn));  break;
    case LayoutType::Vector3D:
      f = Field(FieldIdentifier(name,grid.get_3d_vector_layout(false,CMP,ndims),u,gn));  break;
    default:
      EKAT_ERROR_MSG ("Invalid layout type for this unit test.\n");
  }
  f.allocate_view();

  return f;
}

template<typename Engine>
Field create_field (const std::string& name, const LayoutType lt, const AbstractGrid& grid, Engine& engine) {
  auto f = create_field(name,lt,grid);

  // Use discrete_distribution to get an integer, then use that as exponent for 2^-n.
  // This guarantees numbers that are exactly represented as FP numbers, which ensures
  // the test will produce the expected answer, regardless of how math ops are performed.
  using IPDF = std::discrete_distribution<int>;
  IPDF ipdf ({1,1,1,1,1,1,1,1,1,1});
  auto pdf = [&](Engine& e) {
    return std::pow(2,ipdf(e));
  };
  randomize(f,engine,pdf);

  return f;
}

template<typename RemapperT,typename EngineT>
void run_remapper (const std::shared_ptr<AbstractGrid>& fine_grid,
                   const std::string& map_file,
                   EngineT& engine, int ntests)
{
  const auto& comm = fine_grid->get_comm();
  if (comm.am_i_root()) {
    printf("running tests with map file: %s\n",map_file.c_str());
  }
  auto remapper = std::make_shared<RemapperT>(fine_grid,map_file);

  auto src_grid = remapper->get_src_grid();
  auto tgt_grid = remapper->get_tgt_grid();

  auto bundle_src = create_field("bundle3d_src",LayoutType::Vector3D,*src_grid,engine);
  auto s2d_src   = create_field("s2d_src",LayoutType::Scalar2D,*src_grid,engine);
  auto v2d_src   = create_field("v2d_src",LayoutType::Vector2D,*src_grid,engine);
  auto s3d_src   = create_field("s3d_src",LayoutType::Scalar3D,*src_grid,engine);
  auto v3d_src   = create_field("v3d_src",LayoutType::Vector3D,*src_grid,engine);

  auto bundle_tgt = create_field("bundle3d_tgt",LayoutType::Vector3D,*tgt_grid);
  auto s2d_tgt   = create_field("s2d_tgt",LayoutType::Scalar2D,*tgt_grid);
  auto v2d_tgt   = create_field("v2d_tgt",LayoutType::Vector2D,*tgt_grid);
  auto s3d_tgt   = create_field("s3d_tgt",LayoutType::Scalar3D,*tgt_grid);
  auto v3d_tgt   = create_field("v3d_tgt",LayoutType::Vector3D,*tgt_grid);

  remapper->registration_begins();
  remapper->register_field(s2d_src,s2d_tgt);
  remapper->register_field(v2d_src,v2d_tgt);
  remapper->register_field(s3d_src,s3d_tgt);
  remapper->register_field(v3d_src,v3d_tgt);
  remapper->register_field(bundle_src.get_component(0),bundle_tgt.get_component(0));
  remapper->register_field(bundle_src.get_component(1),bundle_tgt.get_component(1));
  remapper->registration_ends();

  for (int i=0; i<ntests; ++i) {
    if (comm.am_i_root()) {
      printf("  itest = %d\n",i);
    }
    remapper->remap(true);
  }
}

TEST_CASE ("refining_remapper") {

  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = setup_random_test (&comm);

  MPI_Fint fcomm = MPI_Comm_c2f(comm.mpi_comm());
  scorpio::eam_init_pio_subsystem(fcomm);

  bool was_already_inited;
  init_gptl(was_already_inited);
  
  const auto& test_params = ekat::TestSession::get().params;
  const auto& coarsen_map = test_params.at("coarsen_map_file");
  const auto& refine_map  = test_params.at("refine_map_file");
  const int ntests = test_params.find("ntests")==test_params.end() ? 10 : std::stoi(test_params.at("ntests"));

  // Create fine grid.
  const int nlevs = 128;
  const int ngdofs_fine   = scorpio::get_dimlen(coarsen_map,"n_a");
  auto fine_grid   = create_point_grid("fine",ngdofs_fine,nlevs,comm);

  run_remapper<RefiningRemapper>(fine_grid,refine_map,engine,ntests);
  run_remapper<CoarseningRemapper>(fine_grid,coarsen_map,engine,ntests);

  write_timers_to_file (comm,"coarsen-vs-refine.txt");

  // Clean up
  if (not was_already_inited) {
    finalize_gptl();
  }
  scorpio::eam_pio_finalize();
}

} // namespace scream

