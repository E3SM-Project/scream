#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/remap/do_nothing_remapper.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

#include <memory>
#include <numeric>

namespace scream {

MeshFreeGridsManager::
MeshFreeGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p)
 : m_params (p)
 , m_comm   (comm)
{
}

MeshFreeGridsManager::remapper_ptr_type
MeshFreeGridsManager::
do_create_remapper (const grid_ptr_type from_grid,
                    const grid_ptr_type to_grid) const
{
  return std::make_shared<DoNothingRemapper<Real> >(from_grid,to_grid);
}

void MeshFreeGridsManager::
build_grids (const std::set<std::string>& grid_names)
{
  if (grid_names.size()==0) {
    return;
  }

  for (const auto& gn : grid_names) {
    EKAT_REQUIRE_MSG (ekat::contains(supported_grids(),gn),
        "Error! MeshFreeGridsManager only supports 'Point Grid', 'SE Grid DG', and 'SE Grid CG' grids.\n"
        "       Requested grid: " + gn + "\n");
  }

  if (not m_params.isParameter("Reference Grid")) {
    // Set a reference grid.
    if (ekat::contains(grid_names,"Point Grid")) {
      m_params.set<std::string>("Reference Grid","Point Grid");
    } else if (ekat::contains(grid_names,"SE Grid CG")) {
      m_params.set<std::string>("Reference Grid","SE Grid CG");
    } else {
      m_params.set<std::string>("Reference Grid","SE Grid DG");
    }
  }
  std::string ref_grid = m_params.get<std::string>("Reference Grid");

  const bool build_pt    = ekat::contains(grid_names,"Point Grid") || ref_grid=="Point Grid";
  const bool build_se_cg = ekat::contains(grid_names,"SE Grid CG") || ref_grid=="SE Grid CG";
  const bool build_se_dg = ekat::contains(grid_names,"SE Grid DG") || ref_grid=="SE Grid DG";

  const auto& gm_params         = m_params.sublist("Mesh Free");
  const int num_vertical_levels = gm_params.get<int>("Number of Vertical Levels");

  if (build_se_cg) {
    // Assume a "strip" of elements:
    //
    //  *---*---*     *---*
    //  |   |   | ... |   *
    //  *---*---*     *---*
    //
    // partitioned evenly across ranks. Use this config to establish
    // a unique numbering for each element, using the following (np=4 here):
    //
    //   1--5--9-13  13-17-21-25
    //   |  |  |  |   |  |  |  |
    //   2--6-10-14  14-18-22-26
    //   |  |  |  |   |  |  |  |
    //   3--7-11-15  15-19-23-27
    //   |  |  |  |   |  |  |  |
    //   4--8-12-16  16-20-24-28

    const int num_local_elems  = gm_params.get<int>("Number of Local Elements");
    const int num_gp           = gm_params.get<int>("Number of Gauss Points");

    // Set up the degrees of freedom.
    SEGrid::dofs_list_type dofs("", num_local_elems*num_gp*num_gp);
    SEGrid::lid_to_idx_map_type dofs_map("", num_local_elems*num_gp*num_gp, 3);

    auto host_dofs = Kokkos::create_mirror_view(dofs);
    auto host_dofs_map = Kokkos::create_mirror_view(dofs_map);

    // Count unique local dofs. On all elems except the very last one (on rank N),
    // we have num_gp*(num_gp-1) unique dofs;
    int num_elem_unique_dofs = num_gp*(num_gp-1);
    int num_local_unique_dofs = num_local_elems*num_elem_unique_dofs;
    int offset = num_local_unique_dofs*m_comm.rank();

    for (int ie = 0; ie < num_local_elems; ++ie) {
      for (int igp = 0; igp < num_gp-1; ++igp) {
        for (int jgp = 0; jgp < num_gp; ++jgp) {
          int idof = ie*num_gp*num_gp + igp*num_gp + jgp;
          int gid = offset + ie*num_elem_unique_dofs + igp*num_gp + jgp;
          host_dofs(idof) = gid;
          host_dofs_map(idof, 0) = ie;
          host_dofs_map(idof, 1) = igp;
          host_dofs_map(idof, 2) = jgp;
        }
      }
      for (int jgp = 0; jgp < num_gp; ++jgp) {
        int idof = ie*num_gp*num_gp + (num_gp-1)*num_gp + jgp;
        int gid = offset + (ie+1)*num_elem_unique_dofs + jgp;
        host_dofs(idof) = gid;
        host_dofs_map(idof, 0) = ie;
        host_dofs_map(idof, 1) = num_gp-1;
        host_dofs_map(idof, 2) = jgp;
      }
    }

    // Move the data to the device and set the DOFs.
    Kokkos::deep_copy(dofs, host_dofs);
    Kokkos::deep_copy(dofs_map, host_dofs_map);

    // Create the grid, and set the dofs
    auto se_grid = std::make_shared<SEGrid>("SE Grid CG",num_local_elems,num_gp,num_vertical_levels,SEType::CG,m_comm);
    se_grid->setSelfPointer(se_grid);

    se_grid->set_dofs(dofs);
    se_grid->set_lid_to_idx_map(dofs_map);

    m_grids["SE Grid CG"]    = se_grid;
  }
  if (build_se_dg) {
    // Build a set of completely disconnected spectral elements.
    const int num_local_elems  = gm_params.get<int>("Number of Local Elements");
    const int num_gp           = gm_params.get<int>("Number of Gauss Points");

    // Set up the degrees of freedom.
    SEGrid::dofs_list_type dofs("", num_local_elems*num_gp*num_gp);
    SEGrid::lid_to_idx_map_type dofs_map("", num_local_elems*num_gp*num_gp, 3);

    auto host_dofs = Kokkos::create_mirror_view(dofs);
    auto host_dofs_map = Kokkos::create_mirror_view(dofs_map);

    // Count unique local dofs. On all elems except the very last one (on rank N),
    // we have num_gp*(num_gp-1) unique dofs;
    int num_local_dofs = num_local_elems*num_gp*num_gp;
    int offset = num_local_dofs*m_comm.rank();

    for (int ie = 0; ie < num_local_elems; ++ie) {
      for (int igp = 0; igp < num_gp; ++igp) {
        for (int jgp = 0; jgp < num_gp; ++jgp) {
          int idof = ie*num_gp*num_gp + igp*num_gp + jgp;
          int gid = offset + idof;
          host_dofs(idof) = gid;
          host_dofs_map(idof, 0) = ie;
          host_dofs_map(idof, 1) = igp;
          host_dofs_map(idof, 2) = jgp;
        }
      }
    }

    // Move the data to the device and set the DOFs.
    Kokkos::deep_copy(dofs, host_dofs);
    Kokkos::deep_copy(dofs_map, host_dofs_map);

    // Create the grid, and set the dofs
    std::shared_ptr<SEGrid> se_grid;
    se_grid = std::make_shared<SEGrid>("SE Grid DG",num_local_elems,num_gp,num_vertical_levels,SEType::DG,m_comm);
    se_grid->setSelfPointer(se_grid);

    se_grid->set_dofs(dofs);
    se_grid->set_lid_to_idx_map(dofs_map);

    m_grids["SE Grid DG"]    = se_grid;
  }
  if (build_pt) {
    const int num_global_cols  = gm_params.get<int>("Number of Global Columns");
    auto pt_grid = create_point_grid("Point Grid",num_global_cols,num_vertical_levels,m_comm);
    m_grids["Point Grid"] = pt_grid;
  }
}

std::shared_ptr<GridsManager>
create_mesh_free_grids_manager (const ekat::Comm& comm, const int num_local_elems,
                                const int num_gp, const int num_vertical_levels,
                                const int num_global_cols)
{
  ekat::ParameterList gm_params;
  gm_params.sublist("Mesh Free").set("Number of Global Columns",num_global_cols);
  gm_params.sublist("Mesh Free").set("Number of Local Elements",num_local_elems);
  gm_params.sublist("Mesh Free").set("Number of Gauss Points",num_gp);
  gm_params.sublist("Mesh Free").set("Number of Vertical Levels",num_vertical_levels);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  return gm;
}

} // namespace scream
