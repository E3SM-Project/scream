#include "share/grid//abstract_grid.hpp"
#include <Kokkos_CopyViews.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <algorithm>
#include <cstring>
#include <string>

namespace scream
{
// Constructor(s) & Destructor
AbstractGrid::
AbstractGrid (const std::string& name,
              const GridType type,
              const int num_local_dofs,
              const int num_vertical_lev,
              const ekat::Comm& comm)
 : m_type (type)
 , m_name (name)
 , m_num_local_dofs (num_local_dofs)
 , m_num_vert_levs  (num_vertical_lev)
 , m_comm (comm)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (m_num_local_dofs>=0, "Error! Number of local dofs must be non-negative.\n");
  EKAT_REQUIRE_MSG (m_num_vert_levs>=2, "Error! Number of vertical levels must be at least 2.\n");

  m_comm.all_reduce(&m_num_local_dofs,&m_num_global_dofs,1,MPI_SUM);
}

bool AbstractGrid::is_unique () const {
  EKAT_REQUIRE_MSG (m_dofs_set,
      "Error! We cannot establish if this grid is unique before the dofs GIDs are set.\n");

  // Get a copy of gids on host. CAREFUL: do not use create_mirror_view,
  // since it would create a shallow copy on CPU devices, but we need a
  // deep copy, to prevent altering order of gids.
  decltype(m_dofs_gids)::HostMirror dofs_h("",m_dofs_gids.size());
  Kokkos::deep_copy(dofs_h,m_dofs_gids);

  std::sort(dofs_h.data(),dofs_h.data()+m_num_local_dofs);
  auto unique_end = std::unique(dofs_h.data(),dofs_h.data()+m_num_local_dofs);

  int locally_unique = unique_end==(dofs_h.data()+m_num_local_dofs);
  int unique;
  m_comm.all_reduce(&locally_unique,&unique,1,MPI_PROD);
  if (unique==0) {
    return false;
  }

  // Each rank has unique gids locally. Now it's time to verify if they are also globally unique.
  int max_dofs;
  m_comm.all_reduce(&m_num_local_dofs,&max_dofs,1,MPI_MAX);
  std::vector<int> gids(max_dofs);
  int unique_gids = 1;

  for (int pid=0; pid<m_comm.size(); ++pid) {
    // Rank pid broadcasts its gids, everyone else checks if there are duplicates
    if (pid==m_comm.rank()) {
      auto start = dofs_h.data();
      auto end   = start + m_num_local_dofs;
      std::copy(start,end,gids.data());
    }

    int ndofs = m_num_local_dofs;
    m_comm.broadcast(&ndofs,1,pid);
    m_comm.broadcast(gids.data(),ndofs,pid);

    int my_unique_gids = 1;
    if (pid!=m_comm.rank()) {
      // Checking two sorted arrays of length m and n for elements in common is O(m+n) ops.
      int i=0, j=0;
      while (i<m_num_local_dofs && j<ndofs && my_unique_gids==1) {
        if (dofs_h[i]<gids[j]) {
          ++i;
        } else if (dofs_h[i]>gids[j]) {
          ++j;
        } else {
          // Found a match. We can stop here
          my_unique_gids = 0;
          break;
        }
      }
    }
    m_comm.all_reduce(&my_unique_gids,&unique_gids,1,MPI_PROD);
    if (unique_gids==0) {
      break;
    }
  }

  return unique_gids;
}

void AbstractGrid::
set_dofs (const dofs_list_type& dofs)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_dofs_set, "Error! Dofs cannot be re-set, once set.\n");

  EKAT_REQUIRE_MSG (dofs.extent_int(0)==m_num_local_dofs,
      "Error! Wrong size for the input dofs list. It should match the number of local dofs stored in the mesh.\n"
      "       Expected gids list size: " + std::to_string(m_num_local_dofs) + "\n"
      "       Input gids list size   : " + std::to_string(dofs.size()) + "\n");

  m_dofs_gids  = dofs;
  m_dofs_set = true;

#ifndef NDEBUG
  this->check_dofs_list();
#endif
}

void AbstractGrid::
set_lid_to_idx_map (const lid_to_idx_map_type& lid_to_idx)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_lid_to_idx_set, "Error! The lid->idx map cannot be re-set, once set.\n");

  EKAT_REQUIRE_MSG (lid_to_idx.extent_int(0)==m_num_local_dofs &&
                    lid_to_idx.extent_int(1)==get_2d_scalar_layout().rank(),
      "Error! Wrong size(s) for the input lid_to_idx map. They should match (num_local_dofs, 2d layout rank).\n"
      "       Expected sizes: (" + std::to_string(m_num_local_dofs) + "," + std::to_string(get_2d_scalar_layout().rank()) + ")\n"
      "       Input map sizes: (" + std::to_string(lid_to_idx.extent(0)) + "," + std::to_string(lid_to_idx.extent(1)) + ")\n");

  m_lid_to_idx = lid_to_idx;
  m_lid_to_idx_set = true;

#ifndef NDEBUG
  this->check_lid_to_idx_map();
#endif
}

void AbstractGrid::
set_geometry_data (const std::string& name, const geo_view_type& data) {
  // Sanity checks
  EKAT_REQUIRE_MSG (data.extent_int(0)==m_num_local_dofs,
                    "Error! Input geometry data has wrong dimensions.\n");

#ifndef NDEBUG
  this->check_geo_data(name,data);
#endif

  m_geo_views[name] = data;
}

AbstractGrid::gid_type
AbstractGrid::get_global_min_dof_gid () const
{
  EKAT_REQUIRE_MSG (m_dofs_set,
      "Error! You need to set dofs gids before you can compute the global min dof.\n");
  // TODO: we could probably cache these into mutable variables.
  //       But unless we call this method *many* times, it won't matter
  gid_type local_min, global_min;
  auto dofs = get_dofs_gids();
  Kokkos::parallel_reduce(Kokkos::RangePolicy<>(0,get_num_local_dofs()),
      KOKKOS_LAMBDA (const int& i, gid_type& lmin) {
        if (dofs(i) < lmin) {
          lmin = dofs(i);
        }
      },Kokkos::Min<gid_type>(local_min));
  Kokkos::fence();

  m_comm.all_reduce(&local_min,&global_min,1,MPI_MIN);

  return global_min;
}

AbstractGrid::gid_type
AbstractGrid::get_global_max_dof_gid () const
{
  EKAT_REQUIRE_MSG (m_dofs_set,
      "Error! You need to set dofs gids before you can compute the global max dof.\n");
  // TODO: we could probably cache these into mutable variables.
  //       But unless we call this method *many* times, it won't matter
  gid_type local_max, global_max;
  auto dofs = get_dofs_gids();
  Kokkos::parallel_reduce(Kokkos::RangePolicy<>(0,get_num_local_dofs()),
      KOKKOS_LAMBDA (const int& i, gid_type& lmax) {
        if (dofs(i) > lmax) {
          lmax = dofs(i);
        }
      },Kokkos::Max<gid_type>(local_max));
  Kokkos::fence();

  m_comm.all_reduce(&local_max,&global_max,1,MPI_MAX);

  return global_max;
}

} // namespace scream
