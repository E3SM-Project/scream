#include "share/grid/abstract_grid.hpp"

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

  m_comm.all_reduce(&m_num_local_dofs,&m_num_global_dofs,1,MPI_SUM);

  // This grid name is also an alias
  m_aliases.push_back(m_name);
}

void AbstractGrid::add_alias (const std::string& alias)
{
  if (not ekat::contains(m_aliases,alias) and alias!=m_name) {
    m_aliases.push_back(alias);
  }
}

bool AbstractGrid::is_unique () const {
  // Get a copy of gids on host. CAREFUL: do not use the stored dofs,
  // since we need to sort dofs in order to call unique, and we don't
  // want to alter the order of gids in this grid.
  auto dofs_h = m_dofs_gids.clone().get_view<gid_type*>();

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

Field
AbstractGrid::get_dofs_gids () const {
  return m_dofs_gids.get_const();
}

Field
AbstractGrid::get_dofs_gids () {
  return m_dofs_gids;
}

Field
AbstractGrid::get_lid_to_idx_map () const {
  return m_lid_to_idx.get_const();
}

Field
AbstractGrid::get_lid_to_idx_map () {
  return m_lid_to_idx;
}

Field
AbstractGrid::get_geometry_data (const std::string& name) const {
  EKAT_REQUIRE_MSG (has_geometry_data(name),
      "Error! Geometry data '" + name + "' not found.\n");

  return m_geo_fields.at(name).get_const();
}

Field
AbstractGrid::get_geometry_data (const std::string& name) {
  EKAT_REQUIRE_MSG (has_geometry_data(name),
      "Error! Geometry data '" + name + "' not found.\n");

  return m_geo_fields.at(name);
}

Field
AbstractGrid::create_geometry_data (const std::string& name, const FieldLayout& layout,
                                    const ekat::units::Units& units,
                                    const DataType data_type)
{
  EKAT_REQUIRE_MSG (not has_geometry_data(name),
      "Error! Cannot create geometry data, since it already exists.\n"
      "  - grid name: " + this->name() + "\n"
      "  - geo data name: " + name + "\n"
      "  - geo data layout: " + to_string(m_geo_fields.at(name).get_header().get_identifier().get_layout()) + "\n"
      "  - input layout: " + to_string(layout) + "\n");

  FieldIdentifier fid(name,layout,units,this->name(),data_type);
  // Create field and the read only copy as well
  auto& f = m_geo_fields[name] = Field(fid);
  f.allocate_view();
  return f;
}

void
AbstractGrid::set_geometry_data (const Field& f)
{
  EKAT_REQUIRE_MSG (not has_geometry_data(f.name()),
      "Error! Cannot set geometry data, since it already exists.\n"
      "  - grid name: " + this->name() + "\n"
      "  - geo data name: " + f.name() + "\n");

  m_geo_fields[f.name()] = f;
}

AbstractGrid::gid_type
AbstractGrid::get_global_min_dof_gid () const
{
  // TODO: we could probably cache these into mutable variables.
  //       But unless we call this method *many* times, it won't matter
  gid_type local_min, global_min;
  auto dofs = get_dofs_gids().get_view<const gid_type*>();
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
  // TODO: we could probably cache these into mutable variables.
  //       But unless we call this method *many* times, it won't matter
  gid_type local_max, global_max;
  auto dofs = get_dofs_gids().get_view<const gid_type*>();
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

std::list<std::string>
AbstractGrid::get_geometry_data_names () const
{
  std::list<std::string> names;
  for (const auto& it : m_geo_fields) {
    names.push_back(it.first);
  }
  return names;
}

void AbstractGrid::reset_num_vertical_lev (const int num_vertical_lev) {
  m_num_vert_levs = num_vertical_lev;

  // TODO: when the PR storing geo data as Field goes in, you should
  //       invalidate all geo data whose FieldLayout contains LEV/ILEV
}

std::vector<AbstractGrid::gid_type>
AbstractGrid::get_unique_gids () const
{
  // Gather local sizes across all ranks
  std::vector<int> ngids (m_comm.size());
  ngids[m_comm.rank()] = get_num_local_dofs();
  m_comm.all_gather(ngids.data(),1);

  std::vector<int> offsets (m_comm.size()+1,0);
  for (int pid=1; pid<=m_comm.size(); ++pid) {
    offsets[pid] = offsets[pid-1] + ngids[pid-1];
  }
  EKAT_REQUIRE_MSG (offsets[m_comm.size()]==m_num_global_dofs,
      "Error! Something went wrong while computing offsets in AbstractGrid::get_unique_grid.\n");

  // Gather all dofs
  const auto mpi_gid_t = ekat::get_mpi_type<gid_type>();
  std::vector<gid_type> all_gids (m_num_global_dofs);
  auto dofs_gids_h = m_dofs_gids.get_view<const gid_type*,Host>();
  MPI_Allgatherv (dofs_gids_h.data(),m_num_local_dofs,mpi_gid_t,
                  all_gids.data(),ngids.data(),offsets.data(),
                  mpi_gid_t,m_comm.mpi_comm());

  // Figure out unique dofs
  std::vector<gid_type> unique_dofs;
  const auto all_gids_beg = all_gids.begin();
  const auto all_gids_end = all_gids.begin() + offsets[m_comm.rank()];
  const auto my_gids_beg = dofs_gids_h.data();
  const auto my_gids_end = dofs_gids_h.data() + m_num_local_dofs;
  for (auto it=my_gids_beg; it!=my_gids_end; ++it) {
    if (std::find(all_gids_beg,all_gids_end,*it)==all_gids_end) {
      unique_dofs.push_back(*it);
    }
  }

  return unique_dofs;
}

std::vector<int> AbstractGrid::
get_owners (const gid_view_h& gids) const
{
  // In order to ship information around across ranks, it is easier to use
  // an auxiliary grid, where dofs are partitioned across ranks linearly
  // NOTE: we actually don't need the grid itself. We only need to know
  //       what the local number of dofs would be on this rank.
  const int ngdofs = get_num_global_dofs();
  const auto& comm = get_comm();
  int nldofs_linear = ngdofs / comm.size();
  if (comm.rank()<(ngdofs % comm.size())) {
    ++ nldofs_linear;
  }

  // For each pid, compute offsets in the global gids array.
  std::vector<int> offsets(comm.size()+1,0);
  const int ndofs_per_rank = ngdofs / comm.size();
  const int remainder = ngdofs % comm.size();
  for (int pid=1; pid<=comm.size(); ++pid) {
    offsets[pid] = offsets[pid-1] + ndofs_per_rank;
    if ( (pid-1)<remainder ){
      ++offsets[pid];
    }
  }
  EKAT_REQUIRE_MSG (offsets.back()==ngdofs,
      "Error! Something went wrong while calling get_gids_owners.\n"
      "  - grid name: " + this->name() + "\n");

  // Utility lambda: given a GID, retrieve the PID that would own it in a
  // linearly distributed grid, as well as the corresponding LID it would
  // have in that grid on that rank. This is doable without any communication
  // since the gids are partitioned linearly
  auto pid_and_lid = [&] (const gid_type gid) -> std::pair<int,int> {
    auto it = std::upper_bound (offsets.begin(),offsets.end(),gid);
    int pid = std::distance(offsets.begin(),it) - 1;
    int lid = gid - offsets[pid];
    EKAT_REQUIRE_MSG (pid>=0 && pid<comm.size(),
        "Error! Failed to retrieve owner of GID in the linear grid.\n");
    return std::make_pair(pid,lid);
  };

  // The idea is to create a "global" array (partitioned across ranks) of the
  // gids in the linear map, use it to store the owner of each dofs (in the original grid),
  // and finally read that global array for all the input gids.
  struct PidLid {
    int pid;
    int lid;
  };
  std::map<int,std::vector<PidLid>> rma_data;
  std::map<int,std::vector<int>> rma_offsets;
  std::map<int,MPI_Datatype> rma_dtypes;

  auto clear_dtypes = [&] () {
    for (auto& it : rma_dtypes) {
      MPI_Type_free(&it.second);
    }
    rma_dtypes.clear();
  };

  MPI_Win win;
  PidLid* pids_lids_linear;
  MPI_Win_allocate (nldofs_linear*sizeof(MPI_2INT),sizeof(MPI_2INT),
                    MPI_INFO_NULL,comm.mpi_comm(),&pids_lids_linear,&win);
  MPI_Win_fence(0,win);

  // Step 1: each rank loops through its grid gids, and sets owners_linear=rank
  // for all its gids in the grid.
  MPI_Win_fence(0,win);

  //  - 1.a Figure where each local dof specs will be written in the linearly
  //        distributed global array
  auto dofs_gids_h = m_dofs_gids.get_view<const gid_type*,Host>();
  for (int i=0; i<get_num_local_dofs(); ++i) {
    const auto gid = dofs_gids_h[i];
    const auto pidlid = pid_and_lid (gid);
    const auto pid = pidlid.first;
    const auto lid = pidlid.second;
    rma_data[pid].push_back({comm.rank(),i});
    rma_offsets[pid].push_back(lid);
  }

  //  - 1.b: build data types for writing all the data on each tgt rank at once
  for (const auto& it : rma_offsets) {
    auto& dtype = rma_dtypes[it.first];
    std::vector<int> ones(it.second.size(),1);
    MPI_Type_indexed (it.second.size(),ones.data(),it.second.data(),MPI_2INT,&dtype);
    MPI_Type_commit (&dtype);
  }

  //  - 1.c: write on the window
  for (const auto& it : rma_data) {
    const auto pid = it.first;
    const auto dtype = rma_dtypes.at(pid);
    // Note: the dtype already encodes offsets in the remote window, so tgt_disp=0
    MPI_Put (it.second.data(),it.second.size(),MPI_2INT,pid,
             0,1,dtype,win);
  }
  clear_dtypes();
  rma_data.clear();
  rma_offsets.clear();
  MPI_Win_fence(0,win);

  // Step 2: each rank loops over its input gids, and retrieves the owner from the window

  //  - 1.a: Figure out what needs to be read from each rank
  for (size_t i=0; i<gids.size(); ++i) {
    const auto gid = gids[i];
    const auto pidlid = pid_and_lid (gid);
    const auto pid = pidlid.first;
    const auto lid = pidlid.second;
    rma_offsets[pid].push_back(lid);
    rma_data[pid].push_back({-1,-1});
  }

  //  - 1.b: build data types for reading all the data from each tgt rank at once
  for (const auto& it : rma_offsets) {
    auto& dtype = rma_dtypes[it.first];
    std::vector<int> ones(it.second.size(),1);
    MPI_Type_indexed (it.second.size(),ones.data(),it.second.data(),MPI_2INT,&dtype);
    MPI_Type_commit (&dtype);
  }

  //  - 1.c: read from the window
  for (auto& it : rma_data) {
    const auto pid = it.first;
    const auto dtype = rma_dtypes.at(pid);
    // Note: the dtype already encodes offsets in the remote window, so tgt_disp=0
    MPI_Get (it.second.data(),it.second.size(),MPI_2INT,pid,
             0,1,dtype,win);
  }
  clear_dtypes();
  rma_offsets.clear();
  MPI_Win_fence(0,win);

  // Step 3: copy data in rma types into output vector, making sure we keep correct order
  std::vector<int> owners(gids.size(),-1);
  std::map<int,int> curr_data_index;
  for (size_t i=0; i<gids.size(); ++i) {
    const auto gid = gids[i];
    const auto pidlid = pid_and_lid (gid);
    const auto pid = pidlid.first;
    auto it_bool = curr_data_index.emplace(pid,0);
    auto& idx = it_bool.first->second;
    owners[i] = rma_data[pid][idx].pid;
    ++idx;
  }
  rma_data.clear();

  // Clean up
  MPI_Win_free(&win);

  return owners;
}

void AbstractGrid::create_dof_fields (const int scalar2d_layout_rank)
{
  using namespace ShortFieldTagsNames;
  const auto units = ekat::units::Units::nondimensional();

  // The dof gids field is a 1d field, while lid2idx has rank 2.
  // For both, the 1st dim is the num of local dofs. The 2nd dime of
  // lid2idx is the rank of a 2d scalar layout.
  FieldLayout dof_layout({COL},{get_num_local_dofs()});
  FieldLayout lid2idx_layout({COL,CMP},{get_num_local_dofs(),scalar2d_layout_rank});
  m_dofs_gids = Field(FieldIdentifier("gids",dof_layout,units,m_name,DataType::IntType));
  m_lid_to_idx = Field(FieldIdentifier("lid2idx",lid2idx_layout,units,m_name,DataType::IntType));

  m_dofs_gids.allocate_view();
  m_lid_to_idx.allocate_view();
}

void AbstractGrid::copy_data (const AbstractGrid& src, const bool shallow)
{
  if (shallow) {
    m_dofs_gids = src.m_dofs_gids;
  } else {
    m_dofs_gids = src.m_dofs_gids.clone();
  }

  if (shallow) {
    m_lid_to_idx = src.m_lid_to_idx;
  } else {
    m_lid_to_idx = src.m_lid_to_idx.clone();
  }

  for (const auto& name : src.get_geometry_data_names()) {
    if (shallow) {
      m_geo_fields[name] = src.m_geo_fields.at(name);
    } else {
      m_geo_fields[name] = src.m_geo_fields.at(name).clone();
    }
  }
}

} // namespace scream
