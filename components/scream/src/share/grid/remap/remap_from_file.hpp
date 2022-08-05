#ifndef SCREAM_REMAPPER_FROM_FILE_HPP
#define SCREAM_REMAPPER_FROM_FILE_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/scream_types.hpp"

#include "ekat/mpi/ekat_comm.hpp"

namespace scream {

/*
 * A utility structure that stores information about remapping given
 * a source grid and a target grid.
 * 
 * This structure assumes that the remapping is only in the horizontal.
 * This structure assumes that the source and target grid are both defined
 * in the 2D as a single vector of columns.  Thus the remapped value can be
 * represented by the expression:
 *     Y(target_column_j) = sum_(n=0)^(N_j-1) [ w_j(n) * X(source_column_j(n)) ]
 *  where,
 *     target_column_j: is the index in the simulation grid for the j'th DOF.
 *     source_column_j: is the set of column indices in the source data used to map to the j'th DOF.
 *     w_j:             is the set of weights used in mapping to the j'th DOF.
 *     X:               is the source data over all columns.
 *     Y:               is the mapped data on the target mesh for all DOF's.
 *     N_j:             is the number of source columns used to map to the j'th DOF
 */

struct GSMap { // TODO: This can be a different name, using GSMap like CPL for now

  using gid_type = AbstractGrid::gid_type;

  using KT = KokkosTypes<DefaultDevice>;
  template <typename S>
  using view_1d = typename KT::template view_1d<S>;

public:
  virtual ~GSMap() = default;
  GSMap();
  GSMap(const std::string& remap_file_name_,const view_1d<gid_type>& dofs_gids_,const gid_type min_dof_);

  std::string name = "";
  std::string remap_file_name;  // Name of the file where the remap information is stored.

  view_1d<gid_type> unique_source_dof; // A list of unique source dofs, determined from source_grid_dof.  This is used to gather just the source data needed.
  view_1d<gid_type> source_grid_dof;  // Degree of freedom on source mesh
  view_1d<gid_type> target_grid_dof;  // Degree of freedom on target mesh
  view_1d<Real>     map_weight;       // Weight associated with map from source to target
  view_1d<gid_type> dofs_gids;        // Global ID's for dofs on target mesh, defined by the grid.
  // DOF offset for meshes.  For example if dof list is 1-based while the dof indices are 0-based.
  gid_type              source_dof_offset, target_dof_offset;

  void set_remap_indices(
    const ekat::Comm&        comm,
    const std::string&       remap_file_name,
          std::vector<Int>&  seg_dof,
          std::vector<Int>&  seg_start,
          std::vector<Int>&  seg_length);

protected:

  void read_remap_data_nco(
    const ekat::Comm& comm,
    const std::string& remap_file_name
  );

  void set_data_chunksize (
    const ekat::Comm& comm,
    const Int         total_length,
          Int         my_chunksize,
          Int         my_start);

  void construct_data_segments (
    const Int                       my_start,
    const view_1d<Int>::HostMirror& target_dofs,
    std::vector<Int>&               seg_dof,
    std::vector<Int>&               seg_start,
    std::vector<Int>&               seg_length);

  void consolidate_data_segments_on_this_rank(
    const ekat::Comm&        comm,
    std::vector<Int>&        seg_dof,
    std::vector<Int>&        seg_start,
    std::vector<Int>&        seg_length);

  void set_map(
    const std::vector<Int>&               seg_dof,
    const std::vector<Int>&               seg_start,
    const std::vector<Int>&               seg_length);
}; // end struct GSMap  

} // namespace scream

#endif // SCREAM_REMAPPER_FROM_FILE_HPP
