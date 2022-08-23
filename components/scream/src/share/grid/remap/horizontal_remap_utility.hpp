#ifndef EAMXX_HORIZONTAL_REMAP_UTILITY_HPP
#define EAMXX_HORIZONTAL_REMAP_UTILITY_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/scream_types.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {

using gid_type = AbstractGrid::gid_type;

using KT = KokkosTypes<DefaultDevice>;
template <typename S>
using view_1d = typename KT::template view_1d<S>;

template <typename S>
using view_1d_host = typename KT::template view_1d<S>::HostMirror;


/* --------------------------------------------------------------------------------------------- */
/* A simple structure to store remap information for a single target column
 * --------------------------------------------------------------------------------------------- */
struct RemapSegment {

public:
  virtual ~RemapSegment() = default;
  RemapSegment() {};
  RemapSegment(const gid_type dof_, const Int length_) :
    m_dof(dof_),
    m_length(length_)
  {
    source_dofs = view_1d<gid_type>("",m_length);
    source_idx  = view_1d<int>("",m_length);
    weights     = view_1d<Real>("",m_length);
  }

  gid_type          m_dof;       // The degree of freedom (dof) on target mesh this segment maps to.
  view_1d<gid_type> source_dofs; // A view of the dofs on the source grid that map to this target dof.
  view_1d<int>      source_idx;  // Local index in the set of unique columns for source data - set by GSMap
  view_1d<Real>     weights;     // A view of the associated weights in the mapping.
  Int               m_length;    // The size of this segment.

  // Helper function to check a segment is valid.
  void check() {
    EKAT_REQUIRE(source_dofs.extent(0)==m_length);
    EKAT_REQUIRE(weights.extent(0)==m_length);
    EKAT_REQUIRE(source_idx.extent(0)==m_length);
    Real wgt = 0.0;
    Kokkos::parallel_reduce("", m_length, KOKKOS_LAMBDA (const int& ii, Real& lsum) {
      lsum += weights(ii);
    },wgt);
    // The sum of the weights should be = 1.0, note due to round-off error we accept to within a tolerance.
    Real tol = std::numeric_limits<Real>::epsilon() * 10.0;
    EKAT_REQUIRE_MSG(std::abs(wgt-1.0)<tol,"ERROR: checking remap segment for DOF = " + std::to_string(m_dof) + ", total weight = " + std::to_string(wgt) + ".");
  }

  // Function to apply segment to data
  Real apply_segment(
    const view_1d<Real>& source_data)
  {
    Real ret = 0;
    auto source_data_h = Kokkos::create_mirror_view(source_data);
    auto source_idx_h  = Kokkos::create_mirror_view(source_idx);
    auto weights_h     = Kokkos::create_mirror_view(weights);
    Kokkos::deep_copy(source_data_h,source_data);
    Kokkos::deep_copy(source_idx_h,source_idx);
    Kokkos::deep_copy(weights_h   ,weights);
    for (int ii=0; ii< m_length; ii++) {
      int idx = source_idx_h(ii);
      ret += source_data_h(idx)*weights_h(ii);
    }
    return ret;
  }

}; // end struct RemapSegment

/* --------------------------------------------------------------------------------------------- */
/* A structure to store all of the map information 
 * --------------------------------------------------------------------------------------------- */
struct GSMap { // TODO: Following the name in component coupler, could be a different name

public:
  virtual ~GSMap() = default;
  GSMap() {};

/*---------------------------------------------*/
  Int get_num_of_segs() {return map_segments.size();}
  Int get_num_of_dofs() {return m_dofs_gids.size();}
/*---------------------------------------------*/
  RemapSegment get_segment(const gid_type dof) {
    return get_segment_impl(dof);
  }
/*---------------------------------------------*/
  void set_dofs_gids(const view_1d<gid_type>& dofs_gids) {
    REQUIRE(dofs_gids.size()>0);
    m_dofs_gids = view_1d<gid_type>(dofs_gids);
    Kokkos::deep_copy(m_dofs_gids,dofs_gids);
  }
/*---------------------------------------------*/
  std::vector<gid_type> get_unique_dofs() { return unique_dofs; };
/*---------------------------------------------*/
  void add_segment(RemapSegment& segment) {
    add_segment_impl(segment);
  }; 
/*---------------------------------------------*/
  void check() {
    std::vector<bool> found(m_dofs_gids.size(),false);
    auto dofs_gids_h = Kokkos::create_mirror_view(m_dofs_gids);
    Kokkos::deep_copy(dofs_gids_h,m_dofs_gids);
    for (int ii=0; ii<map_segments.size(); ii++) {
      auto& seg = map_segments[ii];
      seg.check();
      for (int jj=0; jj<dofs_gids_h.size(); jj++) {
        if (dofs_gids_h(jj) - m_min_dof == seg.m_dof - source_min_dof) {
          found[jj] = true;
          break;
        }
      }
    }
    REQUIRE(std::find(found.begin(),found.end(),false) == found.end());
    EKAT_REQUIRE_MSG(m_min_dof != -999, "Error in GSMap: m_min_dof has not been set.");
    EKAT_REQUIRE_MSG(source_min_dof != -999, "Error in GSMap: source_min_dof has not been set.");
  };
/*---------------------------------------------*/
  // Useful if needing to grab source data from somewhere, like a data file.
  void set_unique_source_dofs() {
    for (int iseg=0; iseg<get_num_of_segs(); iseg++) {
      const auto& seg = map_segments[iseg];
      const auto& src_dofs = seg.source_dofs;
      const auto& src_dofs_h = Kokkos::create_mirror_view(src_dofs);
      Kokkos::deep_copy(src_dofs_h,src_dofs);
      for (int ii=0; ii<seg.m_length; ii++) {
        auto idx = std::find(unique_dofs.begin(), unique_dofs.end(), src_dofs_h(ii));
        if (idx == unique_dofs.end()) {
          unique_dofs.push_back(src_dofs_h(ii));
        }
      }
    }
    std::sort(unique_dofs.begin(), unique_dofs.end());
    for (int iseg=0; iseg<get_num_of_segs(); iseg++) {
      const auto& seg = map_segments[iseg];
      const auto& src_dofs = seg.source_dofs;
      const auto& src_dofs_h = Kokkos::create_mirror_view(src_dofs);
      auto& src_idx  = seg.source_idx;
      auto  src_idx_h  = Kokkos::create_mirror_view(src_idx);
      for (int ii=0; ii<seg.m_length; ii++) {
        auto idx = std::find(unique_dofs.begin(), unique_dofs.end(), src_dofs_h(ii));
        int idx_ii = idx - unique_dofs.begin();
        src_idx_h(ii) = idx - unique_dofs.begin();
      }
      Kokkos::deep_copy(src_idx,src_idx_h);
    }
  }
/*---------------------------------------------*/
  RemapSegment get_segment_impl(const gid_type dof) {
    RemapSegment ret;
    bool found = false;
    for (int ii=0; ii<map_segments.size(); ii++) {
      if (map_segments[ii].m_dof == dof) {
        found = true;
        ret = map_segments[ii];
        break;
      }
    }
    EKAT_REQUIRE_MSG(found,"Error: GSMap - get_segment, segment for  dof " + std::to_string(dof) + " not found.");
    return ret;
  };
/*---------------------------------------------*/
  void set_segments_from_file(
    const std::string&       remap_file,
    const ekat::Comm&        comm,
    const view_1d<gid_type>& dofs_gids,
    const Int                min_dof)
  {
    // Open remap file and determine the amount of data to be read
    scorpio::register_file(remap_file,scorpio::Read);
    const auto remap_size = scorpio::get_dimlen_c2f(remap_file.c_str(),"n_s"); // Note, here we assume a standard format of col, row, S
    // Distribute responsibility for reading remap data over all ranks
    const int my_rank   = comm.rank();
    const int num_ranks = comm.size();
    // my_chunk will represent the chunk of data this rank will read from file.
    int my_chunk        = remap_size/num_ranks;
    int remainder       = remap_size - (my_chunk*num_ranks);
    if (remainder != 0) {
      my_chunk += my_rank<remainder ? 1 : 0;
    }
    // now determine where this rank start reading the data.
    int* chunks_glob = new int[num_ranks];
    comm.all_gather(&my_chunk,chunks_glob,1);
    int my_start = 0;
    for (int ii=0; ii<my_rank;ii++) {
      my_start += chunks_glob[ii];
    }
    // Check that the total set of chunks covers all the data
    {
      int chunk_check = 0;
      for (int ii=0; ii<num_ranks; ii++) {
        chunk_check += chunks_glob[ii];
      }
      EKAT_REQUIRE_MSG(chunk_check==remap_size,"ERROR: get_remap_indices - Something went wrong distributing remap data among the MPI ranks");
    }
    // Using scream input routines, read remap data from file by chunk
    view_1d<Int>  tgt_col("row",my_chunk); 
    auto tgt_col_h = Kokkos::create_mirror_view(tgt_col);
    std::vector<std::string> vec_of_dims = {"n_s"};
    std::string i_decomp = "Int-n_s";
    scorpio::get_variable(remap_file, "row", "row", vec_of_dims.size(), vec_of_dims, PIO_INT, i_decomp);
    std::vector<int> var_dof(my_chunk);
    std::iota(var_dof.begin(),var_dof.end(),my_start);
    scorpio::set_dof(remap_file,"row",var_dof.size(),var_dof.data());
    scorpio::set_decomp(remap_file);
    scorpio::grid_read_data_array(remap_file,"row",0,tgt_col_h.data()); 
    scorpio::eam_pio_closefile(remap_file);
    // Organize data into sets of target column, start location in data and length of data.
    // At the same time, determine the min_dof for remap column indices.
    std::vector<Int> chunk_dof, chunk_start, chunk_len;
    chunk_dof.push_back(tgt_col_h(0));
    chunk_start.push_back(my_start);
    chunk_len.push_back(1);
    int remap_min_dof = tgt_col_h(0);
    for (int ii=1; ii<my_chunk; ii++) {
      remap_min_dof = std::min(tgt_col_h(ii),remap_min_dof);
      if (tgt_col_h(ii) == chunk_dof.back()) {
        // Then we add one to the length for this chunk.
        chunk_len.back() ++;
      } else {
        // Start a new chunk for a new DOF
        chunk_dof.push_back(tgt_col_h(ii));
        chunk_start.push_back(my_start+ii);
        chunk_len.push_back(1);
      }
    }
    // Pass chunk information among all ranks so they can be consolidated.
    int  num_chunks          = chunk_dof.size();
    int* num_chunks_per_rank = new int[num_ranks];
    int* chunk_displacement  = new int[num_ranks];
    int  total_num_chunks;
    int  global_remap_min_dof;
    comm.all_gather(&num_chunks, num_chunks_per_rank,1);
    comm.all_reduce(&remap_min_dof,&global_remap_min_dof,1,MPI_MIN);
    chunk_displacement[0] = 0;
    total_num_chunks = num_chunks_per_rank[0];
    for (int ii=1; ii<num_ranks; ii++) {
      chunk_displacement[ii] = total_num_chunks;
      total_num_chunks += num_chunks_per_rank[ii];
    }
    int* buff_dof = (int*)calloc(total_num_chunks, sizeof(int));
    int* buff_sta = (int*)calloc(total_num_chunks, sizeof(int));
    int* buff_len = (int*)calloc(total_num_chunks, sizeof(int));
    MPI_Allgatherv(chunk_dof.data(),  chunk_dof.size(),MPI_INT,buff_dof,num_chunks_per_rank,chunk_displacement,MPI_INT,comm.mpi_comm());
    MPI_Allgatherv(chunk_start.data(),chunk_dof.size(),MPI_INT,buff_sta,num_chunks_per_rank,chunk_displacement,MPI_INT,comm.mpi_comm());
    MPI_Allgatherv(chunk_len.data(),  chunk_dof.size(),MPI_INT,buff_len,num_chunks_per_rank,chunk_displacement,MPI_INT,comm.mpi_comm());
    // Now construct and add segments for just the DOF's this rank cares about.
    std::vector<int> seg_dof, seg_start, seg_length;
    var_dof.clear();
    auto dofs_gids_h = Kokkos::create_mirror_view(dofs_gids);
    Kokkos::deep_copy(dofs_gids_h,dofs_gids);
    for (int ii=0; ii<total_num_chunks; ii++) {
      // Search dofs to see if this chunk matches dofs on this rank
      for (int jj=0; jj<dofs_gids_h.extent(0); jj++) {
        if (buff_dof[ii]-global_remap_min_dof == dofs_gids_h(jj)-min_dof) {
          std::vector<int> var_tmp(buff_len[ii]);
          std::iota(var_tmp.begin(),var_tmp.end(),buff_sta[ii]);
          seg_dof.push_back(buff_dof[ii]);
          seg_start.push_back(var_dof.size()); 
          seg_length.push_back(buff_len[ii]);
          var_dof.insert(var_dof.end(),var_tmp.begin(),var_tmp.end());
        } 
      }
    }
    // Now that we know which parts of the remap file this rank cares about we can construct segments
    view_1d<Int>  col("col",var_dof.size()); 
    view_1d<Real> S("S",var_dof.size()); 
    auto col_h = Kokkos::create_mirror_view(col);
    auto S_h = Kokkos::create_mirror_view(S);
    vec_of_dims = {"n_s"};
    i_decomp = "Int-n_s";
    std::string r_decomp = "Real-n_s";
    scorpio::register_file(remap_file,scorpio::Read);
    scorpio::get_variable(remap_file, "col", "col", vec_of_dims.size(), vec_of_dims, PIO_INT, i_decomp);
    scorpio::get_variable(remap_file, "S", "S", vec_of_dims.size(), vec_of_dims, PIO_REAL, r_decomp);
    scorpio::set_dof(remap_file,"col",var_dof.size(),var_dof.data());
    scorpio::set_dof(remap_file,"S",var_dof.size(),var_dof.data());
    scorpio::set_decomp(remap_file);
    scorpio::grid_read_data_array(remap_file,"col",0,col_h.data()); 
    scorpio::grid_read_data_array(remap_file,"S",0,S_h.data()); 
    scorpio::eam_pio_closefile(remap_file);
    Kokkos::deep_copy(col,col_h);
    Kokkos::deep_copy(S,S_h);
    // Construct segments based on data just read from file
    for (int ii=0; ii<seg_dof.size(); ii++) {
      RemapSegment seg(seg_dof[ii],seg_length[ii]);
      int seglength = seg_length[ii];
      int segstart  = seg_start[ii];
      Kokkos::parallel_for("", seglength, KOKKOS_LAMBDA (const int& jj) {
        int idx = segstart + jj;
        seg.source_dofs(jj) = col(idx);
        seg.weights(jj)     = S(idx);
      });
      add_segment(seg);
    }
    m_min_dof = min_dof;
    source_min_dof = global_remap_min_dof;
  } // end set_segments_from_file
/*---------------------------------------------*/
  void add_segment_impl(RemapSegment& segment) {
    // First see if DOF is already in map segments, if so, expand that segment, otherwise just add this one.
    bool seg_found = false;
    RemapSegment seg_match;
    Int          match_loc;
    // Search segments
    for (int ii=0;ii<map_segments.size();ii++) {
      seg_match = map_segments[ii];
      if (seg_match.m_dof == segment.m_dof) {
        // Segment matches
        seg_found = true;
        match_loc = ii;
        break;
      }
    }

    if (seg_found) {
      // Adjust the found segment to include this new segment information.
      Int new_length = seg_match.m_length + segment.m_length;
      RemapSegment new_seg(seg_match.m_dof, new_length);
      auto source_dofs_h = Kokkos::create_mirror_view(new_seg.source_dofs);
      auto weights_h     = Kokkos::create_mirror_view(new_seg.weights);
      Int nsize = seg_match.m_length; 
      Kokkos::parallel_for("", nsize, KOKKOS_LAMBDA (const int& ii) {
        new_seg.source_dofs(ii) = seg_match.source_dofs(ii);
        new_seg.weights(ii)     = seg_match.weights(ii);
      });
      Kokkos::fence();
      nsize = segment.m_length;
      Kokkos::parallel_for("", nsize, KOKKOS_LAMBDA (const int& ii) {
        new_seg.source_dofs(seg_match.m_length+ii) = segment.source_dofs(ii);
        new_seg.weights(seg_match.m_length+ii)     = segment.weights(ii);
      });
      Kokkos::fence();
      map_segments[match_loc] = new_seg;
    } else {
      map_segments.push_back(segment);
    }
  }; // end add_segment_impl
/*---------------------------------------------*/
  void apply_remap(
    const view_1d<Real>& source_data,
          view_1d<Real>& remapped_data)
  {
    auto remap_data_h  = Kokkos::create_mirror_view(remapped_data);
    for (int iseg=0; iseg<get_num_of_segs(); iseg++) {
      auto seg = map_segments[iseg];
      Real remap_val = seg.apply_segment(source_data);
      remap_data_h(iseg) = remap_val;
    }
    Kokkos::deep_copy(remapped_data,remap_data_h);
  }; // end apply_remap
/*---------------------------------------------*/
  gid_type                  m_min_dof = -999;      // The global minimum dof value - needed to offset for array index of data.
  gid_type                  source_min_dof = -999; // Similar to m_min_dof, but the value used in the source dofs.

protected:

/*---------------------------------------------*/
  std::vector<RemapSegment> map_segments;
  std::vector<gid_type>     unique_dofs;
  view_1d<gid_type>         m_dofs_gids;

}; // end struct GSMap

} //namespace scream

#endif // EAMXX_HORIZONTAL_REMAP_UTILITY_HPP
