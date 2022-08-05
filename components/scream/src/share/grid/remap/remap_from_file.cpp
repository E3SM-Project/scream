
#include "share/grid/remap/remap_from_file.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"
#include "ekat/ekat_workspace.hpp"
#include <numeric>

namespace scream
{

/* -----------------------------------------------------------------------
 * Constructor
 * ----------------------------------------------------------------------*/
  GSMap::GSMap(const std::string& remap_file_name_,const view_1d<gid_type>& dofs_gids_,const gid_type min_dof_) : 
    remap_file_name(remap_file_name_),
    dofs_gids(dofs_gids_),
    target_dof_offset(min_dof_)
  {
    // Nothing else to do.
  }

  GSMap::GSMap()
  {
    // Do Nothing
  }
/* -----------------------------------------------------------------------
 * This function will break the remap file data into a set of segments. 
 * Each segment will represent a single degree-of-freedom on the source 
 * mesh, which will start at location n in the contiguous 1-D array 
 * of mapping data and have a certain length. 
 * This approach follows the example used by the MCT coupler to handle 
 * large inter-component maps. 
 * ASSUMPTIONS: 
 * The remap data is cast as a 1-D array of source-to-target indices. 
 * The remap variable "col" represents the degree-of-freedom index for 
 * the source data. 
 * The remap variable "row" represents the degree-of-freedom index for 
 * the target data, i.e. the dof on the simulation grid. 
 * The remap variable "S" represents the weight associated with this 
 * map from source to target. 
 * All three of these variables have length "n_s", defined in the file. 
 * Not required, but it is assumed that the source indices are organized 
 * in such a way that the remapping wieghts for a specific source dof is 
 * contiguous in the 1-D array. 
 * GOAL: 
 * The goal of this function is to have each MPI rank process a chunk 
 * of the remap data and organize it. 
 * Only the "row" variable is read, as we are focused on which segments 
 * of the remap data correspond to which degrees-of-freedom on the target 
 * grid. 
 * The three outputs are vectors of integers which define a specific 
 * segment. 
 * seg_dof: The source degree-of-freedom this segment maps to. 
 * seg_start: The starting element in the 1-D array of remap data for 
 * this segment. 
 * seg_length: The number of remapping weights in this segment. 
 * ------------------------------- 
 * A.S. Donahue (LLNL) July 29, 2022 
 * ----------------------------------------------------------------------*/
void GSMap::
set_remap_indices (
    const ekat::Comm&        comm,
    const std::string&       remap_file_name,
          std::vector<Int>&  seg_dof,
          std::vector<Int>&  seg_start,
          std::vector<Int>&  seg_length)
{
  // At this point the input vectors should all be empty.  We enforce
  // this by clearing them
  seg_dof.clear();
  seg_start.clear();
  seg_length.clear();

  //
}
/* -----------------------------------------------------------------------
 * Helper function to read in a small chunk of the remap data per rank.
 * This function first determines the total size of the data, then sets a
 * fraction of that to be read by this rank, given the overall size of the
 * communicator group.
 * Once the size is known, each rank reads a set of data and organizes the
 * data into a set of segments, each segment has the three properties:
 *   seg_dof: represents the target dof this segment refers to.
 *   seg_str: is the starting index in the global data for this segment
 *   seg_len: is the number of values in this segment.
 * -------------------------------
 *  the "nco" suffix sets assumptions that the remap file has been defined for
 *  use by the nco tool ncremap.  
 * ------------------------------- 
 * A.S. Donahue (LLNL) July 29, 2022 
 * ----------------------------------------------------------------------*/
void GSMap::
read_remap_data_nco (
  const ekat::Comm& comm,
  const std::string& remap_file_name
  )
{
  // NCO based remaps are 1-based in the DOF - so set the appropriate offset
  source_dof_offset = 1;

  // Determine the size of the data to be read on this rank
  scorpio::register_file(remap_file_name,scorpio::Read);
  const Int total_length = scorpio::get_dimlen_c2f(remap_file_name.c_str(),"n_s");
  Int my_chunk_size = 0;
  Int my_start = 0;
  set_data_chunksize(comm,total_length,my_chunk_size,my_start);

  // Read remap data from file
  view_1d<Int>  row("row",my_chunk_size); 
  auto target_dofs_h = Kokkos::create_mirror_view(row);
  std::vector<std::string> vec_of_dims = {"n_s"};
  std::string i_decomp = "Int-n_s";
  scorpio::get_variable(remap_file_name, "row", "row", vec_of_dims.size(), vec_of_dims, PIO_INT, i_decomp);
  std::vector<int> var_dof(my_chunk_size);
  std::iota(var_dof.begin(),var_dof.end(),my_start);
  scorpio::set_dof(remap_file_name,"row",var_dof.size(),var_dof.data());
  scorpio::set_decomp(remap_file_name);
  scorpio::grid_read_data_array(remap_file_name,"row",0,target_dofs_h.data()); 
  scorpio::eam_pio_closefile(remap_file_name);

  // Construct segments from data
  std::vector<Int> seg_dof, seg_start, seg_length; 
  construct_data_segments(my_start, target_dofs_h, seg_dof, seg_start, seg_length);

  // Consolidate segments based on which DOFs are on this rank.  
  consolidate_data_segments_on_this_rank(comm,seg_dof,seg_start,seg_length);
}

/* -----------------------------------------------------------------------*/
void GSMap::
consolidate_data_segments_on_this_rank(
  const ekat::Comm&        comm,
  std::vector<Int>&        seg_dof,
  std::vector<Int>&        seg_start,
  std::vector<Int>&        seg_length
  )
{ 
  // Gather from all ranks the total number of segments and their
  // displacements for an allgatherv
  const int num_of_ranks    = comm.size();
  int  num_of_segs          = seg_dof.size(); // count for passing
  int* num_of_segs_per_rank = new int[num_of_ranks]; // counts for root to recieve
  int* seg_disp             = new int[num_of_ranks]; // displacements
  int  total_num_segs       = 0;
  comm.all_gather(&num_of_segs,num_of_segs_per_rank,1);
  seg_disp[0] = 0;
  total_num_segs = num_of_segs_per_rank[0];
  for (int ii=1;ii<num_of_ranks;ii++) {
    seg_disp[ii]    = seg_disp[ii-1] + num_of_segs_per_rank[ii-1];
    total_num_segs += num_of_segs_per_rank[ii];
  }
  int* buff_dof = (int*)calloc(total_num_segs, sizeof(int));
  int* buff_srt = (int*)calloc(total_num_segs, sizeof(int));
  int* buff_len = (int*)calloc(total_num_segs, sizeof(int));
  MPI_Allgatherv(seg_dof.data(),total_num_segs,MPI_INT,buff_dof,num_of_segs_per_rank,seg_disp,MPI_INT,comm.mpi_comm());
  MPI_Allgatherv(seg_start.data(),total_num_segs,MPI_INT,buff_srt,num_of_segs_per_rank,seg_disp,MPI_INT,comm.mpi_comm());
  MPI_Allgatherv(seg_length.data(),total_num_segs,MPI_INT,buff_len,num_of_segs_per_rank,seg_disp,MPI_INT,comm.mpi_comm());
  // Now that all ranks have all segments, take only the segments this rank cares about.
  // First we clear the output again and repopulate the them with the correct segments.
  for (int ii=0;ii<total_num_segs;ii++) {
    // TODO, the following loop may be faster using Kokkos::parallel_reduce
    for (int jj=0; jj<dofs_gids.extent(0);jj++) {
      if (dofs_gids[jj]-target_dof_offset == buff_dof[ii]-source_dof_offset) {
        // Using offsets to recast the seg_dof to match the target_grid_offset.
        seg_dof.push_back   (buff_dof[ii]-source_dof_offset+target_dof_offset);
        seg_start.push_back (buff_srt[ii]);
        seg_length.push_back(buff_len[ii]);
        break;
      }
    }
  }
}

void GSMap::
construct_data_segments (
  const Int                       my_start,
  const view_1d<Int>::HostMirror& target_dofs,
  std::vector<Int>&               seg_dof,
  std::vector<Int>&               seg_start,
  std::vector<Int>&               seg_length
)
{
  seg_dof.clear();
  seg_start.clear();
  seg_length.clear();
  // Start building individual segments from the data on this rank.
  seg_dof.push_back(target_dofs(0));
  seg_start.push_back(my_start);
  seg_length.push_back(1);
  for (int ii=1;ii<target_dofs.extent(0);ii++) {
    if (target_dofs(ii)==seg_dof.back()) {
      // Then we expand the segment by 1
      seg_length.back() ++;
    } else {
      // Otherwise it is time to start a new segment
      seg_dof.push_back(target_dofs(ii));
      seg_start.push_back(my_start+ii);
      seg_length.push_back(1);
    }
  }
}

/* -----------------------------------------------------------------------*/
void GSMap::
set_data_chunksize (
  const ekat::Comm& comm,
  const Int         total_length,
        Int         my_chunk_size,
        Int         my_start)
{
  // Gather my rank
  const Int my_rank      = comm.rank();
  const Int num_of_ranks = comm.size();
  // Check that length of the data > num_of_ranks
  EKAT_REQUIRE_MSG(num_of_ranks<=total_length,"Error: read_remap_data on file, the data has less values than then the size of the communicator group.");
  // Assign a chunk of data to this rank
  my_chunk_size = total_length/num_of_ranks;
  if (total_length != my_chunk_size*num_of_ranks) {
    // Then there are still extra data values not yet assigned to a rank.  Distribute
    // these by adding one to the chunk size for each rank less than the remainder.
    Int remainder = total_length - (my_chunk_size*num_of_ranks);
    if (comm.rank()<remainder) {
      my_chunk_size+=1;
    }
  }
  // Assign a starting point for my chunk of data
  int* chunk_sizes_global = new int[num_of_ranks];
  comm.all_gather(&my_chunk_size,chunk_sizes_global,1);
  my_start = 0;
  for (int ii=0;ii<my_rank;ii++) {
    my_start += chunk_sizes_global[ii];
  }
  // Sanity check, my_start + my_chunk_size shouldn't be larger than the total size
  EKAT_REQUIRE_MSG(my_start+my_chunk_size < total_length,"Error: read_remap_data on file " + remap_file_name + " for rank " + std::to_string(my_rank) + 
    ".\n my_start of " + std::to_string(my_start) + "\n + my chunksize of " + std::to_string(my_chunk_size) + "\n is larger than the total amount of data, " +
    "which is " + std::to_string(total_length) +".");
  { // Sanity check, sum all chunks and make sure it matches total_length
    int check = 0;
    for (int ii=0;ii<num_of_ranks;ii++) {
      check += chunk_sizes_global[ii];
    }
    EKAT_REQUIRE_MSG(check==total_length,"ERROR: get_remap_indices - Something went wrong distributing remap data among the MPI ranks");
  }

}
/* -----------------------------------------------------------------------*/
void GSMap::
set_map(
  const std::vector<Int>&               seg_dof,
  const std::vector<Int>&               seg_start,
  const std::vector<Int>&               seg_length
)
{
  // Determine the overall length of data on this rank
  Int length = 0;
  for (int ii=0;ii<seg_length.size();ii++) {
    length += seg_length[ii];
  }
  // Allocate GSMap remap views
  source_grid_dof = view_1d<gid_type>("",length); 
  target_grid_dof = view_1d<gid_type>("",length); 
  map_weight      = view_1d<Real>("",length); 
  // Assign mapping values
  // Using segment info build the degree's of freedom in remap file to load on this rank
  std::vector<int> var_dof;
  gid_type idx=0;
  for (int seg=0;seg<seg_length.size();seg++) {
    for (int ii=0;ii<seg_length[seg];ii++) {
      var_dof.push_back(seg_start[seg]+ii);
      target_grid_dof(idx) = seg_dof[seg];
      idx ++;
    }
  }
  // Read just the data needed on this rank.
  scorpio::register_file(remap_file_name,scorpio::Read);
  std::vector<std::string> vec_of_dims = {"n_s"};
  std::string r_decomp = "Real-n_s";
  std::string i_decomp = "Int-n_s";
  scorpio::get_variable(remap_file_name, "S", "S", vec_of_dims.size(), vec_of_dims, PIO_REAL, r_decomp);    // weights
  scorpio::get_variable(remap_file_name, "col", "col", vec_of_dims.size(), vec_of_dims, PIO_INT, i_decomp); // source column id
  scorpio::set_dof(remap_file_name,"S",var_dof.size(),var_dof.data());
  scorpio::set_dof(remap_file_name,"col",var_dof.size(),var_dof.data());
  scorpio::set_decomp(remap_file_name);
  scorpio::grid_read_data_array(remap_file_name,"S",0,  map_weight.data()  ); 
  scorpio::grid_read_data_array(remap_file_name,"col",0,source_grid_dof.data());
  scorpio::eam_pio_closefile(remap_file_name);
  // Shift source column indices to be zero-based and determine the set of unique source columns
  std::vector<gid_type> unique_src;
  for (int ii=0;ii<length;ii++) {
    source_grid_dof(ii) -= source_dof_offset;
    if (std::find(unique_src.begin(),unique_src.end(),source_grid_dof(ii)) != unique_src.end()) {
      unique_src.push_back(source_grid_dof(ii));
    }
  }
  unique_source_dof = view_1d<gid_type>(unique_src.data());
}
/* -----------------------------------------------------------------------*/

} // end namespace scream
