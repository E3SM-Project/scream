#include <catch2/catch.hpp>

#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/vertical_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include <limits>

namespace scream {

constexpr Real FILL_VALUE = std::numeric_limits<float>::max()/1e5;

template<typename ViewT>
typename ViewT::HostMirror
cmvc (const ViewT& v) {
  auto vh = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(vh,v);
  return vh;
}

class VerticalRemapperTester : public VerticalRemapper {
public:
  VerticalRemapperTester (const grid_ptr_type& src_grid,
                          const std::string&   map_file,
                          const Field&         lev_prof,
                          const Field&         ilev_prof,
                          const Real           mask_val)
   : VerticalRemapper(src_grid, map_file, lev_prof, ilev_prof, mask_val)
  {
    // Nothing to do
  }
};

template<typename ViewT>
bool contains (const ViewT& v, const typename ViewT::traits::value_type& entry) {
  const auto vh = cmvc (v);
  const auto beg = vh.data();
  const auto end = vh.data() + vh.size();
  for (auto it=beg; it!=end; ++it) {
    if (*it == entry) {
      return true;
    }
  }
  return false;
}

class CoarseningRemapperTester : public CoarseningRemapper {
public:
  CoarseningRemapperTester (const grid_ptr_type& src_grid,
                            const std::string& map_file,
			    const bool track_mask = false)
   : CoarseningRemapper(src_grid,map_file,track_mask)
  {
    // Nothing to do
  }
  std::vector<gid_t>
  test_triplet_gids (const std::string& map_file) const {
    return CoarseningRemapper::get_my_triplets_gids (map_file,m_src_grid);
  }

  view_1d<int> get_row_offsets () const {
    return m_row_offsets;
  }
  view_1d<int> get_col_lids () const {
    return m_col_lids;
  }
  view_1d<Real> get_weights () const {
    return m_weights;
  }

  grid_ptr_type get_ov_tgt_grid () const {
    return m_ov_tgt_grid;
  }

  view_2d<int>::HostMirror get_send_f_pid_offsets () const {
    return cmvc(m_send_f_pid_offsets);
  }
  view_2d<int>::HostMirror get_recv_f_pid_offsets () const {
    return cmvc(m_recv_f_pid_offsets);
  }

  view_1d<int>::HostMirror get_recv_lids_beg () const {
    return cmvc(m_recv_lids_beg);
  }
  view_1d<int>::HostMirror get_recv_lids_end () const {
    return cmvc(m_recv_lids_end);
  }

  view_2d<int>::HostMirror get_send_lids_pids () const {
    return cmvc(m_send_lids_pids );
  }
  view_2d<int>::HostMirror get_recv_lids_pidpos () const {
    return cmvc(m_recv_lids_pidpos);
  }

  view_1d<int>::HostMirror get_send_pid_lids_start () const {
    return cmvc(m_send_pid_lids_start);
  }

  int gid2lid (const gid_t gid, const grid_ptr_type& grid) const {
    return CoarseningRemapper::gid2lid(gid,grid);
  }
};

template<typename ViewT>
bool view_contains (const ViewT& v, const typename ViewT::traits::value_type& entry) {
  const auto vh = cmvc (v);
  const auto beg = vh.data();
  const auto end = vh.data() + vh.size();
  for (auto it=beg; it!=end; ++it) {
    if (*it == entry) {
      return true;
    }
  }
  return false;
}

void print (const std::string& msg, const ekat::Comm& comm) {
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

// Helper function to create a grid given the number of dof's and a comm group.
std::shared_ptr<AbstractGrid>
build_src_grid(const ekat::Comm& comm, const int nldofs_src, const int nlevs_src) 
{
  auto src_grid = std::make_shared<PointGrid>("src",nldofs_src,nlevs_src,comm);

  auto src_dofs = src_grid->get_dofs_gids();
  auto src_dofs_h = src_dofs.get_view<gid_t*,Host>();
  std::iota(src_dofs_h.data(),src_dofs_h.data()+nldofs_src,nldofs_src*comm.rank());
  src_dofs.sync_to_dev();

  return src_grid;
}

// Helper function to create fields
Field
create_field(const std::string& name, const std::shared_ptr<const AbstractGrid>& grid, const bool twod, const bool vec, const bool mid = false, const int ps = 1, const bool add_mask = false)
{
  constexpr int vec_dim = 3;
  constexpr auto CMP = FieldTag::Component;
  constexpr auto units = ekat::units::Units::nondimensional();
  auto fl = twod
          ? (vec ? grid->get_2d_vector_layout (CMP,vec_dim)
                 : grid->get_2d_scalar_layout ())
          : (vec ? grid->get_3d_vector_layout (mid,CMP,vec_dim)
                 : grid->get_3d_scalar_layout (mid));
  FieldIdentifier fid(name,fl,units,grid->name());
  Field f(fid);
  f.get_header().get_alloc_properties().request_allocation(ps);
  f.allocate_view();

  if (add_mask) {
    // Add a mask to the field
    FieldIdentifier fid_mask(name+"_mask",fl,units,grid->name());
    Field f_mask(fid_mask);
    f_mask.get_header().get_alloc_properties().request_allocation(ps);
    f_mask.allocate_view();
    f_mask.deep_copy(1.0);
    f.get_header().set_extra_data("mask_data",f_mask);
    f.get_header().set_extra_data("mask_value",FILL_VALUE);
  }

  return f;
}

// Helper function to create a remap file
void create_remap_file(const std::string& filename, std::vector<std::int64_t>& dofs,
                       const int na, const int nb, const int ns,
                       const std::vector<Real>& col, const std::vector<Real>& row, const std::vector<Real>& S,
		       const int nlevs, const std::vector<std::int64_t>& dofs_p, const std::vector<Real>& p_tgt) 
{

  scorpio::register_file(filename, scorpio::FileMode::Write);

  scorpio::register_dimension(filename,"n_a", "n_a", na, true);
  scorpio::register_dimension(filename,"n_b", "n_b", nb, true);
  scorpio::register_dimension(filename,"n_s", "n_s", ns, true);
  scorpio::register_dimension(filename,"nlevs","nlevs",nlevs, false);

  scorpio::register_variable(filename,"col","col","none",{"n_s"},"real","int","int-nnz");
  scorpio::register_variable(filename,"row","row","none",{"n_s"},"real","int","int-nnz");
  scorpio::register_variable(filename,"S","S","none",{"n_s"},"real","real","Real-nnz");
  scorpio::register_variable(filename,"p_levs","p_levs","none",{"nlevs"},"real","real","Real-nlevs");

  scorpio::set_dof(filename,"col",dofs.size(),dofs.data());
  scorpio::set_dof(filename,"row",dofs.size(),dofs.data());
  scorpio::set_dof(filename,"S",  dofs.size(),dofs.data());
  scorpio::set_dof(filename,"p_levs",dofs_p.size(),dofs_p.data()); 
  
  scorpio::eam_pio_enddef(filename);

  scorpio::grid_write_data_array(filename,"row",row.data(),ns);
  scorpio::grid_write_data_array(filename,"col",col.data(),ns);
  scorpio::grid_write_data_array(filename,"S",    S.data(),ns);
  scorpio::grid_write_data_array(filename,"p_levs",p_tgt.data(),nlevs);

  scorpio::eam_pio_closefile(filename);
}

TEST_CASE ("vertical_plus_coarsening_remap") {
  using gid_t = AbstractGrid::gid_type;

  // -------------------------------------- //
  //           Init MPI and PIO             //
  // -------------------------------------- //

  ekat::Comm comm(MPI_COMM_WORLD);

  MPI_Fint fcomm = MPI_Comm_c2f(comm.mpi_comm());
  scorpio::eam_init_pio_subsystem(fcomm);

  // -------------------------------------- //
  //           Set grid/map sizes           //
  // -------------------------------------- //

  const int nldofs_src = 10;
  const int nlevs_src  = 2*SCREAM_PACK_SIZE + 2;  // Make sure we check what happens when the vertical extent is a little larger than the max PACK SIZE
  const int nldofs_tgt = 5;
  const int nlevs_tgt  = nlevs_src/2;
  const int ngdofs_src = nldofs_src*comm.size();
  const int ngdofs_tgt = nldofs_tgt*comm.size();
  const int nnz_local  = nldofs_src;
  const int nnz        = nnz_local*comm.size();

  // -------------------------------------- //
  //           Create a map file            //
  // -------------------------------------- //

  print (" -> creating map file ...\n",comm);

  std::string filename = "vertical_plus_coarsening_map_file_np" + std::to_string(comm.size()) + ".nc";
  std::vector<std::int64_t> dofs (nnz_local);
  std::iota(dofs.begin(),dofs.end(),comm.rank()*nnz_local);

  // Create triplets: tgt entry K is the avg of src entries K and K+ngdofs_tgt
  // NOTE: add 1 to row/col indices, since e3sm map files indices are 1-based
  std::vector<Real> col,row,S;
  for (int i=0; i<nldofs_tgt; ++i) {
    row.push_back(1+i+nldofs_tgt*comm.rank());
    col.push_back(1+i+nldofs_tgt*comm.rank());
    S.push_back(0.25);

    row.push_back(1+i+nldofs_tgt*comm.rank());
    col.push_back(1+i+nldofs_tgt*comm.rank() + ngdofs_tgt);
    S.push_back(0.75);
  }

  // Create target pressure levels to be remapped onto
  const Real ptop_tgt = 1;
  const Real pbot_tgt = nlevs_src-2;
  const Real dp_tgt   = (pbot_tgt-ptop_tgt)/(nlevs_tgt-1);
  std::vector<std::int64_t> dofs_p(nlevs_tgt);
  std::iota(dofs_p.begin(),dofs_p.end(),0);
  std::vector<Real> p_tgt;
  for (int ii=0; ii<nlevs_tgt; ++ii) {
    p_tgt.push_back(ptop_tgt + dp_tgt*ii);
  }

  create_remap_file(filename, dofs, ngdofs_src, ngdofs_tgt, nnz, col, row, S, nlevs_tgt, dofs_p, p_tgt);
  print (" -> creating map file ... done!\n",comm);

  // -------------------------------------- //
  //      Build src grid and remapper       //
  // -------------------------------------- //

  print (" -> creating grid and remapper ...\n",comm);

  auto src_grid = build_src_grid(comm, nldofs_src, nlevs_src);

  // We need the source pressure level fields for both p_mid and p_int
  auto pmid_src   = create_field("p_mid",  src_grid, false, false, true,  SCREAM_PACK_SIZE);
  auto pint_src   = create_field("p_int",  src_grid, false, false, false, SCREAM_PACK_SIZE);
  // Set the source pressures
  {
    // By adding 1 to the pbot_tgt and subtrating 1 from ptop_tgt we ensure some masking, which 
    // we also want to check.
    const Real ptop_src = 0.0;
    const Real dp_dx = (pbot_tgt - ptop_src)/nldofs_src;
    auto pmid_v = pmid_src.get_view<Real**,Host>();
    auto pint_v = pint_src.get_view<Real**,Host>();
    for (int ii=0; ii<pmid_v.extent_int(0); ++ii) {
      // For each column, change the surface pressure so that some target values will necessarily be masked
      const Real pbot_src = pbot_tgt - ii*dp_dx;
      const Real dp_src = (pbot_src-ptop_src)/(nlevs_src-1);
      pint_v(ii,0) = ptop_src;
      for (int kk=0; kk<nlevs_src; ++kk) {
        pint_v(ii,kk+1) = pint_v(ii,kk) + dp_src;
        pmid_v(ii,kk)   = 0.5*(pint_v(ii,kk) + pint_v(ii,kk+1));
      }
    }
  }
  pmid_src.sync_to_dev();
  pint_src.sync_to_dev();

  // Vertical remapper
  auto vert_remap = std::make_shared<VerticalRemapperTester>(src_grid,filename,pmid_src,pint_src,FILL_VALUE);
  auto vert_grid = vert_remap->get_tgt_grid();
  REQUIRE(vert_grid->get_num_global_dofs()==src_grid->get_num_global_dofs());
  // Horizontal remapper - coarsening
  auto horz_remap = std::make_shared<CoarseningRemapperTester>(vert_grid,filename,true);
  auto tgt_grid  = horz_remap->get_tgt_grid();
  REQUIRE (tgt_grid->get_num_global_dofs()==ngdofs_tgt);
  print (" -> creating grid and remapper ... done!\n",comm);

  // -------------------------------------- //
  //      Create src/tgt grid fields        //
  // -------------------------------------- //

  print (" -> creating fields ...\n",comm);
  constexpr int vec_dim = 3;

  auto src_s2d   = create_field("s2d",  src_grid,true,false);
  auto src_v2d   = create_field("v2d",  src_grid,true,true);
  auto src_s3d_m = create_field("s3d_m",src_grid,false,false,true, 1);
  auto src_s3d_i = create_field("s3d_i",src_grid,false,false,false,SCREAM_PACK_SIZE);
  auto src_v3d_m = create_field("v3d_m",src_grid,false,true ,true, 1);
  auto src_v3d_i = create_field("v3d_i",src_grid,false,true ,false,SCREAM_PACK_SIZE);

  // Note that once a field has been vertically interpolated onto set pressure levels there is no
  // longer a concept of interface pressure levels, they are all on the same vertical grid which
  // we define as the LEV grid.
  auto vert_s2d   = create_field("s2d",  vert_grid,true,false);
  auto vert_v2d   = create_field("v2d",  vert_grid,true,true);
  auto vert_s3d_m = create_field("s3d_m",vert_grid,false,false,true, 1);
  auto vert_s3d_i = create_field("s3d_i",vert_grid,false,false,true,SCREAM_PACK_SIZE);
  auto vert_v3d_m = create_field("v3d_m",vert_grid,false,true ,true, 1);
  auto vert_v3d_i = create_field("v3d_i",vert_grid,false,true ,true,SCREAM_PACK_SIZE);

  auto tgt_s2d   = create_field("s2d",  tgt_grid,true,false);
  auto tgt_v2d   = create_field("v2d",  tgt_grid,true,true);
  auto tgt_s3d_m = create_field("s3d_m",tgt_grid,false,false,true, 1);
  auto tgt_s3d_i = create_field("s3d_i",tgt_grid,false,false,true,SCREAM_PACK_SIZE);
  auto tgt_v3d_m = create_field("v3d_m",tgt_grid,false,true ,true, 1);
  auto tgt_v3d_i = create_field("v3d_i",tgt_grid,false,true ,true,SCREAM_PACK_SIZE);

  std::vector<Field> src_f = {src_s2d,src_v2d,src_s3d_m,src_s3d_i,src_v3d_m,src_v3d_i};
  std::vector<Field> tgt_f = {tgt_s2d,tgt_v2d,tgt_s3d_m,tgt_s3d_i,tgt_v3d_m,tgt_v3d_i};

  const int nfields = src_f.size();

  std::vector<int> field_col_size (src_f.size());
  std::vector<int> field_col_offset (src_f.size()+1,0);
  for (int i=0; i<nfields; ++i) {
    const auto& f  = src_f[i];  // Doesn't matter if src or tgt
    const auto& fl = f.get_header().get_identifier().get_layout();
    field_col_size[i] = fl.size() / fl.dim(0);
    field_col_offset[i+1] = field_col_offset[i]+field_col_size[i];
  }

  print (" -> creating fields ... done!\n",comm);

  // -------------------------------------- //
  //     Register fields in the remapper    //
  // -------------------------------------- //

  print (" -> registering fields ...\n",comm);
  print (" -> registering fields ... vertical\n",comm);
  vert_remap->registration_begins();
  vert_remap->register_field(src_s2d,  vert_s2d);
  vert_remap->register_field(src_v2d,  vert_v2d);
  vert_remap->register_field(src_s3d_m,vert_s3d_m);
  vert_remap->register_field(src_s3d_i,vert_s3d_i);
  vert_remap->register_field(src_v3d_m,vert_v3d_m);
  vert_remap->register_field(src_v3d_i,vert_v3d_i);
  vert_remap->registration_ends();
  print (" -> registering fields ... horizontal\n",comm);
  horz_remap->registration_begins();
  horz_remap->register_field(vert_s2d,  tgt_s2d);
  horz_remap->register_field(vert_v2d,  tgt_v2d);
  horz_remap->register_field(vert_s3d_m,tgt_s3d_m);
  horz_remap->register_field(vert_s3d_i,tgt_s3d_i);
  horz_remap->register_field(vert_v3d_m,tgt_v3d_m);
  horz_remap->register_field(vert_v3d_i,tgt_v3d_i);
  horz_remap->registration_ends();
  print (" -> registering fields ... done!\n",comm);

  // -------------------------------------- //
  //        Check remapper internals        //
  // -------------------------------------- //

  print (" -> Checking remapper internal state ...\n",comm);

  // Check which triplets are read from map file
  auto src_dofs_h = src_grid->get_dofs_gids().get_view<const gid_t*,Host>();
  auto my_triplets = horz_remap->test_triplet_gids (filename);
  const int num_triplets = my_triplets.size();
  REQUIRE (num_triplets==nnz_local);
  for (int i=0; i<nnz_local; ++i) {
    const auto src_gid = src_dofs_h(i);
    const auto tgt_gid = src_gid % ngdofs_tgt;

    REQUIRE (ekat::contains(my_triplets, 2*tgt_gid + src_gid/ngdofs_tgt));
  }

  // Check overlapped tgt grid
  // NOTE: you need to treat the case of 1 rank separately, since in that case
  //       there are 2 local src dofs impacting the same tgt dof, while with 2+
  //       ranks every local src dof impacts a different tgt dof.
  auto ov_tgt_grid = horz_remap->get_ov_tgt_grid ();
  const int num_loc_ov_tgt_gids = ov_tgt_grid->get_num_local_dofs();
  const int expected_num_loc_ov_tgt_gids = ngdofs_tgt>=nldofs_src ? nldofs_src : ngdofs_tgt;
  REQUIRE (num_loc_ov_tgt_gids==expected_num_loc_ov_tgt_gids);
  const auto ov_gids = ov_tgt_grid->get_dofs_gids().get_view<const gid_t*,Host>();
  for (int i=0; i<num_loc_ov_tgt_gids; ++i) {
    if (comm.size()==1) {
      REQUIRE(ov_gids[i]==i);
    } else {
      const auto src_gid = src_dofs_h[i];
      REQUIRE (view_contains(ov_gids, src_gid % ngdofs_tgt));
    }
  }

  // Check sparse matrix
  auto row_offsets_h = cmvc(horz_remap->get_row_offsets());
  auto col_lids_h    = cmvc(horz_remap->get_col_lids());
  auto weights_h = cmvc(horz_remap->get_weights());
  auto ov_tgt_gids = ov_tgt_grid->get_dofs_gids().get_view<const gid_t*,Host>();
  auto src_gids    = horz_remap->get_src_grid()->get_dofs_gids().get_view<const gid_t*,Host>();

  REQUIRE (col_lids_h.extent_int(0)==nldofs_src);
  REQUIRE (row_offsets_h.extent_int(0)==(num_loc_ov_tgt_gids+1));
  for (int i=0; i<num_loc_ov_tgt_gids; ++i) {
    if (comm.size()==1) {
      REQUIRE (row_offsets_h(i)==(2*i));
    } else {
      REQUIRE (row_offsets_h(i)==i);
    }
  }
  REQUIRE (row_offsets_h(num_loc_ov_tgt_gids)==nldofs_src);

  for (int irow=0; irow<num_loc_ov_tgt_gids; ++irow) {
    const auto row_gid = ov_tgt_gids(irow);
    for (int innz=row_offsets_h(irow); innz<row_offsets_h(irow+1); ++innz) {
      const auto col_lid = col_lids_h(innz);
      const auto col_gid = src_gids(col_lid);
      if (row_gid==col_gid) {
        REQUIRE (weights_h(innz)==0.25);
      } else {
        REQUIRE (weights_h(innz)==0.75);
      }
    }
  }

  // Check internal MPI structures
  const int num_loc_tgt_gids = tgt_grid->get_num_local_dofs();
  const auto tgt_gids = tgt_grid->get_dofs_gids().get_view<const gid_t*,Host>();
  const auto recv_lids_beg = horz_remap->get_recv_lids_beg();
  const auto recv_lids_end = horz_remap->get_recv_lids_end();
  const auto recv_lids_pidpos = horz_remap->get_recv_lids_pidpos();
  // Rank 0 sends everything to itself
  for (int i=0; i<num_loc_tgt_gids; ++i) {
    if (comm.size()==1) {
      // Each tgt dof has one ov_tgt contribution,
      // since the matvec is fully local
      REQUIRE (recv_lids_beg(i)==i);
      REQUIRE (recv_lids_end(i)==(i+1));

      REQUIRE (recv_lids_pidpos(i,0)==comm.rank());
      REQUIRE (recv_lids_pidpos(i,1)==i);
    } else {
      // Each tgt dof has two ov_tgt contributions,
      // since the mat vec happens on two different PIDs.
      REQUIRE (recv_lids_beg(i)==2*i);
      REQUIRE (recv_lids_end(i)==(2*i+2));
      // Figure out where the contributions come from
      const auto src1 = tgt_gids(i);
      const auto src2 = src1 + ngdofs_tgt;
      const auto pid1 = src1 / nldofs_src;
      const auto pid2 = src2 / nldofs_src;
      REQUIRE (recv_lids_pidpos(2*i,0)==pid1);
      REQUIRE (recv_lids_pidpos(2*i+1,0)==pid2);
    }
  }
  print (" -> Checking remapper internal state ... OK!\n",comm);

  // -------------------------------------- //
  //       Generate data for src fields     //
  // -------------------------------------- //

  print (" -> generate src fields data ...\n",comm);
  // Generate data in a deterministic way, so that when we check results,
  // we know a priori what the input data that generated the tgt field's
  // values was, even if that data was off rank.
  // Note, we already check spatially varying data in the direct coarsening
  // remapper test.  To simplify this test we set all source data to a constant
  // wtih respects to horizontal dimension. Thus the vertical interpolation
  // should map to a masked value or a single pressure value.

  const int ntgt_gids = tgt_gids.size();
  const Real pscale = 100.0; //TODO: Make this a random value
  for (const auto& f : src_f) {
    const auto& l = f.get_header().get_identifier().get_layout();
    auto& f_extra = f.get_header().get_extra_data();
    const auto p_src = l.has_tag(FieldTag::LevelMidPoint) ? pmid_src : pint_src;
    auto pres = p_src.get_view<Real**,Host>();
    switch (get_layout_type(l.tags())) {
      case LayoutType::Scalar2D:
      {
        const auto v_src = f.get_view<Real*,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          v_src(i) = src_gids(i);
        }
      } break;
      case LayoutType::Vector2D:
      {
        const auto v_src = f.get_view<Real**,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          for (int j=0; j<vec_dim; ++j) {
            v_src(i,j) = src_gids(i)*vec_dim + j;
        }}
      } break;
      case LayoutType::Scalar3D:
      {
        const int nlevs = l.dims().back();
        const auto v_src = f.get_view<Real**,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          for (int j=0; j<nlevs; ++j) {
            v_src(i,j) = pscale*pres(i,j);
        }}
      } break;
      case LayoutType::Vector3D:
      {
        const int nlevs = l.dims().back();
        const auto v_src = f.get_view<Real***,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          for (int j=0; j<vec_dim; ++j) {
            for (int k=0; k<nlevs; ++k) {
              v_src(i,j,k) = pscale*(j+1)*pres(i,k);
        }}}
      } break;
      default:
        EKAT_ERROR_MSG ("Unexpected layout.\n");
    }
    f.sync_to_dev();
  }
  print (" -> generate src fields data ... done!\n",comm);

  auto combine = [] (const Real lhs, const Real rhs) -> Real {
    return 0.25*lhs + 0.75*rhs;
  };

  // No bwd remap
  REQUIRE_THROWS(horz_remap->remap(false));

  for (int irun=0; irun<5; ++irun) {
    print (" -> run remap ...\n",comm);
    vert_remap->remap(true);
    horz_remap->remap(true);
    print (" -> run remap ... done!\n",comm);

    // -------------------------------------- //
    //          Check remapped fields         //
    // -------------------------------------- //

    print (" -> check tgt fields ...\n",comm);
    // Recall, tgt gid K should be the avg of src gids K and K+ngdofs_tgt,
    // unless we are checking the masked version, which should just be the value of src_gids K+ndofs_tgt because the first value is masked.
    for (size_t ifield=0; ifield<tgt_f.size(); ++ifield) {
      const auto& f = tgt_f[ifield];
      const auto& f_chk = src_f[ifield];
      const auto& l = f.get_header().get_identifier().get_layout();
      const auto& l_chk = f_chk.get_header().get_identifier().get_layout();
      const auto ls = to_string(l);
      const auto p_src = l_chk.has_tag(FieldTag::LevelMidPoint) ? pmid_src : pint_src;
      auto pres_chk = p_src.get_view<Real**,Host>();
      const auto pres = p_tgt;
      std::string dots (25-ls.size(),'.');
      print ("   -> Checking field (" + f.name() + ") with layout " + to_string(l) + " " + dots + "\n",comm);

      f.sync_to_host();

      switch (get_layout_type(l.tags())) {
        case LayoutType::Scalar2D:
        {
          const auto v_tgt = f.get_view<const Real*,Host>();
          for (int i=0; i<ntgt_gids; ++i) {
            const auto gid = tgt_gids(i);
            const auto term1 = gid;
            const auto term2 = gid+ngdofs_tgt;
            REQUIRE ( v_tgt(i) == combine(term1,term2) );
          }
        } break;
        case LayoutType::Vector2D:
        {
          const auto v_tgt = f.get_view<const Real**,Host>();
          for (int i=0; i<ntgt_gids; ++i) {
            const auto gid = tgt_gids(i);
            for (int j=0; j<vec_dim; ++j) {
              const auto term1 = gid*vec_dim+j;
              const auto term2 = (gid+ngdofs_tgt)*vec_dim+j;
              REQUIRE ( v_tgt(i,j)== combine(term1,term2) );
          }}
        } break;
        case LayoutType::Scalar3D:
        {
          const int nlevs = l.dims().back();
          const auto v_tgt = f.get_view<const Real**,Host>();
          for (int i=0; i<ntgt_gids; ++i) {
            const auto gid = tgt_gids(i);
            for (int j=0; j<nlevs; ++j) {
	      // Check if this value should be fully masked
	      // TODO: Need to do a better job of checking if a value should be masked.
	      printf("ASD - (%d, %d): %f,\t%f,\t%f,\t%f,\t%f\n",i,j,v_tgt(i,j),pscale*pres[j],pres_chk(gid,j),pres_chk(gid+ngdofs_tgt,j),pres[j]);
	      if (pres_chk(gid,j)<pres[j] && pres_chk(gid+ngdofs_tgt,j)<pres[j]) {
                REQUIRE ( v_tgt(i,j) == FILL_VALUE );
	      } else {
                REQUIRE ( v_tgt(i,j) == pscale*pres[j] );
	      }
          }}
        } break;
        case LayoutType::Vector3D:
        {
          const int nlevs = l.dims().back();
          const auto v_tgt = f.get_view<const Real***,Host>();
          for (int i=0; i<ntgt_gids; ++i) {
            const auto gid = tgt_gids(i);
            for (int j=0; j<vec_dim; ++j) {
              for (int k=0; k<nlevs; ++k) {
	        // Check if this value should be fully masked
	        printf("ASD - (%d, %d): %f,\t%f,\t%f,\t%f,\t%f\n",i,j,v_tgt(i,j),pscale*pres[j],pres_chk(gid,j),pres_chk(gid+ngdofs_tgt,j),pres[j]);
	        if (pres_chk(gid,k)<pres[k] && pres_chk(gid+ngdofs_tgt,k)<pres[k]) {
                  REQUIRE ( v_tgt(i,j,k) == FILL_VALUE );
	        } else {
                  REQUIRE ( v_tgt(i,j,k) == pscale*pres[k]*(j+1) );
	        }
          }}}
        } break;
        default:
          EKAT_ERROR_MSG ("Unexpected layout.\n");
      }

      print ("   -> Checking field with layout " + to_string(l) + " " + dots + " OK!\n",comm);
    }
    print ("check tgt fields ... done!\n",comm);
  }

  // Clean up scorpio stuff
  scorpio::eam_pio_finalize();
}

} // namespace scream
