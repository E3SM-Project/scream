#include "share/util/scream_vertical_interpolation.hpp"
#include "share/util/scream_universal_constants.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>
#include <ekat/kokkos/ekat_subview_utils.hpp>

namespace scream {

template<int PackSize>
VerticalInterpolation<PackSize>::
VerticalInterpolation (const int ncols, const int nlevs_src, const int nlevs_tgt)
 : VerticalInterpolation(ncols,nlevs_src,nlevs_tgt,constants::DefaultFillValue<Real>::value)
{
  // Nothing to do here
}
template<int PackSize>
VerticalInterpolation<PackSize>::
VerticalInterpolation (const int ncols, const int nlevs_src, const int nlevs_tgt, const Real mask_val)
 : m_lin_interp(ncols,nlevs_src,nlevs_tgt)
{
  // Nothing to do here
}

// =================== SETUP (NON-PACKED) ================== //
template<int PackSize>
void VerticalInterpolation<PackSize>::
setup (const view_2d<const Real>& x_src,
       const view_2d<const Real>& x_tgt)
{
  check_packable (x_src,"src X view");
  check_packable (x_tgt,"tgt X view");

  // Pack inputs, and call packed version
  using PI = ekat::PackInfo<PS>;
  view_2d<const PackT> x_src_p (pack_ptr(x_src.data()),x_src.extent(0),PI::npacks(x_src.extent(1)));
  view_2d<const PackT> x_tgt_p (pack_ptr(x_tgt.data()),x_tgt.extent(0),PI::npacks(x_tgt.extent(1)));
  setup (x_src_p, x_tgt_p);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
setup (const view_2d<const Real>& x_src,
       const view_1d<const Real>& x_tgt)
{
  check_packable (x_src,"src X view");
  check_packable (x_tgt,"tgt X view");

  // Pack inputs, and call packed version
  using PI = ekat::PackInfo<PS>;
  view_2d<const PackT> x_src_p (pack_ptr(x_src.data()),x_src.extent(0),PI::npacks(x_src.extent(1)));
  view_2d<const PackT> x_tgt_p (pack_ptr(x_tgt.data()),PI::npacks(x_tgt.extent(0)));
  setup (x_src_p, x_tgt_p);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
setup (const view_1d<const Real>& x_src,
       const view_2d<const Real>& x_tgt)
{
  check_packable (x_src,"src X view");
  check_packable (x_tgt,"tgt X view");

  // Pack inputs, and call packed version
  using PI = ekat::PackInfo<PS>;
  view_2d<const PackT> x_src_p (pack_ptr(x_src.data()),PI::npacks(x_src.extent(0)));
  view_2d<const PackT> x_tgt_p (pack_ptr(x_tgt.data()),x_tgt.extent(0),PI::npacks(x_tgt.extent(1)));
  setup (x_src_p, x_tgt_p);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
setup (const view_2d<const Real>& x_src,
       const Real x_tgt)
{
  EKAT_REQUIRE_MSG (PackSize==1,
      "Error! When doing vertical interpolation to a *single* vertical level, use PackSize=1.\n");

  check_packable (x_src,"src X view");

  // Pack inputs, and call packed version
  using PI = ekat::PackInfo<PS>;
  view_2d<const PackT> x_src_p (pack_ptr(x_src.data()),x_src.extent(0),PI::npacks(x_src.extent(1)));
  setup (x_src_p, x_tgt);
}

// =================== SETUP (PACKED) ================== //

template<int PackSize>
void VerticalInterpolation<PackSize>::
setup (const view_2d<const PackT>& x_src,
       const view_2d<const PackT>& x_tgt)
{
  check_ncols(x_src,"src X view");
  check_ncols(x_tgt,"tgt X view");
  check_packed_len(x_src,"src X view");
  check_packed_len(x_tgt,"tgt X view");

  const auto policy = m_lin_interp.policy();
  auto vert_interp = m_lin_interp;
  auto li_setup = KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
  
    vert_interp.setup(team,
                      ekat::subview(x_src,icol),
                      ekat::subview(x_tgt,icol));
  };
  Kokkos::parallel_for("vert_interp_setup", policy, li_setup);
  Kokkos::fence();

  m_x_src_2d = x_src;
  m_x_tgt_2d = x_tgt;

  m_xtype = Src2d_Tgt2d;
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
setup (const view_2d<const PackT>& x_src,
       const view_1d<const PackT>& x_tgt)
{
  check_ncols(x_src,"src X view");
  check_packed_len(x_src,"src X view");
  check_packed_len(x_tgt,"tgt X view");

  const auto policy = m_lin_interp.policy();
  auto vert_interp = m_lin_interp;
  auto li_setup = KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
  
    vert_interp.setup(team,
                      ekat::subview(x_src,icol),
                      x_tgt);
  };
  Kokkos::parallel_for("vert_interp_setup", policy, li_setup);
  Kokkos::fence();

  m_x_src_2d = x_src;
  m_x_tgt_1d = x_tgt;

  m_xtype = Src2d_Tgt1d;
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
setup (const view_1d<const PackT>& x_src,
       const view_2d<const PackT>& x_tgt)
{
  check_ncols(x_tgt,"tgt X view");
  check_packed_len(x_src,"src X view");
  check_packed_len(x_tgt,"tgt X view");

  const auto policy = m_lin_interp.policy();
  auto vert_interp = m_lin_interp;
  auto li_setup = KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
  
    vert_interp.setup(team,
                      x_src,
                      ekat::subview(x_tgt,icol));
  };
  Kokkos::parallel_for("vert_interp_setup", policy, li_setup);
  Kokkos::fence();

  m_x_src_1d = x_src;
  m_x_tgt_2d = x_tgt;

  m_xtype = Src1d_Tgt2d;
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
setup (const view_2d<const PackT>& x_src,
       const Real x_tgt)
{
  check_ncols(x_src,"src X view");
  check_packed_len(x_src,"src X view");

  const auto policy = m_lin_interp.policy();
  view_1d<Real> x_tgt_view ("",1);
  Kokkos::deep_copy(x_tgt_view,x_tgt);

  auto vert_interp = m_lin_interp;
  auto li_setup = KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
  
    vert_interp.setup(team,
                      ekat::subview(x_src,icol),
                      x_tgt_view);
  };
  Kokkos::parallel_for("vert_interp_setup", policy, li_setup);
  Kokkos::fence();

  m_x_src_1d = x_src;
  m_x_tgt_0d = x_tgt;

  m_xtype = Src2d_Tgt0d;
}

// =================== INTERPOLATE (NON-PACKED) ================== //

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate (const view_2d<const Real>& y_src,
             const view_2d<      Real>& y_tgt,
             const view_2d<      MaskVal>& mask)
{
  check_packable (y_src, "src Y view");
  check_packable (y_tgt, "tgt Y view");

  // Pack inputs, and call packed version
  using PI = ekat::PackInfo<PS>;
  view_2d<const PackT> y_src_p (pack_ptr(y_src.data()),y_src.extent(0),PI::npacks(y_src.extent(1)));
  view_2d<      PackT> y_tgt_p (pack_ptr(y_tgt.data()),y_tgt.extent(0),PI::npacks(y_tgt.extent(1)));

  view_2d<MaskT> mask_p;
  if (mask.data()!=nullptr) {
    check_packable (mask);
    mask_p = view_2d<MaskT>(mask_ptr(mask.data()),mask.extent(0),PI::npacks(mask.extent(1)));
  }
  interpolate (y_src_p, y_tgt_p, mask_p);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate (const view_3d<const Real>& y_src,
             const view_3d<      Real>& y_tgt,
             const view_2d<      MaskVal>& mask)
{
  check_packable (y_src, "src Y view");
  check_packable (y_tgt, "tgt Y view");

  // Pack inputs, and call packed version
  using PI = ekat::PackInfo<PS>;
  view_3d<const PackT> y_src_p (pack_ptr(y_src.data()),y_src.extent(0),y_src.extent(1),PI::npacks(y_src.extent(2)));
  view_3d<      PackT> y_tgt_p (pack_ptr(y_tgt.data()),y_tgt.extent(0),y_tgt.extent(1),PI::npacks(y_tgt.extent(2)));

  view_2d<MaskT> mask_p;
  if (mask.data()!=nullptr) {
    check_packable (mask);
    mask_p = view_2d<MaskT>(mask_ptr(mask.data()),mask.extent(0),PI::npacks(mask.extent(1)));
  }
  interpolate (y_src_p, y_tgt_p, mask_p);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate (const view_2d<const Real>& y_src,
             const view_1d<      Real>& y_tgt,
             const view_1d<      MaskVal>& mask)
{
  check_packable (y_src, "src Y view");

  // Pack inputs, and call packed version
  using PI = ekat::PackInfo<PS>;
  view_2d<const PackT> y_src_p (pack_ptr(y_src.data()),y_src.extent(0),PI::npacks(y_src.extent(1)));
  view_1d<      PackT> y_tgt_p (pack_ptr(y_tgt.data()),y_tgt.extent(0));

  view_2d<MaskT> mask_p;
  if (mask.data()!=nullptr) {
    check_packable (mask);
    mask_p = view_2d<MaskT>(mask_ptr(mask.data()),mask.extent(0),PI::npacks(mask.extent(2)));
  }
  interpolate (y_src_p, y_tgt_p, mask_p);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate (const view_3d<const Real>& y_src,
             const view_2d<      Real>& y_tgt,
             const view_1d<      MaskVal>& mask)
{
  check_packable (y_src, "src Y view");
  check_packable (y_tgt, "tgt Y view");

  // Pack inputs, and call packed version
  using PI = ekat::PackInfo<PS>;
  view_2d<const PackT> y_src_p (pack_ptr(y_src.data()),y_src.extent(0),y_src.extent(1),PI::npacks(y_src.extent(2)));
  view_1d<      PackT> y_tgt_p (pack_ptr(y_tgt.data()),y_tgt.extent(0),y_tgt.extent(1));

  view_1d<MaskT> mask_p;
  if (mask.data()!=nullptr) {
    check_packable (mask);
    mask_p = view_1d<MaskT>(mask_ptr(mask.data()),mask.extent(0));
  }
  interpolate (y_src_p, y_tgt_p, mask_p);
}

// =================== INTERPOLATE (PACKED) ================== //

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate (const view_2d<const PackT>& y_src,
             const view_2d<      PackT>& y_tgt,
             const view_2d<      MaskT>& mask)
{
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  EKAT_REQUIRE_MSG (m_xtype!=Src2d_Tgt0d,
      "Error! VerticalInterpolation was NOT setup to interp to a single level, but input y views have different rank.\n");

  check_packed_len (y_src, "src Y view");
  check_packed_len (y_tgt, "tgt Y view");
  check_ncols (y_src, "src Y view");
  check_ncols (y_tgt, "tgt Y view");

  m_y_src_2d = y_src;
  m_y_tgt_2d = y_tgt;
  m_mask_2d = mask;

  auto policy = ESU::get_default_team_policy<TagSameDim>(m_ncols,m_lin_interp.km2_pack());
  Kokkos::parallel_for("vert_interp_2d_2d",policy,*this);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate (const view_3d<const PackT>& y_src,
             const view_3d<      PackT>& y_tgt,
             const view_2d<      MaskT>& mask)
{
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  EKAT_REQUIRE_MSG (m_xtype!=Src2d_Tgt0d,
      "Error! VerticalInterpolation was NOT setup to interp to a single level, but input y views have different rank.\n");
  EKAT_REQUIRE_MSG (y_src.extent(1) == y_tgt.extent(1),
      "Error! Cannot setup VerticalInterpolation object due to vec dim incompatibility.\n"
      " - src Y view vec dim extent: " + std::to_string(y_src.extent(1)) + "\n"
      " - tgt Y view vec dim extent: " + std::to_string(y_tgt.extent(1)) + "\n");

  check_packed_len (y_src, "src Y view");
  check_packed_len (y_tgt, "tgt Y view");
  check_ncols (y_src, "src Y view");
  check_ncols (y_tgt, "tgt Y view");

  m_y_src_3d = y_src;
  m_y_tgt_3d = y_tgt;
  m_mask_2d = mask;

  auto policy = ESU::get_default_team_policy<TagSameDim>(m_ncols,m_lin_interp.km2_pack());
  Kokkos::parallel_for("vert_interp_3d_3d",policy,*this);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate (const view_2d<const PackT>& y_src,
             const view_1d<      PackT>& y_tgt,
             const view_1d<      MaskT>& mask)
{
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  EKAT_REQUIRE_MSG (m_xtype==Src2d_Tgt0d,
      "Error! VerticalInterpolation was setup to interp to a single level, but input y views have same rank.\n");

  check_packed_len (y_src, "src Y view");
  check_ncols (y_src, "src Y view");
  check_ncols (y_tgt, "tgt Y view");

  m_y_src_2d = y_src;
  m_y_tgt_1d = y_tgt;
  m_mask_1d = mask;

  auto policy = ESU::get_default_team_policy<TagVertSlice>(m_ncols,m_lin_interp.km2_pack());
  Kokkos::parallel_for("vert_interp_2d_1d",policy,*this);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate (const view_3d<const PackT>& y_src,
             const view_2d<      PackT>& y_tgt,
             const view_1d<      MaskT>& mask)
{
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  EKAT_REQUIRE_MSG (m_xtype==Src2d_Tgt0d,
      "Error! VerticalInterpolation was setup to interp to a single level, but input y views have same rank.\n");
  EKAT_REQUIRE_MSG (y_src.extent(1) == y_tgt.extent(1),
      "Error! Cannot setup VerticalInterpolation object due to vec dim incompatibility.\n"
      " - src Y view vec dim extent: " + std::to_string(y_src.extent(1)) + "\n"
      " - tgt Y view vec dim extent: " + std::to_string(y_tgt.extent(1)) + "\n");

  check_packed_len (y_src, "src Y view");
  check_ncols (y_src, "src Y view");
  check_ncols (y_tgt, "tgt Y view");

  m_y_src_3d = y_src;
  m_y_tgt_2d = y_tgt;
  m_mask_1d = mask;

  auto policy = ESU::get_default_team_policy<TagVertSlice>(m_ncols*y_src.extent(1),m_lin_interp.km2_pack());
  Kokkos::parallel_for("vert_interp_3d_2d",policy,*this);
}

template<int PackSize>
KOKKOS_FORCEINLINE_FUNCTION
void VerticalInterpolation<PackSize>::
operator()(const TagSameDim&, const MemberType& team)
{
  const int icol = team.league_rank() % m_ncols;
  const int ivar = team.league_rank() / m_ncols;

  const bool vector = m_y_src_3d.data()!=nullptr;

  auto x_src_col = m_xtype==Src1d_Tgt2d ? m_x_src_1d : ekat::subview(m_x_src_2d,icol);
  auto x_tgt_col = m_xtype==Src2d_Tgt1d ? m_x_tgt_1d : ekat::subview(m_x_tgt_2d,icol);
  auto y_src_col = vector ? ekat::subview(m_y_src_3d,icol,ivar) : ekat::subview(m_y_src_2d,icol);
  auto y_tgt_col = vector ? ekat::subview(m_y_tgt_3d,icol,ivar) : ekat::subview(m_y_tgt_2d,icol);
  view_1d<MaskT> mask;
  if (m_mask_2d.data()!=nullptr) {
    mask = ekat::subview(m_mask_2d,icol);
  }
  interpolate_1d(team, icol, x_src_col, x_tgt_col, y_src_col, y_tgt_col, mask);
}

template<int PackSize>
KOKKOS_FORCEINLINE_FUNCTION
void VerticalInterpolation<PackSize>::
operator()(const TagVertSlice&, const MemberType& team)
{
  const int icol = team.league_rank() % m_ncols;
  const int ivar = team.league_rank() / m_ncols;

  const bool vector = m_y_src_3d.data()!=nullptr;

  auto  x_src_col = ekat::subview(m_x_src_2d,icol);
  auto  y_src_col = vector ? ekat::subview(m_y_src_3d,icol,ivar) : ekat::subview(m_y_src_2d,icol);
  auto& y_tgt_val = vector ? m_y_tgt_2d(icol,ivar) : m_y_tgt_1d(icol);
  MaskT* mask = nullptr;
  if (m_mask_1d.data()!=nullptr) {
    mask = &m_mask_1d(icol);
  }
  interpolate_0d(team, icol, x_src_col, y_src_col, y_tgt_val, mask);
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate_1d (const MemberType& team,
                const int icol,
                const view_1d<const PackT>& x_src,
                const view_1d<const PackT>& x_tgt,
                const view_1d<const PackT>& y_src,
                const view_1d<      PackT>& y_tgt,
                const view_1d<      MaskT>& mask)
{
  using PI = ekat::PackInfo<PS>;

  m_lin_interp.lin_interp(team, x_src, x_tgt, y_src, y_tgt, icol);

  // If the mask view is valid, compute mask and reset y_tgt if mask=true
  if (mask.data()!=nullptr) {
    auto range = Kokkos::TeamVectorRange(team, y_tgt.extent(0));
    
    auto min = x_src(0)[0];
    auto max = x_src(PI::last_pack_idx(m_nlevs_src))[PI::last_vec_end(m_nlevs_src)];
    auto set_mask = [&](const int k) {
      mask(k) = x_tgt(k)<min or x_tgt(k)>max;
      if (mask(k).any()) {
        y_tgt(k).set(mask(k),m_mask_val);
      }
    };
  }
}

template<int PackSize>
void VerticalInterpolation<PackSize>::
interpolate_0d (const MemberType& team,
                const int icol,
                const view_1d<const PackT>& x_src,
                const view_1d<const PackT>& y_src,
                                    PackT & y_tgt,
                                    MaskT* mask)
{
  using PI = ekat::PackInfo<PS>;

  view_1d<const PackT> x_tgt_v(&m_x_tgt_0d,1);
  view_1d<      PackT> y_tgt_v(&y_tgt,1);

  m_lin_interp.lin_interp(team, x_src, x_tgt_v, y_src, y_tgt_v, icol);

  // If the mask view is valid, compute mask and reset y_tgt if mask=true
  if (mask!=nullptr) {
    auto& m = *mask;
    auto min = x_src(0)[0];
    auto max = x_src(PI::last_pack_idx(m_nlevs_src))[PI::last_vec_end(m_nlevs_src)];
    Kokkos::single(Kokkos::PerTeam(team),[&]{
      m = m_x_tgt_0d<min or m_x_tgt_0d>max;
      if (m.any()) {
        y_tgt.set(m,m_mask_val);
      }
    });
  }
}

// =================== HELPERS ================== //

template<int PackSize>
template<typename ViewT>
void VerticalInterpolation<PackSize>::
check_packable (const ViewT& v, const std::string& view_type)
{
  const int len = v.extent(v.Rank-1);
  EKAT_REQUIRE_MSG (len % PS == 0,
      "Error! Cannot setup VerticalInterpolation object due to pack size incompatibility.\n"
      " - " + view_type + " vertical dim extent: " + std::to_string(len) + "\n"
      " - VerticalInterpolation pack size: " + std::to_string(PS) + "\n");
}

template<int PackSize>
template<typename ViewT>
void VerticalInterpolation<PackSize>::
check_packed_len (const ViewT& v, const std::string& view_type)
{
  const bool is_src = view_type=="src X view" or view_type=="src Y view";
  const int li_npacks = is_src ? m_lin_interp.km1_pack() : m_lin_interp.km2_pack();
  const int len = v.extent(v.Rank-1);
  const std::string type = is_src ? "src" : "tgt";
  EKAT_REQUIRE_MSG (len == li_npacks,
      "Error! Cannot setup VerticalInterpolation object due to number of vertical packs incompatibility.\n"
      " - " + view_type + " num vertical packs: " + std::to_string(len) + "\n"
      " - LinInterp " + type + " num vertical packs: " + std::to_string(li_npacks) + "\n");
}

template<int PackSize>
template<typename ViewT>
void VerticalInterpolation<PackSize>::
check_ncols (const ViewT& v, const std::string& view_type)
{
  EKAT_REQUIRE_MSG (v.extent(0) == m_ncols,
      "Error! Cannot setup VerticalInterpolation object due to ncols incompatibility.\n"
      " - " + view_type + " ncols extent: " + std::to_string(v.extent(0)) + "\n"
      " - ncols used at construction time: " + std::to_string(m_ncols) + "\n");
}

} // namespace scream
