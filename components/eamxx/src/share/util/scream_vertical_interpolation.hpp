#ifndef SCREAM_VERTICAL_INTERPOLATION_HPP
#define SCREAM_VERTICAL_INTERPOLATION_HPP

#include "share/scream_types.hpp"

#include <ekat/util/ekat_lin_interp.hpp>

namespace scream {

/*
 * A lightweight class to wrap the ekat::LinInterp object
 *
 * This class is wrapping the ekat::LinInterp object, with a couple of additions:
 *  - allows to call with either views of Pack or Real: if the latter, we create
 *    a view of Pack obj on the fly, viewing the same data (and checking that
 *    the view extent allows packing)
 *  - if the user provides it, we fill an output view with a mask, expressing
 *    where the output data is masked (i.e., could not be interpolated, so a
 *    pre-defined mask value was set)
 * The interpolation can be setup in 4 different flavors, depending on src/tgt pressure:
 *  - Src2d_Tgt2d: both src and tgt are (ncols,nlevs)
 *  - Src2d_Tgt1d: src is (ncols,nlevs), while tgt is (nlevs)
 *  - Src1d_Tgt2d: src is (nlevs), while tgt is (ncols,nlevs)
 *  - Src2d_Tgt0d: src is (ncols,nlevs), while tgt is a single value for all cols
 * The latter is mainly used by the FieldAtPressureLevel diagnostic, while the second
 * and third are used when saving/loading data to/from a uniform set of pressure levels
 */


template<int PackSize>
class VerticalInterpolation {
public:
  static constexpr int PS = PackSize;

  using PackT = ekat::Pack<Real,PS>;
  using MaskT = ekat::Mask<PS>;
  using MaskVal = typename MaskT::type;
  using ELI = ekat::LinInterp<Real,PackSize>;

  using KT = KokkosTypes<DefaultDevice>;
  using MemberType = typename KT::MemberType;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;
  template<typename T>
  using view_2d = typename KT::template view_2d<T>;
  template<typename T>
  using view_3d = typename KT::template view_3d<T>;
  template<typename T, int N>
  using view_Nd = typename KT::template view_ND<T,N>;

  enum XType {
    Src2d_Tgt2d,  // Both src and tgt are (ncols,nlevs)
    Src1d_Tgt2d,  // Src is (nlevs), while tgt is (ncols,nlevs)
    Src2d_Tgt1d,  // Src is (ncols,nlevs), while tgt is (nlevs)
    Src2d_Tgt0d,  // Src is (ncols,nlevs), while tgt is a single value for all cols
    Unset         // Used to detect if setup method was not called
  };

  VerticalInterpolation (const int ncols, const int nlevs_src, const int nlevs_tgt);
  VerticalInterpolation (const int ncols, const int nlevs_src, const int nlevs_tgt, const Real mask_val);

  void setup (const view_2d<const Real>& x_src,
              const view_2d<const Real>& x_tgt);
  void setup (const view_2d<const Real>& x_src,
              const view_1d<const Real>& x_tgt);
  void setup (const view_1d<const Real>& x_src,
              const view_2d<const Real>& x_tgt);
  void setup (const view_2d<const Real>& x_src,
              const Real x_tgt);

  void setup (const view_2d<const PackT>& x_src,
              const view_2d<const PackT>& x_tgt);
  void setup (const view_2d<const PackT>& x_src,
              const view_1d<const PackT>& x_tgt);
  void setup (const view_1d<const PackT>& x_src,
              const view_2d<const PackT>& x_tgt);
  void setup (const view_2d<const PackT>& x_src,
              const Real x_tgt);

  void interpolate (const view_2d<const Real>& y_src,
                    const view_2d<      Real>& y_tgt,
                    const view_2d<      MaskVal>& mask = view_2d<MaskVal>());
  void interpolate (const view_3d<const Real>& y_src,
                    const view_3d<      Real>& y_tgt,
                    const view_2d<      MaskVal>& mask = view_2d<MaskVal>());
  void interpolate (const view_2d<const Real>& y_src,
                    const view_1d<      Real>& y_tgt,
                    const view_1d<      MaskVal>& mask = view_1d<MaskVal>());
  void interpolate (const view_3d<const Real>& y_src,
                    const view_2d<      Real>& y_tgt,
                    const view_1d<      MaskVal>& mask = view_1d<MaskVal>());

  void interpolate (const view_3d<const PackT>& y_src,
                    const view_3d<      PackT>& y_tgt,
                    const view_2d<      MaskT>& mask = view_2d<MaskT>());
  void interpolate (const view_2d<const PackT>& y_src,
                    const view_2d<      PackT>& y_tgt,
                    const view_2d<      MaskT>& mask = view_2d<MaskT>());
  void interpolate (const view_2d<const PackT>& y_src,
                    const view_1d<      PackT>& y_tgt,
                    const view_1d<      MaskT>& mask = view_1d<MaskT>());
  void interpolate (const view_3d<const PackT>& y_src,
                    const view_2d<      PackT>& y_tgt,
                    const view_1d<      MaskT>& mask = view_1d<MaskT>());

#ifndef KOKKOS_ENABLE_CUDA
protected:
#endif

  void interpolate_1d (const MemberType&   team,
                       const int icol,
                       const view_1d<const PackT>& x_src,
                       const view_1d<const PackT>& x_tgt,
                       const view_1d<const PackT>& y_src,
                       const view_1d<      PackT>& y_tgt,
                       const view_1d<      MaskT>& mask);

  void interpolate_0d (const MemberType&   team,
                       const int icol,
                       const view_1d<const PackT>& x_src,
                       const view_1d<const PackT>& y_src,
                                           PackT & y_tgt,
                                           MaskT* mask);

  struct TagSameDim {};   // y_src and y_tgt have the same layout
  struct TagVertSlice {}; // y_tgt is a vertical slice (at a single level) of y_src

  KOKKOS_FORCEINLINE_FUNCTION
  void operator () (const TagSameDim&, const MemberType& team);
  KOKKOS_FORCEINLINE_FUNCTION
  void operator () (const TagVertSlice&, const MemberType& team);

protected:

  // Helper functions to check that input views to setup/interpolate have compatible extents
  template<typename ViewT>
  void check_packable (const ViewT& v, const std::string& view_type);
  template<typename ViewT>
  void check_packed_len (const ViewT& v, const std::string& view_type);
  template<typename ViewT>
  void check_ncols (const ViewT& v, const std::string& view_type);

  // Helper functions to do the Real->PackT and long->MaskT casts
        MaskT* mask_ptr (      MaskVal* p) { return reinterpret_cast<      MaskT*>(p); }
        PackT* pack_ptr (      Real*    p) { return reinterpret_cast<      PackT*>(p); }
  const PackT* pack_ptr (const Real*    p) { return reinterpret_cast<const PackT*>(p); }

  XType m_xtype = Unset;

  int m_ncols;
  int m_nlevs_src;
  int m_nlevs_tgt;

  Real m_mask_val;

  // When setup is called, we set 1d or 2d for src and tgt, depending on which overload is called
  view_1d<const PackT> m_x_src_1d;
  view_2d<const PackT> m_x_src_2d;
  view_1d<const PackT> m_x_tgt_1d;
  view_2d<const PackT> m_x_tgt_2d;
  PackT m_x_tgt_0d;

  // When interpolate is called, we set 3d/2d for src 3d/2d/1d for tgt/mask, depending on which overload is called
  view_3d<const PackT> m_y_src_3d;
  view_2d<const PackT> m_y_src_2d;
  view_2d<const PackT> m_y_tgt_3d;
  view_2d<const PackT> m_y_tgt_2d;
  view_1d<const PackT> m_y_tgt_1d;
  view_2d<const MaskT> m_mask_2d;
  view_1d<const MaskT> m_mask_1d;
  
  ELI m_lin_interp;
};

} // namespace scream

#include "scream_vertical_interpolation_impl.hpp"

#endif // SCREAM_VERTICAL_INTERPOLATION_HPP
