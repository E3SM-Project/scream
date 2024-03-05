#ifndef SCREAM_BACKWARD_VIEW_1D
#define SCREAM_BACKWARD_VIEW_1D

#include <Kokkos_Core.hpp>

namespace scream
{

template<typename ViewT>
class BwdView1d : public ViewT
{
public:
  static_assert (Kokkos::is_view<ViewT>::value,
                 "Error! Template parameter is not a Kokkos view.");
  static_assert (ViewT::rank==1,
                 "Error! Template parameter ViewT must have rank 1.");

  using value_type = typename ViewT::traits::value_type;

  KOKKOS_INLINE_FUNCTION
  BwdView1d () = default;
  KOKKOS_INLINE_FUNCTION
  BwdView1d (const ViewT& v) : ViewT(v) {}
  KOKKOS_INLINE_FUNCTION
  BwdView1d (const BwdView1d&) = default;

  KOKKOS_INLINE_FUNCTION
  ~BwdView1d () = default;

  KOKKOS_INLINE_FUNCTION
  BwdView1d& operator= (const BwdView1d& src) = default;

  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
  value_type& operator()(I0 i0) {
    return ViewT::operator()(ViewT::extent(0)-i0-1);
  }
  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
  value_type& operator[](I0 i0) {
    return ViewT::operator[](ViewT::extent(0)-i0-1);
  }

  BwdView1d<typename ViewT::HostMirror> host_mirror() const {
    auto fwd_h = Kokkos::create_mirror_view(*this);
    BwdView1d<typename ViewT::HostMirror> bwd_h(fwd_h);
    return BwdView1d<typename ViewT::HostMirror> (fwd_h);
  }
};

} // namespace scream

#endif // SCREAM_BACKWARD_VIEW_1D
