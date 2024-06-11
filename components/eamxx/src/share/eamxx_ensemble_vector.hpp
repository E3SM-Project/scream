#ifndef EAMXX_ENSEMBLE_VECTOR_HPP
#define EAMXX_ENSEMBLE_VECTOR_HPP

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_math.hpp"

#include "scream_config.h"

namespace scream
{

// On device memory, an ensample is a proxy for a strided view
// The reason is that if the field dims are, say, (ncol,nlev),
// to keep the same loops as in the scalar case (and keep
// coalesced access on GPU), we store as (ncol,ens_size,nlev)
// so that an ensemble is a subview(icol,ALL,ilev).
template<typename T, typename Device>
struct Ensemble;

template<typename T>
struct EnsembleView : public Unmanaged<KokkosTypes::view_1d<T>>
{
  static constexpr int N = SCREAM_ENSEMBLE_SIZE;
  using base = Unmanaged<KokkosTypes::view_1d<T>>;

  EnsembleRef (const base& src) : base(src) {}

  KOKKOS_INLINE_FUNCTION
  T mean () const {
    T m = 0;
    for (int i=0; i<N; ++i) {
      m += this->operator()(i);
    }
    return m;
  }

  KOKKOS_INLINE_FUNCTION
  T std () const {
    if constexpr (N==1) {
      return 0;
    } else {
      T s = 0;
      auto m = mean();
      for (int i=0; i<N; ++i) {
        s += pow(this->operator()(i) - m,2)
      }
      return s / (N-1);
    }
  }
};

// TODO: add overloads of op (+-*/), as well as op= and math functions

template<typename T>
struct EnsemblePack<T> : public ekat::Pack<T,SCREAM_ENSEMBLE_SIZE>
{
  static constexpr int N = SCREAM_ENSEMBLE_SIZE;
  using base = ekat::Pack<T,N>;
  using base::Pack;

  KOKKOS_INLINE_FUNCTION
  T mean () const {
    return reduce_sum(*this);
  }

  KOKKOS_INLINE_FUNCTION
  T std () const {
    if constexpr (N==1) {
      return 0;
    } else {
      base diff  = *this - base(mean());
      return reduce_sum(diff*diff) / (N-1);
    }
  }
};

template<typename T, typename DeviceT>
struct EnsembleSwitch
{
  using type = typename std::conditional
                    <
                      ekat::OnGpu<DeviceT::execution_space>::value,
                      EnsembleView<T>,
                      EnsemblePack<T>
                    >::type;
};

// template<typename T, typename DeviceT>
// using Ensemble = typename EnsembleSwitch<T,DeviceT>::type;


template<typename T, typename DeviceT>
struct RefType
{
  using type = typename std::conditional
                    <
                      ekat::OnGpu<DeviceT::execution_space>::value,
                      EnsembleView<T>,
                      EnsemblePack<T>&
                    >::type;
};

template<typename... Props>
class View<En

} // namespace scream

#endif // EAMXX_ENSEMBLE_VECTOR_HPP
