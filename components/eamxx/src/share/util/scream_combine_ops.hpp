#ifndef SCREAM_COMBINE_OPS_HPP
#define SCREAM_COMBINE_OPS_HPP

#include "share/util/scream_universal_constants.hpp"

#include <ekat/ekat_scalar_traits.hpp>
#include <ekat/util/ekat_math_utils.hpp>

#include <type_traits>

namespace scream {

/*
 * Flags used to specify how to handle a variable update.
 * When computing f(x), there are a few ways to combine
 * the result with the value of the variable were we want
 * to store it:
 *    y = alpha*f(x) + beta*y
 *    y = y*f(x)
 *    y = y/f(x)
 * This enum can be used as template arg in some general functions,
 * so that we can write a single f(x), and then combine:
 *    combine<CM>(f(x),y,alpha,beta);
 * The result has zero overhead compared to any specific version,
 * since the if/switch statements involving the combine mode CM
 * are compiled-away by the compiler
 */

enum class CombineMode {
  ScaleUpdate,  // out = beta*out + alpha*in (most generic case)
  Update,       // out = beta*out + in (special case of ScaleUpdate wiht alpha=1)
  ScaleAdd,     // out = out + alpha*in (special case of ScaleUpdate with beta=1)
  ScaleReplace, // out = alpha*in (special case of ScaleUpdate with beta=0)
  Add,          // out = out + in (special case of ScaleUpdate with alpha=beta=1)
  Rescale,      // out = beta*out
  Replace,      // out = in
  Max,          // out = max(in,out)
  Min,          // out = min(in,out)
  Multiply,     // out = out*in
  Divide        // out = out/in
};

// Small helper functions to combine a new value with an old one.
// The template argument help reducing the number of operations
// performed (the switch is resolved at compile time). In the most
// complete form, the function performs
//    result = beta*result + alpha*newVal
// This routine should have no overhead compared to a manual
// update (assuming you call it with the proper CM)
// If UseFill is true, we check newVal and result before deciding how to update.
// A "fill value" is to be interpreted as an "uninitialized" number. So, for instance
//   max(A,fillValue) = A
//   A+fillValue = fillValue

/* Special version of combine that takes a mask into account */
template<CombineMode CM, typename ScalarIn, typename ScalarOut,
         typename CoeffType, bool UseFill>
KOKKOS_FORCEINLINE_FUNCTION
void combine_impl (const ScalarIn& newVal, ScalarOut& result,
                   const CoeffType alpha, const CoeffType beta, const ScalarOut& fill_v = ScalarOut(0));
{
#ifndef NDEBUG
  EKAT_KERNEL_REQUIRE_MSG (needsAlpha<CM>() or alpha==CoeffType(1),
                           "Error! Passing alpha!=1 with a combine mode that does not require alpha.\n");
  EKAT_KERNEL_REQUIRE_MSG (needsBeta<CM>() or beta==CoeffType(0),
                           "Error! Passing beta!=0 with a combine mode that does not require beta.\n");
#endif
  switch (CM) {
    case CombineMode::Replace:
      if constexpr (UseFill) {
        if (newVal!=fill_v) result = newVal;
      } else {
        result = newVal;
      }
      break;
    case CombineMode::Min:
      if constexpr (UseFill) {
        if (newVal!=fill_v)
          if (result!=fill_v)
            result = ekat::impl::min(newVal,result);
          else
            result = newVal;
      } else {
        result = ekat::impl::min(newVal,result);
      }
      break;
    case CombineMode::Max:
      if constexpr (UseFill) {
        if (newVal!=fill_v)
          if (result!=fill_v)
            result = ekat::impl::min(newVal,result);
          else
            result = newVal;
      } else {
        result = ekat::impl::max(newVal,result);
      }
      break;
    case CombineMode::Rescale:
      if constexpr (UseFill) {
        if (result != fill_v) {
          result *= beta;
      } else {
        result *= beta;
      }
      break;
    case CombineMode::ScaleReplace:
      if constexpr (UseFill) {
        if (newVal == fill_v) {
          result = fill_v;
        } else {
          result = alpha*newVal;
        }
      } else {
        result = alpha*newVal;
      }
      break;
    case CombineMode::Update:
      if constexpr (UseFill) {
        if (result == fill_v || newVal == fill_v) {
          result = fill_v;
        } else {
          result *= beta;
          result += newVal;
        }
      } else {
        result *= beta;
        result += newVal;
      }
      break;
    case CombineMode::ScaleUpdate:
      if constexpr (UseFill) {
        if (result == fill_v || newVal == fill_v) {
          result = fill_v;
        } else {
          result *= beta;
          result += alpha*newVal;
        }
      } else {
        result *= beta;
        result += alpha*newVal;
      }
      break;
    case CombineMode::ScaleAdd:
      if constexpr (UseFill) {
        if (result == fill_v || newVal == fill_v) {
          result = fill_v;
        } else {
          result += alpha*newVal;
        }
      } else {
        result += alpha*newVal;
      }
      break;
    case CombineMode::Add:
      if constexpr (UseFill) {
        if (result == fill_v || newVal == fill_v) {
          result = fill_v;
        } else {
          result += newVal;
        }
      } else {
        result += newVal;
      }
      break;
    case CombineMode::Multiply:
      if constexpr (UseFill) {
        if (result == fill_v || newVal == fill_v) {
          result = fill_v;
        } else {
          result *= newVal;
        }
      } else {
        result *= newVal;
      }
      break;
    case CombineMode::Divide:
      if constexpr (UseFill) {
        if (result == fill_v || newVal == fill_v) {
          result = fill_v;
        } else {
          result /= newVal;
        }
      } else {
        result /= newVal;
      }
      break;
      break;
    default:
      EKAT_KERNEL_ERROR_MSG ("Error! Unrecognized/unsupported CombineMode.\n");
  }
}

template<CombineMode CM, typename ScalarIn, typename ScalarOut,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
void combine (const ScalarIn& newVal, ScalarOut& result,
              const CoeffType alpha, const CoeffType beta)
{
  combine_impl<CM,ScalarIn,ScalarOut,CoeffType,false>(newVal,result,alpha,beta);
}

template<CombineMode CM, typename ScalarIn, typename ScalarOut,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
void combine (const ScalarIn& newVal, ScalarOut& result, const ScalarOut& fill_val
              const CoeffType alpha, const CoeffType beta)
{
  combine_impl<CM,ScalarIn,ScalarOut,CoeffType,true>(newVal,result,alpha,beta,fill_val);
}

} // namespace scream

#endif // SCREAM_COMBINE_OPS_HPP
