#ifndef SCREAM_STD_UTILS_HPP
#define SCREAM_STD_UTILS_HPP

#include <algorithm>
#include <memory>

namespace scream {
namespace util {

// This function returns an iterator which is of the same type of c.begin()
template<typename ContainerType, typename T>
auto find (const ContainerType& c, T&& value) -> decltype(c.begin()) {
  return std::find(c.begin(),c.end(),value);
}

template<typename ContainerType, typename T>
bool contains (const ContainerType& c, T&& value) {
  return std::find(c.begin(),c.end(),value) != c.end();
}

} // namespace util

// A set of weak_ptr would not compile, due to the lack of operator<.
// To overcome this, one could add a Compare type to the set template
// arguments. Or, more easily, overload operator<...
// NOTE: this cannot be in the util namespace, or else the compiler won't detect
//       it, when T is declared in the scream namespace

template<typename T>
inline bool operator< (const std::weak_ptr<T>& p, const std::weak_ptr<T>& q)
{
  return p.owner_before(q);
}

} // namespace scream

#endif // SCREAM_STD_UTILS_HPP
