#pragma once

// TODO: remove (適切な場所へ)

#include "risu/prelude.hpp"

#if __has_include(<boost/multiprecision/cpp_int.hpp>)
#  include <boost/multiprecision/cpp_int.hpp>
#endif

#include <limits>

#ifdef __SIZEOF_INT128__
using i128 = __int128;
using u128 = unsigned __int128;
static_assert(sizeof(i128) == 16 && sizeof(u128) == 16, "i128/u128 is not 128 bits");
static_assert(i128(-1) < i128(0) && u128(-1) > u128(0), "i128/u128 signedness is broken");
#else
#  error "__int128 is not available on this target"
#endif

using f80 = long double; // 80 bits or wider
static_assert(std::numeric_limits<f80>::digits >= 64, "f80 has no extended precision");

#if __has_include(<boost/multiprecision/cpp_int.hpp>)
using bigint = boost::multiprecision::cpp_int;
#endif
