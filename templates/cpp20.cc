// SPDX-License-Identifier: MIT
// (c) 2023 TwoSquirrels
// my AtCoder environment: https://github.com/TwoSquirrels/atcoder-env

//#define DEBUG

#ifndef DEBUG
#  pragma GCC optimize("O3,unroll-loops")
#endif // DEBUG

#include <risu/std.hpp>

#if __has_include(<boost/multiprecision/cpp_int.hpp>)
#  include <boost/multiprecision/cpp_int.hpp>
#endif
#if __has_include(<boost/math/special_functions/prime.hpp>)
#  include <boost/math/special_functions/prime.hpp>
#endif
#if __has_include(<boost/multiprecision/miller_rabin.hpp>)
#  include <boost/multiprecision/miller_rabin.hpp>
#endif

#ifndef DEBUG
#  pragma GCC target("avx2,bmi2,popcnt,lzcnt")
#endif // DEBUG

#if __has_include(<atcoder/all>)
#  include <atcoder/all>
#endif

#include <risu/all.hpp>

using namespace std;
namespace rng = std::ranges;
namespace viw = std::ranges::views;
#if __has_include(<atcoder/all>)
using namespace atcoder;
#endif // <atcoder/all>

/// answer

//using mint=modint998244353;

inline auto cp_main() -> str {

  return "";
}
