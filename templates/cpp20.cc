#ifndef DEBUG
#  pragma GCC optimize("O3,unroll-loops")
#endif
#if !defined(DEBUG) || defined(__x86_64__) || defined(__i386__)
#  pragma GCC target("avx2,bmi2,popcnt,lzcnt")
#endif

#include <risu/main.hpp>
#include <risu/all.hpp>

#if __has_include(<atcoder/all>)
#  include <atcoder/all>
#endif

#if __has_include(<boost/multiprecision/cpp_int.hpp>)
#  include <boost/multiprecision/cpp_int.hpp>
#endif
#if __has_include(<boost/math/special_functions/prime.hpp>)
#  include <boost/math/special_functions/prime.hpp>
#endif
#if __has_include(<boost/multiprecision/miller_rabin.hpp>)
#  include <boost/multiprecision/miller_rabin.hpp>
#endif

#include <bits/stdc++.h>

using namespace std;
#if __has_include(<atcoder/all>)
using namespace atcoder;
#endif
using namespace risu;

#pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-variable"
constexpr auto infi = inf<int>();
constexpr auto infl = inf<long long>();
constexpr auto infd = inf<double>();
constexpr auto infld = inf<long double>();
const auto YN = array<string, 2>{ "No", "Yes" };
#pragma GCC diagnostic pop

#define p1 first
#define p2 second
#define eb emplace_back
#define ef emplace_front
#define pb pop_back
#define pf pop_front
#define sz(x) (i64(ssize(x)))
#define reps(i, l, r) for ([[maybe_unused]] const auto i : iota_range(l, r))
#define rep(i, n) reps(i, 0, n)
#define rreps(i, l, r) for ([[maybe_unused]] const auto i : iota_range(l, r) | views::reverse)
#define rrep(i, n) rreps(i, 0, n)
#define each(for_able) for (auto &&for_able##_i : (for_able))
#define all(for_able) (begin(for_able)), (end(for_able))

/// answer

//using mint=modint998244353;
//using mint=f80;

auto risu::main() {
  i64 n=scan;
  out(n);
}
