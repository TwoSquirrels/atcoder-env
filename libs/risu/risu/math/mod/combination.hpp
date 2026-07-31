#pragma once

#include <risu/prelude.hpp>
#include <risu/math/mod/factorial.hpp>

#include <numeric>
#include <type_traits>

template <typename T> auto perm(long long n, long long k) -> T {
  if (k < 0 || k > n) return 0;
  if constexpr (not std::is_integral_v<T>) {
    if (n <= 5000000) return fact<T>(n) * fact<T>(n - k, true);
  }
  auto res = T(1);
  for (auto i = n; i > n - k; --i) res *= i;
  return res;
}

template <typename T> auto comb(long long n, long long k) -> T {
  if (k < 0 || k > n) return 0;
  if (n - k < k) k = n - k;
  if constexpr (std::is_integral_v<T>) {
    auto res = T(1);
    for (auto i = 1LL; i <= k; ++i) {
      const auto num = T(n - k + i), den = T(i);
      const auto g = std::gcd(num, den);
      res /= den / g;
      res *= num / g;
    }
    return res;
  }
  return perm<T>(n, k) * fact<T>(k, true);
}
