#pragma once

#include <risu/prelude.hpp>
#include <risu/math/prime/prime.hpp>

#include <array>
#include <utility>
#include <vector>
#if __has_include(<boost/math/special_functions/prime.hpp>)
#  include <boost/math/special_functions/prime.hpp>
#endif

// thanks to https://perogram.hateblo.jp/entry/2019/03/29/193632
template <bool osa_k = false, typename T> auto factors(T n) -> std::vector<std::pair<T, int>> {
  constexpr auto spf_limit = 1 << 24;
  auto result = std::vector<std::pair<T, int>>();
  if (n < 0) {
    result.emplace_back(-1, 1);
    const auto natural = factors<osa_k>(-n);
    result.insert(result.end(), natural.begin(), natural.end());
  } else if (n == 0) result.emplace_back(0, 1);
  else if (n == 1);
  else if (osa_k && n < spf_limit) {
    static std::array<int, spf_limit> spf;
    if (spf[0] == 0) {
      for (auto i = 0; i < spf_limit; ++i) spf[i] = i % 2 == 0 ? 2 : i % 3 == 0 ? 3 : i;
      spf[0] = 1;
      for (auto p1 = 6; (p1 - 1) * (p1 - 1) <= spf_limit; p1 += 6) {
        if (spf[p1 - 1] == p1 - 1) for (auto i = p1 - 1; i < spf_limit; i += p1 - 1) if (spf[i] == i) spf[i] = p1 - 1;
        if (spf[p1 + 1] == p1 + 1) for (auto i = p1 + 1; i < spf_limit; i += p1 + 1) if (spf[i] == i) spf[i] = p1 + 1;
      }
    }
    while (n != 1) {
      if (not result.empty() && result.back().first == spf[n]) ++result.back().second;
      else result.emplace_back(spf[n], 1);
      n /= spf[n];
    }
  } else {
    auto expo = 0;
    while (n % 2 == 0) { n /= 2; ++expo; }
    if (expo != 0) result.emplace_back(2, expo);
    expo = 0;
    while (n % 3 == 0) { n /= 3; ++expo; }
    if (expo != 0) result.emplace_back(3, expo);
    auto tried = T(3);
#if __has_include(<boost/math/special_functions/prime.hpp>)
    for (auto i = 2; i <= int(boost::math::max_prime); ++i) {
      const auto prime_i = boost::math::prime(i);
      if (prime_i > n / prime_i) break;
      expo = 0;
      while (n % prime_i == 0) { n /= prime_i; ++expo; }
      result.emplace_back(prime_i, expo);
      tried = prime_i;
    }
    if (osa_k && n < spf_limit) {
      const auto rest = factors<osa_k>(n);
      result.insert(result.end(), rest.begin(), rest.end());
      n = 1;
    }
#endif // <boost/math/special_functions/prime.hpp>
    if (is_prime(n)) { result.emplace_back(n, 1); n = 1; }
    for (auto i = (tried + 5) / 6 * 6; (i - 1) * (i - 1) <= n; i += 6) {
      expo = 0;
      while (n % (i - 1) == 0) { n /= (i - 1); ++expo; }
      if (expo != 0) result.emplace_back(i - 1, expo);
      expo = 0;
      while (n % (i + 1) == 0) { n /= (i + 1); ++expo; }
      if (expo != 0) result.emplace_back(i + 1, expo);
      if (osa_k && n < spf_limit) {
        const auto rest = factors<osa_k>(n);
        result.insert(result.end(), rest.begin(), rest.end());
        n = 1;
      }
    }
    if (n != 1) result.emplace_back(n, 1);
  }
  return result;
}

