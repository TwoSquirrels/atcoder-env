#pragma once

#include <risu/prelude.hpp>
#include <risu/math/prime/factors.hpp>

#include <algorithm>
#include <vector>

template <bool osa_k = false, typename T> auto divisors(T n, bool sorted = true) -> std::vector<T> {
  auto result = std::vector<T>(1, 1);
  for (auto [prime, expo] : factors<osa_k>(n)) {
    const auto result_size = int(result.size());
    auto pow_i = T(1);
    for (auto i = 1; i <= expo; ++i) {
      pow_i *= prime;
      for (auto k = 0; k < result_size; ++k) result.emplace_back(result[k] * pow_i);
    }
  }
  if (sorted) std::ranges::sort(result);
  return result;
}

