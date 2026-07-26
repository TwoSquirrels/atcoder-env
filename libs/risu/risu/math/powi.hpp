#pragma once

#include <risu/prelude.hpp>

inline auto powi(int base, int exponent = 2) -> long long {
  auto ans = 1LL;
  for (auto i = exponent; i != 0; i += (i >= 0 ? -1 : 1)) {
    if (i >= 0) ans *= base;
    else ans /= base;
    if (ans == 0) break;
  }
  return ans;
}

