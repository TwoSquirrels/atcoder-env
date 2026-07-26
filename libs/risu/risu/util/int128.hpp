#pragma once

#include <risu/prelude.hpp>

#include <algorithm>
#include <string>

inline auto int128_to_str(__int128_t target) -> std::string {
  auto target_str = std::string();
  auto target_tmp = __uint128_t(target < 0 ? -target : target);
  do {
    target_str += '0' + target_tmp % 10;
    target_tmp /= 10;
  } while (target_tmp != 0);
  if (target < 0) target_str += '-';
  std::ranges::reverse(target_str);
  return target_str;
}

