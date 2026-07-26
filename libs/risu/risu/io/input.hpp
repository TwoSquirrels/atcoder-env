#pragma once

#include <risu/prelude.hpp>
#include <risu/util/traits/is_pair.hpp>
#include <risu/util/traits/istreamable.hpp>
#include <risu/util/traits/iterable.hpp>

#include <iostream>
#include <type_traits>
#ifdef DEBUG
#  include <risu/util/typename.hpp>
#  include <chrono>
#endif // DEBUG

#ifdef DEBUG
inline auto input_total = std::chrono::steady_clock::duration{0};
#endif // DEBUG
template <typename T> inline auto read_stdin(T &&target) -> void {
#ifdef DEBUG
  using namespace std::chrono;
  const auto start = steady_clock::now();
  std::cin >> target;
  const auto end = steady_clock::now();
  input_total += end - start;
#else // DEBUG
  std::cin >> target;
#endif // DEBUG
}
template <typename T> inline auto input(T &&target) -> T {
  using T_V = std::remove_reference_t<T>;
  if constexpr (istreamable_v<T_V>) read_stdin(target);
  else if constexpr (iterable_v<T_V>) for (auto &&target_i : target) input(target_i);
  else if constexpr (is_pair_v<T_V>) {
    input(target.first);
    input(target.second);
  } else if constexpr (std::is_convertible_v<long long, T_V>) {
    auto n = 0LL;
    target = input(n);
  } else {
#ifdef DEBUG
    std::cerr << "[WARN] Type " << get_typename<T_V>()
              << " is invalid, so input is skipped;" << std::endl;
#endif //DEBUG
  }
  return target;
}
// input and initialize
struct Scanner { template <typename T> inline operator T() const { T target; return input(target); } };
inline Scanner scan;

