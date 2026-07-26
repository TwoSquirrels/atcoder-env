#pragma once

#include <risu/prelude.hpp>
#include <risu/util/traits/is_pair.hpp>
#include <risu/util/traits/iterable.hpp>
#include <risu/util/traits/ostreamable.hpp>
#include <risu/util/int128.hpp>

#include <concepts>
#include <iostream>
#if __has_include(<atcoder/modint>)
#  include <atcoder/modint>
#endif

template <typename T> inline auto write_stdout(const T &target, bool flush = false) -> void {
  std::cout << target;
  if (flush) std::cout << std::flush;
}
template <typename T, typename Sep = char> inline auto output(const T &target, Sep separator = ' ', bool flush = false) -> void {
  if constexpr (ostreamable_v<T>) {
    write_stdout(target, flush);
  } else if constexpr (std::convertible_to<T, __int128_t>) {
    write_stdout(int128_to_str(target), flush);
#  if __has_include(<atcoder/modint>)
  } else if constexpr (atcoder::internal::is_modint<T>::value) {
    output(target.val(), separator, flush);
#  endif // <atcoder/modint>
  } else if constexpr (iterable_v<T>) {
    auto separate = false;
    for (const auto &target_i : target) {
      if (separate) write_stdout(separator);
      output(target_i, separator);
      separate = true;
    }
    if (flush) write_stdout("", flush);
  } else if constexpr (is_pair_v<T>) {
    output(target.first, separator);
    write_stdout(separator);
    output(target.second, separator, flush);
  } else {
    write_stdout("<unknown>", flush);
  }
}
template <typename T, typename Sep = char> inline auto outputln(const T &target, Sep separator = ' ', bool flush = false) -> void {
  output(target, separator);
  write_stdout('\n', flush);
}

