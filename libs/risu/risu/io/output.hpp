#pragma once

#if __has_include(<atcoder/modint>)
#  include <atcoder/modint>
#endif

#include <algorithm>
#include <array>
#include <concepts>
#include <ranges>
#include <span>
#include <string_view>
#include <tuple>
#include <utility>

#include <cstdint>
#include <cstdio>
#include <cstring>

namespace risu {

inline constexpr auto LUT4 = ([]() {
  auto lut4 = std::array<char, 40000>{};
  for (const auto i : std::views::iota(0, 10000)) {
    lut4[i * 4 + 0] = static_cast<char>('0' + i / 1000);
    lut4[i * 4 + 1] = static_cast<char>('0' + i / 100 % 10);
    lut4[i * 4 + 2] = static_cast<char>('0' + i / 10 % 10);
    lut4[i * 4 + 3] = static_cast<char>('0' + i % 10);
  }
  return lut4;
})();

// TODO: int128 対応.
class Printer {
public:
  Printer(std::FILE *fp) : fp(fp), pos(0u), isLineEnd(false), isUnbuffered(false) {}
  ~Printer() { operator++(0); }
  auto operator()(const auto &...value) -> Printer & {
    if (std::exchange(isLineEnd, false)) write('\n');
    print(std::tie(value...));
    isLineEnd = true;
    return *this;
  }
  auto operator--(int) -> Printer & {
    isLineEnd = false;
    return *this;
  }
  auto operator++(int) -> void {
    if (std::exchange(isLineEnd, false)) write('\n');
    flush();
  }
private:
  std::FILE *fp;
  std::array<char, 1u << 18> buf;
  std::size_t pos;
  bool isLineEnd;
  bool isUnbuffered;
  auto flush() -> void {
    if (pos == 0u) return; // 不使用時の副作用を防ぐため.
    if (not std::exchange(isUnbuffered, true)) std::setvbuf(fp, nullptr, _IONBF, 0u);
    std::fwrite(buf.data(), sizeof(decltype(buf)::value_type), pos, fp);
    pos = 0u;
  }
  auto write(char c) -> void {
    if (pos == buf.size()) [[unlikely]] flush();
    buf[pos++] = c;
  }
  auto write(std::string_view s) -> void {
    if (pos + s.size() > buf.size()) [[unlikely]] {
      flush();
      if (s.size() >= buf.size()) {
        std::fwrite(s.data(), sizeof(decltype(s)::value_type), s.size(), fp);
        return;
      }
    }
    std::ranges::copy(s, std::span(buf).subspan(pos).begin());
    pos += s.size();
  }
  template <typename UInt> requires (std::unsigned_integral<UInt> && not std::same_as<UInt, char>) auto write(UInt n) {
    // TODO
  }
  // TODO: 整数, コンテナ, / 演算子対応.
  template <typename T> auto print(const T &value) -> void {
    using namespace std;
    if constexpr (is_same_v<T, char>) write(value);
    else if constexpr (is_same_v<T, bool>) write(value ? "Yes" : "No");
    else if constexpr (convertible_to<const T &, string_view>) write(string_view(x));
    else if constexpr (is_floating_point_v<T>) {
      auto t = array<char, 64>();
      const auto [p, ec] = to_chars(t.data(), t.data() + t.size(), value, chars_format::fixed, 15);
      write(string_view(t.data(), p - t.data()));
    } else static_assert(sizeof(T) == 0, "risu::Printer cannot print T.");
  }
};

inline auto out = Printer{stdout};

/*
template <typename T, typename Sep = char> auto output(const T &target, Sep
separator = ' ', bool flush = false) -> void {
  // TODO: 限界高速化.
  using namespace std;
  if constexpr (ostreamable_v<T>) {
    cout << target;
    if (flush) cout << flush;
    return;
  }
  if constexpr (convertible_to<T, __int128_t>) {
    output(int128_to_str(target), flush);
#if __has_include(<atcoder/modint>)
  } else if constexpr (atcoder::internal::is_modint<T>::value) {
    output(target.val(), separator, flush);
#endif
  } else if constexpr (iterable_v<T>) {
    auto separate = false;
    for (const auto &target_i : target) {
      if (separate) output(separator);
      output(target_i, separator);
      separate = true;
    }
    if (flush) output("", flush);
  } else if constexpr (is_pair_v<T>) {
    output(target.first, separator);
    output(separator);
    output(target.second, separator, flush);
  } else if constexpr (is_tuple_v<T>) {
    auto separate = false;
    std::apply([&](const auto &... elems) {
      (([&](const auto &elem) {
        if (separate) output(separator);
        output(elem, separator);
        separate = true;
      })(elems), ...);
    }, target);
    if (flush) output("", flush);
  } else {
    output("<unknown>", flush);
  }
}
template <typename T, typename Sep = char> inline auto outputln(const T &target,
Sep separator = ' ', bool flush = false) -> void { output(target, separator);
  output('\n', separator, flush);
}
*/

}
