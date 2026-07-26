#pragma once

#include <risu/prelude.hpp>
#include <risu/util/traits/is_pair.hpp>
#include <risu/util/traits/is_tuple.hpp>
#include <risu/util/traits/iterable.hpp>
#include <risu/util/traits/ostreamable.hpp>
#include <risu/util/int128.hpp>
#include <risu/util/inf.hpp>
#include <risu/util/typename.hpp>

#include <concepts>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#if __has_include(<atcoder/modint>)
#  include <atcoder/modint>
#endif

template <typename T> auto to_pretty_str(const T &target) -> std::string {
  using namespace std;
  auto str = ""s;
  if constexpr (is_void_v<T>) str += "void"s;
  else if constexpr (is_null_pointer_v<T>) str += "null"s;
  else if constexpr (is_same_v<T, bool>) str += target ? "true"s : "false"s;
  else if constexpr (is_same_v<T, char> || is_same_v<T, char16_t> || is_same_v<T, char32_t> || is_same_v<T, wchar_t>) {
    str += "'"s + target + "'"s;
#if __has_include(<atcoder/modint>)
  } else if constexpr (atcoder::internal::is_modint<T>::value) {
    str += to_string(target.val()) + "(mod)"s;
#endif // <atcoder/modint>
  } else if constexpr (is_arithmetic_v<T>) {
    if constexpr (is_same_v<T, __int128_t>) str += int128_to_str(target);
    else str += to_string(target);
    if constexpr (is_unsigned_v<T>) str += "u"s;
    if constexpr (is_same_v<remove_cv_t<T>, long>) str += "L"s;
    else if constexpr (is_same_v<remove_cv_t<T>, long long>) str += "LL"s;
    else if constexpr (is_same_v<T, __int128_t>) str += "LLL"s;
    if constexpr (std::numeric_limits<T>::is_specialized) {
      if (target >= inf<T>()) str = "inf"s;
      else if constexpr (std::is_signed_v<T>) if (target <= -inf<T>()) str = "-inf"s;
    }
  } else if constexpr (is_pair_v<T>) {
    str += "("s + to_pretty_str(target.first);
    str += ", "s + to_pretty_str(target.second) + ")"s;
  } else if constexpr (is_tuple_v<T>) {
    str += "("s;
    auto separate = false;
    apply([&](const auto &... elems) {
      (([&](const auto &elem) {
        if (separate) str += ", "s;
        str += to_pretty_str(elem);
        separate = true;
      })(elems), ...);
    }, target);
    str += ")"s;
  } else if constexpr (convertible_to<T, string>) {
    str += "\""s + target + "\""s;
  } else if constexpr (is_array_v<T>) {
    str += "["s;
    auto separate = false;
    for (const auto &target_i : target) {
      if (separate) str += ", "s;
      str += to_pretty_str(target_i);
      separate = true;
    }
    str += "]"s;
  } else if constexpr (iterable_v<T>) {
    str += "("s + get_typename<T>(20) + "){"s;
    auto separate = false;
    for (const auto &target_i : target) {
      if (separate) str += ","s;
      str += " "s + to_pretty_str(target_i);
      separate = true;
    }
    if (separate) str += " "s;
    str += "}"s;
  } else if constexpr (ostreamable_v<T>) {
    auto ss = std::ostringstream();
    ss << target;
    str += ss.str();
  } else {
    str += "<"s + get_typename<T>(20);
    str += " ("s + to_string(sizeof(target)) + " byte)>"s;
  }
  return str;
}

