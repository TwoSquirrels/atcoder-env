#pragma once

// TODO: Opus 5 にリファクタさせたら複雑化しすぎたのでシンプルにする.

#include "../util/traits/is_pair.hpp"
#include "../util/traits/is_tuple.hpp"
#include "../util/traits/iterable.hpp"
#include "../util/traits/ostreamable.hpp"
#include "../util/int128.hpp"
#include "../util/inf.hpp"

#if __has_include(<atcoder/modint>)
#  include <atcoder/modint>
#endif

#if __has_include(<cxxabi.h>)
#  include <cxxabi.h>
#endif

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <variant>

#if __has_include(<format>)
#  include <format>
#endif

namespace risu {

inline std::size_t prettify_max_elements = 240;
inline std::size_t prettify_max_depth = 16;
inline std::size_t prettify_max_typename = 20;

template <typename T, typename = void> struct prettifier;

template <typename T> auto prettify_to(std::string &out, const T &target, std::size_t depth = 0) -> void;
template <typename T> auto prettify(const T &target) -> std::string;

namespace internal {

using namespace std;

template <typename T, typename = void> constexpr bool prettifier_defined_v = false;
template <typename T> constexpr bool prettifier_defined_v<T, void_t<decltype(prettifier<T>::apply(declval<const T &>()))>> = true;

template <typename T, typename = void> constexpr bool formattable_v = false;
#ifdef __cpp_lib_format
template <typename T> constexpr bool formattable_v<T, void_t<decltype(declval<formatter<T, char> &>().format(declval<const T &>(), declval<format_context &>()))>> = true;
#endif

// TODO: traits に移動.

template <typename T> constexpr bool char_like_v = is_same_v<T, char> || is_same_v<T, char8_t> || is_same_v<T, char16_t> || is_same_v<T, char32_t> || is_same_v<T, wchar_t>;

template <typename T> struct string_like : false_type {};
template <typename C, typename Tr, typename A> struct string_like<basic_string<C, Tr, A>> : true_type {};
template <typename C, typename Tr> struct string_like<basic_string_view<C, Tr>> : true_type {};
template <typename T> constexpr bool string_like_v = string_like<T>::value;

template <typename T> constexpr bool char_pointer_v = is_pointer_v<T> && char_like_v<remove_cv_t<remove_pointer_t<T>>>;
template <typename T> constexpr bool char_array_v = is_array_v<T> && char_like_v<remove_cv_t<remove_extent_t<T>>>;

template <typename T> struct is_optional : false_type {};
template <typename T> struct is_optional<optional<T>> : true_type {};
template <typename T> constexpr bool is_optional_v = is_optional<T>::value;

template <typename T> struct is_variant : false_type {};
template <typename... Ts> struct is_variant<variant<Ts...>> : true_type {};
template <typename T> constexpr bool is_variant_v = is_variant<T>::value;

template <typename T> struct is_smart_pointer : false_type {};
template <typename T, typename D> struct is_smart_pointer<unique_ptr<T, D>> : true_type {};
template <typename T> struct is_smart_pointer<shared_ptr<T>> : true_type {};
template <typename T> constexpr bool is_smart_pointer_v = is_smart_pointer<T>::value;

template <typename T>
concept map_like = requires {
  typename T::key_type;
  typename T::mapped_type;
};

template <typename T>
concept adapter_like = !iterable_v<T> && requires {
  typename T::container_type;
  typename T::value_type;
};

template <adapter_like A> auto underlying(const A &a) -> const typename A::container_type & {
  struct hack : A {
    static auto get(const A &a) -> const typename A::container_type & { return a.*&hack::c; }
  };
  return hack::get(a);
}

template <typename T> auto raw_typename() -> string {
#ifdef __has_include(<cxxabi.h>)
  auto status = 0;
  const auto free_ = [](char *p) { free(p); };
  const auto demangled = unique_ptr<char, decltype(free_)>(abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status), free_);
  if (status == 0 && demangled) return demangled.get();
#endif
  return typeid(T).name();
}

inline auto simplify_typename(string s) -> string {
  constexpr string_view inline_namespaces[] = {"__cxx11::", "__1::", "__ndk1::"};
  for (const auto &ns : inline_namespaces) {
    for (auto pos = s.find(ns); pos != string::npos; pos = s.find(ns, pos)) s.erase(pos, ns.size());
  }

  constexpr string_view defaults[] = {
    "std::allocator<", "std::less<", "std::greater<", "std::hash<", "std::equal_to<",   "std::char_traits<"
  };
  for (const auto &d : defaults) {
    for (auto from = size_t(0);;) {
      const auto pos = s.find(d, from);
      if (pos == string::npos) break;
      auto begin = pos;
      while (begin > 0 && s[begin - 1] == ' ') --begin;
      if (begin == 0 || s[begin - 1] != ',') { from = pos + 1; continue; }
      --begin;
      auto i = pos + d.size();
      for (auto level = 1; i < s.size() && level > 0; ++i) {
        if (s[i] == '<') ++level;
        else if (s[i] == '>') --level;
      }
      s.erase(begin, i - begin);
      from = begin;
    }
  }

  for (auto pos = s.find(" >"); pos != string::npos; pos = s.find(" >", pos))
    s.erase(pos, 1);

  constexpr pair<string_view, string_view> aliases[] = {
    {"std::basic_string_view<char>", "std::string_view"},
    {"std::basic_string<char>", "std::string"},
  };
  for (const auto &[from, to] : aliases) {
    for (auto pos = s.find(from); pos != string::npos; pos = s.find(from, pos + to.size())) {
      s.replace(pos, from.size(), to);
    }
  }

  return s;
}

template <typename T> auto type_label() -> string {
  auto s = simplify_typename(raw_typename<T>());
  const auto limit = prettify_max_typename;
  if (limit != 0 && limit != string::npos && s.size() > limit)
    s = s.substr(0, limit > 3 ? limit - 3 : 0) + "...";
  return s;
}

inline auto append_hex(string &out, uint_least32_t v) -> void {
  constexpr auto digits = string_view("0123456789abcdef");
  auto buf = string();
  do {
    buf += digits[v & 0xf];
    v >>= 4;
  } while (v);
  out.append(buf.rbegin(), buf.rend());
}

inline auto append_escaped(string &out, char32_t c, char quote, bool raw_high) -> void {
  switch (c) {
    case U'\\': out += "\\\\"; return;
    case U'\n': out += "\\n"; return;
    case U'\r': out += "\\r"; return;
    case U'\t': out += "\\t"; return;
    case U'\0': out += "\\0"; return;
    default: break;
  }
  if (c == static_cast<char32_t>(static_cast<unsigned char>(quote))) {
    out += '\\';
    out += quote;
    return;
  }
  if (0x20 <= c && c != 0x7f && (c < 0x80 || raw_high)) {
    out += static_cast<char>(c);
    return;
  }
  out += "\\u{";
  append_hex(out, static_cast<uint_least32_t>(c));
  out += '}';
}

template <typename C> auto append_quoted(string &out, const C *p, size_t n) -> void {
  constexpr auto raw_high = sizeof(C) == 1;
  out += '"';
  for (auto i = size_t(0); i < n; ++i) {
    append_escaped(out, static_cast<char32_t>(static_cast<make_unsigned_t<C>>(p[i])), '"', raw_high);
  }
  out += '"';
}

template <typename T> auto append_number(string &out, const T &target) -> void {
  if constexpr (is_same_v<T, __int128_t> || is_same_v<T, __uint128_t>) {
    out += int128_to_str(target);
  } else {
#ifdef __cpp_lib_format
    out += format("{}", target);
#else
    out += to_string(target);
#endif
  }
}

template <typename T> auto append_suffix(string &out) -> void {
  if constexpr (is_floating_point_v<T>) {
    if constexpr (is_same_v<T, float>) out += "f";
    else if constexpr (is_same_v<T, long double>) out += "L";
  } else {
    if constexpr (is_unsigned_v<T>) out += "u";
    using S = make_signed_t<T>;
    if constexpr (is_same_v<S, long>) out += "L";
    else if constexpr (is_same_v<S, long long>) out += "LL";
    else if constexpr (is_same_v<S, __int128_t>) out += "LLL";
  }
}

template <typename R> auto append_elements(string &out, const R &range, string_view sep, size_t depth) -> void {
  auto count = size_t(0);
  auto separate = false;
  for (const auto &elem : range) {
    if (count++ == prettify_max_elements) {
      out += sep;
      out += "...";
      if constexpr (requires { size(range); }) {
        const auto total = static_cast<size_t>(size(range));
        if (total >= count) out += " (" + to_string(total - count + 1) + " more)";
      }
      break;
    }
    if (separate) out += sep;
    prettify_to(out, elem, depth + 1);
    separate = true;
  }
}

}

#if __has_include(<atcoder/modint>)
template <typename T>
struct prettifier<T, std::enable_if_t<atcoder::internal::is_modint<T>::value>> {
  static auto apply(const T &target) -> std::string {
    return std::to_string(target.val()) + "(mod " + std::to_string(T::mod()) + ")";
  }
};
#endif

template <typename T>
auto prettify_to(std::string &out, const T &target, std::size_t depth) -> void {
  using namespace std;
  namespace pi = internal;
  using U = remove_cv_t<T>;

  if (depth > prettify_max_depth) { out += "..."; return; }

  if constexpr (pi::prettifier_defined_v<U>) {
    out += prettifier<U>::apply(target);
  }
  else if constexpr (is_null_pointer_v<U>) {
    out += "null";
  } else if constexpr (is_same_v<U, bool>) {
    out += target ? "true" : "false";
  } else if constexpr (pi::char_like_v<U>) {
    out += '\'';
    pi::append_escaped(out, static_cast<char32_t>(static_cast<make_unsigned_t<U>>(target)),
                       '\'', sizeof(U) == 1);
    out += '\'';
  } else if constexpr (is_arithmetic_v<U> || is_same_v<U, __int128_t> || is_same_v<U, __uint128_t>) {
    if constexpr (is_floating_point_v<U>) {
      if (target != target) { out += "nan"; return; }
    }
    if constexpr (numeric_limits<U>::is_specialized) {
      if (target >= inf<U>()) { out += "inf"; return; }
      if constexpr (is_signed_v<U>) {
        if (target <= -inf<U>()) { out += "-inf"; return; }
      }
    }
    const auto head = out.size();
    pi::append_number(out, target);
    if constexpr (is_floating_point_v<U>) {
      if (out.find_first_of(".eE", head) == string::npos) out += ".0";
    }
    pi::append_suffix<U>(out);
  }
  else if constexpr (pi::string_like_v<U>) {
    pi::append_quoted(out, target.data(), target.size());
  } else if constexpr (pi::char_array_v<U>) {
    auto n = extent_v<U>;
    if (n > 0 && target[n - 1] == 0) --n;
    pi::append_quoted(out, static_cast<const remove_extent_t<U> *>(target), n);
  } else if constexpr (pi::char_pointer_v<U>) {
    if (target == nullptr) out += "null";
    else {
      auto n = size_t(0);
      while (target[n] != 0) ++n;
      pi::append_quoted(out, target, n);
    }
  }
  else if constexpr (pi::is_smart_pointer_v<U>) {
    if (!target) out += "null";
    else { out += "*"; prettify_to(out, *target, depth + 1); }
  } else if constexpr (is_pointer_v<U>) {
    if (target == nullptr) out += "null";
    else {
      out += "(" + pi::type_label<U>() + ")0x";
      pi::append_hex(out, static_cast<uintptr_t>(reinterpret_cast<uintptr_t>(target)));
    }
  }
  else if constexpr (pi::is_optional_v<U>) {
    if (!target) out += "none";
    else { out += "some("; prettify_to(out, *target, depth + 1); out += ")"; }
  } else if constexpr (is_same_v<U, monostate>) {
    out += "monostate";
  } else if constexpr (pi::is_variant_v<U>) {
    out += "<" + to_string(target.index()) + ">";
    visit([&](const auto &alt) { prettify_to(out, alt, depth + 1); }, target);
  }
  else if constexpr (is_pair_v<U>) {
    out += "(";
    prettify_to(out, target.first, depth + 1);
    out += ", ";
    prettify_to(out, target.second, depth + 1);
    out += ")";
  } else if constexpr (is_tuple_v<U>) {
    out += "(";
    auto separate = false;
    apply([&](const auto &...elems) {
      (([&](const auto &elem) {
         if (separate) out += ", ";
         prettify_to(out, elem, depth + 1);
         separate = true;
       })(elems), ...);
    }, target);
    out += ")";
  }
  else if constexpr (is_array_v<U>) {
    out += "[";
    pi::append_elements(out, target, ", ", depth);
    out += "]";
  } else if constexpr (pi::map_like<U>) {
    out += "(" + pi::type_label<U>() + "){";
    auto count = size_t(0);
    auto separate = false;
    for (const auto &[key, value] : target) {
      if (count++ == prettify_max_elements) { out += ", ..."; break; }
      out += separate ? ", " : " ";
      prettify_to(out, key, depth + 1);
      out += ": ";
      prettify_to(out, value, depth + 1);
      separate = true;
    }
    out += separate ? " }" : "}";
  } else if constexpr (iterable_v<U>) {
    out += "(" + pi::type_label<U>() + "){";
    const auto head = out.size();
    pi::append_elements(out, target, ", ", depth);
    if (out.size() != head) { out.insert(head, " "); out += " "; }
    out += "}";
  } else if constexpr (pi::adapter_like<U>) {
    out += "(" + pi::type_label<U>() + "){";
    const auto &c = pi::underlying(target);
    const auto head = out.size();
    pi::append_elements(out, c, ", ", depth);
    if (out.size() != head) { out.insert(head, " "); out += " "; }
    out += "}";
  }
  else if constexpr (pi::formattable_v<U>) {
#ifdef __cpp_lib_format
    out += format("{}", target);
#endif
  } else if constexpr (ostreamable_v<U>) {
    auto ss = ostringstream();
    ss << target;
    out += ss.str();
  } else {
    out += "<" + pi::type_label<U>();
    out += " (" + to_string(sizeof(target)) + " byte)>";
  }
}

template <typename T> auto prettify(const T &target) -> std::string {
  auto str = std::string();
  str.reserve(64);
  prettify_to(str, target);
  return str;
}

}
