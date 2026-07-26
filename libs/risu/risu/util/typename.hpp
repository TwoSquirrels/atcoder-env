#pragma once

#include <risu/prelude.hpp>

#include <cstddef>
#include <string>
#include <string_view>
#include <typeinfo>

#if __has_include(<cxxabi.h>)
#  define INCLUDED_CXXABI
#  include <cxxabi.h>
#endif

template <typename T> inline auto get_typename(std::size_t length_limit = std::string::npos) -> std::string {
  auto name = std::string();
#ifdef INCLUDED_CXXABI
  name = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, nullptr);
#else // INCLUDED_CXXABI
  name = typeid(T).name();
#endif // INCLUDED_CXXABI
#ifdef __ANDROID__
  constexpr auto ndk_ns = std::string_view("::__ndk1::");
  for (auto pos = name.find(ndk_ns); pos != std::string::npos; pos = name.find(ndk_ns, pos)) {
    name.replace(pos, ndk_ns.size(), "::");
  }
#endif // __ANDROID__
  if (name.length() > length_limit) name = name.substr(0, length_limit - 3) + "...";
  return name;
}

