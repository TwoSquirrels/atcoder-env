#pragma once

#ifdef DEBUG

#include "prettify.hpp"

#include <chrono>
#include <iostream>
#include <source_location>
#include <string>

namespace risu {
namespace internal {

inline auto debug_total = std::chrono::steady_clock::duration{0};
template <typename F> auto debug_txt_f(F callback, const std::source_location &loc = std::source_location::current()) -> void {
  using namespace std;
  using namespace std::chrono;
  const auto start = steady_clock::now();
  auto msg = "[DEBUG] ("s + loc.file_name() + ":L"s + to_string(loc.line()) + ") "s;
  msg += callback();
  cerr << msg << endl;
  const auto end = steady_clock::now();
  debug_total += end - start;
}

}
}

#  define debug_txt(txt) (::risu::internal::debug_txt_f([&]() { return ::risu::prettify(txt); }))

#else
#  define debug_txt(...)
#endif
