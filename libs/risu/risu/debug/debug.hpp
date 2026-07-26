#pragma once

#include <risu/prelude.hpp>

#ifdef DEBUG

#include <risu/debug/pretty.hpp>

#include <chrono>
#include <iostream>
#include <source_location>
#include <string>

inline auto debug_total = std::chrono::steady_clock::duration{0};
template <typename F>
auto debug_txt_f(F callback, const std::source_location &loc = std::source_location::current()) -> void {
  using namespace std::chrono;
  const auto start = steady_clock::now();
  auto dump_str = std::string("[DEBUG] (") + loc.file_name() + ":L" + std::to_string(loc.line()) + ") ";
  dump_str += callback();
  std::cerr << dump_str << std::endl;
  const auto end = steady_clock::now();
  debug_total += end - start;
}
#  define debug_txt(txt) debug_txt_f([&]() { return to_pretty_str(txt); })

#else // DEBUG
#  define debug_txt(...)
#endif // DEBUG
