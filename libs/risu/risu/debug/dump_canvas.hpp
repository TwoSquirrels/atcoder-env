#pragma once

#include <risu/prelude.hpp>

#ifdef DEBUG

#include <risu/debug/debug.hpp>

#include <cstddef>
#include <source_location>
#include <string>
#include <vector>

inline auto dump_canvas_f(std::string_view label, const std::vector<std::string> &canvas,
                          const std::source_location &loc = std::source_location::current()) -> void {
  debug_txt_f([&]() {
    auto h = int(canvas.size()), w = int(canvas[0].size());
    auto txt = std::string(label) + " (" + std::to_string(h) + " x " + std::to_string(w) + ")\n";
    auto ruler0 = std::string(w, ' '), ruler1 = std::string(w, '-');
    for (auto x = 0; x < w; x += 4) {
      auto x_s = std::to_string(x);
      for (auto i = std::size_t{0}; i < x_s.size(); ++i) ruler0[x - i] = x_s[x_s.size() - i - 1];
      ruler1[x] = '+';
    }
    txt += "   |" + ruler0 + "\n---+" + ruler1;
    for (auto y = 0; y < h; ++y) {
      auto ruler = std::string("   |");
      if (y % 4 == 0) {
        auto y_s = std::to_string(y);
        for (auto i = std::size_t{0}; i < y_s.size(); ++i) ruler[2 - i] = y_s[y_s.size() - i - 1];
        ruler[3] = '+';
      }
      txt += "\n" + ruler + canvas[y];
    }
    return txt;
  }, loc);
}
#  define dump_canvas(canvas) dump_canvas_f((#canvas), (canvas))

#else // DEBUG
#  define dump_canvas(...)
#endif // DEBUG
