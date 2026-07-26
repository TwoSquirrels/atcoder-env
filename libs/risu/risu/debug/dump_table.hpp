#pragma once

#include <risu/prelude.hpp>

#ifdef DEBUG

#include <risu/debug/debug.hpp>
#include <risu/debug/pretty.hpp>
#include <risu/util/chminmax.hpp>

#include <iomanip>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>

// TODO: fix
template <typename T>
auto dump_table_f(std::string label, const std::vector<std::vector<T>> &table,
                  const std::source_location &loc = std::source_location::current()) -> void {
  debug_txt_f([&]() {
    if (table.empty()) return label + ": empty table";
    auto h = int(table.size()), w = int(table[0].empty() ? 0 : table[0].size());
    auto s_table = std::vector(h, std::vector<std::string>(w));
    auto col_width = std::vector<int>(w);
    for (auto x = 0; x < w; ++x) col_width[x] = std::to_string(x).length();
    for (auto y = 0; y < h; ++y) {
      for (auto x = 0; x < w; ++x) {
        s_table[y][x] = to_pretty_str(table[y][x]);
        chmax(col_width[x], (int)s_table[y][x].length());
      }
    }
    auto row_idx_width = int(std::to_string(h - 1).length());
    auto ss = std::stringstream();
    ss << label << " (" << h << " x " << w << ")\n";
    ss << std::string(row_idx_width + 2, ' ');
    for (auto x = 0; x < w; ++x) ss << " " << std::setw(col_width[x]) << x;
    ss << "\n" << std::string(row_idx_width + 2, ' ') << "+";
    for (auto x = 0; x < w; ++x) ss << std::string(col_width[x] + 1, '-') << ((x == w - 1) ? "+" : "");
    ss << "\n";
    for (auto y = 0; y < h; ++y) {
      ss << " " << std::setw(row_idx_width) << y << "|";
      for (auto x = 0; x < w; ++x) ss << " " << std::setw(col_width[x]) << s_table[y][x];
      ss << "\n";
    }
    return ss.str();
  }, loc);
}
#  define dump_table(table) dump_table_f((#table), (table))

#else // DEBUG
#  define dump_table(...)
#endif // DEBUG
