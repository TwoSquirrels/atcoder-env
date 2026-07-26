#pragma once

#include <risu/prelude.hpp>

template <typename Grid> auto rotated(const Grid &grid, int angle = 1) -> Grid {
  angle = (angle % 4 + 4) % 4;
  if (angle == 0 || grid.empty()) return grid;
  using Row = typename Grid::value_type;
  using Cell = typename Row::value_type;
  const auto h = int(grid.size()), w = int(grid[0].size());
  const auto turned = angle % 2 == 1;
  auto result = Grid(turned ? w : h, Row(turned ? h : w, Cell{}));
  for (auto y = 0; y < h; ++y) {
    for (auto x = 0; x < w; ++x) {
      if (angle == 1) result[w - 1 - x][y] = grid[y][x];
      else if (angle == 2) result[h - 1 - y][w - 1 - x] = grid[y][x];
      else result[x][h - 1 - y] = grid[y][x];
    }
  }
  return result;
}
