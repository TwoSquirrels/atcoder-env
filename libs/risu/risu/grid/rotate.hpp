#pragma once

#include <risu/prelude.hpp>

#include <vector>

template <typename T> auto rotate(std::vector<std::vector<T>> grid, int angle = 1) -> std::vector<std::vector<T>> {
  angle %= 4;
  if (angle == 0) return grid;
  auto h = int(grid.size()), w = int(grid[0].size());
  auto rotated = std::vector(w, std::vector<T>(h));
  for (auto y = 0; y < h; ++y) for (auto x = 0; x < w; ++x) rotated[w - 1 - x][y] = grid[y][x];
  return rotate(rotated, angle - 1);
}

