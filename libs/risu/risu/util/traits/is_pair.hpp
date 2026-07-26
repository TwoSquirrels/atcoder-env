#pragma once

#include <risu/prelude.hpp>

#include <utility>

template <typename T> inline constexpr bool is_pair_v = false;
template <typename T, typename U> inline constexpr bool is_pair_v<std::pair<T, U>> = true;

