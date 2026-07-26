#pragma once

#include <risu/prelude.hpp>

#include <tuple>

template <typename T> inline constexpr bool is_tuple_v = false;
template <typename... Types> inline constexpr bool is_tuple_v<std::tuple<Types...>> = true;

