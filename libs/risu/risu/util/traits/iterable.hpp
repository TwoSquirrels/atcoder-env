#pragma once

#include <risu/prelude.hpp>

#include <ranges>

template <typename T> concept iterable_v = std::ranges::range<T>;
