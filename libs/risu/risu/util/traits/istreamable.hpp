#pragma once

#include <risu/prelude.hpp>

#include <iostream>

template <typename T> concept istreamable_v = requires (T a) { std::cin >> a; };
