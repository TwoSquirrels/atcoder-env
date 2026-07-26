#pragma once

#include <risu/prelude.hpp>

#include <iostream>

template <typename T> concept ostreamable_v = requires (T a) { std::cout << a; };
