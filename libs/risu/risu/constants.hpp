#pragma once

#include <risu/prelude.hpp>
#include <risu/util/inf.hpp>

#include <array>
#include <string>

#pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-variable"
constexpr auto infi = inf<int>();
constexpr auto infl = inf<long long>();
constexpr auto infd = inf<double>();
constexpr auto infld = inf<long double>();
const auto YN = std::array<std::string, 2>{ "No", "Yes" };
const auto AB = std::array<std::string, 2>{ "Bob", "Alice" };
const auto FS = std::array<std::string, 2>{ "Second", "First" };
const auto TA = std::array<std::string, 2>{ "Aoki", "Takahashi" };
#pragma GCC diagnostic pop
