#pragma once

#include <risu/prelude.hpp>
#include <risu/util/iota.hpp>

#include <bit>
#include <iterator>
#include <ranges>
#include <string>

// functions
#define tostr to_string
#define p1 first
#define p2 second
#define eb emplace_back
#define ef emplace_front
#define pb pop_back
#define pf pop_front
// cast
#define sz(x) (i64(std::ssize(x)))
#define bit_width(x) (int(std::bit_width(x)))
// repeat
#define reps(i, l, r) for ([[maybe_unused]] const auto i : iota_range(l, r))
#define rep(i, n) reps(i, 0, n)
#define rreps(i, l, r) for ([[maybe_unused]] const auto i : iota_range(l, r) | std::views::reverse)
#define rrep(i, n) rreps(i, 0, n)
#define each(for_able) for (auto &&for_able##_i : (for_able))
// iterate
#define all(for_able) (std::begin(for_able)), (std::end(for_able))
