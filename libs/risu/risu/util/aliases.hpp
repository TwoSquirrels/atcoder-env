#pragma once

#include <functional>
#include <limits>
#include <queue>
#include <vector>

namespace risu {

using i8 = signed char;
using u8 = unsigned char;
static_assert(std::numeric_limits<u8>::digits == 8, "byte is not 8 bits");
using i32 = int;
using u32 = unsigned int;
static_assert(sizeof(i32) == 4 && sizeof(u32) == 4, "i32/u32 is not 32 bits");
using i64 = long long;
using u64 = unsigned long long;
static_assert(sizeof(i64) == 8 && sizeof(u64) == 8, "i64/u64 is not 64 bits");

using f32 = float;
static_assert(std::numeric_limits<f32>::is_iec559 && std::numeric_limits<f32>::digits == 24, "f32 is not IEEE binary32");
using f64 = double;
static_assert(std::numeric_limits<f64>::is_iec559 && std::numeric_limits<f64>::digits == 53, "f64 is not IEEE binary64");

template <typename T, typename Compare = std::less<T>> using heapq = std::priority_queue<T, std::vector<T>, Compare>;

}
