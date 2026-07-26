// 手元ビルドを速くするためのプリコンパイル対象。CMake が -include で読み込む。
// 提出コードには影響しないので、テンプレートと違って pragma は持たない。
#pragma once

#include <risu/std.hpp>

#if __has_include(<boost/multiprecision/cpp_int.hpp>)
#  include <boost/multiprecision/cpp_int.hpp>
#endif
#if __has_include(<boost/math/special_functions/prime.hpp>)
#  include <boost/math/special_functions/prime.hpp>
#endif
#if __has_include(<boost/multiprecision/miller_rabin.hpp>)
#  include <boost/multiprecision/miller_rabin.hpp>
#endif

#if __has_include(<atcoder/all>)
#  include <atcoder/all>
#endif

#include <risu/all.hpp>
