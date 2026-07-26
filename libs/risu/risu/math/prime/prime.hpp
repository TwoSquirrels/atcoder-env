#pragma once

#include <risu/prelude.hpp>

#include <algorithm>
#include <limits>
#if __has_include(<atcoder/math>)
#  include <atcoder/math>
#endif
#if __has_include(<boost/math/special_functions/prime.hpp>)
#  include <boost/math/special_functions/prime.hpp>
#endif
#if __has_include(<boost/multiprecision/miller_rabin.hpp>)
#  include <boost/multiprecision/miller_rabin.hpp>
#endif

template <typename T> auto is_prime(T n) -> bool {
  if (n <= 1) return false;
  if (n == 2 || n == 3) return true;
  if (n % 2 == 0 || n % 3 == 0) return false;
  // miller rabin
#if __has_include(<atcoder/math>)
  if (n <= std::min(1LL << 32, static_cast<long long>(std::numeric_limits<int>::max()))) {
    return atcoder::internal::is_prime_constexpr(n);
  }
#endif // <atcoder/math>
#if __has_include(<boost/multiprecision/miller_rabin.hpp>)
  return boost::multiprecision::miller_rabin_test(n, 25);
#endif // <boost/multiprecision/miller_rabin.hpp>
  // binary search
#if __has_include(<boost/math/special_functions/prime.hpp>)
  if (n <= boost::math::prime(boost::math::max_prime)) {
    auto left = 0, right = int(boost::math::max_prime + 1);
    while (right - left > 1) {
      auto mid = (left + right) / 2;
      (boost::math::prime(mid) <= n ? left : right) = mid;
    }
    return n == boost::math::prime(left);
  }
#endif // <boost/math/special_functions/prime.hpp>
  // trial
  auto tried = T(3);
#if __has_include(<boost/math/special_functions/prime.hpp>)
  for (auto i = 2; i <= int(boost::math::max_prime); ++i) {
    const auto prime_i = boost::math::prime(i);
    if (prime_i > n / prime_i) return true;
    if (n % prime_i == 0) return false;
    tried = prime_i;
  }
#endif // <boost/math/special_functions/prime.hpp>
  for (auto i = (tried + 5) / 6 * 6; (i - 1) * (i - 1) <= n; i += 6) {
    if (n % (i - 1) == 0 || n % (i + 1) == 0) return false;
  }
  return true;
}

