// SPDX-License-Identifier: MIT
// (c) 2023 TwoSquirrels
// my AtCoder environment: https://github.com/TwoSquirrels/atcoder-env

//#define DEBUG

#ifndef INCLUDED_MAIN
#  define INCLUDED_MAIN
#  include __FILE__

//using mint=modint998244353;

inline std::string cp_main() {
  
  return"";
}

#else // INCLUDED_MAIN

// enable debug mode when compiled with -DDEBUG
#ifdef DEBUG
#  define IS_DEBUG (1)
#else
// QCFium method
#  pragma GCC target("avx2")
#  pragma GCC optimize("O3")
#  pragma GCC optimize("unroll-loops")
#  pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#  define IS_DEBUG (0)
#endif

/// includes

// standard libraries
#if __has_include(<bits/stdc++.h>)
#  include <bits/stdc++.h>
#else // <bits/stdc++.h>
// C++17 (clang)
#  include <cassert>
#  include <cfenv>
#  include <cfloat>
#  include <ciso646>
#  include <clocale>
#  include <csetjmp>
#  include <csignal>
#  include <cstdbool>
#  include <cinttypes>
#  include <charconv>
#  include <typeindex>
#  include <any>
#  include <scoped_allocator>
#  include <forward_list>
#  include <list>
#  include <map>
#  include <set>
#  include <valarray>
#  include <variant>
#  include <unordered_map>
#  include <unordered_set>
#  include <queue>
#  include <condition_variable>
#  include <shared_mutex>
#  include <codecvt>
#  include <future>
#  include <regex>
#  include <iostream>
#  include <random>
#  include <ctgmath>
#  include <fstream>
#endif // <bits/stdc++.h>

// boost libraries
#if __has_include(<boost/multiprecision/cpp_int.hpp>)
#  define INCLUDED_BOOST_CPP_INT
#  include <boost/multiprecision/cpp_int.hpp>
#endif
#if __has_include(<boost/math/special_functions/prime.hpp>)
#  define INCLUDED_BOOST_PRIME
#  include <boost/math/special_functions/prime.hpp>
#endif
#if __has_include(<boost/multiprecision/miller_rabin.hpp>)
#  define INCLUDED_BOOST_MILLER_RABIN
#  include <boost/multiprecision/miller_rabin.hpp>
#endif

// atcoder libraries
#if __has_include(<atcoder/all>)
#  define INCLUDED_ACL
#  include <atcoder/all>
#endif

/// sfinae

#if __cplusplus >= 202002L
template <typename T>
concept IsInputable = requires (T a) {
  std::cin >> a;
};
template <typename T>
concept IsOutputable = requires (T a) {
  std::cout << a;
};
template <typename T>
concept IsMutableNumber = requires (T a) {
  a = 0LL;
};
template <typename T>
concept IsPair = requires (T a) {
  a.first;
  a.second;
};
template <typename T>
concept IsEachable = requires (T a) {
  std::begin(a);
  std::end(a);
};
#endif // C++20

/// utils

template <typename T>
constexpr T inf() {
  if (std::numeric_limits<T>::has_infinity) {
    return std::numeric_limits<T>::infinity();
  }
  return std::numeric_limits<T>::max() / 2;
}

template <typename T>
constexpr int sin90(T theta90) {
  if (theta90 % 2 == 0) return 0;
  return theta90 % 4 < 2 ? 1 : -1;
}
template <typename T>
constexpr int cos90(T theta90) {
  return sin90(theta90 + 1);
}
template <typename T>
constexpr int sin45(T theta45) {
  if (theta45 % 4 == 0) return 0;
  return theta45 % 8 < 4 ? 1 : -1;
}
template <typename T>
constexpr int cos45(T theta45) {
  return sin90(theta45 + 2);
}

template <typename T, typename U>
inline bool chmin(T &&a, const U b) {
  const bool compare = a > b;
  if (compare) a = b;
  return compare;
}
template <typename T, typename U>
inline bool chmax(T &&a, const U b) {
  const bool compare = a < b;
  if (compare) a = b;
  return compare;
}

long long int powi(int base, int exponent = 2) {
  long long int ans = 1;
  for (int i = exponent; i != 0; i += (i >= 0 ? -1 : 1)) {
    if (i >= 0) ans *= base;
    else ans /= base;
    if (ans == 0) break;
  }
  return ans;
}

template <typename T>
bool is_prime(T n) {
  if (n <= 1) return false;
  if (n == 2 || n == 3) return true;
  if (n % 2 == 0 || n % 3 == 0) return false;
  // miller rabin
#ifdef INCLUDED_ACL
  if (n <= 1LL << 32 && n <= std::numeric_limits<int>::max()) {
    return atcoder::internal::is_prime_constexpr(n);
  }
#endif // INCLUDED_ACL
#ifdef INCLUDED_BOOST_MILLER_RABIN
  return boost::multiprecision::miller_rabin_test(n, 25);
#endif // INCLUDED_BOOST_MILLER_RABIN
  // binary search
#ifdef INCLUDED_BOOST_PRIME
  if (n <= boost::math::prime(boost::math::max_prime)) {
    int left = 0, right = boost::math::max_prime + 1;
    while (right - left > 1) {
      auto mid = (left + right) / 2;
      (boost::math::prime(mid) <= n ? left : right) = mid;
    }
    return n == boost::math::prime(left);
  }
#endif // INCLUDED_BOOST_PRIME
  // trial
  T tried = 3;
#ifdef INCLUDED_BOOST_PRIME
  for (int i = 2; i <= boost::math::max_prime; i++) {
    const auto prime_i = boost::math::prime(i);
    if (prime_i > n / prime_i) return true;
    if (n % prime_i == 0) {
      return false;
    }
    tried = prime_i;
  }
#endif // INCLUDED_BOOST_PRIM
  for (T i = (tried + 5) / 6 * 6; i * i <= n; i += 6) {
    if (n % (i - 1) == 0 || n % (i + 1) == 0) {
      return false;
    }
  }
  return true;
}

// TODO: divisor enumeration, prime factorization

template <typename T>
T fact(int n) {
  assert(n >= 0);
  static std::vector<T> factorials = { 1 };
  for (int i = factorials.size(); i <= n; i++) {
    factorials.emplace_back(i * factorials[i - 1]);
  }
  return factorials[n];
}
template <typename T>
T perm(int n, int k) {
  return fact<T>(n) / fact<T>(n - k);
}
template <typename T>
T comb(int n, int k) {
  return perm<T>(n, k) / fact<T>(k);
}

inline std::string yn(bool yes) {
  return yes ? "Yes" : "No";
}

template <typename T>
std::string to_pretty_str(T target) {
  return "<TODO>";
}

// input
#ifdef DEBUG
std::chrono::milliseconds input_total_ms{0};
#endif // DEBUG
template <typename T>
inline void read_stdin(T &&target) {
#ifdef DEBUG
  using namespace std::chrono;
  const auto start = system_clock::now();
  std::cin >> target;
  const auto end = system_clock::now();
  input_total_ms += duration_cast<milliseconds>(end - start);
#else // DEBUG
  std::cin >> target;
#endif // DEBUG
}
template <typename T>
inline T input(T &&target) {
#if __cplusplus >= 202002L
  if constexpr (IsInputable<T>) {
    read_stdin(target);
  } else if constexpr (IsEachable<T>) {
    for (auto &&target_i : target) input(target_i);
  } else if constexpr (IsPair<T>) {
    input(target.first);
    input(target.second);
  } else if constexpr (IsMutableNumber<T>) {
    long long int n;
    target = input(n);
  } else {
    // skip
  }
#else // C++20
  read_stdin(target);
#endif // C++20
  return target;
}
// input and initialize
struct Scanner {
  template <typename T>
  inline operator T() const {
    T target;
    return input(target);
  }
} scan;

// output
template <typename T>
inline void write_stdout(T target, bool flush = false) {
  std::cout << target;
  if (flush) std::cout << std::flush;
}
template <typename T, typename Sep = char>
inline void output(T target, Sep separator = ' ', bool flush = false) {
#if __cplusplus >= 202002L
  if constexpr (IsOutputable<T>) {
    write_stdout(target, flush);
  } else if constexpr (std::is_convertible<T, __int128_t>::value) {
    std::string target_str = "";
    __uint128_t target_tmp = target < 0 ? -target : target;
    do {
      target_str += '0' + target_tmp % 10;
      target_tmp /= 10;
    } while (target_tmp != 0);
    if (target < 0) target_str += '-';
    std::reverse(target_str.begin(), target_str.end());
    write_stdout(target_str, flush);
#ifdef INCLUDED_ACL
  } else if constexpr (atcoder::internal::is_modint<T>::value) {
    output(target.val(), separator, flush);
#endif // INCLUDED_ACL
  } else if constexpr (IsEachable<T>) {
    bool separate = false;
    for (const auto target_i : target) {
      if (separate) write_stdout(separator);
      output(target_i, separator);
      separate = true;
    }
    if (flush) write_stdout("", flush);
  } else if constexpr (IsPair<T>) {
    output(target.first, separator);
    write_stdout(separator);
    output(target.second, separator, flush);
  } else {
    write_stdout("<unknown>", flush);
  }
#else // C++20
  write_stdout(target, flush);
#endif // C++20
}
template <typename T, typename Sep = char>
inline void outputln(T target, Sep separator = ' ', bool flush = false) {
  output(target, separator);
  write_stdout('\n', flush);
}

// dump
template <typename... Types>
void dump_stderr(std::string labels,
                 std::tuple<Types...> targets_tupl,
                 int line = -1, std::string file = (__FILE__)) {
  std::cerr << "[DEBUG] ";
  if (line >= 0) {
    std::cerr << "(";
    if (file != (__FILE__)) std::cerr << file << ":";
    std::cerr << "L" << line << ") ";
  }
  std::apply([labels](auto... targets) {
    int i = 0, label_left = 0;
    (([&](auto target) {
      const auto label_len = labels.find(',', label_left) - label_left;
      const auto label =
        std::regex_replace(labels.substr(label_left, label_len),
                           std::regex("^\\s+|\\s+$"), "");
      if (i >= 1) std::cerr << ", ";
      std::cerr << label << ": ";
      std::cerr << /*to_pretty_str*/(target);
      label_left += label_len + 1;
      i++;
    })(targets), ...);
    std::cerr << ";" << std::endl;
  }, targets_tupl);
}
#ifdef DEBUG
#  define dump(...) (dump_stderr((#__VA_ARGS__),                \
                                 std::make_tuple(__VA_ARGS__),  \
                                 (__LINE__), (__FILE__)), true)
#else // DEBUG
#  define dump(...) (false)
#endif // DEBUG

/// main

#if !defined(__INCLUDE_LEVEL__) || __INCLUDE_LEVEL__ <= 1

inline std::string cp_main();

int main() {
  using namespace std;
  using namespace std::chrono;
#  ifdef DEBUG
  cerr << "[INFO] running in debug mode!" << endl;
  const auto start = system_clock::now();
  try {
#  endif // DEBUG
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    cout << fixed << setprecision(8);
    // run!!!
    const auto result = cp_main();
    if (!result.empty()) write_stdout(result);
#  ifdef DEBUG
    write_stdout('\n', true);
  } catch (exception e) {
    write_stdout('\n', true);
    cerr << "[ERROR] " << e.what() << endl;
  }
  const auto end = system_clock::now();
  const auto time_ms =
    duration_cast<milliseconds>(end - start) - input_total_ms;
  cerr << "[INFO] finished in " << time_ms.count() << " ms!" << endl;
#  endif // DEBUG
  return 0;
}

#endif // __INCLUDE_LEVEL__

/// aliases

using i32 = int;
using i64 = long long int;
using i128 = __int128;
using f32 = float;
using f64 = double;
using f128 = long double;
using str = std::string;
template <typename T> using vec = std::vector<T>;
template <typename T> using deq = std::deque<T>;
template <typename T> using list = std::list<T>;
template <typename T, typename Compare = std::less<T>> using p_que =
  std::priority_queue<T, std::vector<T>, Compare>;
template <typename Key, typename Compare = std::less<Key>> using mset =
  std::multiset<Key, Compare>;
template <typename Key, typename T, typename Compare = std::less<Key>> using mmap =
  std::multimap<Key, T, Compare>;
template <typename Key> using u_set = std::unordered_set<Key>;
template <typename Key> using u_mset = std::unordered_multiset<Key>;
template <typename Key, typename T> using u_map = std::unordered_map<Key, T>;
template <typename Key, typename T> using u_mmap = std::unordered_multimap<Key, T>;
template <size_t N> using bset = std::bitset<N>;
#ifdef INCLUDED_CPP_INT
using bigint = boost::multiprecision::cpp_int;
#endif // INCLUDED_CPP_INT

using namespace std;
#if __cplusplus >= 202002L && __has_include(<ranges>)
namespace rng = std::ranges;
namespace viw = std::ranges::views;
#endif // C++20 && <ranges>
#ifdef INCLUDED_ACL
using namespace atcoder;
#endif // INCLUDED_ACL

constexpr auto infi = inf<int>();
constexpr auto infl = inf<long long int>();
constexpr auto infd = inf<double>();
constexpr auto infld = inf<long double>();

// functions
#define tostr to_string
#define dist distance
#define min_e min_element
#define max_e max_element
#define iota_v iota_view
#define p1 first
#define p2 second
#define has contains
#define bgn begin
#define rbgn rbegin
#define eb emplace_back
#define ef emplace_front
#define pb pop_back
#define pf pop_front
// repeat
#define rep(i, n) for (decltype(n) i##_len = (n), i = 0; i < i##_len; i++)
#define reps(i, l, r) for (decltype(r) i##_right = (r), i = (l); i < i##_right; i++)
#define rrep(i, n) for (decltype(n) i##_len = (n), i = i##_len - 1; i >= 0; i--)
#define rreps(i, l, r) for (decltype(l) i##_left = (l), i = (r) - 1; i >= i##_left; i--)
#define each(for_able) for (auto &&for_able##_i : (for_able))
// iterate
#define all(for_able) (std::begin(for_able)), (std::end(for_able))
#define rev(for_able) (std::rbegin(for_able)), (std::rend(for_able))
// lambda
#define pred(x, expression) ([&](auto &&x) -> bool { return (expression); })
#define comp(x, y, expression) ([&](auto &&x, auto &&y) -> bool { return (expression); })

#endif // INCLUDED_MAIN
