// SPDX-License-Identifier: MIT
// (c) 2023 TwoSquirrels
// my AtCoder environment: https://github.com/TwoSquirrels/atcoder-env

//#define DEBUG

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

/// abi
#if __has_include(<cxxabi.h>)
#  define INCLUDED_CXXABI
#  include <cxxabi.h>
#endif // <cxxabi.h>

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

// is_pair
template <typename T> inline constexpr bool is_pair_v = false;
template <typename T, typename U> inline constexpr bool is_pair_v<std::pair<T, U>> = true;

// is_tuple
template <typename T> inline constexpr bool is_tuple_v = false;
template <typename... Types> inline constexpr bool is_tuple_v<std::tuple<Types...>> = true;

// istreamable
#if __cplusplus >= 202002L
template <typename T> concept istreamable_v = requires (T a) { std::cin >> a; };
#else // C++20
template<typename T, typename = void> inline constexpr bool istreamable_v = false;
template<typename T>
inline constexpr bool istreamable_v<T, std::void_t<decltype(std::cin >> std::declval<T&>())>> = true;
#endif // C++20

// ostreamable
#if __cplusplus >= 202002L
template <typename T> concept ostreamable_v = requires (T a) { std::cout << a; };
#else // C++20
template<typename T, typename = void> inline constexpr bool ostreamable_v = false;
template<typename T>
inline constexpr bool ostreamable_v<T, std::void_t<decltype(std::cout << std::declval<T&>())>> = true;
#endif // C++20

// iterable
#if __cplusplus >= 202002L
#  if __has_include(<ranges>)
template <typename T> concept iterable_v = std::ranges::range<T>;
#  else // <ranges>
template <typename T> concept iterable_v = requires (T a) { std::begin(a); std::end(a); };
#  endif // C++20
#else // C++20
template <typename T> inline constexpr bool iterable_v =
  std::is_same_v<decltype(std::begin(std::declval<T>())), decltype(std::end(std::declval<T>()))>;
#endif // C++20

/// utils

template <typename T> inline std::string get_typename(std::size_t length_limit = std::string::npos) {
  std::string name;
#ifdef INCLUDED_CXXABI
  name = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, nullptr);
#else // INCLUDED_CXXABI
  name = typeid(T).name();
#endif // INCLUDED_CXXABI
#ifdef __ANDROID__
  name = std::regex_replace(name, std::regex("::__ndk1::"), "::");
#endif // __ANDROID__
  if (name.length() > length_limit) name = name.substr(0, length_limit - 3) + "...";
  return name;
}

template <typename T> constexpr T inf() {
  if constexpr (std::numeric_limits<T>::has_infinity) return std::numeric_limits<T>::infinity();
  return std::numeric_limits<T>::max() / 2.125L;
}

template <typename T> constexpr int sin90(T theta90) {
  if (theta90 % 2 == 0) return 0;
  return theta90 % 4 < 2 ? 1 : -1;
}
template <typename T> constexpr int cos90(T theta90) {
  return sin90(theta90 + 1);
}
template <typename T> constexpr int sin45(T theta45) {
  if (theta45 % 4 == 0) return 0;
  return theta45 % 8 < 4 ? 1 : -1;
}
template <typename T> constexpr int cos45(T theta45) {
  return sin45(theta45 + 2);
}

template <typename T, typename U> inline bool chmin(T &&a, const U b) {
  const bool compare = a > b;
  if (compare) a = b;
  return compare;
}
template <typename T, typename U> inline bool chmax(T &&a, const U b) {
  const bool compare = a < b;
  if (compare) a = b;
  return compare;
}

long long powi(int base, int exponent = 2) {
  long long ans = 1;
  for (int i = exponent; i != 0; i += (i >= 0 ? -1 : 1)) {
    if (i >= 0) ans *= base;
    else ans /= base;
    if (ans == 0) break;
  }
  return ans;
}

#ifdef INCLUDED_ACL
template <typename T = atcoder::modint>
T mint_inv(T x) {
  constexpr int memo_limit = 1 << 24;
  static std::array<T, memo_limit> memo;
  if (x.val() >= memo_limit) return x.inv();
  if (memo[x.val()] == 0) memo[x.val()] = x.inv();
  return memo[x.val()];
}
#endif // INCLUDED_ACL

template <typename T> bool is_prime(T n) {
  if (n <= 1) return false;
  if (n == 2 || n == 3) return true;
  if (n % 2 == 0 || n % 3 == 0) return false;
  // miller rabin
#ifdef INCLUDED_ACL
  if (n <= std::min(1LL << 32, static_cast<long long>(std::numeric_limits<int>::max()))) {
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
  for (int i = 2; i <= int(boost::math::max_prime); i++) {
    const auto prime_i = boost::math::prime(i);
    if (prime_i > n / prime_i) return true;
    if (n % prime_i == 0) return false;
    tried = prime_i;
  }
#endif // INCLUDED_BOOST_PRIME
  for (T i = (tried + 5) / 6 * 6; (i - 1) * (i - 1) <= n; i += 6) {
    if (n % (i - 1) == 0 || n % (i + 1) == 0) return false;
  }
  return true;
}

// thanks to https://perogram.hateblo.jp/entry/2019/03/29/193632
template <bool osa_k = false, typename T> std::vector<std::pair<T, int>> factors(T n) {
  constexpr int spf_limit = 1 << 24;
  std::vector<std::pair<T, int>> result;
  if (n < 0) {
    result.emplace_back(-1, 1);
    const auto natural = factors<osa_k>(-n);
    result.insert(result.end(), natural.begin(), natural.end());
  } else if (n == 0) result.emplace_back(0, 1);
  else if (n == 1);
  else if (osa_k && n < spf_limit) {
    static std::array<int, spf_limit> spf;
    if (spf[0] == 0) {
      for (int i = 0; i < spf_limit; i++) spf[i] = i % 2 == 0 ? 2 : i % 3 == 0 ? 3 : i;
      spf[0] = 1;
      for (int p1 = 6; (p1 - 1) * (p1 - 1) <= spf_limit; p1 += 6) {
        if (spf[p1 - 1] == p1 - 1) for (int i = p1 - 1; i < spf_limit; i += p1 - 1) if (spf[i] == i) spf[i] = p1 - 1;
        if (spf[p1 + 1] == p1 + 1) for (int i = p1 + 1; i < spf_limit; i += p1 + 1) if (spf[i] == i) spf[i] = p1 + 1;
      }
    }
    while (n != 1) {
      if (!result.empty() && result.back().first == spf[n]) result.back().second++;
      else result.emplace_back(spf[n], 1);
      n /= spf[n];
    }
  } else {
    int expo = 0;
    while (n % 2 == 0) { n /= 2; expo++; }
    if (expo != 0) result.emplace_back(2, expo);
    expo = 0;
    while (n % 3 == 0) { n /= 3; expo++; }
    if (expo != 0) result.emplace_back(3, expo);
    T tried = 3;
#ifdef INCLUDED_BOOST_PRIME
    for (int i = 2; i <= int(boost::math::max_prime); i++) {
      const auto prime_i = boost::math::prime(i);
      if (prime_i > n / prime_i) break;
      expo = 0;
      while (n % prime_i == 0) { n /= prime_i; expo++; }
      result.emplace_back(prime_i, expo);
      tried = prime_i;
    }
    if (osa_k && n < spf_limit) {
      const auto rest = factors<osa_k>(n);
      result.insert(result.end(), rest.begin(), rest.end());
      n = 1;
    }
#endif // INCLUDED_BOOST_PRIME
    if (is_prime(n)) { result.emplace_back(n, 1); n = 1; }
    for (T i = (tried + 5) / 6 * 6; (i - 1) * (i - 1) <= n; i += 6) {
      expo = 0;
      while (n % (i - 1) == 0) { n /= (i - 1); expo++; }
      if (expo != 0) result.emplace_back(i - 1, expo);
      expo = 0;
      while (n % (i + 1) == 0) { n /= (i + 1); expo++; }
      if (expo != 0) result.emplace_back(i + 1, expo);
      if (osa_k && n < spf_limit) {
        const auto rest = factors<osa_k>(n);
        result.insert(result.end(), rest.begin(), rest.end());
        n = 1;
      }
    }
    if (n != 1) result.emplace_back(n, 1);
  }
  return result;
}

template <bool osa_k = false, typename T> std::vector<T> divisors(T n) {
  std::vector<T> result(1, 1);
  for (auto [prime, expo] : factors<osa_k>(n)) {
    const int result_size = result.size();
    T pow_i = 1;
    for (int i = 1; i <= expo; i++) {
      pow_i *= prime;
      for (int k = 0; k < result_size; k++) result.emplace_back(result[k] * pow_i);
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

template <typename T> T fact(int n, bool inv = false) {
  assert(n >= 0);
  static std::vector<std::pair<T, T>> factorials = { { 1, 1 } };
  for (int i = factorials.size(); i <= n; i++) factorials.emplace_back(i * factorials[i - 1].first, 0);
  if (inv && factorials[n].second == 0) factorials[n].second = mint_inv(factorials[n].first);
  return inv ? factorials[n].second : factorials[n].first;
}
template <typename T> T perm(int n, int k) { return fact<T>(n) * fact<T>(n - k, true); }
template <typename T> T comb(int n, int k) { return perm<T>(n, k) * fact<T>(k, true); }

template <typename T> std::vector<std::vector<T>> rotate(std::vector<std::vector<T>> grid, int angle = 1) {
  angle %= 4;
  if (angle == 0) return grid;
  auto h = grid.size(), w = grid[0].size();
  auto rotated = std::vector(w, std::vector<T>(h));
  for (int y = 0; y < h; y++) for (int x = 0; x < w; x++) rotated[w - 1 - x][y] = grid[y][x];
  return rotate(rotated, angle - 1);
}

std::string int128_to_str(__int128_t target) {
  std::string target_str;
  __uint128_t target_tmp = target < 0 ? -target : target;
  do {
    target_str += '0' + target_tmp % 10;
    target_tmp /= 10;
  } while (target_tmp != 0);
  if (target < 0) target_str += '-';
  std::reverse(target_str.begin(), target_str.end());
  return target_str;
}

template <typename T> std::string to_pretty_str(T target) {
  using namespace std;
  string str;
  if constexpr (is_void_v<T>) str += "void"s;
  else if constexpr (is_null_pointer_v<T>) str += "null"s;
  else if constexpr (is_same_v<T, bool>) str += target ? "true"s : "false"s;
  else if constexpr (is_same_v<T, char> || is_same_v<T, char16_t> || is_same_v<T, char32_t> || is_same_v<T, wchar_t>) {
    str += "'"s + target + "'"s;
#ifdef INCLUDED_ACL
    } else if constexpr (atcoder::internal::is_modint<T>::value) {
    str += to_string(target.val()) + "(mod)"s;
#endif // INCLUDED_ACL
  } else if constexpr (is_arithmetic_v<T>) {
    if constexpr (is_same_v<T, __int128_t>) str += int128_to_str(target);
    else str += to_string(target);
    if constexpr (is_unsigned_v<T>) str += "u"s;
    if constexpr (is_same_v<remove_cv_t<T>, long>) str += "L"s;
    else if constexpr (is_same_v<remove_cv_t<T>, long long>) str += "LL"s;
    else if constexpr (is_same_v<T, __int128_t>) str += "LLL"s;
  } else if constexpr (is_pair_v<T>) {
    str += "("s + to_pretty_str(target.first);
    str += ", "s + to_pretty_str(target.second) + ")"s;
  } else if constexpr (is_convertible_v<T, string>) {
    str += "\""s + target + "\""s;
  } else if constexpr (is_array_v<T>) {
    str += "["s;
    bool separate = false;
    for (const auto &target_i : target) {
      if (separate) str += ", "s;
      str += to_pretty_str(target_i);
      separate = true;
    }
    str += "]"s;
  } else if constexpr (iterable_v<T>) {
    str += get_typename<T>(20) + "{"s;
    bool separate = false;
    for (const auto &target_i : target) {
      if (separate) str += ","s;
      str += " "s + to_pretty_str(target_i);
      separate = true;
    }
    if (separate) str += " "s;
    str += "}"s;
  } else {
    str += "<"s + get_typename<T>(20);
    str += " ("s + to_string(sizeof(target)) + " byte)>"s;
  }
  return str;
}

// input
#ifdef DEBUG
std::chrono::milliseconds input_total_ms{0};
#endif // DEBUG
template <typename T> inline void read_stdin(T &&target) {
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
template <typename T> inline T input(T &&target) {
  using T_V = std::remove_reference_t<T>;
  if constexpr (istreamable_v<T_V>) read_stdin(target);
  else if constexpr (iterable_v<T_V>) for (auto &&target_i : target) input(target_i);
  else if constexpr (is_pair_v<T_V>) {
    input(target.first);
    input(target.second);
  } else if constexpr (std::is_convertible_v<long long, T_V>) {
    long long n;
    target = input(n);
  } else {
#ifdef DEBUG
    std::cerr << "[WARN] Type " << get_typename<T_V>()
              << " is invalid, so input is skipped;" << std::endl;
#endif //DEBUG
  }
  return target;
}
// input and initialize
struct Scanner { template <typename T> inline operator T() const { T target; return input(target); } } scan;

// output
template <typename T> inline void write_stdout(T target, bool flush = false) {
  std::cout << target;
  if (flush) std::cout << std::flush;
}
template <typename T, typename Sep = char> inline void output(T target, Sep separator = ' ', bool flush = false) {
  if constexpr (ostreamable_v<T>) {
    write_stdout(target, flush);
  } else if constexpr (std::is_convertible<T, __int128_t>::value) {
    write_stdout(int128_to_str(target), flush);
#  ifdef INCLUDED_ACL
    } else if constexpr (atcoder::internal::is_modint<T>::value) {
    output(target.val(), separator, flush);
#  endif // INCLUDED_ACL
  } else if constexpr (iterable_v<T>) {
    bool separate = false;
    for (const auto target_i : target) {
      if (separate) write_stdout(separator);
      output(target_i, separator);
      separate = true;
    }
    if (flush) write_stdout("", flush);
  } else if constexpr (is_pair_v<T>) {
    output(target.first, separator);
    write_stdout(separator);
    output(target.second, separator, flush);
  } else {
    write_stdout("<unknown>", flush);
  }
}
template <typename T, typename Sep = char> inline void outputln(T target, Sep separator = ' ', bool flush = false) {
  output(target, separator);
  write_stdout('\n', flush);
}

// dump
#ifdef DEBUG
std::chrono::milliseconds dump_total_ms{0};
template <typename... Types>
void dump_stderr(std::string labels, std::tuple<Types...> targets_tupl, int line = -1, std::string file = (__FILE__)) {
  using namespace std::chrono;
  const auto start = system_clock::now();
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
      std::cerr << to_pretty_str(target);
      label_left += label_len + 1;
      i++;
    })(targets), ...);
    std::cerr << ";" << std::endl;
  }, targets_tupl);
  const auto end = system_clock::now();
  dump_total_ms += duration_cast<milliseconds>(end - start);
}
#  define dump(...) (dump_stderr((#__VA_ARGS__), std::make_tuple(__VA_ARGS__), (__LINE__), (__FILE__)), true)
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
  } catch (exception &e) {
    write_stdout('\n', true);
    cerr << "[ERROR] " << e.what() << endl;
  }
  const auto end = system_clock::now();
  const auto time_ms = duration_cast<milliseconds>(end - start) - input_total_ms - dump_total_ms;
  cerr << "[INFO] finished in " << time_ms.count() << " ms!" << endl;
#  endif // DEBUG
  return 0;
}

#endif // __INCLUDE_LEVEL__

/// aliases

using i32 = int; using u32 = unsigned int;
using i64 = long long; using u64 = unsigned long long;
using i128 = __int128; using u128 = unsigned __int128;
using f32 = float; using f64 = double; using f128 = long double;
using str = std::string;
template <typename T> using vec = std::vector<T>;
template <typename T> using deq = std::deque<T>;
template <typename T> using list = std::list<T>;
template <typename T, typename Compare = std::less<T>> using p_que = std::priority_queue<T, std::vector<T>, Compare>;
template <typename Key, typename Compare = std::less<Key>> using mset = std::multiset<Key, Compare>;
template <typename Key, typename T, typename Compare = std::less<Key>> using mmap = std::multimap<Key, T, Compare>;
template <typename Key> using u_set = std::unordered_set<Key>;
template <typename Key> using u_mset = std::unordered_multiset<Key>;
template <typename Key, typename T> using u_map = std::unordered_map<Key, T>;
template <typename Key, typename T> using u_mmap = std::unordered_multimap<Key, T>;
template <size_t N> using bset = std::bitset<N>;
#ifdef INCLUDED_BOOST_CPP_INT
using bigint = boost::multiprecision::cpp_int;
#endif // INCLUDED_BOOST_CPP_INT

using namespace std;
#if __cplusplus >= 202002L && __has_include(<ranges>)
namespace rng = std::ranges;
namespace viw = std::ranges::views;
#endif // C++20 && <ranges>
#ifdef INCLUDED_ACL
using namespace atcoder;
#endif // INCLUDED_ACL

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
constexpr auto infi = inf<int>();
constexpr auto infl = inf<long long>();
constexpr auto infd = inf<double>();
constexpr auto infld = inf<long double>();
const std::array YN = { "No", "Yes" };
const std::array AB = { "Bob", "Alice" };
const std::array FS = { "Second", "First" };
const std::array TA = { "Aoki", "Takahashi" };
#pragma GCC diagnostic pop

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
#define reps(i, l, r) for (std::decay_t<decltype(r)> i##_right = (r), i = (l); i < i##_right; i++)
#define rep(i, n) reps(i, 0, n)
#define rreps(i, l, r) for (std::decay_t<decltype(l)> i##_left = (l), i = (r) - 1; i >= i##_left; i--)
#define rrep(i, n) rreps(i, 0, n)
#define each(for_able) for (auto &&for_able##_i : (for_able))
// iterate
#define all(for_able) (std::begin(for_able)), (std::end(for_able))
#define rev(for_able) (std::rbegin(for_able)), (std::rend(for_able))
// lambda
#define pred(x, expr) ([&](auto &&x) -> bool { return (expr); })
#define comp(x, y, expr) ([&](auto &&x, auto &&y) -> bool { return (expr); })

/// answer

//using mint=modint998244353;

inline str cp_main() {
  
  return "";
}
