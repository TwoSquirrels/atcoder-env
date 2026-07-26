#pragma once

#include <risu/prelude.hpp>

#include <functional>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#if __has_include(<boost/multiprecision/cpp_int.hpp>)
#  include <boost/multiprecision/cpp_int.hpp>
#endif

using i32 = int; using u32 = unsigned int;
using i64 = long long; using u64 = unsigned long long;
using i128 = __int128; using u128 = unsigned __int128;
using f32 = float; using f64 = double; using f80 = long double;
using str = std::string;
template <typename T> using vec = std::vector<T>;
template <typename T, typename Compare = std::less<T>> using p_que = std::priority_queue<T, std::vector<T>, Compare>;
template <typename Key, typename Compare = std::less<Key>> using mset = std::multiset<Key, Compare>;
template <typename Key, typename T, typename Compare = std::less<Key>> using mmap = std::multimap<Key, T, Compare>;
template <typename Key> using u_set = std::unordered_set<Key>;
template <typename Key> using u_mset = std::unordered_multiset<Key>;
template <typename Key, typename T> using u_map = std::unordered_map<Key, T>;
template <typename Key, typename T> using u_mmap = std::unordered_multimap<Key, T>;
#if __has_include(<boost/multiprecision/cpp_int.hpp>)
using bigint = boost::multiprecision::cpp_int;
#endif // <boost/multiprecision/cpp_int.hpp>
