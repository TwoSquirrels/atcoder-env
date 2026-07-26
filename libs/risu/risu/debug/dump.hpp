#pragma once

#include <risu/prelude.hpp>

#ifdef DEBUG

#include <risu/debug/debug.hpp>
#include <risu/debug/pretty.hpp>

#include <cstddef>
#include <source_location>
#include <string>
#include <string_view>
#include <tuple>

// find the next comma that separates macro arguments,
// skipping over brackets and string/char literals
inline auto dump_label_end(std::string_view labels, std::size_t left) -> std::size_t {
  auto depth = 0;
  for (auto i = left; i < labels.size(); ++i) {
    const auto c = labels[i];
    if (c == '(' || c == '[' || c == '{') ++depth;
    else if (c == ')' || c == ']' || c == '}') --depth;
    else if (c == '"' || c == '\'') {
      while (++i < labels.size() && labels[i] != c) {
        if (labels[i] == '\\') ++i;
      }
    } else if (c == ',' && depth == 0) return i;
  }
  return labels.size();
}

inline auto dump_label_trim(std::string_view label) -> std::string_view {
  constexpr auto spaces = std::string_view(" \t\n\r");
  const auto first = label.find_first_not_of(spaces);
  if (first == std::string_view::npos) return {};
  return label.substr(first, label.find_last_not_of(spaces) - first + 1);
}

template <typename... Types>
auto dump_f(std::string_view labels, const std::tuple<Types...> &targets_tupl,
            const std::source_location &loc = std::source_location::current()) -> void {
  debug_txt_f([&]() {
    auto txt = std::string();
    std::apply([&labels, &txt](const auto &... targets) {
      auto i = 0;
      auto label_left = std::size_t{0};
      (([&](const auto &target) {
        const auto label_right = dump_label_end(labels, label_left);
        if (i >= 1) txt += ", ";
        txt += dump_label_trim(labels.substr(label_left, label_right - label_left));
        txt += ": " + to_pretty_str(target);
        label_left = label_right + 1;
        ++i;
      })(targets), ...);
      txt += ";";
    }, targets_tupl);
    return txt;
  }, loc);
}
#  define dump(...) dump_f((#__VA_ARGS__), std::forward_as_tuple(__VA_ARGS__))

#else // DEBUG
#  define dump(...)
#endif // DEBUG
