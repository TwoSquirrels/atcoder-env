#pragma once

#include <risu/prelude.hpp>
#include <risu/io/output.hpp>

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#ifdef DEBUG
#  include <risu/io/input.hpp>
#  include <risu/debug/debug.hpp>
#  include <chrono>
#endif // DEBUG

inline auto cp_main() -> std::string;

#ifndef RISU_NO_MAIN

auto main() -> int {
  using namespace std;
#  ifdef DEBUG
  using namespace std::chrono;
  cerr << "[INFO] running in debug mode!" << endl;
  const auto start = steady_clock::now();
  try {
#  endif // DEBUG
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    cout << fixed << setprecision(12);
    // run!!!
    const auto result = cp_main();
    if (!result.empty()) write_stdout(result);
#  ifdef DEBUG
    write_stdout('\n', true);
  } catch (const exception &e) {
    write_stdout('\n', true);
    cerr << "[ERROR] " << e.what() << endl;
  }
  const auto end = steady_clock::now();
  const auto time_ms = duration_cast<milliseconds>(end - start - input_total - debug_total);
  cerr << "[INFO] finished in " << time_ms.count() << " ms!" << endl;
#  endif // DEBUG
  return 0;
}

#endif // RISU_NO_MAIN
