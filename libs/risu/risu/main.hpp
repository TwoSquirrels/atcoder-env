#pragma once

#include "io/output.hpp"
#ifdef DEBUG
#  include "io/input.hpp"
#  include "debug/debug.hpp"
#endif

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#ifdef DEBUG
#  include <chrono>
#endif

namespace risu {
auto main() -> void;
}

auto main() -> int {
  using namespace std;
  ios_base::sync_with_stdio(false); // so I must not use stdio.
  cout << fixed << setprecision(15);
#ifdef DEBUG
  using namespace std::chrono;
  using namespace risu::internal;
  cerr << "[INFO] running in debug mode!" << endl;
  // TODO: 直近のカバレッジを記録する何かを作る.
  auto exit_code = 0;
  const auto start = steady_clock::now();
  try {
#else
    cin.tie(nullptr); // so I need to flush manually before cin in interactive problems.
#endif
    risu::main();
#ifdef DEBUG
  } catch (const exception &e) {
    cout.flush();
    cerr << "\n[ERROR] " << e.what() << endl;
    exit_code = 1;
  } catch (...) {
    cout.flush();
    cerr << "\n[ERROR] unknown exception" << endl;
    exit_code = 1;
  }
  const auto end = steady_clock::now();
  const auto time_us = duration_cast<microseconds>(end - start - input_total - debug_total);
  cout.flush();
  cerr << "[INFO] finished in " << time_us.count() / 1000.0 << " ms!" << endl;
  return exit_code;
#else
  return 0;
#endif
}
