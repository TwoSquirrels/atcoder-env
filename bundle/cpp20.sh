#!/bin/bash

set -eo pipefail

cd "$(dirname "$0")/../"

risundle main.cc > bundled.cc
if ! g++ -std=gnu++20 -fsyntax-only bundled.cc; then
    echo '[bundle/cpp20.sh] falling back to --no-tree-shaking' >&2
    risundle --no-tree-shaking main.cc > bundled.cc
    g++ -std=gnu++20 -fsyntax-only bundled.cc
fi

echo '// SPDX-FileCopyrightText: 2026 TwoSquirrels'
echo '// SPDX-License-Identifier: MIT'
echo '// My competitive programming environment: https://github.com/TwoSquirrels/atcoder-env'
echo ''
echo '// GCC 13.x-14.1: target pragma before libstdc++ headers triggers an always_inline/target mismatch CE in std::allocator.'
echo '#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ == 13 || (__GNUC__ == 14 && __GNUC_MINOR__ < 2))'
echo '#  include <bits/stdc++.h>'
echo '#endif'
echo ''
cat bundled.cc
