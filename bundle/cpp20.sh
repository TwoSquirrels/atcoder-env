#!/bin/bash

set -eo pipefail

cd "$(dirname "$0")/../"

dedup() { awk '{ print (/^#include/ && seen[$0]++) ? "" : $0 }'; }

risundle main.cc | dedup > submission.cc
if ! g++ -std=gnu++20 -fsyntax-only submission.cc; then
  echo "[bundle/cpp20.sh] falling back to --no-tree-shaking" >&2
  risundle --no-tree-shaking main.cc | dedup > submission.cc
  g++ -std=gnu++20 -fsyntax-only submission.cc
fi
cat submission.cc
