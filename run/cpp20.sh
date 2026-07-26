#!/bin/bash

set -e

cd "$(dirname "$0")/../"

mkdir -p ./cmake-build-debug
cd ./cmake-build-debug/

cmake .. && make && ./atcoder_env
