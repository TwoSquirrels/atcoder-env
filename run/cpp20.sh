#!/bin/bash

cd "$(dirname "$0")/../"

mkdir ./cmake-build-debug
cd ./cmake-build-debug/

cmake .. && make && ./atcoder_env
