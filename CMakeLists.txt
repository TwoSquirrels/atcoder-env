cmake_minimum_required(VERSION 3.28)
project(atcoder_env CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT EXISTS ${CMAKE_SOURCE_DIR}/main.cc)
    configure_file(templates/cpp20.cc ${CMAKE_SOURCE_DIR}/main.cc COPYONLY)
endif()

add_executable(atcoder_env main.cc)

target_include_directories(atcoder_env PRIVATE
    libs/risu
    "$ENV{CPLUS_INCLUDE_PATH}"
)
target_compile_options(atcoder_env PRIVATE
    -Wall -Wextra -g
    "-ffile-prefix-map=${CMAKE_SOURCE_DIR}/="
)
target_compile_definitions(atcoder_env PRIVATE DEBUG)

target_precompile_headers(atcoder_env PRIVATE
    <bits/stdc++.h>
    <atcoder/all>
    <risu/all.hpp>
)
