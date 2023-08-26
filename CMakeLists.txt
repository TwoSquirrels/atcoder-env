cmake_minimum_required(VERSION 3.25)
project(atcoder_env)

set(CMAKE_CXX_STANDARD 20)

add_executable(atcoder_env main.cpp)

include_directories($ENV{CPLUS_INCLUDE_PATH})
target_compile_options(atcoder_env PUBLIC -Wall -Wextra -g -DDEBUG)
