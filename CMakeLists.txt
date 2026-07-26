cmake_minimum_required(VERSION 3.25)
project(atcoder_env CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# main.cc は追跡対象外なので、無い場合はテンプレートから用意する
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/main.cc)
  configure_file(templates/cpp20.cc ${CMAKE_SOURCE_DIR}/main.cc COPYONLY)
endif()

add_executable(atcoder_env main.cc)

# CPLUS_INCLUDE_PATH はシェルの設定なので、GUI から起動した IDE には渡らない
target_include_directories(atcoder_env PRIVATE
  libs/risu
  "$ENV{CPLUS_INCLUDE_PATH}"
  "$ENV{HOME}/ac-library")
# CMake は絶対パスでコンパイルするため、dump の出力する位置情報も絶対パスになる。
# 接頭辞を削ってリポジトリからの相対パスで表示させる。
target_compile_options(atcoder_env PRIVATE
  -Wall -Wextra -g
  "-ffile-prefix-map=${CMAKE_SOURCE_DIR}/=")
target_compile_definitions(atcoder_env PRIVATE DEBUG)

# ヘッダオンリーのライブラリは main.cc を書き換えるたびに全て再解析される。
# 解析済みの状態を保存しておくことで、ライブラリを触らない限り再解析を省く。
target_precompile_headers(atcoder_env PRIVATE "${CMAKE_SOURCE_DIR}/pch.hpp")
