#!/bin/bash

g++ -std=gnu++20 -Wall -Wextra \
    -g -fsanitize=address \
    -DDEBUG "$1" && ./a.out
