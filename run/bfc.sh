#!/bin/bash

cd "$(dirname "$0")/../"

[[ -r BrainFuckCompiler.py ]] || wget 'https://gist.github.com/KartikTalwar/4583510/raw/9a3adf8aeb75ad29ddb6f3431fd04a88ad447f31/BrainFuckCompiler.py'

python2 BrainFuckCompiler.py main.bf main.c && mv main.c main.cc && ./run/cpp20.sh
echo ''
echo '[INFO] Finish excution.'
