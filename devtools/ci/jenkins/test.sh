#!/bin/bash -ex

python devtools/run-ctest.py --job-duration=120 --timeout 600 --in-order $*