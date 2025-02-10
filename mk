#!/bin/bash
#

# Check if command line argument --debug is given
if [ "$1" == "--debug" ]; then
    if [ -d build_dbg ]; then
      rm -rf build_dbg
    fi
       
    meson setup build_dbg --buildtype=debug
    meson compile -C build_dbg -j$(nproc)
    exit 0
fi

if [ -d build ]; then
  rm -rf build
fi
meson setup build
meson compile -C build -j$(nproc)

if [ "$1" == "--run" ]; then
  ./build/src/solver/laneEmden -ne 25 -o 1
fi
