#!/bin/bash
#

# Check if command line argument --debug is given
if [ "$1" == "--debug" ]; then
    if [ -d build_dbg ]; then
      rm -rf build_dbg
    fi
       
    meson setup build_dbg --buildtype=debug
    meson compile -C build_dbg -j$(nproc)

    if [ "$2" == "--run" ]; then
      gdb --args ./build_dbg/src/solver/laneEmden -ne 50 -o 1 -n 1
    fi
    exit 0
fi

if [ -d build ]; then
  rm -rf build
fi
meson setup build
meson compile -C build -j$(nproc)

if [ "$1" == "--run" ]; then
  ./build/src/solver/laneEmden -ne 50 -o 1 -n 1
fi
