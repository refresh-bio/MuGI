#!/bin/bash
export TBBROOT="/mnt/syn2/adanek/MuGI/MuGI/IDX_build/tbb" #
tbb_bin="/mnt/syn2/adanek/MuGI/MuGI/IDX_build/tbb/build/linux_intel64_gcc_cc4.6.2_libc2.14.90_kernel3.6.10_debug" #
if [ -z "$CPATH" ]; then #
    export CPATH="${TBBROOT}/include" #
else #
    export CPATH="${TBBROOT}/include:$CPATH" #
fi #
if [ -z "$LIBRARY_PATH" ]; then #
    export LIBRARY_PATH="${tbb_bin}" #
else #
    export LIBRARY_PATH="${tbb_bin}:$LIBRARY_PATH" #
fi #
if [ -z "$LD_LIBRARY_PATH" ]; then #
    export LD_LIBRARY_PATH="${tbb_bin}" #
else #
    export LD_LIBRARY_PATH="${tbb_bin}:$LD_LIBRARY_PATH" #
fi #
 #
