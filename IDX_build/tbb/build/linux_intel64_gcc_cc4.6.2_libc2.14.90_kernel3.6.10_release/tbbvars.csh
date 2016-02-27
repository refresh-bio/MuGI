#!/bin/csh
setenv TBBROOT "/mnt/syn2/adanek/MuGI/MuGI/IDX_build/tbb" #
setenv tbb_bin "/mnt/syn2/adanek/MuGI/MuGI/IDX_build/tbb/build/linux_intel64_gcc_cc4.6.2_libc2.14.90_kernel3.6.10_release" #
if (! $?CPATH) then #
    setenv CPATH "${TBBROOT}/include" #
else #
    setenv CPATH "${TBBROOT}/include:$CPATH" #
endif #
if (! $?LIBRARY_PATH) then #
    setenv LIBRARY_PATH "${tbb_bin}" #
else #
    setenv LIBRARY_PATH "${tbb_bin}:$LIBRARY_PATH" #
endif #
if (! $?LD_LIBRARY_PATH) then #
    setenv LD_LIBRARY_PATH "${tbb_bin}" #
else #
    setenv LD_LIBRARY_PATH "${tbb_bin}:$LD_LIBRARY_PATH" #
endif #
 #
