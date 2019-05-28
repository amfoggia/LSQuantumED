#!/bin/bash

# --------------------------------------- #
#    TEST: run the specified executable   #
#      and compare the obtained result    #
#          with the reference one.        #
#                                         #
#                                         #
#           ARGS = c++ executable         #
#                  reference file(.ref)   #
#                  output file(.dat)      #
# --------------------------------------- #


EXE=$1; shift
_nprocs=$1; shift
_out_ref=$1; shift
_out_dat=$1; shift
_meson_args=$1

# Run the program
echo "Running for" ${EXE}
mpiexec -np ${_nprocs} ${EXE} ${_meson_args}

# Test if run was successfull
if [[ $? -ne 0 ]]
then
    echo "run failed"
    exit 1
fi

# Compare files
echo "Comparing files"
diff ${_out_ref} ${_out_dat}

# test if diff was successfull
if [[ $? -ne 0 ]]
then
    echo "comparison failed"
    exit 1
fi

rm ${_out_dat}
