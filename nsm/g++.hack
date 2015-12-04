#!/bin/bash

# g++ Hack
#

# in CUDA 6.0 the source code is always the last parameter
SourceFile="${@: -1}"

# get the file extension
Extension=${SourceFile##*.}

if [ "$Extension" == "hpp" ]
then
   StdFlag="-std=c++11"
else
   StdFlag=""
fi

# run now the g++ 4.8 in your own path with the personalized std option
/usr/bin/gcc-4.8 $StdFlag $*
