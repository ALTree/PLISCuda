#!/bin/bash

nvccPath=$(which nvcc)
$nvccPath $* -std=c++11