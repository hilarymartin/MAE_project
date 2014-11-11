#!/bin/bash

rm -f par3
g77 -O2 -I. par3.f -o par3
gfortran -O2 -I. par3.f -o par3

rm -f ped_util
g77 -O2 -I. ped_util.f -o ped_util
gfortran -O2 -I. ped_util.f -o ped_util
