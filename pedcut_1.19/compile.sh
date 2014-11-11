#!/bin/bash

rm -r -f bin
mkdir bin
cd src/PEDIG-cut
rm -f par3
rm -f ped_util
sh make.sh
cd ../PedHunter-cut
rm -f asp
rm -f *.o
make
cd ../..
cp src/PEDIG-cut/par3 bin/.
cp src/PedHunter-cut/asp bin/.
cp src/PedCut_v1.19 bin/.
cp src/pedsum.pl bin/.
chmod +x bin/*
