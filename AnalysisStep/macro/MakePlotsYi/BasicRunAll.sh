#!/bin/bash

./Run 1125 8
./Run 8125 8
./Run 100 8
./Run 101 8
./Run 102 8

for i in -1 0 `seq 100 107`
do
   ./Run $i 7
   ./Run $i 8
done

hadd -f All_Plots8TeV.root Plots8TeV_10[0-7].root
hadd -f All_Plots7TeV.root Plots7TeV_10[0-7].root

./Run 1125 7
./Run 1125 8
./Run 2125 8
./Run 3125 8

./Run 8125 8
./Run 9125 8


cp Plots7TeV_1125.root All_Signal7.root
hadd -f All_Signal8.root Plots8TeV_1125.root Plots8TeV_2125.root Plots8TeV_3125.root
cp Plots8TeV_8125.root All_PSSignal8.root
cp Plots8TeV_9125.root All_SpinTwoMinimalSignal8.root

cp All_Plots8TeV.root All_Signal8.root Plots8TeV_ZX.root Plots8TeV_DATA.root All_SpinTwoMinimalSignal8.root All_PSSignal8.root Plotting/
cp All_Plots7TeV.root All_Signal7.root Plots7TeV_ZX.root Plots7TeV_DATA.root Plotting/


