#! /bin/bash
# example of usage: 
# ./fitAll.sh $PWD/8TeV 2012 2 0 1
# $PWD is important because in that way the absolute path of the root file is put in the combined datacard and the job can be run in batch

# combine the 7 TeV 2D cards: ./fitAll.sh $PWD/7TeV 2011 2 1 0  
# combine the 8 TeV 2D cards: ./fitAll.sh $PWD/8TeV 2012 2 1 0
# combine the 7 + 8 2D cards: ./fitAll.sh $PWD 0 2 1 0

# submit the fits for 7 + 8 TeV: ./fitAll.sh $PWD 0 2 0 1

DIR=$1
yiear=$2
fitDim=$3
combineCards=$4
runFit=$5

pwd=$PWD

if [ $yiear == 2011 ]
    then
    let TeVstr=7
elif [ $yiear == 2012 ]
    then
    let TeVstr=8
else 
    let TeVstr=78
fi

echo "Running for $TeVstr TeV"

masses="`seq 110 1 160` `seq 162 2 290` `seq 295 5 350` `seq 360 10 400` `seq 420 20 1000`"

COMBDIR=$DIR/combcards
mkdir -p $COMBDIR

#first combined the cards 
if [ $combineCards == 1 ]
    then
    if [[ ( $yiear == 2011 ) || ( $yiear == 2012 ) ]]; then
        for m in $masses; do
            echo "Combining card for mass = $m GeV"
            combineCards.py hzz4l_4mu=$DIR/card_$fitDim\D_m$m\_$TeVstr\TeV_4mu.txt \
                hzz4l_4e=$DIR/card_$fitDim\D_m$m\_$TeVstr\TeV_4e.txt \
                hzz4l_2e2mu=$DIR/card_$fitDim\D_m$m\_$TeVstr\TeV_2e2mu.txt \
                > $COMBDIR/card_$fitDim\D_m$m\_$TeVstr\TeV_4l.txt
        done
    else
        for m in $masses; do
             echo "Combining 7 TeV and 8 TeV card for mass = $m GeV"
             combineCards.py \
                 hzz4l_7TeV_4mu=$DIR/7TeV/card_$fitDim\D_m$m\_7TeV_4mu.txt \
                 hzz4l_7TeV_4e=$DIR/7TeV/card_$fitDim\D_m$m\_7TeV_4e.txt \
                 hzz4l_7TeV_2e2mu=$DIR/7TeV/card_$fitDim\D_m$m\_7TeV_2e2mu.txt \
                 hzz4l_8TeV_4mu=$DIR/8TeV/card_$fitDim\D_m$m\_8TeV_4mu.txt \
                 hzz4l_8TeV_4e=$DIR/8TeV/card_$fitDim\D_m$m\_8TeV_4e.txt \
                 hzz4l_8TeV_2e2mu=$DIR/8TeV/card_$fitDim\D_m$m\_8TeV_2e2mu.txt \
                 > $COMBDIR/card_$fitDim\D_m$m\_78TeV_4l.txt
        done
    fi
fi

RMIN=-2
RMAX=2

if [ $runFit == 1 ] 
    then
    for m in $masses; do
        if [[ ( $m < 118 ) || ( ( $m > 160 && $m < 190 ) ) ]]; then
            let RMIN=-5
        fi;
        bsub -q 8nh $PWD/combine-one.sh  -M MaxLikelihoodFit -m $m --rMin $RMIN --rMax 2 $COMBDIR/card_$fitDim\D_m$m\_$TeVstr\TeV_4l.txt ;
    done;
fi;
