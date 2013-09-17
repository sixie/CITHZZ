
###########################################################
# Do Fits For Regression Energy
###########################################################


#53X. Eta bins
# for bin in 0 1 2 3
# do
#     for r9bin in 0 1
#     do
# 	root -l -b -q  CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/scale/FitInclusiveZ.C+\(\"/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_MC_53X.root\",\"mc_EtaBin${bin}_R9Bin${r9bin}\",$bin,$r9bin,1\)
# 	root -l -b -q  CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/scale/FitInclusiveZ.C+\(\"/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_Data_2012.root\",\"data2012_Jan22_EtaBin${bin}_R9Bin${r9bin}\",$bin,$r9bin,0\)
#     done
# done

#42X. Eta bins
for bin in 0 1 2 3
do
    for r9bin in 0 1
    do
	root -l -b -q  CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/scale/FitInclusiveZ.C+\(\"/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_MC_42X.root\",\"mc42X_EtaBin${bin}_R9Bin${r9bin}\",$bin,$r9bin,1\)
	root -l -b -q  CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/scale/FitInclusiveZ.C+\(\"/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_Data_2011.root\",\"data2011_EtaBin${bin}_R9Bin${r9bin}\",$bin,$r9bin,0\)
    done
done

