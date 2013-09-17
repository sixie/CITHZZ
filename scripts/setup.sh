# Packages needed to be checked out


# Combination Tools: need it for roofit signal model
    addpkg   HiggsAnalysis/CombinedLimit   HEAD; 

# EGamma Electron MVA Class
    cvs co -r V00-00-07 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools                    ;
    cvs update -rV00-00-11 EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h ;

# Muon Effective Area Class
    cvs co -r V00-00-10 -d Muon/MuonAnalysisTools UserCode/sixie/Muon/MuonAnalysisTools                      ;

