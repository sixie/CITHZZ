// $Id: rootlogon.C,v 1.1 2012/07/18 19:23:07 sixie Exp $

{

  if (gSystem->Getenv("CMSSW_VERSION")) {

    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I$ROOFITSYS/include/";
    gSystem->SetIncludePath(path.Data());
 
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }
  }

//   gSystem->Load("libFWCoreFWLite.so"); AutoLibraryLoader::enable();
//   gSystem->Load("libDataFormatsFWLite.so");
//   gSystem->Load("libHiggsAnaHZZ4l.so");
//   gSystem->Load("libWWAnalysisAnalysisStep.so");



}
