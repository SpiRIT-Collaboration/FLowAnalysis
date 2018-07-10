TChain *fChain;
//
// How to use command
// RUN=2900,2901,2902.2903 STTPCDIR=/cache/scr/spirit/recoData/20180706/132Sn124Sn/ STVERSION=develop.RC.1603.e284ae9 root Open_runreco.C

void Open_runreco()
{
  fChain = new TChain("cbmsim");
 
  TString rootDir = gSystem -> Getenv("STTPCDIR");
  TString fileversion = gSystem->Getenv("STVERSION");
  TString lRun = gSystem -> Getenv("RUN");
  Int_t   nrun = lRun.Length()/5 + 1;
   
  Int_t ist = 0;
  for(UInt_t j = 0; j < nrun; j++){
    if(ist > lRun.Length() - 2) break;

    TString sRun = lRun(ist, 4);

    Int_t i = 0;
    while(kTRUE && fileversion != ""){

      TString recoFile = Form("run"+sRun+"_s%d.reco."+fileversion+".root",i);
      std::cout << " recoFile " << rootDir+recoFile << std::endl;
      
      if(gSystem->FindFile(rootDir,recoFile)){
	fChain -> Add(recoFile);
      }
      else
	break;
      i++;
      
    }

    ist += 5;
  }
}
