#ifndef OPENFLW_H
#define OPENFLW_H

TChain *rChain;

TString printHeader="";
TString sRun;
TString sSuf;
TString sVer;
TString dVer;
TString oVer;

UInt_t  iRun;
UInt_t  isys;
TString sysName;
Long64_t nEntry;

// Retrivew tree
TClonesArray   *aArray; 
TClonesArray   *aFlowArray;
TClonesArray   *aNLClusterArray;
TClonesArray   *aBDC;
STBDC          *aBeamInfo;
STFlowInfo     *aFlowInfo;
UInt_t          beamPID;

Double_t RPPsi;

void openRunAna(Bool_t bopen=kTRUE);
void OpenChain();
Long64_t SetBranch();
void ShowProcess(Long64_t ievt);
void SaveCanvas(TString fopt = "", Int_t isel=-1);

#include "../flowformat/STFlowInfo.hh"
//#include "../flowformat/STBDC.hh"

#endif
