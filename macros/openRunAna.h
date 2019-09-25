#ifndef OPENFLW_H
#define OPENFLW_H

UInt_t  isys;
TString sysName;

TChain *rChain;

TString printHeader="";
TString aRun;
TString sSuf;
TString sVer;
TString dVer;

// Retrivew tree
TClonesArray *aArray; 
TClonesArray *aFlowArray;
TClonesArray *aNLClusterArray;


#endif
