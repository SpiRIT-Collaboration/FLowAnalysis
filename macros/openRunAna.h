#ifndef OPENFLW_H
#define OPENFLW_H

TChain *rChain;

TString printHeader="";
TString aRun;
TString sSuf;
TString sVer;
TString dVer;
TString oVer;

UInt_t  isys;
TString sysName;

// Retrivew tree
TClonesArray *aArray; 
TClonesArray *aFlowArray;
TClonesArray *aNLClusterArray;

Double_t RPPsi;
#endif
