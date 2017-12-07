#ifndef OPENFLW_H
#define OPENFLW_H

UInt_t  isys[]    = {10, 10, 10, 10};
TString sysName[] = {"132Sn","108Sn","124Sn","112Sn"};

const UInt_t nconfig = 4;
TChain *LChain[nconfig];

TString printHeader="";
TString aRun[nconfig];
TString sDB[nconfig];

#endif
