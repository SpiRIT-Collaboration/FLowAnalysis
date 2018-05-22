#ifndef OPENFLW_H
#define OPENFLW_H

const UInt_t nconfig = 4;
UInt_t  isys[]            = {10, 10, 10, 10};
TString sysName[nconfig]; //{"132Sn","108Sn","124Sn","112Sn"};

TChain *rChain[nconfig];

TString printHeader="";
TString aRun[nconfig];
TString sDB[nconfig];


UInt_t m_bgn = 0;
UInt_t m_end = 1;

#endif
