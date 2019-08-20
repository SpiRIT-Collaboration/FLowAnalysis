#ifndef STSOURCE_HH
#define STSOURCE_HH
#include "FairSource.h"
#include "FairLogger.h"
#include "TChain.h"
#include "TString.h"

class STSource : public FairSource
{
public:
  STSource();
  ~STSource(){};

public:
  static STSource* Instance();

  UInt_t  SetInputChain(TChain *fchain);
  TChain* GetInputChain(TString fcname);

private:
  static STSource* fgRinstance;

  std::map<std::string, TChain*> fInputChain;
 

  FairLogger *fLogger;
};
#endif
