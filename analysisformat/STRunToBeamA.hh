#ifndef STRUNTOBEAMA_H
#define STRUNTOBEAMA_H

namespace STRunToBeamA
{
  // system selection                                                                                                                                               
  // 0: 132Sn + 124Sn : 2841 ~ 3039                                                                                                                                 
  // 1: 108Sn + 112Sn : 2261 ~ 2509                                                                                                                                 
  // 2: 124Sn + 112Sn : 2520 - 2653                                                                                                                                 
  // 3: 112Sn + 124Sn                                                                                                                                              

  const  UInt_t  beam_ID[] = {    0,   1,   2,   3,   4};
  const  UInt_t  beam_A[]  = {  132, 108, 124, 112, 100};
  
  inline UInt_t  GetBeamA(Int_t irun);
  inline TString GetBeamSnA(Int_t irun);

  inline UInt_t GetSystemID(Int_t irun);
};

UInt_t STRunToBeamA::GetBeamA(Int_t irun) 
{
  auto id = GetSystemID(irun);
  return beam_A[id];
}

TString STRunToBeamA::GetBeamSnA(Int_t irun)
{
  auto id = GetSystemID(irun);
  return   Form("%dSn",beam_A[id]);
}
UInt_t STRunToBeamA::GetSystemID(Int_t irun)
{
  UInt_t id = 4;
  if(irun >= 2841 && irun <= 3039)
    id = 0;
  else if(irun >= 2261 && irun <= 2509)
    id = 1;
  else if(irun >= 3059 && irun <= 3184)
    id = 2;
  else if(irun >= 2520 && irun <= 2653)
    id = 3;

  return id;
}

#endif


