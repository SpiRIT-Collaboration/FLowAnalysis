#include "STTriggerArray.hh"

void STTriggerArray::ClearKyotoArray(){

  for(Int_t i = 0; i < 64; i++){
    adch[i] = 0;
    adcl[i] = 0;
    scr[i]  = 0;
    mhit[i] = 0;
    tdclfirst[i]=0;
    tdctfirst[i]=0;
  }
  scORU = 0;
  scORL = 0;
  scOR64 = 0;

  knhit = 0;
  kch.clear();
  kah.clear();
  kal.clear();
  ktL.clear();
  ktT.clear();

  kzpos.clear();
  kxpos.clear();
  kscinum.clear();

}

#if !defined(__CINT__)
ClassImp(STTriggerArray);
#endif
