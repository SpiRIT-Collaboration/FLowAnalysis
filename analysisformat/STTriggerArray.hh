#ifndef STTRIGGERARRAY_H
#define STTRIGGERARRAY_H

#include "vector"
#include <TObject.h>
#include <iostream>


class STTriggerArray : public TObject
{

public:
  Int_t run;
  Int_t evt;
  Long_t time;

  //KYOTO 
  Int_t knhit;
  Int_t adch[64];    // ADC high gain of each ch
  Int_t adcl[64];    // ADC low gain of each ch
  Int_t tdclfirst[64]; // first hit TDC leading
  Int_t tdctfirst[64]; // first hit TDC trailing
  Int_t scr[64];     // Scaler of each ch
  Int_t scORU, scORL, scOR64; // Scaler of each EASIROC, and sum of them
  Int_t mhit[64];    // multi hit information (multi hit TDC)

  // hit paddle information
  std::vector<Int_t> kch;
  std::vector<Int_t> kah;       
  std::vector<Int_t> kal;       
  std::vector<Int_t> ktL;       // MHTDC Leading edge
  std::vector<Int_t> ktT;       // MHTDC Trailing edge
  std::vector<Int_t> kscinum;   // KyotoArray Scinti. number
  std::vector<Double_t> kxpos;
  std::vector<Double_t> kzpos;

  //TriggerBit
  Int_t katnhit;
  std::vector<Double_t> katxpos;
  std::vector<Double_t> katzpos;
  std::vector<Int_t> bitpat;
  Bool_t ftboxmissing;
  unsigned int tboxbitdata; //!
  unsigned int tboxeventdata; //!

  //RPV130
  Int_t rpvnhit;
  std::vector<Int_t> rpvbitpat;

public:
  STTriggerArray(){};
  virtual ~STTriggerArray(){};
  
  void ClearKyotoArray();

  ClassDef(STTriggerArray,1);
};

#endif
