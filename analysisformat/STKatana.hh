#ifndef STKatana_h
#define STKatana_h

#include "vector"
#include <TObject.h>

class STKatanaSignal : public TObject
{
 public:
  int Module;
  int Channel;
  std::vector<int> Amplitude;
  std::vector<int> Tstart;
  std::vector<int> Tstop;  
  std::vector<int> Tmax;
  float Pedestal;
  int Npeak;
  //  std::vector<float> peak_ampl;
  //  std::vector<float> peak_time;
  /* void AddPeak(float ampl, float tstart, float tstop, float tmax); */
  /* void Dump(); */
  STKatanaSignal(){};
  ClassDef(STKatanaSignal,3);
};

class STTriggerBox : public TObject
{
 public:
  unsigned long evnum;
  unsigned long bitpat;
  unsigned long time_stamp;
  int KatanaM;
  int MinBias;
  int VETO;
  int AC;
  int KYOTO;
  int offset;
  std::vector<int> bitpattern;
  //void AddBitPattern(int bit);
  STTriggerBox(){};
  //void Reset();
  //void SetUp();
  ClassDef(STTriggerBox,1);
};


class STKatana : public TObject
{
 public:
  int run_number;
  unsigned long event_number;
  unsigned long event_size;
  unsigned long time_stamp;
  //unsigned long time_stamp1;
  float max_veto;
  int mult; //! 
  std::vector<STKatanaSignal> signal;
  void AddSignal(STKatanaSignal & sig){
    signal.push_back(sig);
    mult++;
  };
  void Reset() {signal.clear(); mult=0; max_veto=0;}; //!   
  // void Dump();
  STKatana();
  ClassDef(STKatana,0);
};

#endif
