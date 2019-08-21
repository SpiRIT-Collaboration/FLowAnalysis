#ifndef KatanaRoot_h
#define KatanaRoot_h

#include "vector"
#include <TObject.h>

class Signal : public TObject
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
  Signal(){};
  ClassDef(Signal,0);
};

class TriggerBox : public TObject
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
  TriggerBox(){};
  //void Reset();
  //void SetUp();
  ClassDef(TriggerBox,2);
};


class KatanaRoot : public TObject
{
 public:
  int run_number;
  unsigned long event_number;
  unsigned long event_size;
  unsigned long time_stamp;
  //unsigned long time_stamp1;
  float max_veto;
  int mult; //! 
  std::vector<Signal> signal;
  void AddSignal(Signal & sig){
    signal.push_back(sig);
    mult++;
  };
  void Reset() {signal.clear(); mult=0; max_veto=0;}; //!   
  // void Dump();
  KatanaRoot();
  ClassDef(KatanaRoot,0);
};

#endif
