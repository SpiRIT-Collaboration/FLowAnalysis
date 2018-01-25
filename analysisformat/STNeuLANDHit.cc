#include "STNeuLANDHit.hh"

ClassImp(STNeuLANDHit)

STNeuLANDHit::STNeuLANDHit()
{
  Clear();
}

void STNeuLANDHit::Clear(Option_t *)
{
  edep = -1;
  tof = -1;
  beta = -1;

  localx = -1; localy = -1; localz = -1;
}
