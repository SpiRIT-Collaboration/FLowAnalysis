#include "STFlowInfo.hh"

STFlowInfo::STFlowInfo()
{
  SnA = 0;
  run = 0;

  Clear();
}

STFlowInfo::STFlowInfo(const STFlowInfo &cp)
{
  Clear();

  run      = cp.run;
  evt      = cp.evt;
  SnA      = cp.SnA;
  beamPID  = cp.beamPID;


  for(UInt_t i = 0; i < 7; i++ )
    ntrack[i] = cp.ntrack[i];

  mtrack0   = cp.mtrack0;
  mtrack1   = cp.mtrack1;
  mtrack2   = cp.mtrack2;
  mtrack3   = cp.mtrack3;
  mtrack4   = cp.mtrack4;
  mtrack5   = cp.mtrack5;
  mtrack6   = cp.mtrack6;

  unitP    = cp.unitP;
  unitP_fc = cp.unitP_fc;
  unitP_rc = cp.unitP_rc;
  unitP_1  = cp.unitP_1;
  unitP_2  = cp.unitP_2;
  unitP_1fc= cp.unitP_1fc;
  unitP_2fc= cp.unitP_2fc;
  mtrack_1 = cp.mtrack_1;
  mtrack_2 = cp.mtrack_1;
  
  for(UInt_t i = 0; i < 3; i++) {
    bsPhi[i]   = cp.bsPhi[i];
    bsPhi_1[i] = cp.bsPhi_1[i];
    bsPhi_2[i] = cp.bsPhi_2[i];
  }

  rpSigma   = cp.rpSigma;
  rpChi[0]  = rpChi[0];
  rpChi[1]  = rpChi[1];
}

STFlowInfo &STFlowInfo::operator=(const STFlowInfo &cp)
{
  if( this != &cp )
    *this = cp;

  return *this;
}


void STFlowInfo::SetNTrack(UInt_t *nval) {
  for(UInt_t i = 0; i < 7; i++ ) {
    ntrack[i] = *(nval + i);
  }
  
  mtrack0 = ntrack[0];
  mtrack1 = ntrack[1];
  mtrack2 = ntrack[2];
  mtrack3 = ntrack[3];
  mtrack4 = ntrack[4];
  mtrack5 = ntrack[5];
  mtrack6 = ntrack[6];

}

void STFlowInfo::SetNTrack(UInt_t nval, UInt_t idx) {

  if(idx < 7) ntrack[idx] = nval;

  switch( idx ){
  case 0:
    mtrack0 = nval;
    break;
  case 1:
    mtrack1 = nval;
    break;
  case 2:
    mtrack2 = nval;
    break;
  case 3:
    mtrack3 = nval;
    break;
  case 4:
    mtrack4 = nval;
    break;
  case 5:
    mtrack5 = nval;
    break;
  case 6:
    mtrack6 = nval;
    break;
  }
}

void STFlowInfo::Clear(){
  
  evt = 0;
  beamPID = 0;
  
  for(UInt_t i = 0; i < 7; i++)
    ntrack[i] = 0;

  mtrack0   = 0;
  mtrack1   = 0;
  mtrack2   = 0;
  mtrack3   = 0;
  mtrack4   = 0;
  mtrack5   = 0;
  mtrack6   = 0;

  unitP_fc   = TVector3(0.,0.,0.);
  unitP_rc   = TVector3(0.,0.,0.);

  unitP      = TVector3(0.,0., 0.);
  unitP_1    = TVector3(0.,0., 0.);
  unitP_2    = TVector3(0.,0., 0.);
  unitP_1fc  = TVector3(0.,0., 0.);
  unitP_2fc  = TVector3(0.,0., 0.);

  mtrack_1 = 0;
  mtrack_2 = 0;

  for(UInt_t i = 0; i < 3; i++) {
    bsPhi[i]  = 0.;
    bsPhi_1[i] = 0.;
    bsPhi_2[i] = 0.;
  }

  rpSigma  = 0.;
  rpChi[0] = 0.;
  rpChi[1] = 0.;
  
}

ClassImp(STFlowInfo);
