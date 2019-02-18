#include "STBDC.hh"


STBDC::STBDC()
{
  run = 0;
  SnA = 0;

  Clear();
}
STBDC::STBDC(const STBDC &cp)
{
  run      = cp.run;
  evt      = cp.evt;
  SnA      = cp.SnA;
  beamPID  = cp.beamPID;

  aoq      = cp.aoq;
  z        = cp.z;
  tof      = cp.tof;
  beta     = cp.beta;
  brho     = cp.brho;
  isGood   = cp.isGood;
  intZ     = cp.intZ;
  intA     = cp.intA;

  bdcax    = cp.bdcax;
  bdcby    = cp.bdcby;
  ProjX    = cp.ProjX;
  ProjY    = cp.ProjY;
  ProjZ    = cp.ProjZ;
  ProjP    = cp.ProjP;
  ProjPX   = cp.ProjPX;
  ProjPY   = cp.ProjPY;
  ProjPZ   = cp.ProjPZ;
  ProjA    = cp.ProjA;
  ProjB    = cp.ProjB;
}

STBDC &STBDC::operator=(const STBDC &cp)
{
  if( this != &cp )
    *this = cp;

  return *this;
}

void STBDC::Clear(){
  evt      = 0;
  beamPID  = 0;
  aoq      = 0.;
  z        = 0.;
  tof      = 0.;
  beta     = 0.;
  brho     = 0.;
  isGood   = 0.;;
  intZ     = 0.;
  intA     = 0.;

  bdcax    = 0.;
  bdcby    = 0.;
  ProjX    = 0.;
  ProjY    = 0.;
  ProjZ    = 0.;
  ProjP    = 0.;
  ProjPX   = 0.;;
  ProjPY   = 0.;;
  ProjPZ   = 0.;;
  ProjA    = 0.;
  ProjB    = 0.;

}

ClassImp(STBDC);
