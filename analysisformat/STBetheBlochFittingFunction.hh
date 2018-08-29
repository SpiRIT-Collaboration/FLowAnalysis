#ifndef BETHEBLOCHFITTER_H
#define BETHEBLOCHFITTER_H

#include <TMath.h>
#include <TF1.h>
#include <TROOT.h>
#include "Math/RootFinder.h"
#include "Math/WrappedTF1.h"


namespace BetheBlochFitter
{
  const Double_t kmpi  = 139.57018;
  const Double_t kmp   = 938.2720813;
  const Double_t kmd   = 1875.612762;
  const Double_t kmt   = 2808.921112;
  const Double_t kmhe3 = 2808.39132;
  const Double_t kmal  = 3727.379378;
  const Double_t kme   = 0.5109989461;

  const Double_t I_Ar  = 285.*1.e-6;
  const Double_t I_CH4 = 131.*1.e-6;


  Double_t fddedx(Double_t *x, Double_t *p);
  
  Double_t fdedx(Double_t Z, Double_t m, Double_t *x, Double_t *p);
  
  inline Double_t fdedx_pim (Double_t *x, Double_t *p) { return fdedx(-1, kmpi,  x, p); }
  inline Double_t fdedx_pi  (Double_t *x, Double_t *p) { return fdedx( 1, kmpi,  x, p); }
  inline Double_t fdedx_p   (Double_t *x, Double_t *p) { return fdedx( 1, kmp,   x, p); }
  inline Double_t fdedx_d   (Double_t *x, Double_t *p) { return fdedx( 1, kmd,   x, p); }
  inline Double_t fdedx_t   (Double_t *x, Double_t *p) { return fdedx( 1, kmt,   x, p); }
  inline Double_t fdedx_he3 (Double_t *x, Double_t *p) { return fdedx( 2, kmhe3, x, p); }
  inline Double_t fdedx_al  (Double_t *x, Double_t *p) { return fdedx( 2, kmal,  x, p); }

  Double_t CalcMass(Double_t* p);

}

#endif
