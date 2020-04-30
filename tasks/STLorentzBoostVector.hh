#ifndef STLORENTZBOOSTVECTOR_HH
#define STLORENTZBOOSTVECTOR_HH
#include <TVector3.h>

namespace STLorentzBoostVector
{
  const Double_t amu  = 931.4940954; //MeV/c2
  const Double_t c    = 299792458.; //m/s 
    
  const TString system[5] = {"(132Sn + 124Sn)", "(108Sn + 112Sn)", "(124Sn + 112Sn)", "(112Sn + 124Sn)", "(p + p)"};
  const Double_t AB[5]    = {  132.           ,  108.            ,   124.           ,   112.           ,  1.      };
  const Double_t mB[5]    = {  131.8906       ,  107.8844964     ,   123.8778449    ,   111.8773895    ,  1.00727646688};
  const Double_t eB_lb[5] = {  268.9          ,  268.9           ,   270.2          ,   270.2          ,  268.9   };
  const Double_t mT[5]    = {  123.8773895    ,  111.8773895     ,   111.8773895    ,   123.8773895    ,  1.00727646688};


  inline TVector3 GetBoostVector(UInt_t sysid);
};

TVector3 STLorentzBoostVector::GetBoostVector(UInt_t sysid = 4)
{
  if( sysid == 5 ) sysid = 0; // For Simulation

  Double_t EkB_lb    = eB_lb[sysid]  * mB[sysid];
  Double_t beamMass  = mB[sysid] * amu;
  Double_t targetMass= mT[sysid] * amu;

  // Beam                                                                                                                                        
  Double_t EB_lb  = EkB_lb + beamMass;
  Double_t PB_lb  = sqrt(EB_lb*EB_lb - beamMass * beamMass);

  TLorentzVector bmVec = TLorentzVector( TVector3(0., 0., PB_lb), EB_lb );
  TLorentzVector tgVec = TLorentzVector( TVector3(0., 0., 0.), targetMass );
  auto totalVec = bmVec + tgVec;

  return totalVec.BoostVector();
}

#endif
