#ifndef STLORENTZBOOSTVECTOR_HH
#define STLORENTZBOOSTVECTOR_HH
#include <TVector3.h>

namespace STLorentzBoostVector
{
  const Double_t amu  = 931.4940954; //MeV/c2
  const Double_t c    = 299792458.; //m/s 
  const Double_t me   = 5.48579909065e-04; //amu 

  const TString system[] = {"(132Sn + 124Sn)", "(108Sn + 112Sn)", "(124Sn + 112Sn)", "(112Sn + 124Sn)", "(p + p)132"   ,"(p + p)108 "  ,"(p + p)124 ","(p + p)112 "};
  const Double_t AB[]    = {  132.           ,  108.            ,   124.           ,   112.           ,  1.            , 1             ,1            ,1};
  const Double_t mB[]    = {  131.917821719  ,  107.911892833   ,   123.905273581  ,   111.904821807  ,  1.00727646688 , 1.00727646688 ,1.00727646688,1.00727646688};
  const Double_t eB_lb[] = {  268.3          ,  268.3           ,   269.8          ,   269.6          ,  268.3         , 268.3         ,269.8        ,269.6        };
  const Double_t mT[]    = {  123.905273581  ,  111.904821807   ,   111.904821807  ,   123.905273581  ,  1.00727646688 , 1.00727646688 ,1.00727646688,1.00727646688};
  const Double_t ne[]    = {    50.          ,   50.            ,    50.           ,    50.           ,  0.            , 0.            ,0.           ,0            };
    

  inline TVector3 GetBoostVector(UInt_t sysid);
};

TVector3 STLorentzBoostVector::GetBoostVector(UInt_t sysid = 4)
{
  
  if( sysid >= 8 ) sysid = 0; // For Simulation

  Double_t beamMass  =(mB[sysid] - ne[sysid]*me);
  Double_t EkB_lb    = eB_lb[sysid]  * beamMass;
  beamMass *= amu;
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
