#include "AMDParticle.hh"

ClassImp(AMDParticle);

AMDParticle::AMDParticle()
: TObject(), 
  fPdg(0), fZ(-2), fN(-1), fCharge(-2), 
  fMomentum(TLorentzVector()), fPosition(TLorentzVector()),
  fEventID(0), fFragmentID(0)
{}

AMDParticle::~AMDParticle()
{}
    
AMDParticle::AMDParticle(Int_t pdg, Int_t z, Int_t n, Int_t chg, 
    TLorentzVector mom, TLorentzVector pos, Int_t eve, Int_t frg)
: fPdg(pdg), fZ(z), fN(n), fCharge(chg), 
  fEventID(eve), fFragmentID(frg)
{
  fMomentum = mom; fPosition = pos;
}

