#ifndef AMDPARTICLE
#define AMDPARTICLE 1

#include "TObject.h"
#include "TLorentzVector.h"

class AMDParticle : public TObject
{
  public:
    AMDParticle();
    virtual ~AMDParticle();
    AMDParticle(Int_t pdg, Int_t z, Int_t n, Int_t chg, 
    	TLorentzVector mom, TLorentzVector pos, Int_t eve, Int_t frg);


    Int_t fPdg;    // pdg code ( pion 211, proton 2212, neutron 2122, deuteron 1000010020, ... 1000000000+Z*10000+    A*10 )
    Int_t fZ;      // # of proton
    Int_t fN;      // # of neutron
    Int_t fCharge; // charge in e unit.
    TLorentzVector fMomentum;  // 4-d momentum at the primary vertex
    TLorentzVector fPosition;  // 4-d position of the primary vertex

    Int_t fEventID; // event # of AMD
    Int_t fFragmentID;  // fragment ID of AMD

    ClassDef(AMDParticle,1);

};

#endif
