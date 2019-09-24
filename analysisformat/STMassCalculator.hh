#ifndef STMassCalculator_h
#define STMassCalculator_h 1

#include "STBBFunction.hh"
#include "TFile.h"
#include "TGraph2D.h"
#include "TVector3.h"
#include "Math/RootFinder.h"
#include "Math/Functor.h"
#include "functional"

class STMassCalculator
{
	public:
		STMassCalculator();
		~STMassCalculator();

		void SetParameter(TString fileName);
		void SetTGraph2D(TGraph2D* g2Par0, TGraph2D* g2Par1);
		Double_t CalcMass(Double_t z,TVector3 mom,Double_t dEdx);
		Double_t CalibdEdx(TVector3 mom,Double_t dEdx);

	private:
		TFile    *fFile;
		TString  fGraphName;
		TGraph2D *g2PIDCalib[2];
		Bool_t   isLoadParameter;
		ROOT::Math::RootFinder finder;
		Double_t *bbPar;
		STBBFunction bbfunc;
		std::function<double(double)> funcBB;
		std::function<double(double)> dfuncBB;
		ROOT::Math::GradFunctor1D *gradfunc1dBB;

};
#endif
