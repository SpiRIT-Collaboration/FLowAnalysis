#include "STMassCalculator.hh"
#include "TError.h"
#include "TSystem.h"

STMassCalculator::STMassCalculator()
	: fFile(nullptr), fGraphName("g2Calib")
{
	gErrorIgnoreLevel = kBreak, kSysError, kFatal;
	finder.SetMethod(ROOT::Math::RootFinder::kGSL_SECANT);
	bbPar = new Double_t[7];
	bbPar[2]=1.; bbPar[3]=0.;
	funcBB  = [&](double x){ return bbfunc.fBBMassFinderEq(&x,bbPar); };
	dfuncBB = [&](double x){ return bbfunc.fBBMassFinderDeriv(&x,bbPar); };
	gradfunc1dBB = new ROOT::Math::GradFunctor1D(funcBB,dfuncBB);

	isLoadParameter = kFALSE;
	std::cout<<"=== STMassCalculator is constructed."<<std::endl;
	std::cout<<"=== Please keep in mind that gErrorIgnoreLevel is changed."<<std::endl;
}

STMassCalculator::~STMassCalculator()
{
	delete bbPar;
	delete gradfunc1dBB;
	delete g2PIDCalib[0];
	delete g2PIDCalib[1];
	delete fFile;
}

void STMassCalculator::SetParameter(TString fileName)
{
	if(!gSystem->IsFileInIncludePath(fileName)){
		std::cout<<"Parameter file is not found in designated name."<<std::endl;
		return;
	}
	fFile = TFile::Open(fileName);
	if(fFile) std::cout<<fileName<<" is loaded."<<std::endl;
	fFile->GetObject(fGraphName+"_132Sn_Proton_Par0",g2PIDCalib[0]);
	fFile->GetObject(fGraphName+"_132Sn_Proton_Par1",g2PIDCalib[1]);
	if(g2PIDCalib[0]!=nullptr&&g2PIDCalib[1]!=nullptr){
		std::cout<<"=== STMassCalculator: "<<g2PIDCalib[0]->GetName()<<std::endl;
		std::cout<<"===               and "<<g2PIDCalib[1]->GetName()<<" were loaded."<<std::endl;
		isLoadParameter=kTRUE;
	}

//	file->Close();
//	delete file;
}

void STMassCalculator::SetTGraph2D(TGraph2D* g2Par0, TGraph2D* g2Par1)
{
	if(isLoadParameter){
		std::cout<<" === STMassCalculator::SetTGraph2D() ==="<<std::endl;
		std::cout<<"      Parameter map is already loaded."<<std::endl;
		return;
	}
	g2PIDCalib[0] = g2Par0;
	g2PIDCalib[1] = g2Par1;
	if(g2PIDCalib[0]!=nullptr&&g2PIDCalib[1]!=nullptr){
		std::cout<<"=== STMassCalculator: "<<g2PIDCalib[0]->GetName()<<std::endl;
		std::cout<<"===               and "<<g2PIDCalib[1]->GetName()<<" were loaded."<<std::endl;
		isLoadParameter=kTRUE;
	}

}

Double_t STMassCalculator::CalcMass(Double_t z, TVector3 mom, Double_t dEdx)
{
	if(!(z==1.||z==2.)||!isLoadParameter) return -1.;
	Double_t mass  = -1.;
	Double_t pitch = TMath::RadToDeg()*TMath::ATan(mom.Py()/mom.Pz());
	Double_t yaw   = TMath::RadToDeg()*TMath::ATan(mom.Px()/mom.Pz());
	bbPar[0] = g2PIDCalib[0]->Interpolate(yaw,pitch);
	bbPar[1] = g2PIDCalib[1]->Interpolate(yaw,pitch);
	bbPar[4] = z;
	bbPar[5] = mom.Mag();
	bbPar[6] = dEdx;
	finder.SetFunction(*gradfunc1dBB,z==1?1500:3500);
	if(finder.Solve())
		mass = finder.Root();
	return mass;
}

Double_t STMassCalculator::CalibdEdx(TVector3 mom,Double_t dEdx)
{
	if(!isLoadParameter) return -1;
	Double_t pitch = TMath::RadToDeg()*TMath::ATan(mom.Py()/mom.Pz());
	Double_t yaw   = TMath::RadToDeg()*TMath::ATan(mom.Px()/mom.Pz());
	return (dEdx-g2PIDCalib[1]->Interpolate(yaw,pitch))/g2PIDCalib[0]->Interpolate(yaw,pitch);
}
