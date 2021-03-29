#ifndef myanalysis_h
#define myanalysis_h 1

#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TEntryList.h"
#include "TDirectory.h"
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TCut.h"

using namespace ROOT;

TString anaPath = "/home/kaneko/analysis/20200529/";
TString recoExt = "develop.1988.bf2b00e";

// variables
std::vector<Int_t> beamID = {0,1,2,3}; 
const Int_t nBeam     = 4;
const Int_t beamA[]   = {132,108,112,124};
const Int_t targetA[] = {124,112,124,112};
const Int_t beamZ[]   = {50,50,50,50};
const Int_t targetZ[] = {50,50,50,50};
const Int_t beamN[]   = {82,58,62,74};
const Int_t targetN[] = {74,62,74,62};
const TString collName[] = {"^{132}Sn+^{124}Sn","^{108}Sn+^{112}Sn","^{112}Sn+^{124}Sn","^{124}Sn+^{112}Sn"};

const Int_t   nParticle  = 6;
const TString pidName[]  = {"Proton","Deuteron","Triton","3He","4He","6He","6Li","7Li"};
const Double_t amu       = 931.478;
const Double_t pidMass[] = { 938.2720813, 1875.612762, 2808.921112, 2808.39132, 3727.379378, 6.0188*amu, 6.0151*amu, 7.016*amu };
//	constexpr Double_t kme   = 0.5109989461;
//
TString obj_name(int beamId=-1, int partId=-1, int multId=-1)
{
	TString name;
	if(beamId!=-1) name = Form("_%dSn",beamA[beamId]);
	if(partId!=-1) name = name+"_"+pidName[partId];
	if(multId!=-1) name = Form(name+"_mbin%d",multId);
	return name;
}


// load TChains and vertex cut files.
TChain       *c[4];
Int_t        run,event;
TLorentzVector *vecAAframe=nullptr;
TLorentzVector *vecNNframe=nullptr;
TLorentzVector *vecBeam=nullptr;
TClonesArray *recoPartArray = nullptr;
TClonesArray *kyotoHit    = nullptr;
TEntryList   *elist[4];
Bool_t       isLoad = kFALSE;
Int_t        multTPC;
Double_t     projA, projB;
Double_t     raveVx, raveVy, raveVz;
Double_t     bdcVx, bdcVy, bdcVz;
Double_t     z,aoq;
Bool_t       sigma30;
Bool_t       isGGClose;

const Int_t minbiasRun[][2] = {{3015,3036},{2467,2485},{2603,2610},{3106,3136}};

void loadAnaTree(int systemID=-1, bool isLoose=false/*loose event selection*/, bool isAll=false/*use full stat.*/)
{
	if(systemID>=0&&systemID<=3){ beamID.clear(); beamID.push_back(systemID); }

	for(auto i: beamID){
		c[i] = new TChain("anaTree");
		c[i]->Add(Form(anaPath+"rootfiles/anaTree%d/run*.root",beamA[i]));

		c[i]->SetBranchAddress("run",&run);
		c[i]->SetBranchAddress("event",&event); 
		
		c[i]->SetBranchAddress("raveVz",&raveVz);
		c[i]->SetBranchAddress("raveVx",&raveVx);
		c[i]->SetBranchAddress("raveVy",&raveVy);
		c[i]->SetBranchAddress("bdcVz",&bdcVz);
		c[i]->SetBranchAddress("bdcVx",&bdcVx);
		c[i]->SetBranchAddress("bdcVy",&bdcVy);
		c[i]->SetBranchAddress("projA",&projA);
		c[i]->SetBranchAddress("projB",&projB);
		c[i]->SetBranchAddress("sigma30",&sigma30);
		c[i]->SetBranchAddress("z",&z);
		c[i]->SetBranchAddress("aoq",&aoq);
		c[i]->SetBranchAddress("isGGClose",&isGGClose);
		c[i]->SetBranchAddress("vecAAframe",&vecAAframe);
		c[i]->SetBranchAddress("vecNNframe",&vecNNframe);
		c[i]->SetBranchAddress("vecBeam",&vecBeam);
		c[i]->SetBranchAddress("kyotoHit",&kyotoHit);
	  c[i]->SetBranchAddress("part",&recoPartArray);
	  c[i]->SetBranchAddress("multTPC",&multTPC);
	}

	auto vertexCutter = TFile::Open(anaPath+"rootfiles/Vertex.root");
	TF1 *f1Vz[4];
	TF1 *f1VxBx[4], *f1VyBy[4];
	//TCut cutRuns = "!(run==2285||run==2848||run==2984||run==2632||run==2601||run==2548||run==2549||run==2652||run==2653)";
	//TCut cutRuns = "!(run==2440||run==2442||run==2453||run==2461)";
	TCut cutRuns = "";

	for(auto i: beamID){
		vertexCutter->GetObject("f1Vz"+obj_name(i),f1Vz[i]);
		vertexCutter->GetObject("f1VxBx"+obj_name(i),f1VxBx[i]);
		vertexCutter->GetObject("f1VyBy"+obj_name(i),f1VyBy[i]);
//		vertexCutter->GetObject(Form("fitVxy_%dSn",beamA[i]),fitVxy[i]);
		Double_t vzPar[4];  f1Vz[i]->GetParameters(vzPar);
		Double_t vxbxPar[3];  f1VxBx[i]->GetParameters(vxbxPar);
		Double_t vybyPar[3];  f1VyBy[i]->GetParameters(vybyPar);
	//	Double_t vxyPar[6]; fitVxy[i]->GetParameters(vxyPar);
		TCut cutMinBias  = Form("!(run>=%d&&run<=%d)",minbiasRun[i][0],minbiasRun[i][1]);
		TCut cutBeam     = Form("beam==%d&&sigma30",beamA[i]);
		TCut cutGGClose  = "!isGGClose";
		TCut cutVertexZ  = Form("TMath::Abs(raveVz-%lf)<=%lf",vzPar[1],vzPar[2]*3.);
		TCut cutVertexXY = "TMath::Abs(raveVx)<=15.&&TMath::Abs(raveVy+205)<=20.";
		TCut cutVxBxCor  = Form("TMath::Abs(raveVx-bdcVx-%lf)<=%lf",vxbxPar[1],vxbxPar[2]*3.);
		TCut cutVyByCor  = Form("TMath::Abs(raveVy-bdcVy-%lf)<=%lf",vybyPar[1],vybyPar[2]*3.);

		TCut cutTotal    = cutRuns+cutMinBias+cutBeam+cutGGClose+cutVertexZ+cutVertexXY+cutVxBxCor+cutVyByCor;
		if(isLoose){
			cutBeam = Form("beam==%d&&TMath::Abs(z-50)<=3&&TMath::Abs(aoq-%lf)<=0.1",beamA[i],(double)beamA[i]/beamZ[i]);
			cutVertexZ  = Form("TMath::Abs(raveVz-%lf)<=%lf",vzPar[1],vzPar[2]*5.);
			cutTotal    = cutRuns+cutMinBias+cutBeam+cutVertexZ+cutVertexXY;
		}
		if(isAll){
			cutTotal = cutRuns+cutMinBias;
		}
	   
		c[i]->Project(Form("elist_%d",i),"",cutTotal,"entrylist");
		gDirectory->GetObject(Form("elist_%d",i),elist[i]);
		c[i]->SetEntryList(elist[i]);
	}

	isLoad = kTRUE;
}


//const Double_t trigEff[]  = { 0.306, 0.356, 0.358, 0.333 };
const Double_t trigEff[]  = { 0.386, 0.385, 0.388, 0.386 };
const Double_t trigB0[]   = { 0.602517, 0.609558 , 0.616405, 0.588172 };

//const Int_t    nMultBin   = 5; // same as FOPI, and no-centrality cut
const Int_t    nMultBin   = 4; 
const TString  centralTitle[] = {"#it{b}_{0}#leq0.15","0.15<#it{b}_{0}#leq0.25","0.25<#it{b}_{0}#leq0.45","trig.bias"};
const Double_t centralCut[] = {0.0225,0.0625,0.2025};
const Double_t b0Cut[] = {0.15,0.25,0.45};
//const Double_t scaledBCut[] = {0.1,0.2,0.4};
//const Double_t scaledBCut[] = {0.2,0.3,0.5};
const Int_t    nBClass  = 4;
//const TString  bTitle[] = { "b_{0}^{MUL}#leq 0.15","0.15 < b_{0}^{MUL}#leq 0.25","0.25 < b_{0}^{MUL}#leq 0.45",
const TString  bTitle[] = { "b_{0}#leq 0.15","0.15 < b_{0}#leq 0.25","0.25 < b_{0}#leq 0.45","trig.bias"};
//const TString  bTitle[] = { "b_{0}^{mult} #leq 0.15","0.15 < b_{0}^{mult} #leq 0.25","0.25 < b_{0}^{mult} #leq 0.45","trig.bias"};


const double pidMRange[][2] = {{500.,1400}, {1400.,2300.}, {2300.,3400.}, {2300.,3300.}, {3300.,4200.}, {5200.,6000.}};

// mass pid analysis
const Int_t nRBin = 30;
const Int_t nRBinBase = 4;
//const Double_t mRange[][2] = {{500.,1400}, {800.,2900}, {1700.,3700.}, {2000.,3800.}, {2600.,4400.}, {4500.,7000.}};
const Double_t mRange[][2] = {{0.,2100}, {700.,3000}, {1500.,4000.}, {1400.,4100.}, {2600.,4600.}, {4500.,7000.}};
const Double_t rRange[][2] = {{100.,1400.}, {200.,2200.}, {300.,2800.}, {300.,1400.}, {400.,1800.}, {1300.,2500.}};
//const Int_t    nMBin[]     = {60,60,60,40,40,30};
const Int_t    nMBin[]     = {80,80,60,60,40,40};

Double_t Gauss(Double_t *x, Double_t *p)
{ return p[0]*TMath::Gaus(x[0],p[1],p[2]); }
Double_t TripleGauss(Double_t *x, Double_t *p)
{ return p[0]*TMath::Gaus(x[0],p[1],p[2])+p[3]*TMath::Gaus(x[0],p[4],p[5])+p[6]*TMath::Gaus(x[0],p[7],p[8]); }
Double_t SigToTotalGauss(Double_t *x, Double_t *p)
{ return Gauss(x,p)/TripleGauss(x,p); }

Double_t TripleVoigt(Double_t *x, Double_t *p)
{ return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[2])+p[3]*TMath::Voigt(x[0]-p[4],p[5],p[5])+p[6]*TMath::Voigt(x[0]-p[7],p[8],p[8]); }
Double_t SigToTotalVoigt(Double_t *x, Double_t *p)
{ return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[2])/TripleVoigt(x,p); }

void dndx(TH1* h, Double_t entries=-1)
{
	if(entries==-1) entries=h->GetEntries();
	if(h->GetDefaultSumw2()) h->Sumw2();
	for(int ibin=0; ibin<=h->GetNbinsX()+1; ++ibin){
		auto cont = h->GetBinContent(ibin);
		auto wBin = h->GetBinWidth(ibin);
		auto err  = h->GetBinError(ibin);
		auto scaleFactor = 1./wBin/entries;
		h->SetBinContent(ibin,cont*scaleFactor);
		h->SetBinError(ibin,err*scaleFactor);
	}
}

void dndxy(TH2* h, Double_t entries=-1)
{
	if(entries==-1) entries=h->GetEntries();
	if(h->GetDefaultSumw2()) h->Sumw2();
	for(int xbin=0; xbin<=h->GetNbinsX()+1; ++xbin)for(int ybin=0; ybin<=h->GetNbinsY()+1; ++ybin){
		auto gbin = h->GetBin(xbin,ybin);  // get global bin number.
		auto cont = h->GetBinContent(gbin);
		auto err  = h->GetBinError(gbin);
		auto wBinX = h->GetXaxis()->GetBinWidth(xbin);
		auto wBinY = h->GetYaxis()->GetBinWidth(ybin);
		auto scaleFactor = 1./wBinX/wBinY/entries;
		h->SetBinContent(gbin,cont*scaleFactor);
		h->SetBinError(gbin,err*scaleFactor);
	}
}

// pbuu analysis
const int nSet = 4;
const TString histNamepBUU[]={
	"132Sn124Sn_Soft","132Sn124Sn_Hard",
	"108Sn112Sn_Soft","108Sn112Sn_Hard"
};
const TString histTitlepBUU[]={
	"^{132}Sn^{124}Sn #gamma=0.5","^{132}Sn^{124}Sn #gamma=1.75",
	"^{108}Sn^{112}Sn #gamma=0.5","^{108}Sn^{112}Sn #gamma=1.75"
};

#endif
