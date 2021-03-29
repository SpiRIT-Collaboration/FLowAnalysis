#include "myAnalysis.h"
#include "ImageProcessing.h"
#include "makeup_functions.h"
#include "DoFlow.h"

using namespace ImageProcessing;

std::vector<TString> embedNames = {"embed132","embed108","embed112","embed124"};
std::vector<TString> partNames  = {"proton","deuteron","triton","he3"};
std::vector<int> embedSys  = {1};
std::vector<int> embedPart = {0,1,2,3};
//std::vector<int> embedPart = {1,2};
//const int nPart=4;

UInt_t  phicutID = 4;
UInt_t  cutndf   = 50;
Float_t cutdist  = 20.;
TCut    phiCut[] = { "abs(mcphi)<=45||abs(mcphi)>=135", //0
		     "abs(mcphi)<=45",                  //1
		     "abs(mcphi)>=135",                 //2
		     "abs(mcphi)>45 && abs(mcphi)<135", //3
		     ""};                               //4
TString phiname[] = { "45or135", "45", "135", "45to135","all"};

TString fileFooder = "_"+phiname[phicutID] + Form("_ndf%d_dis%d",cutndf, (Int_t)cutdist);


const int multThresholds[3][2]={{56,55},{51,50},{39,39}};

TString beamLR = "left"; 

const double invProjA[] = {44.1/1000., 55.2/1000.,45./1000. };
const double yAA[]   = {0.3822, 0.3647, 0.3538, 0.3902};
const double yNN[]   = {0.3696, 0.3697, 0.3705, 0.3706};
const double yBeam[] = {0.7421, 0.7423, 0.7439, 0.7441};
const double uNN[]   = {0.3709, 0.3709, 0.3718, 0.3719};

TString yTitle  = "#it{y}_{0} = #it{y}/#it{y}^{c.o.m.}_{NN}-1";
TString ptTitle = "#it{p_{T}} (GeV/#it{c})";
TString utTitle = "#it{u_{t0}}";
TString yxTitle = "#it{y}_{x}/#it{y}^{c.o.m.}_{NN}";

TString obj_name_eff(int i, int j, int k, int iter)
{
	return (TString)Form(obj_name(i,j,k)+"_iter%d",iter);
}


TString additionalExt   = ".";
//TString additionalExt   = ".EratZ1.";
//TString additionalExt   = ".b0.2.";

TString embedDate = "20200907";
TString embedExt  = ".left.vaprip";


TString embedFileName(int i, int j)
{ 
	TString path = "/home/common/isobe/embedding/"+embedDate+"/treefiles.embed132.1M/";
	return path+"tree_"+partNames[j]+"_"+embedNames[i]+".root";
}

TString reducedFileName = "rootfiles/embedding/reducedEmbed."+embedDate+embedExt+".root"; 

TString filePath = "rootfiles/LCPSpectra/";
TString dataFileName(int i=0,TString momName="va",TString LR="left", Bool_t isELT=kFALSE)
{
	TString dataMomName = momName=="vapri"?"va":momName;
	TString fileName    = "LCPSpectra"+obj_name(i)+"."+dataMomName+"."+LR+".root";
	if(isELT) fileName  = "LCPSpectra"+obj_name(i)+"."+dataMomName+"."+LR+".corrELT.root";
	return filePath+fileName;
}
//
TString filePath2 = "/home/mizuki/SpiRITAnalysis/macros/data/";
TString inputFileName(UInt_t isys=1)
{
  TString dataName = "Acceptance_"+rsys[isys]+"Sn" + fileFooder + ".root";
  return filePath2+dataName;
}
TString histName(UInt_t ipid=0, TString mlabel="_55to80")
{
  TString dataName = "eff_"+lpid[ipid] + mlabel;
  return dataName;
}

/* TString inputFileName(UInt_t isys=0, UInt_t ipid=0, TString version="v52.11.0") */
/* { */
/*   TString dataName = "finYPt_"+rsys[isys]+"Sn_"+fpid[ipid]+"."+ version +".root"; */
/*   return filePath2+dataName; */
/* } */

TString dataSysFileName(int i=0,TString momName="va",TString LR="left")
{
	TString dataMomName  = momName=="vapri"?"va":momName;
	TString fileName    = "LCPSpectra"+obj_name(i)+"."+dataMomName+"."+LR+".systematics.root";
	return filePath+fileName;
}
TString dataRunSysFileName(int i=0,TString momName="va",TString LR="left")
{
	TString dataMomName  = momName=="vapri"?"va":momName;
	TString fileName    = "LCPSpectra"+obj_name(i)+"."+dataMomName+"."+LR+".halfrun.root";
	return filePath+fileName;
}



TString outExtName = embedDate+embedExt;
TString corrFileName(TString momName="va",TString bgAssign="signal")
{ 
	TString fileName   = "UnfoldedLCPSpectra."+momName+"."+outExtName+".root";
	//if(isELT) fileName = "UnfoldedLCPSpectra."+momName+"."+outExtName+".corrELT.root";
	if(bgAssign=="signalV") fileName = "UnfoldedLCPSpectra."+momName+"."+outExtName+".Voigt.root";
	return filePath+fileName; 
};
TString amdCorrFileName(TString momName="va")
{ return filePath+"UnfoldedAMDLCPSpectra."+momName+"."+outExtName+".root"; };
TString artCorrFileName(TString momName="va")
{ return filePath+"UnfoldedArtLCPSpectra."+momName+"."+outExtName+".root"; };
TString corrSysFileName(TString momName="va")
{ return filePath+"UnfoldedLCPSpectra."+momName+"."+outExtName+".systematics.root"; };
TString corrRunSysFileName(TString momName="va")
{ return filePath+"UnfoldedLCPSpectra."+momName+"."+outExtName+".halfrun.root"; };
TString corrMassGateSysFileName(TString momName="va")
{ return filePath+"UnfoldedLCPSpectra."+momName+"."+outExtName+".massgate.root"; };

TString corrFileName_slim(TString comm="")
{ return filePath2+"rootfiles/UnfoldedLCPSpectra"+fileFooder + "_"+ comm+".slim.root"; };


// embed tree branch
TLorentzVector mcmom;
TLorentzVector *mcmomptr=nullptr;
double mcy, mcpt, mcyx, mcut, mcphi;
int    ntrk, nclus, ngrtrk, ndf;
double b0;
// reco tree branch
TLorentzVector recomom;
TLorentzVector *recomomptr=nullptr;
double recoy, recopt, recoyx, recout, recophi, recodist;
void SetBranchAddressMCTree(TTree* mctree)
{
	mctree->SetBranchAddress("mcmom",&mcmomptr);
	mctree->SetBranchAddress("mcy",&mcy);
	mctree->SetBranchAddress("mcyx",&mcyx);
	mctree->SetBranchAddress("mcpt",&mcpt);
	mctree->SetBranchAddress("mcut",&mcut);
	mctree->SetBranchAddress("mcphi",&mcphi);
	mctree->SetBranchAddress("ngrtrk",&ngrtrk);
	mctree->SetBranchAddress("ntrk",&ntrk);
	mctree->SetBranchAddress("b0",&b0);
}

void SetBranchAddressRecoTree(TTree* recotree, TString momName)
{
	recotree->SetBranchAddress(momName+"mom",&recomomptr);
	recotree->SetBranchAddress(momName+"y",&recoy);
	recotree->SetBranchAddress(momName+"yx",&recoyx);
	recotree->SetBranchAddress(momName+"pt",&recopt);
	recotree->SetBranchAddress(momName+"ut",&recout);
	recotree->SetBranchAddress(momName+"phi",&recophi);
	recotree->SetBranchAddress("recodist",&recodist);
	recotree->SetBranchAddress("ndf",&ndf);
	recotree->SetBranchAddress("nclus",&nclus);
	recotree->SetBranchAddress("ngrtrk",&ngrtrk);
	recotree->SetBranchAddress("ntrk",&ntrk);
	recotree->SetBranchAddress("b0",&b0);
	recotree->SetBranchAddress("mcy",&mcy);
	recotree->SetBranchAddress("mcpt",&mcpt);
	recotree->SetBranchAddress("mcut",&mcut);
}

	
/*
const int    nYBin     = 28;
const double yBinWidth = 0.15;
const double yRange    = 2.1;
const int    nPtBin[]  = {28,32,36,36,40,40};
const double ptBinWidth = 0.05;
const double ptRange[] = {1.4,1.6,1.8,1.8,2.,2.};
*/

double getweight(TH2* hw, double y, double pt, TH2* hembed, double yfac=-1., double ptfac=-1.)
{
  if( hw == NULL ) return 0.;

  //  hw->Print();

  double y_w  = y;
  double pt_w = pt;
  if(yfac!=-1)  y_w  = yfac*y;
  if(ptfac!=-1) pt_w = ptfac*pt;
  double ymin  = hw->GetXaxis()->GetXmin();
  double ymax  = hw->GetXaxis()->GetXmax()-0.01*hw->GetXaxis()->GetBinWidth(hw->GetNbinsX());
  y_w = TMath::Range(ymin,ymax,y_w);
  double ptmin = hw->GetYaxis()->GetXmin();
  double ptmax = hw->GetYaxis()->GetXmax()-0.01*hw->GetXaxis()->GetBinWidth(hw->GetNbinsX());
  pt_w = TMath::Range(ptmin,ptmax,pt_w);
  //if(y>=2)std::cout<<" ranged: y"<<y<<", pt"<<pt<<std::endl;
  double econt = hembed->GetBinContent(hembed->GetXaxis()->FindBin(y),hembed->GetYaxis()->FindBin(pt));
  double w = hw->Interpolate(y_w,pt_w)/econt;
  // use reversal at y<0
  if(y<-0.75)   w = hw->Interpolate(-y_w,pt_w)/econt;
  else if(y<0.) w = 0.5*(hw->Interpolate(y_w,pt_w)+hw->Interpolate(-y_w,pt_w))/econt;

  if( std::isinf( w ) ) return 0.;

  return w;
}


TH1D* IntPt(TH2* h, int binLow=0, int binUp=-1, double yLowerLimit=-1.)
{
	if(h==nullptr) return nullptr;
	TAxis *rapAxis = h->GetXaxis(), *ptAxis = h->GetYaxis();
	if(binUp==-1) binUp = h->GetNbinsY();
	TH1D* h1RapProj = new TH1D("h1RapProj","",h->GetNbinsX(),rapAxis->GetXmin(),rapAxis->GetXmax());
	auto yTgtBin = h->GetXaxis()->FindBin(yLowerLimit);
	for(int rapBin = yTgtBin+1; rapBin<=h->GetNbinsX()+1; ++rapBin){
		double contentBuf=0., errorBuf=0.;
		for(int ptBin=binLow; ptBin<=binUp; ++ptBin){
			contentBuf += h->GetBinContent(rapBin,ptBin)*ptAxis->GetBinWidth(ptBin);
			errorBuf   += TMath::Power(h->GetBinError(rapBin,ptBin)*ptAxis->GetBinWidth(ptBin),2);
		}
		h1RapProj->SetBinContent(rapBin,contentBuf);
		h1RapProj->SetBinError(rapBin,TMath::Sqrt(errorBuf));
	}
	return h1RapProj;
}
