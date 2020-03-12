#include "FlowFunctions.h"
#include "../flowformat/STFlowInfo.hh"
#include "openRunAna.C"

Int_t  seltrackID = 4;
UInt_t selReactionPlanef = 1000;
Int_t  seltrack;

TCanvas *cc[12];
UInt_t im = 0;

Double_t xconst[4];
Double_t xmean[4];
Double_t xsig[4];
Double_t yconst[4];
Double_t ymean[4];
Double_t ysig[4];


const UInt_t  nbin = 25;
  
UInt_t ic = -1;
Int_t ntrack[7];

Int_t iVer;

TString flabel;
TString unitpX;
TString unitpY;
TString ltrack;
TString fhead;
TString fOutName ;


Double_t constX;
Double_t meanX ;
Double_t sigX  ;
Double_t constY;
Double_t meanY ;
Double_t sigY  ;

void Initialize();
void SetPsiCorrectionFileHeader(UInt_t isel);
void SaveReCenteringData(UInt_t m);
void Flatten_Psi_ntrack(UInt_t isel = 0);//%% Executable : 
void ReCentering(UInt_t isel = 10, Int_t nmin=0, Int_t nmax=100);

void DoFlattening(UInt_t isel = 10)
{
  gStyle->SetOptStat(0);
  gROOT->Reset();

  UInt_t ichain = 0;

  openRunAna();

  rChain = (TChain*)gROOT->FindObject("rChain");
  if(rChain == NULL) {    
    cout << "rChain is not found " << endl;
    exit(0);
  }
  
  SetPsiCorrectionFileHeader(isel);
  Flatten_Psi_ntrack(isel);

  // later v17
  //  Flatten_Psi_ntrack(10);  //TVector3(unitP->X(),   unitP->Y(),   0.);
  //  Flatten_Psi_ntrack(11);  //TVector3(unitP_1->X(),   unitP_1->Y(),   0.);

}

//________________________________//%% Executable : 
void ReCentering(UInt_t isel = 10, Int_t nmin=0, Int_t nmax=100) //%% Executable : Recentering calibration
{
  constX = 0.;
  meanX  = 0.;
  sigX   = 0.;
  constY = 0.;
  meanY  = 0.;
  sigY   = 0.;

  auto fgX  = new TF1("fgX","gaus",-30,30);;
  auto fgY  = new TF1("fgY","gaus",-30,30);;

  auto hQx  = new TH1D(Form("hQx%d_%d",nmin,isel),"; Qx",100,-20,20);
  auto hQy  = new TH1D(Form("hQy%d_%d",nmin,isel),"; Qy",100,-20,20);

  std::cout << " hQx " << Form("hQx%d_%d",nmin,isel) << std::endl;

  TString htitle = Form("mult >= %d && mult < %d",nmin, nmax);
  hQx->SetName(Form("hQx%d_%d",nmin,isel) );
  hQy->SetName(Form("hQy%d_%d",nmin,isel) );
  hQx->SetTitle(htitle);
  hQy->SetTitle(htitle);

  TString cutdef = Form(ltrack+">=%d&&"+ltrack+"<%d",nmin,nmax); 
  if( isys != 5 ) cutdef += "&&beamPID>0";

  TCut   multcut = TCut(cutdef);
  multcut.Print();

  rChain->Project(Form("hQx%d_%d",nmin,isel), unitpX, multcut);
  rChain->Project(Form("hQy%d_%d",nmin,isel), unitpY, multcut);

  if(hQx->GetEntries() > 0) {

    hQx->Fit("fgX","Q0");
    hQy->Fit("fgY","Q0");
    
    constX= fgX->GetParameter(0);
    meanX = fgX->GetParameter(1);
    sigX  = fgX->GetParameter(2);
    
    constY= fgY->GetParameter(0);
    meanY = fgY->GetParameter(1);
    sigY  = fgY->GetParameter(2);

    LOG(INFO) << " meanX " << meanX << " meanY " << meanY << FairLogger::endl;
  }

  //  delete hQx;
  //  delete hQy;
}


//________________________________//%% Executable : 
void Flatten_Psi_ntrack(UInt_t isel)
{
  std::cout << "flatten_Psi_ntrackbin " << isel << std::endl;

  Initialize();
  
  const UInt_t harm      = 5;
  const UInt_t ntrknbin  = 20;

  // bin setting for multiplicity
  Double_t ntrkbin[ntrknbin+1];
  Double_t ntrk_min = 0;
  //                      0  1  2  3  4  5  6  7  8  9   10  11   12  13
  Double_t ntrk_max[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 100, 50, 100, 50};
  for(UInt_t n = 0; n < ntrknbin+2; n++)
    ntrkbin[n]   = ntrk_max[isel]/ntrknbin * n;

  auto hbaiphi = new TH2D(Form("hbaiphi_%d",isel),  " #Phi before and after; before #Phi [rad]; after #Phi [rad] ", 
			  400,-3.5,3.5,400,-3.5,3.5);
  auto hniphi  = new TH1D(Form("hniphi_%d",isel),   " #Phi no corr.; Azimuthal angle [rad]"  , 200,-3.2,3.2);
  auto hbiphi  = new TH1D(Form("hbiphi_%d",isel),   " #Phi before  ; Azimuthal angle [rad]"  , 200,-3.2,3.2);
  auto haiphi  = new TH1D(Form("haiphi_%d",isel),   " #Phi after   ; Azimuthal angle [rad]"  , 200,-3.2,3.2);

  auto hbntrkiphi = new TH2D(Form("hbntrkiphi_%d",isel)," before ; Number of tracks; #phi"     , 
			     ntrk_max[isel],0,ntrk_max[isel],400,-3.2,3.2);
  auto hantrkiphi = new TH2D(Form("hantrkiphi_%d",isel)," after  ; Number of tracks; #phi"     , 
			     ntrk_max[isel],0,ntrk_max[isel],400,-3.2,3.2);
    
  auto habiphi = new TH2D(Form("habiphi_%d",isel)," ;#Psi_rot; #Psi_fc",200,-3.2,3.2,200,-3.2,3.2);

  //--- retrivewed data
  Double_t bsPhi_1[3];
  Double_t bsPhi_2[3];
  
  TClonesArray *aFlowArray = NULL;
  TClonesArray *aBDCArray  = NULL;

  if( isel >= 10) {
    rChain->SetBranchAddress("STFlow"  ,&aFlowArray);
  }


  // Flattening with a shifting method
  STFlowCorrection *flowcorr[ntrknbin+1];

  for(UInt_t j = 0; j < ntrknbin+1; j++){   

    flowcorr[j] = new STFlowCorrection(rChain, harm, 0); 
    
    if(j < ntrknbin+1){

      ReCentering(isel, ntrkbin[j], ntrkbin[j+1]);

      flowcorr[j]->SetBin_max(0, ntrkbin[j+1]);
    }

    flowcorr[j]->SetBin_min(0, ntrkbin[j]);
	
    Double_t rcX[3];
    rcX[0] = constX;
    rcX[1] = meanX;
    rcX[2] = sigX;
    Double_t rcY[3];
    rcY[0] = constY;
    rcY[1] = meanY;
    rcY[2] = sigY;

    LOG(INFO) << " meanX " << meanX << " meanY " << meanY << FairLogger::endl;

    flowcorr[j]->SetReCenteringParameter("X",rcX);  
    flowcorr[j]->SetReCenteringParameter("Y",rcY);        
  }   


  // process
    
  Int_t nevt = rChain->GetEntries();
  cout << " Number of events " << nevt << endl;

  Int_t icout = 0;
  std::vector<Double_t> ophi;

  nevt=100000;
  for(UInt_t i = 0; i < nevt; i++){
    rChain->GetEntry(i);

    STFlowInfo *aFlow = (STFlowInfo*)aFlowArray->At(0);
    if( aFlow == NULL ) continue;

    if( aFlow->goodEventf == 0 || aFlow->beamPID == 0 ) continue;
    if( aFlow->mtrack4 < 5 ) continue;
      
    Int_t seltrack;
    TVector3 vec;
    
    switch(isel){
    case 10:
      seltrack = aFlow->mtrack4;
      vec      = aFlow->unitP;
      break;
    case 11:
      seltrack = aFlow->mtrack_1;
      vec      = aFlow->unitP_1;
      break;
    case 12:
      seltrack = aFlow->mtrack4;
      vec      = aFlow->unit2P;
      break;
    case 13:
      seltrack = aFlow->mtrack_1;
      vec      = aFlow->unit2P_1;
    }

    UInt_t intrk  = 0;
    UInt_t j = ntrknbin;
    while(1){ 
      if( seltrack >= ntrkbin[j] ){
	intrk = j;
	break;
      }
      j--;
    }
    
    hniphi->Fill(vec.Phi());

   if(intrk==1)
      ophi.push_back(vec.Phi());
    
    if(intrk <= ntrknbin )   
      flowcorr[intrk]->Add(seltrack, vec);
    
    //-----------------------------
  }

  //----------  get corrrection parameters
  for(UInt_t j = 0; j < ntrknbin+1; j++){   

    //    if(flowcorr[j]->GetNPhi() == 0 ) continue;
    UInt_t nphi = flowcorr[j]->ReCenteringFourierCorrection();
    std::cout << " At " << ntrkbin[j]  << " nphi " << nphi <<  std::endl;
      
    vector<Int_t>    mtk   = flowcorr[j]->GetMTrack();
    vector<Double_t> aphi  = flowcorr[j]->GetCorrectedPhi();
    vector<Double_t> bphi  = flowcorr[j]->GetOriginalPhi();
    vector<Double_t> rcphi = flowcorr[j]->GetReCeneringPhi();
    
    cout << "after " << mtk.size()  << endl;
    if(mtk.size() > 0){
      for(UInt_t k = 0; k < (UInt_t)mtk.size(); k++){	  
	hbiphi     ->Fill(rcphi.at(k));
	hbntrkiphi ->Fill(mtk.at(k)   , rcphi.at(k));	  
	
	hantrkiphi ->Fill(mtk.at(k)   , aphi.at(k));
	haiphi     ->Fill(aphi.at(k));	  

	hbaiphi->Fill(abs(bphi.at(k)), abs(rcphi.at(k)));
	  
	if(j == 1)
	  habiphi  ->Fill(ophi.at(k), aphi.at(k));

      }
    }
      
    TString comm1 = sysName + ".v" + dVer + "." + flabel;
    comm1 += Form(".m%d:flatten_Psi_ntrk; ntrack>= %f && ntrack< %f ",
		  j,ntrkbin[j],ntrkbin[j+1]);
 
    TString comm2 = unitpX + " && " + unitpY;
    cout << "saving in ...  " << comm1 << endl;
    flowcorr[j]-> SaveCorrectionFactor(comm1, comm2);    
  }

       
  im++;
  cc[im] = new TCanvas(Form("cc%d",im),Form("cc%d",im),700,500);
  
  hniphi->SetLineColor(2);
  hniphi->Draw("e");
  
  hbiphi->SetLineColor(8);
  hbiphi->Draw("samee");

  haiphi->SetLineColor(4);
  haiphi->Draw("samee");

  auto aLeg = new TLegend(0.75,0.13,0.9,0.3,"");
  aLeg->AddEntry(hniphi,"No Collection ","lp");
  aLeg->AddEntry(hbiphi,"ReCentering","lp");
  aLeg->AddEntry(haiphi,"ReCentering & Shifting","lp");
  
  aLeg->Draw();

  im++;
  cc[im] = new TCanvas(Form("cc%d",im),Form("cc%d",im),700,500);
  
  hbaiphi->Draw("colz");
}

void SaveReCenteringData(UInt_t m)
{
  gSystem->cd("db");

  std::fstream fout;
  fOutName = "ReCent" +  sysName + ".data";
  fout.open(fOutName, std::fstream::out);

  fout << "ReCentering for " << unitpX << " and " << unitpY << endl;
  fout << " Axis  Constatnt       Mean      Signma " << sysName  
       << endl;
  cout << " X: " 
       << setw(10)  << constX 
       << setw(15)  << meanX
       << setw(10)  << sigX
       << endl;

  cout << " Y: " 
       << setw(10)  << constY
       << setw(15)  << meanY
       << setw(10)  << sigY
       << endl;

  fout << " X: " 
       << setw(10)  << constX 
       << setw(15)  << meanX
       << setw(10)  << sigX
       << endl;


  fout << " Y: " 
       << setw(10)  << constY
       << setw(15)  << meanY
       << setw(10)  << sigY
       << endl;

  fout.close();   
  gSystem->cd("..");

}

void SetPsiCorrectionFileHeader(UInt_t isel){
  switch(isel){
  case 10: // total Psi
    unitpX = "unitP.X()";
    unitpY = "unitP.Y()";
    ltrack = "mtrack4";
    flabel = "psi";
    break;
  case 11:
    unitpX = "unitP_1.X()";
    unitpY = "unitP_1.Y()";
    ltrack = "mtrack_1";
    flabel = "subpsi1";
    break;
  case 12:
    unitpX = "unit2P.X()";
    unitpY = "unit2P.Y()";
    ltrack = "mtrack4";
    flabel = "2psi";
    break;
  case 13:
    unitpX = "unit2P_1.X()";
    unitpY = "unit2P_1.Y()";
    ltrack = "mtrack_1";
    flabel = "sub2psi";
    break;
  }
}

void Initialize()
{
  constX= 0.;
  meanX = 0.;
  sigX  = 1.;
  constY= 0.;
  meanY = 0.;
  sigY  = 1.;

}
