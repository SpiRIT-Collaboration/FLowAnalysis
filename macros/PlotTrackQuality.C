#include "openRunAna.C" 
#include "DoFlow.h"

void PlotTrackQuality()
{

  gROOT->Reset();
  openRunAna();

  if(rChain != NULL)
    LOG(INFO) << " DoFLow_adv: System " << isys << "  -> " << sysName << FairLogger::endl;

  else
    exit(0);

  UInt_t isys = 1; //180Sn
  Int_t   embedPartID[]   = {2212,1000010020,1000010030,1000020030};
  TString embedPartName[] = {"H", "2H", "3H", "3He"};
  UInt_t mult[][2] = { {42, 52}, {55, 80}, {40, 55}, {30, 40}, {0, 30}};
  const Int_t mmax = sizeof(mult)/sizeof(UInt_t)/2;
  y_norm = y_cm[isys+6];

  //--------------------------------------------------                                                        
  //--- Booking
  //--------------------------------------------------                                                                                                    
  TFile *outFile = TFile::Open("data/Acceptance_108Sn.root","recreate");

  TH2D *hyuteff[4][mmax];
  TH2D *hypteff[4][mmax];
  TH1I *hmtrack2[4][mmax];
  TH1I *hmtrack6[4][mmax];
  
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)     

  
  for(auto i: ROOT::TSeqI(4))for(auto j: ROOT::TSeqI(mmax) ){

      TString hlabel = Form("hyuteff_"+lpid[i]+"_%dto%d",mult[j][0],mult[j][1]);
      hyuteff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50, 0., 2.5);
      hlabel         = Form("hypteff_"+lpid[i]+"_%dto%d",mult[j][0],mult[j][1]);
      hypteff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50, 0., 2.);
      hlabel         = Form("hmtrack2_"+lpid[i]+"_%dto%d",mult[j][0],mult[j][1]);
      hmtrack2[i][j] = new TH1I(hlabel, hlabel, 80, 0, 80);
      hlabel         = Form("hmtrack6_"+lpid[i]+"_%dto%d",mult[j][0],mult[j][1]);
      hmtrack6[i][j] = new TH1I(hlabel, hlabel, 80, 0, 80);

    }


  //--------------------------------------------------                                                                                                    
  //--- Event Loop                                                                                                                                        
  //--------------------------------------------------                                                                                                    
  Long64_t nEntry = SetBranch();
  for(Long64_t i = 0; i < nEntry; i++){

    ShowProcess(i);

    rChain->GetEntry(i);

    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);
    /// centrality selection                                                                                                                              
    Int_t trackselection = aflow->mtrack2;

    Int_t hmult = -1;
    for( auto j: ROOT::TSeqI(mmax) ) {
      if( trackselection >= mult[j][0] && trackselection < mult[j][1] ) {
	hmult = j;
	break;
      }
    }
    
    if( hmult < 0 ) continue;

    TIter next(aArray);
    STParticle *aPart = NULL;

    //--------------------------------------------------                                                                                                  
    //----- Main loop                                                                                                                                     
    //--------------------------------------------------                                                                                                  
    UInt_t mtk = 0;

    while( (aPart = (STParticle*)next()) ) {
    
      if(aPart->GetGoodTrackFlag() != 1111 || aPart->GetReactionPlaneFlag() == 0 ) continue;

      auto pid   = aPart->GetPID();
      auto pt    = aPart->GetRotatedMomentum().Pt()/1000.;
      auto rapid  = aPart->GetRapiditycm();;
      auto rapidn = rapid / y_norm;
      auto fmass  = aPart->GetMass();
      Double_t u_t0  = aPart->GetRotatedMomentum().Pt()/fmass/u_p;

      Int_t hpid = -1;
      for( auto i: ROOT::TSeqI(4) ){
	if( pid == embedPartID[i] ) {
	  hpid = i;
	  break;
	}
      }

      if( hpid > -1 ) {
	hyuteff[hpid][hmult] -> Fill( rapidn, u_t0);
	hypteff[hpid][hmult] -> Fill( rapidn, pt);
	hmtrack2[hpid][hmult]-> Fill(aflow->mtrack2);
	hmtrack6[hpid][hmult]-> Fill(aflow->mtrack6);
      }
    }
  }


  outFile->Write();

  TCanvas *ccv;
  Int_t id;

  ccv = new TCanvas("ccv2","ccv2",1000,1000);
  ccv->Divide(4,mmax);
  id = 1;

  for( auto i: ROOT::TSeqI(4))for( auto j: ROOT::TSeqI(mmax) ){
      ccv->cd(id); id++;
      hmtrack2[i][j]->Draw();
    }

  ccv = new TCanvas("ccv3","ccv3",1000,1000);
  ccv->Divide(4,mmax);
   id = 1;

  for( auto i: ROOT::TSeqI(4))for( auto j: ROOT::TSeqI(mmax) ){
      ccv->cd(id); id++;
      hmtrack6[i][j]->Draw();
    }

  ccv = new TCanvas("ccv1","ccv1",1000,1000);
  ccv->Divide(4,mmax);
  id = 1;

  for( auto i: ROOT::TSeqI(4))for( auto j: ROOT::TSeqI(mmax) ){
      ccv->cd(id); id++;
      hyuteff[i][j]->Draw("colz");
    }

  ccv = new TCanvas("ccv4","ccv4",1000,1000);
  ccv->Divide(4,mmax);
  id = 1;

  for( auto i: ROOT::TSeqI(4))for( auto j: ROOT::TSeqI(mmax) ){
      ccv->cd(id); id++;
      hypteff[i][j]->Draw("colz");
    }

}

