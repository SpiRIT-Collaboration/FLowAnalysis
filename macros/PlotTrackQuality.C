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
  Int_t   embedPartID[]   = {2212,1000010020,1000010030,1000020030,1000020040};
  TString embedPartName[] = {"H", "2H", "3H", "3He", "4He"};
  const Int_t pidmax = sizeof(embedPartID)/sizeof(Int_t);

  //mtrack2 cut condition
  UInt_t mult[][2] = { {42, 52}, {55, 80}, {40, 55}, {30, 40}, {0, 30}};
  const Int_t mmax = sizeof(mult)/sizeof(UInt_t)/2;
  Bool_t bmult[mmax];

  UInt_t  phicutID = 4;
  UInt_t  cutndf   = 50;
  Float_t cutdist  = 20.;

  TString phiname[] = { "45or135", "45", "135", "45to135","all" };

  y_norm = y_cm[isys+6];
  //--------------------------------------------------                                                        
  //--- Booking
  //--------------------------------------------------     
  //  TString trklabel = "data/Acceptance_108Sn_135.root";               
  //  TString trklabel = "data/Acceptance_108Sn_45.root";                                                        
  //  TString trklabel = "data/Acceptance_108Sn_tight.root";
  //  TString trklabel = "data/Acceptance_108Sn_ndf50.root"; 
  //  TString trklabel = "data/Acceptance_108Sn_45or135_ndf50_dis20.root";           
  //  TString trklabel = "data/Acceptance_108Sn_45or135_ndf20_dis20.root";           
  //  TString trklabel = "data/Acceptance_108Sn_45_ndf50_dis20.root";           
  TString trklabel = "data/Acceptance_108Sn_" + phiname[phicutID];
  trklabel += Form("_ndf%d_dis%d.root",cutndf, (Int_t)cutdist);           

  TFile *outFile = TFile::Open(trklabel,"recreate");
  LOG(INFO) << " output file " << trklabel << FairLogger::endl;

  TH2D *hyuteff[pidmax][mmax];
  TH2D *hypteff[pidmax][mmax];
  TH2D *hypheff[pidmax][mmax];
  TH1I *hmtrack2[pidmax][mmax];
  TH1I *hmtrack6[pidmax][mmax];
  
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)     
  
  for(auto i: ROOT::TSeqI(pidmax))for(auto j: ROOT::TSeqI(mmax) ){

      TString hlabel = Form("hyuteff_"+lpid[i]+"_%dto%d",mult[j][0],mult[j][1]);
      hyuteff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50, 0., 2.5);
      hlabel         = Form("hypteff_"+lpid[i]+"_%dto%d",mult[j][0],mult[j][1]);
      hypteff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50, 0., 2.5);
      hlabel         = Form("hypheff_"+lpid[i]+"_%dto%d",mult[j][0],mult[j][1]);
      hypheff[i][j]  = new TH2D(hlabel, hlabel , 40,-2.,2.,60,-3.15,3.15);

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

    for( auto j: ROOT::TSeqI(mmax) ) {
      if( trackselection >= mult[j][0] && trackselection < mult[j][1] ) {
	bmult[j] = kTRUE;
      }
      else
	bmult[j] = kFALSE;
    }


    TIter next(aArray);
    STParticle *aPart = NULL; 


    //--------------------------------------------------                                                                                                  
    //----- Main loop                                                                                                                                     
    //--------------------------------------------------                                                                                                  
    UInt_t mtk = 0;

    while( (aPart = (STParticle*)next()) ) {
    
      if( aPart->GetGoodTrackFlag() != 1111 ) continue;
      if( aPart->GetReactionPlaneFlag() == 0 ) continue;
      if( aPart->GetNDF() < cutndf ) continue;
      if( aPart->GetDistanceAtVertex() > cutdist ) continue;

      // phi cut
      auto phi    = aPart->GetRotatedMomentum().Phi();
      // w phicut                                                                                                                                              
      switch(phicutID) {
      case 0:
	if( abs(phi) > 45*TMath::DegToRad() && abs(phi) < 135*TMath::DegToRad()) continue;
	break;
      case 1:
	if( abs(phi) > 45*TMath::DegToRad()) continue;
	break;
      case 2:
	if( abs(phi) < 135*TMath::DegToRad() ) continue;
	break;
      case 3:
	if( abs(phi) < 45*TMath::DegToRad() || abs(phi) > 135*TMath::DegToRad()) continue;
	break;
      }


      auto pid   = aPart->GetPID();
      //      auto pid   = aPart->GetPIDTight();
      auto pt    = aPart->GetRotatedMomentum().Pt()/1000.;
      auto rapid  = aPart->GetRapiditycm();;
      auto rapidn = rapid / y_norm;
      auto fmass  = aPart->GetMass();
      Double_t u_t0  = aPart->GetRotatedMomentum().Pt()/fmass/u_p;

      Int_t hpid = -1;
      for( auto i: ROOT::TSeqI(pidmax) ){
	if( pid == embedPartID[i] ) {
	  hpid = i;
	  break;
	}
      }

      if( hpid > -1 ) {
	
	for( auto i : ROOT::TSeqI(mmax) ){
	  if( bmult[i] ) {
	    hyuteff[hpid][i] -> Fill( rapidn, u_t0);
	    hypteff[hpid][i] -> Fill( rapidn, pt);
	    hypheff[hpid][i] -> Fill( rapidn, phi);
	    hmtrack2[hpid][i]-> Fill(aflow->mtrack2);
	    hmtrack6[hpid][i]-> Fill(aflow->mtrack6);
	  }
	}
      }
    }
  }


  outFile->Write();

  TCanvas *ccv;
  Int_t id;
  
  ccv = new TCanvas("ccv2","ccv2",1000,1000);
  ccv->Divide(pidmax,mmax);
  id = 1;

  for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
      ccv->cd(id); id++;
      hmtrack2[i][j]->Draw();
    }

  ccv = new TCanvas("ccv3","ccv3",1000,1000);
  ccv->Divide(4,mmax);
  id = 1;

  for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
      ccv->cd(id); id++;
      hmtrack6[i][j]->Draw();
    }

   ccv = new TCanvas("ccv1","ccv1",1000,1000);
   ccv->Divide(4,mmax);
   id = 1;

   for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
       ccv->cd(id); id++;
       hyuteff[i][j]->Draw("colz");
     }

   ccv = new TCanvas("ccv4","ccv4",1000,1000);
   ccv->Divide(4,mmax);
   id = 1;

   for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
       ccv->cd(id); id++;
       hypheff[i][j]->Draw("colz");
     }
}

