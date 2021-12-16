#include "openRunAna.C" 
#include "DoFlow.h"
#include "../tasks/STLorentzBoostVector.hh"

auto calc_yx = [](TLorentzVector v)
{
  double beta_x = v.Beta()*TMath::Sin(v.Theta())*TMath::Cos(v.Phi());
  return 0.5*TMath::Log((1.+beta_x)/(1.-beta_x));
};


void PlotTrackQuality()
{

  gROOT->Reset();
  openRunAna();

  TDatime dtime;
  TRandom3 grnd(dtime.GetSecond());
  gRandom->SetSeed(dtime.GetSecond());


  if(rChain != NULL)
    LOG(INFO) << " System " << isys << "  -> " << sysName << FairLogger::endl;

  else
    exit(0);

  //  UInt_t isys = 1; //180Sn
  Int_t   embedPartID[]   = {2212,1000010020,1000010030,1000020030,1000020040};
  TString embedPartName[] = {"H", "2H", "3H", "3He", "4He"};
  const Int_t pidmax = sizeof(embedPartID)/sizeof(Int_t);

  //mtrack2 cut condition
  const Int_t mmax = sizeof(trkcut)/sizeof(UInt_t)/2;
  Bool_t btrkcut[mmax];

  //UInt_t  phicutID = 1;
  UInt_t  phicutID = 5;
  UInt_t  cutndf   = 15;
  Float_t cutdist  = 20.;

  TString phiname[] = { "45or135", "45", "135", "45to135","all","-20to30"};

  y_norm = y_cm[isys+6];
  //--------------------------------------------------                                                        
  //--- Booking
  //--------------------------------------------------     
  TString trklabel = "data/Acceptance_"+rsys[isys]+"Sn_" + phiname[phicutID] + ".root";

  TFile *outFile = TFile::Open(trklabel,"recreate");
  LOG(INFO) << " output file " << trklabel << FairLogger::endl;

  TH2D *hyuteff[pidmax][mmax];
  TH2D *hypteff[pidmax][mmax];
  TH2D *hyyxeff[pidmax][mmax];
  TH2D *hyxheff[pidmax][mmax];
  TH2D *hyyteff[pidmax][mmax];
  TH1I *hmtrack1[pidmax][mmax];
  TH1I *hmtrack2[pidmax][mmax];
  TH1I *hmtrack4[pidmax][mmax];

  
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)     
  
  for(auto i: ROOT::TSeqI(pidmax))for(auto j: ROOT::TSeqI(mmax) ){

      TString hlabel = Form("hyuteff_"+lpid[i]+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hyuteff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50, 0., 2.5);
      hlabel         = Form("hypteff_"+lpid[i]+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hypteff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50, 0., 2.5);
      hlabel         = Form("hyyxeff_"+lpid[i]+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hyyxeff[i][j]  = new TH2D(hlabel, hlabel , 40,-2.,2.,40,-2.,2.);
      hlabel         = Form("hyxheff_"+lpid[i]+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hyxheff[i][j]  = new TH2D(hlabel, hlabel , 40,-2.,2.,20,0.,2.);
      hlabel         = Form("hyyteff_"+lpid[i]+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hyyteff[i][j]  = new TH2D(hlabel, hlabel , 40,-2.,2.,20,0.,2.);

      hlabel         = Form("hmtrack1_"+lpid[i]+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hmtrack1[i][j] = new TH1I(hlabel, hlabel, 80, 0, 80);
      hlabel         = Form("hmtrack2_"+lpid[i]+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hmtrack2[i][j] = new TH1I(hlabel, hlabel, 80, 0, 80);
      hlabel         = Form("hmtrack4_"+lpid[i]+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hmtrack4[i][j] = new TH1I(hlabel, hlabel, 80, 0, 80);

    }


  //--------------------------------------------------                                                                                                    
  //--- Event Loop                                                                                                                                        
  //--------------------------------------------------                                                                                                    


  Long64_t nEntry = SetBranch();
  for(Long64_t i = 0; i < nEntry; i++){

    ShowProcess(i);

    Double_t rnd_Phi = 2.*TMath::Pi()*(grnd.Rndm() - 0.5);

    rChain->GetEntry(i);

    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);
    /// centrality selection                                                                                                                              
    Int_t trackselection = aflow->mtrack1;

    UInt_t trkcut_sel = 0;
    for( auto j: ROOT::TSeqI(mmax) ) {
      if( trackselection >= trkcut[j][0] && trackselection < trkcut[j][1] ) {
	btrkcut[j] = kTRUE;
	trkcut_sel++;
      }
      else
	btrkcut[j] = kFALSE;

    }

    if( trkcut_sel == 0 ) continue;

    TIter next(aArray);
    STParticle *aPart = NULL; 

    //--------------------------------------------------                                                                                                  
    //----- Main loop                                                                                                                                     
    //--------------------------------------------------                                                                                                  
    UInt_t mtk = 0;

    while( (aPart = (STParticle*)next()) ) {
    
      if( aPart->GetGoodTrackFlag() != 1111 ) continue;
      //      if( aPart->GetReactionPlaneFlag() == 0 ) continue;

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
      case 5:
	if( phi > 30*TMath::DegToRad() || phi < -20*TMath::DegToRad() ) continue;
	break;
      }


      auto pid   = aPart->GetPID();
      //      auto pid   = aPart->GetPIDTight();
      auto pt    = aPart->GetRotatedMomentum().Pt()/1000.;
      // auto rapid  = aPart->GetRapiditycm();;
      // auto rapidn = rapid / y_norm;

      auto rapid  = aPart->GetRapidity();
      rapidn      = rapid / y_norm - 1;

      auto fmass  = aPart->GetMass();
      Double_t u_t0  = aPart->GetRotatedMomentum().Pt()/fmass/u_p;

      TLorentzVector lrnzVec = aPart->GetLorentzVector();
      TVector3 boostVec = STLorentzBoostVector::GetBoostVector(4+isys);
      lrnzVec.Boost(-boostVec);
      TLorentzVector transvec(lrnzVec.Z(),0,lrnzVec.Pt(),lrnzVec.E());
      TLorentzVector transvecx(lrnzVec.Z(),lrnzVec.Pt()*sin(rnd_Phi), lrnzVec.Pt()*cos(rnd_Phi), lrnzVec.E());
      Double_t yx = transvecx.Rapidity();

      // TLorentzVector transvecx = transvec;
      // transvecx.RotateZ(rnd_Phi);
      // Double_t yx = calc_yx(transvecx);

      Int_t hpid = -1;
      for( auto i: ROOT::TSeqI(pidmax) ){
	if( pid == embedPartID[i] ) {
	  hpid = i;
	  break;
	}
      }

      if( hpid > -1 ) {
	
	for( auto i : ROOT::TSeqI(mmax) ){
	  if( btrkcut[i] ) {
	    hyuteff[hpid][i]  -> Fill( rapidn, u_t0);
	    hypteff[hpid][i]  -> Fill( rapidn, pt);
	    hyyxeff[hpid][i]  -> Fill( rapidn, yx/y_norm);
	    hyxheff[hpid][i]  -> Fill( rapidn, abs(yx)/y_norm);
	    hyyteff[hpid][i]  -> Fill( rapidn, transvec.Rapidity()/y_norm);
	    hmtrack1[hpid][i] -> Fill(aflow->mtrack1);
	    hmtrack2[hpid][i] -> Fill(aflow->mtrack2);
	    hmtrack4[hpid][i] -> Fill(aflow->mtrack4);
	  }
	}
      }
    }
  }


  outFile->Write();

  TCanvas *ccv;
  Int_t id;
  


  if( 0 ) {
    ccv = new TCanvas("ccv1","ccv1",1000,1000);
    ccv->Divide(4,mmax);
    id = 1;
  
    for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
	ccv->cd(id); id++;
	hyyxeff[i][j]->Draw("colz");
      }

    ccv = new TCanvas("ccv2","ccv2",1000,1000);
    ccv->Divide(pidmax,mmax);
    id = 1;

    for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
	ccv->cd(id); id++;
	hmtrack4[i][j]->Draw();
     }

    ccv = new TCanvas("ccv3","ccv3",1000,1000);
    ccv->Divide(4,mmax);
    id = 1;

    for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
	ccv->cd(id); id++;
	hmtrack1[i][j]->Draw();
      }


    ccv = new TCanvas("ccv4","ccv4",1000,1000);
    ccv->Divide(4,mmax);
    id = 1;

    for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
	ccv->cd(id); id++;
	hyxheff[i][j]->Draw("colz");
      }
  }
}

