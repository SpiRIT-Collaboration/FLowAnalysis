#include "openRunAna.C" 
#include "DoAna.h"
#include "../tasks/STLorentzBoostVector.hh"

auto calc_yx = [](TLorentzVector v)
{
  double beta_x = v.Beta()*TMath::Sin(v.Theta())*TMath::Cos(v.Phi());
  return 0.5*TMath::Log((1.+beta_x)/(1.-beta_x));
};


void PlotTrackQuality(UInt_t azmcut=5)
{

  gROOT->Reset();
  openRunAna();

  TDatime dtime;
  TRandom3 grnd(dtime.GetSecond());
  gRandom->SetSeed(dtime.GetSecond());

  std::map< Long64_t, Int_t > PDGtoIDX{{2212, 0},{1000010020,1},{1000010030,2},{1000020030,3},{1000020040,4},{2112,7}};

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
  //  UInt_t  phicutID = 5;
  UInt_t  phicutID = azmcut;
  UInt_t  cutndf   = 15;
  Float_t cutdist  = 20.;

  TString phiname[] = { "45or135", "45", "135", "45to135","all","-30to20"};

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
  TH2D *hypxeff[pidmax][mmax];
  TH2D *hyeteff[pidmax][mmax];
  TH2D *hybgeff[pidmax][mmax];
  TH1I *hmtrack1[pidmax][mmax];
  TH1I *hmtrack2[pidmax][mmax];
  TH1I *hmtrack4[pidmax][mmax];

  TH2D *hypt[pidmax];
  TH2D *hangle[pidmax];
  
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)     
  
  for(auto i: ROOT::TSeqI(pidmax))for(auto j: ROOT::TSeqI(mmax) ){

      TString hlabel = Form("hyuteff_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hyuteff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50, 0., 2.);
      hlabel         = Form("hypteff_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hypteff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50, 0., 2000.);
      hlabel         = Form("hyyxeff_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hyyxeff[i][j]  = new TH2D(hlabel, hlabel , 40,-2.,2.,40,-2.,2.);
      hlabel         = Form("hyxheff_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hyxheff[i][j]  = new TH2D(hlabel, hlabel , 40,-2.,2.,20,0.,2.);
      hlabel         = Form("hyyteff_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hyyteff[i][j]  = new TH2D(hlabel, hlabel , 40,-2.,2.,20,0.,2.);
      hlabel         = Form("hypxeff_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hypxeff[i][j]  = new TH2D(hlabel, hlabel , 40,-2.,2.,60,-1200.,1200.);
      hlabel         = Form("hyeteff_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hyeteff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50,    0., 500.);
      hlabel         = Form("hybgeff_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hybgeff[i][j]  = new TH2D(hlabel, hlabel , 40, -2., 2., 50,    0., 1.2);

      hlabel         = Form("hmtrack1_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hmtrack1[i][j] = new TH1I(hlabel, hlabel, 80, 0, 80);
      hlabel         = Form("hmtrack2_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hmtrack2[i][j] = new TH1I(hlabel, hlabel, 80, 0, 80);
      hlabel         = Form("hmtrack4_"+ncls[i].sName+"_%dto%d",trkcut[j][0],trkcut[j][1]);
      hmtrack4[i][j] = new TH1I(hlabel, hlabel, 80, 0, 80);

    }

  for(auto i: ROOT::TSeqI(pidmax)) {
    TString hlabel     = "hypt_"+ncls[i].sName;
    hypt[i]    = new TH2D(hlabel, hlabel , 50, 0., 1., 50, 0., 2000.);
    hlabel     = "hangle_"+ncls[i].sName;
    hangle[i]  = new TH2D(hlabel, hlabel ,100,0,1.6,100,-3.15,3.15);
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
      UInt_t Lm = trkcut[j][0];
      UInt_t Um = trkcut[j][1];
      if( isys == 0 || isys == 2 ) {
	Lm+=1; Um+=1;
      }
      if( trackselection >= Lm && trackselection < Um ) {
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

      auto pid   = aPart->GetPID();
      Int_t hpid = PDGtoIDX[pid];

      if( hpid > -1 && hpid < 5 ) {
	auto rapid  = aPart->GetRapidity();
	auto pt    = aPart->GetRotatedMomentum().Pt();
	// auto rapid  = aPart->GetRapiditycm();;
	TLorentzVector lrnzVec = aPart->GetLorentzVector();

	if(trackselection>54) {
	  hypt[hpid]    -> Fill( rapid, pt);
	  hangle[hpid]  -> Fill( lrnzVec.Theta(), lrnzVec.Phi());
	}

	// phi cut
	auto phi    = aPart->GetRotatedMomentum().Phi()*TMath::RadToDeg();
	// w phicut                                                                                                                                              
	switch(phicutID) {
	case 0:
	  if( abs(phi) > 45 && abs(phi) < 135) continue;
	  break;
	case 1:
	  if( abs(phi) > 45) continue;
	  break;
	case 2:
	  if( abs(phi) < 135 ) continue;
	  break;
	case 3:
	  if( abs(phi) < 45 || abs(phi) > 135) continue;
	  break;
	case 5:
	  if( phi > 20 || phi < -30 ) continue;
	  break;
	}


	// auto rapidn = rapid / y_norm;
	auto rapidn = rapid / y_norm - 1;
	auto fmass  = aPart->GetMass();
	Double_t u_t0  = aPart->GetRotatedMomentum().Pt()/fmass/u_p;

	TVector3 boostVec = STLorentzBoostVector::GetBoostVector(4+isys);
	lrnzVec.Boost(-boostVec);
	TLorentzVector transvec(lrnzVec.Z(),0,lrnzVec.Pt(),lrnzVec.E());
	TLorentzVector transvecx(lrnzVec.Z(),lrnzVec.Pt()*sin(rnd_Phi), lrnzVec.Pt()*cos(rnd_Phi), lrnzVec.E());
	Double_t yx = transvecx.Rapidity();

	// TLorentzVector transvecx = transvec;
	// transvecx.RotateZ(rnd_Phi);
	// Double_t yx = calc_yx(transvecx);



	auto MassNumber = ncls[hpid].A;

	for( auto i : ROOT::TSeqI(mmax) ){
	  if( btrkcut[i] ) {
	    hyuteff[hpid][i]  -> Fill( rapidn, u_t0);
	    hypteff[hpid][i]  -> Fill( rapidn, pt);
	    hyyxeff[hpid][i]  -> Fill( rapidn, yx/y_norm);
	    hyxheff[hpid][i]  -> Fill( rapidn, abs(yx)/y_norm);
	    hyyteff[hpid][i]  -> Fill( rapidn, transvec.Rapidity()/y_norm);
	    hypxeff[hpid][i]  -> Fill( rapidn, aPart->GetRotatedMomentum().Pt()*cos(rnd_Phi) / MassNumber );
	    hyeteff[hpid][i]  -> Fill( rapidn, (lrnzVec.E()-lrnzVec.M())*sin(lrnzVec.Theta()));
	    hybgeff[hpid][i]  -> Fill( rapidn, lrnzVec.Pt()/lrnzVec.M());
	    hmtrack1[hpid][i] -> Fill(aflow->mtrack1);
	    hmtrack2[hpid][i] -> Fill(aflow->mtrack2);
	    hmtrack4[hpid][i] -> Fill(aflow->mtrack4);
	  }
	}


      }
    }
  }

  outFile->Write();


  auto accfile = TFile::Open("TPCacceptance.root","recreate");
  for( auto i: ROOT::TSeqI(pidmax)) {
    hypt[i] -> Write();
    hangle[i] -> Write();
  }
  accfile->Close();



  TCanvas *ccv;
  Int_t id;
  


  if( 1 ) {
    ccv = new TCanvas("ccv1","ccv1",1000,1000);
    ccv->Divide(mmax,pidmax);
    id = 1;
  
    for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
	ccv->cd(id); id++;
	hyeteff[i][j]->Draw("colz");
      }

    ccv = new TCanvas("ccv1","ccv1",1000,1000);
    ccv->Divide(mmax,pidmax);
    id = 1;
  
    for( auto i: ROOT::TSeqI(pidmax))for( auto j: ROOT::TSeqI(mmax) ){
	ccv->cd(id); id++;
	hybgeff[i][j]->Draw("colz");
      }
  }
}


