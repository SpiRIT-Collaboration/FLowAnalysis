#include "DoRPRes.C" //#include "openRunAna.C" #include "DoFlow.h" #include "SimFunction.C"

//drawing

UInt_t selReactionPlanef = 10000;

TFile* fefffile;

Bool_t bccPsi;

// functions
TString SetupOutputFile(TString fopt);
void PlotPt(UInt_t selid = 2);
 
//-------------------//
void DoFlow_mod(Int_t isel = 2) 
{
  gROOT->Reset();

  openRunAna();

  if(rChain != NULL)     
    LOG(INFO) << " DoFLow_mod: System " << isys << "  -> " << sysName << FairLogger::endl; 

  else
    exit(0);


  gROOT->ProcessLine(".! grep -i void DoFlow_mod.C ");

  // Configuration ================================
  TString su = gSystem -> Getenv("UC");
  if( su != "" ) {
    Ucent = (UInt_t)atoi(su);
    if( Ucent < 0 )
      Ucent = 80;
  }


  su = gSystem -> Getenv("LC");
  if( su != "" ) {
    Lcent = (UInt_t)atoi(su);
    if( Lcent > Ucent || Lcent < 0)
      Lcent = 0;
  }

  LOG(INFO) << " Multiplicity :: " << Lcent << " to " << Ucent << FairLogger::endl;

  //--------------------------------------------------
  TString rpBase = gSystem -> Getenv("RPBS");
  RPBase = rpBase != "" ? atoi(rpBase): 0;

  LOG(INFO) << " Reaction plane base is " << RPBase << FairLogger::endl;

  TString ccPsi = gSystem -> Getenv("CCPSI");
  bccPsi = ccPsi == "" ? kTRUE: kFALSE;
  LOG(INFO) << " Angle dependent correction is ";
  if( bccPsi )
    LOG(INFO) << " !!! ON !!!"  << FairLogger::endl;
  else
    LOG(INFO) << " OFF " << FairLogger::endl;


  //==================================================
  if( isel > -1 )
    LOG(INFO) << " Output Version v" << sVer << "." << oVer << FairLogger::endl;

  if( isel > 0 ) 
    PlotPt((UInt_t)isel);  
  else if( isel == -1 )
    PlotPt(0);  
}

TString SetupOutputFile(TString fopt)
{
  //--------------------------------------------------
  //----- SETUP OUTPUT FILE --------------------------
  //--------------------------------------------------

  gSystem->cd("data");
  TString oVerSub= gSystem -> Getenv("OSUBV");
  TString fName = fopt + oVerSub;

  if( oVerSub == "" ) {
    UInt_t kVer = 0;
    TString checkName;

    while( kTRUE ) {
      checkName = Form(fName + "%d.root", kVer);
      if( !gSystem->FindFile(".",checkName) ) {
        break;
      }
      else
        kVer++;
    }
    fName = Form(fName + "%d.root", kVer);
  }
  else {
    TString checkName = fName + ".root";;
    if( !gROOT->IsBatch() && gSystem->FindFile(".",checkName) ) {
      LOG(INFO) << checkName << " is existing. Do you recreate? (y/n)" << FairLogger::endl;
      TString sAns;
      std::cin >> sAns;
      if(sAns == "N" || sAns == "n") {
	LOG(INFO) << " Retry" << FairLogger::endl;
	exit(0);
      }
    }
    fName += ".root";
  }

  LOG(INFO) << "File " << fName << " will be created. " << FairLogger::endl;

  return fName;
}
//--------------------------------------------------       

//pdt
void PlotPt(UInt_t selid = 2)       //%% Executable :
{
  TDatime beginTime;
  TDatime dtime;

  gStyle->SetOptStat(0);

  TString fHeader = "mod_"+ sysName + "_" + partname[selid]+".v"+sVer+".";
  auto fName = SetupOutputFile( fHeader );
  TFile *fopen = TFile::Open(fName,"recreate");
  
  //------------------------------
  TH1D* hpt[4];
  TH2D* hypt[4];
  TH2D* hthetapt[4];

  Float_t phirange[] = { TMath::Pi()*7/4, TMath::Pi()/4, TMath::Pi()*3/4, TMath::Pi()*5/4, TMath::Pi()*7/4 };  

  for( UInt_t i = 0; i < 4; i++ ) {
    TString label_low = Form("Phi:%4.0f",phirange[i]*TMath::RadToDeg());
    TString label_up  = Form(" to %4.0f",phirange[i+1]*TMath::RadToDeg());

    hpt[i]  = new TH1D(Form("hpt_%d",i) ,label_low + label_up + "; Pt",100,0.,1500.);
    hypt[i] = new TH2D(Form("hypt_%d",i),label_low + label_up + "; Rapidity; Pt", 100,0.,1.2,100.,0.,800);
    hthetapt[i] = new TH2D(Form("hthetapt_%d",i),label_low + label_up + ";#Theta; Pt",100,0.,1.4,100,0.,800.);
  }

  //------------------------------
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)

  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  Long64_t nEntry = SetBranch();
  
  for(Int_t i = 0; i < nEntry; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    /// for reaction plane resolution
    Bool_t bFill = kFALSE;
    Bool_t bRes  = kFALSE;

    beginTime.Copy(dtime);
    if(i%(UInt_t)(nEntry/50) == 0) {
      dtime.Set();
      Int_t ptime = dtime.Get() - beginTime.Get();

      LOG(INFO) << "Processing .... " 
		<< setw(4) << Int_t(((Double_t)i/(Double_t)nEntry)*100.) << " % = "
		<< setw(8) << i << "/"<< nEntry
		<< "--->"
		<< dtime.AsString() << " ---- "
		<< FairLogger::endl;
    }

    /// centrality selection
    if(aflow->mtrack2 > Ucent || aflow->mtrack2 <= Lcent || aflow->mtrack4 < 6) continue;

    //--------------------------------------------------
    //----- Main loop 
    //--------------------------------------------------
    TIter next(aArray);
    STParticle *aPart = NULL;
    UInt_t mtk = 0;
    while( (aPart = (STParticle*)next()) ) {

      mtk++;
      // ---- track quality selection ---
      if( isys != 5 && aPart->GetGoodTrackFlag() != 1111 ) continue;
      //------------------------------

      auto yaw   = aPart->GetYawAngle();
      auto pitch = aPart->GetPitchAngle();
      auto pt    = aPart->GetRotatedMomentum().Pt();
      auto mom   = aPart->GetRotatedMomentum().Mag();
      auto bmass = aPart->GetBBMass();
      auto phi   = aPart->GetRotatedMomentum().Phi();
      auto theta = aPart->GetRotatedMomentum().Theta();
      auto charge= aPart->GetCharge();
      auto rapid = aPart->GetRapidity();

      auto pid   = aPart->GetPID();
      // auto pid   = aPart->GetPIDLoose();  
      //auto pid   = aPart->GetPIDNorm();  

      //----- Particle Selection -----
      Bool_t bpart = kFALSE;
      auto MassNumber = partA[selid];
      if( selid == 8 ) {   //H
	if( pid == partid[2] ) {
	  MassNumber = partA[2];
	  bpart = kTRUE;
	}
	else if( pid == partid[3] ){
	  MassNumber = partA[3];
	  bpart = kTRUE;
	}
	else if( pid == partid[4] ) {
	  MassNumber = partA[4];
	  bpart = kTRUE;
	}
      }
      else if( selid == 9 ){ 
	if( pid == partid[4] && mom > 1000  ) {
	  pid = partid[4];
	  MassNumber = partA[4];
	  bpart = kTRUE;
	}
      }
      else if( selid == 0 && charge == 1 && pid == partid[0] ) //pi+
	bpart = kFALSE;

      else if( selid == 1 && charge == -1 && pid == partid[0] ) //pi-
	bpart = kFALSE;

      else if( selid < 8  && pid == partid[selid] && charge == 1 ){ //p,d,t,3He,4He
	MassNumber = partA[selid];
	bpart = kTRUE;
      }
   
      //@@****************************************
      if( !bpart ) continue;

      if( theta > 45.*TMath::DegToRad() ) continue;

      // if( abs( phi ) > 30.*TMath::DegToRad() ) continue;
      //if( abs( phi ) < 150.*TMath::DegToRad() ) continue;
      //*********************************************

      if( abs(phi) < phirange[1] ) {
	hpt[0]  -> Fill(pt);
	hypt[0] -> Fill(rapid, pt);
	hthetapt[0] -> Fill(theta, pt);
      }
      else if( phi >= phirange[1] && phi < phirange[2] ) {
	hpt[1] -> Fill(pt);
      	hypt[1] -> Fill(rapid, pt);
	hthetapt[1] -> Fill(theta, pt);
      }
      else if( abs(phi) > phirange[2] ) {
	hpt[2] -> Fill(pt);
	hypt[2] -> Fill(rapid, pt);
	hthetapt[2] -> Fill(theta, pt);
      }
      else if( TVector2::Phi_0_2pi(phi) < phirange[4] && TVector2::Phi_0_2pi(phi) >= phirange[3] ) {
	hpt[3] -> Fill(pt);
	hypt[3] -> Fill(rapid, pt);
	hthetapt[3] -> Fill(theta, pt);
      }
    }
  }

  ///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for( UInt_t i = 0; i < 4; i++ ) {
    hpt[i] -> Write();
    hypt[i]-> Write();
    hthetapt[i] -> Write();
  }

  //--------------------------------------------------
  //--- Plotting
  //--------------------------------------------------

  TCanvas* cc;
  UInt_t id=1; UInt_t ic = 0;


  //pt
  cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500); ic++;
  auto hptleg = new TLegend(0.55,0.9,0.55,0.9,"");

  // hpt[1] -> SetLineColor(6);
  // hpt[1] -> SetNormFactor(100);
  // hpt[1] -> Draw();


  hpt[2] -> SetLineColor(4);
  hpt[2] -> SetNormFactor(100);
  hpt[2] -> Draw();

  hpt[0] -> SetLineColor(2);
  hpt[0] -> SetNormFactor(100);
  hpt[0] -> Draw("same");

  // hpt[3] -> SetLineColor(7);
  // hpt[3] -> SetNormFactor(100);
  // hpt[3] -> Draw("same");

  hptleg -> AddEntry(hpt[0],hpt[0]->GetTitle());
  //  hptleg -> AddEntry(hpt[1],hpt[1]->GetTitle());
  hptleg -> AddEntry(hpt[2],hpt[2]->GetTitle());
  //  hptleg -> AddEntry(hpt[3],hpt[3]->GetTitle());
  hptleg -> Draw();

  // yPt
  cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,800); ic++;
  cc->Divide(2,2);

  for( UInt_t i = 0; i < 4; i++ ) {
    cc->cd(i+1);
    hypt[i] -> Draw("colz");
  }

}
    

    
  
    
