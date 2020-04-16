#include "openRunAna.C"
#include "DoFlow.h"
#include "SimFunction.C"
#include "sslib.h"
#include "Complex.C"
#include "poly.C"

//--- <cos[*]dphi> = 
//                    0:<cos>  1:x       2:x^2     3:x^3       4:x^4     5:x^5
Double_t ppar[2][6] = {{-1.,     0.626657, 0.,       -0.09694,   0.02754, -0.002283},  //k=1
		       {-1.,     0.,       0.25,     -0.011414, -0.034726, 0.006815} }; //k=2


//drawing

//multiplicity dependent correction factor
ROOT::Math::Interpolator *itrpvx;
std::vector< Double_t > v1x;
std::vector< Double_t > v1y;
std::vector< Double_t > v1xe;
std::vector< Double_t > v1ye;
std::vector< Double_t > v2x;
std::vector< Double_t > v2y;
std::vector< Double_t > v2xe;
std::vector< Double_t > v2ye;

TDatime beginTime;
TDatime dtime;

// functions
TString  GetPsiRPFileName(UInt_t index = 0);
Double_t GetRPBaseAngle(STFlowInfo *aflow);
Double_t GetError(Double_t x, Double_t y, Double_t xe, Double_t ye);
void     GetResolution();
void     CentralityDependence();            
void     PsiAngleDependence();      
UInt_t   GetPsiRPIndex(Double_t aVal);
UInt_t   Get2PsiRPIndex(Double_t aVal);
UInt_t   GetRPCorrIndex(Double_t mult);
Double_t *GetRPResolutionwChi(Double_t *rpres, TH1D *hphi0_180, TH1D *hphi90_180, UInt_t kk);
Double_t *GetRPResolutionwChi(Double_t *rpres, Double_t chi, Double_t chie, const UInt_t kk);
Double_t *GetRPResolutionwCount(Double_t *rpres, Double_t c0, Double_t c90, UInt_t kk);
Double_t *GetRPResolutionMCos(Double_t *rpres, Double_t mcos, Double_t mcose, UInt_t kk);
void     GetChi2();
void     PsiAngleDependence_old();            //%% Executable :

//--


//-------------------//
void DoRPRes(Int_t isel = 0) 
{
  gROOT->Reset();

  openRunAna(1);

  if(rChain != NULL)     
    LOG(INFO) << " System " << isys << "  -> " << sysName << FairLogger::endl; 

  else
    exit(0);


  gROOT->ProcessLine(".! grep -i void DoRPRes.C ");


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

  //==================================================

  LOG(INFO) << " Output Version v" << oVer << FairLogger::endl;


  if( isel == 0 ) {
    CentralityDependence() ;
    PsiAngleDependence()   ;
  }
  else if( isel == 1 ) {
    //PsiAngleDependence_old()   ;
    PsiAngleDependence()   ;
  }

  else if( isel == 2 ) 
    GetResolution();

  else if( isel == 3 )
    GetChi2();
}

//--------------------------------------------------------------------
void GetChi2()
{
  LOG(INFO) << " GetChi2" << FairLogger::endl;

  auto hcos_sub1ab  = new TH1D("hcos_sub1ab2","cos(a-b)",100,-1.,1.);
  auto hcos2_sub2ab = new TH1D("hcos2_sub2ab","cos2(a-b)",100,-1.,1.);

  auto hdphi10 = new TH1D("hdphi10","dphi0to180",100,0.,TMath::Pi());
  auto hdphi19 = new TH1D("hdphi19","dphi0to90 ",100,0.,TMath::Pi());

  auto hdphi20 = new TH1D("hdphi20","dphi0to180",100,0.,TMath::Pi());
  auto hdphi29 = new TH1D("hdphi29","dphi0to90 ",100,0.,TMath::Pi());

  const UInt_t nybin = 8;
  Float_t ybin[nybin];

  TH1D *hv2obs[nybin];
  for(UInt_t i = 0; i < nybin; i++) {
    ybin[i] = 0.1*(i+1) -0.4; 
    hv2obs[i] = new TH1D(Form("hv2obs_%d",i),"v2^{obs}",200,-1.,1.);
  }

  auto nEntry = SetBranch();  
  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  for(Long64_t i = 0; i < nEntry; i++){

    if(i%(UInt_t)(nEntry/50) == 0 || nEntry < 50) {
      dtime.Set();
      Int_t ptime = dtime.Get() - beginTime.Get();

      LOG(INFO) << "Processing .... " 
		<< setw(4) << Int_t(((Double_t)i/(Double_t)nEntry)*100.) << " % = "
		<< setw(8) << i << "/"<< nEntry
		<< "--->"
		<< dtime.AsString() << " ---- "
		<< FairLogger::endl;
    }

    rChain->GetEntry(i);

    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    if( aflow->mtrack4 < 3 ) continue;

    Double_t subevt_phi =  (aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi() ;
    hcos_sub1ab -> Fill( cos(subevt_phi) );

    hdphi10 -> Fill( abs(TVector2::Phi_mpi_pi(subevt_phi))  );
    if( abs(TVector2::Phi_mpi_pi(subevt_phi)) >= TMath::Pi()/2. )
      hdphi19 -> Fill( abs(TVector2::Phi_mpi_pi(subevt_phi))  );

    Double_t subevt_phi2 = 2.*( (aflow->unit2P_1fc).Phi()/2. - (aflow->unit2P_2fc).Phi()/2. );
    hcos2_sub2ab -> Fill( cos(subevt_phi2) );

    subevt_phi2 = 2.*abs(TVector2::Phi_mpi_pi( (aflow->unit2P_1fc).Phi()/2. - (aflow->unitP_2fc).Phi()));
    subevt_phi2 = abs(TVector2::Phi_mpi_pi(subevt_phi2));

    hdphi20 -> Fill( abs(TVector2::Phi_mpi_pi(subevt_phi2)) );
    if( abs(TVector2::Phi_mpi_pi(subevt_phi2)) <= TMath::Pi()/2. )
      hdphi29 -> Fill( abs(TVector2::Phi_mpi_pi(subevt_phi2)) );

    //--------------------------------------------------
    //----- Main loop 
    //--------------------------------------------------
    TIter next(aArray);
    STParticle *aPart = NULL;
    UInt_t mtk = 0;
    while( (aPart = (STParticle*)next()) ) {

      auto rapid  = aPart->GetRapiditycm();;
      auto pt     = aPart->GetRotatedMomentum().Pt();
      auto dphi   = aPart->GetAzmAngle_wrt_RP();
      auto dphi2  = aPart->GetAzmAngle2_wrt_RP();
      auto px     = pt*cos(dphi); 
      auto py     = pt*sin(dphi); 

      auto v2obs = pow(px/pt, 2) - pow(py/pt,2);

      for(UInt_t i = 0; i < nybin; i++){
	if( rapid < ybin[i] ) { 
	  hv2obs[i] -> Fill(v2obs);
	  break;
	}
      }
    }
  }

  Double_t cos2_res = sqrt(hcos2_sub2ab->GetMean())*2.;

  auto gu_v2 = new TGraph();
  for(UInt_t i = 0; i < nybin; i++ ){
    auto v2 = hv2obs[i]->GetMean()/cos2_res;

    gu_v2->SetPoint(i, ybin[i], v2);

    LOG(INFO) << " Y( " << ybin[i] << ") "
	      << hv2obs[i]->GetMean() << " = " 
	      << v2
	      << FairLogger::endl; 
  }

  Double_t *res = new Double_t[4];
  GetRPResolutionwChi(res, hdphi10, hdphi19, 2);
  LOG(INFO) << " chi_res <cos(dphi)> = " << *res
	    << " <cos(2dphi)> = " << *(res+2)
	    << FairLogger::endl; 

  GetRPResolutionwChi(res, hdphi20, hdphi29, 1);
  LOG(INFO) << " chi_res <cos(2dphi)> = " << *res
	    << FairLogger::endl; 

  id = 1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,700);
  cc->Divide(2,2);

  cc->cd(id); id++;
  hcos2_sub2ab->Draw();
  LOG(INFO) << "<cos2(psi2-psir)> =  sqrt(4* "
	    << hcos2_sub2ab->GetMean() 
	    << ") = "
	    << sqrt(hcos2_sub2ab->GetMean()*4)
 	    << FairLogger::endl; 


  cc->cd(id); id++;
  hcos_sub1ab->Draw();
  LOG(INFO) << "<cos(psi1-psir)> =  sqrt(4* "
	    << hcos_sub1ab->GetMean() 
	    << ") = "
	    << sqrt(hcos_sub1ab->GetMean()*4)
	    << FairLogger::endl; 

  cc->cd(id); id++;
  hdphi10->Draw();
  hdphi19->Draw("same");

  cc->cd(id); id++;
  hdphi20->Draw();
  hdphi29->Draw("same");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  gu_v2->Draw("ALP");

}


void GetResolution()
{
  LOG(INFO) << "GetResolution ... " << FairLogger::endl;

  // Averaged resolution is calculated.

  TH1I *hmult = new TH1I("hmult",";Multiplicity",100,0,100);
  TH1D *hphi0_180  = new TH1D("hphi0_180" , "0to180",100,0.,3.2);
  TH1D *hphi90_180 = new TH1D("hphi90_180","90to180",100,0.,3.2);

  TH1D *hdphi1  = new TH1D("hdphi1",";#Delta(#Psi)" ,100, -TMath::Pi(), TMath::Pi());
  TH1D *hdphi2  = new TH1D("hdphi2",";#Delta(2#Psi)",100, -TMath::Pi(), TMath::Pi());
  TH1D *hdphis  = new TH1D("hdphis",";#Delta(#Psi)" ,100, -TMath::Pi(), TMath::Pi());
  TH1D *hcosd   = new TH1D("hcosd" ,";cos#Delta(#Psi)" ,100, -1., 1.);
  TH1D *hcos2d  = new TH1D("hcos2d",";cos2#Delta(#Psi)" ,100, -1., 1.);

  auto nevt = SetBranch();  
  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  //  nevt = 2000;
  for(Long64_t i = 0; i < nevt; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    if( aflow->mtrack4 < 3 ) continue;

    Double_t subevt_phi  = TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi()-
						(aflow->unitP_2fc).Phi());

    hcosd -> Fill( aflow->cosdPsi);

    //Double_t subevt_phi2 = TVector2::Phi_mpi_pi( (aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi() );
    //Double_t subevt_phi2 = 2.*abs(TVector2::Phi_mpi_pi( (aflow->unit2P_1fc).Phi()/2. - (aflow->unitP_2fc).Phi()));
    Double_t subevt_phi2 = 2.*TVector2::Phi_mpi_pi( (aflow->unit2P_1fc).Phi()/2. - (aflow->unit2P_2fc).Phi()/2.);
    subevt_phi2 = TVector2::Phi_mpi_pi(subevt_phi2);

    hcos2d -> Fill( cos(subevt_phi2) );
      
    hdphi1 -> Fill( subevt_phi );
    hdphi2 -> Fill( subevt_phi2);

    hphi0_180->Fill( abs(subevt_phi) );

    if( abs(subevt_phi) > TMath::Pi()/2. )
      hphi90_180->Fill( abs(subevt_phi) );

    hmult->Fill( aflow->mtrack2 );

  }

  auto  rpres = new Double_t[2]; 
  GetRPResolutionwChi(rpres, hphi0_180, hphi90_180, 1);

  LOG(INFO) << " Ordinal way  <cos (d phi) > = " << *rpres << FairLogger::endl;

  // ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  // hmult->Draw();

  // ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  // hphi0_180->Draw();
  // hphi90_180->SetLineColor(2);
  // hphi90_180->Draw("same");

  // ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  // cc->Divide(1,2);
  // cc->cd(1);
  // hcosd->Draw();
  // cc->cd(2);
  //  hcos2d->Draw();


  LOG(INFO) << " New way --------- >>>  " << FairLogger::endl;


  auto mcos_sub = hcosd->GetMean();
  mcos_sub = sqrt(mcos_sub);
  auto mcos_err = hcosd->GetStdDev()/sqrt((Double_t)hcosd->GetEntries()-1 );

  rpres = GetRPResolutionMCos(rpres, mcos_sub, mcos_err, 1);
  auto chifull1  = sqrt(2) * rpres[0];
  auto chifull1e = sqrt(2) * rpres[1];

  rpres = GetRPResolutionwChi(rpres, chifull1, chifull1e, 1);

  LOG(INFO) << " <cosd > = " << rpres[0] << " +- " << rpres[1] 
	    << " chi_1 = " << chifull1 
	    << FairLogger::endl;

  rpres = GetRPResolutionwChi(rpres, chifull1, chifull1e, 2);

  LOG(INFO) << " <cos2d >m1 = " << rpres[0] << " +- " << rpres[1] 
	    << " chi_1 = " << chifull1 
	    << FairLogger::endl;


  LOG(INFO) << " <cos2(psi_a - psi_b)> ----> " << FairLogger::endl;  
  mcos_sub = hcos2d->GetMean();
  mcos_sub = sqrt(mcos_sub);
  mcos_err = hcos2d->GetStdDev()/sqrt( (Double_t)hcos2d->GetEntries()-1 );

  rpres = GetRPResolutionMCos(rpres,  mcos_sub, mcos_err, 2);
  auto chifull2  = sqrt(2)*rpres[0];
  auto chifull2e = sqrt(2)*rpres[1];

  auto chifull  = (chifull1 + chifull2)/2.;
  auto chifulle = sqrt(pow(chifull1e,2) + pow(chifull2e,2));
  rpres = GetRPResolutionwChi(rpres, chifull, chifulle, 2);

  LOG(INFO) << "###" << FairLogger::endl;
  LOG(INFO) << " <cos 2dsub = " << mcos_sub 
	    << " <cos2d >m2 = " << rpres[0] << " +- " << rpres[1] 
	    << " chi_2 = "      << chifull << " +- " << chifulle
	    << FairLogger::endl;
  
}
//--------------------------------------------------
//@@@@@
void PsiAngleDependence()            //%% Executable :
{
  LOG(INFO) << " PsiAngleDependence .... " << FairLogger::endl;

  gROOT->Reset();
  gROOT->cd();


  TString fName = GetPsiRPFileName(0);
  auto GraphSave = new TFile(fName,"recreate");

  TH1I *hmult  = new TH1I("hmult" ,"multiplicity",100,0,100);
  TH1I *hmult1 = new TH1I("hmult1","multiplicity",100,0,100);
  TH1I *hmult2 = new TH1I("hmult2","multiplicity",100,0,100);
  auto hpsindx = new TH2D("hpsindx",";#Psi ;index ",100,-TMath::Pi(),TMath::Pi(),20,0.,20.);

  TH1D *hpsi1[npsi];
  TH1D *hpsi2[npsi];
  TH1D *hcos1_sub1ab[npsi];
  TH1D *hcos2_sub2ab[npsi];
  TH1D *hdpsi_sub1ab[npsi];
  TH1D *hdpsi_sub2ab[npsi];


  for(UInt_t k = 0; k < npsi; k++){
    hpsi1[k]        = new TH1D(Form("hpsi1_%d",k),"",100,-3.15,3.15);
    hpsi2[k]        = new TH1D(Form("hpsi2_%d",k),"",100,-3.15,3.15);
    hcos1_sub1ab[k] = new TH1D(Form("hcos1_sub1ab_%d",k),"<cos(phi)>",100,-1.,1.);
    hcos2_sub2ab[k] = new TH1D(Form("hcos2_sub2ab_%d",k),"<cos(2phi)>",100,-1.,1.);
    hdpsi_sub1ab[k]   = new TH1D(Form("hdpsi_sub1ab_%d",k),"#Delta #Psi",100.,-TMath::Pi(), TMath::Pi());
    hdpsi_sub2ab[k]   = new TH1D(Form("hdpsi_sub2ab_%d",k),"#Delta #Psi",100.,-TMath::Pi(), TMath::Pi());
  }

  auto hiphi = new TH1I("hiphi","hiphi",15,0,15);


  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  Long64_t nEntry = SetBranch();

  //  nEntry = 1000;
  for(Int_t i = 0; i < nEntry; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    /// for reaction plane resolution  
    Bool_t bFill = kFALSE;
    Bool_t bRes  = kFALSE;

    if(i%(UInt_t)(nEntry/50) == 0 || nEntry < 50) {
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

    hmult ->Fill( aflow->mtrack4 );
    hmult1->Fill( aflow->mtrack1 );
    hmult2->Fill( aflow->mtrack2 );

    bRes = kTRUE; //@@@ 
    TIter next(aArray);
    STParticle *aPart = NULL;

    auto RPangle = GetRPBaseAngle(aflow);
    UInt_t iphi  = GetPsiRPIndex(RPangle);

    hpsi1[iphi] ->Fill( RPangle );

    //    Double_t subevt_phi1 = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi()));
    Double_t subevt_phi1 = TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi());
    hcos1_sub1ab[iphi] -> Fill( cos(subevt_phi1) );
    hdpsi_sub1ab[iphi] -> Fill( subevt_phi1 );

    hpsi2[iphi] ->Fill( RPangle );
    Double_t subevt_phi2 = 2.*abs(TVector2::Phi_mpi_pi( (aflow->unit2P_1fc).Phi()/2. - (aflow->unit2P_2fc).Phi()/2.));
    hcos2_sub2ab[iphi] -> Fill(cos( subevt_phi2 ) );
    hdpsi_sub2ab[iphi] -> Fill( TVector2::Phi_mpi_pi(2.* subevt_phi2 ) );

  }


  //-----------
  auto gv_psi1 = new TGraphErrors();
  gv_psi1->SetName("gv_psi1");
  gv_psi1->SetTitle("; #psi; <cos(#Delta #Psi)>");

  auto gv_dpsi1 = new TGraphErrors();
  gv_dpsi1->SetName("gv_dpsi1");
  gv_dpsi1->SetTitle("; #psi; <cos(#Delta #Psi)>");

  auto gv_psi2 = new TGraphErrors();
  gv_psi2->SetName("gv_psi2");
  gv_psi2->SetTitle("; #psi; <cos(#Delta 2#Psi)>");
   

  UInt_t id = 1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1000);
  cc->Divide(2, npsi);

  UInt_t ip = 0; UInt_t ip2 = 0;
  Double_t *rpres = new Double_t[2];

  UInt_t nst[] = {6, 0};
  for(UInt_t j = 0; j < 2; j++) {
    for(UInt_t i = nst[j]; i < nst[j]+6; i++) {
      
      if( hpsi1[i]->GetEntries() < 5 ) continue;

      Double_t psi = hpsi1[i]->GetMean();

      Double_t meancos_sub = hcos1_sub1ab[i]->GetMean();
      Double_t meancos_err = hcos1_sub1ab[i]->GetStdDev()/sqrt((Double_t)(hcos1_sub1ab[i]->GetEntries()-1));
      meancos_sub  = sqrt( meancos_sub );
      meancos_err  = 0.5 * meancos_err/meancos_sub ;

      LOG(INFO) << " <cos(dpsi)>[" << i << "] = " << meancos_sub << " +- " << meancos_err << FairLogger::endl;
      
      rpres = GetRPResolutionMCos(rpres, meancos_sub, meancos_err, 1);
      Double_t chifull  = sqrt(2) * rpres[0];
      Double_t chifulle = sqrt(2) * rpres[1];

      rpres = GetRPResolutionwChi(rpres, chifull, chifulle, 1);

      gv_psi1->SetPoint(ip, psi, rpres[0]);
      gv_psi1->SetPointError(ip, 0., rpres[1]);

      
      auto nc90 = hdpsi_sub1ab[i] -> Integral(0,25) + hdpsi_sub1ab[i] -> Integral(76, 100);
      auto nc0  = hdpsi_sub1ab[i] -> GetEntries();
      
      rpres = GetRPResolutionwCount(rpres, nc0, nc90, 1);
      gv_dpsi1->SetPoint(ip, psi, rpres[0]);
      gv_dpsi1->SetPointError(ip, 0., rpres[1]);

      ip++;

      
      if( hpsi2[i]->GetEntries() > 0 ){
	rpres = GetRPResolutionwChi(rpres, chifull, chifulle, 2);

	psi = hpsi2[i]->GetMean();
	gv_psi2->SetPoint(ip2, psi, rpres[0]);
	gv_psi2->SetPointError(ip2, 0., rpres[1]);

	ip2++;
      }

      cc->cd(id); id++;
      hdpsi_sub1ab[i]->Draw();


      cc->cd(id); id++;
      hpsi2[i]->Draw();


    }
  }

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));  
  gv_psi1->SetLineColor(2);
  gv_psi1->Draw("ALP");

  gv_dpsi1->SetLineColor(4);
  gv_dpsi1->Draw("same");


  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));  
  gv_psi2->SetLineColor(2);
  gv_psi2->Draw("ALP");


  gv_psi1->Write();
  gv_psi2->Write();
  gv_dpsi1->Write();

  GraphSave->Write();

  LOG(INFO) << GraphSave->GetName() << " is created. " << FairLogger::endl;

} // void PsiAngleDependence() 




void PsiAngleDependence_old()            //%% Executable :
{
  LOG(INFO) << " PsiAngleDependence_old .... " << FairLogger::endl;

  gROOT->Reset();
  gROOT->cd();


  TString fName = GetPsiRPFileName(1);
  auto GraphSave = new TFile(fName,"recreate");

  TH1I *hmult  = new TH1I("hmult" ,"multiplicity",100,0,100);
  TH1I *hmult1 = new TH1I("hmult1","multiplicity",100,0,100);
  TH1I *hmult2 = new TH1I("hmult2","multiplicity",100,0,100);
  auto hpsindx = new TH2D("hpsindx",";#Psi ;index ",100,-TMath::Pi(),TMath::Pi(),20,0.,20.);

  TH1I *hphibin[npsi];
  TH1D *hphi0_180[npsi];
  TH1D *hphi90_180[npsi];
  TH1D *h2phi0_180[npsi];
  TH1D *h2phi90_180[npsi];
  TH1D *hpsi1[npsi];
  TH1D *hpsi2[npsi];
  TH1D *hcos1_sub1ab[npsi];
  TH1D *hcos2_sub2ab[npsi];

  for(UInt_t k = 0; k < npsi; k++){
    TString htitle = Form("hphi0_180_%d",k);
    hphi0_180[k]  = new TH1D(htitle, "",100,0.,3.2);
    htitle = Form("hphi90_180_%d",k);
    hphi90_180[k] = new TH1D(htitle,"",100,0.,3.2);
    hpsi1[k] = new TH1D(Form("hpsi1_%d",k),"",100,-3.15,3.15);

    htitle = Form("h2phi0_180_%d",k);
    h2phi0_180[k]  = new TH1D(htitle, "",100,0.,3.2);
    htitle = Form("h2phi90_180_%d",k);
    h2phi90_180[k] = new TH1D(htitle,"",100,0.,3.2);
    hpsi2[k] = new TH1D(Form("hpsi2_%d",k),"",100,-3.15,3.15);

    hcos1_sub1ab[k] = new TH1D(Form("hcos1_sub1ab_%d",k),"<cos(phi)>",100,-1.,1.);
    hcos2_sub2ab[k] = new TH1D(Form("hcos2_sub2ab_%d",k),"<cos(2phi)>",100,-1.,1.);

  }

  auto hiphi = new TH1I("hiphi","hiphi",15,0,15);


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

    if(i%(UInt_t)(nEntry/50) == 0 || nEntry < 50) {
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

    hmult ->Fill( aflow->mtrack4 );
    hmult1->Fill( aflow->mtrack1 );
    hmult2->Fill( aflow->mtrack2 );

    bRes = kTRUE; //@@@ 
    TIter next(aArray);
    STParticle *aPart = NULL;

    auto RPangle = GetRPBaseAngle(aflow);
    UInt_t iphi  = GetPsiRPIndex(RPangle);

    hiphi->Fill( iphi );
    hpsi1[iphi] ->Fill( RPangle );

    Double_t subevt_phi1 = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi()));
    hcos1_sub1ab[iphi] -> Fill(cos( subevt_phi1 ) );
    hphi0_180[iphi] ->Fill(subevt_phi1);

    if( subevt_phi1 > TMath::Pi()/2. )
      hphi90_180[iphi]->Fill(subevt_phi1);

    //iphi = Get2PsiRPIndex( aflow->unit2P_fc.Phi()/2. );
    
    hpsindx->Fill(aflow->unit2P_fc.Phi()/2., (Double_t)iphi);

    Double_t subevt_phi2 = 2.*abs(TVector2::Phi_mpi_pi( (aflow->unit2P_1fc).Phi()/2. - (aflow->unit2P_2fc).Phi()/2.));
    hcos2_sub2ab[iphi] -> Fill(cos( subevt_phi2 ) );

    subevt_phi2 = 2.*abs(TVector2::Phi_mpi_pi( (aflow->unit2P_1fc).Phi()/2. - (aflow->unitP_2fc).Phi()));
    subevt_phi2 = abs(TVector2::Phi_mpi_pi(subevt_phi2));


    //    hpsi2[iphi] ->Fill( aflow->unit2P_fc.Phi()/2. );
    hpsi2[iphi] ->Fill( RPangle );

    h2phi0_180[iphi] ->Fill(subevt_phi2);

    if( subevt_phi2 < TMath::Pi()/2. )
      h2phi90_180[iphi]->Fill(subevt_phi2);
  }

  auto gv_psi1 = new TGraphErrors();
  gv_psi1->SetName("gv_psi1");
  gv_psi1->SetTitle("; #psi; <cos(#Delta #Psi)>");

  auto gv_psi2 = new TGraphErrors();
  gv_psi2->SetName("gv_psi2");
  gv_psi2->SetTitle("; #psi; <cos(#Delta 2#Psi)>");
   
  auto gv_cospsi1 = new TGraphErrors();
  gv_cospsi1->SetName("gv_cospsi1");

  auto gv_cospsi2 = new TGraphErrors();
  gv_cospsi2->SetName("gv_cospsi2");

  UInt_t id = 1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1000);
  cc->Divide(2,npsi);

  UInt_t ip = 0; UInt_t ip2 = 0;
  Double_t *rpres = new Double_t[4];

  UInt_t nst[] = {6, 0};
  for(UInt_t j = 0; j < 2; j++) {
    for(UInt_t i = nst[j]; i < nst[j]+6; i++) {

      if( hphi0_180[i]->GetEntries() < 5 ) continue;

      GetRPResolutionwChi(rpres, hphi0_180[i], hphi90_180[i], 1);
    
      Double_t psi = hpsi1[i]->GetMean();
      gv_psi1->SetPoint(ip, psi, rpres[0]);
      gv_psi1->SetPointError(ip, 0., rpres[1]);

      Double_t mcos = hcos1_sub1ab[i]->GetMean();
      Double_t scos = 1./mcos* abs( hcos1_sub1ab[i]->GetStdDev()/mcos/(Double_t)hcos1_sub1ab[i]->GetEntries() );
      gv_cospsi1->SetPoint(ip, psi, sqrt(mcos*2.) );
      gv_cospsi1->SetPointError(ip, 0., scos);
      ip++;
      
      GetRPResolutionwChi(rpres, hphi0_180[i], hphi90_180[i], 2);

      if( hpsi2[i]->GetEntries() > 0 ){
	psi = hpsi2[i]->GetMean();
	gv_psi2->SetPoint(ip2, psi, rpres[0]);
	gv_psi2->SetPointError(ip2, 0., rpres[1]);

	mcos = hcos2_sub2ab[i]->GetMean();
	scos = 1./mcos* abs( hcos2_sub2ab[i]->GetStdDev()/mcos/(Double_t)hcos2_sub2ab[i]->GetEntries() );
	gv_cospsi2->SetPoint(ip2, psi, -sqrt(mcos*2.) );
	gv_cospsi2->SetPointError(ip2, 0, scos );

	ip2++;
      }


      h2phi0_180[i]->Write();

      cc->cd(id); id++;
      h2phi0_180[i]->Draw();
      h2phi90_180[i]->SetLineColor(2);
      h2phi90_180[i]->Draw("same");

      cc->cd(id); id++;
      hpsi2[i]->Draw();

    }
  }

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));  
  gv_psi1->SetLineColor(2);
  gv_psi1->Draw("ALP");
  gv_cospsi1->SetLineColor(4);
  gv_cospsi1->Draw("same");


  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));  
  gv_psi2->SetLineColor(2);
  gv_psi2->Draw("ALP");

  //  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));  
  gv_cospsi2->SetLineColor(4);
  gv_cospsi2->Draw("same");



  hmult->Write();
  hmult1->Write();
  hmult2->Write();
  hiphi->Write();
  gv_psi1->Write();
  gv_psi2->Write();
  gv_cospsi1->Write();
  gv_cospsi2->Write();

  LOG(INFO) << GraphSave->GetName() << " is created. " << FairLogger::endl;

} // void PsiAngleDependence() 

void CentralityDependence()            //%% Executable :
{
  LOG(INFO) << " CentralityDependence .... " << FairLogger::endl;

  gROOT->Reset();
  gROOT->cd();

  TString fName = "data/mlt_"+ sysName + ".v"+ oVer +".root";

  auto GraphSave = new TFile(fName,"recreate");

  
  TH1I *hmult  = new TH1I("hmult" ,"multiplicity",100,0,mrange[0]);
  TH1I *hmult1 = new TH1I("hmult1","multiplicity",100,0,mrange[0]);
  TH1I *hmult2 = new TH1I("hmult2","multiplicity",100,0,mrange[0]);
  TH1I *hmultbin1[mbin];
  TH1I *hmultbin4[mbin];
  TH1D *hphi0_180[mbin];
  TH1D *hphi90_180[mbin];

  TH1D *h2phi0_180[mbin];
  TH1D *h2phi90_180[mbin];

  for(UInt_t k = 0; k < mbin; k++){
    TString htitle = Form("hmultbin4_%d",k);
    hmultbin4[k]  = new TH1I(htitle,"",100,0,mrange[0]);
    htitle = Form("hmultbin1_%d",k);
    hmultbin1[k]  = new TH1I(htitle,"",100,0,mrange[0]);

    htitle = Form("hphi0_180_%d",k);
    hphi0_180[k]  = new TH1D(htitle,"",100,0.,3.2);
    htitle = Form("hphi90_180_%d",k);
    hphi90_180[k] = new TH1D(htitle,"",100,0.,3.2);

    htitle = Form("h2phi0_180_%d",k);
    h2phi0_180[k]  = new TH1D(htitle,"",100,0.,3.2);
    htitle = Form("h2phi90_180_%d",k);
    h2phi90_180[k] = new TH1D(htitle,"",100,0.,3.2);
  }

  Long64_t nEntry = SetBranch();

  nEntry=10;
  for(Long64_t i = 0; i < nEntry; i++) {

    rChain->GetEntry(i);

    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);
    if( aflow == NULL || aflow->mtrack4 <= 5 ) continue;
    hmult  -> Fill( aflow->mtrack4 );
    hmult1 -> Fill( aflow->mtrack1 );
    hmult2 -> Fill( aflow->mtrack2 );

    UInt_t ik = 0;
    while( ik < mbin ) {
      if( aflow->mtrack4 > mrange[ik] ) break;
      ik++;
    }
    
    hmultbin1[ik] -> Fill( aflow->mtrack2 );
    hmultbin4[ik] -> Fill( aflow->mtrack4 );

    Double_t subdphi = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi()));

    hphi0_180[ik]->Fill( subdphi );
    if( subdphi > TMath::Pi()/2. )
      hphi90_180[ik]->Fill( subdphi );

    subdphi = abs(TVector2::Phi_mpi_pi(2.*((aflow->unit2P_1fc).Phi()/2. - (aflow->unit2P_2fc).Phi()/2.)));

    h2phi0_180[ik]->Fill( subdphi );
    if( subdphi > TMath::Pi()/2. )
      h2phi90_180[ik]->Fill( subdphi );

  }
   

  auto gv_mcos1 = new TGraphErrors();
  gv_mcos1->SetName("gv_mcos1");
  gv_mcos1->SetTitle(";Multiplicity; <cos(#Delta #Psi)>");
  auto gv_mcos1m1 = new TGraphErrors();
  gv_mcos1m1->SetName("gv_mcos1m1");
  gv_mcos1m1->SetTitle(";Multiplicity; <cos(#Delta #Psi)>");


  auto gv_mcos2 = new TGraphErrors();
  gv_mcos2->SetName("gv_mcos2");
  gv_mcos2->SetTitle(";Multiplicity; <cos(2#Delta #Psi)>");
  auto gv_mcos2m1 = new TGraphErrors();
  gv_mcos2m1->SetName("gv_mcos2m1");
  gv_mcos2m1->SetTitle(";Multiplicity; <cos(2#Delta #Psi)>");

  auto cc80 = new TCanvas("cc80","cc80",500,1000);
  cc80->Divide(2,10);
   
  UInt_t id = 1;
  UInt_t ip = 0;
  for(UInt_t i = 0; i < mbin; i++) {
    if( hphi0_180[i]->GetEntries() < 5 ) continue;
    cc80->cd(id); id++;
    hphi0_180[i] ->Draw();
    hphi90_180[i]->Draw("same");

    cc80->cd(id); id++;
    hmultbin4[i]->Draw();
    
    Double_t *rpres = new Double_t[4];

    GetRPResolutionwChi(rpres, hphi0_180[i], hphi90_180[i], 1);
    LOG(INFO) << " Multiplicity " << hmultbin4[i]->GetMean() << " > " << mrange[i]
	      << " k=1  <cos(Phi)> = " << rpres[0] << " +- " << rpres[1] 
	      << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
	      << FairLogger::endl;

    GetRPResolutionwChi(rpres, h2phi0_180[i], h2phi90_180[i], 2);
    LOG(INFO) << " Multiplicity " << hmultbin4[i]->GetMean() << " > " << mrange[i]
	      << " k=2 <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
	      << FairLogger::endl;
     
    if( hphi90_180[i]->GetEntries() > 15 && !std::isnan(rpres[0]) && rpres[1] < 0.1) {
      gv_mcos1->SetPoint(ip, hmultbin4[i]->GetMean(), rpres[0]);
      gv_mcos1->SetPointError(ip, hmultbin4[i]->GetStdDev()/sqrt((Double_t)hmultbin4[i]->GetEntries()), rpres[1]);

      gv_mcos2->SetPoint(ip, hmultbin4[i]->GetMean(), rpres[2]);
      gv_mcos2->SetPointError(ip, hmultbin4[i]->GetStdDev()/sqrt((Double_t)hmultbin4[i]->GetEntries()), rpres[3]);

      gv_mcos1m1->SetPoint(ip, hmultbin1[i]->GetMean(), rpres[0]);
      gv_mcos1m1->SetPointError(ip, hmultbin1[i]->GetStdDev()/sqrt((Double_t)hmultbin1[i]->GetEntries()), rpres[1]);
      gv_mcos2m1->SetPoint(ip, hmultbin1[i]->GetMean(), rpres[2]);
      gv_mcos2m1->SetPointError(ip, hmultbin1[i]->GetStdDev()/sqrt((Double_t)hmultbin1[i]->GetEntries()), rpres[3]);
      ip++;
    }
  }


  auto cc81 = new TCanvas("cc81","cc81");
  gv_mcos1m1->Draw("ALP");

  auto cc82 = new TCanvas("cc82","cc82");
  gv_mcos2m1->Draw("ALP");

  gv_mcos1->Write();
  gv_mcos2->Write();  
  gv_mcos1m1->Write();
  gv_mcos2m1->Write();  
  GraphSave->Write();
  hmult->Write();
  hmult1->Write();
  hmult2->Write();

  LOG(INFO) <<  GraphSave->GetName() << " is created. " << FairLogger::endl;
 
}

///######################################################################
/// sub functions
///######################################################################
Double_t  GetError(Double_t x, Double_t y, Double_t xe, Double_t ye)
{
  if( y == 0 ) return 0.;
  Double_t xc = abs(x/y);
  return  xc * sqrt(pow(xe/x,2) + pow(ye/y,2));
}

Double_t GetRPBaseAngle(STFlowInfo *aflow) {
  auto RPangle   = aflow->unitP_fc.Phi();

  if( isys == 5 && RPBase == 2)
    RPangle   = RPPsi;
  else if( RPBase == 1 )
    RPangle   = aflow->unitP.Phi();      

  return RPangle;
}


//**************************************************
UInt_t GetPsiRPIndex(Double_t aVal)
{
  return UInt_t( TVector2::Phi_0_2pi(aVal)/dphi);
}

UInt_t Get2PsiRPIndex(Double_t aVal)
{
  //  return UInt_t( TVector2::Phi_0_2pi(aVal - TMath::Pi()/2.)/dphi);
  return UInt_t( TVector2::Phi_0_2pi(aVal)/dphi);
}

//**************************************************
UInt_t GetRPCorrIndex(Double_t aVal) //???
{
  UInt_t xn = v1x.size() - 1;


  if( xn == 1 )
    return 0;

  if( aVal <= v1x.at( v1x.size() - 1 ) ) {
    xn = (UInt_t)itrpvx->Eval( (Double_t)aVal ) ;
  
    if( abs(v1x.at(xn) - aVal) > abs(v1x.at(xn+1) - aVal) )
      xn += 1;
  }

  return xn;
}

//**************************************************
TString GetPsiRPFileName(UInt_t index = 0)
{
  TString fn;
  if( index == 1 )
    fn = "data/dpsi_"+ sysName + ".v" + oVer;
  else
    fn = "data/cpsi_"+ sysName + ".v" + oVer;


  if( RPBase < 3 && RPBase > 0)
    fn += Form("_%d", RPBase);

  fn += Form(".m%02dto%02d.root",Lcent,Ucent);

  return fn;
}

//**************************************************
Bool_t SetRPResolution()
{
  itrpvx = new ROOT::Math::Interpolator(20, ROOT::Math::Interpolation::kPOLYNOMIAL);

  TString fname = "data/mlt_"+ sysName + ".v" + sVer + ".root";
  TFile *fOpen = TFile::Open(fname);
  if( fOpen == NULL ) {
    LOG(ERROR) << "Please do it doflow -2 " << FairLogger::endl;
    return kFALSE;
  }
  
  LOG(INFO) << fname << " is opened. " << FairLogger::endl;

  auto hgv_mcos1 = (TGraphErrors*)fOpen->Get("gv_mcos1");
  auto hgv_mcos2 = (TGraphErrors*)fOpen->Get("gv_mcos2");

  std::vector< Double_t > ix;
    
  Double_t x, y, xe, ye;
  UInt_t k = 0;
  ix.push_back(k); 
  v1x.push_back(0.); v1xe.push_back(0.);
  v1y.push_back(0.); v1ye.push_back(0.);
  k++;
  for( Int_t i = (Int_t)hgv_mcos1->GetN()-1; i > -1; i-- ) {
    hgv_mcos1->GetPoint(i, x, y);
    xe = hgv_mcos1->GetErrorX(i);
    ye = hgv_mcos1->GetErrorY(i);


    ix.push_back(k); 
    v1x.push_back(x);
    v1y.push_back(y);
    v1xe.push_back(xe);
    v1ye.push_back(ye);

    k++;
  }
  if( v1x.size() > 1 )
    itrpvx->SetData(v1x,ix);


  v2x.push_back(0.); v2xe.push_back(0.);
  v2y.push_back(0.); v2ye.push_back(0.);

  k = 0;
  for( Int_t i = (Int_t)hgv_mcos2->GetN()-1; i > -1; i-- ) {
    hgv_mcos2->GetPoint(i, x, y);
    xe = hgv_mcos2->GetErrorX(i);
    ye = hgv_mcos2->GetErrorY(i);


    v2x.push_back(x);
    v2y.push_back(y);
    v2xe.push_back(xe);
    v2ye.push_back(ye);

    LOG(DEBUG) << "v2 resolution "<< k << " th " << v2x.at(k) << " vs " << v2y.at(k) << " +- " << v2ye.at(k) << FairLogger::endl; 
    k++;
  }

  fOpen->Close();

  return kTRUE;
}

//**************************************************
Double_t *GetRPResolutionwCount(Double_t *rpres, Double_t c0, Double_t c90, UInt_t kk) 
{

  if( rpres == NULL ) return rpres;

  if( c0 == 0 ) {
    LOG(INFO) << " 0 entry " << endl;
    rpres[0] = 1.;
    rpres[1] = 0.;
    return rpres;
  }


  Double_t mr   = c90/c0;
  Double_t mre2 = 1./c0; 

  Double_t chi   = sqrt(abs(-4.* log(2.* mr)));
  Double_t dchin = sqrt(-4.* log(2.*(mr - 1./sqrt(c90))));
  Double_t chie2 = pow( chi - dchin, 2 );
  Double_t chie = sqrt(chie2);


  Double_t par1[] = {0.626657, -0.09694, 0.02754, -0.002283};
  Double_t par2[] = {0.25,  -0.011414, -0.034726, 0.006815};

  //v1
  if( kFALSE ) {
    // //v1
    if( kk == 1 )
      rpres[0] = par1[0]*chi + par1[1]*pow(chi,3) + par1[2]*pow(chi,4) + par1[3]*pow(chi,5);  
    // //v2
    else if( kk == 2 )
      rpres[0] = par2[0]*pow(chi,2) + par2[1]*pow(chi,3) + par2[2]*pow(chi,4) + par2[3]*pow(chi,5);
  }
  else {
    if( kk == 1 )
      rpres[0] = sqrt(TMath::Pi())/(2.*sqrt(2))*chi*exp(-chi*chi/4.)*
	(ROOT::Math::cyl_bessel_i(0,chi*chi/4.)+ROOT::Math::cyl_bessel_i(1,chi*chi/4.));

    else if( kk == 2 )
      rpres[0] = sqrt(TMath::Pi())/(2.*sqrt(2))*chi*exp(-chi*chi/4.)*
	(ROOT::Math::cyl_bessel_i((kk-1.)/2.,chi*chi/4.)+ROOT::Math::cyl_bessel_i((kk+1.)/2.,chi*chi/4.));
  }

  Double_t err2 = 0;
  if( kk == 1 ) {
    err2 = mre2*pow(par1[0] + 3.*par1[1]*pow(chi,2) + 4.*par1[2]*pow(chi,3) + 5.*par1[3]*pow(chi,4),2);
    rpres[1] = sqrt( err2 );
  }
  else if( kk == 2 ) {
    err2     = mre2*pow( par2[0] + 2.*par2[1]*chi + 4.*par2[2]*pow(chi,3) + 5.*par2[3]*pow(chi,4),2);
    rpres[1] = sqrt( err2 );
  }

  
  LOG(INFO) << " Getting correction factor wCont ############" ;
  if( kk == 1 ) {
    LOG(INFO) << " k = " << kk << " : " 
	      << std::setw(6)  << c90 << "/" << c0 << " = " << mr << " chi = " << chi 
	      << " <cos(dphi)> "
	      << std::setw(14) << rpres[0] << " +- " << std::setw(10) << rpres[1]
	      << FairLogger::endl;
  }
  else if( kk == 2) {
    LOG(INFO) << " <cos(2dphi)> "
	      << std::setw(14) << rpres[2] << " +- " << std::setw(10) << rpres[3]
	      << FairLogger::endl; 
  }
  if( std::isnan( rpres[0] ) || std::isnan( rpres[1] ) ) {
    rpres[0] = 1.;
    rpres[1] = 0.;
    LOG(ERROR) << rpres[0] << FairLogger::endl;
  }

  return rpres;
} //Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180)
  

Double_t *GetRPResolutionwChi(Double_t *rpres, TH1D *hphi0_180, TH1D *hphi90_180, UInt_t kk) 
{
  if( rpres == NULL ) return rpres;

  hphi90_180->SetLineColor(2);

  Double_t m0 = hphi0_180->GetEntries();
  Double_t m1 = hphi90_180->GetEntries();

  return GetRPResolutionwCount(rpres, m0, m1, kk);

} //Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180)

Double_t *GetRPResolutionMCos(Double_t *rpres, Double_t mcos, Double_t mcose, UInt_t kk )            //%% Executable : 
{
  if( rpres == NULL )       return rpres;
  if( kk != 1 && kk != 2)   return rpres;


  LOG(INFO) << " GetRPReosolutionMCos()..... " << FairLogger::endl;

  rpres[0] = 1.;
  rpres[1] = 0.;

  Complex xx[5];
  Complex aa[6];

  ppar[kk-1][0] = -mcos;

  for(UInt_t i = 0; i < 6; i++ )
    aa[i] = tocomplex( ppar[kk-1][5-i], 0.);

  UInt_t ks = kk;
  dka(aa, xx, 5, 0.001, 100);

  std::vector<Double_t> xsolve;
  std::vector<Double_t> xerr;

  for(UInt_t i = 0; i < 5; i++ ) {
    if( abs(xx[i].i) < 0.01 && xx[i].r > 0 && xx[i].r < 3.) 
      {
	xsolve.push_back( xx[i].r );
	
	Double_t err = 0;
	for(UInt_t j = 1; j < 5; j++ )
	  err += pow(ppar[ks-1][j]*mcose,2);

	xerr.push_back(sqrt(err));
      }
  }


  if( xsolve.size() == 1 ) 
    {
      LOG(INFO) << xsolve.at(0) << " and " << xerr.at(0)  << FairLogger::endl;

      *rpres     = xsolve.at(0);
      *(rpres+1) = xerr.at(0);
    }

    else if( xsolve.size() > 1 ) 
      {
	LOG(ERROR) << " Resolution is not solved Please check it size: " << xsolve.size() << FairLogger::endl;

	for(UInt_t i = 0; i < (UInt_t)xsolve.size(); i++ ) {
	    LOG(INFO) << " chi_sub = " << xsolve.at(i)
		      << FairLogger::endl;    
	}
      }
  return rpres;
}

//====>>>
Double_t *GetRPResolutionwChi(Double_t *rpres, Double_t chi, Double_t chie, const UInt_t k)            //%% Executable : 
{
  if( rpres == NULL ) return rpres;

  LOG(INFO) << " GetRPResolutionwChi(double)..... " << FairLogger::endl;

  rpres[0] = 1.;
  rpres[1] = 0.;

  Float_t kk = (Float_t)k;


  cout << " wChi  k = " << k << endl; 

  if( k == 1 || k == 2 ) 
    {
      rpres[0] = sqrt(TMath::Pi())/(2.*sqrt(2))*chi*exp(-chi*chi/4.)*
	(ROOT::Math::cyl_bessel_i((kk-1.)/2.,chi*chi/4.)+ROOT::Math::cyl_bessel_i((kk+1.)/2.,chi*chi/4.));

      Double_t rrr = 0.; 
      for( UInt_t i = 1; i < 6; i++ ) 
       	rrr += ppar[k-1][i]*pow(chi,i);
      
      rpres[0] = rrr;

      Double_t err2 = 0.;
      for( Double_t i = 6; i > 1; i-- ) 
	err2 += pow( (i-1.) * ppar[k-1][(Int_t)i-1] * pow(chi, i-2), 2 );
      err2 = pow(chie,2)*err2;

      rpres[1] = sqrt( err2 );
    }

  LOG(INFO) << " Getting correction factor for k = " << k 
	    << std::setw(6)  << " chi = " << chi << " +- " << chie 
   	    << std::setw(14) << " res = " << rpres[0] << " +- " << std::setw(10) << rpres[1]
   	    << FairLogger::endl;

  if( std::isnan( rpres[0] ) || std::isnan( rpres[1] ) ) {
    rpres[0] = 1.;
    rpres[1] = 0.;
    LOG(INFO) << rpres[0] << FairLogger::endl;
  }

  return rpres;
} //Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180)


