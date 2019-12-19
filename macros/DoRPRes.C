#include "openRunAna.C"
#include "DoFlow.h"
#include "SimFunction.C"


auto *fcos1 = new TF1("fcos1","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());

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

//PSi_lab dependent correction factor
std::vector< Double_t > v1psix;
std::vector< Double_t > v1psiy;
std::vector< Double_t > v1psixe;
std::vector< Double_t > v1psiye;
std::vector< Double_t > v2psix;
std::vector< Double_t > v2psiy;
std::vector< Double_t > v2psixe;
std::vector< Double_t > v2psiye;


// functions
TString  GetPsiRPFileName();
Double_t GetRPBaseAngle(STFlowInfo *aflow);
Double_t GetError(Double_t x, Double_t y, Double_t xe, Double_t ye);
void     GetResolution();
void     CentralityDependence();            
void     PsiAngleDependence();      
UInt_t   GetPsiRPIndex(Double_t aVal);
UInt_t   GetRPCorrIndex(Double_t mult);
Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180);
//--


//-------------------//
void DoRPRes(Int_t isel = 0) 
{
  gROOT->Reset();

  openRunAna();

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
    PsiAngleDependence()   ;
  }

  else if( isel == 2 ) 
    GetResolution();
}

void GetResolution()
{
  // Averaged resolution is calculated.

  TH1I *hmult = new TH1I("hmult",";Multiplicity",100,0,100);
  TH1D *hphi0_180  = new TH1D("hphi0_180" , "0to180",100,0.,3.2);
  TH1D *hphi90_180 = new TH1D("hphi90_180","90to180",100,0.,3.2);

  Int_t nevt = SetBranch();  
  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  for(Int_t i = 0; i < nevt; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    if( aflow->mtrack4 < 3 ) continue;

    Double_t subevt_phi = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi()-
						   (aflow->unitP_2fc).Phi()));
      
    hphi0_180->Fill( subevt_phi );
    if( subevt_phi > TMath::Pi()/2. )
      hphi90_180->Fill( subevt_phi );

    hmult->Fill( aflow->mtrack2 );

  }

  auto reaolution = GetRPResolutionwChi(hphi0_180, hphi90_180);

  LOG(INFO) << " <cos (d phi) > = " << *reaolution << FairLogger::endl;


  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hmult->Draw();

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hphi0_180->Draw();
  hphi90_180->SetLineColor(2);
  hphi90_180->Draw("same");
  
}

void PsiAngleDependence()            //%% Executable :
{
  LOG(INFO) << " PsiAngleDependence .... " << FairLogger::endl;

  gROOT->Reset();
  gROOT->cd();

  TDatime beginTime;
  TDatime dtime;

  TString fName = GetPsiRPFileName();
  auto GraphSave = new TFile(fName,"recreate");

  TH1I *hmult  = new TH1I("hmult" ,"multiplicity",100,0,100);
  TH1I *hmult1 = new TH1I("hmult1","multiplicity",100,0,100);
  TH1I *hmult2 = new TH1I("hmult2","multiplicity",100,0,100);

  TH1I *hphibin[npsi];
  TH1D *hphi0_180[npsi];
  TH1D *hphi90_180[npsi];
  TH1D *h2phi0_180[npsi];
  TH1D *h2phi90_180[npsi];
  TH1D *hpsi[npsi];


  for(UInt_t k = 0; k < npsi; k++){
    TString htitle = Form("hphi0_180_%d",k);
    hphi0_180[k]  = new TH1D(htitle, "",100,0.,3.2);
    htitle = Form("hphi90_180_%d",k);
    hphi90_180[k] = new TH1D(htitle,"",100,0.,3.2);
    hpsi[k] = new TH1D(Form("hpsi%d",k),"",100,-3.15,3.15);

    htitle = Form("h2phi0_180_%d",k);
    h2phi0_180[k]  = new TH1D(htitle, "",100,0.,3.2);
    htitle = Form("h2phi90_180_%d",k);
    h2phi90_180[k] = new TH1D(htitle,"",100,0.,3.2);
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


    hmult ->Fill( aflow->mtrack4 );
    hmult1->Fill( aflow->mtrack1 );
    hmult2->Fill( aflow->mtrack2 );

    bRes = kTRUE; //@@@ 
    TIter next(aArray);
    STParticle *aPart = NULL;


    auto RPangle = GetRPBaseAngle(aflow);
    UInt_t iphi  = GetPsiRPIndex(RPangle);

    hpsi[iphi] ->Fill( TVector2::Phi_mpi_pi(RPangle) );
    hiphi->Fill( iphi );

    Double_t subdphi = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi()));
    hphi0_180[iphi] ->Fill(subdphi);

    if( subdphi > TMath::Pi()/2. )
      hphi90_180[iphi]->Fill(subdphi);

    subdphi = abs(TVector2::Phi_mpi_pi(2.*((aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi())));
    h2phi0_180[iphi] ->Fill(subdphi);

    if( subdphi > TMath::Pi()/2. )
      h2phi90_180[iphi]->Fill(subdphi);

  }

  auto gv_psi1 = new TGraphErrors();
  gv_psi1->SetName("gv_psi1");
  gv_psi1->SetTitle("; #psi; <cos(#Delta #Psi)>");

  auto gv_psi2 = new TGraphErrors();
  gv_psi2->SetName("gv_psi2");
  gv_psi2->SetTitle("; #psi; <cos(#Delta 2#Psi)>");
   

  UInt_t id = 1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1000);
  cc->Divide(2,npsi);

  UInt_t ip = 0;
  Double_t *rpres = new Double_t[4];

  UInt_t nst[] = {6, 0};
  for(UInt_t j = 0; j < 2; j++) {
    for(UInt_t i = nst[j]; i < nst[j]+6; i++) {
      cc->cd(id); id++;
      hphi0_180[i]->Draw();
      hphi90_180[i]->SetLineColor(2);
      hphi90_180[i]->Draw("same");

      cc->cd(id); id++;
      hpsi[i]->Draw();

      if( hphi0_180[i]->GetEntries() < 5 ) continue;

      rpres = GetRPResolutionwChi(hphi0_180[i], hphi90_180[i]);
    
      Double_t psi = hpsi[i]->GetMean();
      gv_psi1->SetPoint(ip, psi, rpres[0]);
      gv_psi1->SetPointError(ip, 0., rpres[1]);
      

      //      rpres = GetRPResolutionwChi(h2phi0_180[i], h2phi90_180[i]);

      gv_psi2->SetPoint(ip, psi, rpres[2]);
      gv_psi2->SetPointError(ip, 0., rpres[3]);
      ip++;

    }
  }

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));  
  gv_psi1->SetMarkerStyle(2);
  gv_psi1->Draw("ALP");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));  
  gv_psi2->SetMarkerStyle(2);
  gv_psi2->Draw("ALP");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hiphi->Draw();


  hmult->Write();
  hmult1->Write();
  hmult2->Write();
  hiphi->Write();
  gv_psi1->Write();
  gv_psi2->Write();

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

    subdphi = abs(TVector2::Phi_mpi_pi(2.*((aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi())));

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

    rpres = GetRPResolutionwChi(hphi0_180[i], hphi90_180[i]);
    LOG(INFO) << " Multiplicity " << hmultbin4[i]->GetMean() << " > " << mrange[i]
	      << " <cos(Phi)> = " << rpres[0] << " +- " << rpres[1] 
	      << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
	      << FairLogger::endl;

    // rpres = GetRPResolutionwChi(h2phi0_180[i], h2phi90_180[i]);
    // LOG(INFO) << " Multiplicity " << hmultbin4[i]->GetMean() << " > " << mrange[i]
    // 	      << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
    // 	      << FairLogger::endl;
     
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
  
  
  // if( abs( aVal ) > TMath::Pi() ) return (UInt_t)v1psix.size();

  // for(auto i = 1; i < (UInt_t)v1psix.size(); i++ ) {

  //   if( aVal < v1psix.at(i) ) {
  //     return i - 1;
  //   }

  // }

  // return (UInt_t)v1psix.size() - 1;
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
TString GetPsiRPFileName()
{
  //  TString fn = "data/bpsi_"+ sysName + ".v" + oVer;
  TString fn = "data/cpsi_"+ sysName + ".v" + oVer;
  //  TString fn = "data/dpsi_"+ sysName + ".v" + oVer;
  //  TString fn = "data/epsi_"+ sysName + ".v" + oVer;
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
Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180)            //%% Executable : 
{
  Double_t *rpres = new Double_t[4];
  hphi90_180->SetLineColor(2);

  Double_t m0 = hphi0_180->GetEntries();
  Double_t m1 = hphi90_180->GetEntries();

  if( m0 == 0 ) {
    LOG(INFO) << " No entry in hphi0_180 " << endl;
    rpres[0] = 1.;
    rpres[1] = 0.;
    rpres[2] = 1.;
    rpres[3] = 0.;
    return rpres;
  }

  Double_t mr    = m1/m0;
  Double_t mre2  = 1./m0;

  Double_t chi   = sqrt(-4.* log(2.* mr));
  Double_t dchin = sqrt(-4.* log(2.*(mr-1./sqrt(m1))));
  Double_t chie2 = pow( chi - dchin, 2 );
  Double_t chie = sqrt(chie2);

  Double_t par1[] = {0.626657, -0.09694, 0.02754, -0.002283};
  Double_t par2[] = {0.25,  -0.011414, -0.034726, 0.006815};

  //v1
  if( kFALSE ) {
    // //v1
    rpres[0] = par1[0]*chi + par1[1]*pow(chi,3) + par1[2]*pow(chi,4) + par1[3]*pow(chi,5);  
    // //v2
    rpres[2] = par2[0]*pow(chi,2) + par2[1]*pow(chi,3) + par2[2]*pow(chi,4) + par2[3]*pow(chi,5);
  }
  else {
    rpres[0] = sqrt(TMath::Pi())/(2.*sqrt(2))*chi*exp(-chi*chi/4.)*
      (ROOT::Math::cyl_bessel_i(0,chi*chi/4.)+ROOT::Math::cyl_bessel_i(1,chi*chi/4.));
    rpres[2] = sqrt(TMath::Pi())/(2.*sqrt(2))*chi*exp(-chi*chi/4.)*
      (ROOT::Math::cyl_bessel_i(0.5,chi*chi/4.)+ROOT::Math::cyl_bessel_i(1.5,chi*chi/4.));
  }

  Double_t err2 = mre2*pow(par1[0] + 3.*par1[1]*pow(chi,2) + 4.*par1[2]*pow(chi,3) + 5.*par1[3]*pow(chi,4),2);
  rpres[1] = sqrt( err2 );
  err2     = mre2*pow( par2[0] + 2.*par2[1]*chi + 4.*par2[2]*pow(chi,3) + 5.*par2[3]*pow(chi,4),2);
  rpres[3] = sqrt( err2 );
 
  LOG(INFO) << " Getting correction factor " 
	    << std::setw(6)  << m1 << "/" << m0 << " = " << mr << " chi = " << chi 
	    << " v1 "
   	    << std::setw(14) << rpres[0] << " +- " << std::setw(10) << rpres[1]
   	    << " v2 "
   	    << std::setw(14) << rpres[2] << " +- " << std::setw(10) << rpres[3]
   	    << FairLogger::endl;

  if( std::isnan( rpres[0] ) || std::isnan( rpres[1] ) ) {
    rpres[0] = 1.;
    rpres[1] = 0.;
    rpres[2] = 1.;
    rpres[3] = 0.;
    LOG(INFO) << rpres[0] << FairLogger::endl;
  }


  return rpres;
} //Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180)


