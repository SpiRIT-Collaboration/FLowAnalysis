#include "FlowFunctions.h"
#include "openFlw.C"
#include "GetLorentzBoost.C"


auto *fcos1 = new TF1("fcos1","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());

Int_t  seltrackID = 4;
UInt_t selReactionPlanef = 10;
Int_t  seltrack;

//drawing
TCanvas *cc[20];
const UInt_t nsys = 4;
const UInt_t nprt = 5;
TString  iopt[]     = {"","same","same","same","same", "same"};
UInt_t   imark[]    = {20, 21, 22, 23, 25};
  
UInt_t ic = -1;
const Int_t nbinx = 30;

// rapidity
Double_t rapid_max = 0.5;
Double_t rapid_min = 0.2;
const Double_t ycm[]      = {0.382453,  0.364873, 0.390302, 0.354066, 0.371326};
const Double_t ybeam_cm[] = {0.360199,  0.377779, 0.354065, 0.390301, 0.371326};
TLorentzVector* bmVec;
TLorentzVector* tgVec;

const UInt_t nspec = 5;
const UInt_t  nbin = 9;

//UInt_y   y_nbin  = 100;
Double_t y_min[] = {0., 0., 0., 0., 0.};
Double_t y_max[] = {1.2, 1.2, 1.2, 1.2, 1.2};
Double_t y_bin[] = {0.05, 0.05, 0.01, 0.01, 0.01};
Double_t y_binx  =  0.1;

Double_t yrange1[] = { -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.5};
const UInt_t ybin1 = sizeof(yrange1)/sizeof(Double_t);

Double_t yrange2[] = {-0.2, -0.05,  0.05, 0.2, 0.35, 0.5};
const UInt_t ybin2 = sizeof(yrange2)/sizeof(Double_t);

// TCutG* ganglecut = new TCutG();
// ganglecut->SetPoint(0,   0., 3.2);
// ganglecut->SetPoint(1, 1.43, 3.2);
// ganglecut->SetPoint(2, 0.94, 2.8);
// ganglecut->SetPoint(3, 1.26,2.44);

TString  partname[] = {"pi-","pi+","proton","deuteron" ,"triton", "3He", "4He", "neutron"};
UInt_t   partid[]   = {211,    211,    2212, 1000010020, 1000010030, 1000020030, 1000020040, 2112};

Double_t cutbmass[] = {191.1, 191.1, 1165.9, 2249.9, 3475.2};
Double_t cutlmass[] = {0.,      0. ,    0.,     0. , 2500.};


// Retrivew tree
Int_t ntrack[7];
TClonesArray *aArray; // = new TClonesArray("STParticle",100);
TClonesArray *aNLClusterArray;// = new TClonesArray("STNueLANDCluster",100);
Double_t aoq;
Double_t z;
Int_t    snbm;
Double_t ProjA;
Double_t ProjB;
Int_t    mtrack_1;
Int_t    mtrack_2;
TVector3 *unitP_fc    = NULL;
TVector3 *unitP_rc    = NULL;
TVector2 *unitP_1     = NULL;
TVector2 *unitP_2     = NULL;
TVector2 *unitP_lang  = NULL;
TBranch  *bunitP_fc;
TBranch  *bunitP_rc;
TBranch  *bunitP_1;
TBranch  *bunitP_2;
TBranch  *brpphi   = 0;
TBranch  *biphi    = 0;
TBranch  *bdeltphi = 0;
TBranch  *bunitP_lang = 0;
TVector2 *eisott;
TVector2 *eisotg;
TVector2 *eisobm;


Int_t   iVer;
TString sVer;

// functions
void     GetFittingParameters(TH1D &h1, Double_t pp[6]);
void     GetFittingParameters(TH1D &h1, Double_t pp[6], Double_t corr[2]);
Double_t GetRapidity_cm(TVector3 p, Double_t mass, TVector3 bvec);
UInt_t   SetBranch(UInt_t m=0);

void     PlotPtDependence(UInt_t selid);     
Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180);
void     PlotNeuLANDv1v2();
void     PlotNeuLANDCosv1v2();
void     PlotNeuLANDProperty(UInt_t iout=0);
void     PlotCosPtDependence(UInt_t selid);

Double_t  GetError(Double_t x, Double_t y, Double_t xe, Double_t ye)
{
  if( y == 0 && ye == 0 ) return 0.;
  Double_t xc = abs(x/y);
  return  xc * sqrt(pow(xe/x,2) + pow(ye/y,2));
}

//Double_t rapoffset[] = {0., 0.018, -0.008, 0.028};
Double_t rapoffset[] = {0.01113, -0.00645, 0.01898, -0.01726};
//UInt_t   cent[]     = {70, 40, 35, 20, 0};    //12.6, 20.9, 61.6, 14.5 %
UInt_t   cent[]     = {70, 35, 28, 0};        //30.7, 33.9, 43.5 %
Double_t isobin[]   = {1., 0.8, 0.7, 0.5, 0};
UInt_t   npb = 30;

UInt_t  icent = 0;
UInt_t  jcent = 4;
UInt_t iaz = 0;

//--------- main ----------//
void calcFlw(UInt_t isel = 0) 
{
  gROOT->Reset();


  openFlw();

  for(UInt_t i = 0; i < 4; i++){
    if(rChain[i] != NULL) {    
      std::cout << " System " << i << " "  << isys[i] << "  -> " << sysName[isys[i]] << std::endl; 
    }
  }

  if(rChain[0] == NULL)
    exit(0);

  gROOT->ProcessLine(".! grep -i void calcFlw.C | grep '//%%'");


  // Multiplicity cut ================================
  TString su = gSystem -> Getenv("UC");
  if( su != "" ) {
    UInt_t ucent = (UInt_t)atoi(su);
    if( ucent < sizeof(cent)/sizeof(UInt_t)-1)
      icent = ucent;
  }


  su = gSystem -> Getenv("LC");
  if( su != "" ) {
    UInt_t lcent = (UInt_t)atoi(su);
    if(  lcent < sizeof(cent)/sizeof(UInt_t))
      jcent = lcent;
  }

  std::cout << " Multiplicity :: " << cent[icent] << " to " << cent[jcent] << std::endl;

  // azimuthal cut ===================================
  TString az = gSystem -> Getenv("AZ");
  if( az != "" ) 
    iaz = (UInt_t)atoi(az);

  std::cout << " Azimuthal angle cut " << iaz << std::endl;
  //==================================================


  //  PlotNeuLANDProperty(1);
  //  PlotPtDependence(isel);

  if( isel > 0 && isel < 5)
    PlotCosPtDependence(isel);  
  else if ( isel == 5 )
    PlotNeuLANDCosv1v2(); 
  else if ( isel > 5 && isel <=  13)
    PlotCosPtDependence(isel);  
    

  //  PlotNeuLANDv1v2();
}
//**************************************************
//**************************************************
//pdt
void PlotCosPtDependence(UInt_t selid = 2)       //%% Executable :
{

  std::cout << "PlotCosPtDependence(" << selid << ")" << std::endl;

  gStyle->SetOptStat(0);

  auto cutfile = new TFile("data/RegionCut.root");
  TCutG *goodThetaPhi = (TCutG*)cutfile->Get("goodThetaPhi");
  cutfile->Close();    

  gSystem->cd("data");


  TString fName = Form("cosYPT_n%2dto%d_",cent[icent],cent[jcent]) + sysName[isys[0]] + "yoffp_" + Form("az%d",iaz) + partname[selid];
  //  TString fName = Form("cosYPT_n%2dto%d_",cent[icent],cent[jcent]) + sysName[isys[0]] + "_" + partname[selid];
  //  TString fName = Form("cosYPT_iso%2dto%d_",(UInt_t)isobin[icent],(UInt_t)isobin[jcent]) + sysName[isys[0]] + "_" + partname[selid]+".root";


  auto GraphSave = new TFile(fName+".root","recreate");
  std::cout << "File " << GraphSave->GetName() << " is created. " << std::endl;

  auto hphi0_180  = new TH1D("hphi0_180" ,"#Phi from  0 to 180; #Phi",100,0.,3.2);
  auto hphi90_180 = new TH1D("hphi90_180","#Phi from 90 to 180; #Phi",100,0.,3.2);

  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  Int_t RPflag = 0;
  if(selid >= 2) RPflag = 10;

  cout << " Particle " << partname[selid] << endl;
  
  // PT binning
  Double_t pt_max = 800.;
  UInt_t   nbin1  = 10; //16
  UInt_t   nbin2  = 8; //10
  Double_t dpt1   = pt_max/(Double_t)nbin1;
  Double_t dpt2   = pt_max/(Double_t)nbin2;

  std::cout << " Rapidity binning " << ybin1 << std::endl;
  TString rangeLabel1[ybin1];
  for(UInt_t i = 0; i < ybin1; i++ ){
    if( i == 0 )
      rangeLabel1[0] = Form(" y < %5.2f ",yrange1[0]);
    else if ( i == ybin1 -1 )
      rangeLabel1[i] = Form("%5.2f <= y "    ,yrange1[ybin1-1]);
    else 
      rangeLabel1[i] = Form("%5.2f <= y < %5.2f",yrange1[i-1],yrange1[i]);
  }
  TString rangeLabel2[ybin2];
  for(UInt_t i = 0; i < ybin2; i++ ){
    if( i == 0 )
      rangeLabel2[0] = Form(" y < %5.2f ",yrange2[0]);
    else if ( i == ybin1 -1 )
      rangeLabel2[i] = Form("%5.2f <= y "    ,yrange2[ybin2-1]);
    else 
      rangeLabel2[i] = Form("%5.2f <= y < %5.2f",yrange2[i-1],yrange2[i]);
  }

  // booking
  auto hmass   = new TH2D("hmass",   ";P/Q; Mass [MeV/c^2]"     ,200,  0.,2500., 200, 0.,7000);
  auto hmassy  = new TH2D("hmassy",  "; Y_{cm}; Mass [MeV/C^2]" ,200, -0.4, 0.5, 200, 0.,7000);
  auto hpy     = new TH2D("hpy",     "; Y_{cm}; P/Q [MeV/c]"    ,200, -0.4, 0.5, 200, 0.,2500.);
  auto hptmass = new TH2D("hptmass", "; Pt [MeV/c]; Mass [MeV/c^2]",200, 0., 800,200, 0.,7000);
  auto hazm    = new TH1D("hazm",    "; #phi"                   ,100, -3.2, 3.2);
  auto hmidv2  = new TH1D("hmidv2",  "; #Delta#phi"             ,100, -3.2, 3.2);

  //  TString hlabel = (TString)Form("iso %2d to %d ; Y_{cm}; Pt [MeV/c]",(UInt_t)isobin[icent],(UInt_t)isobin[jcent]);
  TString hlabel = (TString)Form("ntrack[4] %2d to %d ; Y_{cm}; Pt [MeV/c]",(UInt_t)cent[icent],(UInt_t)cent[jcent]);
  auto hyptacp = new TH2D("hyptacp", hlabel ,200, -0.4, 0.5, 200, 0., 800);
  auto heisoratio = new TH1D("heisoratio",";EISO_beam Y",100,0.,1.);
  auto hntrack   = new TH1I("hntrack",  "; Number of good track", 60, 0, 60);


  TH1D *hyphi1[ybin1];
  TH1D *hyphi2[ybin2];
  TH1D *hyptphi1[ybin1][nbin1];
  TH1D *hyptphi2[ybin2][nbin2];

  std::vector< std::vector< Double_t > > cosv1x(ybin1);
  std::vector< std::vector< Double_t > > cosv2x(ybin2);
  Double_t  cosv1[ybin1];
  Double_t  cosv2[ybin2];
  Double_t  sinv1[ybin1];
  Double_t  sinv2[ybin2];

  std::vector< std::vector< std::vector< Double_t> > > cosv1ptx(ybin1);
  Double_t  cosv1pt[ybin1][nbin1];
  Double_t  sinv1pt[ybin1][nbin1];
  for( UInt_t i = 0; i < ybin1; i++) {
    cosv1[i]   = 0.;
    sinv1[i]   = 0.;
    cosv1x[i].clear();
    cosv1ptx[i].resize(nbin1);

    hyphi1[i] = new TH1D( Form("hyphi1_%d",i),rangeLabel1[i]+"#Delta #phi"   , npb, -3.15, 3.15);

    for( UInt_t j = 0; j < nbin1; j++ ){
      cosv1pt[i][j] = 0.;
      sinv1pt[i][j] = 0.;

      hyptphi1[i][j] = new TH1D( Form("hyptphi1_%d%d",i,j),rangeLabel1[i]+"#Delta #phi"   , npb, -3.15, 3.15);
    }
  }

  std::vector< std::vector< std::vector< Double_t> > > cosv2ptx(ybin2);
  Double_t  cosv2pt[ybin2][nbin2];
  Double_t  sinv2pt[ybin2][nbin2];
  for( UInt_t i = 0; i < ybin2; i++) {
    cosv2[i]   = 0.;
    sinv2[i]   = 0.;
    cosv2x[i].clear();
    cosv2ptx[i].resize(nbin2);
    hyphi2[i] = new TH1D( Form("hyphi2_%d",i),rangeLabel2[i]+"2x #Delta #phi", npb, -3.15, 3.15);

    for( UInt_t j = 0; j < nbin; j++ ){
      cosv2pt[i][j] = 0.;
      sinv2pt[i][j] = 0.;
      hyptphi2[i][j] = new TH1D( Form("hyptphi2_%d%d",i,j),rangeLabel2[i]+"2x #Delta #phi"   , npb, -3.15, 3.15);
    }
  }
    
    
  // Lorentz Transform    
  TVector3 boostVec = LorentzBoost(4);
  
  Int_t nevt = SetBranch();
  cout << " Number of events " << nevt << endl;
  
  for(Int_t i = 0; i < nevt; i++){
    rChain[0]->GetEntry(i);
      
    /// good event selection
    if(ntrack[2] == 0) continue;

    /// centrality selection
    if(ntrack[4] > cent[icent] || ntrack[4] <= cent[jcent] ) continue;
    Double_t isoratio = (eisobm->Y()+eisotg->Y())/(eisobm->Mod()+eisotg->Mod());
    //if( isoratio > isobin[icent] || isoratio <= isobin[jcent] ) continue;

    heisoratio->Fill( isoratio );
    hntrack->Fill( ntrack[4] );


    //    if( unitP_fc->Phi() < 0. || unitP_fc->Phi() > 50.*TMath::Pi()/180.) continue; 
    //    if( abs(unitP_fc->Phi()) < 130.*TMath::Pi()/180. ) continue;

    Double_t subevt_phi = abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi())); 
    hphi0_180->Fill( subevt_phi );
    if( subevt_phi > TMath::Pi()/2. )
      hphi90_180->Fill( subevt_phi );


    TIter next(aArray);
    STParticle *aPart = NULL;

    while( (aPart = (STParticle*)next()) ) {
	
      auto rpf   = aPart->GetReactionPlaneFlag();
      auto pid   = aPart->GetPID();
      auto bmass = aPart->GetBBMass();
      auto phi   = aPart->GetRotatedMomentum().Phi();
      auto theta = aPart->GetRotatedMomentum().Theta();
      
      if( pid == partid[selid] ) { //&& 

	if( iaz == 13 ){
	  if( !goodThetaPhi->IsInside(theta, phi) ) continue;
	}
	else if( iaz == 12){
	  if( abs(phi) < 135.*TMath::DegToRad() &&  abs(phi) >  45.*TMath::DegToRad() ) continue;
	}
	else if( iaz == 1 || iaz == 10 || iaz == 11 ) {
	  if( abs(phi) >  45.*TMath::DegToRad() ) continue;
	}
	else if( iaz== 2 || iaz == 20 || iaz == 21 ) {
	  if( abs(phi) < 135.*TMath::DegToRad() ) continue;
	}
	else if( iaz == 3 || iaz == 30 || iaz == 31 ) {
	  if( abs(phi) > 135.*TMath::DegToRad() ||  abs(phi) <  45.*TMath::DegToRad() ) continue;
	}

	if( iaz == 10 || iaz == 20 || iaz == 30 ){
	  if( phi > 0 ) continue;
	}
	else if( iaz == 11 || iaz == 21 || iaz == 31 ) {
	  if( phi < 0 ) continue;
	}

	hazm -> Fill(phi);



	auto pt    = aPart->GetRotatedMomentum().Pt();
	auto dphi  = aPart->GetAzmAngle_wrt_RP();
	auto rapid = aPart->GetRapiditycm();;
	
	
	rapid -= rapoffset[isys[0]];

	hmass->Fill( aPart->GetRotatedMomentum().Mag(), bmass);
	
	hyptacp->Fill( rapid, pt );
	hmassy ->Fill( rapid, bmass );
	hpy    ->Fill( rapid, aPart->GetRotatedMomentum().Mag() );
	hptmass->Fill( pt,   bmass );


	UInt_t irapid = ybin1 - 1;
	for( UInt_t i = 0; i < ybin1; i++){
	  if(rapid < yrange1[i]){
	    irapid = i;
	    break;
	  }
	}
	
	hyphi1[irapid]->Fill(dphi);
	cosv1x[irapid].push_back( rapid );
	cosv1[irapid] += cos(dphi);
	sinv1[irapid] += sin(dphi);

	UInt_t ipt = nbin1 - 1;
	for(UInt_t i = 0; i < nbin1; i++){
	  if( pt < dpt1*(i+1)) {
	    ipt = i;
	    break;
	  }
	}

	hyptphi1[irapid][ipt]->Fill(dphi);
	cosv1ptx[irapid][ipt].push_back( pt );
	cosv1pt[irapid][ipt] += cos(dphi);
	sinv1pt[irapid][ipt] += sin(dphi);


	irapid = ybin2 - 1;

	for( UInt_t i = 0; i < ybin2; i++){
	  if(rapid < yrange2[i]){
	    irapid = i;
	    break;
	  }
	}

	hyphi2[irapid]->Fill(TVector2::Phi_mpi_pi(2.*dphi));
	cosv2x[irapid].push_back( rapid );
	cosv2[irapid] += cos(2.*dphi);
	sinv2[irapid] += sin(2.*dphi);

	ipt = nbin2 - 1;
	for(UInt_t i = 0; i < nbin2; i++){
	  if( pt < dpt2*(i+1)) {
	    ipt = i;
	    break;
	  }
	}

	hyptphi2[irapid][ipt]->Fill(TVector2::Phi_mpi_pi(2.*dphi));
	cosv2ptx[irapid][ipt].push_back( pt );
	cosv2pt[irapid][ipt] += cos(2.*dphi);
	sinv2pt[irapid][ipt] += sin(2.*dphi);


      }
    }
  }
      
    
  TGraphErrors *gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  TGraphErrors *gv_v2 = new TGraphErrors();
  gv_v2->SetName("gv_v2");
  
  TGraphErrors *gPt_v1[ybin1];
  TGraphErrors *gPt_v2[ybin1];
    
  for(UInt_t kn = 0; kn < ybin1 ; kn++){      
    gPt_v1[kn] = new TGraphErrors();
    gPt_v1[kn]->SetName((TString)Form("gPt_v1%d",kn));
    TString sname = partname[selid]+"; Pt [MeV/c]; v1";
    gPt_v1[kn]->SetTitle(sname);
  }
  
  for(UInt_t kn = 0; kn < ybin2 ; kn++){      
    gPt_v2[kn] = new TGraphErrors();
    gPt_v2[kn]->SetName((TString)Form("gPt_v2%d",kn));
    TString sname = partname[selid]+"; Pt [MeV/c]; v2";
    gPt_v2[kn]->SetTitle(sname);
  }


  Double_t *rpres = new Double_t[4];
  rpres = GetRPResolutionwChi(hphi0_180, hphi90_180);
  std::cout << " <cos(Phi)> = " << rpres[0] << " +- " << rpres[1] 
	    << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
	    << std::endl;


  std::cout << " ---- Resutls ---------------------" << std::endl;
  
  UInt_t kl = 0;
  UInt_t id1 = 0;
  UInt_t id2 = 0;
  Double_t para[6];

    
  for(UInt_t kn = 0; kn < ybin1; kn++){
      
    if( cosv1x[kn].size() == 0 ) continue;

    //v1 rapidity dependence
    Double_t rapm  = TMath::Mean(cosv1x[kn].begin(), cosv1x[kn].end());
    Double_t rape  = TMath::StdDev(cosv1x[kn].begin(), cosv1x[kn].end());
      
    Double_t yv1  = cosv1[kn] / (Double_t)cosv1x[kn].size();
    Double_t yv1e = sinv1[kn] / (Double_t)cosv1x[kn].size();

    Double_t yv1c  = yv1/rpres[0];
    Double_t yv1ce = GetError(yv1, rpres[0], yv1e, rpres[1]);

    gv_v1->SetPoint( kl,      rapm, yv1c);
    gv_v1->SetPointError( kl, rape, yv1ce);
    kl++;

    // pt dependence 
    UInt_t il = 0; 
    for(UInt_t jn = 0; jn < nbin1; jn++){

      if( cosv1ptx[kn][jn].size() > 0 ) {	

	Double_t ptc  = TMath::Mean(cosv1ptx[kn][jn].begin(), cosv1ptx[kn][jn].end());
	Double_t ptce = TMath::StdDev(cosv1ptx[kn][jn].begin(), cosv1ptx[kn][jn].end());

	Double_t ptv1  = cosv1pt[kn][jn] / (Double_t)cosv1ptx[kn][jn].size();
	Double_t ptv1e = sinv1pt[kn][jn] / (Double_t)cosv1ptx[kn][jn].size();

	yv1c  = ptv1/rpres[0];
	yv1ce = GetError(ptv1, rpres[0], ptv1e, rpres[1]);

	if( kn == 0 ){
	  cout << jn << " erro " << sinv1pt[kn][jn] << " yvc1  " << yv1ce <<  endl;
	  cout << sinv1pt[kn][jn] << " --  " << (Double_t)cosv1ptx[kn][jn].size() << endl;
	}


	gPt_v1[kn]->SetPoint(il, ptc, yv1c);
	gPt_v1[kn]->SetPointError(il, ptce, yv1ce);
	    
	il++;
      }
    }
    gPt_v1[kn]->Write();
  }

  kl = 0;
  for(UInt_t kn = 0; kn < ybin2; kn++){

    if( cosv2x[kn].size() == 0 ) continue;

    Double_t rapm  = TMath::Mean(cosv2x[kn].begin(), cosv2x[kn].end());
    Double_t rape  = TMath::StdDev(cosv2x[kn].begin(), cosv2x[kn].end());

    Double_t yv2  = cosv2[kn] / (Double_t)cosv2x[kn].size();
    Double_t yv2e = sinv2[kn] / (Double_t)cosv2x[kn].size();

    Double_t yv2c  = yv2 / rpres[2];
    Double_t yv2ce = GetError(yv2, rpres[2], yv2e, rpres[3]);


    gv_v2->SetPoint( kl,    rapm, yv2c);
    gv_v2->SetPointError( kl, rape, yv2ce);      
    kl++;

    UInt_t il = 0; 
    for(UInt_t jn = 0; jn < nbin2; jn++){
	
      if( cosv2ptx[kn][jn].size() > 0 ){

	Double_t ptc  = TMath::Mean(cosv2ptx[kn][jn].begin(), cosv2ptx[kn][jn].end());
        Double_t ptce = TMath::StdDev(cosv2ptx[kn][jn].begin(), cosv2ptx[kn][jn].end());

	Double_t ptv2  = cosv2pt[kn][jn] / (Double_t)cosv2ptx[kn][jn].size();
        Double_t ptv2e = sinv2pt[kn][jn] / (Double_t)cosv2ptx[kn][jn].size();

	yv2c  = ptv2/ rpres[2];
	yv2ce = GetError(ptv2, rpres[2], ptv2e, rpres[3]);

	gPt_v2[kn]->SetPoint(il, ptc, yv2c);
	gPt_v2[kn]->SetPointError(il, ptce, yv2ce);
	  
	il++;
      }
    }
    gPt_v2[kn]->Write();
  }

  gv_v1->Write();
  gv_v2->Write();


  //====================
  //plotting
  ic++; id=1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  cc[ic]->Divide(2,2);

  cc[ic]->cd(id); id++;
  hmass->Draw("colz");

  cc[ic]->cd(id); id++;
  hntrack->Draw();

  cc[ic]->cd(id); id++;
  hpy->Draw("colz");

  // cc[ic]->cd(id); id++;
  // hptmass->Draw("colz");

  cc[ic]->cd(id); id++;
  hazm->Draw();
  //  heisoratio->Draw("colz");  
  

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hyptacp->Draw("colz");
  


  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  cc[ic]->Divide(ybin1,1);
  for(UInt_t kn = 0; kn < ybin1; kn++){
    cc[ic]->cd(id); id++;
    gPt_v1[kn]->Draw("ALP");
  }

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1600,500);
  cc[ic]->Divide(ybin2,1);
  for(UInt_t kn = 0; kn < ybin2; kn++){
    cc[ic]->cd(id); id++;
    gPt_v2[kn]->Draw("ALP");
  }


  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v1->Draw("ALP");

  auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v2->Draw("ALP");

  ic++; 
  cc[ic]   = new TCanvas("dyphi","dphi1 y bin",500,1200);
  cc[ic]->Divide(2, ybin1);

  cc[ic+1]   = new TCanvas("dphi1","dphi1 y bin and pt bin",2000,1200);
  cc[ic+1]->Divide(nbin1, ybin1);

  cc[ic+2] = new TCanvas("dphi2","dphi2 ybin and pt bin",1400,1200);
  cc[ic+2]->Divide(nbin2, ybin2);

  id = 1;
  for(UInt_t i = 0; i < ybin1; i++) {
    cc[ic]->cd(2*(i+1)-1);
    hyphi1[i]->SetNormFactor(npb);
    hyphi1[i]->Draw("e");

    for(UInt_t j = 0; j < nbin1; j++){
      cc[ic+1]->cd(id);    id++;

      
      if(hyptphi1[i][j]->GetEntries() == 0) continue;
      
      hyptphi1[i][j]->SetNormFactor(npb);
      hyptphi1[i][j]->Draw("e");
    }
  } 

  id = 1;
  for(UInt_t i = 0; i < ybin2; i++) {
    cc[ic]->cd(2*(i+1));
    hyphi2[i]->SetNormFactor(60);
    hyphi2[i]->Draw("e");

    for(UInt_t j = 0; j < nbin2; j++){
      cc[ic+2]->cd(id); id++;
      
      if(hyptphi2[i][j]->GetEntries() == 0) continue;
      
      hyptphi2[i][j]->SetNormFactor(npb);
      hyptphi2[i][j]->Draw("e");

    }
  } 
  ic += 3;

  hphi0_180->Write();
  hyptacp->Write();
  heisoratio->Write();
  hntrack->Write();

  gSystem->cd("..");

  SaveCanvas(fName);

  //  std::cout << " TOP " << hisoratio->GetEntries()/


  //  ic++; 
  //  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

  //  gROOT->cd();
}

//neutron
void PlotNeuLANDCosv1v2()               //%% Executable :
{

  std::cout << "-----> PlotNeuLANDCosv1v2(); ----->" << std::endl;

  gStyle->SetOptStat(0);

  gSystem->cd("data");

  TString fName = Form("cosYPT_n%2dto%d_",cent[icent],cent[jcent]) + sysName[isys[0]] + "yoffp_neutron";
  //  TString fName = Form("cosYPT_n%2dto%d_",cent[icent],cent[jcent]) + sysName[isys[0]] + "_neutron";
  auto GraphSave = new TFile(fName+".root","recreate");
  std::cout << "File " << GraphSave->GetName() << " is created. " << std::endl;

  auto hphi0_180  = new TH1D("hphi0_180" ,"#Phi from  0 to 180; #Phi",100,0.,3.2);
  auto hphi90_180 = new TH1D("hphi90_180","#Phi from 90 to 180; #Phi",100,0.,3.2);

  // PT binning
  Double_t pt_max = 800.;
  UInt_t   nbin1  = 10; //16
  UInt_t   nbin2  = 8; //10
  Double_t dpt1   = pt_max/(Double_t)nbin1;
  Double_t dpt2   = pt_max/(Double_t)nbin2;

  std::cout << " Rapidity binning " << ybin1 << std::endl;
  TString rangeLabel1[ybin1];
  for(UInt_t i = 0; i < ybin1; i++ ){
    if( i == 0 )
      rangeLabel1[0] = Form(" y < %5.2f ",yrange1[0]);
    else if ( i == ybin1 -1 )
      rangeLabel1[i] = Form("%5.2f <= y "    ,yrange1[ybin1-1]);
    else 
      rangeLabel1[i] = Form("%5.2f <= y < %5.2f",yrange1[i-1],yrange1[i]);
  }
  TString rangeLabel2[ybin2];
  for(UInt_t i = 0; i < ybin2; i++ ){
    if( i == 0 )
      rangeLabel2[0] = Form(" y < %5.2f ",yrange2[0]);
    else if ( i == ybin1 -1 )
      rangeLabel2[i] = Form("%5.2f <= y "    ,yrange2[ybin2-1]);
    else 
      rangeLabel2[i] = Form("%5.2f <= y < %5.2f",yrange2[i-1],yrange2[i]);
  }

  // booking
  TString hlabel = (TString)Form("ntrack[4] %2d to %d ; Y_{cm}; Pt [MeV/c]",(UInt_t)cent[icent],(UInt_t)cent[jcent]);
  auto hyptacp = new TH2D("hyptacp", hlabel ,200, -0.4, 0.5, 200, 0., 800);
  auto heisoratio = new TH1D("heisoratio",";EISO_beam Y",100,0.,1.);
  auto hntrack   = new TH1I("hntrack",  "; Number of good track", 60, 0, 60);


  TH1D *hyphi1[ybin1];
  TH1D *hyphi2[ybin2];
  TH1D *hyptphi1[ybin1][nbin1];
  TH1D *hyptphi2[ybin2][nbin2];

  std::vector< std::vector< Double_t > > cosv1x(ybin1);
  std::vector< std::vector< Double_t > > cosv2x(ybin2);
  Double_t  cosv1[ybin1];
  Double_t  cosv2[ybin2];
  Double_t  sinv1[ybin1];
  Double_t  sinv2[ybin2];

  std::vector< std::vector< std::vector< Double_t> > > cosv1ptx(ybin1);
  Double_t  cosv1pt[ybin1][nbin1];
  Double_t  sinv1pt[ybin1][nbin1];
  for( UInt_t i = 0; i < ybin1; i++) {
    cosv1[i]   = 0.;
    sinv1[i]   = 0.;
    cosv1x[i].clear();
    cosv1ptx[i].resize(nbin1);

    hyphi1[i] = new TH1D( Form("hyphi1_%d",i),rangeLabel1[i]+"#Delta #phi"   , npb, -3.15, 3.15);

    for( UInt_t j = 0; j < nbin1; j++ ){
      cosv1pt[i][j] = 0.;
      sinv1pt[i][j] = 0.;

      hyptphi1[i][j] = new TH1D( Form("hyptphi1_%d%d",i,j),rangeLabel1[i]+"#Delta #phi"   , npb, -3.15, 3.15);
    }
  }

  std::vector< std::vector< std::vector< Double_t> > > cosv2ptx(ybin2);
  Double_t  cosv2pt[ybin2][nbin2];
  Double_t  sinv2pt[ybin2][nbin2];
  for( UInt_t i = 0; i < ybin2; i++) {
    cosv2[i]   = 0.;
    sinv2[i]   = 0.;
    cosv2x[i].clear();
    cosv2ptx[i].resize(nbin2);
    hyphi2[i] = new TH1D( Form("hyphi2_%d",i),rangeLabel2[i]+"2x #Delta #phi", npb/2, 0., 3.15);

    for( UInt_t j = 0; j < nbin; j++ ){
      cosv2pt[i][j] = 0.;
      sinv2pt[i][j] = 0.;
      hyptphi2[i][j] = new TH1D( Form("hyptphi2_%d%d",i,j),rangeLabel2[i]+"2x #Delta #phi"   , npb/2, 0., 3.15);
    }
  }
    
    
  // Lorentz Transform    
  TVector3 boostVec = LorentzBoost(4);
  
  UInt_t vid = 0;
  if( snbm == 108 )
    vid = 1;


  Int_t nevt = SetBranch(0);
  cout << " Number of events " << nevt << endl;
  
  for(Int_t i = 0; i < nevt; i++){

    rChain[0]->GetEntry(i);
      
    /// good event selection
    if(ntrack[2] == 0) continue;

    /// centrality selection
    if(ntrack[4] > cent[icent] || ntrack[4] <= cent[jcent] ) continue;
    Double_t isoratio = (eisobm->Y()+eisotg->Y())/(eisobm->Mod()+eisotg->Mod());
    //if( isoratio > isobin[icent] || isoratio <= isobin[jcent] ) continue;

    heisoratio->Fill( isoratio );
    hntrack->Fill( ntrack[4] );

    Double_t subevt_phi = abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi())); 
    hphi0_180->Fill( subevt_phi );
    if( subevt_phi > TMath::Pi()/2. )
      hphi90_180->Fill( subevt_phi );

      
    TIter next(aNLClusterArray);
    STNeuLANDCluster* aNLClust = NULL;

    while( (aNLClust = (STNeuLANDCluster*)next())) {
	
      aNLClust->SetBeamAngle(ProjA/1000., ProjB/1000.);
      auto pid      = aNLClust->GetPID();
      auto phi      = aNLClust->GetP().Phi();
      auto pt       = aNLClust->GetP().Pt();
      auto rapid    = GetRapidity_cm( aNLClust->GetP(), 939.5731, -boostVec);
      auto veto_all = aNLClust->GetVetoHitAll(vid);
      auto veto_bar = aNLClust->GetVetoHitOne(vid);
      auto veto_mid = aNLClust->GetVetoHitMid(vid);

      rapid -= rapoffset[isys[0]];

      auto dphi     = TVector2::Phi_mpi_pi(phi - unitP_fc->Phi());
      
      if( pid == 2112 && veto_all == 0  ) { //&& 

	hyptacp->Fill( rapid, pt );

	UInt_t irapid = ybin1 - 1;
	for( UInt_t i = 0; i < ybin1; i++){
	  if(rapid < yrange1[i]){
	    irapid = i;
	    break;
	  }
	}
	
	hyphi1[irapid]->Fill(dphi);
	cosv1x[irapid].push_back( rapid );
	cosv1[irapid] += cos(dphi);
	sinv1[irapid] += sin(dphi);

	UInt_t ipt = nbin1 - 1;
	for(UInt_t i = 0; i < nbin1; i++){
	  if( pt < dpt1*(i+1)) {
	    ipt = i;
	    break;
	  }
	}

	hyptphi1[irapid][ipt]->Fill(dphi);
	cosv1ptx[irapid][ipt].push_back( pt );
	cosv1pt[irapid][ipt] += cos(dphi);
	sinv1pt[irapid][ipt] += sin(dphi);

	  
	irapid = ybin2 - 1;

	for( UInt_t i = 0; i < ybin2; i++){
	  if(rapid < yrange2[i]){
	    irapid = i;
	    break;
	  }
	}

	hyphi2[irapid]->Fill(abs(TVector2::Phi_mpi_pi(2.*dphi)));
	cosv2x[irapid].push_back( rapid );
	cosv2[irapid] += cos(2.*dphi);
	sinv2[irapid] += sin(2.*dphi);

	ipt = nbin2 - 1;
	for(UInt_t i = 0; i < nbin2; i++){
	  if( pt < dpt2*(i+1)) {
	    ipt = i;
	    break;
	  }
	}

	hyptphi2[irapid][ipt]->Fill(abs(TVector2::Phi_mpi_pi(2.*dphi)));
	cosv2ptx[irapid][ipt].push_back( pt );
	cosv2pt[irapid][ipt] += cos(2.*dphi);
	sinv2pt[irapid][ipt] += sin(2.*dphi);

      }
    }
  }
      
    
  TGraphErrors *gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  TGraphErrors *gv_v2 = new TGraphErrors();
  gv_v2->SetName("gv_v2");
  
  TGraphErrors *gPt_v1[ybin1];
  TGraphErrors *gPt_v2[ybin1];
    
  for(UInt_t kn = 0; kn < ybin1 ; kn++){      
    gPt_v1[kn] = new TGraphErrors();
    gPt_v1[kn]->SetName((TString)Form("gPt_v1%d",kn));
    TString sname = "neut ; Pt [MeV/c]; v1";
    gPt_v1[kn]->SetTitle(sname);
  }
  
  for(UInt_t kn = 0; kn < ybin2 ; kn++){      
    gPt_v2[kn] = new TGraphErrors();
    gPt_v2[kn]->SetName((TString)Form("gPt_v2%d",kn));
    TString sname = "neut ; Pt [MeV/c]; v2";
    gPt_v2[kn]->SetTitle(sname);
  }


  Double_t *rpres = new Double_t[4];
  rpres = GetRPResolutionwChi(hphi0_180, hphi90_180);
  std::cout << " <cos(Phi)> = " << rpres[0] << " +- " << rpres[1] 
	    << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
	    << std::endl;


  std::cout << " ---- Resutls ---------------------" << std::endl;
  
  UInt_t kl = 0;
  UInt_t id1 = 0;
  UInt_t id2 = 0;
  Double_t para[6];

    
  for(UInt_t kn = 0; kn < ybin1; kn++){
      
    if( cosv1x[kn].size() == 0 ) continue;

    Double_t rapm  = TMath::Mean(cosv1x[kn].begin(), cosv1x[kn].end());
    Double_t rape  = TMath::StdDev(cosv1x[kn].begin(), cosv1x[kn].end());
      
    Double_t yv1  = cosv1[kn] / (Double_t)cosv1x[kn].size();
    Double_t yv1e = sinv1[kn] / (Double_t)cosv1x[kn].size();

    Double_t yv1c  = yv1/rpres[0];
    Double_t yv1ce = GetError(yv1, rpres[0], yv1e, rpres[1]);

    gv_v1->SetPoint( kl,      rapm, yv1c);
    gv_v1->SetPointError( kl, rape, yv1ce);
    kl++;

    // pt dependence 
    UInt_t il = 0; 
    for(UInt_t jn = 0; jn < nbin1; jn++){

      if( cosv1ptx[kn][jn].size() > 0 ) {	

	Double_t ptc  = TMath::Mean(cosv1ptx[kn][jn].begin(), cosv1ptx[kn][jn].end());
	Double_t ptce = TMath::StdDev(cosv1ptx[kn][jn].begin(), cosv1ptx[kn][jn].end());

	Double_t ptv1  = cosv1pt[kn][jn] / (Double_t)cosv1ptx[kn][jn].size();
	Double_t ptv1e = sinv1pt[kn][jn] / (Double_t)cosv1ptx[kn][jn].size();

	yv1c  = ptv1/rpres[0];
	yv1ce = GetError(ptv1, rpres[0], ptv1e, rpres[1]);

	if( kn == 0 ){
	  cout << jn << " erro " << sinv1pt[kn][jn] << " yvc1  " << yv1ce <<  endl;
	  cout << sinv1pt[kn][jn] << " --  " << (Double_t)cosv1ptx[kn][jn].size() << endl;
	}


	gPt_v1[kn]->SetPoint(il, ptc, yv1c);
	gPt_v1[kn]->SetPointError(il, ptce, yv1ce);
	    
	il++;
      }
    }
    gPt_v1[kn]->Write();
  }

  kl = 0;
  for(UInt_t kn = 0; kn < ybin2; kn++){

    if( cosv2x[kn].size() == 0 ) continue;

    Double_t rapm  = TMath::Mean(cosv2x[kn].begin(), cosv2x[kn].end());
    Double_t rape  = TMath::StdDev(cosv2x[kn].begin(), cosv2x[kn].end());

    Double_t yv2  = cosv2[kn] / (Double_t)cosv2x[kn].size();
    Double_t yv2e = sinv2[kn] / (Double_t)cosv2x[kn].size();

    Double_t yv2c  = yv2 / rpres[2];
    Double_t yv2ce = GetError(yv2, rpres[2], yv2e, rpres[3]);


    gv_v2->SetPoint( kl,    rapm, yv2c);
    gv_v2->SetPointError( kl, rape, yv2ce);      
    kl++;

    UInt_t il = 0; 
    for(UInt_t jn = 0; jn < nbin2; jn++){
	
      if( cosv2ptx[kn][jn].size() > 0 ){

	Double_t ptc  = TMath::Mean(cosv2ptx[kn][jn].begin(), cosv2ptx[kn][jn].end());
        Double_t ptce = TMath::StdDev(cosv2ptx[kn][jn].begin(), cosv2ptx[kn][jn].end());

	Double_t ptv2  = cosv2pt[kn][jn] / (Double_t)cosv2ptx[kn][jn].size();
        Double_t ptv2e = sinv2pt[kn][jn] / (Double_t)cosv2ptx[kn][jn].size();

	yv2c  = ptv2/ rpres[2];
	yv2ce = GetError(ptv2, rpres[2], ptv2e, rpres[3]);

	gPt_v2[kn]->SetPoint(il, ptc, yv2c);
	gPt_v2[kn]->SetPointError(il, ptce, yv2ce);
	  
	il++;
      }
    }
    gPt_v2[kn]->Write();
  }

  gv_v1->Write();
  gv_v2->Write();



  //plotting
  ic++; id=1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  cc[ic]->Divide(1,2);

  cc[ic]->cd(id); id++;
  hntrack->Draw("colz");

  cc[ic]->cd(id); id++;
  heisoratio->Draw("colz");  
  

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hyptacp->Draw("colz");
  


  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  cc[ic]->Divide(ybin1,1);
  for(UInt_t kn = 0; kn < ybin1; kn++){
    cc[ic]->cd(id); id++;
    gPt_v1[kn]->Draw("ALP");
  }

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1600,500);
  cc[ic]->Divide(ybin2,1);
  for(UInt_t kn = 0; kn < ybin2; kn++){
    cc[ic]->cd(id); id++;
    gPt_v2[kn]->Draw("ALP");
  }


  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v1->Draw("ALP");

  auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v2->Draw("ALP");

  ic++; 
  cc[ic]   = new TCanvas("dyphi","dphi1 y bin",500,1200);
  cc[ic]->Divide(2, ybin1);

  cc[ic+1]   = new TCanvas("dphi1","dphi1 y bin and pt bin",2000,1200);
  cc[ic+1]->Divide(nbin1, ybin1);

  cc[ic+2] = new TCanvas("dphi2","dphi2 ybin and pt bin",1600,1200);
  cc[ic+2]->Divide(nbin2, ybin2);

  id = 1;
  for(UInt_t i = 0; i < ybin1; i++) {
    cc[ic]->cd(2*(i+1)-1);
    hyphi1[i]->SetNormFactor(npb);
    hyphi1[i]->Draw("e");

    for(UInt_t j = 0; j < nbin1; j++){
      cc[ic+1]->cd(id);    id++;

      
      if(hyptphi1[i][j]->GetEntries() == 0) continue;
      
      hyptphi1[i][j]->SetNormFactor(npb);
      hyptphi1[i][j]->Draw("e");
    }
  } 

  id = 1;
  for(UInt_t i = 0; i < ybin2; i++) {
    cc[ic]->cd(2*(i+1));
    hyphi2[i]->SetNormFactor(npb/2);
    hyphi2[i]->Draw("e");

    for(UInt_t j = 0; j < nbin2; j++){
      cc[ic+2]->cd(id); id++;
      
      if(hyptphi2[i][j]->GetEntries() == 0) continue;
      
      hyptphi2[i][j]->SetNormFactor(npb/2);
      hyptphi2[i][j]->Draw("e");

    }
  } 
  ic += 3;

  hphi0_180->Write();
  hyptacp->Write();
  heisoratio->Write();
  hntrack->Write();

  gSystem->cd("..");

  SaveCanvas(fName);

  //  std::cout << " TOP " << hisoratio->GetEntries()/


  //  ic++; 
  //  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

  //  gROOT->cd();
}


//**************************************************
//**************************************************

void PlotPtDependence(UInt_t selid = 2)       //%% Executable :
{

  std::cout << "----->  PlotPtDependence(" << selid << ")" << std::cout;

  gStyle->SetOptStat(0);

  gSystem->cd("data");

  //  TString fName = Form("YPT_n%2dto%d_",cent[icent],cent[jcent]) + sysName[isys[0]] + "_" + partname[selid]+"_1.root";
  TString fName = Form("YPT_iso%2dto%d_",(UInt_t)isobin[icent]*10,(UInt_t)isobin[jcent]*10) + sysName[isys[0]] + "_" + partname[selid]+".root";
  auto GraphSave = new TFile(fName,"recreate");
  std::cout << "File " << fName << " is created. " << std::endl;

  auto hphi0_180  = new TH1D("hphi0_180" ,"#Phi from  0 to 180; #Phi",100,0.,3.2);
  auto hphi90_180 = new TH1D("hphi90_180","#Phi from 90 to 180; #Phi",100,0.,3.2);

  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  Int_t RPflag = 0;
  if(selid >= 2) RPflag = 10;

  cout << " Particle " << partname[selid] << endl;
  
  // Rapidity binning
  TArrayD arr_range1(ybin1, yrange1);
  TArrayD arr_range2(ybin2, yrange2);

  // PT binning
  Double_t pt_max = 800.;
  UInt_t   nbin1  = 10; //16
  UInt_t   nbin2  = 8; //10
  Double_t dpt1   = pt_max/(Double_t)nbin1;
  Double_t dpt2   = pt_max/(Double_t)nbin2;


  std::cout << " Rapidity binning " << ybin1 << std::endl;
  TString rangeLabel[ybin1];
  for(UInt_t i = 0; i < ybin1; i++ ){
    if( i == 0 )
      rangeLabel[0] = Form(" y < %5.2f ",yrange1[0]);
    else if ( i == ybin1 -1 )
      rangeLabel[i] = Form("%5.2f <= y "    ,yrange1[ybin1-1]);
    else 
      rangeLabel[i] = Form("%5.2f <= y < %5.2f",yrange1[i-1],yrange1[i]);
  }

  
  auto hmass   = new TH2D("hmass",   ";P/Q; Mass [MeV/c^2]"     ,200,  0.,2500., 200, 0.,7000);
  // auto hyptacp = new TH2D("hyptacp", (TString)Form("m %2d to %d ; Y_{cm}; Pt [MeV/c]",cent[icent],cent[jcent]) ,200, -0.4, 0.5, 200, 0., 800);   
  auto hyptacp = new TH2D("hyptacp", (TString)Form("iso %2d to %d ; Y_{cm}; Pt [MeV/c]",
						   (UInt_t)isobin[icent],(UInt_t)isobin[jcent]) ,200, -0.4, 0.5, 200, 0., 800);
  auto hmassy  = new TH2D("hmassy",  "; Y_{cm}; Mass [MeV/C^2]",200, -0.4, 0.5, 200, 0.,7000);
  auto hpy     = new TH2D("hpy",     "; Y_{cm}; P/Q [MeV/c]"   ,200, -0.4, 0.5, 200, 0.,2500.);
  auto hptmass = new TH2D("hptmass", "; Pt [MeV/c]; Mass [MeV/c^2]",200, 0., 800,200, 0.,7000);
  auto hisoratio = new TH1D("hisoratio",";EISO_beam Y",100,0.,1.);

  TH2D *hypt1[ybin1];
  TH2D *hypt2[ybin1];
  TH2D *hyphi1[ybin1];
  TH2D *hyphi2[ybin2];
  TH2D *hyptphi1[ybin1][nbin1];
  TH2D *hyptphi2[ybin2][nbin2];
    
  // Lorentz Transform    
  TVector3 boostVec = LorentzBoost(4);
  
  for(UInt_t kn = 0; kn < ybin1; kn++){ 

    TString sname = rangeLabel[kn];
    hypt1[kn]   = new TH2D((TString)Form("hypt1_%d",kn),   sname+"; Rapidity; Pt [MeV/c]"  , 200, -0.4, 0.5, 200, 0., 800);
    hyphi1[kn]  = new TH2D((TString)Form("hyphi1_%d",kn),  sname+"; #Phi"   , 100, -3.2,  3.2, 200, 0., 800);

    for(UInt_t pn = 0; pn < nbin1; pn++)
      hyptphi1[kn][pn] = new TH2D((TString)Form("hyptphi1_%d%d",kn,pn), sname+"; #Delta #Phi" , 100, -3.2, 3.2, 200, 0., 1000.); 
  }


  for(UInt_t kn = 0; kn < ybin2; kn++){ 
    TString sname = rangeLabel[kn];
    hypt2[kn]   = new TH2D((TString)Form("hypt2_%d",kn),   sname+"; Rapidity; Pt [MeV/c]"  , 200, -0.4, 0.5, 200, 0., 800);
    hyphi2[kn]  = new TH2D((TString)Form("hyphi2_%d",kn),  sname+";2*#Phi"  ,  50,  0.,   3.2, 200, 0., 800);
    
    for(UInt_t pn = 0; pn < nbin2; pn++)
      hyptphi2[kn][pn] = new TH2D((TString)Form("hyptphi2_%d%d",kn,pn),"pt; #Delta 2*#Phi", 50, 0., 3.2,   200, 0., 1000.); 
  }

  
  Int_t nevt = SetBranch();
  cout << " Number of events " << nevt << endl;
  
  for(Int_t i = 0; i < nevt; i++){
    rChain[0]->GetEntry(i);
      
    if(ntrack[2] == 0) continue;
    //    if(ntrack[4] > cent[icent] || ntrack[4] <= cent[jcent] ) continue;
    Double_t isoratio = eisobm->Y()/eisobm->Mod();
    if( isoratio > isobin[icent] || isoratio <= isobin[jcent] ) continue;

    hisoratio->Fill( isoratio );

    Double_t subevt_phi = abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi())); 
    hphi0_180->Fill( subevt_phi );
    if( subevt_phi > TMath::Pi()/2. )
      hphi90_180->Fill( subevt_phi );


      
    TIter next(aArray);
    STParticle *aPart = NULL;

    while( (aPart = (STParticle*)next()) ) {
	
      auto rpf   = aPart->GetReactionPlaneFlag();
      auto pid   = aPart->GetPID();
      auto bmass = aPart->GetBBMass();
      
      if( pid == partid[selid] ) { //&& 
      //      if(aPart->GetCharge() == pcharge  && aPart->GetIndividualRPAngle() > -9 && bmass <= cutbmass[selid] && bmass >= cutlmass[selid]){

	//	if( pid > 2000 && (rpf == 110 || rpf == 210 ) ) continue;

	auto pt    = aPart->GetRotatedMomentum().Pt();
	auto dphi  = aPart->GetAzmAngle_wrt_RP();
	auto rapid = aPart->GetRapiditycm();;
	

	hmass->Fill( aPart->GetRotatedMomentum().Mag(), bmass);
	
	hyptacp->Fill( rapid, pt );
	hmassy ->Fill( rapid, bmass );
	hpy    ->Fill( rapid, aPart->GetRotatedMomentum().Mag() );
	hptmass->Fill( pt,   bmass );

	UInt_t irapid = ybin1 - 1;
	for( UInt_t i = 0; i < ybin1; i++){
	  if(rapid < yrange1[i]){
	    irapid = i;
	    break;
	  }
	}
	
	hypt1[irapid]->Fill(rapid, pt);
	hyphi1[irapid]->Fill(dphi, pt);
	UInt_t ipt = nbin1 - 1;
	for(UInt_t i = 0; i < nbin1; i++){
	  if( pt < dpt1*(i+1)) {
	    ipt = i;
	    break;
	  }
	}
	hyptphi1[irapid][ipt]->Fill(dphi, pt);
	

	  
	irapid = ybin2 - 1;
	for( UInt_t i = 0; i < ybin2; i++){
	  if(rapid < yrange2[i]){
	    irapid = i;
	    break;
	  }
	}

	hypt2[irapid]->Fill(rapid, pt);
	hyphi2[irapid]->Fill(abs( TVector2::Phi_mpi_pi(2.*dphi)), pt );
	
	ipt = nbin2 - 1;
	for(UInt_t i = 0; i < nbin2; i++){
	  if( pt < dpt2*(i+1)) {
	    ipt = i;
	    break;
	  }
	}
	hyptphi2[irapid][ipt]->Fill(abs( TVector2::Phi_mpi_pi(2.*dphi) ), pt);
      }
    }
  }
    
  
    
  TGraphErrors *gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  TGraphErrors *gv_v2 = new TGraphErrors();
  gv_v2->SetName("gv_v2");
  
  TGraphErrors *gPt_v1[ybin1];
  TGraphErrors *gPt_v2[ybin1];
    
  for(UInt_t kn = 0; kn < ybin1 ; kn++){      
    gPt_v1[kn] = new TGraphErrors();
    gPt_v1[kn]->SetName((TString)Form("gPt_v1%d",kn));
    TString sname = partname[selid]+"; Pt [MeV/c]; v1";
    gPt_v1[kn]->SetTitle(sname);
  }
  
  for(UInt_t kn = 0; kn < ybin2 ; kn++){      
    gPt_v2[kn] = new TGraphErrors();
    gPt_v2[kn]->SetName((TString)Form("gPt_v2%d",kn));
    TString sname = partname[selid]+"; Pt [MeV/c]; v2";
    gPt_v2[kn]->SetTitle(sname);
  }


  Double_t *rpres = new Double_t[4];
  rpres = GetRPResolutionwChi(hphi0_180, hphi90_180);
  std::cout << " <cos(Phi)> = " << rpres[0] << " +- " << rpres[1] 
	    << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
	    << std::endl;

  ic++; 
  cc[ic]   = new TCanvas("dphi1","dphi1",2000,1200);
  cc[ic]->Divide(nbin1, ybin1);
  
  cc[ic+1] = new TCanvas("dphi2","dphi2",1400,1200);
  cc[ic+1]->Divide(nbin2, ybin2);
  

  std::cout << " ---- Resutls ---------------------" << std::endl;
  
  UInt_t kl = 0;
  UInt_t id1 = 0;
  UInt_t id2 = 0;
  Double_t para[6];

  cc[ic+2] = new TCanvas(Form("cc%d",ic+2),Form("cc%d",ic+2), 1000, 1200);
  cc[ic+2]->Divide(2, ybin1);
    
  for(UInt_t kn = 0; kn < ybin1; kn++){
      
    Double_t rapm  = hypt1[kn]->GetMean(1);
    Double_t rape  = hypt1[kn]->GetStdDev(1);

    // v1 vs rapidity
    cc[ic+2]->cd(2*(kn+1)-1);
    auto hyphi = (TH1D*)hyphi1[kn]->ProjectionX(Form("hyphiv1_%d",kn),0,-1,"eo");
    
    Double_t corr[2]={rpres[0], rpres[1]};
    GetFittingParameters(*hyphi, para, corr);
    //    hyphi->SetMaximum(1.5);
    //    hyphi->SetMinimum(0.5);
      
    gv_v1->SetPoint( kl,     rapm, para[4]);
    gv_v1->SetPointError( kl, rape, para[5]);
    kl++;

    // pt dependence 
    UInt_t il = 0; 
    for(UInt_t jn = 0; jn < nbin1; jn++){

      id1++; cc[ic]->cd(id1); 

      if( hyptphi1[kn][jn]->GetEntries() > 0 ) {	

	Double_t ptc  = hyptphi1[kn][jn]->GetMean(2);
	Double_t ptce = hyptphi1[kn][jn]->GetStdDev(2);
	    
	auto hypt = (TH1D*)hyptphi1[kn][jn]->ProjectionX();
	GetFittingParameters(*hypt, para, corr);
	    

	gPt_v1[kn]->SetPoint(il, ptc, para[4]);
	gPt_v1[kn]->SetPointError(il, ptce, para[5]);
	    
	il++;
      }
    }
    gPt_v1[kn]->Write();
  }

  kl = 0;
  for(UInt_t kn = 0; kn < ybin2; kn++){

    Double_t rapm  = hypt2[kn]->GetMean(1);
    Double_t rape  = hypt2[kn]->GetStdDev(1);

    // v2 vs rapidity
    cc[ic+2]->cd(2*(kn+1));

    auto hyphii = (TH1D*)hyphi2[kn]->ProjectionX(Form("hyphiv2_%d",kn),0,-1,"eo");
      
    Double_t corr[2] = {rpres[2], rpres[3]};
    GetFittingParameters(*hyphii, para, corr);
    //    hyphii->SetMaximum(1.1);
    //    hyphii->SetMinimum(0.9);

    gv_v2->SetPoint( kl,    rapm, para[4]);
    gv_v2->SetPointError( kl, rape, para[5]);      
    kl++;

    UInt_t il = 0; 
    for(UInt_t jn = 0; jn < nbin2; jn++){
	
      id2++; cc[ic+1]->cd(id2); 
      if( hyptphi2[kn][jn]->GetEntries() > 0 ){

	Double_t ptc  = hyptphi2[kn][jn]->GetMean(2);
	Double_t ptce = hyptphi2[kn][jn]->GetStdDev(2);

	auto hypt = (TH1D*)hyptphi2[kn][jn]->ProjectionX();

	GetFittingParameters(*hypt, para, corr);
	//	hypt->SetMaximum(1.1);
	//	hypt->SetMinimum(0.9);

	gPt_v2[kn]->SetPoint(il, ptc, para[4]);
	gPt_v2[kn]->SetPointError(il, ptce, para[5]);
	  
	il++;
      }
    }
    gPt_v2[kn]->Write();
  }

  gv_v1->Write();
  gv_v2->Write();

  ic+=3;


  //plotting
  ic++; id=1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  cc[ic]->Divide(2,2);

  cc[ic]->cd(id); id++;
  hmass->Draw("colz");

  cc[ic]->cd(id); id++;
  hmassy->Draw("colz");

  cc[ic]->cd(id); id++;
  hpy->Draw("colz");

  // cc[ic]->cd(id); id++;
  // hptmass->Draw("colz");

  cc[ic]->cd(id); id++;
  hisoratio->Draw("colz");  
  

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hyptacp->Draw("colz");
  


  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
  cc[ic]->Divide(ybin1,1);
  for(UInt_t kn = 0; kn < ybin1; kn++){
    cc[ic]->cd(id); id++;
    gPt_v1[kn]->Draw("ALP");
  }

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
  cc[ic]->Divide(ybin2,1);
  for(UInt_t kn = 0; kn < ybin2; kn++){
    cc[ic]->cd(id); id++;
    gPt_v2[kn]->Draw("ALP");
  }

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,400);
  cc[ic]->Divide(ybin1,1);
  for(UInt_t kn = 0; kn < ybin1; kn++){
    cc[ic]->cd(id); id++;
    hypt1[kn]->Draw("colz");
  }

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,400);
  cc[ic]->Divide(ybin2,1);
  for(UInt_t kn = 0; kn < ybin2; kn++){
    cc[ic]->cd(id); id++;
    hypt2[kn]->Draw("colz");
  }

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hypt1[0]->Draw("colz");
  for(UInt_t kn = 1; kn < ybin1; kn++)
    hypt1[kn]->Draw("same colz");

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hypt2[0]->Draw("colz");
  for(UInt_t kn = 1; kn < ybin2; kn++)
    hypt2[kn]->Draw("same colz");


  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v1->Draw("ALP");

  auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v2->Draw("ALP");

  hphi0_180->Write();
  hyptacp->Write();
  GraphSave->WriteObject(&arr_range1, "v1_yrange");
  GraphSave->WriteObject(&arr_range2, "v2_yrange");
  
  gSystem->cd("..");


  ic++; 
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

  //  gROOT->cd();
}

//**************************************************
//**************************************************

Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180)            //%% Executable : 
{
  hphi90_180->SetLineColor(2);

  Double_t m0 = hphi0_180->GetEntries();
  Double_t m1 = hphi90_180->GetEntries();

  Double_t chi = sqrt(-2.* log(2.* m1/m0));
  
  Double_t m01e = m1/m0*sqrt(1./m1+1./m0);
  Double_t chie = chi - sqrt(-2.* log(2. * (m1/m0+m01e)));

  Double_t *rpres = new Double_t[4];

  rpres[0] = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*chi;
  rpres[1] = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*(chi+chie) - rpres[0];
  rpres[2] = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi,2);
  rpres[3] = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi+chie,2) - rpres[2];

  // ic++;
  // cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  // hphi0_180->Draw();
  // hphi90_180->Draw("same");

  //  Double_t para[6];
  //  GetFittingParameters( *hphi0_180, para);

  //  cout << std::setw(14) << para[1] << " +- " << std::setw(10) << para[3]
  //       << std::endl;

  return rpres;
}


void PlotT3He()         //
{

  // cluster ratio is in progress.

  gStyle->SetOptStat(0);

  //  TVector3 boostVec = LorentzBoost(4);

  auto haccpTri    = new TH2D("haccpTri"   , "Tri; Rapidity ; Pt [MeV/c]",100, -0.4, 0.6, 100,  0.,800.);
  auto haccp3He    = new TH2D("haccp3He"   , "3He; Rapidity ; Pt [MeV/c]",100, -0.4, 0.6, 100,  0.,800.);
  auto haccp4He    = new TH2D("haccp4He"   , "4He; Rapidity ; Pt [MeV/c]",100, -0.4, 0.6, 100,  0.,800.);

  auto hmT3He      = new TH2I("hmT3He"     , " Multiplicity; Number of Trition; Number of 3He",15, 0, 15, 15, 0, 15); 

  auto hphiDlt = new TH1D("hphiDlt",";#Delta( #phi_{tri} - #phi_{3He}",100,-3.2,3.2);
  auto hphiTri = new TH1D("hphiTri","Triton ;#phi",100,-3.2,3.2);
  auto hphi3He = new TH1D("hphi3He","3He ;#phi",100,-3.2,3.2);
  auto hphiCll = new TH2D("hphiCll","; #phi_{tri}; #phi_{3He}",100,-3.2,3.2, 100,-3.2,3.2);

  std::vector< TVector3 > momTri;
  std::vector< TVector3 > mom3He;
  std::vector< TVector3 > mom4He;


  Int_t nevt = SetBranch(0);
  for(Int_t i = 0; i < nevt; i++){
    rChain[0]->GetEntry(i);

    TIter next(aArray);
    STParticle *aPart = NULL;
    momTri.clear();
    mom3He.clear();
    mom4He.clear();


    auto Psi = unitP_fc->Phi();
    TVector2 momVTri(0.,0.);
    TVector2 momV3He(0.,0.);
  
    while( (aPart = (STParticle*)next()) ) {

      auto chr = aPart->GetCharge();
      auto rpf = aPart->GetReactionPlaneFlag();
      auto pid = aPart->GetPID();
      auto mom = aPart->GetRotatedMomentum();
      auto pt  = mom.Pt();
      auto rapid = aPart->GetRapiditycm();
      
      mom.RotateY(-Psi);
      //      TVector2 momPt = mom

      // STRecoTrack *atrack = (STRecoTrack*)aPart->GetRecoTrack();
      // Double_t nclst = (Double_t)atrack->GetClusterIDArray()->size();
      // Double_t cclst = (Double_t)db->GetClusterNum(chr, mom.Theta(), mom.Phi(), mom.Mag());
      // cout << " nnumber of cluster " << nclst
      // 	   << " expected cluster " << cclst
      // 	//	   << " ratio " << nclst/cclust
      // 	   << endl;

      
      if( pid == partid[4]) {
	haccpTri -> Fill(rapid, pt);

	if( rapid > 0 ) {
	  hphiTri -> Fill(aPart->GetAzmAngle_wrt_RP());
	  momVTri += mom.XYvector();
	}
	else {
 	  hphiTri -> Fill(TVector2::Phi_mpi_pi( aPart->GetAzmAngle_wrt_RP() - TMath::Pi() ) );
	  momVTri += mom.XYvector().Rotate(TMath::Pi());; 
	}
	momTri.push_back(mom);
      }
      else if( pid == partid[5]) {
	haccp3He -> Fill(rapid, pt);

	if( rapid > 0 ) {
	  hphi3He -> Fill(aPart->GetAzmAngle_wrt_RP());
	  momV3He += mom.XYvector();
	}
	else {
 	  hphi3He -> Fill(TVector2::Phi_mpi_pi( aPart->GetAzmAngle_wrt_RP() - TMath::Pi() ) );
	  momV3He += mom.XYvector().Rotate(TMath::Pi());; 
	}

	mom3He.push_back(mom);
      }
      else if( pid == partid[6]) {
	haccp4He -> Fill(rapid, pt);

	mom4He.push_back(mom);
      }
    }

    hmT3He->Fill( momTri.size(), mom3He.size() );

    if( momTri.size() > 0 && mom3He.size() ) {
      Double_t phitri = TVector2::Phi_mpi_pi(momVTri.Phi());
      Double_t phi3he = TVector2::Phi_mpi_pi(momV3He.Phi()); 

      hphiTri->Fill(phitri);
      hphi3He->Fill(phi3he);
      
      auto dphi = momVTri.DeltaPhi(momV3He);
      hphiDlt->Fill(dphi);
      
      hphiCll->Fill( Psi, dphi );
    }
  }

  ic++;
  auto ccv = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  ccv->Divide(1,3);
  ccv->cd(1);
  haccpTri->Draw("colz");
  ccv->cd(2);
  haccp3He->Draw("colz");
  ccv->cd(3);
  haccp4He->Draw("colz");

  ic++;
  ccv = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  ccv->Divide(2,2);
  ccv->cd(1);
  hmT3He->Draw("colz");

  ccv->cd(2);
  hphiDlt->Draw();

  ccv->cd(3);
  hphiTri->Draw();

  ccv->cd(4);
  hphi3He->Draw();

  ic++;
  ccv = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hphiCll->Draw("colz");

  
  gROOT->cd();
}


void JPS()
{
  UInt_t m = 0;

  auto hphir0 = new TH1D("hphir0","fRapidity<0.1; #Delta#phi; dN/d#phi",100,-180.,180);
  rChain[m]->Project("hphi0","fdeltphi","fRapidity<0.1");

  auto hphir1 = new TH1D("hphir1","Rapidity(0.36 ~ 0.46); #Delta#phi; dN/d#phi",100,-180.,180);
  rChain[m]->Project("hphi1","fdeltphi","fRapidity>0.36&&fRapidity<0.46");

  auto hphir2 = new TH1D("hphir2","fRapidity>0.7; #Delta#phi; dN/d#phi",100,-180.,180);
  rChain[m]->Project("hphi2","fdeltphi","fRapidity>0.7");



}


void FlatteningCheck()            
{
  //----- Parametres
  Int_t ntrack[7];

  //----- Booking
  TH2D *hphitheta[4];
  TH2D *hphimtrck[4];

  for(Int_t m = m_bgn; m < m_end; m++){
    TString hname = Form("hphitheta%d",m);
    hphitheta[m] = new TH2D(hname, sysName[isys[m]]+"; #Theta ; #Phi",100,0,0.8, 100,-3.2, 3.2);

    hname = Form("hphimtrck%d",m);
    hphimtrck[m] = new TH2D(hname, sysName[isys[m]]+"; Number of Track ; #Phi",60,0,60, 100,-3.2, 3.2);
  }


  //----- Filling
  for(Int_t m = m_bgn; m < m_end; m++){


    Int_t nEntry = SetBranch(m);

    for(Int_t i = 0; i < nEntry; i++){
      aArray->Clear();

      rChain[m]->GetEntry(i);

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	// auto pid   = aPart->GetPID();
	// auto charg = aPart->GetCharge();
	// auto rapid = aPart->GetRapidity();
	// auto vp    = aPart->GetFlattenMomentum();
	// auto dltphi= aPart->GetAzmAngle_wrt_RP();;
	auto phi   = aPart->GetIndividualRPAngle();
	auto theta = aPart->GetRotatedMomentum().Theta();
	auto flag  = aPart->GetReactionPlaneFlag();

	//	if(flag > 110 ){
	if(flag >= selReactionPlanef ){
	  hphitheta[m]->Fill( theta, phi );
	  hphimtrck[m]->Fill( ntrack[4], phi ); 
	}
      }
    }
  }


  //----- Drawing 

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  cc[ic]->Divide(2,m_end);

  UInt_t id = 1;

  for(Int_t m = m_bgn; m < m_end; m++){
    cc[ic]->cd(id); id++;
    hphitheta[m]->Draw("colz");

    cc[ic]->cd(id); id++;
    hphimtrck[m]->Draw("colz");
  }

}



void PlotSubEvent(Double_t ml, Double_t mu)   
{

  Double_t mlt[] = {ml, mu};
  //    Double_t mlt[2] = {0., 8.};
  //Double_t mlt[2] = {8., 16.};
  // Double_t mlt[2] = {16., 24.};
  // Double_t mlt[2] = {24., 32.};
  // Double_t mlt[2] = {32., 40.};
  // Double_t mlt[2] = {40., 100.};

  TCut mcrot = Form("ntrack[4]>%f&&ntrack[4]<%f",mlt[0]*2.,mlt[1]*2.);
  TCut mc1r  = Form("mtrack_1>%f&&mtrack_1<%f"  ,mlt[0],mlt[1]);
  TCut mc2r  = Form("mtrack_2>%f&&mtrack_2<%f"  ,mlt[0],mlt[1]);

  TString sname;
  
  sname = mcrot.GetTitle();
  auto *hrotx = new TH1D("hrotx","All   "+sname+ ";Qx" ,100,-12.,12.);
  auto *h1rx  = new TH1D("h1rx", "sub_1 "+sname+ ";Qx" ,100,-12.,12.);
  auto *h2rx  = new TH1D("h2rx", "sub_2 "+sname+ ";Qx" ,100,-12.,12.);
  auto *hroty = new TH1D("hroty","All   "+sname+ ";Qy" ,100,-12.,12.);
  auto *h1ry  = new TH1D("h1ry", "sub_1 "+sname+ ";Qy" ,100,-12.,12.);
  auto *h2ry  = new TH1D("h2ry", "sub_2 "+sname+ ";Qy" ,100,-12.,12.);
 

  rChain[0]->Project("hrotx","unitP2_rot.X()",mcrot);
  rChain[0]->Project("h1rx" ,"unitP_1r.X()"  ,mc1r);
  rChain[0]->Project("h2rx" ,"unitP_2r.X()"  ,mc2r);
					      
  rChain[0]->Project("hroty","unitP2_rot.Y()",mcrot);
  rChain[0]->Project("h1ry" ,"unitP_1r.Y()"  ,mc1r);
  rChain[0]->Project("h2ry" ,"unitP_2r.Y()"  ,mc2r);


  //----- Drawing                                                                                                                          
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  hrotx->SetLineColor(2);
  hrotx->SetNormFactor(1);
  h1rx ->SetLineColor(4);
  h1rx ->SetNormFactor(1);
  h2rx ->SetLineColor(6);
  h2rx ->SetNormFactor(1);

  h1rx->Draw();
  h2rx->Draw("same");
  hrotx->Draw("same");

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  hroty->SetLineColor(2);
  hroty->SetNormFactor(1);
  h1ry ->SetLineColor(4);
  h1ry ->SetNormFactor(1);
  h2ry ->SetLineColor(6);
  h2ry ->SetNormFactor(1);

  h1ry->Draw();
  h2ry->Draw("same");
  hroty->Draw("same");
}



void PlotNeuLANDPsi()                         //%%
{
  TH2D *hACCp[4];
  TH2D *hnlACCp[4];
  TH2D *hnlACCn[4];
  TH1D *hfcrn[4][4];
  TH1D *hfcrp[4][4];
  TH1D *hfcrnd[4][4];

  TCut ncCut[4];
  ncCut[0]="";
  ncCut[1]="ncPID==2112";
  ncCut[2]=ncCut[1]&&"ncRapidity<=0.37";
  ncCut[3]=ncCut[1]&&"ncRapidity>0.37";

  TCut pcCut[4];
  pcCut[0]="";
  pcCut[1]="ncPID==2212";
  pcCut[2]=pcCut[1]&&"ncRapidity<=0.37";
  pcCut[3]=pcCut[1]&&"ncRapidity>0.37";

  //----- Booking            
  TString hname;
  for(Int_t m = m_bgn; m < m_end; m++){
    hname = Form("hnlACCn%d",m);
    hnlACCn[m] = new TH2D(hname, hname,  200, 0., 1., 200., 0., 500.);
    hnlACCn[m]->SetMarkerColor(4);
    rChain[m]->Project(hname,"ncP.Pt():ncRapidity","ncPID==2112");

    hname = Form("hnlACCp%d",m);
    hnlACCp[m] = new TH2D(hname, hname,  200, 0., 1., 200., 0., 500.);
    hnlACCp[m]->SetMarkerColor(2);
    rChain[m]->Project(hname,"ncP.Pt():ncRapidity","ncPID==2212");

    hname = Form("hACCp%d",m);
    hACCp[m] = new TH2D(hname, hname,  200, 0., 1., 200., 0., 500.);
    hACCp[m]->SetTitle("; Rapidity; Pt [GeV/c]");
    rChain[m]->Project(hname,"fRotatedP3.Pt():fRapidity","fPID==2212&&fReactionPlanef>100");

    for(UInt_t k = 0; k < 4; k++){
      hname = Form("hfcrn%d%d",m,k);
      hfcrn[m][k] = new TH1D(hname, hname+ncCut[k].GetTitle(), 30, -3.15, 3.15);
      rChain[m]->Project(hname,"TVector2::Phi_mpi_pi(unitP_fc.Phi()-ncP.Phi())",ncCut[k]);

      hname = Form("hfcrnd%d%d",m,k);
      hfcrnd[m][k] = new TH1D(hname, hname+ncCut[k].GetTitle(), 60, -180., 180.);
      rChain[m]->Project(hname,"TVector2::Phi_mpi_pi(unitP_fc.Phi()-ncP.Phi())*180./TMath::Pi()",ncCut[k]);

      hname = Form("hfcrp%d%d",m,k);
      hfcrp[m][k] = new TH1D(hname, hname+pcCut[k].GetTitle(), 60,  -180., 180.);
      rChain[m]->Project(hname,"TVector2::Phi_mpi_pi(unitP_fc.Phi()-ncP.Phi())*180./TMath::Pi()",pcCut[k]);
    }

    if(kFALSE){
      ic++; id = 1;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

      hACCp[m]->Draw("colz");
      hnlACCn[m]->Draw("same");

      auto PRLabel = new TLatex(0.7,50,"TPC proton");
      auto NLLabel = new TLatex(0.2,100,"Neutron");
      NLLabel->SetTextColor(0);
      PRLabel->SetTextColor(0);
      NLLabel->Draw();
      PRLabel->Draw();


      ic++; id = 1;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
      cc[ic]->Divide(4,2);

      for(UInt_t k = 0; k < 4; k++){
	cc[ic]->cd(id); id++;
	hfcrn[m][k]->SetLineColor(4);
	hfcrn[m][k]->Draw("e");
      }
      for(UInt_t k = 0; k < 4; k++){
	cc[ic]->cd(id); id++;
	hfcrp[m][k]->SetLineColor(2);
	hfcrp[m][k]->Draw("e");
      }
    }

    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    hfcrnd[m][2]->SetLineColor(2);
    hfcrnd[m][2]->SetMarkerColor(2);
    hfcrnd[m][2]->SetNormFactor(60);
    hfcrnd[m][2]->SetMarkerStyle(20);
    hfcrnd[m][2]->SetTitle(";#Psi;1/NdN/d#Psi");
    hfcrnd[m][3]->SetLineColor(4);
    hfcrnd[m][3]->SetMarkerColor(4);
    hfcrnd[m][3]->SetNormFactor(60);
    hfcrnd[m][3]->SetMarkerStyle(20);
    hfcrnd[m][3]->SetTitle(";#Psi;1/NdN/d#Psi");
    hfcrp[m][3]->SetLineColor(8);
    hfcrp[m][3]->SetNormFactor(60);
    hfcrp[m][3]->SetMarkerColor(8);
    hfcrp[m][3]->SetMarkerStyle(20);;

    hfcrp[m][3]->Draw("e");
    hfcrnd[m][3]->Draw("samee");
    hfcrnd[m][2]->Draw("samee");



    auto aLeg = new TLegend(0.65,0.7,0.9,0.9,"");
    aLeg->AddEntry(hfcrnd[m][2] ,"Neutron hit y < y_cm","lp");
    aLeg->AddEntry(hfcrnd[m][3] ,"Neutron hit y > y_cm","lp");
    aLeg->AddEntry(hfcrp[m][3] ,"Proton on NL y > y_cm","lp");

    aLeg->Draw();


  }
}


//--------------------------------------------//%% Executable : 
void PlotNeuLANDProperty(UInt_t iout)           //%% Executable : 
{
  std::cout << " Executing : PlotNeuLANDProperty(" << iout << ")" << std::endl;

  //----- Parametres                                                                                                                       
  UInt_t nybin  = 100;
  TFile* hout;
  UInt_t total_neut = 0;
 
  //----- Canvas
  if( iout == 1 ){
    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  }

  //----- Booking
  auto hmult  = new TH1I("hmult0","Neutron; Multiplicity ", 14,0.,14);
  auto hmultc = new TH1I("hmultc","Charged particles; Multiplicity ", 8,0,8);
  auto hmultnc= new TH2I("hmultnc"," ;Neutron Multiplicity; Charged Particle Multiplicity",14,0,14,8,0,8);
  auto hmulttn= new TH2I("hmulttn"," ;TPC multiplicity; Neutron Multiplicity",65,0,65,14,0,14);
  auto hmulttc= new TH2I("hmulttc"," ;TPC multiplicity; NeuLAND Charged Particle Multiplicity",65,0,65,8,0,8);

  auto haccp  = new TH2D("haccp0","; Rapidity ; Pt [MeV/c]",nybin, -0.4, 0.5, 100,   0., 800);

  //----- Event loop
  for(Int_t m = m_bgn; m < m_end; m++){

    total_neut = 0;
    TVector3 boostVec = LorentzBoost(4);

  //----- Output file
    if( iout == 1 ) {
      TString fName = "NLdbProp.cm" + sysName[isys[m]] + ".root";
      gSystem->cd("data");
      hout = new TFile(fName,"recreate");
    }
    
    hmult->Reset();
    haccp->Reset();
    hmultc->Reset();

    rChain[m]->SetBranchAddress("ntrack",ntrack);
    rChain[m]->SetBranchAddress("STNeuLANDCluster", &aNLClusterArray);
    Int_t nevt = rChain[m]->GetEntries();

    cout << " Number of events " << nevt << endl;
    
    UInt_t vID = 0;
    if( snbm == 108 )
      vID = 1;

    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      UInt_t ncharge = 0;
      UInt_t nneut = 0;
      TIter next(aNLClusterArray);
      STNeuLANDCluster* aNLClust = NULL;

      while( (aNLClust = (STNeuLANDCluster*)next()) ){

	aNLClust->SetBeamAngle(ProjA/1000., ProjB/1000.);
        auto pid      = aNLClust->GetPID();
        auto pt       = aNLClust->GetP().Pt();
        auto veto_all = aNLClust->GetVetoHitAll(vID);
        auto veto_bar = aNLClust->GetVetoHitOne(vID);
	
	if( pid == 2112 && aNLClust->GetTOF()>40 ){
	  auto rapidity = GetRapidity_cm( aNLClust->GetP(), 939.5731, -boostVec);
	  haccp->Fill(rapidity, pt);
	  nneut++;
	  total_neut++;
	}

	else if( pid != 2112 && pid > 0 && veto_all==1 )
	  ncharge++;

      }
      
      hmult->Fill(nneut);
      hmultc->Fill(ncharge);
      hmultnc->Fill(nneut,ncharge);

      hmulttn->Fill(ntrack[4], nneut);
      hmulttc->Fill(ntrack[4], ncharge);
    }
  
    
    //    hmult->Draw();

    if(iout == 1){
      hmult->Write();
      haccp->Write();

      hout->Close();
      gSystem->cd("..");
    }

    gStyle->SetOptStat(0);
    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    haccp->Draw("colz");


    ic++;
    auto cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    cc->SetLogy();
    hmult->Draw();

    ic++;
    cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    cc->SetLogy();
    hmultc->Draw();

    ic++;
    cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    cc->SetLogz();
    hmultnc->Draw("colz");

    ic++;
    cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    cc->SetLogz();
    hmulttn->Draw("colz");

    ic++;
    cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    cc->SetLogz();
    hmulttc->Draw("colz");


    std::cout << " Number of neutron " << total_neut << " / " << hmult->GetEntries() 
	      << " = " << Float_t(total_neut/hmult->GetEntries())
	      << " : " << nevt << std::endl;
  }
}



void PlotNeuLANDv1v2()                        //%% Executable : 
{
  std::cout << "-----> PlotNeuLANDv1v2() ----->" << std::endl;

  //----- Parametres                                                                                                                       
  TFile *hout;

  auto hphi0_180  = new TH1D("hphi0_180" ,"#Phi from  0 to 180; #Phi",100,0.,3.2);
  auto hphi90_180 = new TH1D("hphi90_180","#Phi from 90 to 180; #Phi",100,0.,3.2);

  auto hisoratio = new TH1D("hisoratio",";EISO_beam Y",100,0.,1.);


  // PT binning                                                                                                                            
  Double_t pt_max = 800.;
  UInt_t nbin1 = 10; //16                                                                                                                  
  UInt_t nbin2 = 8; //10    
  Double_t dpt1 = pt_max/(Double_t)nbin1;
  Double_t dpt2 = pt_max/(Double_t)nbin2;

  std::cout << " Rapidity binning " << ybin1 << std::endl;
  TString rangeLabel[ybin1];
  for(UInt_t i = 0; i < ybin1; i++ ){
    if( i == 0 )
      rangeLabel[0] = Form(" y < %f ",yrange1[0]);
    else if ( i == ybin1 -1 )
      rangeLabel[i] = Form("%f <= y "    ,yrange1[ybin1-1]);
    else
      rangeLabel[i] = Form("%f <= y < %f",yrange1[i-1],yrange1[i]);
  }

  auto hyptacp = new TH2D("hyptacp", "; Y_{cm}; Pt [MeV/c]"    ,200, -0.4, 0.5, 200, 0., 800);

  TH2D *hypt[ybin1];
  TH2D *hypt2[ybin1];
  TH2D *hyphi1[ybin1];
  TH2D *hyphi2[ybin2];
  TH2D *hyptphi1[ybin1][nbin1];
  TH2D *hyptphi2[ybin2][nbin2];


  // Lorentz Transform                                                                                                                    
  TVector3 boostVec = LorentzBoost(4);

  for(UInt_t kn = 0; kn < ybin1; kn++){

    TString sname = rangeLabel[kn];
    hypt[kn]    = new TH2D((TString)Form("hypt_%d",kn),    sname+"; Rapidity; Pt [MeV/c]"  , 200, -0.4, 0.5, 200, 0., 800);
    hyphi1[kn]  = new TH2D((TString)Form("hyphi1_%d",kn),  sname+"; #Phi"   , 100, -3.2,  3.2, 200, 0., 800);
    
    for(UInt_t pn = 0; pn < nbin1; pn++)
      hyptphi1[kn][pn] = new TH2D((TString)Form("hyptphi1_%d%d",kn,pn),"pt; #Delta #Phi" , 100, -3.2, 3.2, 200, 0., 1000.);
  }

  for(UInt_t kn = 0; kn < ybin2; kn++){
    TString sname = rangeLabel[kn];
    hypt2[kn]   = new TH2D((TString)Form("hypt2_%d",kn),   sname+"; Rapidity; Pt [MeV/c]"  , 200, -0.4, 0.5, 200, 0., 800);
    hyphi2[kn]  = new TH2D((TString)Form("hyphi2_%d",kn),  sname+";2*#Phi"  ,  50,  0.,   3.2, 200, 0., 800);

    for(UInt_t pn = 0; pn < nbin2; pn++)
      hyptphi2[kn][pn] = new TH2D((TString)Form("hyptphi2_%d%d",kn,pn),"pt; #Delta 2*#Phi", 50, 0., 3.2,   200, 0., 1000.);
  }


    
  UInt_t vid = 0;
  if( snbm == 108 )
    vid = 1;


  Int_t nevt = SetBranch(0);

  cout << " Number of events " << nevt << endl;

  for(Int_t i = 0; i < nevt; i++){
    rChain[0]->GetEntry(i);

    //    if(ntrack[4] > cent[icent] || ntrack[4] <= cent[jcent] ) continue;
    Double_t isoratio = eisobm->Y()/eisobm->Mod();
    if( isoratio > isobin[icent] || isoratio <= isobin[jcent] ) continue;
    hisoratio->Fill(isoratio);

    // Reaction Plane resolution evaluated with subevents
    Double_t subevt_phi = abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi()));
    hphi0_180->Fill( subevt_phi );
    if( subevt_phi > TMath::Pi()/2. )
      hphi90_180->Fill( subevt_phi );


    TIter next(aNLClusterArray);
    STNeuLANDCluster* aNLClust = NULL;
    
    while( (aNLClust = (STNeuLANDCluster*)next()) ){
      
      aNLClust->SetBeamAngle(ProjA/1000., ProjB/1000.);
      
      auto pid      = aNLClust->GetPID();
      auto phi      = aNLClust->GetP().Phi(); 
      auto pt       = aNLClust->GetP().Pt();
      
      auto rapid    = GetRapidity_cm( aNLClust->GetP(), 939.5731, -boostVec);
      auto veto_all = aNLClust->GetVetoHitAll(vid);
      auto veto_bar = aNLClust->GetVetoHitOne(vid);
      auto veto_mid = aNLClust->GetVetoHitMid(vid);
      

      auto dphi     = TVector2::Phi_mpi_pi(unitP_fc->Phi() - phi);
  	
      if( pid == 2112 && veto_all == 0){ // neutron
	//   	if( pid == 2112 ){ // neutron


	hyptacp->Fill( rapid, pt );

	UInt_t irapid = ybin1 - 1;
	for( UInt_t i = 0; i < ybin1; i++){
	  if(rapid < yrange1[i]){
	    irapid = i;
	    break;
	  }
	}
	
	hypt[irapid]->Fill(rapid, pt);
	hyphi1[irapid]->Fill(dphi, pt);
	UInt_t ipt = nbin1 - 1;
	for(UInt_t i = 0; i < nbin1; i++){
	  if( pt < dpt1*(i+1)) {
	    ipt = i;
	    break;
	  }
	}

	hyptphi1[irapid][ipt]->Fill(dphi, pt);
	
	irapid = ybin2 - 1;
	for( UInt_t i = 0; i < ybin2; i++){
	  if(rapid < yrange2[i]){
	    irapid = i;
	    break;
	  }
	}
          
	ipt = nbin2 - 1;
	for(UInt_t i = 0; i < nbin2; i++){
	  if( pt < dpt2*(i+1)) {
	    ipt = i;
	    break;
	  }
	}
	hypt2[irapid]->Fill(rapid, pt);
	hyphi2[irapid]->Fill(abs( TVector2::Phi_mpi_pi(2.*dphi) ), pt);
	hyptphi2[irapid][ipt]->Fill(abs( TVector2::Phi_mpi_pi(2.*dphi) ), pt);
      }
    }
  }
  


  *hyptphi1[0][0] = *hyptphi1[0][0] + *hyptphi1[0][1];
  hyptphi1[0][1]->Reset();

  *hyptphi1[1][1] = *hyptphi1[1][1] + *hyptphi1[1][2];
  hyptphi1[1][2]->Reset();

  *hyptphi1[2][1] = *hyptphi1[2][1] + *hyptphi1[2][3];
  hyptphi1[2][3]->Reset();

  *hyptphi1[3][1] = *hyptphi1[3][1] + *hyptphi1[3][2];
  hyptphi1[3][2]->Reset();

  *hyptphi1[3][3] = *hyptphi1[3][3] + *hyptphi1[3][4];
  hyptphi1[3][4]->Reset();
    
  *hyptphi1[4][5] = *hyptphi1[4][5] + *hyptphi1[4][6];
  hyptphi1[4][6]->Reset();

  *hyptphi1[5][2] = *hyptphi1[5][2] + *hyptphi1[5][3];
  hyptphi1[5][3]->Reset();

  *hyptphi1[5][6] = *hyptphi1[5][6] + *hyptphi1[5][7];
  hyptphi1[5][7]->Reset();
  *hyptphi1[5][6] = *hyptphi1[5][6] + *hyptphi1[5][8];
  hyptphi1[5][8]->Reset();

  *hyptphi1[6][3] = *hyptphi1[6][3] + *hyptphi1[6][4];
  hyptphi1[6][4]->Reset();

  *hyptphi1[6][7] = *hyptphi1[6][7] + *hyptphi1[6][8];
  hyptphi1[6][8]->Reset();
  *hyptphi1[6][7] = *hyptphi1[6][7] + *hyptphi1[6][9];
  hyptphi1[6][9]->Reset();

  *hyptphi1[7][4] = *hyptphi1[7][4] + *hyptphi1[7][5];
  hyptphi1[7][5]->Reset();



  //v2 pt
  *hyptphi2[1][1] = *hyptphi2[1][1] + *hyptphi2[1][2];
  hyptphi2[1][2]->Reset();

  *hyptphi2[2][3] = *hyptphi2[2][3] + *hyptphi2[2][4];
  hyptphi2[2][4]->Reset();

  *hyptphi2[3][5] = *hyptphi2[3][5] + *hyptphi2[3][6];
  hyptphi2[3][6]->Reset();

  *hyptphi2[4][2] = *hyptphi2[4][2] + *hyptphi2[4][3];
  hyptphi2[4][3]->Reset();

  *hyptphi2[4][6] = *hyptphi2[4][6] + *hyptphi2[4][7];
  hyptphi2[4][7]->Reset();




  // output    
  //************************************************** 
  std::cout << " ---- Resutls ---------------------" << std::endl;

  TString fName = Form("YPT_n%2dto%d_",cent[icent],cent[jcent]) + sysName[isys[0]] + "_neutron.root";
  gSystem->cd("data");
  auto GraphSave = new TFile(fName,"recreate");
  std::cout << " File " << fName << " is created." << std::endl;


  TGraphErrors *gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  TGraphErrors *gv_v2 = new TGraphErrors();
  gv_v2->SetName("gv_v2");

  TGraphErrors *gPt_v1[ybin1];
  TGraphErrors *gPt_v2[ybin1];

  for(UInt_t kn = 0; kn < ybin1 ; kn++){      
    gPt_v1[kn] = new TGraphErrors();
    gPt_v1[kn]->SetName((TString)Form("gPt_v1%d",kn));
    TString sname = "Neutron; Pt [MeV/c]; v1";
    gPt_v1[kn]->SetTitle(sname);
  }

  for(UInt_t kn = 0; kn < ybin2 ; kn++){      
    gPt_v2[kn] = new TGraphErrors();
    gPt_v2[kn]->SetName((TString)Form("gPt_v2%d",kn));
    TString sname = "Neutron; Pt [MeV/c]; v2";
    gPt_v2[kn]->SetTitle(sname);
  }


  Double_t *rpres = new Double_t[4];
  rpres = GetRPResolutionwChi(hphi0_180, hphi90_180);

  ic++; 
  cc[ic]   = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hyptacp->Draw("colz");

  ic++; 
  cc[ic]   = new TCanvas("dphi1","dphi1",2000,1200);
  cc[ic]->Divide(nbin1, ybin1);

  cc[ic+1] = new TCanvas("dphi2","dphi2",1400,1200);
  cc[ic+1]->Divide(nbin2, ybin2);

  cc[ic+2] = new TCanvas(Form("cc%d",ic+2),Form("cc%d",ic+2), 1000, 1200);
  cc[ic+2]->Divide(2, ybin1);

  UInt_t kl = 0;
  UInt_t id1 = 0;
  UInt_t id2 = 0;
  Double_t para[6];
    
  for(UInt_t kn = 0; kn < ybin1; kn++){
    // v1 vs rapidity
    cc[ic+2]->cd(2*(kn+1)-1);

    if(  hypt[kn]->GetEntries() < 0 ) continue;
      
    Double_t rapm  = hypt[kn]->GetMean(1);
    Double_t rape  = hypt[kn]->GetStdDev(1);


    if( hyphi1[kn]->GetEntries() < 0 ) continue;
    auto hyphi = (TH1D*)hyphi1[kn]->ProjectionX(Form("hyphiv1_%d",kn),0,-1,"eo");

    Double_t corr[2]= {rpres[0], rpres[1]};
    GetFittingParameters(*hyphi, para, corr);
    // hyphi->SetMaximum(1.5);
    // hyphi->SetMinimum(0.5);
      
    if( !std::isnan(para[4]) && !std::isinf(para[4] )) {
      gv_v1->SetPoint( kl,     rapm, para[4]);
      gv_v1->SetPointError( kl, rape, para[5]);
      kl++;
    }

    // pt dependence 
    UInt_t il = 0; 
    for(UInt_t jn = 0; jn < nbin1; jn++){

      id1++; cc[ic]->cd(id1); 

      if( hyptphi1[kn][jn]->GetEntries() > 0 ) {	

	Double_t ptc  = hyptphi1[kn][jn]->GetMean(2);
	Double_t ptce = hyptphi1[kn][jn]->GetStdDev(2);
	    
	auto hypt = (TH1D*)hyptphi1[kn][jn]->ProjectionX();
	GetFittingParameters(*hypt, para, corr);
	    

	if( !std::isnan(para[4]) && !std::isinf(para[4] )) {
	  gPt_v1[kn]->SetPoint(il, ptc, para[4]);
	  gPt_v1[kn]->SetPointError(il, ptce, para[5]);
	    
	  il++;
	}
      }
    }
    gPt_v1[kn]->Write();
  }

  kl = 0;
  for(UInt_t kn = 0; kn < ybin2; kn++){

    // v2 vs rapidity
    cc[ic+2]->cd(2*(kn+1));

    if(  hypt2[kn]->GetEntries() < 0 ) continue;

    Double_t rapm  = hypt2[kn]->GetMean(1);
    Double_t rape  = hypt2[kn]->GetStdDev(1);

    if( hyphi2[kn]->GetEntries() < 0 ) continue;
    auto hyphii = (TH1D*)hyphi2[kn]->ProjectionX(Form("hyphiv2_%d",kn),0,-1,"eo");
      
    Double_t corr[2] = {rpres[2], rpres[3]};
    GetFittingParameters(*hyphii, para, corr);
    //hyphii->SetMaximum(1.1);
    //hyphii->SetMinimum(0.9);

    if( !std::isnan(para[4]) && !std::isinf(para[4] )) {
      gv_v2->SetPoint( kl,    rapm, para[4]);
      gv_v2->SetPointError( kl, rape, para[5]);
      kl++;
    }
    

    UInt_t il = 0; 
    for(UInt_t jn = 0; jn < nbin2; jn++){
	
      id2++; cc[ic+1]->cd(id2); 
      if( hyptphi2[kn][jn]->GetEntries() > 0 ){

	Double_t ptc  = hyptphi2[kn][jn]->GetMean(2);
	Double_t ptce = hyptphi2[kn][jn]->GetStdDev(2);

	auto hypt = (TH1D*)hyptphi2[kn][jn]->ProjectionX();

	GetFittingParameters(*hypt, para, corr);

	if( !std::isnan(para[4]) && !std::isinf(para[4] )) {
	  gPt_v2[kn]->SetPoint(il, ptc, para[4]);
	  gPt_v2[kn]->SetPointError(il, ptce, para[5]);
	    
	  il++;
	}
      }
    }
    gPt_v2[kn]->Write();
  }


  ic+=3;


  //plotting
  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
  cc[ic]->Divide(ybin1,1);
  for(UInt_t kn = 0; kn < ybin1; kn++){
    cc[ic]->cd(id); id++;
    gPt_v1[kn]->Draw("ALP");
  }

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
  cc[ic]->Divide(ybin2,1);
  for(UInt_t kn = 0; kn < ybin2; kn++){
    cc[ic]->cd(id); id++;
    gPt_v2[kn]->Draw("ALP");
  }

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,400);
  cc[ic]->Divide(ybin2,1);
  for(UInt_t kn = 0; kn < ybin2; kn++){
    cc[ic]->cd(id); id++;
    hypt2[kn]->Draw("colz");
  }

  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,400);
  cc[ic]->Divide(ybin1,1);
  for(UInt_t kn = 0; kn < ybin1; kn++){
    cc[ic]->cd(id); id++;
    hypt[kn]->Draw("colz");
  }


  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v1->Draw("ALP");

  auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();


  ic++; id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v2->Draw("ALP");


  gv_v1->Write();
  gv_v2->Write();

  
  gSystem->cd("..");
}

void GetFittingParameters(TH1D &h1, Double_t pp[6])
{
  for(UInt_t i = 0; i < 6; i++)
    pp[i] = 0.;


  Double_t nbin = h1.GetNbinsX();
  Double_t scf  = h1.GetEntries();
  h1.SetNormFactor(nbin);
  // h1.SetMaximum(1.2);
  // h1.SetMinimum(0.8);
  fcos1->SetParameter(0, 1.);
  fcos1->SetParameter(1, 0.1);

  h1.Fit("fcos1","","",-3.1,3.1);

  pp[0] = fcos1->GetParameter(0);
  pp[1] = fcos1->GetParameter(1);
  pp[2] = fcos1->GetParError(0);
  pp[3] = fcos1->GetParError(1);
  pp[3] = pp[1]*sqrt(pp[2]/pp[0]*pp[2]/pp[0] + pp[3]/pp[1]*pp[3]/pp[1]);
  pp[1] = pp[1]/pp[0];
  
  

  h1.SetNormFactor(-1);
}
void GetFittingParameters(TH1D &h1, Double_t pp[6], Double_t corr[2])
{
  GetFittingParameters(h1,pp);
  
  if( corr[0] != 0 && corr[1] != 0 ) {
    pp[4] = pp[1]/corr[0];
    pp[5] = sqrt( pow(pp[4],2)*( pow(pp[3]/pp[1],2) + pow(corr[1]/corr[0], 2) ));
  }
  else {
    pp[4] = 0.;
    pp[5] = 0.;
  }
}


void PlotISO()                   
{
  //----- Parametres

  //----- Booking
  for(Int_t m = m_bgn; m < m_end; m++){
    auto hisobmx   = new TH1D(Form("hisobmx_%d",m), "ISO Longitudianl at Target", 100,0., 20000);
    auto hisobmy   = new TH1D(Form("hisobmy_%d",m), "ISO Longitudianl at Target", 100,0., 20000);
    auto hisobmphi = new TH1D(Form("hisobmphi_%d",m), "ISO Longitudianl at Target", 100,0., 20000);
  }

  //----- Filling
  for(Int_t m = m_bgn; m < m_end; m++){
    Int_t nEntry = rChain[m]->GetEntries();

    for(Int_t i = 0; i < nEntry; i++){
      rChain[m]->GetEntry(i);
    }
  }
  //----- Drawing 
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,1000);
  cc[ic]->Divide(3,m_end);

  UInt_t id = 1;
  for(Int_t m = m_bgn; m < m_end; m++){
  }
}


void DetectorBias()
{
  UInt_t m = 0;

  SetBranch(m);
  Int_t nEntry = rChain[m]->GetEntries();


  Double_t  cosphi = 0.;
  Double_t  sinphi = 0.;
  Double_t  cos2phi= 0.;
  Double_t  sin2phi= 0.;

  Double_t count = 0.;
  for(Int_t i = 0; i < nEntry; i++){
    rChain[m]->GetEntry(i);

    if(ntrack[4] > cent[icent] || ntrack[4] <= cent[jcent] ) continue;      
      
    count++;

    cosphi += cos(unitP_fc->Phi());      
    sinphi += sin(unitP_fc->Phi());

    cos2phi += cos(2.*unitP_fc->Phi());      
    sin2phi += sin(2.*unitP_fc->Phi());

  }
  
  std::cout << " <cos Phi> " << cosphi/count 
	    << " <sin phi> " << sinphi/count << std::endl;

  std::cout << " <cos 2Phi> " << cos2phi/count 
	    << " <sin 2phi> " << sin2phi/count << std::endl;


}

void CentralityDependence()            //%% Executable :
{
  //  UInt_t mrange[] = {70, 50, 45, 40, 35, 25, 20, 0};
  UInt_t mrange[] = {70, 35, 28, 0};

  const UInt_t mbin = sizeof(mrange)/sizeof(UInt_t) - 1;
  TH1D *hphi0_180;
  TH1D *hphi90_180;
  TH1I *hmult;

  auto gv_mcos1 = new TGraphErrors();
  gv_mcos1->SetTitle(";Multiplicity; <cos(#Psi)>");
  auto gv_mcos2 = new TGraphErrors();
  gv_mcos2->SetTitle(";Multiplicity; <cos(2#Psi)>");


  auto cc80 = new TCanvas("cc80","cc80",500,1000);
  cc80->Divide(2,7);
  
  UInt_t id = 1;
  for(UInt_t i = 0; i < mbin; i++) {
    cc80->cd(id); id++;
    TCut hcut = Form("ntrack[4]>%u&&ntrack[4]<=%u",mrange[i+1],mrange[i]);

    TString htitle = Form("hphi0_%d",i);
    hphi0_180  = new TH1D(htitle, "",100,0.,3.2);
    rChain[0]->Project(htitle,"abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi()))",hcut);
    hphi0_180->Draw();

    htitle = Form("hphi90_%d",i);
    hphi90_180 = new TH1D(htitle,"",100,0.,3.2);
    TCut phicut = "abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi()))>1.5707963";
    rChain[0]->Project(htitle,"abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi()))",hcut&&phicut);
    hphi90_180->Draw("same");

    hmult = new TH1I(Form("hmult_%u",i),"",70,0.,70.);
    rChain[0]->Project(Form("hmult_%u",i),"ntrack[4]",hcut);
    cc80->cd(id); id++;
    hmult->Draw();
    

     Double_t *rpres = new Double_t[4];
     rpres = GetRPResolutionwChi(hphi0_180, hphi90_180);
     std::cout << mrange[i+1] << " ~ " << mrange[i] 
     	      << " <cos(Phi)> = " << rpres[0] << " +- " << rpres[1] 
     	      << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
     	      << std::endl;
     
     hmult->GetMean();

     gv_mcos1->SetPoint(i, hmult->GetMean(), rpres[0]);
     gv_mcos1->SetPointError(i, hmult->GetStdDev(), rpres[1]);

     gv_mcos2->SetPoint(i, hmult->GetMean(), rpres[2]);
     gv_mcos2->SetPointError(i, hmult->GetStdDev(), rpres[3]);
  }

  hmult = new TH1I("hmult",";Multiplicity",70,0,70);
  rChain[0]->Project("hmult","ntrack[4]");


  auto cc79 = new TCanvas("cc79","cc79");
  hmult->Draw("");

  auto cc81 = new TCanvas("cc81","cc81");
  gv_mcos1->Draw("ALP");

  auto cc82 = new TCanvas("cc82","cc82");
  gv_mcos2->Draw("ALP");
}



//--------------------------------------------//%% Executable : 
void Template()                   
{
  //----- Parametres

  //----- Booking
  for(Int_t m = m_bgn; m < m_end; m++){

  }

  //----- Filling
  for(Int_t m = m_bgn; m < m_end; m++){
    Int_t nEntry = rChain[m]->GetEntries();

    for(Int_t i = 0; i < nEntry; i++){
      rChain[m]->GetEntry(i);
    }
  }
  //----- Drawing 
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,1000);
  cc[ic]->Divide(3,m_end);

  UInt_t id = 1;
  for(Int_t m = m_bgn; m < m_end; m++){
  }
}


UInt_t SetBranch(UInt_t m=0)
{

  if(aArray != NULL)
    aArray->Clear();

  if(rChain[m] == NULL) {
    std::cout << " no file is loaded " << std::endl;
    return 0;
  }

  rChain[m]->SetBranchAddress("STParticle",&aArray);
  rChain[m]->SetBranchAddress("ntrack",   ntrack);
  rChain[m]->SetBranchAddress("aoq",&aoq);
  rChain[m]->SetBranchAddress("z",&z);
  rChain[m]->SetBranchAddress("snbm",&snbm);
  rChain[m]->SetBranchAddress("ProjA",&ProjA);
  rChain[m]->SetBranchAddress("ProjB",&ProjB);
  rChain[m]->SetBranchAddress("unitP_fc"  ,&unitP_fc,&bunitP_fc);
  rChain[m]->SetBranchAddress("unitP_rc"  ,&unitP_rc,&bunitP_rc);
  rChain[m]->SetBranchAddress("unitP_1"   ,&unitP_1,&bunitP_1);
  rChain[m]->SetBranchAddress("unitP_2"   ,&unitP_2,&bunitP_2);
  rChain[m]->SetBranchAddress("mtrack_1"  ,&mtrack_1);    
  rChain[m]->SetBranchAddress("mtrack_2"  ,&mtrack_2);
  rChain[m]->SetBranchAddress("unitP_lang",&unitP_lang,&bunitP_lang);
  rChain[m]->SetBranchAddress("eisott"    ,&eisott);
  rChain[m]->SetBranchAddress("eisobm"    ,&eisobm);
  rChain[m]->SetBranchAddress("eisotg"    ,&eisotg);
  rChain[m]->SetBranchAddress("STNeuLANDCluster", &aNLClusterArray);

  
  SetupDB();

  return rChain[m]->GetEntries();

}

