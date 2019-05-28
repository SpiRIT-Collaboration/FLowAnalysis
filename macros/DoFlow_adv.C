#include "openRunAna.C"
#include "DoFlow.h"

auto *fcos1 = new TF1("fcos1","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());

//drawing

UInt_t selReactionPlanef = 10000;

// Retrivew tree
TClonesArray *aArray; 
TClonesArray *aFlowArray;
TClonesArray *aNLClusterArray;

ROOT::Math::Interpolator *itrpvx;

std::vector< Double_t > v1x;
std::vector< Double_t > v1y;
std::vector< Double_t > v1xe;
std::vector< Double_t > v1ye;

std::vector< Double_t > v2x;
std::vector< Double_t > v2y;
std::vector< Double_t > v2xe;
std::vector< Double_t > v2ye;



// functions
void     GetFittingParameters(TH1D &h1, Double_t pp[6]);
void     GetFittingParameters(TH1D &h1, Double_t pp[6], Double_t corr[2]);
Double_t GetRapidity_cm(TVector3 p, Double_t mass, TVector3 bvec);
UInt_t   SetBranch();
void     PlotCosPtDependence(UInt_t selid);
TString  SetupOutputFile(TString fopt);

Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180);
Bool_t   SetRPResolution();
Double_t *GetRPResolution(UInt_t vn, UInt_t mult);
UInt_t   GetV1RapidityIndex(Double_t val);
UInt_t   GetV2RapidityIndex(Double_t val);
UInt_t   GetV1PtIndex(Double_t val);
UInt_t   GetV2PtIndex(Double_t val);
UInt_t    GetMultiplicityIndex(UInt_t mult);

Double_t  GetError(Double_t x, Double_t y, Double_t xe, Double_t ye)
{
  if( y == 0 && ye == 0 ) return 0.;
  Double_t xc = abs(x/y);
  return  xc * sqrt(pow(xe/x,2) + pow(ye/y,2));
}

//Double_t rapoffset[] = {0., 0.018, -0.008, 0.028};
Double_t rapoffset[] = {0.01113, -0.00645, 0.01898, -0.01726};
UInt_t   npb = 30;


//-------------------//
void DoFlow_adv(UInt_t isel = 2) 
{
  gROOT->Reset();

  openRunAna();

  if(rChain != NULL)     
    std::cout << " System " << isys << "  -> " << sysName << std::endl; 

  else
    exit(0);


  gROOT->ProcessLine(".! grep -i void DoFlow.C ");


  // Multiplicity cut ================================
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

  std::cout << " Multiplicity :: " << Lcent << " to " << Ucent << std::endl;
  //==================================================

  //  PlotCosPtDependence(isel);  


}


//pdt
void PlotCosPtDependence(UInt_t selid = 2)       //%% Executable :
{
  gStyle->SetOptStat(0);

  if(SetRPResolution())
    std::cout << " RP resolution was ready " << std::endl;

  std::cout << "PlotCosPtDependence(" << selid << ")" << std::endl;
  dpt1 = pt_max/(Double_t)ptbin1;
  dpt2 = pt_max/(Double_t)ptbin2;

  // auto cutfile = new TFile("db/RegionCut.root");
  // TCutG *goodThetaPhi = (TCutG*)cutfile->Get("goodThetaPhi");
  // cutfile->Close();    


  TString fHeader = "advYPt_"+ sysName + "_" + partname[selid]+".v"+sVer+".";
  auto fName = SetupOutputFile( fHeader );

  auto GraphSave  = new TFile(fName,"recreate");

  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  std::cout << " Rapidity binning " << ybin1 << std::endl;
  TString rangeLabel1[ybin1];
  for(UInt_t i = 0; i < ybin1; i++ ){
    if( i == 0 )
      rangeLabel1[0] = Form(" y < %5.2f ",yrange1[0]);
    else if ( i == ybin1 - 2 )
      rangeLabel1[i] = Form("%5.2f <= y "    ,yrange1[i]);
    else 
      rangeLabel1[i] = Form("%5.2f <= y < %5.2f",yrange1[i-1],yrange1[i]);
  }
  TString rangeLabel2[ybin2];
  for(UInt_t i = 0; i < ybin2; i++ ){
    if( i == 0 )
      rangeLabel2[0] = Form(" y < %5.2f ",yrange2[0]);
    else if ( i == ybin2 - 2 )
      rangeLabel2[i] = Form("%5.2f <= y "    ,yrange2[i]);
    else 
      rangeLabel2[i] = Form("%5.2f <= y < %5.2f",yrange2[i-1],yrange2[i]);
  }

  //------------------------------
  std::cout << " v1 BIN y " << ybin1 << " pt " << ptbin1 << " mult " << mbin << std::endl;
  std::cout << " v2 BIN y " << ybin2 << " pt " << ptbin2 << " mult " << mbin << std::endl;

  
  TH1I *hmult = new TH1I("hmult",";Multiplicity",80,0,80);
  TH1F *hdydptcos1[ybin1][ptbin1][mbin];
  TH1F *hdydptcos2[ybin1][ptbin1][mbin];
  TH1D *hdydpt1[ybin1][ptbin1];
  TH1D *hdy1[ybin1];
  TH1D *hdydpt2[ybin1][ptbin1];
  TH1D *hdy2[ybin1];

  //-----  booking
  for( UInt_t iy = 0; iy < ybin1; iy++ ) {
    hdy1[iy] = new TH1D(Form("hdy1_y%d",iy),rangeLabel1[iy]+";Rapidity", 500,-0.5,0.8);

    for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {
      hdydpt1[iy][ipt] = new TH1D(Form("hdydpt1_y%dpt%d",iy,ipt),rangeLabel1[iy]+";Pt", 500,0., 1300);

      for( UInt_t im = 0; im < mbin; im++ ) 
	hdydptcos1[iy][ipt][im] = new TH1F(Form("hdydptcos1_y%dpt%dm%d",iy,ipt,im),"; cos(#phi - #Psi)" , 100., -1., 1.);
    }
  }

  for( UInt_t iy = 0; iy < ybin2; iy++ ) {
    hdy2[iy] = new TH1D(Form("hdy2_y%d",iy),rangeLabel2[iy]+";Rapidity", 500,-0.5,0.8);

    for( UInt_t ipt = 0; ipt < ptbin2; ipt++ ) {
      hdydpt2[iy][ipt] = new TH1D(Form("hdydpt2_y%dpt%d",iy,ipt),rangeLabel2[iy]+"; Pt", 500,0., 1300);

      for( UInt_t im = 0; im < mbin; im++ ) {
	hdydptcos2[iy][ipt][im] = new TH1F(Form("hdydptcos2_y%dpt%dm%d",iy,ipt,im),"; cos2(#phi - #Psi)", 100.,-1.,1.);
      }
    }
  }

  //------------------------------


  Int_t nevt = SetBranch();
  
  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  for(Int_t i = 0; i < nevt; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    /// for reaction plane resolution
    Bool_t bFill = kFALSE;


    /// centrality selection
    if(aflow->mtrack4 > Ucent || aflow->mtrack4 <= Lcent ) continue;

    hmult->Fill( aflow->mtrack4 );

    TIter next(aArray);
    STParticle *aPart = NULL;

    //--------------------------------------------------
    //----- Main loop 
    //--------------------------------------------------
    while( (aPart = (STParticle*)next()) ) {
      
      auto rpf   = aPart->GetReactionPlaneFlag();
      auto pid   = aPart->GetPID();
      auto yaw   = aPart->GetYawAngle();
      auto pitch = aPart->GetPitchAngle();
      auto ndf   = aPart->GetNDF();


      //------------------------------
      //----- Particle Selection -----
      //      if( pid == partid[selid] && rpf > 10000 ){ //&& yaw < 0 && pitch < 0) {
      if( (selid == 8 && (pid == partid[2] || pid == partid[3] || pid == partid[4]) &&
	   yaw > 0) ||
	  (pid == partid[selid]  && yaw > 0 ) ) {
      
	bFill = kTRUE;
      
	auto pt    = aPart->GetRotatedMomentum().Pt();
	auto dphi  = aPart->GetAzmAngle_wrt_RP();
	auto rapid = aPart->GetRapiditycm();;	
	auto bmass = aPart->GetBBMass();
	auto phi   = aPart->GetRotatedMomentum().Phi();
	auto theta = aPart->GetRotatedMomentum().Theta();

	UInt_t irapid1 = GetV1RapidityIndex(rapid);
	UInt_t ipt1    = GetV1PtIndex(pt);

	UInt_t irapid2 = GetV2RapidityIndex(rapid);
	UInt_t ipt2    = GetV2PtIndex(pt);
	UInt_t im      = GetMultiplicityIndex(aflow->mtrack4);
	

	//	cout << " im " << im << " " << ipt2 << " " << irapid2 << " " << rapid << " " << pt << " " << dphi <<  endl;

	hdy1[irapid1]->Fill( rapid );
	hdy2[irapid2]->Fill( rapid );

	hdydpt1[irapid1][ipt1] -> Fill( pt );
	hdydpt2[irapid2][ipt2] -> Fill( pt );

	hdydptcos1[irapid1][ipt1][im]->Fill( cos(dphi) );
	hdydptcos2[irapid2][ipt2][im]->Fill( cos(2.*dphi) );

      }
    }
  }
  //--------------------------------------------------
  //--- Enf of event Loop
  //--------------------------------------------------
  std::cout << " End of Event Loop " << std::endl;

  TGraphErrors *gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  TGraphErrors *gv_v2 = new TGraphErrors();
  gv_v2->SetName("gv_v2");

  TGraphErrors *gPt_v1[ybin1];
  TGraphErrors *gPt_v2[ybin2];

  for(UInt_t kn = 0; kn < ybin1 ; kn++){
    gPt_v1[kn] = new TGraphErrors();
    gPt_v1[kn]->SetName((TString)Form("gPt_v1%d",kn));
    TString sname = partname[selid]+" "+rangeLabel1[kn]+"; Pt [MeV/c]; v1";
    gPt_v1[kn]->SetTitle(sname);
  }

  for(UInt_t kn = 0; kn < ybin2 ; kn++){
    gPt_v2[kn] = new TGraphErrors();
    gPt_v2[kn]->SetName((TString)Form("gPt_v2%d",kn));
    TString sname = partname[selid]+" "+rangeLabel2[kn]+"; Pt [MeV/c]; v2";
    gPt_v2[kn]->SetTitle(sname);
  }

  //-------------------- v1
  UInt_t iyv = 0;
  for( UInt_t iy = 0; iy < ybin1; iy++ ) {
    Double_t v_y   = 0.;
    Double_t v_ye2 = 0.;

    UInt_t iptv = 0;
    for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {
      
      Double_t ptn    = (Double_t)hdydpt1[iy][ipt]->GetEntries();
      Double_t ptmean = hdydpt1[iy][ipt]->GetMean();
      Double_t pte    = hdydpt1[iy][ipt]->GetStdDev()/sqrt(ptn);

      Double_t v_pt    = 0;
      Double_t vn_pt   = 0;
      Double_t v_pte2  = 0.;
      Double_t imm     = 0.;
      for( UInt_t im = 0; im < mbin; im++ ){
	
	Double_t *rpres = new Double_t[3];
	rpres = GetRPResolution(1, im);

	Double_t vn = (Double_t)hdydptcos1[iy][ipt][im]->GetEntries();
	Double_t vm = hdydptcos1[iy][ipt][im]->GetMean();
	Double_t ve = hdydptcos1[iy][ipt][im]->GetStdDev()/sqrt(vn);

	if( rpres[1] > 0 && vn > 0 ) {
	  Double_t vc = vm/rpres[1] * vn;
	  v_pt  += vc;
	  vn_pt += vn;

	  Double_t vce2   = pow(vm/rpres[1], 2) * ( pow(ve/vm,2) + pow(rpres[2]/rpres[1],2) );

	  v_pte2  += vce2*pow(vn,2);
	}
      }

      Double_t v1  = v_pt / vn_pt;
      Double_t v1e = sqrt(v_pte2)/vn_pt;

      if( !std::isnan(ptmean) && !std::isnan(pte) ) {
	gPt_v1[iy]->SetPoint(iptv, ptmean, v1);
	gPt_v1[iy]->SetPointError(iptv, pte, v1e);
	iptv++;

	v_y   += v1*ptn;
	v_ye2 += pow(v1e*ptn,2); 

      }
    }

    Double_t yn    = (Double_t)hdy1[iy]->GetEntries();
    Double_t v1_y  = v_y / yn;
    Double_t v1_ye = sqrt(v_ye2)/yn;
    
    if( !std::isnan(v1_y) && !std::isnan(v1_ye) ) {
      Double_t ymean = hdy1[iy]->GetMean();
      Double_t ystdv = hdy1[iy]->GetStdDev();
    
      gv_v1->SetPoint( iyv, ymean, v1_y);
      gv_v1->SetPointError( iyv, ystdv, v1_ye);
      iyv++;
    }

    gPt_v1[iy]->Write();
    //    gv_v1->SetPoint(iv1y, ymean, 0.);
  }

  //--------------------v2
  iyv = 0;
  for( UInt_t iy = 0; iy < ybin2; iy++ ) {
    Double_t v_y   = 0.;
    Double_t v_ye2 = 0.;

    UInt_t iptv = 0;

    for( UInt_t ipt = 0; ipt < ptbin2; ipt++ ){
    
      Double_t ptn    = (Double_t)hdydpt2[iy][ipt]->GetEntries();
      Double_t ptmean = hdydpt2[iy][ipt]->GetMean();
      Double_t pte    = hdydpt2[iy][ipt]->GetStdDev()/sqrt(ptn);

      Double_t v_pt    = 0;
      Double_t vn_pt   = 0;
      Double_t v_pte2  = 0.;

      for( UInt_t im = 0; im < mbin; im++ ){
	
	Double_t *rpres = new Double_t[3];
	rpres = GetRPResolution(2, im);
	
	Double_t vn = (Double_t)hdydptcos2[iy][ipt][im]->GetEntries();
	Double_t vm = hdydptcos2[iy][ipt][im]->GetMean();
	Double_t ve = hdydptcos2[iy][ipt][im]->GetStdDev()/sqrt(vn);
	
	if( rpres[1] > 0 && vn > 0 ) {
	  Double_t vc = vm/rpres[1] * vn;
	  v_pt  += vc;
	  vn_pt += vn;

	  Double_t vce2   = pow(vm/rpres[1], 2) * ( pow(ve/vm,2) + pow(rpres[2]/rpres[1],2) );
	  v_pte2  += vce2*pow(vn,2);
	}
      }

      Double_t v2  = v_pt / vn_pt;
      Double_t v2e = sqrt(v_pte2)/vn_pt;

      if( !std::isnan(ptmean) && !std::isnan(pte) ) {
	gPt_v2[iy]->SetPoint(iptv, ptmean, v2);
	gPt_v2[iy]->SetPointError(iptv, pte, v2e);
	iptv++;


	v_y   += v2*ptn;
        v_ye2 += pow(v2e*ptn,2);

      }
    }

    Double_t yn    = (Double_t)hdy2[iy]->GetEntries();
    Double_t v2_y  = v_y / yn;
    Double_t v2_ye = sqrt(v_ye2)/yn;


    if( !std::isnan(v2_y) && !std::isnan(v2_ye) ) {
      Double_t ymean = hdy2[iy]->GetMean();
      Double_t ystdv = hdy2[iy]->GetStdDev();
      
      gv_v2->SetPoint( iyv, ymean, v2_y);
      gv_v2->SetPointError( iyv, ystdv, v2_ye);
      iyv++;
    }

    gPt_v2[iy]->Write();
  }

  gv_v1->Write();
  gv_v2->Write();
  hmult->Write();
      
  //--------------------------------------------------
  //--- Ploting
  //--------------------------------------------------
  id=1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  cc->Divide(10);

  for( UInt_t iy = 0; iy < 10; iy++ ) {
    //    for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {

      cc->cd(id); id++;
      gPt_v1[iy] -> SetMarkerStyle(20);
      gPt_v1[iy] -> SetMarkerColor(4);
      gPt_v1[iy] -> Draw("ALP");

      //    }
  }

  id=1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  cc->Divide(7);

  for( UInt_t iy = 0; iy < 7; iy++ ) {
    //    for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {

      cc->cd(id); id++;
      gPt_v2[iy] -> SetMarkerStyle(20);
      gPt_v2[iy] -> SetMarkerColor(4);
      gPt_v2[iy] -> Draw("ALP");

      //    }
  }


  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v1->Draw("ALP");

  auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v2->Draw("ALP");



}


UInt_t GetV1RapidityIndex(Double_t val)
{
  UInt_t irapid = ybin1 - 2;
  for( UInt_t i = 0; i < ybin1; i++){
    if( val < yrange1[i]){
      irapid = i;
      break;
    }
  }
  return irapid;
}

UInt_t GetV2RapidityIndex(Double_t val)
{
  UInt_t irapid = ybin2 - 2;
  for( UInt_t i = 0; i < ybin2; i++){
    if( val < yrange2[i]){
      irapid = i;
      break;
    }
  }
  return irapid;
}

UInt_t GetV1PtIndex(Double_t val)
{
  UInt_t ipt = ptbin1 - 1;
  for(UInt_t i = 0; i < ptbin1; i++){
    if( val < dpt1*(i+1)) {
      ipt = i;
      break;
    }
  }
  return ipt;
}
UInt_t GetV2PtIndex(Double_t val)
{
  UInt_t ipt = ptbin2 - 1;
  for(UInt_t i = 0; i < ptbin2; i++){
    if( val < dpt2*(i+1)) {
      ipt = i;
      break;
    }
  }
  return ipt;
}




//**************************************************
UInt_t GetMultiplicityIndex(UInt_t mult)
{
  UInt_t xn = v1x.size() - 1;
  if( mult <= v1x.at( v1x.size() - 1 ) ) {
    xn = (UInt_t)itrpvx->Eval( (Double_t)mult ) ;
  
    if( abs(v1x.at(xn) - mult) > abs(v1x.at(xn+1) - mult) )
      xn += 1;
  }
  return xn;
}

Double_t *GetRPResolution(UInt_t vn, UInt_t xn)
{

  Double_t *rpcor = new Double_t[3];
  rpcor[0] = 0.;
  rpcor[1] = 0.;
  rpcor[2] = 0.;

  if( xn < v1x.size() ) {

    if( vn == 1 ){

      rpcor[0] = v1x.at(xn);   // multiplicity
      rpcor[1] = v1y.at(xn);   // v1 corr 
      rpcor[2] = v1ye.at(xn);  // v1 corr error
    }
    else {
      rpcor[0] = v2x.at(xn);
      rpcor[1] = v2y.at(xn);
      rpcor[2] = v2ye.at(xn);
    }
  }
  return rpcor;

}
//**************************************************
Bool_t SetRPResolution()
{
  itrpvx = new ROOT::Math::Interpolator(20, ROOT::Math::Interpolation::kPOLYNOMIAL);

  TString fname = "data/mlt_"+ sysName + ".v25.0.5.root";
  TFile *fOpen = TFile::Open(fname);
  if( fOpen == NULL ) kFALSE;

  auto hgv_mcos1 = (TGraphErrors*)fOpen->Get("gv_mcos1");
  auto hgv_mcos2 = (TGraphErrors*)fOpen->Get("gv_mcos2");

  std::vector< Double_t > ix;
    
  Double_t x, y, xe, ye;
  UInt_t k = 0;
  ix.push_back(k); 
  v1x.push_back(0.); v1xe.push_back(0.);
  v1y.push_back(0.); v1ye.push_back(0.);
  k++;
  for( UInt_t i = (UInt_t)hgv_mcos1->GetN()-1; i > 0; i-- ) {
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
  itrpvx->SetData(v1x,ix);


  v2x.push_back(0.); v2xe.push_back(0.);
  v2y.push_back(0.); v2ye.push_back(0.);

  k = 0;
  for( UInt_t i = (UInt_t)hgv_mcos2->GetN()-1; i > 0; i-- ) {
    hgv_mcos2->GetPoint(i, x, y);
    xe = hgv_mcos2->GetErrorX(i);
    ye = hgv_mcos2->GetErrorY(i);

    v2x.push_back(x);
    v2y.push_back(y);
    v2xe.push_back(xe);
    v2ye.push_back(ye);

    cout << k << " th " << v2x.at(k) << " vs " << v2y.at(k) << " +- " << v2ye.at(k) << endl; 
    k++;
  }

  fOpen->Close();

  return kTRUE;
}

//**************************************************
Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180)            //%% Executable : 
{
  hphi90_180->SetLineColor(2);

  Double_t m0 = hphi0_180->GetEntries();
  Double_t m1 = hphi90_180->GetEntries();

  if( m0 == 0 ) {
    std::cout << " No entry in hphi0_180 " << endl;
    return 0;
  }

  //  Double_t chi = sqrt(-2.* log(2.* m1/m0));
  Double_t mr    = m1/m0;
  Double_t mre2  = m1/m0 * (1./m0 + 1./m1);
  Double_t chi   = sqrt(-4.* log(2.* mr));
  Double_t chie2 = -4.*log(2.*mr)/pow(mr,2)* mre2;

  Double_t *rpres = new Double_t[4];
  //v1
  if( chi > 0. ) {

    //v1
    rpres[0] = 0.626657*chi - 0.09694*pow(chi,3) + 0.02754 * pow(chi,4) - 0.002283*pow(chi,5);
  
    Double_t err2 = pow(0.626657,2) * chie2;
    err2 += pow(0.09694,2) *  6.* pow(chi,4) * chie2;
    err2 += pow(0.02754,2) * 16.* pow(chi,6) * chie2;
    err2 += pow(0.002283,2)* 25.* pow(chi,8) * chie2;
    rpres[1] = sqrt( err2 );

    //v2
    rpres[2] = 0.25*pow(chi,2) - 0.011414*pow(chi,3) - 0.034726*pow(chi,4) + 0.006815*pow(chi,5);
    err2  = pow(0.25,2) * chie2;
    err2 += pow(0.011414,2) *  9.* pow(chi,4) * chie2;
    err2 += pow(0.034726,2) * 16.* pow(chi,6) * chie2;
    err2 += pow(0.006815,2) * 25.* pow(chi,5) * chie2;
    rpres[3] = sqrt( err2 );
  }

  else {
    Double_t chie = sqrt(chie2);
    rpres[0] = sqrt(TMath::Pi())/(2.*sqrt(2))*chi*exp(-chi*chi/4.)*
      (ROOT::Math::cyl_bessel_i(0,chi*chi/4.)+ROOT::Math::cyl_bessel_i(1,chi*chi/4.));
    rpres[1] = 0.;
    rpres[2] = sqrt(TMath::Pi())/(2.*sqrt(2))*chi*exp(-chi*chi/4.)*
      (ROOT::Math::cyl_bessel_i(0.5,chi*chi/4.)+ROOT::Math::cyl_bessel_i(1.5,chi*chi/4.));
    rpres[3] = 0.;
  }

  std::cout << " Resolution : v1 " 
	    << std::setw(14) << rpres[0] << " +- " << std::setw(10) << rpres[1]
	    << " v2 "
	    << std::setw(14) << rpres[2] << " +- " << std::setw(10) << rpres[3]
	    << std::endl;

  return rpres;
}


void FlatteningCheck()            
{
  //----- Parametres
  Int_t ntrack[7];

  //----- Booking
  TH2D *hphitheta;
  TH2D *hphimtrck;

  hphitheta = new TH2D("hphitheta", sysName+"; #Theta ; #Phi",100,0,0.8, 100,-3.2, 3.2);
  hphimtrck = new TH2D("hphimtrck", sysName+"; Number of Track ; #Phi",60,0,60, 100,-3.2, 3.2);

  //----- Filling
  Int_t nEntry = SetBranch();

  for(Int_t i = 0; i < nEntry; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);
    if( aflow == NULL ) continue;
    

    TIter next(aArray);
    STParticle *aPart = NULL;
    
    while( (aPart = (STParticle*)next()) ) {

      auto phiv  = aPart->GetIndividualRPVector();
      auto theta = aPart->GetRotatedMomentum().Theta();
      auto flag  = aPart->GetReactionPlaneFlag();

      if(flag >= selReactionPlanef ){
	hphitheta->Fill( theta, phiv.Phi() );
	hphimtrck->Fill( aflow->mtrack4, phiv.Phi() ); 
      }
    }
  }


  //----- Drawing 
  ic++;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hphitheta->Draw("colz");

  ic++;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hphimtrck->Draw("colz");

}



void PlotSubEvent(Double_t ml=30, Double_t mu=80)   
{

  Double_t mlt[] = {ml, mu};
  //    Double_t mlt[2] = {0., 8.};
  //Double_t mlt[2] = {8., 16.};
  // Double_t mlt[2] = {16., 24.};
  // Double_t mlt[2] = {24., 32.};
  // Double_t mlt[2] = {32., 40.};
  // Double_t mlt[2] = {40., 100.};

  TCut mcrot = Form("aFlow->mtrack4>%f&&aFlow->mtrack4<%f",mlt[0]*2.,mlt[1]*2.);
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
 

  rChain->Project("hrotx","unitP2_fc.X()",mcrot);
  rChain->Project("h1rx" ,"unitP_1fc.X()"  ,mc1r);
  rChain->Project("h2rx" ,"unitP_2fc.X()"  ,mc2r);
					      
  rChain->Project("hroty","unitP2_fc.Y()",mcrot);
  rChain->Project("h1ry" ,"unitP_1fc.Y()"  ,mc1r);
  rChain->Project("h2ry" ,"unitP_2fc.Y()"  ,mc2r);


  //----- Drawing       
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  hrotx->SetLineColor(2);
  hrotx->SetNormFactor(1);
  h1rx ->SetLineColor(4);
  h1rx ->SetNormFactor(1);
  h2rx ->SetLineColor(6);
  h2rx ->SetNormFactor(1);

  h1rx->Draw();
  h2rx->Draw("same");
  hrotx->Draw("same");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

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

void DetectorBias()
{
  UInt_t m = 0;

  SetBranch();
  Int_t nEntry = rChain->GetEntries();


  Double_t  cosphi = 0.;
  Double_t  sinphi = 0.;
  Double_t  cos2phi= 0.;
  Double_t  sin2phi= 0.;

  Double_t count = 0.;
  for(Int_t i = 0; i < nEntry; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    
    if(aflow->mtrack4 > Ucent || aflow->mtrack4 <= Lcent ) continue;
    count++;

    cosphi += cos( (aflow->unitP_fc).Phi());      
    sinphi += sin( (aflow->unitP_fc).Phi());

    cos2phi += cos(2.*(aflow->unitP_fc).Phi());      
    sin2phi += sin(2.*(aflow->unitP_fc).Phi());

  }
  
  std::cout << " <cos Phi> " << cosphi/count 
	    << " <sin phi> " << sinphi/count << std::endl;

  std::cout << " <cos 2Phi> " << cos2phi/count 
	    << " <sin 2phi> " << sin2phi/count << std::endl;


}

void CentralityDependence()            //%% Executable :
{

  TString fHeader = "mlt_"+ sysName + ".v"+sVer+".";
  auto fName = SetupOutputFile( fHeader );

  auto GraphSave = new TFile(fName,"recreate");

  
  TH1I *hmult = new TH1I("hmult","multiplicity",80,0,80);
  TH1I *hmultbin[mbin];
  TH1D *hphi0_180[mbin];
  TH1D *hphi90_180[mbin];

  for(UInt_t k = 0; k < mbin; k++){
    TString htitle = Form("hmultbin_%d",k);
    hmultbin[k]  = new TH1I(htitle,"",80,0,80);

    htitle = Form("hphi0_180_%d",k);
    hphi0_180[k]  = new TH1D(htitle, "",100,0.,3.2);
    htitle = Form("hphi90_180_%d",k);
    hphi90_180[k] = new TH1D(htitle,"",100,0.,3.2);
  }


  Long64_t nEntry = SetBranch();

  for(Long64_t i = 0; i < nEntry; i++) {

    rChain->GetEntry(i);

    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);
    if( aflow == NULL || aflow->mtrack4 == 0 ) continue;
    hmult -> Fill( aflow->mtrack4 );

    UInt_t ik = 0;
    while( ik < mbin ){
      if( aflow->mtrack4 > mrange[ik] ) break;
      ik++;
    }
    
    //    cout << ik << " " << aflow->mtrack4 << endl;

    hmultbin[ik] -> Fill( aflow->mtrack4 );

    Double_t subdphi = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi()));
    //    Double_t subdphi = abs(TVector2::Phi_mpi_pi((aflow->unitP_1).Phi() - (aflow->unitP_2).Phi()));
    hphi0_180[ik]->Fill( subdphi );
    if( subdphi > TMath::Pi()/2. )
      hphi90_180[ik]->Fill( subdphi );
  }
   

  auto gv_mcos1 = new TGraphErrors();
  gv_mcos1->SetName("gv_mcos1");
  gv_mcos1->SetTitle(";Multiplicity; <cos(#Psi)>");

  auto gv_mcos2 = new TGraphErrors();
  gv_mcos2->SetName("gv_mcos2");
  gv_mcos2->SetTitle(";Multiplicity; <cos(2#Psi)>");

  auto cc80 = new TCanvas("cc80","cc80",500,1000);
  cc80->Divide(2,mbin);
   
  UInt_t id = 1;
  UInt_t ip = 0;
  for(UInt_t i = 0; i < mbin; i++) {
    cc80->cd(id); id++;

    if( hphi0_180[i]->GetEntries() < 5 ) continue;

    hphi0_180[i] ->Draw();
    hphi90_180[i]->Draw("same");

    cc80->cd(id); id++;
    hmultbin[i]->Draw();
    

    Double_t *rpres = new Double_t[4];
    rpres = GetRPResolutionwChi(hphi0_180[i], hphi90_180[i]);
    std::cout << mrange[i+1] << " ~ " << mrange[i] 
     	      << " <cos(Phi)> = " << rpres[0] << " +- " << rpres[1] 
     	      << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
     	      << std::endl;
     
    if( hphi90_180[i]->GetEntries() > 15 && !std::isnan(rpres[0])) {
      gv_mcos1->SetPoint(ip, hmultbin[i]->GetMean(), rpres[0]);
      gv_mcos1->SetPointError(ip, hmultbin[i]->GetStdDev(), rpres[1]);

      gv_mcos2->SetPoint(ip, hmultbin[i]->GetMean(), rpres[2]);
      gv_mcos2->SetPointError(ip, hmultbin[i]->GetStdDev(), rpres[3]);
      ip++;
    }
  }
  auto cc81 = new TCanvas("cc81","cc81");
  gv_mcos1->Draw("ALP");

  auto cc82 = new TCanvas("cc82","cc82");
  gv_mcos2->Draw("ALP");

  gv_mcos1->Write();
  gv_mcos2->Write();  
  GraphSave->Write();

}


UInt_t SetBranch()
{
  if(rChain == NULL) {
    std::cout << " no file is loaded " << std::endl;
    return 0;
  }

  std::cout << " Nentry ->  " << rChain->GetEntries() << std::endl;

  if(aArray != NULL)
    aArray->Clear();

  rChain->SetBranchAddress("STParticle",&aArray);
  rChain->SetBranchAddress("STFlow"    ,&aFlowArray);

  return rChain->GetEntries();
}


TString SetupOutputFile(TString fopt)
{
  //--------------------------------------------------
  //----- SETUP OUTPUT FILE --------------------------
  //--------------------------------------------------
  gSystem->cd("data");
  TString oVer = gSystem->Getenv("OUTVER");
  TString fName = fopt + oVer;

  if( oVer == "" ) {
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
      std::cout << checkName << " is existing. Do you recreate? (y/n)" << std::endl;
      TString sAns;
      std::cin >> sAns;
      if(sAns == "N" || sAns == "n") {
	std::cout << " Retry" << std::endl;
	exit(0);
      }
    }
    fName += ".root";
  }

  std::cout << "File " << fName << " is created. " << std::endl;  

  return fName;
  //--------------------------------------------------
}
