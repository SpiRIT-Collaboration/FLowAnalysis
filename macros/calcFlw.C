#include "FlowFunctions.h"
#include "openFlw.C"

Int_t  seltrackID = 4;
UInt_t selReactionPlanef = 10;
Int_t  seltrack;

TCanvas *cc[12];

//const Double_t ycm      = 0.388568; // 132Sn + 124Sn before Nov.15 201 //278MeV/u
//                         0: 132Sn + 124Sn,  1: 108Sn + 112Sn  //270MeV/u
const Double_t ycm[]      = {0.383198,  0.365587}; 
const Double_t ybeam_cm[] = {0.36599 ,  0.378516};
const Double_t ybeam_lb[] = {0.749188,  0.744103};
UInt_t sys[]              = {10, 10, 10, 10};


const UInt_t nspec = 5;
const UInt_t  nbin = 9;

//UInt_y   y_nbin  = 100;
Double_t y_min[] = {0., 0., 0., 0., 0.};
Double_t y_max[] = {0.8, 0.8, 0.8, 0.8, 0.8};
Double_t y_bin[] = {0.05, 0.05, 0.01, 0.01, 0.01};
Double_t y_binx  = 0.1;

Double_t pt_prmax  =  800.;
Double_t pt_dtmax  = 1000.;
Double_t pt_trmax  = 1400.;
Double_t pt_pimax  = 300.;

Double_t pt_prmin  =  0.;
Double_t pt_dtmin  =  0.;
Double_t pt_trmin  =  0.;
Double_t pt_pimin  =  0.;

UInt_t   pt_nbin   = 100;

Double_t pt_min   = 0;
Double_t pt_max   = 1100;
Double_t pt_dbin  = (pt_max - pt_min)/(Double_t)(pt_nbin-1); 

TString  partname[] = {"pi-","pi+","proton","deuteron","triton"};
UInt_t   partid[]   = {211, 211, 2212, 1000010020, 1000010030};
UInt_t   icol[]     = {7,5,2,4,3};
TString  iopt[]     = {"","same","same","same","same"};
UInt_t   imark[]    = {23, 33, 20, 21, 22};
  
UInt_t ic = -1;

const Int_t nbinx = 30;
TChain *rChain[4];
Int_t ntrack[7];
auto aArray = new TClonesArray("STParticle",100);


UInt_t m_bgn = 0;
UInt_t m_end = 1;

// histogram
TH1D* hptpm[4][nbin];
TH1D* hptpp[4][nbin];
TH1D* hptpr[4][nbin];
TH1D* hptdt[4][nbin];
TH1D* hpttr[4][nbin];


const UInt_t nphi = 30;
std::vector< std::vector < std::vector< Double_t > > > labphi;
std::vector< std::vector < std::vector< Double_t > > > subcos;

Int_t   iVer;
TString sVer;
Int_t   mtrack;

// pt dependence

Double_t rapid_max = 0.4;
Double_t rapid_min = 0.2;
UInt_t pxbooking();
UInt_t DrawCenterLine(TMultiGraph *mg);
// Double_t *vMean(vector<Double_t> &vec);
// Double_t *vn(UInt_t hm, vector<Double_t> &vphi);

void plotv1v2(UInt_t selid=2);
void PtDependece(UInt_t hrm=1);
void YDependece(UInt_t hrm=1);
void PhiYbin();
void PhiYbinf();
void PxDistribution(UInt_t nplot=0);
void hpt_plot();
void meanPx();
void dndy();
void CalculateResolution(UInt_t m, Double_t Phi, Double_t phi_sub);
void SaveRPResolution(Int_t sn, UInt_t size, Double_t *x, Double_t *xe, Double_t *y, Double_t *ye);
void ShiftingCorrection(STParticle *apar);


void calcFlw()
{
  gROOT->Reset();

  UInt_t ichain = 0;

  openFlw();

  for(UInt_t i = 0; i < 4; i++){
    rChain[ichain] = (TChain*)gROOT->FindObject(Form("rChain%d",i));
    if(rChain[ichain] != NULL) {    
      sys[ichain] = GetSystem(ichain);
      std::cout << " System " << ichain << " "  << sys[ichain] << "  -> " << sysName[sys[ichain]] << std::endl; 
      ichain++;
    }
  }

  if(rChain[0] == NULL)
    exit(0);

  m_end = ichain;

  std::cout << " ichain " << ichain << " m_end " << m_end << std::endl;
  
  gROOT->ProcessLine(".! grep -i void calcFlw.C | grep '//%%'");
}


void dndy()                       //%% Executable : Make plots of dNdy for p, d, t, pi+-
{

  //----- booking
  TH1D* hrap[4][5];
  TH1D* hnpart[4][5];
  TH1D* hmtrack[4];

  auto aLeg0 = new TLegend(0.7,0.7,0.9 ,0.9,"");
  auto aLeg1 = new TLegend(0.7,0.7,0.9 ,0.9,"");

  for(Int_t m = m_bgn; m < m_end; m++){

    for(UInt_t i = 0; i < 5; i++){

      UInt_t y_nbin = UInt_t((y_max[i] - y_min[i])/y_bin[i]);
      TString hname = Form("hrap%d_%d",m,i);
      TString htitle= partname[i] + "; Rapidity; dN/dy";
      hrap[m][i] = new TH1D(hname, htitle, y_nbin, y_min[i], y_max[i]);
      hrap[m][i] ->SetLineColor(icol[sys[m]]);

      hname  = Form("hnpart%d_%d",m,i);
      htitle = partname[i]+" ; Multiplicity";
      hnpart[m][i] = new TH1D(hname, htitle, 25,0,25);
      hnpart[m][i]->SetLineColor(icol[m]);
    }

    hmtrack[m] = new TH1D(Form("hmtrack%d",m),"Number of good tracks; Multiplicity",80, 0, 80);
    hmtrack[m] -> SetLineColor(icol[sys[m]]);
    
    aLeg0->AddEntry(hrap[m][4],sysName[sys[m]],"lp");
    aLeg1->AddEntry(hnpart[m][4],sysName[sys[m]],"lp");
  }

  //------------------------

  for(Int_t m = m_bgn; m < m_end; m++){

    aArray->Clear();
    mtrack = 0;

    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("ntrack",   ntrack);

    Int_t nEntry = rChain[m]->GetEntries();
    
    for(Int_t i = 0; i < nEntry; i++){
      rChain[m]->GetEntry(i);

      TIter next(aArray);
      STParticle *aPart = NULL;

      hmtrack[m]->Fill(mtrack);
      UInt_t npart[5] = {0,0,0,0,0};

      while( (aPart = (STParticle*)next()) ) {

	auto flag  = aPart->GetReactionPlaneFlag();
	auto pid   = aPart->GetPID();
	auto charg = aPart->GetCharge();
	auto p     = aPart->GetRotatedMomentum();
	auto rapid = aPart->GetRapidity();
	auto theta = aPart->GetRotatedMomentum().Theta();

	if(pid == partid[0] ){
	  if( charg < 0){
	    hrap[m][0]->Fill(rapid);
	    npart[0]++;
	  }
	  else {
	    hrap[m][1]->Fill(rapid);
	    npart[1]++;
	  }
	}
	
	if(flag > 1000){
	  for(UInt_t i = 2; i < 5; i++){
	    if(pid == partid[i]){
	      hrap[m][i]->Fill(rapid);
	      npart[i]++;
	    }
	  }
	}
      }
      for(UInt_t i = 0; i < 5; i++)                                                                                                       
        hnpart[m][i]->Fill(npart[i]);
    }
    

    std::cout << " -------------------- " << std::endl;
    for(UInt_t i = 0; i < 5; i++) {

      //----- Normalizing
      auto scl = 1./nEntry * (Double_t)hrap[m][i]->GetNbinsX() / (y_max[i] - y_min[i]);
      hrap[m][i]->Scale(scl);

      scl = 1./nEntry;
      hnpart[m][i]->Scale(scl);
      auto mean   = hnpart[m][i]->GetMean();
      auto meaner = hnpart[m][i]->GetMeanError();
      std::cout << setw(12) << partname[i] << " : " << mean << " +- " << meaner << std::endl;

    }
    hmtrack[m]->Scale(1./nEntry);
    std::cout << " -------------------- " << std::endl;
  }


  //----- Drawing
  //----- canvas
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1400,500);
  cc[ic]->Divide(3,1);

  const UInt_t nsys = m_end - 1;


  Double_t hmax[nsys];


  // plot p, d, and t
  for(UInt_t ip = 2; ip < 5; ip++){
    cc[ic]->cd(ip-1); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hrap[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    cout << "ip : " << ip  << " gmax " << gmax << endl;
    
    hrap[0][ip]->SetMaximum(gmax*1.1);

    for(UInt_t m = m_bgn; m < m_end; m++) 
      hrap[m][ip]->Draw(iopt[m]); 
     
    if(ip == 4)
      aLeg0->Draw();
  }


  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
  cc[ic]->Divide(2,1);

 
  for(UInt_t ip = 0; ip < 2; ip++){
    cc[ic]->cd(ip+1); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hrap[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    hrap[0][ip]->SetMaximum(gmax*1.1);
    for(UInt_t m = m_bgn; m < m_end; m++) 
      hrap[m][ip]->Draw(iopt[m]); 
  

    if(ip == 1)
      aLeg1->Draw();
  }


  // //----cc1
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,800);
  cc[ic]->Divide(3,2);
  
  for(UInt_t ip = 2; ip < 5; ip++){
    cc[ic]->cd(ip-1); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hnpart[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    hnpart[0][ip]->SetMaximum(gmax*1.1);
    for(UInt_t m = m_bgn; m < m_end; m++) 
      hnpart[m][ip]->Draw(iopt[m]);

    if(ip == 4)
      aLeg1->Draw();
  }

  io = 0;
  for(UInt_t ip = 0; ip < 2; ip++){
    cc[ic]->cd(ip+4); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hnpart[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    hnpart[0][ip]->SetMaximum(gmax*1.1);
    for(UInt_t m = m_bgn; m < m_end; m++)
      hnpart[m][ip]->Draw(iopt[m]);
    

    if(ip == 1)
      aLeg1->Draw();
  }
  
}


//________________________________//%%
void meanPx()                     //%% Executable : Make  plots of <px> vs rapidity 
{
  PxDistribution(1);

  TGraphErrors  *gpr[4];
  TGraphErrors  *gdt[4];
  TGraphErrors  *gtr[4];
  TGraphErrors  *gpm[4];
  TGraphErrors  *gpp[4];

  auto mgpdt = new TMultiGraph();
  mgpdt->SetTitle("p,d and t; Rapidity_lab; <px>/A [MeV/c]");
  auto aLeg0 = new TLegend(0.1,0.7,0.35,0.9,"");

  auto mgpi  = new TMultiGraph();
  mgpi->SetTitle("#pi^{+-}; Rapidity_lab; <px> [MeV/c]");
  auto aLeg1 = new TLegend(0.1,0.7,0.35,0.9,"");

  for(Int_t m = m_bgn; m < m_end; m++){
   
    Double_t rap[nbin];
    Double_t rape[nbin];
 
    Double_t mptpm[nbin];
    Double_t mptpp[nbin];
    Double_t mptpr[nbin];
    Double_t mptdt[nbin];
    Double_t mpttr[nbin];

    Double_t mptpme[nbin];
    Double_t mptppe[nbin];
    Double_t mptpre[nbin];
    Double_t mptdte[nbin];
    Double_t mpttre[nbin];

    UInt_t   npr = 0;
    UInt_t   ndt = 0;
    UInt_t   ntr = 0;
    UInt_t   npm = 0;
    UInt_t   npp = 0;


    // get mean 
    for(UInt_t k = 0; k < nbin; k++){
      rap[k]   = y_min[0] + y_binx*(k+0.5);
      rape[k]  = y_binx*0.5;

      if(hptpm[m][k]->GetEntries() > 0){
	mptpm[npm]  = hptpm[m][k]->GetMean(); 
	mptpme[npm] = hptpm[m][k]->GetMeanError(); 
	npm++;
      }

      if(hptpp[m][k]->GetEntries() > 0){
	mptpp[npp]  = hptpp[m][k]->GetMean();
	mptppe[npp] = hptpp[m][k]->GetMeanError();
	npp++;
      }

      if(hptpr[m][k]->GetEntries() > 0){
	mptpr[npr]  = hptpr[m][k]->GetMean();
	mptpre[npr] = hptpr[m][k]->GetMeanError();
	npr++;
      }

      if(hptdt[m][k]->GetEntries() > 0){
	mptdt[ndt]  = hptdt[m][k]->GetMean()/2.;
	mptdte[ndt] = hptdt[m][k]->GetMeanError();
	ndt++;
      }

      if(hpttr[m][k]->GetEntries() > 0){
	mpttr[ntr] = hpttr[m][k]->GetMean()/3.;
	mpttre[ntr] = hpttr[m][k]->GetMeanError();
	ntr++;
      }
    }


    gpr[m] = new TGraphErrors(npr, rap, mptpr, rape, mptpre);
    gdt[m] = new TGraphErrors(ndt, rap, mptdt, rape, mptdte);
    gtr[m] = new TGraphErrors(ntr, rap, mpttr, rape, mpttre);
    gpm[m] = new TGraphErrors(npm, rap, mptpm, rape, mptpme);
    gpp[m] = new TGraphErrors(npp, rap, mptpp, rape, mptppe);
   
  

    TString atitle = Form("gpr%d",m);
    gpr[m]->SetName(atitle);
    gpr[m]->SetLineColor(icol[sys[m]]);
    gpr[m]->SetMarkerStyle(imark[sys[m]]);
    gpr[m]->SetMarkerColor(icol[sys[m]]);
    gpr[m]->SetTitle("Proton; y_lab; <Px> (MeV/c)");
    mgpdt->Add(gpr[m],"lp");
    aLeg1->AddEntry(gpr[m],"proton   "+sysName[sys[m]],"lp");

    atitle = Form("gdt%d",m);
    gdt[m]->SetName(atitle);
    gdt[m]->SetLineColor(icol[sys[m]]);
    gdt[m]->SetMarkerStyle(21);
    gdt[m]->SetMarkerColor(icol[sys[m]]);
    gdt[m]->SetTitle("Deuteron; y_lab; <Px> (MeV/c)");
    mgpdt->Add(gdt[m],"lp");
    aLeg1->AddEntry(gdt[m],"deuteron "+sysName[sys[m]],"lp");

    atitle = Form("gtr%d",m);
    gtr[m]->SetName(atitle);
    gtr[m]->SetLineColor(icol[sys[m]]);
    gtr[m]->SetMarkerStyle(22);
    gtr[m]->SetMarkerColor(icol[sys[m]]);
    gtr[m]->SetTitle("Triton; y_lab; <Px> (MeV/c)");
    mgpdt->Add(gtr[m],"lp");
    aLeg1->AddEntry(gtr[m],"trition  "+sysName[sys[m]],"lp");

    atitle = Form("gpm%d",m);
    gpm[m]->SetName(atitle);
    gpm[m]->SetLineColor(icol[sys[m]]);
    gpm[m]->SetMarkerStyle(23);
    gpm[m]->SetMarkerColor(icol[sys[m]]);
    gpm[m]->SetTitle("pi-; y_lab; <Px> (MeV/c)");
    mgpi->Add(gpm[m],"ip");
    aLeg0->AddEntry(gpm[m],"#pi^{-}   "+sysName[sys[m]],"lp");

    atitle = Form("gpp%d",m);
    gpp[m]->SetName(atitle);
    gpp[m]->SetLineColor(icol[sys[m]]);
    gpp[m]->SetMarkerStyle(33);
    gpp[m]->SetMarkerColor(icol[sys[m]]);
    gpp[m]->SetTitle("pi+; y_lab; <Px> (MeV/c)");
    mgpi->Add(gpp[m],"ip");
    mgpi->Add(gpp[m],"ip");
    aLeg0->AddEntry(gpp[m],"#pi^{+}   "+sysName[sys[m]],"lp");
  }

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mgpdt -> Draw("a");
  aLeg1->Draw();


  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mgpi  ->Draw("a");
  aLeg0->Draw();
}

UInt_t pxbooking()  // used by meanPx()
{
  //----- booking

  for(Int_t m = m_bgn; m < m_end; m++){

    for(UInt_t i = 0; i < nbin; i++){

      Double_t yL = y_min[0] + y_binx*i;
      Double_t yU = yL + y_binx;

      TString hname = Form("hptpr%d%d",m,i);
      TString htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpr[m][i] = new TH1D(hname,"Proton :  "+htitle,pt_nbin, pt_prmin, pt_prmax);
      hptpr[m][i] ->SetLineColor(icol[sys[m]]);

      hname = Form("hptdt%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptdt[m][i] = new TH1D(hname,"Deuteron :"+htitle,pt_nbin, pt_dtmin, pt_dtmax);
      hptdt[m][i] ->SetLineColor(icol[sys[m]]);

      hname = Form("hpttr%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hpttr[m][i] = new TH1D(hname,"Triton :  "+htitle,pt_nbin, pt_trmin, pt_trmax);
      hpttr[m][i] ->SetLineColor(icol[sys[m]]);

      hname = Form("hptpm%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpm[m][i] = new TH1D(hname,"Pi- :     "+htitle,pt_nbin, pt_pimin, pt_pimax);
      hptpm[m][i] ->SetLineColor(icol[sys[m]]);

      hname = Form("hptpp%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpp[m][i] = new TH1D(hname,"Pi+ :     "+htitle,pt_nbin, pt_pimin, pt_pimax);
      hptpp[m][i] ->SetLineColor(icol[sys[m]]);

    }
  }
  return 1;
}



void PxDistribution(UInt_t nplot)  // used by meanPx()
{
  pt_prmin  =  -pt_prmax;
  pt_dtmin  =  -pt_dtmax;
  pt_trmin  =  -pt_trmax;
  pt_pimin  =  -pt_pimax;

  pxbooking();
  
  for(Int_t m = m_bgn; m < m_end; m++){

    
    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("ntrack" ,ntrack);


    Int_t nEntry = rChain[m]->GetEntries();
    cout << " Number of events " << nEntry << endl;


    for(Int_t i = 0; i < nEntry; i++){
      aArray->Clear();

      rChain[m]->GetEntry(i);

      if(ntrack[2] == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	if( aPart->GetReactionPlaneFlag() == 0)
	  continue;

	auto pid   = aPart->GetPID();
	auto charg = aPart->GetCharge();
	auto rapid = aPart->GetRapidity();
	//	auto vp    = aPart->GetFlattenMomentum();
	auto vp    = aPart->GetRotatedMomentum();
	auto dltphi= aPart->GetAzmAngle_wrt_RP();;
	vp.SetPhi(dltphi);

	auto px    = vp.X();
	//	auto px    = vp.Pt();

	//rapid = (rapid - ycm[sys[m]]) / ybeam_cm[sys[m]];

	for(UInt_t k = 0; k < nbin; k++){
	  
	  if(rapid < y_min[0] + y_binx * (Double_t)(k+1)) {

	    if(pid == partid[0] && charg < 0)
	      hptpm[m][k]->Fill(px);
	    else if(pid == partid[1] && charg > 0)
	      hptpp[m][k] ->Fill(px);
	    else if(pid == partid[2])
	      hptpr[m][k] ->Fill(px);
	    else if(pid == partid[3])
	      hptdt[m][k] ->Fill(px);
	    else if(pid == partid[4])
	      hpttr[m][k] ->Fill(px);

	    break;
	    
	  }
	}

      }
    }

  }
  
  if(nplot == 0) hpt_plot(); 
  
}
 

void hpt_plot()
{

  for(Int_t m = m_bgn; m < m_end; m++){

    Int_t nEntry = rChain[m]->GetEntries();

    auto sclp  = 1./nEntry * (Double_t)pt_nbin / (pt_prmax - pt_prmin) * y_binx;
    auto scld  = 1./nEntry * (Double_t)pt_nbin / (pt_dtmax - pt_dtmin) * y_binx;
    auto sclt  = 1./nEntry * (Double_t)pt_nbin / (pt_trmax - pt_trmin) * y_binx;
    auto sclpi = 1./nEntry * (Double_t)pt_nbin / (pt_pimax - pt_pimin) * y_binx;

    for(UInt_t j = 0; j < nbin; j++){
      hptpr[m][j]->Scale(sclp);
    }

    ic = 0;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else {
      cc[ic]->cd();

      for(UInt_t j = 0; j < nbin; j++){
	cc[ic]->cd(j+1);
	hptpr[1][j]->Draw(iopt[0]);
	hptpr[0][j]->Draw(iopt[1]);

      }
    }

    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else
      cc[ic]->cd();

    for(UInt_t j = 0; j < nbin; j++){
      cc[ic]->cd(j+1);
      hptdt[m][j]->Scale(scld);
      hptdt[m][j]->Draw(iopt[m]);
    }

    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else
      cc[ic]->cd();

    for(UInt_t j = 0; j < nbin; j++){
      cc[ic]->cd(j+1);
      hpttr[m][j]->Scale(sclt);
      hpttr[m][j]->Draw(iopt[m]);
    }

    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else
      cc[ic]->cd();

    for(UInt_t j = 0; j < nbin; j++){
      cc[ic]->cd(j+1);
      hptpm[m][j]->Scale(sclpi);
      hptpm[m][j]->Draw(iopt[m]);
    }

    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else
      cc[ic]->cd();

    for(UInt_t j = 0; j < nbin; j++){
      cc[ic]->cd(j+1);
      hptpp[m][j]->Scale(sclpi);
      hptpp[m][j]->Draw(iopt[m]);
    }
  }
}


void comparev1v2(UInt_t m = 0)                 //%% Executable : compare v1 and v2 among p/d/t
{
  TH1D *hphi[5][nbin][4];
  TH1D *h2phi[5][nbin][4];
  TH1D *hcos1phi[5][nbin][4];
  Int_t nphi    = 30;
  Double_t dphi = 2./(Double_t)nphi;

  TString hname;
  
  for(Int_t m = m_bgn; m < m_end; m++){

    for(UInt_t pn = 0; pn < 5; pn++){
      for(UInt_t xbin = 0; xbin < nbin; xbin++){
	hname = Form("hphi%d%d_%d",pn,xbin,m);
	hphi[pn][xbin][m]     = new TH1D(hname,partname[pn],60,-3.15,3.15);

	hname = Form("h2phi%d%d_%d",pn,xbin,m);
	h2phi[pn][xbin][m] = new TH1D(hname,partname[pn]   ,nphi,-3.15,3.15);
	
	hname = Form("hcos1phi%d%d_%d",pn,xbin,m);
	hcos1phi[pn][xbin][m] = new TH1D(hname,partname[pn],nphi,-1.,1.);

      }
    }

    vector< vector< vector<Double_t> > >  bphi(5);
    vector< vector< vector<Double_t> > > rpdbin(5);
    for(UInt_t j = 0; j < 5; j++){
      bphi[j].resize(nbin);
      rpdbin[j].resize(nbin);
    }
     
    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("ntrack" ,ntrack);

    Int_t nEntry = rChain[m]->GetEntries();
    cout << " Number of events " << nEntry << endl;

    for(Int_t i = 0; i < nEntry; i++){
      aArray->Clear();
      rChain[m]->GetEntry(i);

      if(ntrack[2] == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	if( aPart->GetReactionPlaneFlag() == 0)
          continue;

        auto pid   = aPart->GetPID();
        auto charg = aPart->GetCharge();
        auto rapid = aPart->GetRapidity();
        auto vp    = aPart->GetRotatedMomentum();
	auto dltphi= aPart->GetAzmAngle_wrt_RP();;
        vp.SetPhi(dltphi);

        auto px    = vp.X();

	UInt_t xbin = nbin + 1;
	for(UInt_t k = 0; k < nbin; k++){
          if(rapid < y_min[0] + y_binx * (Double_t)(k+1) && rapid < y_max[0]){ 
	    xbin = k;
	    break;
	  }
	}
	
	UInt_t selid = 5;
	if(pid == partid[0] && charg < 0)
	  selid = 0;
	else if(pid == partid[1] && charg > 0)
	  selid = 1;
	else if(pid == partid[2])
	  selid = 2;
	else if(pid == partid[3])
	  selid = 3;
	else if(pid == partid[4])
	  selid = 4;

	
	if(selid < 5 && xbin < nbin){
	  hphi[selid][xbin][m]->Fill(dltphi);
	  hcos1phi[selid][xbin][m]->Fill(cos(dltphi));
	  h2phi[selid][xbin][m]->Fill(2.*dltphi);

	  rpdbin[selid][xbin].push_back(rapid);
	  bphi[selid][xbin].push_back(dltphi);
	}
      }
    }

    Double_t rap[5][nbin];
    Double_t rape[5][nbin];

    Double_t v1[5][nbin];
    Double_t v1e[5][nbin];
    Double_t v2[5][nbin];
    Double_t v2e[5][nbin];
    Int_t   npar[5];

    TGraphErrors *gv_v1[5];
    TGraphErrors *gv_v2[5];


    auto aLeg_v1 = new TLegend(0.1,0.7,0.35,0.9,"v1");
    auto mg_v1 = new TMultiGraph();
    mg_v1->SetTitle("; Rapidity; v1");
    
    auto aLeg_v2 = new TLegend(0.1,0.7,0.35,0.9,"v2");
    auto mg_v2 = new TMultiGraph();
    mg_v2->SetTitle("; Rapidity; v2");

    for(UInt_t k = 0; k < 5; k++) npar[k] = 0;

    for(UInt_t pn = 2; pn < 5; pn++){

      for(UInt_t xbin = 0; xbin < nbin; xbin++){
	auto *getx = new Double_t[2];
	auto *gety = new Double_t[3];

	if( rpdbin[pn][xbin].size() > 1 ){
	  getx = vMean( rpdbin[pn][xbin] );

	  if(getx[0] > -900 ) {
	    gety = vn( 1, bphi[pn][xbin] );
	  
	    rap[pn][xbin]  = getx[0];
	    rape[pn][xbin] = getx[1];
	    
	    v1[pn][xbin]   = gety[0];
	    v1e[pn][xbin]  = gety[1];
	  
	    gety = vn( 2, bphi[pn][xbin] );

	    v2[pn][xbin]   = gety[0];
	    v2e[pn][xbin]   = gety[1];
	  }
	  npar[pn]++;

	  std::cout << setw(5) << pn << " : " << xbin 
		    << setw(12) << "  " << rap[pn][xbin]
		    <<  " v1 " << setw(12) << v1[pn][xbin] << " +- " << v1e[pn][xbin]  
		    <<  " v2 " << setw(12) << v2[pn][xbin] << " +- " << v2e[pn][xbin]  << std::endl;
	}
      }


      //%%%%%%%%%%%%%%%%%%%%%%%%%%%
      gv_v1[pn] = new TGraphErrors(npar[pn], rap[pn], v1[pn], rape[pn], v1e[pn]);
      gv_v2[pn] = new TGraphErrors(npar[pn], rap[pn], v2[pn], rape[pn], v1e[pn]);
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%

      gv_v1[pn]->SetTitle(partname[pn]);
      gv_v1[pn]->SetLineColor(icol[pn]);
      gv_v1[pn]->SetMarkerStyle(imark[pn]);
      gv_v1[pn]->SetMarkerColor(icol[pn]);
      mg_v1->Add(gv_v1[pn]);
      aLeg_v1->AddEntry(gv_v1[pn], partname[pn],"lp");

      gv_v2[pn]->SetTitle(partname[pn]);
      gv_v2[pn]->SetLineColor(icol[pn]);
      gv_v2[pn]->SetMarkerStyle(imark[pn]);
      gv_v2[pn]->SetMarkerColor(icol[pn]);
      mg_v2->Add(gv_v2[pn]);
      aLeg_v2->AddEntry(gv_v2[pn],partname[pn],"lp");
    }

    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    mg_v1->Draw("ALP");
    aLeg_v1->Draw();

    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    mg_v2->Draw("ALP");
    aLeg_v2->Draw();


    
    ic++; 
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,1000);
    cc[ic]->Divide(3,7); id = 1;
    
    for(UInt_t xbin = 0; xbin < 8; xbin++){
      for(UInt_t pn = 2; pn < 5; pn++){
	cc[ic]->cd(id); id++;
	hname = Form(partname[pn]+" Rapidity %f; #delta#phi; v1",rap[pn][xbin]); 
	hphi[pn][xbin][m]->SetTitle(hname);
	hphi[pn][xbin][m]->SetTitleSize(0.1);
	hphi[pn][xbin][m]->Draw();
      }
    }

    ic++; 
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,1000);
    cc[ic]->Divide(3,7); id = 1;
    for(UInt_t xbin = 0; xbin < 8; xbin++){
      for(UInt_t pn = 2; pn < 5; pn++){
	cc[ic]->cd(id); id++;
	hname = Form(partname[pn]+" Rapidity %f; #delta#phi; v2",rap[pn][xbin]); 
	hcos1phi[pn][xbin][m]->SetTitle(hname);
	hcos1phi[pn][xbin][m]->Draw();
      }
    }
    ic++; 
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,1000);
    cc[ic]->Divide(3,7); id = 1;
    for(UInt_t xbin = 0; xbin < 8; xbin++){
      for(UInt_t pn = 2; pn < 5; pn++){
	cc[ic]->cd(id); id++;
	hname = Form(partname[pn]+"Rapidity %f",rap[pn][xbin]); 
	h2phi[pn][xbin][m]->SetTitle(hname);
	h2phi[pn][xbin][m]->SetTitleSize(0.1);
	h2phi[pn][xbin][m]->Draw();
      }
    }
  }
}


void plotv1v2(UInt_t selid=2)     //%% Executable : v1 and v2 as a function of rapidity
{
  if(selid > 4 ) return;

  gStyle->SetGridColor(7);
  gStyle->SetGridStyle(1);

  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  Int_t RPflag = 0;
  if(selid >= 2) RPflag = 10;

  cout << " Particle " << partname[selid] << endl;

  auto haccp = new TH2D("haccp",partname[selid]+"; Rapidity ; Pt [MeV/c]",100,0.,1.5,100,  0.,800.);
  auto haccx = new TH2D("haccx",partname[selid]+"; Rapidity ; Px [MeV/c]",100,0.,1.5,100,-500,500.);
  

  for(Int_t m = m_bgn; m < m_end; m++){

    auto ymin = y_min[0];

    aArray->Clear();
    
    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("ntrack" ,ntrack);

    vector< vector<Double_t> > bphi(nbin);
    vector< vector<Double_t> > rpdbin(nbin);


    for(UInt_t k = 0; k < nbin; k++){
      bphi[k].clear();
      rpdbin[k].clear();
    }


    Int_t nevt = rChain[m]->GetEntries();
    cout << " Number of events " << nevt << endl;


    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      if(ntrack[2] == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto pid = aPart->GetPID();

	if(pid == partid[selid] && aPart->GetCharge() == pcharge ){

	  auto rapid = aPart-> GetRapidity();

	  haccp -> Fill(rapid, aPart->GetRotatedMomentum().Pt());
	  haccx -> Fill(rapid, aPart->GetRotatedMomentum().Pt()*sin(aPart->GetAzmAngle_wrt_RP()));

	  
	  if(aPart->GetIndividualRPAngle() > -9 ) {
	  
	    for(UInt_t k = 0; k < nbin; k++){

	      if(rapid < ymin + y_binx*(Double_t)(k+1)) {
		bphi[k].push_back(aPart->GetAzmAngle_wrt_RP());
		rpdbin[k].push_back(rapid);
	
		//		cout << "RP angle " << aPart->GetIndividualRPAngle() << " rapidity " << rapid  << " " << k << endl;
		break;
	      }
	    }
	  }
	}
      }
    }

    Double_t xval[nbin];
    Double_t xvale[nbin];
    Double_t yval1[nbin];
    Double_t yval1e[nbin];
    Double_t yval2[nbin];
    Double_t yval2e[nbin];

    UInt_t jncount = 0;
    Int_t  jnfirst = -1;
    for(UInt_t jn = 0; jn < nbin; jn++){

      Double_t *getx = new Double_t[2];
      Double_t *gety = new Double_t[2];
      getx   = vMean(rpdbin[jn]);

      if(getx[0] > -900 && bphi[jn].size() > 0){

	gety   = vn(1, bphi[jn]);

	xval[jn]  = getx[0];
	yval1[jn] = gety[0]; ///rp_res[m];

	xvale[jn]  = getx[1];
	yval1e[jn] = gety[1];


	gety    = vn(2, bphi[jn]);
	yval2[jn] = gety[0]; ///rp_res[m];
	yval2e[jn]= gety[1];
	
	jncount++;
	if(jnfirst == -1)
	  jnfirst = jn;
      }
      else {
	xval[jn] = -999.;

      }
    }

    //print
    std::cout << " ---- Resutls ---------------------" << std::endl;
    for(UInt_t i = 0; i < nbin; i++)
      std::cout << setw(5) << i << " : " << setw(12)
		<< xval[i] << " +-" << xvale[i] 
		<<  " v1 " << setw(12) << yval1[i] << " +- " << yval1e[i] << "  w  " << bphi[i].size() << std::endl;  
    std::cout << " ----------------------------------" << std::endl;
    for(UInt_t i = 0; i < nbin; i++)
      std::cout << setw(5) << i << " : " << setw(12)
		<< xval[i] << " +-"  << xvale[i]
		<<  " v2 " << setw(12) << yval2[i] << " +- " << yval2e[i] << "  w  " << bphi[i].size() << std::endl;  

    cout << " jncount " << jncount << " jnfirst " << jnfirst << endl;


    ///    const UInt_t nsl = 9;
    const UInt_t nsl = jncount;
    Double_t rap[nsl];
    Double_t rape[nsl];
    Double_t v1[nsl];
    Double_t v1e[nsl];
    Double_t v2[nsl];
    Double_t v2e[nsl];

    for(UInt_t isl = 0; isl < nsl; isl++){
      UInt_t iisl = jnfirst + isl;
      rap[isl]  = xval[iisl];
      rape[isl] = xvale[iisl];
      
      v1[isl]   = yval1[iisl];
      v1e[isl]  = yval1e[iisl];
      v2[isl]   = yval2[iisl];
      v2e[isl]  = yval2e[iisl];

    }

    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

    auto gv_v1 = new TGraphErrors(nsl, rap, v1, rape, v1e);
    gv_v1->SetName("gv_v1");

    // gv_v1->SetMaximum(0.1);
    // gv_v1->SetMinimum(-0.1);
    gv_v1->SetLineColor(4);
    gv_v1->SetMarkerStyle(20);
    gv_v1->SetMarkerColor(4);
    gv_v1->SetTitle(partname[selid]+"; Ycm/Ycm_beam ; v1(a.u.)");
   
    gv_v1->Draw("ALP");

    auto xmin = gv_v1->GetYaxis()->GetXmin();
    auto xmax = gv_v1->GetYaxis()->GetXmax();
    auto aLineX = new TLine(0, xmin, 0,  xmax);
    aLineX->SetLineColor(1);
    aLineX->SetLineStyle(3);
    aLineX->Draw();

    xmin = gv_v1->GetXaxis()->GetXmin();
    xmax = gv_v1->GetXaxis()->GetXmax();
    auto aLineY = new TLine(xmin, 0, xmax, 0);
    aLineY->SetLineColor(1);
    aLineY->SetLineStyle(3);
    aLineY->Draw();


    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

    
    auto gv_v2 = new TGraphErrors(nsl, rap, v2, rape, v2e);
    gv_v2->SetName("gv_v2");

    // gv_v2->SetMaximum(0.015);
    // gv_v2->SetMinimum(-0.015);
    gv_v2->SetLineColor(2);
    gv_v2->SetMarkerStyle(20);
    gv_v2->SetMarkerColor(2);
    gv_v2->SetTitle(partname[selid]+"; Ycm/Ycm_beam ; v2(a.u.)");

    gv_v2->Draw("ALP");

  
    UInt_t id = 1;
    ic++;    
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);


    //    cc[ic]->Divide(1,2);
    
    //    cc[ic]->cd(id); id++;
    haccx -> Draw("colz");
    
    //    cc[ic]->cd(id); 
    //    haccx -> Draw("colz");

  }
}

UInt_t DrawCenterLine(TMultiGraph *mg)
{
  // Center Line
  auto xmin = mg->GetYaxis()->GetXmin();
  auto xmax = mg->GetYaxis()->GetXmax();
  auto aLineX = new TLine(0, xmin, 0,  xmax);
  aLineX->SetLineColor(1);
  aLineX->SetLineStyle(3);
  aLineX->Draw();

  xmin = mg->GetXaxis()->GetXmin();
  xmax = mg->GetXaxis()->GetXmax();
  auto aLineY = new TLine(xmin, 0, xmax, 0);
  aLineY->SetLineColor(1);
  aLineY->SetLineStyle(3);
  aLineY->Draw();

  return 1;
}




void check(UInt_t selid=2)
{
  if(selid > 4 ) return;


  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  Int_t RPflag = 0;
  if(selid >= 2) RPflag = 10;


  cout << " Particle " << partname[selid] << endl;

  auto h1 = new TH2D("h1"," ycm/ycm_beam ; #Delta#phi ",100,-1.1, 1.5,100,-3.2, 3.2);
  auto h2 = new TH2D("h2"," ycm/ycm_beam ; #Delta#phi ",100,-1.1, 1.5,100,-3.2, 3.2);

  auto h3 = new TH2D("h3"," #theta ; #Delta#phi ",200,0.,1.6,200,-3.2, 3.2);
  auto h4 = new TH2D("h4"," #theta ; #Delta#phi ",200,0.,1.6,200,-3.2, 3.2);


  Double_t rp_res[2] = {0.107183, 0.0003226};
  Double_t rp_rese[2]= {0.69 ,0.71};

  Double_t aoq;
  Double_t z;

  for(Int_t m = m_bgn; m < m_end; m++){

    aArray->Clear();
    
    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("ntrack",ntrack);
    rChain[m]->SetBranchAddress("aoq",&aoq);
    rChain[m]->SetBranchAddress("z",&z);
    rChain[m]->SetBranchAddress("ntrack" ,&ntrack);

    vector< vector<Double_t> > bphi  (nbin);
    vector< vector<Double_t> > rpdbin(nbin);



    Int_t nevt = rChain[m]->GetEntries();
    cout << " Number of events " << nevt << endl;


    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      if(ntrack[2] == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto pid = aPart->GetPID();
	auto rapid = aPart-> GetRapidity();
	rapid = (rapid - ycm[sys[m]]) / ybeam_cm[sys[m]];

	if(aPart->GetReactionPlaneFlag() >= RPflag){

	  h1 -> Fill(rapid, aPart->GetFlattenMomentum().Phi());
	  h3 -> Fill(aPart->GetFlattenMomentum().Theta(), aPart->GetFlattenMomentum().Phi());

	  if(pid == partid[selid] && aPart->GetCharge() == pcharge){

	    h2 -> Fill(rapid, aPart->GetFlattenMomentum().Phi());
	    h4 -> Fill(aPart->GetFlattenMomentum().Theta(), aPart->GetFlattenMomentum().Phi());
	  }
	}
      }
    }

    UInt_t id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

    cc[ic]->Divide(1,2);
    id = 1;
    cc[ic]->cd(id); id++;
    h1 -> Draw("colz");
    
    cc[ic]->cd(id); 
    h2 -> Draw("colz");

    ic++;

    id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->Divide(1,2);
    
    cc[ic]->cd(id); id++;
    h3 -> Draw("colz");
    
    cc[ic]->cd(id); 
    h4 -> Draw("colz");

    ic++;
  }
}


void PtDependece(UInt_t hrm = 1)
{

  cout << " PtDependence v"<<hrm << endl;

  UInt_t   selid[3] = {2,3,4};

  for(Int_t m = m_bgn; m < m_end; m++){

    aArray->Clear();
    rChain[m]->SetBranchAddress("STParticle",&aArray);


    vector< vector< vector<Double_t> > > bphi;
    vector< vector< vector<Double_t> > > ptbin;
    bphi.resize(3);
    ptbin.resize(3);

    for(UInt_t im = 0; im < 3; im++){
      bphi[im] .resize(pt_nbin);
      ptbin[im].resize(pt_nbin);
    }


    Int_t nevt = rChain[m]->GetEntries();
    cout << " Number of events " << nevt << endl;

    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto pid   = aPart->GetPID();
	auto rapid = aPart-> GetRapidity();
      	rapid = (rapid - ycm[sys[m]])/ybeam_cm[sys[m]];

	if(  aPart->GetBestTrackFlag() > 0 && (rapid < rapid_max && rapid > rapid_min) ){

	  auto pt = aPart->GetRotatedMomentum().Pt();

	  UInt_t k = 0; 
	  while( k < pt_nbin ){
	    if( pt < pt_min + pt_dbin*(Double_t)(k+1) )
	      break;
	    
	    k++;
	  }

	  //	  cout << " K " << k << " " << pt << endl;

	  if(k < pt_nbin){
	    for(UInt_t j = 0; j < 3; j++){
	      if( pid == partid[selid[j]] ){
		bphi[j][k].push_back(aPart->GetAzmAngle_wrt_RP());
		ptbin[j][k].push_back(pt);
		
	      }
	    }
	  }
	}
      }
    }

    Double_t xval [3][pt_nbin];
    Double_t xvale[3][pt_nbin];
    Double_t yval [3][pt_nbin];
    Double_t yvale[3][pt_nbin];

    for(UInt_t j = 0; j < 3; j++) {
      for(UInt_t jn = 0; jn < pt_nbin; jn++) {

	auto getx = vMean(ptbin[j][jn]);
	auto gety = vn(hrm, bphi[j][jn]);

	xval [j][jn] = getx[0];
	xvale[j][jn] = getx[1];

	yval [j][jn] = gety[0];
	yvale[j][jn] = gety[1];

      }
    }
    

    const UInt_t nsl = 11;
    UInt_t iisl = 0;
    auto gr0 = new TGraphErrors(nsl, &xval[0][iisl], &yval[0][iisl], &xvale[0][iisl], &yvale[0][iisl]);
    gr0->SetName("gr0");
    gr0->SetLineColor(2);
    gr0->SetMarkerStyle(20);
    gr0->SetMarkerColor(2);

    auto gr1 = new TGraphErrors(nsl, &xval[1][iisl], &yval[1][iisl], &xvale[1][iisl], &yvale[1][iisl]);
    gr1->SetName("gr1");
    gr1->SetLineColor(4);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerColor(4);

    auto gr2 = new TGraphErrors(nsl, &xval[2][iisl], &yval[2][iisl], &xvale[2][iisl], &yvale[2][iisl]);
    gr2->SetName("gr2");
    gr2->SetLineColor(8);
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerColor(8);


    auto mg = new TMultiGraph();
    TString aTitle = Form("; Pt[MeV/c]; v%d(a.u.)",hrm); 
    mg->SetTitle(aTitle);
    mg->SetName("mg");
    mg->Add(gr0,"lp");
    mg->Add(gr1,"lp");
    mg->Add(gr2,"lp");

    // Legend                                                                                                                               

    TString sconf = Form("%4.2f < y < %4.2f",rapid_min,rapid_max);
    auto aLeg = new TLegend(0.1,0.7,0.35,0.9,sconf);
    aLeg->AddEntry(gr0,partname[selid[0]],"lp");
    aLeg->AddEntry(gr1,partname[selid[1]],"lp");
    aLeg->AddEntry(gr2,partname[selid[2]],"lp");


    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    //    gr0->Draw();

    mg  ->Draw("a");
    aLeg->Draw();

    ic++;
  }
}


void YDependece(UInt_t hrm=1)
{
  cout << " YDependece v"<< hrm  << endl;

  UInt_t   selid[3] = {2,3,4};

  for(Int_t m = m_bgn; m < m_end; m++){

    aArray->Clear();
    rChain[m]->SetBranchAddress("STParticle",&aArray);

    vector< vector< vector<Double_t> > > bphi;
    vector< vector< vector<Double_t> > > xbin;
    bphi.resize(3);
    xbin.resize(3);

    for(UInt_t im = 0; im < 3; im++){
      bphi[im].resize(nbin);
      xbin[im].resize(nbin);
    }


    Int_t nevt = rChain[m]->GetEntries();
    cout << " Number of events " << nevt << endl;

    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto pid   = aPart->GetPID();
	auto rapid = aPart-> GetRapidity();
	rapid = (rapid - ycm[sys[m]])/ybeam_cm[sys[m]];

	if(  aPart->GetBestTrackFlag() > 0  ){

	  UInt_t k = 0; 
	  while( k < nbin ){
	    if( rapid < y_min[0] + y_binx*(Double_t)(k+1) )
	      break;
	    
	    k++;
	  }

	  if(k < nbin){
	    for(UInt_t j = 0; j < 3; j++){
	      if( pid == partid[selid[j]] ){
		bphi[j][k].push_back(aPart->GetAzmAngle_wrt_RP());
		xbin[j][k].push_back(rapid);
		
	      }
	    }
	  }
	}
      }
    }

    Double_t xval [4][nbin];
    Double_t xvale[4][nbin];
    Double_t yval [4][nbin];
    Double_t yvale[4][nbin];

    for(UInt_t j = 0; j < 3; j++) {
      for(UInt_t jn = 0; jn < nbin; jn++) {

	auto getx = vMean(xbin[j][jn]);
	auto gety = vn(hrm, bphi[j][jn]);

	xval [j][jn] = getx[0];
	xvale[j][jn] = getx[1];

	yval [j][jn] = gety[0];
	yvale[j][jn] = gety[1];

      }
    }

    Double_t yratio [2][nbin];
    Double_t yratioe[2][nbin];

    // reversed
    for(UInt_t jn = 0; jn < nbin; jn++) {
      xval [3][jn] = -xval[0][jn];
      xvale[3][jn] = xvale[0][jn];

      yval [3][jn] =  yval[0][jn] * pow(-1, hrm);
      yvale[3][jn] = yvale[0][jn];

      if(yval[0][jn] != 0){
	yratio [0][jn] = abs(yval[1][jn] / yval[0][jn]);
	yratioe[0][jn] = pow(yvale[1][jn],2)/pow(yval[1][jn],2) + pow(yvale[0][jn],2)/pow(yval[0][jn],2) ;
	yratioe[0][jn] = abs(yratio [0][jn]) * sqrt(yratioe[0][jn]);

	// cout << " yratio 0 " << jn << " " << setw(10)
	//      << yval[1][jn] << " / " << yval[0][jn] << " = " << yratio [0][jn]
	//      << " +- " << yratioe[0][jn] << " " << xval[0][jn] <<endl;
      }

      if(yval[0][jn] != 0){
	yratio [1][jn] = abs(yval[2][jn] / yval[0][jn]);
	yratioe[1][jn] = pow(yvale[2][jn],2)/pow(yval[2][jn],2) + pow(yvale[0][jn],2)/pow(yval[0][jn],2) ;
	yratioe[1][jn] = abs(yratio [1][jn]) * sqrt(yratioe[1][jn]);

	// cout << xval[0][jn] << " " <<  xval[2][jn] << endl; 
	// cout << " yratio 1 " << jn << " : " << yratio[1][jn] << " +- " << yratioe[1][jn] << endl;
      }
    }
    

    const UInt_t nsl = 10;
    UInt_t iisl = 2;
    auto gr0 = new TGraphErrors(nsl, &xval[0][iisl], &yval[0][iisl], &xvale[0][iisl], &yvale[0][iisl]);
    gr0->SetName("gr0");
    gr0->SetLineColor(2);
    gr0->SetMarkerStyle(20);
    gr0->SetMarkerColor(2);

    auto gr1 = new TGraphErrors(nsl, &xval[1][iisl], &yval[1][iisl], &xvale[1][iisl], &yvale[1][iisl]);
    gr1->SetName("gr1");
    gr1->SetLineColor(4);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerColor(4);

    auto gr2 = new TGraphErrors(nsl, &xval[2][iisl], &yval[2][iisl], &xvale[2][iisl], &yvale[2][iisl]);
    gr2->SetName("gr2");
    gr2->SetLineColor(8);
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerColor(8);

    auto *rv_gr0 = new TGraphErrors(nsl, &xval[3][iisl], &yval[3][iisl], &xvale[3][iisl], &yvale[3][iisl]);
    rv_gr0->SetName("rv_gv0");
    rv_gr0->SetLineColor(2);
    rv_gr0->SetMarkerStyle(24);
    rv_gr0->SetMarkerColor(2);
    rv_gr0->SetLineStyle(3);


    auto mg = new TMultiGraph();
    TString aTitle = Form("; Ycm/Ycm_beam; v%d(a.u.)",hrm);
    mg->SetTitle(aTitle);
    mg->SetName("mg");
    mg->Add(gr0,"lp");
    mg->Add(gr1,"lp");
    mg->Add(gr2,"lp");
    mg->Add(rv_gr0,"lp");

    if(hrm == 2)
      mg->SetMaximum(0.001);

    // Legend                                                                                                                               

    TString sconf = ""; //Form("%4.2f < y < %4.2f",rapid_min,rapid_max);
    auto aLeg = new TLegend(0.1,0.7,0.35,0.9,sconf);
    aLeg->AddEntry(gr0,partname[selid[0]],"lp");
    aLeg->AddEntry(gr1,partname[selid[1]],"lp");
    aLeg->AddEntry(gr2,partname[selid[2]],"lp");
    aLeg->AddEntry(rv_gr0,partname[selid[0]]+" reversed","lp");

    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    //    gr0->Draw();

    mg  ->Draw("a");
    aLeg->Draw();

    DrawCenterLine(mg);

    iisl = 2;
    auto gt0 = new TGraphErrors(nsl, &xval[0][iisl], &yratio[0][iisl], &xvale[0][iisl], &yratioe[0][iisl]);
    gt0->SetName("gt0");
    gt0->SetLineColor(4);
    gt0->SetMarkerStyle(21);
    gt0->SetMarkerColor(4);

    auto gt1 = new TGraphErrors(nsl, &xval[0][iisl], &yratio[1][iisl], &xvale[0][iisl], &yratioe[1][iisl]);
    gt1->SetName("gt0");
    gt1->SetLineColor(8);
    gt1->SetMarkerStyle(22);
    gt1->SetMarkerColor(8);

    auto mgr = new TMultiGraph();
    aTitle = Form("; Ycm/Ycm_beam; v%d /v%d(proton)",hrm, hrm);
    mgr->SetTitle(aTitle);
    mgr->SetName("mgr");
    mgr->Add(gt0,"lp");
    mgr->Add(gt1,"lp");

    auto aLegr = new TLegend(0.1,0.7,0.35,0.9,"");
    aLegr->AddEntry(gt0,partname[selid[1]],"lp");
    aLegr->AddEntry(gt1,partname[selid[2]],"lp");


    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    mgr->Draw("a");
    aLegr->Draw();


  }
}


void GetRPResolution()            //%% Executable : Plot Phi and subevent, Phi_A and Phi_B correlation
{
  labphi.resize(4);
  subcos.resize(4);

  // Booking
  TH1D *hrpphi[4];
  TH2D *hdltphi[4];
  TH2D *hsubphi[4];
  TGraphErrors *gscos[4];

  for(UInt_t m = m_bgn; m < m_end; m++){

    TString hname = Form("hrpphi%d",m);
    hrpphi[m] = new TH1D(hname,sysName[sys[m]]+";#Phi ",60, -3.2, 3.2);
    hrpphi[m] -> SetLineColor(icol[sys[m]]);

    hname = Form("hdltphi%d",m);
    hdltphi[m] = new TH2D(hname,sysName[sys[m]]+";#Phi ; #Phi_A - #Phi_B",60,-3.2,3.2, 60,-3.2,3.2);
    hdltphi[m]-> SetLineColor(icol[sys[m]]);

    hname = Form("hsubphi%d",m);
    hsubphi[m] = new TH2D(hname,sysName[sys[m]]+";#Phi_A ; Phi_B",60,-3.2,3.2, 60,-3.2,3.2);

  }
    
  // Retreview

  for(Int_t m = m_bgn; m < m_end; m++){
    
    Double_t dphi[nphi];
    Double_t dphie[nphi];
    Double_t scos[nphi];
    Double_t scose[nphi];

    Int_t    snbm;
    Int_t    mtrack;
    Int_t    mtrack_1;
    Int_t    mtrack_2;

    TVector2 *unitP_lang  = NULL;
    TVector2 *unitP_1     = NULL;
    TVector2 *unitP_2     = NULL;

    TBranch  *bunitP_lang;
    TBranch  *bunitP_1;
    TBranch  *bunitP_2;
    TBranch  *brpphi   = 0;
    TBranch  *biphi    = 0;
    TBranch  *bdeltphi = 0;


    labphi[m].resize(nphi);
    subcos[m].resize(nphi);
    cout << " iphi size " << labphi.size() << endl;

    rChain[m]->SetBranchAddress("snbm",&snbm);
    rChain[m]->SetBranchAddress("unitP_lang",&unitP_lang,&bunitP_lang);
    rChain[m]->SetBranchAddress("unitP_1"   ,&unitP_1,&bunitP_1);
    rChain[m]->SetBranchAddress("unitP_2"   ,&unitP_2,&bunitP_2);
    rChain[m]->SetBranchAddress("mtrack"    ,&mtrack);
    rChain[m]->SetBranchAddress("mtrack_1"  ,&mtrack_1);
    rChain[m]->SetBranchAddress("mtrack_2"  ,&mtrack_2);

    Int_t nEntry = rChain[m]->GetEntries();

  
    for(Int_t i = 0; i < nEntry; i++){
      rChain[m]->GetEntry(i);

      if(mtrack_2>0){
	hrpphi[m]  -> Fill(TVector2::Phi_mpi_pi(unitP_lang->Phi()));
	hdltphi[m] -> Fill(TVector2::Phi_mpi_pi(unitP_lang->Phi()), TVector2::Phi_mpi_pi(unitP_1->Phi() - unitP_2->Phi()));
	hsubphi[m] -> Fill(TVector2::Phi_mpi_pi(unitP_1->Phi()), TVector2::Phi_mpi_pi(unitP_2->Phi()) );

	CalculateResolution(m, unitP_lang->Phi(), TVector2::Phi_mpi_pi(unitP_1->Phi() - unitP_2->Phi()) );

      }
    }


    Double_t *getx = new Double_t[2];
    Double_t *gety = new Double_t[2];

    UInt_t j = 0;
    for(UInt_t k = 0; k < nphi; k++){
      getx  = vMean( labphi[m][k] );
      gety  = vMean( subcos[m][k] );
      if(getx[0] > -999) {
	dphi[j] = getx[0];
	dphie[j]= getx[1];
	scos[j] = gety[0];
	scose[j]= gety[1]/sqrt( subcos[m][k].size());
	j++;
      }
    }
    gscos[m] = new TGraphErrors(j, dphi, scos, dphie, scose);
    gscos[m]->SetTitle(sysName[sys[m]]+";#Phi; <cos(#phi_A-#phi_B)>");

    SaveRPResolution(m, nphi, dphi, dphie, scos, scose);
  }


  // Draw
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1500, 400*m_end);
  cc[ic]->Divide(4,m_end);
  
  UInt_t id = 1;
  
  for(Int_t m = m_bgn; m < m_end; m++){
    cc[ic]->cd(id); id++;
    hrpphi[m] -> Draw();


    cc[ic]->cd(id); id++;
    hsubphi[m]-> Draw("colz");

    cc[ic]->cd(id); id++;
    hdltphi[m]-> Draw("colz");

    cc[ic]->cd(id); id++;
    gscos[m]->Draw();
  }  

}

void SaveRPResolution(Int_t m, UInt_t size, Double_t *x, Double_t *xe, Double_t *y, Double_t *ye)
{
  std::fstream fout;


  TString fName = "RP" +  sysName[sys[m]] + ".data";

  gSystem->cd("db");

  fout.open(fName, std::fstream::out);

  for(UInt_t i = 0; i < size; i++)
    fout << std::setw(4) << i 
	 << std::setw(10) << x[i] << " +- " << std::setw(10) << xe[i]
	 << std::setw(10) << y[i] << " +- " << std::setw(10) << ye[i] << std::endl;
    
  gSystem->cd("..");

  fout.close();
}


void CalculateResolution(UInt_t m, Double_t Phi, Double_t phi_sub)
{
  Double_t minphi =  -1.*TMath::Pi();
  Double_t dltphi =  -2.*minphi/nphi;

  Phi = TVector2::Phi_mpi_pi( Phi );


  Double_t iphi = 0;
  for(UInt_t iphi = 0; iphi < (UInt_t)nphi; iphi++){
    
    if( Phi < minphi + dltphi*(iphi+1)) {
      labphi[m][iphi].push_back(Phi);
      subcos[m][iphi].push_back(cos(phi_sub) );
      
      //cout << " m " << m
      //   << " iphi " << iphi << " " << Phi 
      //   << endl;
      
      break;

    
    }
  }
  
}

//________________________________//%%
void FlatteningCheck()            //%% Executable : 
{
  //----- Parametres
  Int_t ntrack[7];

  //----- Booking
  TH2D *hphitheta[4];
  TH2D *hphimtrck[4];

  for(Int_t m = m_bgn; m < m_end; m++){
    TString hname = Form("hphitheta%d",m);
    hphitheta[m] = new TH2D(hname, sysName[sys[m]]+"; #Theta ; #Phi",100,0,0.8, 100,-3.2, 3.2);

    hname = Form("hphimtrck%d",m);
    hphimtrck[m] = new TH2D(hname, sysName[sys[m]]+"; Number of Track ; #Phi",35,0,35, 100,-3.2, 3.2);
  }


  //----- Filling
  for(Int_t m = m_bgn; m < m_end; m++){


    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("ntrack",ntrack);

    Int_t nEntry = rChain[m]->GetEntries();

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
	auto phi   = aPart->GetFlattenMomentum().Phi();
	auto theta = aPart->GetFlattenMomentum().Theta();
	auto flag  = aPart->GetReactionPlaneFlag();

	//	if(flag > 110 ){
	if(flag >= selReactionPlanef ){
	  hphitheta[m]->Fill( theta, phi );
	  hphimtrck[m]->Fill( ntrack[5], phi ); 
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

//________________________________//%% Executable : 
void Phi()                        //%% Executable : 
{
  //----- Parametres
  Int_t    mtrack;
  Int_t    mtrack_1;
  Int_t    mtrack_2;

  TVector2 *unitP_lang  = NULL;
  TVector2 *unitP_1     = NULL;
  TVector2 *unitP_2     = NULL;

  TBranch  *bunitP_lang;
  TBranch  *bunitP_1;
  TBranch  *bunitP_2;
  TBranch  *brpphi   = 0;
  TBranch  *biphi    = 0;
  TBranch  *bdeltphi = 0;

  //----- Booking
  for(Int_t m = m_bgn; m < m_end; m++){
 

  }

  //----- Filling
  for(Int_t m = m_bgn; m < m_end; m++){
    rChain[m]->SetBranchAddress("unitP_lang",&unitP_lang,&bunitP_lang);
    rChain[m]->SetBranchAddress("unitP_1"   ,&unitP_1,&bunitP_1);
    rChain[m]->SetBranchAddress("unitP_2"   ,&unitP_2,&bunitP_2);
    rChain[m]->SetBranchAddress("mtrack"    ,&mtrack);
    rChain[m]->SetBranchAddress("mtrack_1"  ,&mtrack_1);
    rChain[m]->SetBranchAddress("mtrack_2"  ,&mtrack_2);

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

//________________________________//%% Executable :                        
void CorrectedPsi()                        //%%
{
  auto hrt  = new TH1D("hrt" ,"Original; #Psi",60,-3.15,3.15);
  auto hrc  = new TH1D("hrc" ,"ReCentering; #Psi",60,-3.15,3.15);
  auto hfc  = new TH1D("hfc" ,"RC + Shifting; #Psi",60,-3.15,3.15);


  rChain[0]->Project("hrt" , "TVector2::Phi_mpi_pi(unitP2_rot.Phi())");
  rChain[0]->Project("hrc" , "unitP_rc.Phi()");
  rChain[0]->Project("hfc" , "unitP_fc.Phi()");


  //----- Drawing                                                                                                                          
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  hrt->SetLineColor(2);
  hrc ->SetLineColor(8);
  hfc ->SetLineColor(4);

  hrt->Draw("e");
  hrc->Draw("samee");
  hfc->Draw("samee");


  auto aLeg = new TLegend(0.15,0.7,0.45,0.9,"");
  aLeg->AddEntry(hrt ,"No Collection ","lp");
  aLeg->AddEntry(hrc ,"ReCentering","lp");
  aLeg->AddEntry(hfc ,"ReCentering & Shifting","lp");

  aLeg->Draw();
}


//________________________________//%% Executable : 
void PlotSubEvent(Double_t ml, Double_t mu)               //%%
{

  Double_t mlt[] = {ml, mu};
  //    Double_t mlt[2] = {0., 8.};
  //Double_t mlt[2] = {8., 16.};
  // Double_t mlt[2] = {16., 24.};
  // Double_t mlt[2] = {24., 32.};
  // Double_t mlt[2] = {32., 40.};
  // Double_t mlt[2] = {40., 100.};

  TCut mcrot = Form("ntrack[4]>%f&&ntrack[4]<%f",mlt[0],mlt[1]);
  TCut mc1r  = Form("mtrack_1>%f&&mtrack_1<%f",mlt[0],mlt[1]);
  TCut mc2r  = Form("mtrack_2>%f&&mtrack_2<%f",mlt[0],mlt[1]);

  TString sname;
  
  sname = mcrot.GetTitle();
  auto *hrotx = new TH1D("hrotx","All   "+sname+";Qx" ,100,-12.,12.);
  auto *h1rx  = new TH1D("h1rx", "sub_1 "+sname+ ";Qx" ,100,-12.,12.);
  auto *h2rx  = new TH1D("h2rx", "sub_2 "+sname+ ";Qx" ,100,-12.,12.);
  auto *hroty = new TH1D("hroty","All   "+sname+";Qy",100,-12.,12.);
  auto *h1ry  = new TH1D("h1ry", "sub_1 "+sname+ ";Qy" ,100,-12.,12.);
  auto *h2ry  = new TH1D("h2ry", "sub_2 "+sname+ ";Qy" ,100,-12.,12.);
 

  rChain[0]->Project("hrotx","unitP2_rot.X()",mcrot);
  rChain[0]->Project("h1rx","unitP_1r.X()"   ,mc1r);
  rChain[0]->Project("h2rx","unitP_2r.X()"   ,mc2r);
					      
  rChain[0]->Project("hroty","unitP2_rot.Y()",mcrot);
  rChain[0]->Project("h1ry","unitP_1r.Y()"   ,mc1r);
  rChain[0]->Project("h2ry","unitP_2r.Y()"   ,mc2r);


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


//________________________________//%% Executable : 
void PlotNeuLANDPsi()             //%%
{
  TH2D *hACCp[4];
  TH2D *hnlACCp[4];
  TH2D *hnlACCn[4];
  TH1D *hfcrn[4][4];
  TH1D *hfcrp[4][4];

  TCut ncCut[4];
  ncCut[0]="";
  ncCut[1]="ncPID==2112";
  ncCut[2]=ncCut[1]&&"ncRapidity<=0.4";
  ncCut[3]=ncCut[1]&&"ncRapidity>0.4";

  TCut pcCut[4];
  pcCut[0]="";
  pcCut[1]="ncPID==2212";
  pcCut[2]=pcCut[1]&&"ncRapidity<=0.4";
  pcCut[3]=pcCut[1]&&"ncRapidity>0.4";

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
    hACCp[m]->SetMarkerColor(2);
    rChain[m]->Project(hname,"fRotatedP3.Pt():fRapidity","fPID==2212&&fgoodtrackf==1");

    for(UInt_t k = 0; k < 4; k++){
      hname = Form("hfcrn%d%d",m,k);
      hfcrn[m][k] = new TH1D(hname, hname+ncCut[k].GetTitle(), 60, -3.15, 3.15);
      rChain[m]->Project(hname,"unitP_fc.Phi()",ncCut[k]);

      hname = Form("hfcrp%d%d",m,k);
      hfcrp[m][k] = new TH1D(hname, hname+pcCut[k].GetTitle(), 60, -3.15, 3.15);
      rChain[m]->Project(hname,"unitP_fc.Phi()",pcCut[k]);
    }


    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->Divide(2,2);

    cc[ic]->cd(id); id++; 
    hACCp[m]->Draw("colz");
    hnlACCn[m]->Draw("same");
    hnlACCp[m]->Draw("same");

    id++;
    cc[ic]->cd(id); id++;
    hnlACCn[m]->Draw("colz");

    cc[ic]->cd(id); id++;
    hnlACCp[m]->Draw("colz");


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
}



//________________________________//%% Executable : 
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


