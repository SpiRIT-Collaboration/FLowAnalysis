#include "STFlowCorrection.hh"
#include "FlowFunctions.h"
#include "openRComp.C"

TCanvas *cc[12];

//const Double_t ycm      = 0.388568; // 132Sn + 124Sn before Nov.15 201 //278MeV/u
//                         0: 132Sn + 124Sn,  1: 108Sn + 112Sn  //270MeV/u
const Double_t ycm[]      = {0.383198,  0.365587}; 
const Double_t ybeam_cm[] = {0.36599 ,  0.378516};
const Double_t ybeam_lb[] = {0.749188,  0.744103};
UInt_t sys[2]   = {0, 0};

//TString sysName[2] = {"132Sn+124Sn","108Sn+112Sn"};
TString sysName[2] = {"132Sn","108Sn"};


const UInt_t nspec = 5;
const UInt_t  nbin = 16;

//UInt_y   y_nbin  = 100;
Double_t y_min[] = {0., 0., 0., 0., 0.};
Double_t y_max[] = {1.8, 1.8, 1.2, 1., 0.8};
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
UInt_t   icol[2]    = {4,2};
TString  iopt[2]    = {"","same"};
  
UInt_t ic = 0;

const Int_t nbinx = 30;

TChain *rChain[2];
UInt_t m_bgn = 0;
UInt_t m_end = 1;

// histogram
TH1D* hptpm[2][nbin];
TH1D* hptpp[2][nbin];
TH1D* hptpr[2][nbin];
TH1D* hptdt[2][nbin];
TH1D* hpttr[2][nbin];

TMultiGraph *mg;
TMultiGraph *mgr;


Int_t mtrack;
auto aArray = new TClonesArray("STParticle",100);

// pt dependence

Double_t rapid_max = 0.4;
Double_t rapid_min = 0.2;
UInt_t pxbooking();
UInt_t DrawCenterLine();
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

void calcRComp()
{
  gROOT->Reset();

  UInt_t ichain = 0;

  //  gROOT->Macro("openRComp.C");
  openRComp();

  rChain[ichain] = (TChain*)gROOT->FindObject(Form("rChain%d",ichain));
  if(rChain[ichain] != NULL) {    
    sys[ichain] = GetSystem(ichain);
    std::cout << " System" << ichain << " "  << sysName[sys[ichain]] << std::endl; 
    ichain++;
  }
  

  rChain[ichain] = (TChain*)gROOT->FindObject(Form("rChain%d",ichain));
  if(rChain[ichain] == NULL) 
    ichain = 0;
  else{
    sys[ichain] = GetSystem(ichain);
    std::cout << " System" << ichain << " "  << sysName[sys[ichain]] << std::endl; 
  }

  if(rChain[0] == NULL && rChain[1] == NULL)
    exit(0);

  


  cout << " ichain " << ichain <<endl;

  m_end = ichain+1;
  
  gROOT->ProcessLine(".! grep -i void calcRComp.C | grep '//%%'");

  //  meanPx();
  //  dndy();
}


void dndy()                       //%% Executable : Make plots of dNdy for p, d, t, pi+-
{

  //----- booking
  TH1D* hrap[2][5];
  TH1D* hnpart[2][5];
  TH1D* hmtrack[2];

  for(Int_t m = m_bgn; m < m_end; m++){

    for(UInt_t i = 0; i < 5; i++){

      UInt_t y_nbin = UInt_t((y_max[i] - y_min[i])/y_bin[i]);
      TString hname = Form("hrap%d_%d",m,i);
      TString htitle= partname[i] + "; Rapidity; dN/dy";
      hrap[m][i] = new TH1D(hname, htitle, y_nbin, y_min[i], y_max[i]);
      hrap[m][i] ->SetLineColor(icol[m]);

      hname  = Form("hnpart%d_%d",m,i);
      htitle = partname[i]+" ; Multiplicity";
      hnpart[m][i] = new TH1D(hname, htitle, 15,0,15);
      hnpart[m][i]->SetLineColor(icol[m]);
    }

    hmtrack[m] = new TH1D(Form("hmtrack%d",m),"Number of good tracks; Multiplicity",50, 0, 50);
    hmtrack[m] -> SetLineColor(icol[m]);
    
  }

  //------------------------

  for(Int_t m = m_bgn; m < m_end; m++){

    aArray->Clear();
    mtrack = 0;

    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("mtrack" ,&mtrack);

    Int_t nEntry = rChain[m]->GetEntries();
    
    for(Int_t i = 0; i < nEntry; i++){
      rChain[m]->GetEntry(i);

      TIter next(aArray);
      STParticle *aPart = NULL;

      hmtrack[m]->Fill(mtrack);
      UInt_t npart[5] = {0,0,0,0,0};

      while( (aPart = (STParticle*)next()) ) {

	auto pid   = aPart->GetPID();
	auto charg = aPart->GetCharge();
	auto p     = aPart->GetRotatedMomentum();
	auto rapid = aPart->GetRapidity();
	auto theta = aPart->GetRotatedMomentum().Theta();

	if(theta < 0.8){
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
  ic = 0;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,800);
  cc[ic]->Divide(3,2);

  //  auto aLeg0 = new TLegend(0.7,0.7,0.9,0.9,"");
  auto aLeg0 = new TLegend(0.1,0.7,0.35,0.9,"");
  if(hrap[0][4] != NULL) aLeg0->AddEntry(hrap[0][4],sysName[0],"lp");
  if(hrap[1][4] != NULL) aLeg0->AddEntry(hrap[1][4],sysName[1],"lp");
    
  id = 1;
  cc[ic]->cd(1); 
  UInt_t io = 0;
  if(hrap[1][2] != NULL) {hrap[1][2]->Draw(iopt[io]); io++;} 
  if(hrap[0][2] != NULL) {hrap[0][2]->Draw(iopt[io]); io++;}


  cc[ic]->cd(2); 
  io = 0;
  if(hrap[0][3] != NULL) {hrap[0][3]->Draw(iopt[io]); io++;}
  if(hrap[1][3] != NULL) {hrap[1][3]->Draw(iopt[io]); io++;} 

  cc[ic]->cd(3); 
  io = 0;
  if(hrap[0][4] != NULL) {hrap[0][4]->Draw(iopt[io]); io++;}
  if(hrap[1][4] != NULL) {hrap[1][4]->Draw(iopt[io]); io++;}

  aLeg0->Draw();

  cc[ic]->cd(4);
  io = 0;
  if(hrap[0][0] != NULL) {hrap[0][0]->Draw(iopt[io]); io++;}
  if(hrap[1][0] != NULL) {hrap[1][0]->Draw(iopt[io]); io++;}

  cc[ic]->cd(5);
  io = 0;
  if(hrap[1][1] != NULL) {hrap[1][1]->Draw(iopt[io]); io++;}
  if(hrap[0][1] != NULL) {hrap[0][1]->Draw(iopt[io]); io++;}

  //----cc1
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,800);
  cc[ic]->Divide(3,2);
  // 
  auto aLeg1 = new TLegend(0.7,0.7,0.9,0.9,"");
  aLeg1->AddEntry(hnpart[0][4],sysName[0],"lp");
  aLeg1->AddEntry(hnpart[1][4],sysName[1],"lp");

  cc[ic]->cd(1); 
  hnpart[0][2]->Draw(iopt[0]);
  hnpart[1][2]->Draw(iopt[1]);

  cc[ic]->cd(2);
  hnpart[0][3]->Draw(iopt[0]);
  hnpart[1][3]->Draw(iopt[1]);

  cc[ic]->cd(3); 
  hnpart[1][4]->Draw(iopt[0]);
  hnpart[0][4]->Draw(iopt[1]);
  aLeg1->Draw();

  cc[ic]->cd(4);
  cc[ic]->cd(4)->SetLogy();
  hnpart[0][0]->Draw(iopt[0]);
  hnpart[1][0]->Draw(iopt[1]);

  cc[ic]->cd(5);
  cc[ic]->cd(5)->SetLogy();
  hnpart[0][1]->Draw(iopt[0]);
  hnpart[1][1]->Draw(iopt[1]);

  cc[ic]->cd(6);
  hmtrack[0]->Draw();
  hmtrack[1]->Draw(iopt[1]);

}

//________________________________//%%
void meanPx()                     //%% Executable : Make  plots of <px> vs rapidity 
{
  PxDistribution(1);

  TGraphErrors  *gpr[2];
  TGraphErrors  *gdt[2];
  TGraphErrors  *gtr[2];
  TGraphErrors  *gpm[2];
  TGraphErrors  *gpp[2];

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
   
    ic = 5;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,800);
      cc[ic]->Divide(3,2);
    }

    UInt_t id = 1;

    cc[ic]->cd(id);id++;
    TString atitle = Form("gpr%d",m);
    gpr[m]->SetName(atitle);
    gpr[m]->SetLineColor(icol[m]);
    gpr[m]->SetMarkerStyle(20+m);
    gpr[m]->SetMarkerColor(icol[m]);
    gpr[m]->SetTitle("Proton; y_lab; <Px> (MeV/c)");
    gpr[m]->Draw(iopt[m]);


    cc[ic]->cd(id);id++;
    atitle = Form("gdt%d",m);
    gdt[m]->SetName(atitle);
    gdt[m]->SetLineColor(icol[m]);
    gdt[m]->SetMarkerStyle(20+m);
    gdt[m]->SetMarkerColor(icol[m]);
    gdt[m]->SetTitle("Deuteron; y_lab; <Px> (MeV/c)");
    gdt[m]->Draw(iopt[m]);

    cc[ic]->cd(id);id++;
    atitle = Form("gtr%d",m);
    gtr[m]->SetName(atitle);
    gtr[m]->SetLineColor(icol[m]);
    gtr[m]->SetMarkerStyle(20+m);
    gtr[m]->SetMarkerColor(icol[m]);
    gtr[m]->SetTitle("Triton; y_lab; <Px> (MeV/c)");
    gtr[m]->Draw(iopt[m]);

    cc[ic]->cd(id);id++;
    atitle = Form("gpm%d",m);
    gpm[m]->SetName(atitle);
    gpm[m]->SetLineColor(icol[m]);
    gpm[m]->SetMarkerStyle(20+m);
    gpm[m]->SetMarkerColor(icol[m]);
    gpm[m]->SetTitle("pi-; y_lab; <Px> (MeV/c)");
    gpm[m]->Draw(iopt[m]);

    cc[ic]->cd(id);id++;
    atitle = Form("gpp%d",m);
    gpp[m]->SetName(atitle);
    gpp[m]->SetLineColor(icol[m]);
    gpp[m]->SetMarkerStyle(20+m);
    gpp[m]->SetMarkerColor(icol[m]);
    gpp[m]->SetTitle("pi+; y_lab; <Px> (MeV/c)");
    gpp[m]->Draw(iopt[m]);
  }


  mg = new TMultiGraph();
  mg->SetTitle("; Rapidity_lab; <px>/A [MeV/c]");
  gpr[0]->SetLineColor(2);
  gpr[0]->SetMarkerStyle(20);
  gpr[0]->SetMarkerColor(2);

  gpr[1]->SetLineColor(2);
  gpr[1]->SetMarkerStyle(24);
  gpr[1]->SetMarkerColor(2);
  gpr[1]->SetLineStyle(3);

  gdt[0]->SetLineColor(4);
  gdt[0]->SetMarkerStyle(21);
  gdt[0]->SetMarkerColor(4);
  gdt[1]->SetLineColor(4);
  gdt[1]->SetMarkerStyle(25);
  gdt[1]->SetMarkerColor(4);
  gdt[1]->SetLineStyle(3);

  gtr[0]->SetLineColor(8);
  gtr[0]->SetMarkerStyle(22);
  gtr[0]->SetMarkerColor(8);
  gtr[1]->SetLineColor(8);
  gtr[1]->SetMarkerStyle(26);
  gtr[1]->SetMarkerColor(8);
  gtr[1]->SetLineStyle(3);
  
  mg->Add(gpr[0],"lp");
  mg->Add(gdt[0],"lp");
  mg->Add(gtr[0],"lp");
  mg->Add(gpr[1],"lp");
  mg->Add(gdt[1],"lp");
  mg->Add(gtr[1],"lp");


  auto aLeg = new TLegend(0.1,0.7,0.35,0.9,"");
  aLeg->AddEntry(gpr[0],"proton   "+sysName[0],"lp");
  aLeg->AddEntry(gdt[0],"deuteron "+sysName[0],"lp");
  aLeg->AddEntry(gtr[0],"trition  "+sysName[0],"lp");
  aLeg->AddEntry(gpr[1],"proton   "+sysName[1],"lp");
  aLeg->AddEntry(gdt[1],"deuteron "+sysName[1],"lp");
  aLeg->AddEntry(gtr[1],"trition  "+sysName[1],"lp");

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mg  ->Draw("a");
  aLeg->Draw();

}

UInt_t pxbooking()
{
  //----- booking

  for(Int_t m = m_bgn; m < m_end; m++){

    for(UInt_t i = 0; i < nbin; i++){

      Double_t yL = y_min[0] + y_binx*i;
      Double_t yU = yL + y_binx;

      TString hname = Form("hptpr%d%d",m,i);
      TString htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpr[m][i] = new TH1D(hname,"Proton :  "+htitle,pt_nbin, pt_prmin, pt_prmax);
      hptpr[m][i] ->SetLineColor(icol[m]);

      hname = Form("hptdt%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptdt[m][i] = new TH1D(hname,"Deuteron :"+htitle,pt_nbin, pt_dtmin, pt_dtmax);
      hptdt[m][i] ->SetLineColor(icol[m]);

      hname = Form("hpttr%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hpttr[m][i] = new TH1D(hname,"Triton :  "+htitle,pt_nbin, pt_trmin, pt_trmax);
      hpttr[m][i] ->SetLineColor(icol[m]);

      hname = Form("hptpm%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpm[m][i] = new TH1D(hname,"Pi- :     "+htitle,pt_nbin, pt_pimin, pt_pimax);
      hptpm[m][i] ->SetLineColor(icol[m]);

      hname = Form("hptpp%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpp[m][i] = new TH1D(hname,"Pi+ :     "+htitle,pt_nbin, pt_pimin, pt_pimax);
      hptpp[m][i] ->SetLineColor(icol[m]);

    }
  }
  return 1;
}



void PxDistribution(UInt_t nplot)
{
  pt_prmin  =  -pt_prmax;
  pt_dtmin  =  -pt_dtmax;
  pt_trmin  =  -pt_trmax;
  pt_pimin  =  -pt_pimax;

  pxbooking();
  
  for(Int_t m = m_bgn; m < m_end; m++){

    aArray->Clear();
    mtrack = 0;
    
    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("mtrack" ,&mtrack);


    Int_t nEntry = rChain[m]->GetEntries();
    cout << " Number of events " << nEntry << endl;


    for(Int_t i = 0; i < nEntry; i++){
      rChain[m]->GetEntry(i);

      if(mtrack == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto pid   = aPart->GetPID();
	auto charg = aPart->GetCharge();
	auto rapid = aPart->GetRapidity();
	auto vp    = aPart->GetFlattenMomentum();
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







  //old-------------------------  
void PhiYbinf() // old
{

  TH1D *hphi_prt[nbin];
  TH1D *hphi_dtr[nbin];
  TH1D *hphi_trt[nbin];

  TString hname;
  TString htitle;
  for(UInt_t i = 0; i < nbin; i++){
    htitle = Form("%f5.2 < ycm/ycm_beam < %f5.2 ;[Rad]",y_min[0]+i*y_binx, y_min[0]+(i+1.)*y_binx); 

    hname  = Form("hphi_prt%d",i);
    hphi_prt[i] = new TH1D(hname,"Proton"+htitle ,15, 0.,3.1);
    hname  = Form("hphi_dtr%d",i);
    hphi_dtr[i] = new TH1D(hname,"Deutron"+htitle,15, 0.,3.1);
    hname  = Form("hphi_trt%d",i);
    hphi_trt[i] = new TH1D(hname,"Triton"+htitle ,15, 0.,3.1);
  }
  Double_t nfct = 15./3.1;


  Int_t mtrack;
  auto aArray = new TClonesArray("STParticle",100);

  for(Int_t m = m_bgn; m < m_end; m++){

    aArray->Clear();
    
    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("mtrack" ,&mtrack);


    Int_t nevt = rChain[m]->GetEntries();
    cout << " Number of events " << nevt << endl;


    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      if(mtrack == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto pid = aPart->GetPID();
	auto rapid = aPart-> GetRapidity();
	rapid = (rapid - ycm[sys[m]]) / ybeam_cm[sys[m]];

	for(UInt_t k = 0; k < nbin; k++){
	  
	  if(rapid < y_min[0] + y_binx * (Double_t)(k+1)) {
	    if(pid == partid[2])
	      hphi_prt[k]->Fill(abs(aPart->GetAzmAngle_wrt_RP()));
	    else if(pid == partid[3])
	      hphi_dtr[k]->Fill(abs(aPart->GetAzmAngle_wrt_RP()));
	    else if(pid == partid[4])
	      hphi_trt[k]->Fill(abs(aPart->GetAzmAngle_wrt_RP()));

	    break;
	    
	  }
	}
      }
    }   

    gStyle -> SetTitleFontSize(0.1);

    UInt_t ic = 0;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,1200);
    cc[ic]->Divide(3,16);

    UInt_t id = 1;
    for(UInt_t i = 0; i < nbin; i++){
      cc[ic]->cd(id); id++; 
      auto N = hphi_prt[i]->GetEntries();
      hphi_prt[i]->Scale(nfct/(Double_t)N);
      hphi_prt[i]->Draw();

      cc[ic]->cd(id); id++; 
      hphi_dtr[i]->Scale(nfct/(Double_t)N);
      hphi_dtr[i]->Draw();

      cc[ic]->cd(id); id++; 
      hphi_trt[i]->Scale(nfct/(Double_t)N);
      hphi_trt[i]->Draw();
    }
  }
}

void PhiYbin() // old
{

  TH1D *hphi_prt[nbin/2];
  TH1D *hphi_dtr[nbin/2];
  TH1D *hphi_trt[nbin/2];

  TString hname;
  TString htitle;
  for(UInt_t i = 0; i < nbin/2; i++){
    htitle = Form("%f5.2 < ycm/ycm_beam < %f5.2 ;[Rad]",y_min[0]+i*y_binx*2,y_min[0]+(i+1.)*y_binx*2); 

    hname  = Form("hphi_prt%d",i);
    hphi_prt[i] = new TH1D(hname,"Proton"+htitle ,15, 0.,3.1);
    hname  = Form("hphi_dtr%d",i);
    hphi_dtr[i] = new TH1D(hname,"Deutron"+htitle,15, 0.,3.1);
    hname  = Form("hphi_trt%d",i);
    hphi_trt[i] = new TH1D(hname,"Triton"+htitle ,15, 0.,3.1);
  }
  Double_t nfct = 15./3.1;


  Int_t mtrack;
  auto aArray = new TClonesArray("STParticle",100);

  for(Int_t m = m_bgn; m < m_end; m++){

    aArray->Clear();
    
    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("mtrack" ,&mtrack);


    Int_t nevt = rChain[m]->GetEntries();
    cout << " Number of events " << nevt << endl;


    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      if(mtrack == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto pid = aPart->GetPID();
	auto rapid = aPart-> GetRapidity();
	rapid = (rapid - ycm[sys[m]]) / ybeam_cm[sys[m]];

	for(UInt_t k = 0; k < nbin/2; k++){
	  
	  if(rapid < y_min[0] + y_binx*2.*(Double_t)(k+1)) {

	    if(pid == partid[2])
	      hphi_prt[k]->Fill(abs(aPart->GetAzmAngle_wrt_RP()));
	    else if(pid == partid[3])
	      hphi_dtr[k]->Fill(abs(aPart->GetAzmAngle_wrt_RP()));
	    else if(pid == partid[4])
	      hphi_trt[k]->Fill(abs(aPart->GetAzmAngle_wrt_RP()));

	    break;
	    
	  }
	}
      }
    }   

    gStyle -> SetTitleFontSize(0.1);

    UInt_t ic = 0;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,1200);
    cc[ic]->Divide(3,8);

    
    UInt_t id = 1;
    for(UInt_t i = 0; i < 8; i++){
      cc[ic]->cd(id); id++; 
      auto N = hphi_prt[i]->GetEntries();
      hphi_prt[i]->Scale(nfct/(Double_t)N);
      hphi_prt[i]->Draw();

      cc[ic]->cd(id); id++; 
      hphi_dtr[i]->Scale(nfct/(Double_t)N);
      hphi_dtr[i]->Draw();

      cc[ic]->cd(id); id++; 
      hphi_trt[i]->Scale(nfct/(Double_t)N);
      hphi_trt[i]->Draw();
    }
  }
}

void plotv1v2(UInt_t selid=2)     //%% Executable :   v1 and v2 as a function of rapidity
{
  if(selid > 4 ) return;

  gStyle->SetGridColor(7);
  gStyle->SetGridStyle(1);


  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  Int_t RPflag = 0;
  if(selid >= 2) RPflag = 10;


  cout << " Particle " << partname[selid] << endl;

  auto haccp = new TH2D("haccp",partname[selid]+";ycm/ycm_beam ; Pt [MeV/c]",100,-1.1, 1.5,100,0.,800.);
  auto haccx = new TH2D("haccx",partname[selid]+";ycm/ycm_beam ; Px [MeV/c]",100,-1.1, 1.5,100,-500,500.);


  Double_t rp_res[2] = {0.107183, 0.0003226};
  Double_t rp_rese[2]= {0.69 ,0.71};



  Int_t ntrack[7];
  Int_t mtrack;
  
  auto aArray = new TClonesArray("STParticle",100);

  for(Int_t m = m_bgn; m < m_end; m++){

    auto ymin = (y_min[0] - ycm[sys[m]])/ybeam_cm[sys[m]];

    aArray->Clear();
    
    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("ntrack",ntrack);
    rChain[m]->SetBranchAddress("mtrack" ,&mtrack);

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

      if(mtrack == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto pid = aPart->GetPID();

	if(pid == partid[selid] && aPart->GetCharge() == pcharge ){

	  auto rapid = aPart-> GetRapidity();

	  rapid = (rapid - ycm[sys[m]]) / ybeam_cm[sys[m]];

	  haccp -> Fill(rapid, aPart->GetRotatedMomentum().Pt());
	  haccx -> Fill(rapid, aPart->GetFlattenMomentum().X());

	  
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
    UInt_t jnfirst = 0;
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
	if(jnfirst == 0)
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
		<< xval[i] << " +-" <<  " v1 " << setw(12) << yval1[i] << " +- " << yval1e[i] << "  w  " << bphi[i].size() << std::endl;  
    std::cout << " ----------------------------------" << std::endl;
    for(UInt_t i = 0; i < nbin; i++)
      std::cout << setw(5) << i << " : " << setw(12)
		<< xval[i] << " +-" <<  " v2 " << setw(12) << yval2[i] << " +- " << yval2e[i] << "  w  " << bphi[i].size() << std::endl;  

    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    ic++;


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

    //cc[ic]->cd(id); id++;
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


    
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    ic++;

    
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
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);


    //    cc[ic]->Divide(1,2);
    
    //    cc[ic]->cd(id); id++;
    haccp -> Draw("colz");
    
    //    cc[ic]->cd(id); 
    //    haccx -> Draw("colz");

    ic++;
  }
}

UInt_t DrawCenterLine()
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


void RPresolution()
{

  UInt_t   nbin = 10;
  Double_t dbin = 6.4/(Double_t)nbin;

  Int_t    mtrack;
  Int_t    mtrack_1;
  Int_t    mtrack_2;
  Int_t ntrack[7];
  Double_t aoq;
  Double_t z;

  TVector2 *unitP_lang  =NULL;
  TVector2 *unitP_1=NULL;
  TVector2 *unitP_2=NULL;

  TBranch  *bunitP_lang;
  TBranch  *bunitP_1;
  TBranch  *bunitP_2;
  TBranch *brpphi=0;
  TBranch *biphi=0;
  TBranch *bdeltphi=0;


  for(Int_t m = m_bgn; m < m_end; m++){

    rChain[m]->SetBranchAddress("ntrack",ntrack);
    rChain[m]->SetBranchAddress("aoq",&aoq);
    rChain[m]->SetBranchAddress("z",&z);

    rChain[m]->SetBranchAddress("unitP_lang",&unitP_lang,&bunitP_lang);
    rChain[m]->SetBranchAddress("unitP_1"   ,&unitP_1,&bunitP_1);
    rChain[m]->SetBranchAddress("unitP_2"   ,&unitP_2,&bunitP_2);
    rChain[m]->SetBranchAddress("mtrack"    ,&mtrack);
    rChain[m]->SetBranchAddress("mtrack_1"  ,&mtrack_1);
    rChain[m]->SetBranchAddress("mtrack_2"  ,&mtrack_2);



    auto nevt = rChain[m]->GetEntries();


    vector< Double_t > alldphi;
    vector< vector<Double_t> > dphi(nbin);
    vector< vector<Double_t> > xphi(nbin);
    
    for(auto i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      if(mtrack_1 < 1 || mtrack_2 < 1) continue;

      alldphi.push_back( TVector2::Phi_mpi_pi( abs( unitP_1->Phi() - unitP_2->Phi() ) ) );

      for(auto k = 0; k < nbin; k++){
	if( unitP_lang->Phi() < dbin*(Double_t)(k+1) ){
	  dphi[k].push_back( TVector2::Phi_mpi_pi( abs( unitP_1->Phi() - unitP_2->Phi() ) ) );
	  xphi[k].push_back( unitP_lang->Phi()  );
	  break;
	}
      }
    }

    Double_t *allcos = new Double_t[2];
    allcos = vn(1, alldphi);
    
    cout << " Resolution " << allcos[0]/sqrt(2.) << " +- " << allcos[1] << endl;


    Double_t xval[nbin];
    Double_t xvale[nbin];
    Double_t yval1[nbin];
    Double_t yval1e[nbin];
    Double_t yval2[nbin];
    Double_t yval2e[nbin];

    for(UInt_t j = 0; j < nbin; j++){
      Double_t *getx = new Double_t[2];
      Double_t *gety = new Double_t[2];

      getx   = vMean(xphi[j]);
      gety   = vn(1, dphi[j]);

      xval[j]  = getx[0];
      yval1[j] = gety[0]/sqrt(2.);

      xvale[j]  = getx[1];
      yval1e[j] = gety[1];

    }


    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    ic++;

    auto *gv_v1 = new TGraphErrors(nbin, xval, yval1, xvale, yval1e);
    gv_v1->SetName("gv_v1");


    gv_v1->SetLineColor(4);
    gv_v1->SetMarkerStyle(20);
    gv_v1->SetMarkerColor(4);
    gv_v1->SetTitle("Subevents correlation ; #Phi ; 1/sqrt(2)*< cos(#Delta#Phi) >");
   
    gv_v1->Draw("ALP");

  }
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
  Int_t ntrack[7];
  Int_t mtrack;
  auto aArray = new TClonesArray("STParticle",100);

  for(Int_t m = m_bgn; m < m_end; m++){

    aArray->Clear();
    
    rChain[m]->SetBranchAddress("STParticle",&aArray);
    rChain[m]->SetBranchAddress("ntrack",ntrack);
    rChain[m]->SetBranchAddress("aoq",&aoq);
    rChain[m]->SetBranchAddress("z",&z);
    rChain[m]->SetBranchAddress("mtrack" ,&mtrack);

    vector< vector<Double_t> > bphi  (nbin);
    vector< vector<Double_t> > rpdbin(nbin);



    Int_t nevt = rChain[m]->GetEntries();
    cout << " Number of events " << nevt << endl;


    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      if(mtrack == 0) continue;

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

  auto aArray = new TClonesArray("STParticle",100);

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


    mg = new TMultiGraph();
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

  auto aArray = new TClonesArray("STParticle",100);

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


    mg = new TMultiGraph();
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

    DrawCenterLine();

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

    mgr = new TMultiGraph();
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

