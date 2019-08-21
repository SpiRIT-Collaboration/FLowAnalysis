#include "DoFlow.h"
void DoIntegralv1(UInt_t selid);
void DoIntegralv2(UInt_t selid);

Double_t ptmax = 2000.;
Double_t ptbin = 1000.;


void DoIntegralv(UInt_t selv=2, UInt_t selid=2)
{
  if( selv == 1 )
    DoIntegralv1(selid);
  else
    DoIntegralv2(selid);
}

void DoIntegralv1(UInt_t selid=2)
{
  TCanvas *cc;
  UInt_t  ic = -1;

  
  TString fdata = "data/advYPt_132Sn_"+partname[selid]+".v29.1.24.root";
  std::cout << fdata << std::endl;

  auto fopen = TFile::Open(fdata);
  auto fcorr = TFile::Open("data/correctedPt.phi60.mult5_55.20190625.root");

  TMultiGraph  *v1pt   = new TMultiGraph("v1pt",";pt,v1");
  TGraphErrors *v1corr = new TGraphErrors();
  v1corr->SetName("v1corr");

  auto gv_v1 = (TGraphErrors*)fopen->Get("gv_v1");
  if( gv_v1 == NULL ) exit(0);

  
  UInt_t iyplot = 0;
  for(UInt_t iy = 0; iy < 12; iy++ ) {

    TString v1name = Form("gPt_v1%d",iy);
    TString cchist = Form("h1PtCor_MN1_132Sn_%dH_py%d",(UInt_t)partA[selid],iy);

    auto gr1 = (TGraphErrors*)fopen->Get(v1name);
    auto hpt = (TH1D*)        fcorr->Get(cchist);
  
    if( gr1 == NULL || hpt == NULL || hpt->GetEntries()==0 ) continue;
    v1pt->Add(gr1,"lp");


    TGraph *v1eval      = new TGraph();
    TGraph *v1deriv     = new TGraph();
    v1eval ->SetName(Form("v1eval%d" ,iy));
    v1deriv->SetName(Form("v1deriv%d",iy));

    TH1D *hptc = new TH1D(Form("hptc_%d",iy),"",800,0.,ptmax);
    
    
    Double_t dpt   = ptmax/ptbin;
    Double_t v1tot = 0.;
    Double_t tot   = 0.;

    for(UInt_t ipt = 0; ipt < (UInt_t)ptbin; ipt++) {
      Double_t ptat = dpt*(ipt+1);

      auto v1p = gr1->Eval(ptat); 
      if( ptat > 650 )
	v1p = gr1->Eval(650.);

      auto v1 = v1p * hpt->Interpolate(ptat);
      
      tot   += hpt->Interpolate(ptat);
      v1tot += v1;
      
      hptc->Fill(ptat, v1tot);
      
      v1eval->SetPoint(ipt, ptat, v1p);
      
      if( ptat > 300)
	v1deriv->SetPoint(ipt, ptat, v1tot/tot);
      
    }

    
    Double_t y, x;
    gv_v1->GetPoint(iyplot, x, y);

    Double_t v1c = 0.;
    if( x <= 0 )
      v1c = v1deriv->GetYaxis()->GetXmin();
    else
      v1c = v1deriv->GetYaxis()->GetXmax();

    v1corr->SetPoint(iyplot, x, v1c);
    iyplot++;


    std::cout << iy << " -> " << " original " << x << " : " << y ;
    std::cout << " average " << v1tot/tot << " = " << v1tot << " / " << tot << endl;
    std::cout << " min  " << v1c << endl;


    UInt_t id = 1;
    ic++;   cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1000);
    cc->Divide(1,4);

    cc->cd(id); id++;
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColor(4);
    gr1->SetLineColor(4);
    gr1->Draw("ALP");
    v1deriv->Draw("same");

    cc->cd(id); id++;
    hpt->Draw("");

    cc->cd(id); id++;
    v1deriv->Draw("ALP");


    cc->cd(id); 
    hptc->Draw();

  }

  ic++;   cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
  v1corr->SetMarkerStyle(20);
  v1corr->SetMarkerColor(4);
  v1corr->SetLineColor(4);
  v1corr->Draw("ALP");

  ic++;   cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  v1pt->Draw("ALP");


  fdata.ReplaceAll(".root",".v1c.root");
  auto fout = TFile::Open(fdata,"recreate");
  v1corr->Write();

}


void DoIntegralv2(UInt_t selid=2)
{
  TCanvas *cc;
  UInt_t  ic = -1;

  TString fdata = "data/advYPt_132Sn_"+partname[selid]+".v29.1.24.root";
  auto fopen = TFile::Open(fdata);
  auto fcorr = TFile::Open("data/correctedPt.phi60.mult5_55.20190625.root");

  TMultiGraph  *v2pt   = new TMultiGraph("v2pt",";pt,v2");
  TGraphErrors *v2corr = new TGraphErrors();
  v2corr->SetName("v2corr");

  auto gv_v2 = (TGraphErrors*)fopen->Get("gv_v2");
  if( gv_v2 == NULL ) exit(0);

  
  UInt_t iyplot = 0;
  for(UInt_t iy = 0; iy < 8; iy++ ) {

    TString v2name = Form("gPt_v2%d",iy);
    TString cchist = Form("h1PtCor_MN2_132Sn_%dH_py%d",(UInt_t)partA[selid],iy);

    auto gr1 = (TGraphErrors*)fopen->Get(v2name);
    auto hpt = (TH1D*)        fcorr->Get(cchist);
  
    if( gr1 == NULL || hpt == NULL || hpt->GetEntries()==0 ) continue;
    //    gr1->SetMarkerSize(1.);
    //    gr1->SetMarkerColor(kRed+iy);
    //    gr1->SetLineColor(kRed+iy);
    v2pt->Add(gr1,"lp");



    TGraph *v2eval      = new TGraph();
    TGraph *v2deriv     = new TGraph();
    v2eval ->SetName(Form("v2eval%d" ,iy));
    v2deriv->SetName(Form("v2deriv%d",iy));

    TH1D *hptc = new TH1D(Form("hptc_%d",iy),"",800,0.,ptmax);
    
    
    Double_t dpt   = ptmax/ptbin;
    Double_t v1tot = 0.;
    Double_t tot   = 0.;

    for(UInt_t ipt = 0; ipt < (UInt_t)ptbin; ipt++) {
      Double_t ptat = dpt*(ipt+1);

      auto v1p = gr1->Eval(ptat); 
      if( ptat > 650 )
	v1p = gr1->Eval(650.);

      auto v1 = v1p * hpt->Interpolate(ptat);
      
      tot   += hpt->Interpolate(ptat);
      v1tot += v1;
      
      hptc->Fill(ptat, v1tot);
      
      v2eval->SetPoint(ipt, ptat, v1p);
      
      if( ptat > 300)
	v2deriv->SetPoint(ipt, ptat, v1tot/tot);
      
    }

    auto minv2 = v2deriv->GetYaxis()->GetXmin();
  
    
    Double_t y, x;
    gv_v2->GetPoint(iyplot, x, y);
    v2corr->SetPoint(iyplot, x, minv2);
    iyplot++;


    std::cout << iy << " -> " << " original " << x << " : " << y ;
    std::cout << " average " << v1tot/tot << " = " << v1tot << " / " << tot << endl;
    std::cout << " min  " << minv2 << endl;


    UInt_t id = 1;
    ic++;   cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1000);
    cc->Divide(1,4);

    cc->cd(id); id++;
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColor(4);
    gr1->SetLineColor(4);
    gr1->Draw("ALP");
    v2deriv->Draw("same");

    cc->cd(id); id++;
    hpt->Draw("");

    cc->cd(id); id++;
    v2deriv->Draw("ALP");


    cc->cd(id); 
    hptc->Draw();

  }

  ic++;   cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
  v2corr->SetMarkerStyle(20);
  v2corr->SetMarkerColor(4);
  v2corr->SetLineColor(4);
  v2corr->Draw("ALP");

  ic++;   cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  v2pt->Draw("ALP");

  fdata.ReplaceAll(".root",".v2c.root");
  auto fout = TFile::Open(fdata,"recreate");
  //  gROOT->cd();
  v2corr->Write();

}
