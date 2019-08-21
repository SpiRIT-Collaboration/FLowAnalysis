#include "SetStyle.C"

void PlotOMEG2019()
{
  SetStyle();

  TString pname = "proton";
  TString pname1= "triton";

  auto _file0  = TFile::Open("data/advYPt_132Sn_"+pname+".v29.1.24.v2c.root");
  auto _file4  = TFile::Open("data/advYPt_132Sn_"+pname+".v29.1.24.v1c.root");
  auto _file10 = TFile::Open("data/advYPt_132Sn_"+pname1+".v29.1.24.v2c.root");
  auto _file14 = TFile::Open("data/advYPt_132Sn_"+pname1+".v29.1.24.v1c.root");
  auto _file3 = TFile::Open("data/advYPt_132Sn_"+pname+".v29.1.24.root");
  auto _file1 = TFile::Open("../../pBUU/data/pBUU2sn132_sn124_energy270_gamma0.50_b5_withCluster."+pname+".root");
  auto _file2 = TFile::Open("../../pBUU/data/pBUU2sn132_sn124_energy270_gamma1.75_b5_withCluster."+pname+".root");


  auto v2corr  = (TGraphErrors*)_file0 ->Get("v2corr");
  auto v2corr1 = (TGraphErrors*)_file10->Get("v2corr");
  auto gv_v2g5 = (TGraphErrors*)_file1->Get("gv_v2");
  auto gv_v2g1 = (TGraphErrors*)_file2->Get("gv_v2");
  auto gv_v2   = (TGraphErrors*)_file3->Get("gv_v2");

  auto v1corr  = (TGraphErrors*)_file4 ->Get("v1corr");
  auto v1corr1 = (TGraphErrors*)_file14->Get("v1corr");
  auto gv_v1g5 = (TGraphErrors*)_file1->Get("gv_v1");
  auto gv_v1g1 = (TGraphErrors*)_file2->Get("gv_v1");
  auto gv_v1   = (TGraphErrors*)_file3->Get("gv_v1");

  if( v2corr == NULL )
    std::cout << " v2corr is not found." << std::endl;

  gv_v2->SetMarkerStyle(20);
  gv_v2->SetMarkerColor(4);
  gv_v2->SetLineColor(4);

  gv_v2g5->SetMarkerStyle(20);
  gv_v2g5->SetMarkerColor(4);
  gv_v2g5->SetLineColor(4);

  gv_v2g1->SetMarkerStyle(20);
  gv_v2g1->SetMarkerColor(8);
  gv_v2g1->SetLineColor(8);

  v2corr->SetMarkerStyle(20);
  v2corr->SetMarkerColor(2);
  v2corr->SetLineColor(2);

  v2corr1->SetMarkerStyle(20);
  v2corr1->SetMarkerColor(4);
  v2corr1->SetLineColor(4);

  TMultiGraph *mv2 = new TMultiGraph("mv2","; y_{cm}/y_{beam}; v2");
  //mv2->Add(gv_v2, "LP");
  // mv2->Add(gv_v2g5, "LP");
  // mv2->Add(gv_v2g1, "LP");
  mv2->Add(v2corr , "LP");
  mv2->Add(v2corr1, "LP");

  auto leg = new TLegend(0.23, 0.73, 0.54, 0.88,"");
  leg->AddEntry(v2corr,pname+"");
  leg->AddEntry(v2corr1,pname1+"");
  //leg->AddEntry(gv_v2,"w/o corr.");
  // leg->AddEntry(gv_v2g5,"pBUU g0.5");
  // leg->AddEntry(gv_v2g1,"pBUU g1.75");

  TCanvas *cc;
  UInt_t ic = -1;
  ic++;   cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  
  mv2->GetYaxis()->SetRangeUser(-0.1, 0.05);
  mv2->GetXaxis()->SetRangeUser(-1.2,1.2);
  mv2->Draw("ALP");
  leg->Draw();


  //------------ v1
  if( v1corr == NULL )
    std::cout << " v1corr is not found." << std::endl;

  gv_v1->SetMarkerStyle(20);
  gv_v1->SetMarkerColor(4);
  gv_v1->SetLineColor(4);

  gv_v1g5->SetMarkerStyle(20);
  gv_v1g5->SetMarkerColor(4);
  gv_v1g5->SetLineColor(4);

  gv_v1g1->SetMarkerStyle(20);
  gv_v1g1->SetMarkerColor(8);
  gv_v1g1->SetLineColor(8);

  v1corr->SetMarkerStyle(20);
  v1corr->SetMarkerColor(2);
  v1corr->SetLineColor(2);

  v1corr1->SetMarkerStyle(20);
  v1corr1->SetMarkerColor(4);
  v1corr1->SetLineColor(4);

  TMultiGraph *mv1 = new TMultiGraph("mv1","; y_{cm}/y_{beam}; v1");
  //  mv1->Add(gv_v1, "LP");
  // mv1->Add(gv_v1g5, "LP");
  // mv1->Add(gv_v1g1, "LP");
  mv1->Add(v1corr , "LP");
  mv1->Add(v1corr1, "LP");

  auto leg1 = new TLegend(0.23, 0.73, 0.54, 0.88,"");
  leg1->AddEntry(v1corr, pname+"");
  leg1->AddEntry(v1corr1, pname1+"");
  //  leg1->AddEntry(gv_v1,"w/o corr.");
  //  leg1->AddEntry(gv_v1g5,"pBUU g0.5");
  //leg1->AddEntry(gv_v1g1,"pBUU g1.75");

  ic++;   cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

  mv1->GetYaxis()->SetRangeUser(-0.45, 0.45);
  mv1->GetXaxis()->SetRangeUser(-1.2,1.2);
  mv1->Draw("ALP");
  leg1->Draw();

}
