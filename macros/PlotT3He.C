{
  TFile *_file10 = TFile::Open("data/advYPt_132Sn_triton.v37.1.3.root");
  TFile *_file11 = TFile::Open("data/advYPt_132Sn_3He.v37.1.3.root");

  TFile *_file20 = TFile::Open("data/advYPt_108Sn_triton.v37.1.3.root");
  TFile *_file21 = TFile::Open("data/advYPt_108Sn_3He.v37.1.3.root");


  auto gv1_132trit = (TGraphErrors*)_file10->Get("gu_v2");
  auto gv1_1323He  = (TGraphErrors*)_file11->Get("gu_v2");

  auto gv1_108trit = (TGraphErrors*)_file20->Get("gu_v2");
  auto gv1_1083He  = (TGraphErrors*)_file21->Get("gu_v2");


  auto gv1_132sub = new TGraphErrors();
  gv1_132sub->SetName("gv1_132Tsub3He");

  auto gv1_108sub = new TGraphErrors();
  gv1_108sub->SetName("gv1_108Tsub3He");

  Double_t xt,yt,xh,yh;
  for(UInt_t ix = 0; ix < 7; ix++){
    
    gv1_132trit->GetPoint(ix, xt, yt);
    gv1_1323He ->GetPoint(ix, xh, yh);
    Double_t yte = gv1_132trit->GetErrorY(ix);
    Double_t yhe = gv1_1323He ->GetErrorY(ix);
    Double_t ye  = sqrt(yte*yte + yhe*yhe);

    gv1_132sub->SetPoint(ix, xt, yt-yh);
    gv1_132sub->SetPointError(ix, 0, ye);


    gv1_108trit->GetPoint(ix, xt, yt);
    gv1_1083He ->GetPoint(ix, xh, yh);
    yte = gv1_108trit->GetErrorY(ix);
    yhe = gv1_1083He ->GetErrorY(ix);
    ye  = sqrt(yte*yte + yhe*yhe);

    gv1_108sub->SetPoint(ix, xt, yt-yh);
    gv1_108sub->SetPointError(ix, 0, ye);
  }

  UInt_t ic = -1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

  gv1_132sub->SetMarkerStyle(20);
  gv1_132sub->SetMarkerColor(2);
  gv1_132sub->SetLineColor(2);

  gv1_108sub->SetMarkerStyle(20);
  gv1_108sub->SetMarkerColor(4);
  gv1_108sub->SetLineColor(4);

  //auto mgr = new TMultiGraph("mgr",";y_{cm}/y_{cm};v1(t-^3He)");
  auto mgr = new TMultiGraph("mgr",";y_{cm}/y_{cm};v2(t-^3He)");

  mgr->Add(gv1_132sub);
  mgr->Add(gv1_108sub);

  mgr->Draw("ALP");

}
