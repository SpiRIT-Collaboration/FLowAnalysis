{
  auto gPt_v11=(TGraphErrors*)_file0->Get("gPt_v11");
  auto c1 = new TCanvas("c1","c1",400,450);
  //  c1->SetWindowSize(400,450);
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.05);


  gPt_v11->GetXaxis()->SetLabelSize(0.04);
  gPt_v11->GetYaxis()->SetLabelSize(0.04);
  gPt_v11->GetXaxis()->SetTitleSize(0.05);
  gPt_v11->GetYaxis()->SetTitleSize(0.05);
  gPt_v11->GetYaxis()->SetTitleOffset(1.2);

  gPt_v11->GetYaxis()->SetRangeUser(-0.26,0.06);
  gPt_v11->Draw("ALP");
}
