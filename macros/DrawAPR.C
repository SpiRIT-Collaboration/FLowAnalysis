void DrawAPR()
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(22);

  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleXSize(0.04);
  gStyle->SetTitleXOffset(1.2);

  gStyle->SetTitleYSize(0.04);
  gStyle->SetTitleYOffset(1.1);

  gStyle->SetLabelSize(0.04);
  

  auto afile = new TFile("data/APR2018.root");
  
  TF1  *ffv1   = (TF1*)afile->Get("fv1");
  TF1  *ffv2   = (TF1*)afile->Get("fv2");
  TF1  *ffv1p  = (TF1*)afile->Get("fv1p");

  TH1D *hphi0 = (TH1D*)afile->Get("hphi0");
  hphi0->SetTitle("; #Delta(#phi - #Psi); 2#pi/N dN/d(#phi - #Psi)");

  TH1D *hphi1 = (TH1D*)afile->Get("hphi1");
  hphi1->SetTitle("; #Delta(#phi - #Psi); 2#pi/N dN/d(#phi - #Psi)");

  TH1D *hphi2 = (TH1D*)afile->Get("hphi2");
  hphi2->SetTitle("; #Delta(#phi - #Psi); 2#pi/N dN/d(#phi - #Psi)");

  auto lb0 = new TLatex(1.6, 1.17,"(a)");
  lb0->SetTextSize(0.08);
  auto lb1 = new TLatex(2.1, 1.031,"(b)");
  lb1->SetTextSize(0.08);
  auto lb2 = new TLatex(2.08, 1.17,"(c)");
  lb2->SetTextSize(0.08);

  hphi0->GetXaxis()->SetLabelSize(0.05);
  hphi0->GetYaxis()->SetLabelSize(0.05);
  hphi0->GetXaxis()->SetTitleOffset(1.);
  hphi0->GetXaxis()->SetTitleSize(0.05);
  hphi0->GetYaxis()->SetTitleOffset(1.5);
  hphi0->GetYaxis()->SetTitleSize(0.05);

  hphi1->GetXaxis()->SetLabelSize(0.05);
  hphi1->GetYaxis()->SetLabelSize(0.05);
  hphi1->GetXaxis()->SetTitleOffset(1.);
  hphi1->GetXaxis()->SetTitleSize(0.05);
  hphi1->GetYaxis()->SetTitleOffset(1.5);
  hphi1->GetYaxis()->SetTitleSize(0.05);

  hphi2->GetXaxis()->SetLabelSize(0.05);
  hphi2->GetYaxis()->SetLabelSize(0.05);
  hphi2->GetXaxis()->SetTitleOffset(1.);
  hphi2->GetXaxis()->SetTitleSize(0.05);
  hphi2->GetYaxis()->SetTitleOffset(1.5);
  hphi2->GetYaxis()->SetTitleSize(0.05);


  auto c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(3,1);

  c1->SetBorderSize(1);
  c1->GetPad(1)->SetRightMargin(0.01);
  c1->GetPad(1)->SetLeftMargin(0.16);
  c1->GetPad(2)->SetRightMargin(0.01);
  c1->GetPad(2)->SetLeftMargin(0.16);
  c1->GetPad(3)->SetRightMargin(0.01);
  c1->GetPad(3)->SetLeftMargin(0.16);

  c1->cd(1);
  hphi0->Draw("e");
  ffv1->Draw("same");
  lb0->Draw();


  c1->cd(2);
  hphi1->Draw("e");
  ffv2->Draw("same");
  lb1->Draw();

  c1->cd(3);
  hphi2->Draw("e");
  ffv1p->Draw("same");
  lb2->Draw();


}
