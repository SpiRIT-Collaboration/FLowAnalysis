void SetStyle()
{

  gStyle -> SetOptStat(0);

  //gStyle -> SetTitleFont(1);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleXSize(0.04);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYSize(0.04);
  gStyle->SetTitleYOffset(1.1);
  gStyle->SetLabelSize(0.04);

  gStyle -> SetPadRightMargin(0.08);
  gStyle -> SetPadTopMargin(0.08);
  gStyle -> SetPadBottomMargin(0.15);
  gStyle -> SetPadLeftMargin(0.15);
  gStyle -> SetTitleOffset(1.2, "x");
  gStyle -> SetTitleOffset(1.2, "y");
  gStyle -> SetTitleSize(0.06, "x");
  gStyle -> SetTitleSize(0.06, "y");
  gStyle -> SetTitleSize(0.06, "z");
  gStyle -> SetLabelSize(0.06, "x");
  gStyle -> SetLabelSize(0.06, "y");
  gStyle -> SetLabelSize(0.06, "z");
  gStyle -> SetTextFont(42);
  gStyle -> SetLegendTextSize(0.06);
  gStyle -> SetLegendBorderSize(0);

}
