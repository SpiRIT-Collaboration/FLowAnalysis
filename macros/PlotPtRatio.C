void PlotPtRatio()
{
  Double_t x0[7], y0[7], ex0[7], ey0[7];
  x0[0]=65.909;  y0[0]=-0.00141385; ex0[0]=23.6858; ey0[0]=5.26389e-05;
  x0[1]=153.372; y0[1]=0.000662766; ex0[1]=28.5175; ey0[1]=0.00265336;
  x0[2]=249.824; y0[2]=-0.00545246; ex0[2]=28.6952; ey0[2]=5.6e-05;
  x0[3]=347.548; y0[3]=-0.0229501;  ex0[3]=28.6582; ey0[3]=0.000510946;
  x0[4]=445.792; y0[4]=-0.0348913;  ex0[4]=28.584;  ey0[4]=0.00281607;
  x0[5]=544.918; y0[5]=-0.0478054;  ex0[5]=28.5746; ey0[5]=0.00292702;
  x0[6]=645.34;  y0[6]=-0.0610387;  ex0[6]=28.7195; ey0[6]=0.00251856;

  Double_t x1[7], y1[7], ex1[7], ey1[7];
  x1[0]=65.6738; y1[0]=-0.0027583;  ex1[0]=23.8151; ey1[0]=0.00125086;
  x1[1]=153.089; y1[1]=-0.00414966; ex1[1]=28.5323; ey1[1]=0.00324138;
  x1[2]=249.5;   y1[2]=-0.00850036; ex1[2]=28.713;  ey1[2]=0.000430481;
  x1[3]=347.372; y1[3]=-0.0181119;  ex1[3]=28.7082; ey1[3]=0.000792238;
  x1[4]=445.625; y1[4]=-0.0316037;  ex1[4]=28.5943; ey1[4]=0.00486763;
  x1[5]=544.644; y1[5]=-0.0489098;  ex1[5]=28.5095; ey1[5]=0.00287252;
  x1[6]=644.999; y1[6]=-0.0543365;  ex1[6]=28.7195; ey1[6]=0.00433803;


  auto gPt_v22ratio = new TGraphErrors();
  gPt_v22ratio->SetName("gPt_v22ratio");
  gPt_v22ratio->SetTitle(";Pt; v2(108Sn/132Sn)");

  for( UInt_t i = 0; i < 7; i++ ) {

    auto yr = y1[i]/ y0[i];
      
    gPt_v22ratio->SetPoint(i, x0[i], yr);
  
  }

  auto cc0 = new TCanvas("cc0","cc0");
  gPt_v22ratio->Draw("ALP");
}

//  LocalWords:  x0 y0 ey0 26389e 6e
