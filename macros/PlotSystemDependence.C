Color_t icol[] = {  kRed, kBlue, kSpring, kMagenta, kOrange, kViolet};

void PlotSystemDependence()
{
  Int_t    ax1[12];
  Double_t ay1[12];
  Double_t aex1[12];
  Double_t aey1[12];
  Int_t    ax2[12];
  Double_t ay2[12];
  Double_t aex2[12];
  Double_t aey2[12];

  //void data_v39.1.0()
  ax1[0]=0;   ay1[0]=0.339771;   aex1[0]=0;  aey1[0]=0.00314579;
  ax1[1]=0;   ay1[1]=0.498374;   aex1[1]=0;  aey1[1]=0.00384451;
  ax1[2]=0;   ay1[2]=0.718733;   aex1[2]=0;  aey1[2]=0.00499022;
  ax1[3]=0;   ay1[3]=0.76137;    aex1[3]=0;  aey1[3]=0.00749204;

  ax1[4]=2;   ay1[4]=0.332628;   aex1[4]=0;  aey1[4]=0.00326629;
  ax1[5]=2;   ay1[5]=0.468602;   aex1[5]=0;  aey1[5]=0.00431556;
  ax1[6]=2;   ay1[6]=0.664346;   aex1[6]=0;  aey1[6]=0.00633387;
  ax1[7]=2;   ay1[7]=0.734197;   aex1[7]=0;  aey1[7]=0.0082982;

  ax1[8]=1;   ay1[8]=0.330547;   aex1[8]=0;  aey1[8]=0.00381516;
  ax1[9]=1;   ay1[9]=0.465106;   aex1[9]=0;  aey1[9]=0.00468351;
  ax1[10]=1;  ay1[10]=0.671668;  aex1[10]=0; aey1[10]=0.00639089;
  ax1[11]=1;  ay1[11]=0.730247;  aex1[11]=0; aey1[11]=0.00896295;

  ax2[0]=0;   ay2[0]=-0.0615601; aex2[0]=0;  aey2[0]=0.00870416;
  ax2[1]=0;   ay2[1]=-0.0942898; aex2[1]=0;  aey2[1]=0.00893068;
  ax2[2]=0;   ay2[2]=-0.126753;  aex2[2]=0;  aey2[2]=0.00967359;
  ax2[3]=0;   ay2[3]=-0.123328;  aex2[3]=0;  aey2[3]=0.00998096;

  ax2[4]=2;   ay2[4]=-0.0477709; aex2[4]=0;  aey2[4]=0.00998066;
  ax2[5]=2;   ay2[5]=-0.079231;  aex2[5]=0;  aey2[5]=0.0101252;
  ax2[6]=2;   ay2[6]=-0.112404;  aex2[6]=0;  aey2[6]=0.0117281;
  ax2[7]=2;   ay2[7]=-0.134334;  aex2[7]=0;  aey2[7]=0.0122794;

  ax2[8]=1;   ay2[8]=-0.0597437; aex2[8]=0;  aey2[8]=0.0100378;
  ax2[9]=1;   ay2[9]=-0.0870866; aex2[9]=0;  aey2[9]=0.0100385;
  ax2[10]=1;  ay2[10]=-0.102265; aex2[10]=0; aey2[10]=0.0107724;
  ax2[11]=1;  ay2[11]=-0.113521; aex2[11]=0; aey2[11]=0.0113059;
  //--------------


  TString  sysname[] = {"132SN","108Sn","112Sn"};
  Double_t sysdlt[]  = {0.22,    0.09,    0.15};
  Double_t sysA[]    = {256.,    220.,    236};
  TString  lpart[]   = {"1H",   "2H",    "3H",  "3He"};

  auto mv1 = new TMultiGraph();
  auto mv2 = new TMultiGraph();
  auto lv1 = new TLegend(0.6, 0.2, 0.7, 0.5);
  auto lv2 = new TLegend(0.2, 0.2, 0.3, 0.5);


  TGraphErrors *grv1[4];
  TGraphErrors *grv2[4];

  for(UInt_t is = 0; is < 4; is++){

    grv1[is] = new TGraphErrors();
    grv2[is] = new TGraphErrors();

    for(UInt_t ip = 0; ip < 3; ip++) {

      UInt_t ik = 4*ip+is;
      grv1[is]->SetPoint(ip, sysdlt[ax1[ik]], ay1[ik]);
      grv1[is]->SetPointError(ip, 0., aey1[ik]);

      grv2[is]->SetPoint(ip, sysdlt[ax2[ik]], -ay2[ik]);
      grv2[is]->SetPointError(ip, 0., aey2[ik]);

     
    }      

    grv1[is]->SetMarkerStyle(20);
    grv2[is]->SetMarkerStyle(20);
    grv1[is]->SetMarkerColor(icol[is]);
    grv2[is]->SetMarkerColor(icol[is]);
    grv1[is]->SetLineColor(icol[is]);
    grv2[is]->SetLineColor(icol[is]);

    mv1->Add(grv1[is]);
    mv2->Add(grv2[is]);
      
    lv1->AddEntry(grv1[is], lpart[is], "LP");
    lv2->AddEntry(grv2[is], lpart[is], "LP");

  }


  TCanvas *cc;

  UInt_t ic = 0;
  cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;
  mv1->SetTitle("; (n-p)/A; v1 slope");
  mv1->Draw("AP");
  lv1->Draw();

  cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;
  mv2->SetTitle("; (n-p)/A; -v2 max");
  mv2->Draw("AP");
  lv2->Draw();
}

void PlotSystemDependence_part()
{
  Double_t ax1[12];
  Double_t ay1[12];
  Double_t aex1[12];
  Double_t aey1[12];
  Double_t ax2[12];
  Double_t ay2[12];
  Double_t aex2[12];
  Double_t aey2[12];

  //void data_v39.1.0()
  ax1[0]=0;   ay1[0]=0.339771;   aex1[0]=0;  aey1[0]=0.00314579;
  ax1[1]=0;   ay1[1]=0.498374;   aex1[1]=0;  aey1[1]=0.00384451;
  ax1[2]=0;   ay1[2]=0.718733;   aex1[2]=0;  aey1[2]=0.00499022;
  ax1[3]=0;   ay1[3]=0.76137;    aex1[3]=0;  aey1[3]=0.00749204;

  ax1[4]=3;   ay1[4]=0.332628;   aex1[4]=0;  aey1[4]=0.00326629;
  ax1[5]=3;   ay1[5]=0.468602;   aex1[5]=0;  aey1[5]=0.00431556;
  ax1[6]=3;   ay1[6]=0.664346;   aex1[6]=0;  aey1[6]=0.00633387;
  ax1[7]=3;   ay1[7]=0.734197;   aex1[7]=0;  aey1[7]=0.0082982;

  ax1[8]=1;   ay1[8]=0.330547;   aex1[8]=0;  aey1[8]=0.00381516;
  ax1[9]=1;   ay1[9]=0.465106;   aex1[9]=0;  aey1[9]=0.00468351;
  ax1[10]=1;  ay1[10]=0.671668;  aex1[10]=0; aey1[10]=0.00639089;
  ax1[11]=1;  ay1[11]=0.730247;  aex1[11]=0; aey1[11]=0.00896295;

  ax2[0]=0;   ay2[0]=-0.0615601; aex2[0]=0;  aey2[0]=0.00870416;
  ax2[1]=0;   ay2[1]=-0.0942898; aex2[1]=0;  aey2[1]=0.00893068;
  ax2[2]=0;   ay2[2]=-0.126753;  aex2[2]=0;  aey2[2]=0.00967359;
  ax2[3]=0;   ay2[3]=-0.123328;  aex2[3]=0;  aey2[3]=0.00998096;

  ax2[4]=3;   ay2[4]=-0.0477709; aex2[4]=0;  aey2[4]=0.00998066;
  ax2[5]=3;   ay2[5]=-0.079231;  aex2[5]=0;  aey2[5]=0.0101252;
  ax2[6]=3;   ay2[6]=-0.112404;  aex2[6]=0;  aey2[6]=0.0117281;
  ax2[7]=3;   ay2[7]=-0.134334;  aex2[7]=0;  aey2[7]=0.0122794;

  ax2[8]=1;   ay2[8]=-0.0597437; aex2[8]=0;  aey2[8]=0.0100378;
  ax2[9]=1;   ay2[9]=-0.0870866; aex2[9]=0;  aey2[9]=0.0100385;
  ax2[10]=1;  ay2[10]=-0.102265; aex2[10]=0; aey2[10]=0.0107724;
  ax2[11]=1;  ay2[11]=-0.113521; aex2[11]=0; aey2[11]=0.0113059;
  //--------------


  TString  sysname[] = {"132SN","108Sn","112Sn"};
  Double_t sysdlt[]  = {0.22,    0.09,    0.15};
  Double_t sysA[]    = {256.,    220.,    236};
  TString  lpart[]   = {"1H",   "2H",    "3He"};

  auto mv1 = new TMultiGraph();
  auto mv2 = new TMultiGraph();
  auto lv1 = new TLegend(0.6, 0.2, 0.8, 0.3);
  auto lv2 = new TLegend(0.2, 0.2, 0.5, 0.3);


  TGraphErrors *grv1[3];
  TGraphErrors *grv2[3];

  UInt_t ik = 0;
  for(UInt_t is = 0; is < 3; is++){

    grv1[is] = new TGraphErrors();
    grv2[is] = new TGraphErrors();

    for(UInt_t ip = 0; ip < 4; ip++) {

      grv1[is]->SetPoint(ip, ip, ay1[ik]);
      grv1[is]->SetPointError(ip, 0., aey1[ik]);

      grv2[is]->SetPoint(ip, ip, ay2[ik]);
      grv2[is]->SetPointError(ip, 0., aey2[ik]);

      ik++;
    }      

    grv1[is]->SetMarkerStyle(20);
    grv2[is]->SetMarkerStyle(20);
    grv1[is]->SetMarkerColor(is+2);
    grv2[is]->SetMarkerColor(is+2);
    grv1[is]->SetLineColor(is+2);
    grv2[is]->SetLineColor(is+2);

    mv1->Add(grv1[is]);
    mv2->Add(grv2[is]);
      
    lv1->AddEntry(grv1[is], sysname[is], "lp");
    lv2->AddEntry(grv2[is], sysname[is], "lp");

  }


  TCanvas *cc;

  UInt_t ic = 0;
  cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;
  mv1->SetTitle("; A(Beam+Target); v1");
  mv1->Draw("AP");
  lv1->Draw();

  cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;
  mv2->SetTitle("; A(Beam+Target); v2");
  mv2->Draw("AP");
  lv2->Draw();
}


