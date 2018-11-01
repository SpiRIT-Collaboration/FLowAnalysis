#include <algorithm>  
#include "SetStyle.C"
TCanvas *cc[10];
auto fv1v2 = new TF1("fv1v2","1+2.*[0]*cos(x-[2]) + 2.*[1]*cos(2.* (x-[2]))",-3.14, 3.14); 
auto fgaus = new TF1("fgaus","gaus",0.,100.);
std::vector< TVector2 > rndV;
TVector2 vcsum;

UInt_t itry = 100;

STBootStrap *bstrap  = new STBootStrap(1);
STBootStrap *bstrapD = new STBootStrap(1);

TH1D *hphidist;
TH2D *hpxy;
TGraphErrors *gstrap;
TGraphErrors *gsdiv;
TGraphErrors *gsdiff;
TGraphErrors *gscomp;

TH1D *hphi;
TH1D *hphio;
TH1D *hstd;

TH1D *hmdphi[20];
TH1D *hmdphio[20];

TH1D *hbsPhiDist;

Bool_t first = kTRUE;

void  SingleBootStrap(UInt_t np, Double_t val);
void  RndBootStrap(UInt_t np, UInt_t ith);
void  RPResolution(UInt_t np, UInt_t ntry);
void  DoubleBootStrap();

void BootStrapTest()
{
  fv1v2->SetParameter(0, 0.);
  fv1v2->SetParameter(1, -0.05);

  //  RPResolution();


  gROOT->ProcessLine(".! grep -i void BootStrapTest.C ");

  //  fv1v2->SetLineColor(2);
  //  fv1v2->Draw("LP");

  DoubleBootStrap();
}


void  Rndv1v2Distribution(UInt_t np = 50)
{
  TH1F *hphidist;
  //  if( hphidist == NULL ) {
    TFile *fin = TFile::Open("FlwRUN2841m113_rf.v7.0.Phi.root");
    //fin->ls();
    if( !fin->IsOpen() ) return;
    hphidist = (TH1F*)fin->Get("hphi");
    Double_t entry = hphidist->GetEntries();
    Double_t nbin  = (Double_t)hphidist->GetNbinsX()/(hphidist->GetXaxis()->GetXmax()-hphidist->GetXaxis()->GetXmin());
    //  }


  // auto hexpphi  = new TH1D("hexpphi","",100,-3.2, 3.2);
  // auto hexpphic = new TH1D("hexpphic","",100,-3.2, 3.2);

  //fin->Close();
  if( hphidist == NULL ) return;

  TDatime dtime;
  TRandom3 grnd(dtime.GetSecond());
  gRandom->SetSeed(dtime.GetSecond());
  TRandom3 ggr(dtime.GetMinute());
  
  rndV.clear();
  vcsum = TVector2(0, 0);

  // TGraph *gr  = new TGraph();
  // TGraph *grd = new TGraph();

  UInt_t i = 0; UInt_t k = 0;
  while( i < np ) {
    auto rnd = ggr.Rndm();
    auto phi = fv1v2->GetRandom();

    auto xbin = hphidist->FindBin(phi);
    auto rate = hphidist->GetBinContent(xbin)/entry*nbin;

    //    hexpphic->Fill(phi);
    //    cout << " phi " << phi << " " << xbin << " " << rate << " : " << rnd << endl;
    
    // gr->SetPoint(k, phi, rate); 
    // grd->SetPoint(k, phi, rnd);
    // k++;
    
    if( rnd > rate ) continue;
    i++;
    
    //    cout << " --> " << phi << " " << xbin << " " << rate << endl;

    //    hexpphi->Fill(phi);

    rndV.push_back( TVector2( cos(phi), sin(phi) ) );

    vcsum += TVector2( cos(phi), sin(phi) );
  }

  // cc[0] = new TCanvas("cc0","cc0");
  // hexpphi->Draw();
  // //  hphidist->Draw("same");

  // cc[1] = new TCanvas("cc1","cc1");
  // hexpphic->Draw();

  // cc[2] = new TCanvas("cc2","cc2");
  // gr->SetMarkerStyle(20);
  // gr->Draw("AP");

  // cc[3] = new TCanvas("cc3","cc3");
  // grd->SetMarkerStyle(20);
  //  grd->Draw("AP");
}

void RndBootStrap(UInt_t np = 50, UInt_t ith = 0)
{

  Rndv1v2Distribution(np);
  
  if( ith == 0 ) {
    if( hphidist == NULL)
      hphidist = new TH1D("hphidist","", 100, -3.14, 3.14);
    else
      hphidist->Reset();

    if( hpxy == NULL)
      hpxy     = new TH2D("hpxy","",100,-1.2, 1.2, 100, -1.2, 1.2);
    else
      hpxy->Reset();
  }


  bstrap->clear();


  for(UInt_t i = 0; i < np; i++) {
      
    auto phi = rndV.at(i);

    TVector2 vec = rndV.at(i); 
      
    if( ith == 0 )
      hphidist->Fill( TVector2::Phi_mpi_pi(vec.Phi()) );


    if( ith == 0 )
      hpxy->Fill( vec.X(), vec.Y());

    bstrap->Add(vec);

    bstrapD->Add( TVector2::Phi_mpi_pi( vec.Phi() ) );

  }

  bstrap->BootStrapping(itry);

  bstrapD->BootStrapping(itry);

  hbsPhiDist = bstrap->GetBSPhiDistribution();
}

void DetectorBias(UInt_t np = 50, UInt_t ntry = 1000)
{
  for( UInt_t n = 0; n < ntry; n++ ) {
    

    
  } 
  

}



void RPResolution(UInt_t np = 50, UInt_t ntry = 1000)
{
  TDatime dtime;
  TRandom3 grnd(0);
  TRandom3 grr(0);
  gRandom->SetSeed(dtime.GetSecond());

  auto hPhiAll  = new TH1D("hPhiAll","" ,100, -3.14, 3.14);
  auto hPhiAllo = new TH1D("hPhiAllo","",100, -3.14, 3.14);
  auto hPhiGen  = new TH1D("hPhiGen","" ,100, -3.14, 3.14);
  
  auto hAoGen   = new TH2D("hAoGen","Ordinary vs Generated", 100, -3.14, 3.14, 100, -3.14, 3.14);
  auto hAlGen   = new TH2D("hAlGen","BootStrap vs Generated", 100, -3.14, 3.14, 100, -3.14, 3.14);

  auto hdAoGen  = new TH2D("hdAoGen","Diff Ordinary vs Generated; R.P. Angle[Rad]",  100, -3.14, 3.14, 100, -2.5, 2.5);
  auto hdAlGen  = new TH2D("hdAlGen","Diff BootStrap vs Generated;R.P. Angle[Rad]", 100, -3.14, 3.14, 100, -2.5, 2.5);
  auto hdAlAo   = new TH2D("hdAlAo","Diff BootStrap vs Ordinary;  R.P. Angle[Rad]", 100, -3.14, 3.14, 100, -2.5, 2.5);


  std::cout << " Number of particles " << np
	    << " ; Number of BootStrapping " << itry
	    << " ; Statistics " << ntry
	    << std::endl;

  for( UInt_t n = 0; n < ntry; n++ ){

    if( n%(ntry/10) == 0 ) std::cout << "Processing.. "<< n << " " << Double_t(n)/Double_t(ntry)*100. << "%" << std::endl;

    Double_t initPhi = (2.*grr.Rndm() - 1) * TMath::Pi();
    fv1v2->SetParameter(2, initPhi);
    hPhiGen->Fill(initPhi);

    RndBootStrap(np, n+1);
    hPhiAll->Fill(bstrap->GetMean());
		  
    hPhiAllo->Fill(bstrap->GetOrdinaryMean());

    hAoGen->Fill(initPhi, bstrap->GetOrdinaryMean());
    hAlGen->Fill(initPhi, bstrap->GetMean());

    hdAoGen->Fill(initPhi, TVector2::Phi_mpi_pi(bstrap->GetOrdinaryMean()-initPhi));
    hdAlGen->Fill(initPhi, TVector2::Phi_mpi_pi(bstrap->GetMean()-initPhi));

    hdAlAo->Fill(initPhi, TVector2::Phi_mpi_pi(bstrap->GetMean()-bstrap->GetOrdinaryMean()));

  }    

  cc[0] = new TCanvas("cc0","cc0");
  hPhiAll->Draw("");

  cc[1] = new TCanvas("cc1","cc1");
  hAoGen->Draw("colz");

  cc[2] = new TCanvas("cc2","cc2");
  hAlGen->Draw("colz");

  cc[3] = new TCanvas("cc3","cc3");
  hdAoGen->Draw("colz");

  cc[4] = new TCanvas("cc4","cc4");
  hdAlGen->Draw("colz");

  cc[5] = new TCanvas("cc5","cc5");
  hdAlAo->Draw("colz");

  cc[6] = new TCanvas("cc6","cc6");
  hdAoGen->ProjectionY("hdAoGen_py");
  TH1D *hdAoGen_py = (TH1D*)gROOT->Get("hdAoGen_py");
  hdAoGen_py->Draw();

  cc[7] = new TCanvas("cc7","cc7");
  hdAlGen->ProjectionY("hdAlGen_py");
  TH1D *hdAlGen_py = (TH1D*)gROOT->Get("hdAlGen_py");
  hdAlGen_py->Draw();

}


void MultiplicityDependenceOrdinary(Double_t val = 0., UInt_t mtry = 1000)
{
  fv1v2->SetParameter(2, val);

  auto gmdphi = new TGraphErrors();
  auto gmdphio= new TGraphErrors();
  
  auto hmoddiff = new TH2D("hmoddiff",";sum.Mod; Diff",100, 0., 0.8, 100, -3.2, 3.2);
 
  const UInt_t ntry = 1;
  for(UInt_t n = 0; n < ntry; n++){
    UInt_t mult = (n+1)*30;
    
    hmdphi[n]  = new TH1D(Form("hmdphi%d",n) ,Form("mult %d",mult), 100, -3.2, 3.2);
    
    std::vector< Double_t > buf_diff;
    for(UInt_t k = 0; k < mtry; k++){
      Rndv1v2Distribution(mult);
      
      hmdphi[n]->Fill(TVector2::Phi_mpi_pi( vcsum.Phi() - val) );

      hmoddiff->Fill( vcsum.Mod()/(Double_t)mult, TVector2::Phi_mpi_pi( vcsum.Phi() - val) );
    }

  }

  

  
  cc[0] = new TCanvas("cc0","cc0");
  hmdphi[0]->Draw();

  cc[1] = new TCanvas("cc1","cc1");
  hmoddiff->Draw("colz");

}

void MultiplicityDependence(Double_t val = 0., UInt_t mtry = 100)
{
  fv1v2->SetParameter(2, val);

  const UInt_t ntry = 6;
  auto gmdphi = new TGraphErrors();
  auto gmdphio= new TGraphErrors();
  
  auto hmoddiff = new TH2D("hmoddiff",";sum.Mod; Diff",100, 0., 0.8, 100, -3.2, 3.2);

  TH1D *herror[ntry];
  TGraph *hmoder[ntry];

  for(UInt_t n = 0; n < ntry; n++){
    UInt_t mult = (n+1)*5;
    
    hmdphi[n]  = new TH1D(Form("hmdphi%d",n) ,Form("mult %d",mult), 100, -3.2, 3.2);
    herror[n]  = new TH1D(Form("herror%d",n) ,Form("mult %d",mult), 100, 0., 3.2);
    hmoder[n]  = new TGraph();
    hmoder[n]->SetName(Form("hmoder%d",n));
    hmoder[n]->SetTitle(Form("mult %d",mult));
    
    std::vector< Double_t > buf_diff;
    for(UInt_t k = 0; k < mtry; k++){
      RndBootStrap(mult, 1);
      
      buf_diff.push_back( TVector2::Phi_mpi_pi(bstrap->GetMean() - val) );
      hmdphi[n]->Fill(TVector2::Phi_mpi_pi( bstrap->GetMean() - val) );

      hmoddiff->Fill(bstrap->GetOrdinarySum().Mod()/(Double_t)mult, TVector2::Phi_mpi_pi(bstrap->GetMean() - val) );


      herror[n]->Fill(bstrap->GetError());
      hmoder[n]->SetPoint(k, bstrap->GetOrdinarySum().Mod()/(Double_t)mult, bstrap->GetError());
    }

    Double_t diff = TMath::Mean( buf_diff.begin(), buf_diff.end() );
    Double_t diffs= TMath::StdDev( buf_diff.begin(), buf_diff.end() );

    cout << mtry << " n " << setw(3) <<  n << " : " << setw(15) << diff << " +- " << setw(15) <<  diffs << endl;

    gmdphi->SetPoint(n, (Double_t)mult, diff);
    gmdphi->SetPointError(n, 0., diffs);

  }
  
  cc[0] = new TCanvas("cc0","cc0");
  gmdphi->SetTitle((TString)Form("%d try",mtry)+"; Multiplicity; Deviation");
  gmdphi->SetMarkerStyle(21);
  gmdphi->SetMarkerColor(4);
  gmdphi->Draw("AP");
  hmdphi[0]->Draw();


  cc[1] = new TCanvas("cc1","cc1",1200,1000);
  cc[1]->Divide(2,3);
  for(UInt_t i = 0; i < ntry; i++) {
    cc[1]->cd(i+1);
    hmdphi[i]->Draw();
  }

  cc[2] = new TCanvas("cc2","cc2");
  hmoddiff->Draw("colz");

  cc[3] = new TCanvas("cc3","cc3",1200,1000);
  cc[3]->Divide(2,3);
  for(UInt_t i = 0; i < ntry; i++) {
    cc[3]->cd(i+1);
    herror[i]->Draw();
  }

  cc[4] = new TCanvas("cc4","cc4",1200,1000);
  cc[4]->Divide(2,3);
  for(UInt_t i = 0; i < ntry; i++) {
    cc[4]->cd(i+1);
    hmoder[i]->SetMarkerStyle(20);
    hmoder[i]->Draw("AP");
  }

}


void SingleBootStrap(UInt_t np = 50, Double_t val = 0.)
{
  Double_t gphi = val;
  fv1v2->SetParameter(2, gphi);
  //  cc[9] = new TCanvas("cc9","cc9");

  RndBootStrap(np, 0);

  Double_t sphi = bstrap->GetOrdinaryMean();
  Double_t phie = bstrap->GetOrdinaryStdDev();


  gstrap = new TGraphErrors();
  gstrap->SetTitle(";Number of BootStrapping; Reaction Plane Orientaion [rad]");

  gsdiv  = new TGraphErrors();
  gsdiff = new TGraphErrors();


  if (hphi == NULL)
    hphi  = new TH1D("hphi" ,"BootStrap method",100, sphi-2.5, sphi+2.5);
  else
    hphi->Reset();


  for(UInt_t it = 1; it < itry; it++){

    if( it%100 == 0 ) std::cout << "Processing.. "<< it << " " << Double_t(it)/Double_t(itry)*100. << std::endl;

    gstrap->SetPoint(it-1, (Double_t)it, bstrap->GetResidualMean(it));
    gstrap->SetPointError(it-1, 0., bstrap->GetResidualStdDev(it));

    hphi->Fill( bstrap->GetResidualMean(it) );

    auto diff = bstrap->GetResidualStdDev(it);
    gsdiv->SetPoint(it-1, (Double_t)it, diff);

    gsdiff->SetPoint(it-1, (Double_t)it, bstrap->GetResidualMean(it)-sphi);
  }



  cout << " BST mean " << bstrap->GetMean() << " +- " << bstrap->GetError() << endl;
  cout << " vsum Mag " << (bstrap->GetOrdinarySum()).Mod() << endl;
  cout << " vsum Phi " << sphi << " +- " << bstrap->GetOrdinaryStdDev() << endl;
  cout << " deviation " << abs(bstrap->GetMean()-val)/bstrap->GetStdDev() << endl;

  cout << " ---- > Double < ----- " << endl;
  cout << " mean " << bstrapD->GetMean() << " +- " << bstrap->GetError() << endl;


  Double_t yp = hphi->GetMaximum();
  // gscomp->SetPoint(0, val, yp);
  // gscomp->SetPointError(0, 0., 0.);

  auto arrOrdn = new TArrow(sphi-bstrap->GetOrdinaryStdDev(), yp/2., sphi+bstrap->GetOrdinaryStdDev(), yp/2., 0.01, "<>");
  arrOrdn->SetLineColor(2);
  auto arrBS   = new TArrow(bstrap->GetMean()-bstrap->GetError(), yp, bstrap->GetMean()+bstrap->GetError(), yp, 0.01, "<>");
  arrBS->SetLineColor(4);



  // plotting histgrams
  if( cc[0] == NULL ) {
    cc[0] = new TCanvas("cc0","cc0",1200,1000);
  }

  cc[0]->Divide(2,3);

  UInt_t id = 1;
  cc[0]->cd(id); id++;
  fv1v2->Draw("LP");

  //  cc[5] = new TCanvas("cc5","cc5");
  cc[0]->cd(id); id++;
  hphidist->Draw();


  cc[0]->cd(id); id++;
  gstrap->SetMarkerStyle(20);
  gstrap->SetMarkerColor(4);
  gstrap->Draw("ap");

  auto aLine = new TLine(0., sphi, gstrap->GetXaxis()->GetXmax(), sphi);
  aLine->SetLineColor(2);
  aLine->Draw();


  //  cc[3] = new TCanvas("cc3","cc3");
  cc[0]->cd(id); id++;
  hpxy->SetMarkerStyle(20);
  hpxy->SetMarkerSize(0.5);
  hpxy->SetMarkerColor(2);
  hpxy->Draw();

  Double_t x = (bstrap->GetOrdinarySum()).X()/(Double_t)np;
  Double_t y = (bstrap->GetOrdinarySum()).Y()/(Double_t)np;
  auto vline = new TArrow(0.,0., x, y, 0.01, ">");
  vline->SetLineColor(4);
  vline->Draw();
  
  x = cos( bstrap->GetMean() );
  y = sin( bstrap->GetMean() );
  auto bline = new TArrow(0.,0., x, y, 0.01, ">");
  bline->SetLineColor(2);
  bline->Draw();
  auto kline = new TLine(0.,0., cos(gphi), sin(gphi));
  kline->SetLineColor(8);
  kline->Draw();


  //  cc[1] = new TCanvas("cc1","cc1");
  cc[0]->cd(id); id++;
  hphi->SetFillColor(4);
  hphi->SetLineColor(4);
  hphi->Draw();
  auto oline = new TLine(sphi, 0., sphi, hphi->GetMaximum());
  oline->SetLineColor(2);
  oline->Draw();
  auto gline = new TLine(gphi, 0., gphi, hphi->GetMaximum());
  gline->SetLineColor(8);
  gline->Draw();
  
  arrOrdn->Draw();
  arrBS->Draw();

}

void SaveCanvas(TString str = "")
{

  Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
  for(Int_t i = 0; i < iCanvas; i++) {
    TString scname = "bs" + str + Form("%d.png",i);
    gROOT->GetListOfCanvases()->At(i)->SaveAs(scname);
  }
}




void DoubleBootStrap()
{
  SetStyle();

  std::vector< Double_t > data0{
    3.11, 8.88, 9.26, 18.36, 18.43, 19.27, 24.58, 25.13, 26.24, 36.37, 38.64, 39.16, 50.39, 52.75, 54.80,
    10.81, 12.69, 13.78, 19.50, 19.54, 20.16, 26.26, 27.65, 28.06, 41.02, 42.97, 44.08, 59.07, 61.22, 70.32,
    15.23, 20.59, 28.08, 44.67, 82.70,
    15.62, 17.00, 22.22, 23.04, 28.38, 32.03, 45.40, 46.69, 85.76, 86.37,
    17.39, 24.47, 34.98, 48.65, 93.34 };

  std::vector< Double_t > data1{3.12, 0.00, 1.57, 19.67, 0.22, 2.20};

  std::vector< Double_t > data2{
    142, 232, 175, 50, 197.5, 146.5, 149.4, 155, 705, 1850,
    132.5, 200, 362, 215, 260, 307, 116.7, 449.9, 266, 244.9, 66.407, 166, 290, 164.95, 375,
    244.95, 335, 210.95, 1370, 265, 256, 296, 148.5, 335, 987.5,
    324.5, 222, 225, 215.5, 179.8, 217, 684.5, 257, 570, 270, 252.95, 507, 330, 149.95, 190 };

  std::vector< Double_t > data3{ 43, 45, 74, 79, 84, 88, 100, 100, 109, 113, 139, 144, 191, 198, 522, 598, 
      53, 56, 80, 80, 89, 91, 101, 102, 114, 118, 145, 147, 211, 214, 56, 57, 81, 81, 91, 92, 
      102, 102, 121, 123, 156, 162, 243, 249, 58, 81, 92, 103, 126, 174, 329, 66, 82, 97, 104, 128, 178, 380, 
      67, 73, 83, 83, 99, 99, 107, 108, 137, 138, 179, 184, 403, 511};

  std::vector< Double_t > data6{3.13797, -1.9299, -1.18035, 2.7284, -1.48609, -0.0786593, 2.8077, 1.33971, 
      0.209839, 1.31363, 1.43447, 0.84654, -0.990638, 1.71389, 0.102992, -1.89179, -0.128132, -0.568639, 
      -1.54816, -1.60144, -2.90241, -0.890058, -1.72326, 2.70165, 0.420397, 2.36085, 0.885057, -0.00727319, 
      0.31821, -1.80072, -1.09973, -2.24691, -2.65091, 0.788468, 1.22745, 0.72817, 1.15998, -2.37638, 
      1.07555, -2.31623, -2.16261, 0.0125769, -1.63575, -1.14374, -2.4979, -2.17052, 0.24836, -2.50429, -1.12503, 2.31213};

  std::vector< Double_t > data7{3.13818, -2.10651, -1.44706, 2.75549, -1.72148, -0.37903, 2.83222, 1.20491, 
      -0.0763654, 1.17399, 1.31745, 0.628329, -1.27235, 1.64885, -0.189652, -2.07395, -0.429916, -0.870439, 
      -1.7762, -1.82293, -2.93077, -1.17826, -1.92895, 2.72929, 0.151065, 2.38053, 0.672539, -0.305086, 
      0.0399924, -1.99582, -1.37325, -2.37492, -2.71541, 0.561993, 1.07201, 0.493534, 0.992499, -2.48385, 
      0.893496, -2.43326, -2.30386, -0.284414, -1.8529, -1.41362, -2.58612, -2.31054, -0.0351751, -2.59151, -1.39648, 2.32838};


  std::vector< Double_t > data8{ 2.99426, 2.67382, -0.135867, 2.86401, -2.33724, 1.19059, 2.32359, 1.95302, -2.30459};//  ; 2.69684

  //5.4028730 * 5.5671996
  std::vector< Double_t > data9{ 0.0564980, 0.2115968, -0.663450, -0.634802, -1.148398, -2.998520, 2.6465250, 
      2.9963681, -2.541403, -2.943897, -2.076007, -1.654803, 1.5584259, -0.985322, -2.362052, -2.920972, -2.693225, 0.5024560};

  //0.3734517 * 0.3366785 
  std::vector< Double_t > data10{0.110158, 5.92771, 1.17038, 0.897871, 4.38574, 4.2925, 3.86337, 3.50714, 6.27399, 1.89132, 0.385921, 5.1921, 5.25916, 3.10276, 5.40767, 0.558457};


  Double_t genmean = 50.;
  fgaus->SetParameter(0,100.);
  fgaus->SetParameter(1,genmean);
  fgaus->SetParameter(2,3.);

  std::vector< Double_t > data4;
  for(UInt_t j = 0; j < 200; j++ )
    data4.push_back( fgaus->GetRandom() );


  std::vector< Double_t > data5{61, 88, 89, 89, 90, 92, 93, 94, 98, 98, 101, 102, 105, 108, 109, 113, 114, 115, 120, 138.};



  std::vector< Double_t > data = data10; 
  genmean  = 0.;

  // ----------------------------------------------------------------------//
  UInt_t dataSize = (UInt_t)data.size();

  //  std::sort( data.begin(), data.end() );
  
  auto horiginal = new TH1D("horiginal","data",100, 0., 6.4);
  auto hboot     = new TH1D("hboot","; #Psi^{A}_{BS}" ,100, 0.,9.);
  auto hqq       = new TGraph();
  auto hbsqq     = new TGraph();
  auto horig     = new TGraph(); horig->SetTitle("; cos(#phi_{k}); sin(#phi_{k})");
  auto hbootstd = new TGraphErrors(); hbootstd->SetTitle(";Number of bootstrapping; #Psi^{A}_{BS}");


  Double_t sum = 0.;

  for( UInt_t i = 0; i < dataSize; i++ ){

    sum += data[i];
    //bstrap->Add( data[i] );
    bstrap->Add( TVector2(cos(data[i]), sin(data[i])) );
  

    horiginal->Fill( data[i] );
    horig->SetPoint(i, cos(data[i]), sin(data[i]) );
  }

  Double_t average = sum/(Double_t)dataSize;

  Double_t dev = 0;
  for( UInt_t i = 0; i < dataSize; i++ ) 
    dev += pow(data[i] - average, 2);

  Double_t dev2 = dev/(Double_t)(dataSize - 1);

  Double_t stddev = TMath::StdDev( data.begin(), data.end() );

  for( UInt_t i = 0; i < (UInt_t)dataSize; i++) {
    Double_t percent = (Double_t)i / (Double_t)dataSize;
    Double_t z_score = ROOT::Math::gaussian_quantile(percent, stddev)/stddev;
    
    if( !std::isnan(z_score) && !std::isinf(z_score) )
      hqq->SetPoint(i, data.at(i), z_score );
  }


  //---------- BootStrap
  
  UInt_t ntry = 10000;
  bstrap->BootStrapping(ntry);

  std::vector< Double_t > replace = bstrap->GetReplaceVector();

  std::sort(replace.begin(), replace.end());

  sum = 0.;
  for( UInt_t i = 0; i < (UInt_t)replace.size(); i++) {
    hboot->Fill(replace.at(i)+bstrap->GetOrdinarySum().Phi());

    sum += replace.at(i);

    hbootstd->SetPoint(i, i, bstrap->GetResidualMean(i) + bstrap->GetOrdinarySum().Phi());
    hbootstd->SetPointError(i, 0, bstrap->GetResidualStdDev(i));
  }
   
  average = sum/(Double_t)replace.size();
  Double_t bsdev = 0.;
  for( UInt_t i = 0; i < (UInt_t)replace.size(); i++)  
    bsdev += pow( replace.at(i) - average , 2);

  Double_t bsdev2 = bsdev/(Double_t)(replace.size() - 1);


  cout << " ----------- " << endl;
  cout << " Ordinary --> " 
       << " Mean        " << bstrap->GetOrdinaryMean()
       << " s           " << bstrap->GetOrdinaryStdDev()
       << " s/sqrt(n)   " << bstrap->GetOrdinaryStdDev()/sqrt( (Double_t)ntry )
       << endl;

  cout << " ---------------------- " << endl;

  Double_t bsaverage = bstrap->GetMean();
  Double_t bsstddev  = bstrap->GetStdDev();
  
  cout << " Bootstrap --> " 
       << " Mean " << bstrap->GetMean()
       << " s  "   << bstrap->GetStdDev()
       << " CL low "<< bstrap->GetCLLow() 
       << " CL up  "<< bstrap->GetCLUp() 
       << " +-  +  "<< bstrap->GetError()
       << " replaceing size " << replace.size()

       << endl;


  

  
  cout << " Deviation ++++++++++++++++++++ " << endl;

  
  cout << " Ordinary  : diff  " << bstrap->GetOrdinaryMean() - genmean
       << " ->  " << abs(bstrap->GetOrdinaryMean() - genmean) / bstrap->GetOrdinaryStdDev()
       << endl;

  
  cout << " BootStrap : 95%CL " << bstrap->GetMean() 
       << " +-  " << bstrap->GetError()
       << " ---> "
       << bstrap->GetMean() - bstrap->GetError() << " ~ " 
       << bstrap->GetMean() + bstrap->GetError()
       << endl;


  // QQ plot

  for( UInt_t i = 0; i < (UInt_t)replace.size(); i++) {
    
    Double_t percent = (Double_t)i / (Double_t)replace.size();
    Double_t z_score = ROOT::Math::gaussian_quantile(percent, bsstddev)/bsstddev;

    if( !std::isnan(z_score) && !std::isinf(z_score) )
      hbsqq->SetPoint(i,  replace.at(i), z_score );

    //    cout << z_score << " vs " << replace.at(i) << endl;
  }


  gStyle->SetOptStat(0);

  UInt_t id = 1;
  cc[0] = new TCanvas("cc0","cc0", 1200, 500);
  cc[0]->Divide(3,1);

  cc[0]->cd(id); id++;
  //  horiginal->Draw();
  horig->SetMarkerStyle(20);
  horig->SetMarkerSize(0.5);
  horig->SetMarkerColor(4);
  horig->Draw("AP");
  auto oline = new TArrow(0.,0.,(bstrap->GetOrdinarySum()).Unit().X(),(bstrap->GetOrdinarySum()).Unit().Y(),0.01,">");
  oline->SetLineColor(2);
  oline->Draw();



  cc[0]->cd(id); id++;
  hbootstd->SetLineColor(7);
  hbootstd->SetMarkerStyle(20);  
  hbootstd->SetMarkerSize(0.5);  
  hbootstd->SetMarkerColor(2);
  hbootstd->Draw("ALP");

  //  cc[1] = new TCanvas("cc1","cc1");
  cc[0]->cd(id); id++;
  hboot->Draw();
  

  // cc[0]->cd(id); id++; id++;
  // fgaus->Draw("lp");

  // id = 1;
  // cc[1] = new TCanvas("cc1","cc1", 1000, 500);
  // cc[1]->Divide(2,1);

  // cc[1]->cd(id); id++;
  // hqq->SetMarkerStyle(20);
  // hqq->SetMarkerSize(0.5);
  // hqq->SetMarkerColor(2);
  // hqq->Draw("AP");

  // cc[1]->cd(id); id++;
  // hbsqq->SetMarkerStyle(20);
  // hbsqq->SetMarkerSize(0.2);
  // hbsqq->SetMarkerColor(2);
  //  hbsqq->Draw("AP");

}


void mathdist()
{
  auto gph1 = new TGraph();

  for(UInt_t i = 0; i < 500; i++ ){

    Double_t x = (Double_t)i * 1;
    //Double_t y = ROOT::Math::normal_pdf(x, 1.);
    //    Double_t y = ROOT::Math::gaussian_quantile(x, 2.);

    Double_t y = ROOT::Math::poisson_pdf(x, 148);

    cout << "x " << x << " y " << y << endl;

    
    if( !std::isnan(y) && !std::isinf(y) ) 
      gph1->SetPoint(i-1, x, y);
  }

  gph1->Draw();

}


void mytest()
{
  std::vector< Double_t > data{61, 61, 61, 88, 89, 89, 90, 93, 93, 94, 102, 105, 108, 109, 109, 114, 115, 115, 120, 138};

  std::sort(data.begin(), data.end());
  
  cout << " mean " << TMath::Mean(data.begin(), data.end()) << endl;
  
  std::vector< Double_t >::iterator  itrmedian = data.begin() + data.size()/2 - 1;

  Double_t median = Double_t(*itrmedian);
  if( data.size()%2 == 0) median = ( median + Double_t(*(itrmedian+1)) )/2.;
			    
  cout << " median " << median << " " << Double_t(*itrmedian) << " + " << Double_t(*(itrmedian+1)) << endl;
  

  
  for( std::vector< Double_t >::iterator it = data.begin(); it != data.end(); it++ )
    cout << it-data.begin() << " : " << *it << endl;

}
