TCanvas *cc[6];
auto fv1v2 = new TF1("fv1v2","1+2.*[0]*cos(x) + 2.*[1]*cos(2.* x)",-3.14, 3.14); 
UInt_t itry = 1500;

STBootStrap *bstrap = new STBootStrap(1);

TH1D *hphidist;
TH2D *hpxy;
TGraphErrors *gstrap;
TGraphErrors *gsdiv;
TGraphErrors *gsdiff;
TH1D *hphi;
TH1D *hphio;
TH1D *hstd;

TH1D *hbsPhiDist;

void  SingleBootStrap(UInt_t np = 100);
void  RndBootStrap(UInt_t np = 100, UInt_t ith = 0);
void  RPResolution();

void BootStrapTest()
{
  fv1v2->SetParameter(0, 0.1);
  fv1v2->SetParameter(1, 0.);

  //  RPResolution();
}

void RndBootStrap(UInt_t np = 100, UInt_t ith = 0)
{

  TDatime dtime;
  Int_t itime = dtime.Get();
  TRandom3 gRandom(0);
  gRandom.SetSeed(itime);

  if( ith == 0 ) {
    hphidist = new TH1D("hphidist","", 100, -3.1, 3.1);
    hpxy     = new TH2D("hpxy","",100,-1.2, 1.2, 100, -1.2, 1.2);
  }


  UInt_t npart = np;

  bstrap->clear();
  TVector2 vcsum = TVector2(0, 0);

  for(UInt_t i = 0; i < npart; i++) {
      
    auto phi = fv1v2->GetRandom();
      
    if( ith == 0 )
      hphidist->Fill(phi);

    TVector2 vec = TVector2( cos(phi), sin(phi) ).Unit(); 

    if( ith == 0 )
      hpxy->Fill(cos(phi), sin(phi));

    bstrap->Add(vec);

    vcsum += vec;
  }

  auto vsphi = TVector2::Phi_mpi_pi( vcsum.Phi() );
  bstrap->BootStrapping(itry);

  hbsPhiDist = bstrap->GetBSPhiDistribution();
}

void SingleBootStrap(UInt_t np = 100)
{

  RndBootStrap(np, 0);

  Double_t vsphi = bstrap->GetOriginalPhi();

  gstrap = new TGraphErrors();
  gstrap->SetTitle(";Number of Trial; Reaction Plane Orientaion [rad]");

  gsdiv  = new TGraphErrors();
  gsdiff = new TGraphErrors();

  hphi  = new TH1D("hphi" ,"BootStrap method",100,vsphi-0.3, vsphi+0.3);
  hphio = new TH1D("hphio","Original phi"    ,100,vsphi-0.3, vsphi+0.3);
  hstd  = new TH1D("hstd" ,"BootStrap method standard deviation",100,0.,1.2);

  for(UInt_t it = 1; it < itry; it++){

    if( it%100 == 0 ) std::cout << "Processing.. "<< it << " " << Double_t(it)/Double_t(itry)*100. << std::endl;

    gstrap->SetPoint(it-1, (Double_t)it, bstrap->GetResidualMean(it));
    gstrap->SetPointError(it-1, 0., bstrap->GetResidualStdDev(it));

    hphi->Fill( bstrap->GetResidualMean(it) );
    hstd->Fill( bstrap->GetResidualStdDev(it) );

    hphio->Fill( vsphi );
    //    auto diff = bstrap->GetResidualMean(it) - bstrap->GetResidualMean(it-1);
    auto diff = bstrap->GetResidualStdDev(it);
    gsdiv->SetPoint(it-1, (Double_t)it, diff);

    gsdiff->SetPoint(it-1, (Double_t)it, bstrap->GetResidualMean(it)-vsphi);
    //   gsdiff->SetPointError(it-1, 0., bstrap->GetResidualStdDev(it) );
  }
  
  cout << " BST mean " << bstrap->GetMean() << endl;
  cout << " vsum Phi " << vsphi << endl;


  cc[0] = new TCanvas("cc0","cc0");
  gstrap->SetMarkerStyle(20);
  gstrap->SetMarkerColor(2);
  gstrap->Draw("alp");

  auto aLine = new TLine(0., vsphi, gstrap->GetXaxis()->GetXmax(), vsphi);
  aLine->Draw();

  cc[1] = new TCanvas("cc1","cc1");
  hphi->Draw();
  hphio->SetLineColor(2);
  hphio->Draw("same");

  // cc[2] = new TCanvas("cc2","cc2");
  // hstd->Draw();

  cc[3] = new TCanvas("cc3","cc3");
  hpxy->SetMarkerStyle(20);
  hpxy->SetMarkerSize(0.5);
  hpxy->SetMarkerColor(2);
  hpxy->Draw();

  Double_t x = (bstrap->GetVectorSum()).Unit().X();
  Double_t y = (bstrap->GetVectorSum()).Unit().Y();
  auto vline = new TLine(0.,0., x, y);
  vline->SetLineColor(4);
  vline->Draw();

  // auto f1 = new TF1("f1","[0]/sqrt(x)",0,10000);

  // gStyle->SetOptFit(1111);
  // gsdiv->Fit("f1");

  // gsdiv->SetMinimum(0);
  // gsdiv->Draw("ALP");

   cc[4] = new TCanvas("cc4","cc4");
   hbsPhiDist->Draw();


   // gsdiff->Draw("alp");

  cc[5] = new TCanvas("cc5","cc5");
  hphidist->Draw();

}

void SaveCanvas(TString str = "")
{

  Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
  for(Int_t i = 0; i < iCanvas; i++) {
    TString scname = "bs" + str + Form("%d.png",i);
    gROOT->GetListOfCanvases()->At(i)->SaveAs(scname);
  }
}


void RPResolution()
{
  
  auto hPhiAll  = new TH1D("hPhiAll","" ,100, -3.14, 3.14);
  auto hPhiAllo = new TH1D("hPhiAllo","",100, -3.14, 3.14);
  
  
  for( UInt_t n = 0; n < 1000; n++ ){
    
    RndBootStrap(100, n+1);
    hPhiAll->Fill(bstrap->GetMean());
		  
    hPhiAllo->Fill(bstrap->GetOriginalPhi());
  }    

  cc[0] = new TCanvas("cc0","cc0");
  hPhiAll->Draw("");

  cc[1] = new TCanvas("cc1","cc1");
  hPhiAllo->Draw("");
  
}
