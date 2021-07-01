#include "DoFlow.h"

void DrawUt(TH2D *h2, UInt_t sel, TString opt="", Color_t scol=1);

void PlotUt()
{
  TH2D* h245;
  TH2D* h2135;

   TFile* infile1 = TFile::Open("data/finYPt_108Sn_3He.v52.15.98.root");
   h245 = (TH2D*)infile1->Get("hyutacpcrr");  
   h245->SetName("hyutacpcrr45");
   DrawUt(h245,2,"",2);

   h245 = (TH2D*)infile1->Get("hyutacp");
   h245->SetName("hyutacp45");
   DrawUt(h245,2,"same",8);

   TFile *infile2 = TFile::Open("data/finYPt_108Sn_3He.v52.15.99.root");
   h2135 = (TH2D*)infile2->Get("hyutacp");
   h2135->SetName("hyutacp135");
   DrawUt(h2135,2,"same",7);

   h2135 = (TH2D*)infile2->Get("hyutacpcrr");  
   h2135->SetName("hyutacpcrr135");
   DrawUt(h2135,2,"same",4);

  // TFile *infile1 = TFile::Open("data/rootfiles/UnfoldedLCPSpectra_45_40to55.slim.root");
  // TH2D* h2c = (TH2D*)infile1->Get("h2UtYEff_108Sn_3He_mbin0_iter1");
  // h2c->SetName("h2UtYEff_108Sn_3He_mbin0_iter1_45");
  // DrawUt(h2c,1);

  // TFile *infile2 = TFile::Open("data/rootfiles/UnfoldedLCPSpectra_135_40to55.slim.root");
  // h2c = (TH2D*)infile2->Get("h2UtYEff_108Sn_3He_mbin0_iter1");
  // h2c->SetName("h2UtYEff_108Sn_3He_mbin0_iter1_135");
  // DrawUt(h2c,1,"same",2);
}


void DrawUt(TH2D *h2, UInt_t sel, TString opt, Color_t scol)
{
  TString hname = h2->GetName();
  cout << " DrawUt -> " << hname << endl;

  if( opt != "same" ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    h2->Draw("colz");
  }
    
  Double_t *yrange = &yrange1nrm[1];
  UInt_t ybin = ybin1;

  if( sel == 2 ) {
    yrange = &yrange2nrm[1];
    ybin = ybin2;
  }


  TH1D *hacp;
  if( opt != "same" ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1600,1200);
    cc->Divide(4, ybin/4);
  }

  for( UInt_t i = 0; i < ybin; i++ ) {
    cc->cd(i+1);

    auto firstbin = h2->GetXaxis()->FindBin(*(yrange+i));
    auto lastbin  = h2->GetXaxis()->FindBin(*(yrange+i+1));

    //    cout << " ybin " << i << " " << *(yrange+i) << " y " << firstbin << " ~ " << *(yrange+i+1) << endl;

    hacp = (TH1D*) h2->ProjectionY(hname+Form("_%d_%d",sel,i), firstbin, lastbin);
    hacp -> SetTitle(Form("%3.2f < y < %3.2f" , *(yrange+i), *(yrange+i+1))); 
    if( hacp->GetEntries() > 0) {
      //      hacp -> SetNormFactor(50);
      //      hacp -> Scale(1./hacp->GetEntries() );

      hacp->SetLineColor(scol);
      
      hacp -> Draw(opt);
    }
    else 
      cout << " NO entry in ybin " << i << endl;

  }
}

