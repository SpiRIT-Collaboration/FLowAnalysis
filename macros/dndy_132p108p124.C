#include "openFlw.h"



void dndy_132p108p124()
{
  gStyle->SetOptTitle(0);

  TCanvas *cc[5];
  UInt_t   ic = 0;
  UInt_t   icol[]     = {4,2,3,7};
  TString  iopt[]     = {"","same","same","same"};
  
  const UInt_t nsys = 3;
  UInt_t m_bgn = 0;
  UInt_t m_end = nsys;


  auto ifile = new TFile("dndy_132p108p124.root");

  TH1D* hrap[nsys][5];
  TH1D* hnpart[nsys][5];
  TH1D* hmtrack[nsys];

  auto aLeg0 = new TLegend(0.12,0.7,0.42,0.9,"");
  auto aLeg1 = new TLegend(0.7 ,0.7,0.9 ,0.9,"");

  for(Int_t m = m_bgn; m < m_end; m++){

    for(UInt_t i = 0; i < 5; i++){

      TString hname  = Form("hrap%d_%d",m,i);
      hrap[m][i] = (TH1D*)gROOT->FindObject(hname);
      hrap[m][i] -> SetLineColor(icol[m]);
      //      hrap[m][i] ->SetFillStyle(3444);
      //  hrap[m][i] ->SetFillColor(icol[m]);

      hname  = Form("hnpart%d_%d",m,i);
      hnpart[m][i] = (TH1D*)gROOT->FindObject(hname);
      hnpart[m][i] -> SetLineColor(icol[m]);
      //      hnpart[m][i] ->SetFillStyle(3444);
      //hnpart[m][i] ->SetFillColor(icol[m]);
    }

    TString hname  = Form("hmtrack%d",m);
    hmtrack[m] = (TH1D*)gROOT->FindObject(hname);
    hmtrack[m] -> SetLineColor(icol[m]);
    
    aLeg0->AddEntry(hrap[m][4]  ,sysName[m],"lp");
    aLeg1->AddEntry(hnpart[m][4],sysName[m],"lp");
  }

  //------------------------


  //----- Drawing
  //----- canvas
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1400,500);
  cc[ic]->Divide(3,1);

  Double_t hmax[nsys];

  for(UInt_t ip = 2; ip < 5; ip++){
    cc[ic]->cd(ip-1); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hrap[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    cout << "ip : " << ip  << " gmax " << gmax << endl;
    
    hrap[0][ip]->SetMaximum(gmax*1.1);

    for(UInt_t m = m_bgn; m < m_end; m++) 
      hrap[m][ip]->Draw(iopt[m]); 
     

    if(ip == 4)
      aLeg0->Draw();

  }

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
  cc[ic]->Divide(2,1);

 
  for(UInt_t ip = 0; ip < 2; ip++){
    cc[ic]->cd(ip+1); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hrap[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    hrap[0][ip]->SetMaximum(gmax*1.1);
    for(UInt_t m = m_bgn; m < m_end; m++) 
      hrap[m][ip]->Draw(iopt[m]); 
  

    if(ip == 1)
      aLeg1->Draw();
  }


  // //----cc1
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,800);
  cc[ic]->Divide(3,2);
  
  for(UInt_t ip = 2; ip < 5; ip++){
    cc[ic]->cd(ip-1); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hnpart[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    hnpart[0][ip]->SetMaximum(gmax*1.1);
    for(UInt_t m = m_bgn; m < m_end; m++) 
      hnpart[m][ip]->Draw(iopt[m]);

    if(ip == 4)
      aLeg1->Draw();
  }

  io = 0;
  for(UInt_t ip = 0; ip < 2; ip++){
    cc[ic]->cd(ip+4); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hnpart[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    hnpart[0][ip]->SetMaximum(gmax*1.1);
    for(UInt_t m = m_bgn; m < m_end; m++)
      hnpart[m][ip]->Draw(iopt[m]);
    

    if(ip == 1)
      aLeg1->Draw();
  }

}


