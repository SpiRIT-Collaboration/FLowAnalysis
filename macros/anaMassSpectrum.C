#include "myFunction.h"
using namespace BetheBlochFitter;

void anaMassSpectrum()
{
  auto chain = new TChain("anaTree");
  chain->Add("rootfiles/anaTree/run29*.root");
  //chain->Add("rootfiles/anaTree/run290*.root");
  //chain->Add("rootfiles/anaTree/run2901.root");
  //chain->Add("rootfiles/anaTree/run2902.root");
  Int_t runNo;
  TClonesArray* part = nullptr;
  chain->SetBranchAddress("run",&runNo);
  chain->SetBranchAddress("part",&part);

  TCut cutVertex = "TMath::Abs(vertex.Z()+12.8)<=3.&&TMath::Abs(vertex.X())<=15.&&TMath::Abs(vertex.Y()+226.06)<=20.";
  chain->Project("list","",cutVertex,"entrylist");
  auto list = (TEntryList*)gDirectory->FindObject("list");
  chain->SetEntryList(list);

  Double_t fitterPara[5];
  auto fitFile = new TFile("rootfiles/cut/BBFitter.root");
  auto fit = (TF1 *) fitFile -> Get("fit_proton");

  fitterPara[0] = fit -> GetParameter(0);
  fitterPara[1] = fit -> GetParameter(1);

  auto out = new TFile("rootfiles/mass.root","recreate");
  TTree* tree = new TTree("mass","");
  Double_t dedx, mom, mass;
  Int_t q;
  tree->Branch("run",&runNo,"run/I");
  tree->Branch("q",&q,"q/I");
  tree->Branch("dedx",&dedx,"dedx/D");
  tree->Branch("mom",&mom,"mom/D");
  tree->Branch("mass",&mass,"mass/D");

  TH1* h1Mass[4];
  Int_t nEvent[4]={};
  for(Int_t i=0; i<4; i++)
    h1Mass[i] = new TH1D(Form("h1Mass%d",i),"",500,0,6000);

  Int_t iTree=0;
  for(auto i: ROOT::TSeqL(list->GetN())){
    auto treeEntry  = list->GetEntryAndTree(i,iTree);
    auto chainEntry = treeEntry + chain->GetTreeOffset()[iTree];
    chain->GetEntry(chainEntry);

    Int_t id = -1;
    if(runNo>2200&&runNo<2510) id = 0;
    if(runNo>2520&&runNo<2660) id = 1;
    if(runNo>2800&&runNo<3040) id = 2;
    if(runNo>3000&&runNo<3300) id = 3;
    if(id==-1)continue;
    nEvent[id]++;

    for(auto j: ROOT::TSeqL(part->GetEntries())){
      auto t = (myParticle*)part->At(j);

      if(t->GetDist()>10||t->GetNCluster()<20) continue;
      //if(t->GetDist()>20||t->GetNCluster()<20||t->GetNCluster()/t->GetNTheoCluster()<0.7) continue;

      dedx = t->GetdEdx();
      q = t->GetPolarity();
      mom = t->GetRigidity().Mag();
      fitterPara[2] = q*mom;
      fitterPara[3] = q;
      fitterPara[4] = dedx;
      mass = CalcMass(fitterPara);

      h1Mass[id]->Fill(mass);

      tree->Fill();

    }



  }

  for(Int_t i=0; i<4; i++){
    h1Mass[i]->SetLineColor(i+3);
    h1Mass[i]->Scale(1./nEvent[i]);
  }

  out->Write();



}
