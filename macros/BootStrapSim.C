#include <algorithm>  
#include "SetStyle.C"
TClonesArray *aArray;
TClonesArray *aFlowArray;


TCanvas *cc;
UInt_t ic = -1;

void  BootStrapSim()
{
  auto rChain = new TChain("cbmsim");
  rChain->Add("data/run2841_BTt.v40.0.root");
  rChain->Add("data/run2843_BTt.v40.0.root");

  if( rChain == NULL) return;

  rChain->SetBranchAddress("STParticle",&aArray);
  rChain->SetBranchAddress("STFlow"    ,&aFlowArray);


  auto nEntry = rChain->GetEntries();
  cout << " entry " << nEntry << endl;


  auto bs_unitP = new STBootStrap();

  auto gr_boots = new TGraphErrors();
  auto hbsPsi_o  = new TH1D("hbsPsi_o","; #Psi(original)",100,-3.14,3.14);
  auto hbsPsi_b  = new TH1D("hbsPsi_b","; #Psi(bootstrap)",100,-3.14,3.14);

  
  for(Long64_t i = 0; i < nEntry; i++){

    rChain->GetEntry(i);


    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);
    if( aflow == NULL ) continue;

    UInt_t mtrk = 0;
    bs_unitP->clear();

    TIter next(aArray);
    STParticle *aPart = NULL;

    while( (aPart = (STParticle*)next()) ) {

      if( aPart->GetReactionPlaneFlag()%2 == 1 ) {
    	
    	Double_t wt = aPart->GetRPWeight();
    	TVector2 ptr= aPart->GetRotatedPt().Unit();
    	
	bs_unitP->Add(wt * ptr);
    	mtrk++;
      }
    }
    

    if( mtrk > 0 ) {
      bs_unitP->BootStrapping(1000);    

      hbsPsi_b->Fill(TVector2::Phi_mpi_pi(bs_unitP->GetMean()));
      hbsPsi_o->Fill(aflow->unitP.Phi());

      
  }

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hbsPsi_b->Draw();

  hbsPsi_o->SetLineColor(4);
  hbsPsi_o->Draw("same");
  
}

