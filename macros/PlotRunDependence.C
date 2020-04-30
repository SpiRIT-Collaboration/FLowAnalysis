#include "openRunAna.C"

//RUN={$RNF132} root PlotRunDependence.C

void PlotRunDependence()
{
  gROOT->Reset();
  openRunAna();

  if(rChain != NULL)
    LOG(INFO) << " DoFLow_adv: System " << isys << "  -> " << sysName << FairLogger::endl;
  else
    exit(0);

  TCanvas *cc;
  UInt_t ic = 0;

  auto hdedxRun = new TH2D("hdedxRun",";Run; dEdx",115,0,115,  200,0., 150);

  Long64_t nEntry = SetBranch();

  //--------------------------------------------------                                                              
  //--- Event Loop                                                                                                  
  //--------------------------------------------------                                                              
  Int_t irun = -1;
  UInt_t   srun = 0;
  for(Long64_t i = 0; i < nEntry; i++){

    ShowProcess(i);

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);
    auto run  = aflow->run;
    
    if( run != srun ) {
      irun++;
      srun = run;
      hdedxRun->GetXaxis()->SetBinLabel(irun, Form("%04d", srun));
    }


    TIter next(aArray);
    STParticle *aPart = NULL;

    //--------------------------------------------------                                                            
    //----- Main loop                                                                                               
    //--------------------------------------------------                                                            
    UInt_t mtk = 0;
    while( (aPart = (STParticle*)next()) ) {

      if( aPart->GetPID() == 2212 && aPart->GetGoodTrackFlag()==1111 && aPart->GetP()>500 && aPart->GetP()<700) {
	auto dedx = aPart->GetdEdx();

	hdedxRun->Fill((Double_t)irun, dedx);
      }
    }
  }
    
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
    hdedxRun->Draw("colz");
}
