
void SaveCanvas(TString fopt = "", Int_t isel=-1)
{
  TString printHeader = "cp_" + fopt;

  if(isel > -1)
    gROOT->GetListOfCanvases()->At(isel)->SaveAs(printHeader+Form("_%d",isel)+".png");

  else {
    
    TString mdir = printHeader + gSystem->Now().AsString();
    gSystem->mkdir(mdir);
    gSystem->cd(mdir);
    
    Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
    for(Int_t i = 0; i < iCanvas; i++)
      gROOT->GetListOfCanvases()->At(i)->SaveAs(printHeader+Form("_%d",i)+".png");

    gSystem->cd("..");
    std::cout << " Figures were saved in " << mdir << std::endl;

  }
}


