#ifndef ImageProcessing_h
#define ImageProcessing_h 1
namespace ImageProcessing
{
	TH2D* hSmoothened;

	TH2D* BoxBlur(TH2*);
	TH2D* GaussianBlur(TH2*,Int_t,Double_t);
	TH2D* MedianBlur(TH2*,Int_t);
//	TH2D* BilateralBlur(TH2*,Int_t,Double_t);
}

TH2D* ImageProcessing::BoxBlur(TH2* h)
{
	if(!h) return nullptr;

	hSmoothened = (TH2D*)h->Clone(Form("%s_BoxBlur",h->GetName()));

	return hSmoothened;

}

TH2D* ImageProcessing::GaussianBlur(TH2* h, Int_t kernelSize=3, Double_t sigma=1./TMath::Sqrt(2.*TMath::Log(2)))
{
  if(!h||kernelSize%2!=1){
    std::cout<<" GaussianBlur didn't work. h is nullptr or even kernelSize."<<std::endl;
    return nullptr;
  }

  hSmoothened = (TH2D*)h->Clone(Form("%s_GaussBlur",h->GetName()));
	
  Double_t totW  = 0.;
  auto w    = new Double_t*[kernelSize];
  auto cont = new Double_t*[kernelSize];
  auto err  = new Double_t*[kernelSize];
  for(auto xxbin: ROOT::TSeqI(kernelSize)){
    w[xxbin]    = new Double_t[kernelSize];
    cont[xxbin] = new Double_t[kernelSize];
    err[xxbin]  = new Double_t[kernelSize];
    for(auto yybin: ROOT::TSeqI(kernelSize))
      w[xxbin][yybin] = TMath::Gaus(xxbin-(kernelSize-1)/2,0,sigma,kTRUE)*TMath::Gaus(yybin-(kernelSize-1)/2,0.,sigma,kTRUE);
		
  }
  for(int xbin=0; xbin<=h->GetNbinsX()+1; ++xbin)for(int ybin=0; ybin<=h->GetNbinsY()+1; ++ybin){	// loop all over cells of input TH2.
      totW = 0.;
      for(auto xxbin: ROOT::TSeqI(kernelSize))for(auto yybin: ROOT::TSeqI(kernelSize)){	// make weights for designated kernel size.
	  int kx = xbin-(kernelSize-1)/2+xxbin;
	  int ky = ybin-(kernelSize-1)/2+yybin;
	  cont[xxbin][yybin]=0; err[xxbin][yybin]=0;
	  if(kx>0 && kx<h->GetNbinsX()+1 && ky>0 && ky<h->GetNbinsY()+1){
	    auto gbin = h->GetBin(kx,ky);  // get global bin number.
	    cont[xxbin][yybin] = h->GetBinContent(gbin);
	    err[xxbin][yybin]  = h->GetBinError(gbin);
	    totW += TMath::Gaus(xxbin-(kernelSize-1)/2,0,sigma,kTRUE)*TMath::Gaus(yybin-(kernelSize-1)/2,0.,sigma,kTRUE);
	  }
	}
      double scont = 0, serr = 0.;
      for(auto xxbin: ROOT::TSeqI(kernelSize))for(auto yybin: ROOT::TSeqI(kernelSize)){
	  scont += cont[xxbin][yybin]*w[xxbin][yybin]/totW;
	  serr  += err[xxbin][yybin]*w[xxbin][yybin]/totW;
	}
      auto gbin = hSmoothened->GetBin(xbin,ybin);  // get global bin number.
      hSmoothened->SetBinContent(gbin,scont);
      hSmoothened->SetBinError(gbin,serr);
    }
	
  for(auto xxbin: ROOT::TSeqI(kernelSize)){
    delete[] w[xxbin];	delete[] cont[xxbin];   delete[] err[xxbin];
  }
  delete[] w; delete[] cont; delete[] err;

  return hSmoothened;
}

TH2D* ImageProcessing::MedianBlur(TH2* h, Int_t kernelSize=3)
{
	if(!h) return nullptr;

	hSmoothened = (TH2D*)h->Clone(Form("%s_MedianBlur",h->GetName()));

	return hSmoothened;
}


/*
TH2D* ImageProcessing::BilateralBlur(TH2* h, Int_t kernelSize=3, Double_t sigma)
{
	if(!h||kernelSize%2!=1) return nullptr;

	hSmoothened = (TH2D*)h->Clone(Form("%s_BilateralBlur",h->GetName()));

	return hSmoothened;

}
*/


#endif
