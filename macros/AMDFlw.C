//#include "../analysisformat/AMDParticle.hh"

TChain *fChain = NULL;

TH2D *h2[10];
TH2D *hangle;
TCanvas *cc[20];
TCanvas *ccv;
UInt_t ic  = 0;
TString pfname[] = {"pi-","pi+","proton","deuteron" ,"triton"};
TString pname[] = {"pi-","pi+","prot","deut" ,"trit","neut"};
UInt_t PID[] = {211,    211,    2212, 1000010020, 1000010030,  2112};

TString printHeader;
TString beamHeader;

const UInt_t nybin = 8;
std::vector< std::vector< Double_t > >vrapd(nybin);
TH1D* hphi1[nybin];
TH1D* hphi2[nybin];

void SaveCanvas(TString fopt = "", Int_t isel=-1)
{
  if(isel > -1)
    gROOT->GetListOfCanvases()->At(isel)->SaveAs(printHeader+fopt+Form("_%d",isel)+".png");

  else {
    Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
    for(Int_t i = 0; i < iCanvas; i++)
      gROOT->GetListOfCanvases()->At(i)->SaveAs(printHeader+fopt+Form("_%d",i)+".png");
  }
}

void PlotCosv1v2(UInt_t ipid = 2)
{
  gStyle->SetOptStat(0);

  printHeader += pname[ipid];

  TFile *froot = new TFile("data/"+printHeader+".root","recreate");

  auto hyptacp = new TH2D("hyptacp", beamHeader+":"+pname[ipid]+";Rapidity; Pt[MeV/c]" ,100, -0.65, 0.65, 100, 0., 1000);
  auto hphi    = new TH1D("hphi"   , "phi", 100, -3.2, 3.2);


  Double_t yrange1[] = {-0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35,  0.5};
  const UInt_t ybin1 = sizeof(yrange1)/sizeof(Double_t);

  Double_t yrange2[] = {-0.3, -0.1, -0.05,  0.05, 0.2, 0.3, 0.5};
  const UInt_t ybin2 = sizeof(yrange2)/sizeof(Double_t);

  std::vector< std::vector< Double_t > > cosv1x(ybin1);
  std::vector< std::vector< Double_t > > cosv2x(ybin2);
  Double_t  cosv1[ybin1];
  Double_t  cosv2[ybin2];
  Double_t  sinv1[ybin1];
  Double_t  sinv2[ybin2];

  for( UInt_t i = 0; i < ybin1; i++) {
    cosv1[i]   = 0.;
    sinv1[i]   = 0.;
    cosv1x[i].clear();
  }

  for( UInt_t i = 0; i < ybin2; i++) {
    cosv2[i]   = 0.;
    sinv2[i]   = 0.;
    cosv2x[i].clear();
  }



  TClonesArray *aAMDArray = NULL;
  fChain->SetBranchAddress("AMDParticle",&aAMDArray);


  Int_t nEntry = fChain -> GetEntries();

  for(Int_t i = 0; i < nEntry; i++){

    fChain -> GetEntry(i);

    Int_t ntrack = aAMDArray->GetEntries();

    TIter next(aAMDArray);
    
    AMDParticle* aAMDParticle = NULL;

    while( (aAMDParticle = (AMDParticle*)next()) ) {

      if( aAMDParticle->fPdg == PID[ipid] ) {

	// if( aAMDParticle->fMomentum.Theta() > 0.75 &&
	//     abs(aAMDParticle->fMomentum.Phi()) < 2 && abs(aAMDParticle->fMomentum.Phi()) > 1) continue;

	Double_t rapid = aAMDParticle->fMomentum.Rapidity();
	Double_t dphi  = aAMDParticle->fMomentum.Phi();


	hyptacp->Fill(rapid, aAMDParticle->fMomentum.Pt());

        UInt_t irapid = ybin1 - 1;
        for( UInt_t i = 0; i < ybin1; i++){
          if(rapid < yrange1[i]){
            irapid = i;
            break;
          }
        }

        cosv1x[irapid].push_back( rapid );
        cosv1[irapid] += cos(dphi);
        sinv1[irapid] += sin(dphi);

	irapid = ybin2 - 1;
	for( UInt_t i = 0; i < ybin2; i++){
          if(rapid < yrange2[i]){
            irapid = i;
            break;
          }
        }

	if( irapid == 5 )
	  hphi->Fill(2.*dphi);

        cosv2x[irapid].push_back( rapid );
        cosv2[irapid] += cos(2.*dphi);
        sinv2[irapid] += sin(2.*dphi);
      }
    }
  }


  TGraphErrors *gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  gv_v1->SetTitle(pname[ipid]+";Rapidity; v1");
  TGraphErrors *gv_v2 = new TGraphErrors();
  gv_v2->SetName("gv_v2");
  gv_v2->SetTitle(pname[ipid]+";Rapidity; v2");

  UInt_t kl = 0;
  for(UInt_t kn = 0; kn < ybin1; kn++){

    if( cosv1x[kn].size() == 0 ) continue;

    Double_t rapm  = TMath::Mean(cosv1x[kn].begin(), cosv1x[kn].end());
    Double_t rape  = TMath::StdDev(cosv1x[kn].begin(), cosv1x[kn].end());

    Double_t yv1  = cosv1[kn] / (Double_t)cosv1x[kn].size();
    Double_t yv1e = sinv1[kn] / (Double_t)cosv1x[kn].size();

    Double_t yv1c  = yv1;
    Double_t yv1ce = abs(yv1e);

    if( !std::isnan( yv1c ) && !std::isinf( yv1c ) ) {
      gv_v1->SetPoint( kl,      rapm, yv1c);
      gv_v1->SetPointError( kl, rape, yv1ce);
      kl++;
    }
  }

  kl = 0;
  for(UInt_t kn = 0; kn < ybin2; kn++){

    if( cosv2x[kn].size() == 0 ) continue;

    Double_t rapm  = TMath::Mean(cosv2x[kn].begin(), cosv2x[kn].end());
    Double_t rape  = TMath::StdDev(cosv2x[kn].begin(), cosv2x[kn].end());

    Double_t yv2  = cosv2[kn] / (Double_t)cosv2x[kn].size();
    Double_t yv2e = sinv2[kn] / (Double_t)(cosv2x[kn].size());

    Double_t yv2c  = yv2;
    Double_t yv2ce = abs(yv2e);

    if( !std::isnan( yv2c ) && !std::isinf( yv2c ) ) {
      gv_v2->SetPoint( kl,    rapm, yv2c);
      gv_v2->SetPointError( kl, rape, yv2ce);
      kl++;

      cout << " cos " << cosv2[kn] << " & " << sinv2[kn] << " size " << kn << " " << cosv2x[kn].size() 
	   << " rapidty " << rapm << " + " << rape 
	   << " <cos> " << yv2c << " +- " << yv2ce
	   << endl;
      

    }
  }



  ic++; 
  ccv = new TCanvas(Form("ccv%d",ic),Form("ccv%d",ic));
  hyptacp->Draw("colz");

  ic++;
  ccv = new TCanvas(Form("ccv%d",ic),Form("ccv%d",ic));
  gv_v1->Draw("ALP");

  auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();

  ic++; 
  ccv = new TCanvas(Form("ccv%d",ic),Form("ccv%d",ic));
  gv_v2->Draw("ALP");

  gv_v1->Write();
  gv_v2->Write();
  froot->Write();

}



void AMDFlw(UInt_t isel = 0, UInt_t ipid = 2)
{
  TString fname[4];

  switch(isel) {
  case 0:
    fname[0] = "132Sn124Sn270AMeV_cluster_SLy4-L108_bimp01.root";
    fname[1] = "132Sn124Sn270AMeV_cluster_SLy4-L108_bimp13.root";
    fname[2] = "132Sn124Sn270AMeV_cluster_SLy4-L108_bimp35.root";
    fname[3] = "132Sn124Sn270AMeV_cluster_SLy4-L108_bimp57.root";
    break;

  case 1:
    fname[0] = "132Sn124Sn270AMeV_cluster_SLy4_bimp01.root";
    fname[1] = "132Sn124Sn270AMeV_cluster_SLy4_bimp13.root";
    fname[2] = "132Sn124Sn270AMeV_cluster_SLy4_bimp35.root";
    fname[3] = "132Sn124Sn270AMeV_cluster_SLy4_bimp57.root";
    break;

  case 2:
    fname[0] = "132Sn124Sn270AMeV_nocluster_SLy4_bimp01.root";
    fname[1] = "132Sn124Sn270AMeV_nocluster_SLy4_bimp13.root";
    fname[2] = "132Sn124Sn270AMeV_nocluster_SLy4_bimp35.root";
    fname[3] = "132Sn124Sn270AMeV_nocluster_SLy4_bimp57.root";
    break;

  case 3:
    fname[0] = "132Sn124Sn270AMeV_nocluster_SLy4-L108_bimp01.root";
    fname[1] = "132Sn124Sn270AMeV_nocluster_SLy4-L108_bimp13.root";
    fname[2] = "132Sn124Sn270AMeV_nocluster_SLy4-L108_bimp35.root";
    fname[3] = "132Sn124Sn270AMeV_nocluster_SLy4-L108_bimp57.root";
    break;

  case 4:
    fname[0] = "108Sn112Sn270AMeV_nocluster_SLy4-L108_bimp01.root";
    fname[1] = "108Sn112Sn270AMeV_nocluster_SLy4-L108_bimp13.root";
    fname[2] = "108Sn112Sn270AMeV_nocluster_SLy4-L108_bimp35.root";
    fname[3] = "108Sn112Sn270AMeV_nocluster_SLy4-L108_bimp57.root";
    break;

  case 5:
    fname[0] = "108Sn112Sn270AMeV_nocluster_SLy4_bimp01.root";
    fname[1] = "108Sn112Sn270AMeV_nocluster_SLy4_bimp13.root";
    fname[2] = "108Sn112Sn270AMeV_nocluster_SLy4_bimp35.root";
    fname[3] = "108Sn112Sn270AMeV_nocluster_SLy4_bimp57.root";
    break;

  case 6:
    fname[0] = "108Sn112Sn270AMeV_cluster_SLy4-L108_bimp01.root";
    fname[1] = "108Sn112Sn270AMeV_cluster_SLy4-L108_bimp13.root";
    fname[2] = "108Sn112Sn270AMeV_cluster_SLy4-L108_bimp35.root";
    fname[3] = "108Sn112Sn270AMeV_cluster_SLy4-L108_bimp57.root";
    break;

  case 7:       
    fname[0] = "108Sn112Sn270AMeV_cluster_SLy4_bimp01.root";
    fname[1] = "108Sn112Sn270AMeV_cluster_SLy4_bimp13.root";
    fname[2] = "108Sn112Sn270AMeV_cluster_SLy4_bimp35.root";
    fname[3] = "108Sn112Sn270AMeV_cluster_SLy4_bimp57.root";
  }                    


  beamHeader  = fname[0](0,5);
  printHeader = beamHeader+"-"+fname[0](18, fname[0].Length()-30 )+"-";

  fChain = new TChain("amdTree");

  TString ffdir = "/cache/scr/spirit/mizuki/SpiRITAnalysis/macros/data/AMD/";
  for(UInt_t i = 0; i < 3; i++){
    if( gSystem->FindFile(ffdir, fname[i]) ) {
	fChain->Add(fname[i]);
	std::cout << " File : " << fname[i] << " is added. " << std::endl;
    }
  }

  // Booking(ipid);
  // Analysis(ipid);
  // getv1v2(ipid);


  PlotCosv1v2(ipid);

  //  SaveCanvas("_"+pname[ipid]);
}




