//#include "../analysisformat/AMDParticle.hh"

TChain *fChain = NULL;

TH2D *h2[10];
TH2D *hangle;
TCanvas *cc[20];
UInt_t ic  = 0;
TString pname[] = {"prt","neut"};
TString pfname[] = {"proton","neutron"};
UInt_t PID[] = {2212, 2112};

TString printHeader;

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

void Analysis(UInt_t ipid)
{


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

	Double_t rapidity = 0.5 * log( ( aAMDParticle->fMomentum.E() + aAMDParticle->fMomentum.Pz() ) /
				       ( aAMDParticle->fMomentum.E() - aAMDParticle->fMomentum.Pz() ) ); 

	h2[0]->Fill( rapidity , aAMDParticle->fMomentum.Pt());
	h2[1]->Fill( aAMDParticle->fMomentum.X(), aAMDParticle->fMomentum.Y());
	h2[2]->Fill( rapidity , aAMDParticle->fMomentum.X());

	h2[3]->Fill( aAMDParticle->fMomentum.Phi(), rapidity );
	h2[4]->Fill( TVector2::Phi_mpi_pi(2.*aAMDParticle->fMomentum.Phi()), rapidity );

	hangle->Fill( aAMDParticle->fMomentum.Theta(), aAMDParticle->fMomentum.Phi() );

	for(UInt_t i = 0; i < nybin; i++) {
	  Double_t upbin = h2[3]->GetYaxis()->GetBinUpEdge(i);
	  if( rapidity < upbin ) {
	    vrapd[i].push_back(rapidity);
	    hphi1[i]->Fill( aAMDParticle->fMomentum.Phi() );
	    hphi2[i]->Fill( TVector2::Phi_mpi_pi(2.*aAMDParticle->fMomentum.Phi()) );
	    break;
	  }
	}
      }
    }
  }

  // for(UInt_t i = 0; i < 4; i++){
  //   cc[i] = new TCanvas(Form("cc%d",i),Form("cc%d",i)); 
  //   h2[i]->Draw("colz");
  // }

}

void getv1v2(UInt_t ipid)
{

  TH1D *hydphi1[nybin];
  TH1D *hydphi2[nybin];
  Double_t nbin1 = (Double_t)h2[3]->GetXaxis()->GetNbins();
  Double_t nbin2 = (Double_t)h2[4]->GetXaxis()->GetNbins();

  auto *fcos1 = new TF1("fcos1","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());

  auto gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  gv_v1->SetTitle(pfname[ipid]+"; Rapidity ; v1 ");

  auto gv_v2 = new TGraphErrors();
  gv_v2->SetName("gv_v2");
  gv_v2->SetTitle(pfname[ipid]+"; Rapidity ; v2 ");

  
  cc[ic]   = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); 
  cc[ic+1] = new TCanvas(Form("cc%d",ic+1),Form("cc%d",ic+1)); 

  cc[ic]  ->Divide(2, nybin/2);
  cc[ic+1]->Divide(2, nybin/2);

  UInt_t id  = 1;
  UInt_t idd = 1;


  UInt_t il = 0;
  for(UInt_t jn = 0; jn < nybin; jn++){
    
    if( vrapd[jn].size() == 0 ) continue;

    Double_t rpd  = TMath::Mean(vrapd[jn].begin(), vrapd[jn].end());
    Double_t rpde = TMath::StdDev(vrapd[jn].begin(), vrapd[jn].end());
    
    
    cc[ic]->cd(id); id++;
    //    hydphi1[jn] = (TH1D*)h2[3]->ProjectionX((TString)Form("hydphi1_%d",jn+1),jn+1, jn+1,"eo");
    hydphi1[jn] = hphi1[jn];
    
    Double_t scf = hydphi1[jn]->GetEntries();
    hydphi1[jn]->Scale(nbin1/scf);
    hydphi1[jn]->GetXaxis()->SetRange(2, 99);
    
    hydphi1[jn]->Fit("fcos1","","",-3.,3.);
    Double_t p0   = fcos1->GetParameter(0);
    Double_t p1   = fcos1->GetParameter(1);
    Double_t p0e  = fcos1->GetParError(0);
    Double_t p1e  = fcos1->GetParError(1);
    Double_t v1   = p1;
    Double_t v1e  = p1e; //sqrt( pow(v1,2)*( pow(p1e/p1,2) + pow(p0e/p0,2) )) ;
      
    gv_v1->SetPoint(il, rpd, v1);
    gv_v1->SetPointError(il, rpde, v1e);
      
    std::cout << setw(5) << jn << " w  c : " << setw(12)
	      << rpd  << " +- " << rpde 
	      <<  " v1 " << setw(12) << v1 << " +- " << setw(10) << v1e << std::endl;
    std::cout << " ----------------------------------" << std::endl;



    cc[ic+1]->cd(idd); idd++;
    //    hydphi2[jn] = (TH1D*)h2[4]->ProjectionX((TString)Form("hydphi2_%d",jn+1),jn+1, jn+1,"eo");
    hydphi2[jn] = hphi2[jn];

    scf = hydphi2[jn]->GetEntries();
    hydphi2[jn]->Scale(nbin2/scf);
    hydphi2[jn]->GetXaxis()->SetRange(2, 99);

    
    hydphi2[jn]->Fit("fcos1","","",-3.,3.);
    p0   = fcos1->GetParameter(0);
    p1   = fcos1->GetParameter(1);
    p0e  = fcos1->GetParError(0);
    p1e  = fcos1->GetParError(1);
    Double_t v2   = 0.5*p1/p0;
    Double_t v2e  = 0.5*sqrt( pow(v2,2)*( pow(p1e/p1,2) + pow(p0e/p0,2) )) ;
      
    gv_v2->SetPoint(il, rpd, v2);
    gv_v2->SetPointError(il, rpde, v2e);

    std::cout << setw(5) << jn << " w  c : " << setw(12)
	      << rpd  << " +- " << rpde 
	      <<  " v2 " << setw(12) << v2 << " +- " << setw(10) << v2e << std::endl;
    std::cout << " ==================================" << std::endl;

    il++;
      
  }

  ic += 2;
  // cc[ic]   = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;
  // h2[3]->Draw("colz");

  // cc[ic]   = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;
  // h2[4]->Draw("colz");


  cc[ic]   = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;
  gv_v1->Draw("ALP");

  gv_v1->Write();
  auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();



  cc[ic]   = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;
  gv_v2->Draw("ALP");

  gv_v2->Write();
  auto LineV2 = new TLine(0.,gv_v2->GetYaxis()->GetXmin(), 0., gv_v2->GetYaxis()->GetXmax());
  auto LineH2 = new TLine(gv_v2->GetXaxis()->GetXmin(),    0., gv_v2->GetXaxis()->GetXmax(), 0.);
  LineV2->SetLineStyle(3);
  LineH2->SetLineStyle(3);
  LineV2->Draw();
  LineH2->Draw();

  cc[ic]   = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;
  hangle->Draw("colz");
}



void Booking(UInt_t ipid)
{
  TFile *froot = new TFile("data/"+printHeader+"_ts"+pname[ipid]+".root","recreate");

  h2[0] = new TH2D("h2_0",pfname[ipid]+"; Rapidity; Pt", 100, -0.2,   1, 100,    0, 500);
  h2[1] = new TH2D("h2_1",pfname[ipid]+"; px; py"      , 100, -500, 500, 100, -500, 500);
  h2[2] = new TH2D("h2_2",pfname[ipid]+"; Rapidity; px", 100, -0.2,   1, 100, -500, 500);
  h2[3] = new TH2D("h2_3",pfname[ipid]+"; #Delta(#phi); Rapidity",100, -3.2, 3.2, nybin, -0.4, 0.45);
  h2[4] = new TH2D("h2_4",pfname[ipid]+"; #Delta(2x#phi); Rapidity",100, -3.2, 3.2, nybin, -0.4, 0.45);

  hangle    = new TH2D("hangle","; #theta; #phi", 100, 0., 1.6, 100, -3.2, 3.2);


  for(UInt_t i = 0 ; i < nybin; i++ ){
    hphi1[i] = new TH1D(Form("hphi1_%d",i),Form("hphi1_%d",i), 100, -3.2, 3.2);
    hphi2[i] = new TH1D(Form("hphi2_%d",i),Form("hphi2_%d",i), 100,  0. , 3.2);
  }

}



void amdFlow(UInt_t isel = 0, UInt_t ipid = 0)
{
  TString fdir;
  TString fname[4];

  TString ffdir = "/data/Q18393/kaneko/amd/rootfiles/132Sn124Sn270AMeV/cluster/SLy4-L108/";
  ffdir = "data/";

  switch(isel) {
  case 0:
    fdir = "/data/Q18393/kaneko/amd/rootfiles/132Sn124Sn270AMeV/cluster/SLy4-L108/";
    fname[0] = "132Sn124Sn270AMeV_cluster_SLy4-L108_bimp01.root";
    fname[1] = "132Sn124Sn270AMeV_cluster_SLy4-L108_bimp13.root";
    fname[2] = "132Sn124Sn270AMeV_cluster_SLy4-L108_bimp35.root";
    fname[3] = "132Sn124Sn270AMeV_cluster_SLy4-L108_bimp57.root";
    break;

  case 1:
    fdir = "/data/Q18393/kaneko/amd/rootfiles/132Sn124Sn270AMeV/cluster/SLy4/";
    fname[0] = "132Sn124Sn270AMeV_cluster_SLy4_bimp01.root";
    fname[1] = "132Sn124Sn270AMeV_cluster_SLy4_bimp13.root";
    fname[2] = "132Sn124Sn270AMeV_cluster_SLy4_bimp35.root";
    fname[3] = "132Sn124Sn270AMeV_cluster_SLy4_bimp57.root";
    break;

  case 2:
    fdir = "/data/Q18393/kaneko/amd/rootfiles/132Sn124Sn270AMeV/nocluster/SLy4/";
    fname[0] = "132Sn124Sn270AMeV_nocluster_SLy4_bimp01.root";
    fname[1] = "132Sn124Sn270AMeV_nocluster_SLy4_bimp13.root";
    fname[2] = "132Sn124Sn270AMeV_nocluster_SLy4_bimp35.root";
    fname[3] = "132Sn124Sn270AMeV_nocluster_SLy4_bimp57.root";
    break;

  case 3:
    fdir = "/data/Q18393/kaneko/amd/rootfiles/132Sn124Sn270AMeV/nocluster/SLy4-L108/";
    fname[0] =  "132Sn124Sn270AMeV_nocluster_SLy4-L108_bimp01.root";
    fname[1] =  "132Sn124Sn270AMeV_nocluster_SLy4-L108_bimp13.root";
    fname[2] =  "132Sn124Sn270AMeV_nocluster_SLy4-L108_bimp35.root";
    fname[3] =  "132Sn124Sn270AMeV_nocluster_SLy4-L108_bimp57.root";
    break;

  case 4:
    fdir = "/data/Q18393/kaneko/amd/rootfiles/108Sn112Sn270AMeV/nocluster/SLy4-L108/";
    fname[0] =  "108Sn112Sn270AMeV_nocluster_SLy4-L108_bimp01.root";
    fname[1] =  "108Sn112Sn270AMeV_nocluster_SLy4-L108_bimp13.root";
    fname[2] =  "108Sn112Sn270AMeV_nocluster_SLy4-L108_bimp35.root";
    fname[3] =  "108Sn112Sn270AMeV_nocluster_SLy4-L108_bimp57.root";
    break;

  case 5:
    fdir = "/data/Q18393/kaneko/amd/rootfiles/108Sn112Sn270AMeV/nocluster/SLy4/";
    fname[0] =  "108Sn112Sn270AMeV_nocluster_SLy4_bimp01.root";
    fname[1] = 	"108Sn112Sn270AMeV_nocluster_SLy4_bimp13.root";
    fname[2] = 	"108Sn112Sn270AMeV_nocluster_SLy4_bimp35.root";
    fname[3] = 	"108Sn112Sn270AMeV_nocluster_SLy4_bimp57.root";
    break;

  case 6:
    fdir = "/data/Q18393/kaneko/amd/rootfiles/108Sn112Sn270AMeV/cluster/SLy4-L108/";
    fname[0] =  "108Sn112Sn270AMeV_cluster_SLy4-L108_bimp01.root";
    fname[1] =	"108Sn112Sn270AMeV_cluster_SLy4-L108_bimp13.root";
    fname[2] =	"108Sn112Sn270AMeV_cluster_SLy4-L108_bimp35.root";
    fname[3] =	"108Sn112Sn270AMeV_cluster_SLy4-L108_bimp57.root";
    break;

  case 7:
    fdir = "/data/Q18393/kaneko/amd/rootfiles/108Sn112Sn270AMeV/cluster/SLy4/";
    fname[0] =  "108Sn112Sn270AMeV_cluster_SLy4_bimp01.root";
    fname[1] =  "108Sn112Sn270AMeV_cluster_SLy4_bimp13.root";
    fname[2] =  "108Sn112Sn270AMeV_cluster_SLy4_bimp35.root";
    fname[3] =  "108Sn112Sn270AMeV_cluster_SLy4_bimp57.root";
  }                    

  printHeader = fname[0](0,fname[0].Length()-12);

  fChain = new TChain("amdTree");
  fChain->Add(ffdir+fname[0]);
  fChain->Add(ffdir+fname[1]);
  fChain->Add(ffdir+fname[2]);
  fChain->Add(ffdir+fname[3]);


  Booking(ipid);
  Analysis(ipid);
  getv1v2(ipid);

  SaveCanvas("_"+pname[ipid]);
}




