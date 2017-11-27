TCanvas *cc[12];
TChain  *rChain[2];

const Int_t nbinx = 30;
Double_t dxbin = 1./(Double_t)nbinx; 

STPID::PID selectPID = STPID::kNon; // kPion, kProton, kDeuteron, kTriton, k3He, k4He                                                     

TString printName;
TF1 *f1;

UInt_t m_bgn = 0;
UInt_t m_end = 1;


Int_t ichain = 0;
TCut proton  = "fPID==2212";
TCut deuteron= "fPID==1000010020";
TCut triton  = "fPID==1000010030";

void drawRComp()
{
  gROOT->Reset();
  
  gROOT->Macro("openRComp.C");

  rChain[ichain] = (TChain*)gROOT->FindObject(Form("rChain%d",ichain));
  if(rChain[ichain] != NULL) ichain++;

  rChain[ichain] = (TChain*)gROOT->FindObject(Form("rChain%d",ichain));
  if(rChain[ichain] == NULL) ichain = 0;
  
  if(rChain[0] == NULL && rChain[1] == NULL)
    exit(0);
  
  cout << " ichain " << ichain <<endl;

  m_end = ichain+1;

}


void NuSYM()
{
  gROOT->cd();
  gROOT->Clear();



  Int_t ic = 0;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  auto hRPphi    = new TH1D("hRPphi","; R. P.#Phi",120,-3.2,3.2);
  rChain[0] -> Draw("TVector2::Phi_mpi_pi(unitP_lang.Phi())>>hRPphi","ntrack[4]>2");
  hRPphi->SetLineColor(4);


  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  auto hsubphi  = new TH2D("hsubphi","SubEvnt Corr",120,-3.2,3.2,120,-3.2,3.2);
  rChain[0] -> Draw("TVector2::Phi_mpi_pi(unitP_1.Phi()-unitP_2.Phi()):TVector2::Phi_mpi_pi(unitP_lang.Phi())>>hsubphi",
		    "mtrack_1>0&&mtrack_2>0","colz");


  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  auto hsubcorr0 = new TH1D("hsubcorr0","#Phi_A - #Phi_B", 120,-3.2,3.2);
  rChain[0] -> Draw("TVector2::Phi_mpi_pi(unitP_1.Phi()-unitP_2.Phi())>>hsubcorr0","mtrack_1>0&&mtrack_2>0");
  


}


void deltphi()
{
  const UInt_t nbin = 3;


  TH1D *hyphi[nbin][2];

  Double_t ybin = (0.65 - 0.15)/(Double_t)nbin;

  Double_t ymin[3] = {0.15, 0.32, 0.48};
  Double_t ymax[3] = {0.32, 0.48, 0.65};

  gStyle -> SetOptStat(0);
  gStyle -> SetPadBottomMargin(0.10);
  gStyle -> SetPadLeftMargin(0.11);
  gStyle -> SetPadRightMargin(0.12);
  gStyle -> SetTitleFontSize(0.06);

  UInt_t ic = 0;
  
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),500,800);
  cc[ic]->Divide(1,3);


  for(UInt_t m = m_bgn; m < m_end; m++){
    hyphi[0][m] = new TH1D(Form("hyphi0_%d",m),"0.15 < y < 0.32 ;#Delta(#phi - #PhiR.P.); counts",20,0.,3.1);
    hyphi[1][m] = new TH1D(Form("hyphi1_%d",m),"0.32 < y < 0.48 ;#Delta(#phi - #PhiR.P.); counts",20,0.,3.1);
    hyphi[2][m] = new TH1D(Form("hyphi2_%d",m),"0.48 < y < 0.65 ;#Delta(#phi - #PhiR.P.); counts",20,0.,3.1);
    
    for(UInt_t i = 0; i < nbin; i++){

      cc[ic]->cd(i+1); 

      hyphi[i][m]->GetXaxis()->SetLabelSize(0.05);
      hyphi[i][m]->GetXaxis()->SetTitleSize(0.05);
      hyphi[i][m]->GetYaxis()->SetTitleSize(0.05);

      hyphi[i][m]->SetMarkerStyle(20);
      hyphi[i][m]->SetMarkerColor(2+m);
      hyphi[i][m]->SetLineColor(2+m);


      TString hname = Form("hyphi%d_%d",i,m);
      
      TString opt;
      if(m == 0) opt = "e";
      else       opt = "same e";
	
      rChain[m] -> Draw("abs(fdeltphi)>>"+hname,Form("mtrack>2&&fPID==2212&&fRapidity>%f&&fRapidity<%f",ymin[i],ymax[i]),opt);


    }
  }

  cc[ic]->cd(1);
  auto aLeg = new TLegend(0.15,0.7,0.3,0.9,"");
  aLeg->AddEntry(hyphi[0][0],"REAL ","lp");
  aLeg->AddEntry(hyphi[0][1],"MIXed","lp");

  aLeg->Draw();

}

void phi()
{
  gROOT->cd();

  //plot <phi> pseudorapidity 
  TH2D *hphiy[2];
  TH2D *hprtphiy[2];

  TH1D *hiphi[2];
  TH1D *huphi[2];
  TH1D *hsphi1[2];
  TH1D *hsphi2[2];

  TH2D *hphitheta[2];
  TH2D *hsubphi[2];
  TH2D *hphimtk[2];

  TH2D *hphith[2];
  TH2D *hprtphith[2];

  TH2D *hdtrphiy[2];
  TH2D *htrtphiy[2];


  Bool_t idv = kTRUE;
  TString opt = "";
  for(UInt_t m = m_bgn; m < m_end; m++){

    hphitheta[m]= new TH2D(Form("hphitheta%d",m),"All; #phi ; #theta",100,0.,1.5,100,-3.2,3.2);
    hiphi[m]    = new TH1D(Form("hiphi%d",m),"All; #varphi",100,-3.2,3.2);
    huphi[m]    = new TH1D(Form("huphi%d",m),"R. P.; #Phi",100,-3.2,3.2);
    hsphi1[m]   = new TH1D(Form("hsphi1%d",m),"SubEvent 1; #Phi",100,-3.2,3.2);
    hsphi2[m]   = new TH1D(Form("hsphi2%d",m),"SubEvent 2; #Phi",100,-3.2,3.2);
    hsubphi[m]  = new TH2D(Form("hsubphi%d",m),"SubEvnt Corr",100,-3.2,3.2,100,-3.2,3.2);
    //    hphiy[m]    = new TH2D(Form("hphiy%d",m),"All ; pseudorapidity; #Delta#varphi_i",nbinx,0,3,60,-3.2,3.2);
    hphiy[m]    = new TH2D(Form("hphiy%d",m),   " All    ; rapidity; #Delta#varphi_i",nbinx,-1.,1.,60,-3.2,3.2);
    hprtphiy[m] = new TH2D(Form("hprtphiy%d",m),"Proton  ; rapidity; #Delta#varphi_i",nbinx,-1.,1.,60,-3.2,3.2);
    hdtrphiy[m] = new TH2D(Form("hdtrphiy%d",m),"Deutron ; rapidity; #Delta#varphi_i",nbinx,-1.,1.,60,-3.2,3.2);
    htrtphiy[m] = new TH2D(Form("htrtphiy%d",m),"Triton  ; rapidity; #Delta#varphi_i",nbinx,-1.,1.,60,-3.2,3.2);

    hphimtk[m]  = new TH2D(Form("hphimtk%d",m),"All phi vs ntrack ",60,0,60,100,-3.2,3.2);
    hphith[m]   = new TH2D(Form("hphith%d",m),"All ; #theta; #Delta#verphi",100,0.,1.6,100,-3.2,3.2);
    hprtphith[m]= new TH2D(Form("hprtphith%d",m),"Proton ; #theta; #Delta#verphi",100,0.,1.6,100,-3.2,3.2);


    UInt_t ic = 0;
    UInt_t col = 4;
    UInt_t id = 1;
    
    if( m == 1) {
      ic  = 0;
      opt = "same";
      col = 2;
      idv  = kFALSE;
    }
    else { 
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
      cc[ic]->Divide(2,2);
    }

    TCut gtrack =  "ntrack[4]>2&&fReactionPlanef>10";

    id = 1;
    cc[ic]->cd(id); id++;
    TString hname = Form("hiphi%d",m);
    rChain[m] -> Draw("frphi>>"+hname,gtrack,opt);
    hiphi[m]->SetLineColor(col);

    cc[ic]->cd(id); id++;
    hname = Form("huphi%d",m);
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(unitP_lang.Phi())>>"+hname,"ntrack[4]>2",opt);
    huphi[m]->SetLineColor(col);


    cc[ic]->cd(id); id++;
    hname = Form("hsphi1%d",m);
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(unitP_1.Phi())>>"+hname,"mtrack_1>1",opt);
    hsphi1[m]->SetLineColor(col);

    cc[ic]->cd(id); id++;
    hname = Form("hsphi2%d",m);
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(unitP_2.Phi())>>"+hname,"mtrack_2>1",opt);
    hsphi2[m]->SetLineColor(col);

    cout << " m " << m << " ic " << ic << endl;


    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),800,1000);
    cc[ic]->Divide(2,2);


    id = 1;
    cc[ic]->cd(id); id++;

    auto y_cm    = 0.388568;
    auto ybm_cm  = 0.36599;
    TString rpd = Form("(fRapidity - %f)/%f >>",y_cm,ybm_cm);

    hname = Form("hphiy%d",m);
    //    rChain[m] -> Draw("TVector2::Phi_mpi_pi(frphi-frpphi):fpsudoRapidity>>"+hname,gtrack,"colz");
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(frphi-frpphi):"+rpd+hname,gtrack,"colz");

    cc[ic]->cd(id); id++;
    hname = Form("hprtphiy%d",m);
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(frphi-frpphi):"+rpd+hname,gtrack*"fPID==2212","colz");

    cc[ic]->cd(id); id++;
    hname = Form("hdtrphiy%d",m);
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(frphi-frpphi):"+rpd+hname,gtrack*"fPID==1000010020","colz");

    cc[ic]->cd(id); id++;
    hname = Form("htrtphiy%d",m);
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(frphi-frpphi):"+rpd+hname,gtrack*"fPID==1000010030","colz");



    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    hname = Form("hsubphi%d",m);
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(unitP_1.Phi()-unitP_2.Phi()):TVector2::Phi_mpi_pi(unitP_lang.Phi())>>"
		      +hname, "mtrack_1>0&&mtrack_2>0","colz");

    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->Divide(1,2);

    cc[ic]->cd(id); id++;
    hname = Form("hphith%d",m);
    rChain[m]->Draw("fdeltphi:fRotatedP3.Theta()>>"+hname,gtrack,"colz");

    cc[ic]->cd(id); id++;
    hname = Form("hprtphith%d",m);
    rChain[m]->Draw("fdeltphi:fRotatedP3.Theta()>>"+hname,gtrack*"fPID==2212","colz");


  }


}



void PtDistribution()
{
  gROOT->Reset();


  Double_t y_min = 0;
  Double_t y_bin = 0.1;
  const UInt_t  nbin = 16;

  Double_t pt_pmax =  800.;
  Double_t pt_dmax = 1000.;
  Double_t pt_tmax = 1400.;
  Double_t pt_pimax = 300.;  
  UInt_t   pt_nbin  = 100;


  TH1D* hptpm[2][nbin];
  TH1D* hptpp[2][nbin];
  TH1D* hptprt[2][nbin];
  TH1D* hptdtr[2][nbin];
  TH1D* hpttrt[2][nbin];

  UInt_t icol[2]={4,2};
  TString iopt[2]={"","same"};

 

  UInt_t ic = -1;   

  for(UInt_t m = m_bgn; m < m_end; m++){
 
    auto nEntry = (Double_t)rChain[m]->GetEntries();

    for(UInt_t i = 0; i < nbin; i++){
    
      Double_t yL = y_min + y_bin*i;
      Double_t yU = yL + y_bin;

      TString hname = Form("hptprt%d%d",m,i);
      TString htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptprt[m][i] = new TH1D(hname,"Proton : "+htitle,pt_nbin, 0., pt_pmax);
      rChain[m]->Project(hname,"fRotatedP3.Pt()","fPID==2212&&"+htitle);

      hname = Form("hptdtr%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptdtr[m][i] = new TH1D(hname,"Deuteron : "+htitle,pt_nbin, 0., pt_dmax);
      rChain[m]->Project(hname,"fRotatedP3.Pt()","fPID==1000010020&&"+htitle);

      hname = Form("hpttrt%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hpttrt[m][i] = new TH1D(hname,"Triton : "+htitle,pt_nbin, 0., pt_tmax);
      rChain[m]->Project(hname,"fRotatedP3.Pt()","fPID==1000010030&&"+htitle);

      hname = Form("hptpm%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpm[m][i] = new TH1D(hname,"Pi- : "+htitle,pt_nbin, 0., pt_pimax);
      rChain[m]->Project(hname,"fRotatedP3.Pt()","fPID==211&&fChar<0&&"+htitle);

      hname = Form("hptpp%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpp[m][i] = new TH1D(hname,"Pi+ : "+htitle,pt_nbin, 0., pt_pimax);
      rChain[m]->Project(hname,"fRotatedP3.Pt()","fPID==211&&fChar>0&&"+htitle);

    }
  

    auto sclp  = 1./nEntry * (Double_t)pt_nbin / pt_pmax  * y_bin;
    auto scld  = 1./nEntry * (Double_t)pt_nbin / pt_dmax  * y_bin;
    auto sclt  = 1./nEntry * (Double_t)pt_nbin / pt_tmax  * y_bin;
    auto sclpi = 1./nEntry * (Double_t)pt_nbin / pt_pimax * y_bin;


    for(UInt_t j = 0; j < nbin; j++){
      hptprt[m][j]->Scale(sclp);
      hptprt[m][j]->SetLineColor(icol[m]);
    }

    ic = 0;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else {
      cc[ic]->cd();
 
      for(UInt_t j = 0; j < nbin; j++){
	cc[ic]->cd(j+1);
	hptprt[1][j]->Draw(iopt[0]);
	hptprt[0][j]->Draw(iopt[1]);

      }
    }


    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else
      cc[ic]->cd();

    for(UInt_t j = 0; j < nbin; j++){
      cc[ic]->cd(j+1);
      hptdtr[m][j]->Scale(scld);
      hptdtr[m][j]->SetLineColor(icol[m]);
      hptdtr[m][j]->Draw(iopt[m]);
    }


    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else
      cc[ic]->cd();

    for(UInt_t j = 0; j < nbin; j++){
      cc[ic]->cd(j+1);
      hpttrt[m][j]->Scale(sclt);
      hpttrt[m][j]->SetLineColor(icol[m]);
      hpttrt[m][j]->Draw(iopt[m]);
    }

    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else
      cc[ic]->cd();

    for(UInt_t j = 0; j < nbin; j++){
      cc[ic]->cd(j+1);
      hptpm[m][j]->Scale(sclpi);
      hptpm[m][j]->SetLineColor(icol[m]);
      hptpm[m][j]->Draw(iopt[m]);
    }

    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,1500);
      cc[ic]->Divide(2,nbin/2);
    }
    else
      cc[ic]->cd();

    for(UInt_t j = 0; j < nbin; j++){
      cc[ic]->cd(j+1);
      hptpp[m][j]->Scale(sclpi);
      hptpp[m][j]->SetLineColor(icol[m]);
      hptpp[m][j]->Draw(iopt[m]);
    }
  }

}


void RP()
{
  gROOT->Reset();

  TH1D *h1[2];
  TH1D *h2[2];
  TH1D *h3[2];
  TH2D *h4[2];

  TH1D *h5[2];
  TH1D *h6[2];
  TH2D *h7[2];
  TH1D *h8[2];

  TH2D *h9[2];  
  TH2D *h10[2];

  UInt_t ic = 0;
  UInt_t id = 1;

  TString hname;

  TCut gevt = "mtrack>2&&mtrack_1>1&&mtrack_2>1";

  for(UInt_t m = m_bgn; m < m_end; m++){


    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->Divide(1,2);    
    id = 1;

    cc[ic]->cd(id); id++;
    hname = Form("h9%d",m);
    h9[m] = new TH2D(hname,"RP corrected; mtrack; #Phi",100,0,100,100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_lang.Phi()):mtrack>>"+hname,gevt,"colz");

    cc[ic]->cd(id); id++;
    hname = Form("h10%d",m);
    h10[m] = new TH2D(hname,"RP Not corrected; mtrack; #Phi",100,0,100,100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_rot.Phi()):mtrack>>"+hname,gevt,"colz");

    ic++; 
    continue;



    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->Divide(2,2);

    cc[ic]->cd(id); id++;
    hname = Form("h2%d",m);
    h2[m] = new TH1D(hname,"RP Sub A Not corrected; #Phi",100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_1r.Phi())>>"+hname,gevt,"colz");
    
    cc[ic]->cd(id); id++;
    hname = Form("h3%d",m);
    h3[m] = new TH1D(hname,"RP Sub B Not corrected; #Phi",100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_2r.Phi())>>"+hname,gevt,"colz");

    cc[ic]->cd(id); id++;
    hname = Form("h4%d",m);
    h4[m] = new TH2D(hname,"RP Sub B Not corrected; #Phi_A - #Phi_B; #Phi",100,-3.2,-3.2,100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_1r.Phi()-unitP_2r.Phi()):TVector2::Phi_mpi_pi(unitP_rot.Phi())>>"+hname,gevt,"colz");
    

    cc[ic]->cd(id); id++;
    hname = Form("h1%d",m);
    h1[m] = new TH1D(hname,"RP Not corrected; #Phi",100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_rot.Phi())>>"+hname,gevt,"colz");


    ic++; id=1;

    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->Divide(2,2);


    cc[ic]->cd(id); id++;
    hname = Form("h5%d",m);
    h5[m] = new TH1D(hname,"RP Sub A corrected; #Phi",100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_1.Phi())>>"+hname,gevt,"colz");
    
    cc[ic]->cd(id); id++;
    hname = Form("h6%d",m);
    h6[m] = new TH1D(hname,"RP Sub B corrected; #Phi",100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_2.Phi())>>"+hname,gevt,"colz");

    cc[ic]->cd(id); id++;
    hname = Form("h7%d",m);
    h7[m] = new TH2D(hname,"RP corrected; #Phi_A - #Phi_B; #Phi",100,-3.2,-3.2,100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_1.Phi()-unitP_2.Phi()):TVector2::Phi_mpi_pi(unitP_lang.Phi())>>"+hname,gevt,"colz");
    

    cc[ic]->cd(id); id++;
    hname = Form("h8%d",m);
    h8[m] = new TH1D(hname,"RP corrected; #Phi",100,-3.2,-3.2);
    rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_lang.Phi())>>"+hname,gevt,"colz");


  }
}


void phi_mtk()
{
  gROOT->Reset();
  
  TH1D *hPHI[2][10];
  TH1D *hphia[2][10];
  TH1D *hphib[2][10];
  TH1I *hmtk[2];

  TString opt = "";
  for(UInt_t m = m_bgn; m < m_end; m++){


    Bool_t idv = kFALSE;
    UInt_t ic = m;

    if( m == 1) {
      ic = 0;
      opt = "same";
    }
    else {
      idv = kTRUE;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,1000);
      cc[ic]->Divide(3,6);
    }

    UInt_t col = 4;
    if( m == 1)
      col = 2;


    Int_t mtrk_binw = 10;
    UInt_t id = 1;

    for(UInt_t i = 0; i < 6; i++){
      Int_t mtrk_l = mtrk_binw * i;
      Int_t mtrk_u = mtrk_l + mtrk_binw;

      TCut gtrack = "mtrack>2&&fReactionPlanef>=10";       
      TCut mtcut  = Form("mtrack>%d&&mtrack<=%d",mtrk_l,mtrk_u);
      TCut hcut   = gtrack&&mtcut;

      TString hname = Form("hPHI%d%d",m,i);
      hPHI[m][i] = new TH1D(hname,hcut,100,-3.2,3.2);
      hPHI[m][i]->SetLineColor(col);
      cc[ic]->cd(id); id++;
      rChain[m]->Draw("TVector2::Phi_mpi_pi(unitP_lang.Phi())>>"+hname,hcut,opt);


      hname = Form("hphia%d%d",m,i);
      hphia[m][i] = new TH1D(hname,hcut,100,-3.2,3.2);
      cc[ic]->cd(id); id++;
      rChain[m]->Draw("frphi>>"+hname,hcut,opt);
      Double_t val = hphia[m][i]->GetMaximum();
      //      hphia[m][i]->SetMinimum(val*0.9);
      //      hphia[m][i]->SetMaximum(val*1.5);
      hphia[m][i]->Draw(opt+" e");
      hphia[m][i]->SetLineColor(col);

      hname = Form("hphib%d%d",m,i);
      hphib[m][i] = new TH1D(hname,hcut,100,-3.2,3.2);
      hphib[m][i]->SetLineColor(col);
      cc[ic]->cd(id); id++;
      rChain[m]->Draw("fphi>>"+hname,hcut,opt);


    }

    ic++;
    if( idv ) {
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
      //      cc[ic]->Divide(2,2);
    }

    id = 1;


    cc[ic]->cd();
    TString hname = Form("hmtk%d",m);

    hmtk[m] = new TH1I(hname,"multiplicity ",60,0,60);
    rChain[m]->Draw("mtrack>>"+hname,"",opt);
    hmtk[m]->SetLineColor(col);
    hmtk[m]->Draw(opt);

    cout << " m " << m << " ic " << ic <<  " id " << id << " opt " << opt << endl;
    



  }
}


void RPphi()
{
  gROOT->cd();
  //  TCut myselect = "mtrack>25";
  //  TCut myselect = "mtrack<25&&mtrack>15";
  //  TCut myselect = "mtrack<=15";
  TCut myselect = "mtrack<=10";

  //plot <phi> pseudorapidity 
  TH1D *huphi[2];
  TH1D *hsubphi[2];
  TLegend *aLeg[2];

  UInt_t ic = 0;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);


  for(UInt_t m = m_bgn; m < m_end; m++){

    huphi[m]    = new TH1D(Form("huphi%d",m),    "R.P.; #Phi"                 ,100,-3.2,3.2);
    hsubphi[m]  = new TH1D(Form("hsubphi%d",m),"Sub Events Correlation; #Phi_A - #Phi_B; counts",100,-3.1,3.1);


    UInt_t col = 2+m;
    UInt_t id = 1;

    TString opt = "e";    
    if(m == 1) 
      opt = "same e";

    
    ic = 0;
    cc[ic]->cd();
    huphi[m]->SetLineColor(col);
    //    huphi[m]->SetMaximum(16000);
    huphi[m]->SetMinimum(0);
    TString hname = Form("huphi%d",m);
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(unitP_lang.Phi())>>"+hname,myselect,opt);

    


    ic++;
    cc[ic]->cd();
    hname = Form("hsubphi%d",m);
    rChain[m] -> Draw("TVector2::Phi_mpi_pi(unitP_1.Phi()-unitP_2.Phi())>>"+hname,myselect&&"mtrack_1>0&&mtrack_2>0",opt);
    hsubphi[m]->SetMarkerStyle(20);
    hsubphi[m]->SetMarkerColor(col);
    hsubphi[m]->SetLineColor(col);

  }

  

  ic = 0;
  cc[ic]->cd(); ic++;
  aLeg[0] = new TLegend(0.15,0.15,0.3,0.3,"");
  aLeg[0]->AddEntry(huphi[0],"REAL ","lp");
  aLeg[0]->AddEntry(huphi[1],"MIXed ","lp");
  aLeg[0]->Draw();

  cc[ic]->cd(); 
  aLeg[1] = new TLegend(0.15,0.75,0.3,0.9,"");
  aLeg[1]->AddEntry(hsubphi[0],"REAL ","lp");
  aLeg[1]->AddEntry(hsubphi[1],"MIXed ","lp");
  aLeg[1]->Draw();


}

//*************************
void PID()
{
  gROOT->cd();

  bool drawPIDLine = true;
  Int_t nbins = 200;
  Int_t p1 = -500;
  Int_t p2 = 2500;
  Int_t dedx1 = 0;
  Int_t dedx2 = 800;

  TString tLabel[2]={"108Sn","132Sm"};

  TH2D *histPID[2];
  TH2D *histacc[2];
  TH1D *histmom[2];
  TH1D *histrapd[2];
  TH1D *histdedx[2];

  Int_t ic = -1;
  for(UInt_t m = m_bgn; m < m_end; m++){
    UInt_t id = 1;
    TString hname = Form("histPID%d",m);
    histPID[m] = new TH2D(hname,";p/Z (MeV/c); dEdx (ADC/mm)",nbins,p1,p2,nbins,dedx1,dedx2);
    rChain[m]->Project(hname,"fdEdx:fP","fPID==2212");

    hname = Form("histmom%d",m);
    histmom[m] = new TH1D(hname,";p (MeV/c)",nbins,p1,p2);
    rChain[m]->Project(hname,"fP","fPID==2212");

    hname = Form("histrapd%d",m);
    histrapd[m] = new TH1D(hname,";Rapidity",nbins,0.,1.);
    rChain[m]->Project(hname,"fRapidity","fPID==2212");

    hname = Form("histdedx%d",m);
    histdedx[m] = new TH1D(hname,";dEdx (ADC/mm)",nbins,0.,100.);
    rChain[m]->Project(hname,"fdEdx","fPID==2212");

    hname = Form("histacc%d",m);
    histacc[m] = new TH2D(hname,"Acceptance; Rapidity; pt(MeV/C)",200,0,1.,200,0,1000);
    rChain[m]->Project(hname,"fRotatedP3.Pt():fRapidity","fPID==2212");


    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->Divide(2,2);
       
    cc[ic]->cd(id); id++;
    histPID[m]->Draw("colz");

    cc[ic]->cd(id); id++;
    histacc[m]->Draw("colz");

    cc[ic]->cd(id); id++;
    histdedx[m]->Draw();

    cc[ic]->cd(id); id++;
    histrapd[m]->Draw();

  

    if( m == 1 ){
      TH2D *histaccDiv = new TH2D( (*histacc[0])/(*histacc[1]) );
  
      ic++;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
      histaccDiv->Draw("colz");
      cc[ic]->SetLogz();
    }
  }
}


void CalibCheck()
{
  gROOT->cd();


  TH2D* histc[4];

  histc[0] = new TH2D("histc0","108Sn Proton dEdx vs RUN",500,2260,2510,200,10,400);
  histc[1] = new TH2D("histc1","132Sn Proton dEdx vs RUN",500,2840,3040,200,10,400);
  histc[2] = new TH2D("histc2","108Sn Proton dEdx vs RUN",20,2260,2280,200,10,400);



  cc[0] = new TCanvas("cc0","cc0",700,500);

  TCut momcut = "fP>800&&fP<1000";
  rChain[0]->Draw("fdEdx:irun>>histc0",momcut,"colz");


  cc[1] = new TCanvas("cc1","cc1",700,500);
  rChain[1]->Draw("fdEdx:irun>>histc1",momcut,"colz");
  

  cc[2] = new TCanvas("cc2","cc2",700,500);
  rChain[0]->Draw("fdEdx:irun>>histc2",momcut,"colz");

  
  
  
  

}


void acceptance()
{
  gROOT->Reset();
  
  TH2D *haccprt[2];
  TH2D *haccdtr[2];
  TH2D *hacctrt[2];


  UInt_t ic = -1;
  for(UInt_t m = m_bgn; m < m_end; m++){

    ic++;
    UInt_t id = 1;

    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->Divide(2,2);

    TString hname = Form("haccptr%d",m);
    haccprt[m] = new TH2D(hname,"Proton acc",200,0,1.2,200,0,1000);
    
    cc[ic]->cd(id); id++;
    rChain[m]->Draw("fRotatedP3.Pt():fRapidity>>"+hname,"fPID==2212","colz");

    hname = Form("haccdtr%d",m);
    haccdtr[m] = new TH2D(hname,"Deuteron acc",200,0,1.,200,0,1000);
    
    
    cc[ic]->cd(id); id++;
    rChain[m]->Draw("fRotatedP3.Pt():fRapidity>>"+hname,"fPID==1000010020","colz");

    hname = Form("hacctrt%d",m);
    hacctrt[m] = new TH2D(hname,"Trition acc",200,0,0.8,200,0,1000);
    
    cc[ic]->cd(id); id++;
    rChain[m]->Draw("fRotatedP3.Pt():fRapidity>>"+hname,"fPID==1000010030","colz");
    

  }
}

void dndy()
{
  gROOT->cd();
  UInt_t ic = -1;

  bool  drawPIDLine = true;
  Int_t nbins = 100;
  auto y_min  = 0.;
  auto y_max  = 1.2;
  auto y_tmin  = 0.;
  auto y_tmax  = 0.8;
  auto y_pimin  = 0.2;
  auto y_pimax  = 1.8;

  auto y_cm    = 0.388568;
  auto ybm_cm  = 0.36599;
  //  TString rpd = Form("(fRapidity - %f)/%f",y_cm,ybm_cm);
  TString rpd = "fRapidity";


  TH1D *histrapdpm[2];
  TH1D *histrapdpp[2];
  TH1D *histrapdp[2];
  TH1D *histrapdd[2];
  TH1D *histrapdt[2];

  TH2D *histaccpm[2];
  TH2D *histaccpp[2];
  TH2D *histaccp[2];
  TH2D *histaccd[2];
  TH2D *histacct[2];

  Double_t pmnum[2];
  Double_t ppnum[2];
  Double_t pnum[2];
  Double_t dnum[2];
  Double_t tnum[2];

  UInt_t icol[2]={4,2};
  TString iopt[2]={"","same"};

  for(UInt_t m = m_bgn; m < m_end; m++){

    auto nEntry = (Double_t)rChain[m]->GetEntries();


    UInt_t id = 1;
    TString hname = Form("histrapdp%d",m);
    histrapdp[m] = new TH1D(hname,"Proton;ylab; dNdy",nbins, y_min, y_max);
    rChain[m]->Project(hname,rpd,"fPID==2212&&ftheta<0.8");

    hname = Form("histrapdd%d",m);
    histrapdd[m] = new TH1D(hname,"Deuteron;ylab; dNdy",nbins, y_min, y_max);
    rChain[m]->Project(hname,rpd,"fPID==1000010020&&ftheta<0.8");

    hname = Form("histrapdt%d",m);
    histrapdt[m] = new TH1D(hname,"Trition;ylab; dNdy",nbins, y_tmin, y_tmax);
    rChain[m]->Project(hname,rpd,"fPID==1000010030&&ftheta<0.8");

    hname = Form("histrapdpm%d",m);
    histrapdpm[m] = new TH1D(hname,"Pi-;ylab; dNdy",nbins, y_pimin, y_pimax);
    rChain[m]->Project(hname,rpd,"fPID==211&&fChar<0&&ftheta<0.8");

    hname = Form("histrapdpp%d",m);
    histrapdpp[m] = new TH1D(hname,"Pi+;ylab; dNdy",nbins, y_pimin, y_pimax);
    rChain[m]->Project(hname,rpd,"fPID==211&&fChar>0&&ftheta<0.8");


    hname = Form("histaccpm%d",m);
    histaccpm[m] = new TH2D(hname,"Pi-; Rapidity; pt(MeV/C)",200,  y_pimin, y_pimax, 200,0,300);
    rChain[m]->Project(hname,"fRotatedP3.Pt():fRapidity","fPID==211&&fChar<0&&ftheta<0.8");

    hname = Form("histaccpp%d",m);
    histaccpp[m] = new TH2D(hname,"Pi+; Rapidity; pt(MeV/C)",200, y_pimin, y_pimax, 200,0,300);
    rChain[m]->Project(hname,"fRotatedP3.Pt():fRapidity","fPID==211&&fChar>0&&ftheta<0.8");

    hname = Form("histaccp%d",m);
    histaccp[m] = new TH2D(hname,"Proton; Rapidity; pt(MeV/C)",200,y_min, y_max,200,0,1000);
    rChain[m]->Project(hname,"fRotatedP3.Pt():fRapidity","fPID==2212&&ftheta<0.8");

    hname = Form("histaccd%d",m);
    histaccd[m] = new TH2D(hname,"Deuteron; Rapidity; pt(MeV/C)",200,y_min, y_max,200,0,1000);
    rChain[m]->Project(hname,"fRotatedP3.Pt():fRapidity","fPID==1000010020&&ftheta<0.8");

    hname = Form("histacct%d",m);
    histacct[m] = new TH2D(hname,"Triton; Rapidity; pt(MeV/C)",200,y_tmin, y_tmax,200,0,1000);
    rChain[m]->Project(hname,"fRotatedP3.Pt():fRapidity","fPID==1000010030&&ftheta<0.8");


    auto scl  = 1./nEntry * (Double_t)nbins / (y_max - y_min) ;
    auto sclt = 1./nEntry * (Double_t)nbins / (y_tmax - y_tmin) ;
    auto sclpi= 1./nEntry * (Double_t)nbins / (y_pimax - y_pimin) ;

    ic = 0;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    }      
    else
      cc[ic]->cd();
    histrapdpp[m]->Scale(sclpi);
    histrapdpp[m]->SetLineColor(icol[m]);
    histrapdpp[m]->Draw(iopt[m]);
    ppnum[m] = histrapdpp[m]->Integral()/nbins;

    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    }      
    else
      cc[ic]->cd();
    histrapdpm[m]->Scale(sclpi);
    histrapdpm[m]->SetLineColor(icol[m]);
    histrapdpm[m]->Draw(iopt[m]);
    pmnum[m] = histrapdpm[m]->Integral()/nbins;


    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
      histrapdp[0]->SetMaximum(10.);
    }      
    else
      cc[ic]->cd();
    histrapdp[m]->Scale(scl);
    histrapdp[m]->SetLineColor(icol[m]);
    histrapdp[m]->Draw(iopt[m]);
    pnum[m] = histrapdp[m]->Integral()/nbins;


    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    }      
    else
      cc[ic]->cd();
    histrapdd[m]->Scale(scl);
    histrapdd[m]->SetLineColor(icol[m]);
    histrapdd[m]->Draw(iopt[m]);
    dnum[m] = histrapdd[m]->Integral()/nbins;

    ic++;
    if(m == 0){
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    }      
    else
      cc[ic]->cd();
    histrapdt[m]->Scale(sclt);
    histrapdt[m]->SetLineColor(icol[m]);
    histrapdt[m]->Draw(iopt[m]);
    tnum[m] = histrapdt[m]->Integral()/nbins;

    // if(m == 0)
    //   ic++;
    // else
    //   ic = 
    // cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    // histaccp[m]->Draw("colz");

    // ic++;
    // cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    // histaccd[m]->Draw("colz");
    
    // ic++;
    // cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    // histacct[m]->Draw("colz");

    // ic++;
    // cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    // histaccpm[m]->Draw("colz");

    // ic++;
    // cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    // histaccpp[m]->Draw("colz");


    if( m == 1 ){
      ic = 5;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
      TH1D *histrapdpmDiv = new TH1D( (*histrapdpm[0])/(*histrapdpm[1]) );
      histrapdpmDiv->SetTitle("Pi- ; ylab;  dNdy(132Sn/ 108Sn)");
      histrapdpmDiv->Draw();

      ic++;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
      TH1D *histrapdppDiv = new TH1D( (*histrapdpp[0])/(*histrapdpp[1]) );
      histrapdppDiv->SetTitle("Pi+ ; ylab;  dNdy(132Sn/ 108Sn)");
      histrapdppDiv->Draw();

      ic++;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
      TH1D *histrapdpDiv = new TH1D( (*histrapdp[0])/(*histrapdp[1]) );
      histrapdpDiv->SetTitle("Proton ; ylab;  dNdy(132Sn/ 108Sn)");
      histrapdpDiv->Draw();

      ic++;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
      TH1D *histrapddDiv = new TH1D( (*histrapdd[0])/(*histrapdd[1]) );
      histrapddDiv->SetTitle("Deuteron; ylab;  dNdy(132Sn/ 108Sn)");
      histrapddDiv->Draw();

      ic++;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
      TH1D *histrapdtDiv = new TH1D( (*histrapdt[0])/(*histrapdt[1]) );
      histrapdtDiv->SetTitle("Triton ; ylab;  dNdy(132Sn/ 108Sn)");
      histrapdtDiv->Draw();

    }

    cout << m << endl;
    // cout << "number of p" << " >> " << pnum[m] << " : " << dnum[m] << " : " << tnum[m] << " : " << pnum[m]+dnum[m]+tnum[m] << endl;
    // cout << "number of n" << " >> " << "        : " << dnum[m] << " : " << tnum[m]*2<<" : " << dnum[m]+tnum[m]*2 << endl;
    // cout << " Ratio (n/p) " << (dnum[m]+tnum[m]*2)/(pnum[m]+dnum[m]+tnum[m]) << endl;

    cout << "Ratio (pi- / pi+)" << pmnum[m] << " : " << ppnum[m] << endl;
  }  

  cout << " 108/132 : p " << (pnum[0]+dnum[0]+tnum[0])/(pnum[1]+dnum[1]+tnum[1]) << " n " << (dnum[0]+tnum[0]*2)/(dnum[1]+tnum[1]*2) << endl; 

}
