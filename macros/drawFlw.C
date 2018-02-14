#include "openFlw.C"

TCanvas *cc[12];
TChain  *rChain[2];
Int_t ic = -1;


const Int_t nbinx = 30;
Double_t dxbin = 1./(Double_t)nbinx; 

STPID::PID selectPID = STPID::kNon; // kPion, kProton, kDeuteron, kTriton, k3He, k4He                                                     

TString printName;

UInt_t m_bgn = 0;
UInt_t m_end = 1;


Int_t ichain = 0;
TCut proton  = "fPID==2212";
TCut deuteron= "fPID==1000010020";
TCut triton  = "fPID==1000010030";

void LoadNLgCut();

void drawFlw()
{
  gROOT->Reset();
  
  openFlw();

  rChain[ichain] = (TChain*)gROOT->FindObject(Form("rChain%d",ichain));
  if(rChain[ichain] != NULL) ichain++;

  rChain[ichain] = (TChain*)gROOT->FindObject(Form("rChain%d",ichain));
  if(rChain[ichain] == NULL) ichain = 0;
  
  if(rChain[0] == NULL && rChain[1] == NULL)
    exit(0);
  
  cout << " ichain " << ichain <<endl;

  m_end = ichain+1;

  gROOT->ProcessLine(".! grep -i void drawFlw.C |grep //%%");
}




void deltphi()   //%%  
{
  const UInt_t nbin = 3;


  TH1D *hyphi[nbin][4];

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

void phi()       //%%
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



void PtDistribution() //%
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


void RP()        //%%
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


void phi_mtk()   //%%
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


void RPphi()     //%%
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
void PID()       //%%
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


void CalibCheck()//%%
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


void acceptance()//%%
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


void dndy()      //%%
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


                  
//________________________________//%% Executable : 
void NeuLAND_Comtamination()      //%% Executable : 
{
  LoadNLgCut();

  TCut deltE = "ncedep<70&&ncedep>50";

  auto rChain0 = rChain[0];

  auto htof_all0 = new TH1D("htof_all0","TOF veto_all=0",150,0.,150.);  // pure neutron spectrum
  auto htof_all1 = new TH1D("htof_all1","TOF veto_all=1",150,0.,150.);  // all charged particle + more neutron contami.
  auto htof_bar0 = new TH1D("htof_bar0","TOF veto_bar=0",150,0.,150.);
  auto htof_bar1 = new TH1D("htof_bar1","TOF veto_bar=1",150,0.,150.);
  auto htof_mid0 = new TH1D("htof_mid0","TOF veto_mid=0",150,0.,150.);
  auto htof_mid1 = new TH1D("htof_mid1","TOF veto_mid=1",150,0.,150.);
  auto htof_loose0 = new TH1D("htof_loose0","TOF veto_loose=0",150,0.,150.);  // more neutron + more proton contami
  auto htof_loose1 = new TH1D("htof_loose1","TOF veto_loose=1",150,0.,150.);  // 


  auto htofE_all0   = new TH2D("htofE_all0" ,"TOF vs Edep veto_all==0",200,0.,500.,200,0.,150.);
  auto htofE_all1   = new TH2D("htofE_all1" ,"TOF vs Edep veto_all==1",200,0.,500.,200,0.,150.);
  auto htofE_all0p  = new TH2D("htofE_all0p","p TOF vs Edep veto_all==0",200,0.,500.,200,0.,150.);
  auto htofE_all1p  = new TH2D("htofE_all1p","p TOF vs Edep veto_all==1",200,0.,500.,200,0.,150.);
  auto htofE_all0d  = new TH2D("htofE_all0d","d TOF vs Edep veto_all==0",200,0.,500.,200,0.,150.);
  auto htofE_all1d  = new TH2D("htofE_all1d","d TOF vs Edep veto_all==1",200,0.,500.,200,0.,150.);
  auto htofE_all0t  = new TH2D("htofE_all0t","t TOF vs Edep veto_all==0",200,0.,500.,200,0.,150.);
  auto htofE_all1t  = new TH2D("htofE_all1t","t TOF vs Edep veto_all==1",200,0.,500.,200,0.,150.);
  auto htofE_all0tn = new TH2D("htofE_all0tn","t TOF vs Edep veto_all==0",200,0.,500.,200,0.,150.);
  auto htofE_all1tn = new TH2D("htofE_all1tn","t TOF vs Edep veto_all==1",200,0.,500.,200,0.,150.);

  auto htofE_bar1   = new TH2D("htofE_bar1","TOF vs Edep veto_bar==1",200,0.,500.,200,0.,150.);
  auto htofE_bar0   = new TH2D("htofE_bar0","TOF vs Edep veto_bar==0",200,0.,500.,200,0.,150.);
  auto htofE_bar0d  = new TH2D("htofE_bar0d","TOF vs Edep veto_bar==0",200,0.,500.,200,0.,150.);
  auto htofE_bar0t  = new TH2D("htofE_bar0t","TOF vs Edep veto_bar==0",200,0.,500.,200,0.,150.);
  auto htofE_mid1   = new TH2D("htofE_mid1","TOF vs Edep veto_mid==1",200,0.,500.,200,0.,150.);
  auto htofE_mid0   = new TH2D("htofE_mid0 "," TOF vs Edep veto_mid==0",200,0.,500.,200,0.,150.);
  auto htofE_mid0p  = new TH2D("htofE_mid0p","p TOF vs Edep veto_mid==0",200,0.,500.,200,0.,150.);
  auto htofE_mid0d  = new TH2D("htofE_mid0d","d TOF vs Edep veto_mid==0",200,0.,500.,200,0.,150.);
  auto htofE_mid0t  = new TH2D("htofE_mid0t","t TOF vs Edep veto_mid==0",200,0.,500.,200,0.,150.);

  auto htofE_loose1 = new TH2D("htofE_loose1 ","TOF vs Edep veto_loose==1",200,0.,500.,200,0.,150.);
  auto htofE_loose0 = new TH2D("htofE_loose0 ","TOF vs Edep veto_loose==0",200,0.,500.,200,0.,150.); 
  auto htofE_loose0p= new TH2D("htofE_loose0p","TOF vs Edep p veto_loose==0",200,0.,500.,200,0.,150.);
  auto htofE_loose0d= new TH2D("htofE_loose0d","TOF vs Edep d veto_loose==0",200,0.,500.,200,0.,150.);
  auto htofE_loose0t= new TH2D("htofE_loose0t","TOF vs Edep t veto_loose==0",200,0.,500.,200,0.,150.);

  auto htofE_alln   = new TH2D("htofE_alln","TOF vs Edep veto_all==0",400,0.,500.,400,0.,300.);
  auto htofE_barn   = new TH2D("htofE_barn","TOF vs Edep veto_bar==0",400,0.,500.,400,0.,300.);
  auto htofE_midn   = new TH2D("htofE_midn","TOF vs Edep veto_mid==0",400,0.,500.,400,0.,300.);
  auto htofE_loosen = new TH2D("htofE_loosen","TOF vs Edep veto_loose==0",400,0.,500.,400,0.,300.);

  auto htofE_loosep = new TH2D("htofEd_loosep","TOF vs Edep veto_loose==0 and prot",100,0.,500.,100,0.,150.);
  auto htofE_loosed = new TH2D("htofEd_loosed","TOF vs Edep veto_loose==0 and deut",100,0.,500.,100,0.,150.);
  auto htofE_looset = new TH2D("htofEd_looset","TOF vs Edep veto_loose==0 and trit",100,0.,500.,100,0.,150.);



  rChain0->Project("htof_all1","nctof",deltE&&"ncveto_all==1");
  rChain0->Project("htof_all0","nctof",deltE&&"ncveto_all==0");
  rChain0->Project("htof_bar1","nctof",deltE&&"ncveto_bar==1");
  rChain0->Project("htof_bar0","nctof",deltE&&"ncveto_bar==0");
  rChain0->Project("htof_mid1","nctof",deltE&&"ncveto_mid==1");
  rChain0->Project("htof_mid0","nctof",deltE&&"ncveto_mid==0");
  rChain0->Project("htof_loose1","nctof",deltE&&"ncveto_loose==1");
  rChain0->Project("htof_loose0","nctof",deltE&&"ncveto_loose==0");


  rChain0->Project("htofE_all1"  ,"nctof:ncedep","ncveto_all==1");
  rChain0->Project("htofE_all0"  ,"nctof:ncedep","ncveto_all==0");
  rChain0->Project("htofE_all1p" ,"nctof:ncedep","gcutNLProton&&ncveto_all==1");
  rChain0->Project("htofE_all0p" ,"nctof:ncedep","gcutNLProton&&ncveto_all==0");
  rChain0->Project("htofE_all1d" ,"nctof:ncedep","gcutNLDeuteron&&ncveto_all==1");
  rChain0->Project("htofE_all0d" ,"nctof:ncedep","gcutNLDeuteron&&ncveto_all==0");
  rChain0->Project("htofE_all1t" ,"nctof:ncedep","gcutNLTriton&&ncveto_all==1");
  rChain0->Project("htofE_all0t" ,"nctof:ncedep","gcutNLTriton&&ncveto_all==0");
  rChain0->Project("htofE_all1tn","nctof:ncedep","!gcutNLNeutron&&gcutNLTriton&&ncveto_all==1");
  rChain0->Project("htofE_all0tn","nctof:ncedep","!gcutNLNeutron&&gcutNLTriton&&ncveto_all==0");

  rChain0->Project("htofE_bar1"  ,"nctof:ncedep","ncveto_bar==1");
  rChain0->Project("htofE_bar0"  ,"nctof:ncedep","ncveto_bar==0");
  rChain0->Project("htofE_mid1"  ,"nctof:ncedep","ncveto_mid==1");
  rChain0->Project("htofE_mid0"  ,"nctof:ncedep","ncveto_mid==0");

  rChain0->Project("htofE_bar0d"  ,"nctof:ncedep","!gcutNLNeutron&&gcutNLDeuteron&&ncveto_bar==0");
  rChain0->Project("htofE_bar0t"  ,"nctof:ncedep","!gcutNLNeutron&&gcutNLTriton&&ncveto_bar==0");

  rChain0->Project("htofE_mid0p" ,"nctof:ncedep","gcutNLProton&&ncveto_mid==0");
  rChain0->Project("htofE_mid0d" ,"nctof:ncedep","!gcutNLNeutron&&gcutNLDeuteron&&ncveto_mid==0");
  rChain0->Project("htofE_mid0t" ,"nctof:ncedep","!gcutNLNeutron&&gcutNLTriton&&ncveto_mid==0");

  rChain0->Project("htofE_loose1" ,"nctof:ncedep","ncveto_loose==1");
  rChain0->Project("htofE_loose0" ,"nctof:ncedep","ncveto_loose==0");
  rChain0->Project("htofE_loose0p","nctof:ncedep","gcutNLProton&&ncveto_loose==0");
  rChain0->Project("htofE_loose0d","nctof:ncedep","!gcutNLNeutron&&gcutNLDeuteron&&ncveto_loose==0");
  rChain0->Project("htofE_loose0t","nctof:ncedep","!gcutNLNeutron&&gcutNLTriton&&ncveto_loose==0");

  rChain0->Project("htofE_loosen","nctof:ncedep","gcutNLNeutron&&ncveto_loose==0");
  rChain0->Project("htofE_alln"  ,"nctof:ncedep","gcutNLNeutron&&ncveto_all==0");
  rChain0->Project("htofE_barn"  ,"nctof:ncedep","gcutNLNeutron&&ncveto_bar==0");
  rChain0->Project("htofE_midn"  ,"nctof:ncedep","gcutNLNeutron&&ncveto_mid==0");

  rChain0->Project("htofE_loosep","nctof:ncedep","gcutNLProton&&ncveto_loose==0");
  rChain0->Project("htofE_loosed","nctof:ncedep","!gcutNLNeutron&&gcutNLDeuteron&&ncveto_loose==0");
  rChain0->Project("htofE_looset","nctof:ncedep","!gcutNLNeutron&&gcutNLTriton&&ncveto_loose==0");


  // Normalization factor
  Double_t nall0   = htof_all0  ->Integral(45.,60.);
  Double_t nall1   = htof_all1  ->Integral(45.,60.);
  Double_t nmid0   = htof_mid0  ->Integral(45.,60.);
  Double_t nmid1   = htof_mid1  ->Integral(45.,60.);
  Double_t nbar0   = htof_bar0  ->Integral(45.,60.);
  Double_t nbar1   = htof_bar1  ->Integral(45.,60.);
  Double_t nloose0 = htof_loose0->Integral(45.,60.);
  Double_t nloose1 = htof_loose1->Integral(45.,60.);
  
  // Charged particle spectrum after neutron subtraction
  auto htof_purecharge = new TH1D((*htof_all1) - nall1/nall0*(*htof_all0));
  htof_purecharge->SetName("htof_purecharge");

  cout << " All normalization " << nall1/nall0 << endl;

  // all
  Double_t fct = 0.792;
  TH2D *mtofE_all0 = new TH2D( fct * (*htofE_all0) );
  TH2D *stofE_all1 = new TH2D( (*htofE_all1) - (*mtofE_all0) );
  stofE_all1->SetName("stofE_all1");

  TH2D *mtofE_all0p = new TH2D( fct * (*htofE_all0p) );
  TH2D *stofE_all1p = new TH2D( (*htofE_all1p) - (*mtofE_all0p) );
  stofE_all1p->SetName("stofE_all1p");

  TH2D *mtofE_all0d = new TH2D( fct * (*htofE_all0d) );
  TH2D *stofE_all1d = new TH2D( (*htofE_all1d) - (*mtofE_all0d) );
  stofE_all1d->SetName("stofE_all1d");

  TH2D *mtofE_all0t = new TH2D( fct * (*htofE_all0t) );
  TH2D *stofE_all1t = new TH2D( (*htofE_all1t) - (*mtofE_all0t) );
  TH2D *stofE_all0t = new TH2D( (*htofE_all0t) - (*mtofE_all0t) );
  stofE_all1t->SetName("stofE_all1t");

  cout << "ALL  Total number " << endl;
  cout << "proton   " << stofE_all1p->GetEntries() << endl;
  cout << "deuteron " << stofE_all1d->GetEntries() << endl;
  cout << "triton   " << stofE_all1t->GetEntries() << endl;
  
  Double_t rpt = stofE_all1p->GetEntries() / stofE_all1t->GetEntries();
  Double_t rdt = stofE_all1d->GetEntries() / stofE_all1t->GetEntries();
  cout << " p/t " << rpt << endl;
  cout << " d/t " << rdt << endl;

  cout << " ---- ratio  ----" << endl;
  cout << " all $ " << endl;
  Double_t rtt = htofE_all0tn->GetEntries()/htofE_all1tn->GetEntries(); 
  cout << " rato t_0/t_1(all) " << htofE_all0tn->GetEntries() << " / " << htofE_all1tn->GetEntries()  
       << " = " <<  rtt << endl;
  cout << " n = " << htofE_alln->GetEntries() << endl;
  cout << " p/n = " << stofE_all1p->GetEntries()*rtt / htofE_alln->GetEntries() << endl;
  
  // mid
  cout << " bar $" << endl;
  rtt = htofE_bar0t->GetEntries()/htofE_all1tn->GetEntries();
  cout << " rato t_0/t_1(all) " <<  rtt << endl;
  cout << " P_cont = " << stofE_all1p->GetEntries()*rtt  << endl;
  cout << " n = " << htofE_barn->GetEntries() << endl;
  cout << " p/n = " << stofE_all1p->GetEntries()*rtt / htofE_barn->GetEntries() << endl;

  // mid
  cout << " mid $" << endl;
  rtt = htofE_mid0t->GetEntries()/htofE_all1tn->GetEntries();
  cout << " rato t_0/t_1(all) " <<  rtt << endl;
  cout << " P_cont = " << stofE_all1p->GetEntries()*rtt  << endl;
  cout << " n = " << htofE_midn->GetEntries() << endl;
  cout << " p/n = " << stofE_all1p->GetEntries()*rtt / htofE_midn->GetEntries() << endl;

  // loose
  cout << " Loose $ " << endl;
  rtt = htofE_loose0t->GetEntries()/htofE_all1tn->GetEntries();
  cout << " rato t_0/t_1(all) " <<  rtt << endl;
  cout << " P_cont = " << stofE_all1p->GetEntries()*rtt  << endl;
  cout << " n = " << htofE_loosen->GetEntries() << endl;
  cout << " p_cont/n = " << stofE_all1p->GetEntries()*rtt / htofE_loosen->GetEntries() << endl;



  // 1d 
  cout << " ------ 1D -------------- " << endl;
  Double_t pcut[] = {67.7, 86.8};
  Double_t dcut[] = {84.4 ,106.};
  Double_t tcut[] = {101. ,119.};
  Double_t pure_p = htof_purecharge->Integral(pcut[0], pcut[1]);
  Double_t pure_d = htof_purecharge->Integral(dcut[0], dcut[1]);
  Double_t pure_t = htof_purecharge->Integral(tcut[0], tcut[1]);

  Double_t pure_pt = pure_p / pure_t;
  Double_t pure_dt = pure_d / pure_t;


  // Charged in loose cut
  auto htof_loose01  = new TH1D( nloose1/nall0 * (*htof_all0) );
  htof_loose01->SetName("htof_loose01");

  auto htof_loosesub = new TH1D( (*htof_loose1) - nloose1/nall0*(*htof_all0));
  htof_loosesub->SetName("htof_loosesub");

  auto htof_all_loose = new TH1D( nloose0/nall0 * (*htof_all0));
  htof_all_loose->SetName("htof_all_loose");
  cout << " neut all0 / loose0 " << nall0/nloose0 << endl;

  auto htof_all_loosesub = new TH1D( (*htof_loose0) - nloose0/nall0 * (*htof_all0) );
  htof_all_loosesub->SetName("htof_all_loosesub");

  

  Double_t a_los_p = htof_all_loosesub->Integral(pcut[0], pcut[1]);
  Double_t a_los_d = htof_all_loosesub->Integral(dcut[0], dcut[1]);
  Double_t a_los_t = htof_all_loosesub->Integral(tcut[0], tcut[1]);

  cout << " tof_all_loosesub scale :"  << nall0/nloose0  << endl;
  cout << " Number of proton   : " << a_los_p << " / " << pure_p << " = " << a_los_p/pure_p << endl;
  cout << " Number of deuteron : " << a_los_d << " / " << pure_d << " = " << a_los_d/pure_d << endl;
  cout << " Number of triton   : " << a_los_t << " / " << pure_t << " = " << a_los_t/pure_t << endl;


  cout << " Number of Neutron " << endl;
  cout << " loose " << htofE_loosen->GetEntries() << endl;
  cout << " all   " << htofE_alln->GetEntries() << endl;
  cout << " bar   " << htofE_barn->GetEntries() << endl;
  cout << " mid   " << htofE_midn->GetEntries() << endl;

  cout << "p contami in n with loose " << a_los_p/htofE_loosen->GetEntries() *100. << " %" << endl;


  // plot
  if(kFALSE){
    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

    htof_loose1 ->Draw();

    htof_loose01->SetLineColor(2);
    htof_loose01->SetFillStyle(3004);
    htof_loose01->SetFillColor(2);
    htof_loose01->Draw("same");

    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    htof_loosesub->Draw();
  


    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->SetLogy();

    htof_all_loose->SetLineColor(2);
    //htof_all_loose->SetFillStyle(3004);
    htof_all_loose->SetFillColor(2);
    htof_all_loose->Draw();
    //  htof_all0->Draw("same");
    htof_loose0->SetLineColor(4);
    htof_loose0->Draw("same");
    auto prgnl = new TLine(pcut[0],0, pcut[0], htof_all_loose->GetMaximum());
    auto prgnr = new TLine(pcut[1],0, pcut[1], htof_all_loose->GetMaximum());
    auto drgnl = new TLine(dcut[0],0, dcut[0], htof_all_loose->GetMaximum());
    auto drgnr = new TLine(dcut[1],0, dcut[1], htof_all_loose->GetMaximum());
    auto trgnl = new TLine(tcut[0],0, tcut[0], htof_all_loose->GetMaximum());
    auto trgnr = new TLine(tcut[1],0, tcut[1], htof_all_loose->GetMaximum());
  
    prgnl->SetLineColor(8);
    prgnl->Draw("AL");
    prgnr->SetLineColor(8);
    prgnr->Draw("AL");

    drgnl->SetLineColor(8);
    drgnl->Draw("AL");
    drgnr->SetLineColor(8);
    drgnr->Draw("AL");

    trgnl->SetLineColor(8);
    trgnl->Draw("AL");
    trgnr->SetLineColor(8);
    trgnr->Draw("AL");
  
    htof_purecharge->SetLineColor(8);
    htof_purecharge->Draw("same");


    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    htof_all_loosesub->Draw();

    prgnl->SetLineColor(8);
    prgnl->SetY2(htof_all_loosesub->GetMaximum());
    prgnl->Draw("AL");
    prgnr->SetLineColor(8);
    prgnr->SetY2(htof_all_loosesub->GetMaximum());
    prgnr->Draw("AL");

    drgnl->SetLineColor(8);
    drgnl->SetY2(htof_all_loosesub->GetMaximum());
    drgnl->Draw("AL");
    drgnr->SetLineColor(8);
    drgnr->SetY2(htof_all_loosesub->GetMaximum());
    drgnr->Draw("AL");

    trgnl->SetLineColor(8);
    trgnl->SetY2(htof_all_loosesub->GetMaximum());
    trgnl->Draw("AL");
    trgnr->SetLineColor(8);
    trgnr->SetY2(htof_all_loosesub->GetMaximum());
    trgnr->Draw("AL");
  }

  if(kFALSE){
    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc[ic]->Divide(2,2);


    UInt_t id = 1;
    cc[ic]->cd(id); 
    cc[ic]->GetPad(id)->SetLogz(); id++;

    htofE_all0->Draw("colz");
    cc[ic]->cd(id);
    cc[ic]->GetPad(id)->SetLogz(); id++;
    htofE_all1->Draw("colz");

    cc[ic]->cd(id); 
    cc[ic]->GetPad(id)->SetLogz(); id++;
    htofE_loose0->Draw("colz");

    cc[ic]->cd(id); 
    cc[ic]->GetPad(id)->SetLogz();id++;
    htofE_loose1->Draw("colz");
  }

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  htof_purecharge->Draw();

  ic++; 
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  cc[ic]->SetLogz();
  stofE_all1->Draw("colz");

  ic++; 
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  cc[ic]->SetLogz();
  stofE_all1p->Draw("colz");
  ic++; 
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  cc[ic]->SetLogz();
  stofE_all1d->Draw("colz");
  ic++; 
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  cc[ic]->SetLogz();
  stofE_all1t->Draw("colz");
}

void LoadNLgCut()
{
  TFile *_file0 = TFile::Open("data/gcutNLNeutron.root");
  TCutG *gcutNLNeutron =(TCutG*)_file0->Get("gcutNLNeutron");
  _file0->Close();
  //  gcutNLNeutron->Print();

  TFile *_file1 = TFile::Open("data/gcutNLProton.root");
  TCutG *gcutNLProton =(TCutG*)_file1->Get("gcutNLProton");
  _file1->Close();
  //  gcutNLProton->Print();

  TFile *_file2 = TFile::Open("data/gcutNLDeuteron.root");
  TCutG *gcutNLDeuteron =(TCutG*)_file2->Get("gcutNLDeuteron");
  _file2->Close();
  //  gcutNLDeuteron->Print();

  TFile *_file3 = TFile::Open("data/gcutNLTriton.root");
  TCutG *gcutNLTriton =(TCutG*)_file3->Get("gcutNLTriton");
  _file3->Close();
  //  gcutNLTriton->Print();

}

//  LocalWords:  TH2D htofE
