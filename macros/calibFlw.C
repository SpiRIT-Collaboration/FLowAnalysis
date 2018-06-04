#include "FlowFunctions.h"
#include "openFlw.C"

Int_t  seltrackID = 4;
UInt_t selReactionPlanef = 1000;
Int_t  seltrack;

TCanvas *cc[12];

UInt_t sys[]              = {10, 10, 10, 10};

Double_t xconst[4];
Double_t xmean[4];
Double_t xsig[4];
Double_t yconst[4];
Double_t ymean[4];
Double_t ysig[4];


const UInt_t  nbin = 16;
  
UInt_t ic = -1;
Int_t ntrack[7];
auto aArray = new TClonesArray("STParticle",100);

Int_t iVer;
TString sVer;
Int_t mtrack;

TString unitpX;
TString unitpY;
TString fhead;
TString fOutName ;
Double_t constX= 0.;
Double_t meanX = 0.;
Double_t sigX  = 1.;
Double_t constY= 0.;
Double_t meanY = 0.;
Double_t sigY  = 1.;

// pt dependence
void SetPsiCorrectionFileHeader(UInt_t isel);
void SetEnvironment();
void ShiftingCorrection(STParticle *apar);
void LoadReCenteringCorrection(UInt_t m);
void SaveReCenteringData(UInt_t m);
void Flatten_Psi_ntrackthetabin(UInt_t isel = 0);//%% Executable : 
void calibFlw()
{
  gROOT->Reset();

  SetEnvironment();

  UInt_t ichain = 0;

  openFlw();

  for(UInt_t i = 0; i < 4; i++){
    rChain[ichain] = (TChain*)gROOT->FindObject(Form("rChain%d",i));
    if(rChain[ichain] != NULL) {    
      sys[ichain] = GetSystem(ichain);
      std::cout << " System " << ichain << " "  << sys[ichain] << "  -> " << sysName[sys[ichain]] << std::endl; 
      ichain++;
    }
  }

  if(rChain[0] == NULL)
    exit(0);

  m_end = ichain;

  std::cout << " ichain " << ichain << " m_end " << m_end << std::endl;
  
  gROOT->ProcessLine(".! grep -i void calibFlw.C | grep '//%%'");

  Flatten_Psi_ntrackthetabin(2);
  Flatten_Psi_ntrackthetabin(4);
}

//________________________________//%% Executable : 
void ReCentering(UInt_t isel = 2, Int_t nmin=0, Int_t nmax=100) //%% Executable : Recentering calibration
{
  TH1D *hQx[4];
  TH1D *hQy[4];
  TF1  *fgX[2];
  TF1  *fgY[2];


  ic++; UInt_t id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700, 500*m_end);
  cc[ic]->Divide(2,m_end);

  for(UInt_t m = m_bgn; m < m_end; m++){

    fgX[m]  = new TF1(Form("fgX%d_%d",nmin,m),"gaus",-30,30);;
    fgY[m]  = new TF1(Form("fgY%d_%d",nmin,m),"gaus",-30,30);;

    hQx[m] = new TH1D(Form("hQx%d_%d",nmin,m),"; Qx",100,-30,30);
    hQy[m] = new TH1D(Form("hQy%d_%d",nmin,m),"; Qy",100,-30,30);

    TString htitle = Form("mult >= %d && mult < %d",nmin, nmax);
    hQx[m]->SetTitle(htitle);
    hQy[m]->SetTitle(htitle);

    SetPsiCorrectionFileHeader(isel);

    
    TString cutdef;
    if(isel <= 3)
      cutdef = Form("ntrack[%d]>=%d && ntrack[%d]<%d",seltrackID,nmin,seltrackID,nmax);
    else
      cutdef = Form("mtrack_1>=%d   && mtrack_1<%d",nmin,nmax);

    TCut multcut = TCut(cutdef);

    rChain[m]->Project(Form("hQx%d_%d",nmin,m), unitpX, multcut);
    rChain[m]->Project(Form("hQy%d_%d",nmin,m), unitpY, multcut);

    cc[ic]->cd(id); id++;
    hQx[m]->Fit(Form("fgX%d_%d",nmin,m),"Q0");

    cc[ic]->cd(id); id++;
    hQy[m]->Fit(Form("fgY%d_%d",nmin,m),"Q0");

    constX= fgX[m]->GetParameter(0);
    meanX = fgX[m]->GetParameter(1);
    sigX  = fgX[m]->GetParameter(2);

    constY= fgY[m]->GetParameter(0);
    meanY = fgY[m]->GetParameter(1);
    sigY  = fgY[m]->GetParameter(2);

    //SaveReCenteringData(m);

  }

}

void SaveReCenteringData(UInt_t m)
{
  gSystem->cd("db");

  std::fstream fout;
  fOutName = "ReCent" +  sysName[sys[m]] + ".data";
  fout.open(fOutName, std::fstream::out);

  fout << "ReCentering for " << unitpX << " and " << unitpY << endl;
  fout << " Axis  Constatnt       Mean      Signma " << sysName[sys[m]]  
       << endl;
  cout << " X: " 
       << setw(10)  << constX 
       << setw(15)  << meanX
       << setw(10)  << sigX
       << endl;

  cout << " Y: " 
       << setw(10)  << constY
       << setw(15)  << meanY
       << setw(10)  << sigY
       << endl;

  fout << " X: " 
       << setw(10)  << constX 
       << setw(15)  << meanX
       << setw(10)  << sigX
       << endl;


  fout << " Y: " 
       << setw(10)  << constY
       << setw(15)  << meanY
       << setw(10)  << sigY
       << endl;

  fout.close();   
  gSystem->cd("..");

}

void SetPsiCorrectionFileHeader(UInt_t isel)
{
  fhead = "Psi";
  switch(isel){
  case 0:
    unitpX = "unitP_rot.X()";
    unitpY = "unitP_rot.Y()";
    fhead += "rt";
    break;
  case 1:
    unitpX = "unitP_ave.X()";
    unitpY = "unitP_ave.Y()";
    fhead += "av";
    break;
  case 2:
    unitpX = "unitP2_rot.X()";
    unitpY = "unitP2_rot.Y()";
    fhead += "2rt";
    break;
  case 3:
    unitpX = "unitP2_ave.X()";
    unitpY = "unitP2_ave.Y()";
    fhead += "2av";
    break;
  case 4:
    unitpX = "unitP_1r.X()";
    unitpY = "unitP_1r.Y()";
    fhead += "s1r";
    break;
  case 5:
    unitpX = "unitP_2r.X()";
    unitpY = "unitP_2r.Y()";
    fhead += "s2r";
    break;
  }
}


//________________________________//%% Executable : 
void Flatten_Psi_ntrackthetabin(UInt_t isel)
{
  std::cout << "flatten_Psi_ntrackbin" << std::endl;

  std::cout << "From " << m_bgn << " to " << m_end << std::endl;

  const UInt_t harm      = 5;
  const UInt_t thetanbin = 0;
  const UInt_t ntrknbin  = 5;
  
  // bin setting for theta
  Double_t thetabin[thetanbin+1];
  Double_t theta_min = 0.;
  Double_t theta_max = TMath::Pi()/2.;

  for(UInt_t n = 0; n < thetanbin+2; n++){
    if(thetanbin != 0) 
      thetabin[n]    = theta_max/(Double_t)thetanbin * (Double_t)n;

    else {
      if(n == 0)    thetabin[n]    = -1.;
      else          thetabin[n]    = TMath::Pi()/2.;
    }
  }

  // bin setting for multiplicity
  Double_t ntrkbin[ntrknbin+1];
  Double_t ntrk_min = 0;
  Double_t ntrk_max[] = {40, 40, 40, 40, 22, 22};
  for(UInt_t n = 0; n < ntrknbin+2; n++)
    ntrkbin[n]   = ntrk_max[isel]/ntrknbin * n;

  TH2D *hbaiphi[2];
  TH1D *hniphi[2];
  TH1D *hbiphi[2];
  TH1D *haiphi[2];
  TH2D *hbthetaiphi[2];
  TH2D *hathetaiphi[2];
  TH2D *hbntrkiphi[2];
  TH2D *hantrkiphi[2];
  TH2D *habiphi[2];
 
  UInt_t im = 0;


  for(UInt_t m = m_bgn; m < m_end; m++){

    hbaiphi[m] = new TH2D(Form("hbaiphi%d",m),  " #Phi before and after; before #Phi [rad]; after #Phi [rad] ", 
			  400,-3.5,3.5,400,-3.5,3.5);
    hniphi[m]  = new TH1D(Form("hniphi%d",m),   " #Phi no corr.; Azimuthal angle [rad]"  , 200,-3.2,3.2);
    hbiphi[m]  = new TH1D(Form("hbiphi%d",m),   " #Phi before  ; Azimuthal angle [rad]"  , 200,-3.2,3.2);
    haiphi[m]  = new TH1D(Form("haiphi%d",m),   " #Phi after   ; Azimuthal angle [rad]"  , 200,-3.2,3.2);

    hbthetaiphi[m]= new TH2D(Form("hbthetaiphi%d",m), " before ; theta; #phi;  "           , 200,0.,theta_max, 400,-3.2,3.2); 
    hathetaiphi[m]= new TH2D(Form("hathetaiphi%d",m), " after  ; theta; #phi;  "           , 200,0.,theta_max, 400,-3.2,3.2); 

    hbntrkiphi[m] = new TH2D(Form("hbntrkiphi%d",m)," before ; Number of tracks; #phi"     , 
			     ntrk_max[isel],0,ntrk_max[isel],400,-3.2,3.2);
    hantrkiphi[m] = new TH2D(Form("hantrkiphi%d",m)," after  ; Number of tracks; #phi"     , 
			     ntrk_max[isel],0,ntrk_max[isel],400,-3.2,3.2);
    
    habiphi[m] = new TH2D(Form("habiphi%d",m)," ;#Psi_rot; #Psi_fc",200,-3.2,3.2,200,-3.2,3.2);

    auto unitP_ave  = new TVector3();
    auto unitP_rot  = new TVector3();
    auto unitP2_ave = new TVector2();
    auto unitP2_rot = new TVector2();
    auto unitP_1r   = new TVector2();
    auto unitP_2r   = new TVector2();
    UInt_t mtrack_1;
    UInt_t mtrack_2;


    rChain[m]->SetBranchAddress("ntrack",ntrack);
    rChain[m]->SetBranchAddress("unitP_ave",&unitP_ave);
    rChain[m]->SetBranchAddress("unitP_rot",&unitP_rot);

    if( isel >= 2) {
      rChain[m]->SetBranchAddress("unitP2_ave",&unitP2_ave);
      rChain[m]->SetBranchAddress("unitP2_rot",&unitP2_rot);
      rChain[m]->SetBranchAddress("unitP_1r"  ,&unitP_1r);
      rChain[m]->SetBranchAddress("unitP_2r"  ,&unitP_2r);
      rChain[m]->SetBranchAddress("mtrack_1"  ,&mtrack_1);
      rChain[m]->SetBranchAddress("mtrack_2"  ,&mtrack_2);
    }

    // Flattening with a shifting method
    STFlowCorrection *flowcorr[ntrknbin+1][thetanbin+1];

    for(UInt_t j = 0; j < ntrknbin+1; j++){   
      for(UInt_t i = 0; i < thetanbin+1; i++)  { 

	flowcorr[j][i] = new STFlowCorrection(rChain[m], harm, m); 

	if(j < ntrknbin+1 && i < thetanbin+1){

	  ReCentering(isel, ntrkbin[j], ntrkbin[j+1]);

	  flowcorr[j][i]->SetBin_max(0, ntrkbin[j+1]);
	  flowcorr[j][i]->SetBin_max(1, thetabin[i+1]);
	}
	flowcorr[j][i]->SetBin_min(0, ntrkbin[j]);
	flowcorr[j][i]->SetBin_min(1, thetabin[i]);

	
	Double_t rcX[3];
	rcX[0] = constX;
	rcX[1] = meanX;
	rcX[2] = sigX;
	Double_t rcY[3];
	rcY[0] = constY;
	rcY[1] = meanY;
	rcY[2] = sigY;

	flowcorr[j][i]->SetReCenteringParameter("X",rcX);  
	flowcorr[j][i]->SetReCenteringParameter("Y",rcY);  
	
      }   
    }
    
    Int_t nevt = rChain[m]->GetEntries();
    cout << " Number of events " << nevt << endl;

    Int_t icout = 0;
    std::vector<Double_t> ophi;


    for(UInt_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      UInt_t intrk  = 0;

      UInt_t j = ntrknbin;

      Int_t seltrack;
      if( isel < 4 )        seltrack = ntrack[seltrackID];
      else if( isel == 4 )  seltrack = mtrack_1;
      else if( isel == 5 )  seltrack = mtrack_2;

      while(1){ 
	if( seltrack >= ntrkbin[j] ){
	  intrk = j;
	  break;
	}
	j--;
      }

      // if(intrk == 1){
      // 	cout << " ntrack " << seltrack
      // 	     << " bin " << intrk
      // 	     << " ntrkbin[" << j << "] " << ntrkbin[j]
      // 	     << endl;
      // }

      //------------------------
      UInt_t itheta = 0;

      TVector3 vec = *unitP_rot;
      if(isel == 1) vec = *unitP_ave;
      else if(isel == 2) vec = TVector3(unitP2_rot->X(), unitP2_rot->Y(), 0.);
      else if(isel == 3) vec = TVector3(unitP2_ave->X(), unitP2_ave->Y(), 0.);
      else if(isel == 4) vec = TVector3(unitP_1r->X(),   unitP_1r->Y(),   0.); 
      else if(isel == 5) vec = TVector3(unitP_2r->X(),   unitP_2r->Y(),   0.); 

      hniphi[m]->Fill(vec.Phi());

      // vec.SetX( (vec.X()-meanX)/sigX );
      // vec.SetY( (vec.Y()-meanY)/sigY );


      Double_t theta = vec.Theta();
      UInt_t k = thetanbin;
      while(1){ 
	if( theta >= thetabin[k] ){
	  itheta = k;
	  break;
	}
	k--;
      }

      if(intrk==1)
	ophi.push_back(vec.Phi());

      if(intrk <= ntrknbin && itheta <= thetanbin) {
	//	flowcorr[intrk][itheta]->Add(seltrack, phi,theta);

	flowcorr[intrk][itheta]->Add(seltrack, vec);
	
      }

      //-----------------------------
    }

    //----------  get corrrection parameters
    for(UInt_t j = 0; j < ntrknbin+1; j++){   

      for(UInt_t i = 0; i < thetanbin+1; i++) {   
	
	UInt_t nphi = flowcorr[j][i]->ReCenteringFourierCorrection();
	std::cout << " At " << ntrkbin[j] << " : " << thetabin[i] << " nphi " << nphi <<  std::endl;

	vector<Int_t>    mtk   = flowcorr[j][i]->GetMTrack();
	vector<Double_t> aphi  = flowcorr[j][i]->GetCorrectedPhi();
	vector<Double_t> atheta= flowcorr[j][i]->GetTheta();
	vector<Double_t> bphi  = flowcorr[j][i]->GetOriginalPhi();
	vector<Double_t> rcphi = flowcorr[j][i]->GetReCeneringPhi();

	if( aphi.size() != atheta.size() ){
	  std::cout << " size of pair doesn't match " << aphi.size() << " : " << atheta.size() << std::endl;
	  continue;
	}
	
	cout << "after " << mtk.size()  << endl;
	if(mtk.size() > 0){
	  for(UInt_t k = 0; k < (UInt_t)mtk.size(); k++){	  
	    hbiphi[m]     ->Fill(rcphi.at(k));
	    hbthetaiphi[m]->Fill(atheta.at(k), rcphi.at(k));	  
	    hbntrkiphi[m] ->Fill(mtk.at(k)   , rcphi.at(k));	  

	    hathetaiphi[m]->Fill(atheta.at(k), aphi.at(k));
	    hantrkiphi[m] ->Fill(mtk.at(k)   , aphi.at(k));
	    haiphi[m]     ->Fill(aphi.at(k));	  
	    
	    if(j == 1)
	      habiphi[m]  ->Fill(ophi.at(k), aphi.at(k));
	  }
	}

	TString comm1 = Form(fhead+"cv%d.m%dn%d:flatten_Psi_ntrkthetabin; ntrack>= %f && ntrack< %f theta>= %f && theta< %f",
			     iVer,j,i,ntrkbin[j],ntrkbin[j+1],thetabin[i],thetabin[i+1]);

	TString comm2 = unitpX + " && " + unitpY;
	cout << "save " << comm1 << endl;

	

	flowcorr[j][i]-> SaveCorrectionFactor(comm1, comm2);    
      }
    }

    if(kFALSE){  
      im++;
      cc[im] = new TCanvas(Form("cc%d",im),Form("cc%d",im),700,500);
      cc[im]->Divide(2,2);
    
      UInt_t iv = 1;
      cc[im]->cd(iv); iv++;
      if(hbthetaiphi[m])  hbthetaiphi[m]->Draw("colz");
    
      cc[im]->cd(iv); iv++;
      if(hbntrkiphi[m])   hbntrkiphi[m]->Draw("colz");

      cc[im]->cd(iv); iv++;
      if(hathetaiphi[m])  hathetaiphi[m]->Draw("colz");

      cc[im]->cd(iv); iv++;
      if(hantrkiphi[m])   hantrkiphi[m]->Draw("colz");
    }
    
    cc[im]->cd(1);
   
    im++;
    cc[im] = new TCanvas(Form("cc%d",im),Form("cc%d",im),700,500);


    hniphi[m]->SetLineColor(2);
    hniphi[m]->Draw("e");

    hbiphi[m]->SetLineColor(8);
    hbiphi[m]->Draw("samee");

    haiphi[m]->SetLineColor(4);
    haiphi[m]->Draw("samee");

    auto aLeg = new TLegend(0.15,0.7,0.3,0.9,"");
    aLeg->AddEntry(hniphi[m],"No Collection ","lp");
    aLeg->AddEntry(hbiphi[m],"ReCentering","lp");
    aLeg->AddEntry(haiphi[m],"ReCentering & Shifting","lp");

    aLeg->Draw();

    // im++;
    // cc[im] = new TCanvas(Form("cc%d",im),Form("cc%d",im),700,500);
    // habiphi[m]->Draw("colz");


    im++;
    if(m == 0){
      cc[im] = new TCanvas(Form("cc%d",im),Form("cc%d",im),700,500);
      cc[im]->Divide(2,2);
    
      iv = 1;
      cc[im]->cd(iv); iv++;
    
      auto hvphi  = new TH1D("hvphi"  ,"phi"   ,100,-3.2,3.2);
      auto hvthet = new TH1D("hvtheta","theta" ,100,0.,1.4);
      auto hvmtk  = new TH1I("hvmtk"  ,"mtrack", 60,0,60);
    
      vector<Double_t>::iterator itr;
      vector<Int_t>::iterator   iitr;

      cout << " dbase " << flowcorr[0][0]->GetFileName() << endl;
      cout << " m " << flowcorr[0][0]->GetBin_min(0) << " ~ " << flowcorr[0][0]->GetBin_max(0) 
	   << " t " << flowcorr[0][0]->GetBin_min(1) << " ~ " << flowcorr[0][0]->GetBin_max(1)
	   << endl;

      cout << " ------------" << endl;
      cout << " m " << ntrkbin[0]  << " ~ " << ntrkbin[1]
	   << " t " << thetabin[0] << " ~ " << thetabin[1]
	   << endl;

      vector<Double_t> vec1 =  flowcorr[0][0]->GetOriginalPhi();
      for(itr=vec1.begin(); itr!=vec1.end(); itr++)      
	hvphi->Fill(*itr);
      vec1.clear();

      hvphi->Draw();


      
      cc[im]->cd(iv); iv++;
      vec1 =  flowcorr[0][0]->GetTheta();
      for(itr=vec1.begin(); itr!=vec1.end(); itr++)      
	hvthet->Fill(*itr);
    
      hvthet->Draw();
    
      cc[im]->cd(iv); iv++;
      vector<Int_t> vec2 =  flowcorr[0][0]->GetMTrack();
      for(iitr = vec2.begin(); iitr != vec2.end(); iitr++)
	hvmtk->Fill(*iitr);

      hvmtk->Draw();

    }

    ///---- debug=----
    std::vector<Double_t> bfv = flowcorr[0][0]->GetOriginalPhi();
    std::vector<Double_t> crv = flowcorr[0][0]->GetCorrectedPhi();
    
    // for(UInt_t i = 0; i < (UInt_t)bfv.size(); i++) {
    //   if( bfv.at(i) >= 2.49 && bfv.at(i) < 2.5) 
    // 	cout << i << "th " << bfv.at(i) << " -> " << crv.at(i) << endl;
    // }
  }
}

//________________________________//%% Executable : 
void Flatten_iphi_ntrkthetabin()  //%% Not used  : 
{
  std::cout << "flatten_iphi_thetamtkbin" << std::endl;

  std::cout << "From " << m_bgn << " to " << m_end << std::endl;

  const UInt_t harm = 20;
  
  const UInt_t thetanbin = 40;
  Double_t thetabin[thetanbin+1];
  Double_t theta_min = 0.;
  Double_t theta_max = 1.4;
  for(UInt_t n = 0; n < thetanbin+2; n++)
    thetabin[n]    = theta_max/(Double_t)thetanbin * (Double_t)n;

  const UInt_t ntrknbin=5;
  Double_t ntrkbin[ntrknbin+1];
  Double_t ntrk_min = 0;
  Double_t ntrk_max = 40;
  for(UInt_t n = 0; n < ntrknbin+2; n++)
    ntrkbin[n]   = ntrk_max/ntrknbin * n;

  TH2D *hbaiphi[2];
  TH1D *hbiphi[2];
  TH1D *haiphi[2];
  TH2D *hbthetaiphi[2];
  TH2D *hathetaiphi[2];
  TH2D *hbntrkiphi[2];
  TH2D *hantrkiphi[2];
 
  UInt_t im = 0;

  for(UInt_t m = m_bgn; m < m_end; m++){

    //    OutputTree(rChain[m]);

    hbaiphi[m] = new TH2D(Form("hbaiphi%d",m),  " #phi_{i} before and after; before #phi_{i} [rad]; after #phi_{i} [rad] ", 
			  400,-3.5,3.5,400,-3.5,3.5);
    hbiphi[m]  = new TH1D(Form("hbiphi%d",m),   " #phi_{i} before; Azimuthal angle [rad]"  , 400,-3.2,3.2);
    haiphi[m]  = new TH1D(Form("haiphi%d",m),   " #phi_{i} after ; Azimuthal angle [rad]"  , 400,-3.2,3.2);
    hbthetaiphi[m]= new TH2D(Form("hbthetaiphi%d",m), " before ; theta; #phi;  "           , 200,0.,1.6, 400,-3.2,3.2); 
    hathetaiphi[m]= new TH2D(Form("hathetaiphi%d",m), " after  ; theta; #phi;  "           , 200,0.,1.6, 400,-3.2,3.2); 

    hbntrkiphi[m] = new TH2D(Form("hbntrkiphi%d",m)," before ; Number of tracks; #phi"     , 40,0,40,400,-3.2,3.2);
    hantrkiphi[m] = new TH2D(Form("hantrkiphi%d",m)," after  ; Number of tracks; #phi"     , 40,0,40,400,-3.2,3.2);

    p_org = new TClonesArray("TVector3", 100);
    p_rot = new TClonesArray("TVector3", 100);

    rChain[m]->SetBranchAddress("ntrack",ntrack);
    rChain[m]->SetBranchAddress("STParticle",&aArray);

    // Flattening with a shifting method
    STFlowCorrection *flowcorr[ntrknbin+1][thetanbin+1];

    for(UInt_t j = 0; j < ntrknbin+1; j++){   
      for(UInt_t i = 0; i < thetanbin+1; i++)  { 

	flowcorr[j][i] = new STFlowCorrection(rChain[m], harm, m); 

	flowcorr[j][i]->SetBin_min(0, ntrkbin[j]);
	flowcorr[j][i]->SetBin_min(1, thetabin[i]);
	
	if(j < ntrknbin+1 && i < thetanbin+1){
	  flowcorr[j][i]->SetBin_max(0, ntrkbin[j+1]);
	  flowcorr[j][i]->SetBin_max(1, thetabin[i+1]);
	}
      }   
    }
    
    Int_t nevt = rChain[m]->GetEntries();
    cout << " Number of events " << nevt << endl;

    Int_t icout = 0;

    //    cout << " at " << thetabin[thetanbin-1] << " " << thetabin[thetanbin] << endl;
    
    for(UInt_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      UInt_t intrk  = 0;

      UInt_t j = ntrknbin;

      //     Int_t seltrack = ntrack[4];
      Int_t seltrack = ntrack[5];

      while(1){ 
	if( seltrack >= ntrkbin[j] ){
	  intrk = j;
	  break;
	}
	j--;
      }

      //------------------------
      UInt_t itheta = 0;
      TIter next(aArray);
      STParticle *aPart1 = NULL;
      
      while( (aPart1 = (STParticle*)next()) ){

	if( aPart1->GetReactionPlaneFlag() == 20 ){

	  ShiftingCorrection(aPart1);
	  
	  Double_t phi   = aPart1->GetRotatedMomentum().Phi();
	  Double_t theta = aPart1->GetRotatedMomentum().Theta();
	  Double_t P     = aPart1->GetRotatedMomentum().Mag();

	  UInt_t k = thetanbin;
	  while(1){ 
	    if( theta >= thetabin[k] ){
	      itheta = k;
	      break;
	    }
	    k--;
	  }

	  hbthetaiphi[m]->Fill(theta     , phi);	  
	  hbntrkiphi[m] ->Fill(seltrack  , phi);	  

	  if(intrk <= ntrknbin && itheta <= thetanbin) {
	    flowcorr[intrk][itheta]->Add(seltrack, phi,theta);

	    aPart1->SetFlattenBinID(intrk, itheta);

	  }
	}
      
	//      cflw->Fill();
      }
      //-----------------------------
    }

    //----------  get corrrection parameters
    for(UInt_t j = 0; j < ntrknbin+1; j++){   

      for(UInt_t i = 0; i < thetanbin+1; i++) {   
	
	UInt_t nphi = flowcorr[j][i]->FourierCorrection();
	std::cout << " At " << ntrkbin[j] << " : " << thetabin[i] << std::endl;

	// if(nphi == 0) {
	//   std::cout << " no data is stored " << std::endl;
	//   continue;
	// }

	vector<Int_t>    mtk   = flowcorr[j][i]->GetMTrack();
	vector<Double_t> aphi  = flowcorr[j][i]->GetCorrectedPhi();
	vector<Double_t> atheta= flowcorr[j][i]->GetTheta();

	// std::cout << " size of pair  " << aphi.size() << " : " << atheta.size() << " / " << nphi << std::endl;
	// std::cout << " mtrack mean " << flowcorr[j][i]->GetMTrackMean() << " theta mean " << flowcorr[j][i]->GetThetaMean() << endl;

	if( aphi.size() != atheta.size() ){
	  std::cout << " size of pair doesn't match " << aphi.size() << " : " << atheta.size() << std::endl;
	  continue;
	}
	
	for(UInt_t k = 0; k < (UInt_t)mtk.size(); k++){	  
	  hathetaiphi[m]->Fill(atheta.at(k), aphi.at(k));
	  hantrkiphi[m] ->Fill(mtk.at(k)   , aphi.at(k));
	  haiphi[m]     ->Fill(aphi.at(k));	  
	}
	
	TString comm = Form("phicv%d.m%dn%d:flatten_iphi_ntrkthetabin; mtrack> %f && mtrack< %f theta> %f && theta< %f",
			    iVer,j,i,ntrkbin[j],ntrkbin[j+1],thetabin[i],thetabin[i+1]);
	cout << "save " << comm << endl;
	flowcorr[j][i]-> SaveCorrectionFactor(comm);    
      }
    }
  

    cc[im] = new TCanvas(Form("cc%d",im),Form("cc%d",im),700,1000);
    cc[im]->Divide(2,3);
    
    UInt_t iv = 1;
    cc[im]->cd(iv); iv++;
    if(hbthetaiphi[m])  hbthetaiphi[m]->Draw("colz");
    
    cc[im]->cd(iv); iv++;
    if(hbntrkiphi[m])   hbntrkiphi[m]->Draw("colz");

    cc[im]->cd(iv); iv++;
    if(hathetaiphi[m])  hathetaiphi[m]->Draw("colz");

    cc[im]->cd(iv); iv++;
    if(hantrkiphi[m])   hantrkiphi[m]->Draw("colz");
    
    cc[im]->cd(iv); iv++;
    if(haiphi[m])       haiphi[m]->Draw();
    
    cc[im]->cd(1);
   
    im++;


    // if(m == 0){
    //   cc[im] = new TCanvas(Form("cc%d",im),Form("cc%d",im),700,500);
    //   cc[im]->Divide(2,2);
    
    //   iv = 1;
    //   cc[im]->cd(iv); iv++;
    
    //   auto hvphi  = new TH1D("hvphi"  ,"phi"   ,100,-3.2,3.2);
    //   auto hvthet = new TH1D("hvtheta","theta" ,100,0.,1.4);
    //   auto hvmtk  = new TH1I("hvmtk"  ,"mtrack", 60,0,60);
    
    //   vector<Double_t>::iterator itr;
    //   vector<Int_t>::iterator   iitr;

    //   cout << " dbase " << flowcorr[3][3]->GetFileName() << endl;
    //   cout << " m " << flowcorr[3][3]->GetBin_min(0) << " ~ " << flowcorr[3][3]->GetBin_max(0) 
    // 	   << " t " << flowcorr[3][3]->GetBin_min(1) << " ~ " << flowcorr[3][3]->GetBin_max(1)
    // 	   << endl;

    //   cout << " ------------" << endl;
    //   cout << " m " << ntrkbin[2]  << " ~ " << ntrkbin[3]
    // 	   << " t " << thetabin[2] << " ~ " << thetabin[3]
    // 	   << endl;

    //   vector<Double_t> vec1 =  flowcorr[3][3]->GetOriginalPhi();
    //   for(itr=vec1.begin(); itr!=vec1.end(); itr++)      
    // 	hvphi->Fill(*itr);
    //   vec1.clear();

    //   hvphi->Draw();


    //   cc[im]->cd(iv); iv++;
    //   vec1 =  flowcorr[3][3]->GetTheta();
    //   for(itr=vec1.begin(); itr!=vec1.end(); itr++)      
    // 	hvthet->Fill(*itr);
    
    //   hvthet->Draw();
    
    //   cc[im]->cd(iv); iv++;
    //   vector<Int_t> vec2 =  flowcorr[3][3]->GetMTrack();
    //   for(iitr = vec2.begin(); iitr != vec2.end(); iitr++)
    // 	hvmtk->Fill(*iitr);

    //   hvmtk->Draw();

    // }
    
    
    // fout->Close();
    // delete cflw;
    // delete fout;
  }
}

void ShiftingCorrection(STParticle *apt)
{
  
  //  TVector2::Phi_mpi_pi(atan(/));
}

void LoadReCenteringCorrection(UInt_t m)
{

  TString sysname = sysName[sys[m]];
  TString fName = "ReCent" + sysname + ".data";

  gSystem->cd("db");
  std::fstream fin;
  fin.open(fName, std::fstream::in);
  if(fin == NULL) {
    std::cout << "A file " << fName << " was not found " << std::endl;
    exit(0);
  }
  
  TString sget;
  //  fin >> sget;
  
  while(!fin.eof()){
    fin >> sget;
    
    //    cout << " sget " << sget << endl;
    
    if(sget == "X:"){
      fin >> sget;
      xconst[m] = atof(sget);
      
      fin >> sget;
      xmean[m]  = atof(sget);
      
      fin >> sget;
      xsig[m]   = atof(sget);
    }

    if(sget == "Y:"){
      fin >> sget;
      yconst[m] = atof(sget);
      
      fin >> sget;
      ymean[m]  = atof(sget);
      
      fin >> sget;
      ysig[m]   = atof(sget);
    }
  }

  fin.close();
  gSystem->cd("..");



  cout << " X: const " << xconst[m] << " Mean " << xmean[m] << " sigma " << xsig[m] << endl;
  cout << " Y: const " << yconst[m] << " Mean " << ymean[m] << " sigma " << ysig[m] << endl;
}

void SetEnvironment()
{
  sVer = gSystem -> Getenv("VER");  // Version ID                                                                                                       
  if( sVer == "") {
    cout << " Please type " << endl;
    cout << "$ RUN0={####} DB0=#.# VER=# root flw_process3.C " << endl;
    exit(0);
  }

  iVer = atoi(sVer);

  cout << " VER " << iVer << endl;
}

