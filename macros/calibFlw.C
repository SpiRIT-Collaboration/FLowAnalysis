#include "FlowFunctions.h"
#include "openFlw.C"

Int_t  seltrackID = 4;
UInt_t selReactionPlanef = 1000;
Int_t  seltrack;

TCanvas *cc[12];
UInt_t im = 0;

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

Int_t iVer;
TString sVer;
Int_t mtrack;

TString unitpX;
TString unitpY;
TString fhead;
TString fOutName ;

Double_t constX;
Double_t meanX ;
Double_t sigX  ;
Double_t constY;
Double_t meanY ;
Double_t sigY  ;

void Initialize();
void SetPsiCorrectionFileHeader(UInt_t isel);
void SetEnvironment();
void LoadReCenteringCorrection(UInt_t m);
void SaveReCenteringData(UInt_t m);
void Flatten_Psi_ntrackthetabin(UInt_t isel = 0);//%% Executable : 

void calibFlw()
{
  gStyle->SetOptStat(0);
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

  m_end = 1;

  std::cout << " ichain " << ichain << " m_end " << m_end << std::endl;
  
  gROOT->ProcessLine(".! grep -i void calibFlw.C | grep '//%%'");

  Flatten_Psi_ntrackthetabin(2);  //TVector3(unitP2_rot->X(), unitP2_rot->Y(), 0.);
  Flatten_Psi_ntrackthetabin(4);  //TVector3(unitP_1r->X(),   unitP_1r->Y(),   0.);
  //  Flatten_Psi_ntrackthetabin(6);  //TVector3(unitP_1r->Mod()*cos(bsPhi_1[0]), unitP_1r->Mod()*sin(bsPhi_1[0]), 0.);
}



//________________________________//%% Executable : 
void ReCentering(UInt_t isel = 2, Int_t nmin=0, Int_t nmax=100) //%% Executable : Recentering calibration
{
  auto fgX  = new TF1("fgX","gaus",-30,30);;
  auto fgY  = new TF1("fgY","gaus",-30,30);;

  auto hQx  = new TH1D(Form("hQx%d_%d",nmin,isel),"; Qx",100,-20,20);
  auto hQy  = new TH1D(Form("hQy%d_%d",nmin,isel),"; Qy",100,-20,20);

  TString htitle = Form("mult >= %d && mult < %d",nmin, nmax);
  hQx->SetTitle(htitle);
  hQy->SetTitle(htitle);


  TString cutdef;
  if(isel <= 3)
    cutdef = Form("ntrack[%d]>=%d && ntrack[%d]<%d",seltrackID,nmin,seltrackID,nmax);
  else
    cutdef = Form("mtrack_1>=%d   && mtrack_1<%d",nmin,nmax);

  TCut multcut = TCut(cutdef);

  rChain[0]->Project(Form("hQx%d_%d",nmin,isel), unitpX, multcut);
  rChain[0]->Project(Form("hQy%d_%d",nmin,isel), unitpY, multcut);


  hQx->Fit("fgX","Q0");
  hQy->Fit("fgY","Q0");

  constX= fgX->GetParameter(0);
  meanX = fgX->GetParameter(1);
  sigX  = fgX->GetParameter(2);

  constY= fgY->GetParameter(0);
  meanY = fgY->GetParameter(1);
  sigY  = fgY->GetParameter(2);

  delete hQx;
  delete hQy;
}


//________________________________//%% Executable : 
void Flatten_Psi_ntrackthetabin(UInt_t isel)
{
  std::cout << "flatten_Psi_ntrackbin " << isel << std::endl;

  std::cout << "From " << m_bgn << " to " << m_end << std::endl;

  Initialize();
  SetPsiCorrectionFileHeader(isel);
  

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
  Double_t ntrk_max[] = {60, 60, 60, 60, 30, 30, 30, 30};
  for(UInt_t n = 0; n < ntrknbin+2; n++)
    ntrkbin[n]   = ntrk_max[isel]/ntrknbin * n;

 

  auto hbaiphi = new TH2D(Form("hbaiphi_%d",isel),  " #Phi before and after; before #Phi [rad]; after #Phi [rad] ", 
			  400,-3.5,3.5,400,-3.5,3.5);
  auto hniphi  = new TH1D(Form("hniphi_%d",isel),   " #Phi no corr.; Azimuthal angle [rad]"  , 200,-3.2,3.2);
  auto hbiphi  = new TH1D(Form("hbiphi_%d",isel),   " #Phi before  ; Azimuthal angle [rad]"  , 200,-3.2,3.2);
  auto haiphi  = new TH1D(Form("haiphi_%d",isel),   " #Phi after   ; Azimuthal angle [rad]"  , 200,-3.2,3.2);

  auto hbthetaiphi= new TH2D(Form("hbthetaiphi_%d",isel), " before ; theta; #phi;  "           , 200,0.,theta_max, 400,-3.2,3.2); 
  auto hathetaiphi= new TH2D(Form("hathetaiphi_%d",isel), " after  ; theta; #phi;  "           , 200,0.,theta_max, 400,-3.2,3.2); 

  auto hbntrkiphi = new TH2D(Form("hbntrkiphi_%d",isel)," before ; Number of tracks; #phi"     , 
			     ntrk_max[isel],0,ntrk_max[isel],400,-3.2,3.2);
  auto hantrkiphi = new TH2D(Form("hantrkiphi_%d",isel)," after  ; Number of tracks; #phi"     , 
			     ntrk_max[isel],0,ntrk_max[isel],400,-3.2,3.2);
    
  auto habiphi = new TH2D(Form("habiphi_%d",isel)," ;#Psi_rot; #Psi_fc",200,-3.2,3.2,200,-3.2,3.2);

  auto unitP_ave  = new TVector3();
  auto unitP_rot  = new TVector3();
  auto unitP2_ave = new TVector2();
  auto unitP2_rot = new TVector2();
  auto unitP_1r   = new TVector2();
  auto unitP_2r   = new TVector2();
  UInt_t mtrack_1;
  UInt_t mtrack_2;
  Double_t bsPhi_1[3];
  Double_t bsPhi_2[3];


  rChain[0]->SetBranchAddress("ntrack",ntrack);
  rChain[0]->SetBranchAddress("unitP_ave",&unitP_ave);
  rChain[0]->SetBranchAddress("unitP_rot",&unitP_rot);


  if( isel >= 2) {
    rChain[0]->SetBranchAddress("unitP2_ave",&unitP2_ave);
    rChain[0]->SetBranchAddress("unitP2_rot",&unitP2_rot);
    rChain[0]->SetBranchAddress("unitP_1r"  ,&unitP_1r);
    rChain[0]->SetBranchAddress("unitP_2r"  ,&unitP_2r);
    rChain[0]->SetBranchAddress("mtrack_1"  ,&mtrack_1);
    rChain[0]->SetBranchAddress("mtrack_2"  ,&mtrack_2);
    rChain[0]->SetBranchAddress("bsPhi_1"   ,bsPhi_1);
    rChain[0]->SetBranchAddress("bsPhi_2"   ,bsPhi_2);
  }

  // Flattening with a shifting method
  STFlowCorrection *flowcorr[ntrknbin+1][thetanbin+1];

  for(UInt_t j = 0; j < ntrknbin+1; j++){   
    for(UInt_t i = 0; i < thetanbin+1; i++)  { 

      flowcorr[j][i] = new STFlowCorrection(rChain[0], harm, 0); 

      if(j < ntrknbin+1 && i < thetanbin+1){

	if(isel != 6)
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
    
  Int_t nevt = rChain[0]->GetEntries();
  cout << " Number of events " << nevt << endl;

  Int_t icout = 0;
  std::vector<Double_t> ophi;


  for(UInt_t i = 0; i < nevt; i++){
    rChain[0]->GetEntry(i);

    UInt_t intrk  = 0;

    UInt_t j = ntrknbin;
    
    Int_t seltrack;
    if( isel < 4 )        seltrack = ntrack[seltrackID];
    else if( isel == 4 )  seltrack = mtrack_1;
    else if( isel == 5 )  seltrack = mtrack_2;
    else if( isel == 6 )  seltrack = mtrack_1;
    else if( isel == 7 )  seltrack = mtrack_2;


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

    // see -> SetPsiCorrectionFileHeader()
    TVector3 vec = *unitP_rot;
    if(isel == 1) vec = *unitP_ave;
    else if(isel == 2) vec = TVector3(unitP2_rot->X(), unitP2_rot->Y(), 0.);
    else if(isel == 3) vec = TVector3(unitP2_ave->X(), unitP2_ave->Y(), 0.);
    else if(isel == 4) vec = TVector3(unitP_1r->X(),   unitP_1r->Y(),   0.); 
    else if(isel == 5) vec = TVector3(unitP_2r->X(),   unitP_2r->Y(),   0.); 
    else if(isel == 6) vec = TVector3(unitP_1r->Mod()*cos(bsPhi_1[0]), unitP_1r->Mod()*sin(bsPhi_1[0]), 0.);
    else if(isel == 7) vec = TVector3(unitP_2r->Mod()*cos(bsPhi_2[0]), unitP_2r->Mod()*sin(bsPhi_2[0]), 0.);

    hniphi->Fill(vec.Phi());

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
	  hbiphi     ->Fill(rcphi.at(k));
	  hbthetaiphi->Fill(atheta.at(k), rcphi.at(k));	  
	  hbntrkiphi ->Fill(mtk.at(k)   , rcphi.at(k));	  
	  
	  hathetaiphi->Fill(atheta.at(k), aphi.at(k));
	  hantrkiphi ->Fill(mtk.at(k)   , aphi.at(k));
	  haiphi     ->Fill(aphi.at(k));	  
	  
	  if(j == 1)
	    habiphi  ->Fill(ophi.at(k), aphi.at(k));
	}
      }
      
      // finename :: comm1(fhead+"cv%d.m%dn%d);
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
    if(hbthetaiphi)  hbthetaiphi->Draw("colz");
    
    cc[im]->cd(iv); iv++;
    if(hbntrkiphi)   hbntrkiphi->Draw("colz");
    
    cc[im]->cd(iv); iv++;
    if(hathetaiphi)  hathetaiphi->Draw("colz");
    
    cc[im]->cd(iv); iv++;
    if(hantrkiphi)   hantrkiphi->Draw("colz");
  }
    
   
  im++;
  cc[im] = new TCanvas(Form("cc%d",im),Form("cc%d",im),700,500);
  
  hniphi->SetLineColor(2);
  hniphi->Draw("e");
  
  hbiphi->SetLineColor(8);
  hbiphi->Draw("samee");

  haiphi->SetLineColor(4);
  haiphi->Draw("samee");

  auto aLeg = new TLegend(0.75,0.13,0.9,0.3,"");
  aLeg->AddEntry(hniphi,"No Collection ","lp");
  aLeg->AddEntry(hbiphi,"ReCentering","lp");
  aLeg->AddEntry(haiphi,"ReCentering & Shifting","lp");
  
  aLeg->Draw();
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
    cout << "$ RUN0={$RNF###} DB0=\"_rf.v#.#\" VER=# root calibFlw.C " << endl;
    exit(0);
  }

  iVer = atoi(sVer);

  cout << " VER " << iVer << endl;
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

void SetPsiCorrectionFileHeader(UInt_t isel){
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
  case 6:
    unitpX = "unitP_1r.Mod()*cos(bsPhi_1[0])";
    unitpY = "unitP_1r.Mod()*sin(bsPhi_1[0])";
    fhead += "bs_1";
    break;
  case 7:
    unitpX = "unitP_2r.Mod()*cos(bsPhi_2[0])";
    unitpY = "unitP_2r.Mod()*sin(bsPhi_2[0])";
    fhead += "bs_2";
    break;
  }


}

void Initialize()
{
  constX= 0.;
  meanX = 0.;
  sigX  = 1.;
  constY= 0.;
  meanY = 0.;
  sigY  = 1.;

}
