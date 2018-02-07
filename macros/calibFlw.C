#include "FlowFunctions.h"
#include "openFlw.C"

Int_t  seltrackID = 4;
UInt_t selReactionPlanef = 10;
Int_t  seltrack;

TCanvas *cc[12];

UInt_t sys[]              = {10, 10, 10, 10};
const UInt_t  nbin = 16;
  
UInt_t ic = -1;
TChain *rChain[4];

UInt_t m_bgn = 0;
UInt_t m_end = 1;
Int_t ntrack[7];
auto aArray = new TClonesArray("STParticle",100);

Int_t iVer;
TString sVer;
Int_t mtrack;

// pt dependence


void ShiftingCorrection(STParticle *apar);


void calibFlw()
{
  gROOT->Reset();

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
}


//________________________________//%% Executable : 
void ReCentering(UInt_t isel = 0) //%% Executable : Recentering calibration
{
  TH1D *hQx[4];
  TH1D *hQy[4];
  TF1  *fg[4];


  UInt_t ic = 0; UInt_t id = 1;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700, 500*m_end);
  cc[ic]->Divide(2,m_end);

  for(UInt_t m = m_bgn; m < m_end; m++){
    fg[m]  = new TF1(Form("fg_%d",m),"gaus",-30,30);;


    hQx[m] = new TH1D(Form("hQx%d",m),"; Qx",100,-10,10);
    hQy[m] = new TH1D(Form("hQy%d",m),"; Qy",100,-10,10);

    TString fName ;

    if(isel == 0){
      rChain[m]->Project(Form("hQx%d",m), "unitP_ave.X()");
      rChain[m]->Project(Form("hQy%d",m), "unitP_ave.Y()");
      fName = "ReCent" +  sysName[sys[m]] + ".data";
    }
    else {
      rChain[m]->Project(Form("hQx%d",m), "unitP_rot.X()");
      rChain[m]->Project(Form("hQy%d",m), "unitP_rot.Y()");
      fName = "rotReCent" +  sysName[sys[m]] + ".data";
    }

    gSystem->cd("db");

    std::fstream fout;  
    fout.open(fName, std::fstream::out);
    fout << " Axis  Constatnt       Mean      Signma " << sysName[sys[m]]  
	 << endl;

    cc[ic]->cd(id); id++;
    hQx[m]->Fit(Form("fg_%d",m));

    fout << " X: " 
	 << setw(10)  << fg[m]->GetParameter(0) 
	 << setw(15)  << fg[m]->GetParameter(1)
	 << setw(10)  << fg[m]->GetParameter(2)
	 << endl;

    cc[ic]->cd(id); id++;
    hQy[m]->Fit(Form("fg_%d",m));

    fout << " Y: " 
	 << setw(10)  << fg[m]->GetParameter(0) 
	 << setw(15)  << fg[m]->GetParameter(1)
	 << setw(10)  << fg[m]->GetParameter(2)
	 << endl;

    fout.close();   
    gSystem->cd("..");

  }
}

//________________________________//%% Executable : 
void flatten_iphi_mtrkthetabin()  //%% Executable : 
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

  const UInt_t mtrknbin=5;
  Double_t mtrkbin[mtrknbin+1];
  Double_t mtrk_min = 0;
  Double_t mtrk_max = 40;
  for(UInt_t n = 0; n < mtrknbin+2; n++)
    mtrkbin[n]   = mtrk_max/mtrknbin * n;

  TH2D *hbaiphi[2];
  TH1D *hbiphi[2];
  TH1D *haiphi[2];
  TH2D *hbthetaiphi[2];
  TH2D *hathetaiphi[2];
  TH2D *hbmtrkiphi[2];
  TH2D *hamtrkiphi[2];
 
  UInt_t im = 0;

  for(UInt_t m = m_bgn; m < m_end; m++){

    //    OutputTree(rChain[m]);

    hbaiphi[m] = new TH2D(Form("hbaiphi%d",m),  " #phi_{i} before and after; before #phi_{i} [rad]; after #phi_{i} [rad] ", 
			  400,-3.5,3.5,400,-3.5,3.5);
    hbiphi[m]  = new TH1D(Form("hbiphi%d",m),   " #phi_{i} before; Azimuthal angle [rad]"  , 400,-3.2,3.2);
    haiphi[m]  = new TH1D(Form("haiphi%d",m),   " #phi_{i} after ; Azimuthal angle [rad]"  , 400,-3.2,3.2);
    hbthetaiphi[m]= new TH2D(Form("hbthetaiphi%d",m), " before ; theta; #phi;  "           , 200,0.,1.6, 400,-3.2,3.2); 
    hathetaiphi[m]= new TH2D(Form("hathetaiphi%d",m), " after  ; theta; #phi;  "           , 200,0.,1.6, 400,-3.2,3.2); 

    hbmtrkiphi[m] = new TH2D(Form("hbmtrkiphi%d",m)," before ; Number of tracks; #phi"     , 40,0,40,400,-3.2,3.2);
    hamtrkiphi[m] = new TH2D(Form("hamtrkiphi%d",m)," after  ; Number of tracks; #phi"     , 40,0,40,400,-3.2,3.2);

    p_org = new TClonesArray("TVector3", 100);
    p_rot = new TClonesArray("TVector3", 100);

    rChain[m]->SetBranchAddress("ntrack",ntrack);
    rChain[m]->SetBranchAddress("STParticle",&aArray);

    // Flattening with a shifting method
    STFlowCorrection *flowcorr[mtrknbin+1][thetanbin+1];

    for(UInt_t j = 0; j < mtrknbin+1; j++){   
      for(UInt_t i = 0; i < thetanbin+1; i++)  { 

	flowcorr[j][i] = new STFlowCorrection(rChain[m], harm, m); 

	flowcorr[j][i]->SetBin_min(0, mtrkbin[j]);
	flowcorr[j][i]->SetBin_min(1, thetabin[i]);
	
	if(j < mtrknbin+1 && i < thetanbin+1){
	  flowcorr[j][i]->SetBin_max(0, mtrkbin[j+1]);
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

      UInt_t imtrk  = 0;

      UInt_t j = mtrknbin;

      //     Int_t seltrack = ntrack[4];
      Int_t seltrack = ntrack[5];

      while(1){ 
	if( seltrack >= mtrkbin[j] ){
	  imtrk = j;
	  break;
	}
	j--;
      }

      UInt_t itheta = 0;
      TIter next(aArray);
      STParticle *aPart1 = NULL;
      
      while( (aPart1 = (STParticle*)next()) ){


	//	if( aPart1->GetReactionPlaneFlag() == 10 ){
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

	  hbthetaiphi[m]->Fill(theta      , phi);	  
	  hbmtrkiphi[m] ->Fill(seltrack  , phi);	  

	  if(imtrk <= mtrknbin && itheta <= thetanbin) {
	    flowcorr[imtrk][itheta]->Add(seltrack, phi,theta);

	    aPart1->SetFlattenBinID(imtrk, itheta);

	  }
	}
      
	//      cflw->Fill();
      }
    }
    //----------  get corrrection parameters
    for(UInt_t j = 0; j < mtrknbin+1; j++){   

      for(UInt_t i = 0; i < thetanbin+1; i++) {   
	
	UInt_t nphi = flowcorr[j][i]->FourierCorrection();
	std::cout << " At " << mtrkbin[j] << " : " << thetabin[i] << std::endl;

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
	  hamtrkiphi[m] ->Fill(mtk.at(k)   , aphi.at(k));
	  haiphi[m]     ->Fill(aphi.at(k));	  
	}
	
	TString comm = Form("tmpcv%d.m%dn%d:flatten_iphi_mtrkthetabin; mtrack> %f && mtrack< %f theta> %f && theta< %f",
			    iVer,j,i,mtrkbin[j],mtrkbin[j+1],thetabin[i],thetabin[i+1]);
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
    if(hbmtrkiphi[m])   hbmtrkiphi[m]->Draw("colz");

    cc[im]->cd(iv); iv++;
    if(hathetaiphi[m])  hathetaiphi[m]->Draw("colz");

    cc[im]->cd(iv); iv++;
    if(hamtrkiphi[m])   hamtrkiphi[m]->Draw("colz");
    
    cc[im]->cd(iv); iv++;
    if(haiphi[m])       haiphi[m]->Draw();
    
    cc[im]->cd(1);
   
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

      cout << " dbase " << flowcorr[3][3]->GetFileName() << endl;
      cout << " m " << flowcorr[3][3]->GetBin_min(0) << " ~ " << flowcorr[3][3]->GetBin_max(0) 
	   << " t " << flowcorr[3][3]->GetBin_min(1) << " ~ " << flowcorr[3][3]->GetBin_max(1)
	   << endl;

      cout << " ------------" << endl;
      cout << " m " << mtrkbin[2]  << " ~ " << mtrkbin[3]
	   << " t " << thetabin[2] << " ~ " << thetabin[3]
	   << endl;

      vector<Double_t> vec1 =  flowcorr[3][3]->GetOriginalPhi();
      for(itr=vec1.begin(); itr!=vec1.end(); itr++)      
	hvphi->Fill(*itr);
      vec1.clear();

      hvphi->Draw();


      cc[im]->cd(iv); iv++;
      vec1 =  flowcorr[3][3]->GetTheta();
      for(itr=vec1.begin(); itr!=vec1.end(); itr++)      
	hvthet->Fill(*itr);
    
      hvthet->Draw();
    
      cc[im]->cd(iv); iv++;
      vector<Int_t> vec2 =  flowcorr[3][3]->GetMTrack();
      for(iitr = vec2.begin(); iitr != vec2.end(); iitr++)
	hvmtk->Fill(*iitr);

      hvmtk->Draw();

    }
    
    
    // fout->Close();
    // delete cflw;
    // delete fout;
  }
}


