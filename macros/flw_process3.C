TCanvas *cc[6];

const Int_t nbinx = 30;

TChain *rChain[2];

TH2D *h2_r;
TH2D *h2_m;
TH1D *h2s_r[nbinx];
TH1D *h2s_m[nbinx];

TH2D *h2cos_r;
TH2D *h2cos_m;
TH1D *h2coss_r[nbinx];
TH1D *h2coss_m[nbinx];

TH2D *h2cos2_r;
TH2D *h2cos2_m;
TH1D *h2cos2s_r[nbinx];
TH1D *h2cos2s_m[nbinx];

TH1D *h1phimid_r;
TH1D *h1phimid_m;
TH1D *h1phimid_rm;
TF1  *f1;

TFile *fout;
TTree *cflw;



Int_t    ntrack[7];
Int_t    mtrack;
Int_t    mtrack_1;
Int_t    mtrack_2;

TClonesArray *aParticleArray;

Int_t             iRun=0;
Double_t          z = 0.;
Double_t          aoq = 0.;
vector<Double_t> *iphi=0;
vector<Double_t> *rpphi=0;
vector<Double_t> *rapid=0;
vector<Double_t> *prapid=0;
vector<Int_t>    *pid=0;
vector<Double_t> *pz=0;
vector<Double_t> *px=0;
vector<Double_t> *py=0;
vector<Double_t> *deltphi=0;
TClonesArray     *p_org=0;
TClonesArray     *p_rot=0;


TVector2 *unitP  =NULL;
TVector2 *unitP_1=NULL;
TVector2 *unitP_2=NULL;

TBranch  *bunitP;
TBranch  *bunitP_1;
TBranch  *bunitP_2;
TBranch *brpphi=0;
TBranch *biphi=0;
TBranch *brapid=0;
TBranch *bprapid=0;
TBranch *bpid=0;
TBranch *bpz=0;
TBranch *bpx=0;
TBranch *bpy=0;
TBranch *bdeltphi=0;

UInt_t m_bgn = 0;
UInt_t m_end = 1;
Int_t  ichain = 0;


void flatten_iphi_mtrkthetabin();
void OutputTree(TChain *rCh);
void flatten_Subevent();

void calcFlattenParameter()
{
  gROOT->Reset();

  gROOT->Macro("openRComp.C");

  rChain[ichain] = (TChain*)gROOT->FindObject(Form("rChain%d",ichain));
  if(rChain[ichain] == NULL ) exit(0);


  flatten_iphi_mtrkthetabin();

}

void flatten_iphi_mtrkthetabin()
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
  Double_t mtrk_max = 30;
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

    hbmtrkiphi[m] = new TH2D(Form("hbmtrkiphi%d",m)," before ; Number of tracks; #phi"     , 100,0,100,400,-3.2,3.2);
    hamtrkiphi[m] = new TH2D(Form("hamtrkiphi%d",m)," after  ; Number of tracks; #phi"     , 100,0,100,400,-3.2,3.2);

    aParticleArray = new TClonesArray("STParticle",100);

    p_org = new TClonesArray("TVector3", 100);
    p_rot = new TClonesArray("TVector3", 100);

    rChain[m]->SetBranchAddress("irun",&iRun);
    rChain[m]->SetBranchAddress("z",&z);
    rChain[m]->SetBranchAddress("aoq",&aoq);
    rChain[m]->SetBranchAddress("mtrack",&mtrack);
    rChain[m]->SetBranchAddress("ntrack",ntrack);

    rChain[m]->SetBranchAddress("STParticle",&aParticleArray);


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
      while(1){ 
	if( mtrack >= mtrkbin[j] ){
	  imtrk = j;
	  break;
	}
	j--;
      }

      UInt_t itheta = 0;
      TIter next(aParticleArray);
      STParticle *aPart1 = NULL;
      
      while( (aPart1 = (STParticle*)next()) ){

	// if(aPart1->GetReactionPlaneFlag() >= 11 && aPart1->GetReactionPlaneFlag() <= 13){
	// if(aPart1->GetReactionPlaneFlag() >= 10 && aPart1->GetRotatedMomentum().Mag()<2500){
	if(aPart1->GetRotatedMomentum().Mag() > 0){

	  
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

	  hbthetaiphi[m]->Fill(theta , phi);	  
	  hbmtrkiphi[m] ->Fill(mtrack  , phi);	  

	  if(imtrk <= mtrknbin && itheta <= thetanbin) {
	    flowcorr[imtrk][itheta]->Add(mtrack,phi,theta);

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
	
	TString comm = Form("m%dn%dmtktheta:flatten_iphi_mtrkthetabin; mtrack> %f && mtrack< %f theta> %f && theta< %f",
			    j,i,mtrkbin[j],mtrkbin[j+1],thetabin[i],thetabin[i+1]);
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

void OutputTree(TChain *rCh)
{

  TString fchain = ((TFile*)rCh->GetListOfFiles()->At(0))->GetFile()->GetTitle();
  
  Ssiz_t bgn = fchain.First("_");
  TString ssName = fchain(bgn-7,bgn+12);
  
  TString foutname = "../data/"+ssName+".c.root";

  fout = new TFile(foutname,"recreate");
  cflw  = new TTree("cflw","Flattened ");


  cout << "Output file is " << foutname << endl;


  //-- output                                                                                                       
  cflw->Branch("irun",&iRun,"irun/I");
  cflw->Branch("aoq",&aoq,"aoq/D");
  cflw->Branch("z",&z,"z/D");
  //  cflw->Branch("STVertex",&vertexArray);
  cflw->Branch("STParticle",&aParticleArray);
  cflw->Branch("ntrack",ntrack,"ntrack[7]/I");



}
