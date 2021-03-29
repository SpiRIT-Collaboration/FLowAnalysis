#include "unfoldLCPSpectra.h"

void makeReducedTree();
void unfold(TString momName="vapri", TString bgAssign="signal");
void unfold_slim(TString momName="vapri", TString bgAssign="signal");

int nIter=1;

Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u) 


// multiplicity cut
UInt_t trksel=0;
UInt_t trkcut[][2] = {{42, 52},
		      {48, 60},
		      {33, 48},
		      {24, 33},
		      { 0, 24}};
TString trklabel = Form("%02dto%02d",trkcut[trksel][0],trkcut[trksel][1]);

void unfoldLCPSpectra()
{
  //makeReducedTree();
  unfold_slim("vapri");
}

/**/
std::vector<int> mostCentralBins={0};
std::vector<int> iterBins={0,nIter-1};
std::vector<int> endIter={nIter-1};

const Double_t momCut[] = { 2000., 4000., 4200. };
/**/

TString nameId(int beam=-1, int pid=-1, int mbin=-1, int iter=-1, int ybin=-1)
{
	TString      name = beam!=-1?Form("_%dSn",beamA[beam]):"";
	if(pid!=-1)  name = name+"_"+pidName[pid];
	if(mbin!=-1) name = Form(name+"_mbin%d",mbin);
	if(iter!=-1) name = Form(name+"_iter%d",iter);
	if(ybin!=-1) name = Form(name+"_ybin%d",ybin);
	return name;
}

void makeReducedTree()
{
  auto WorkDir = gSystem->pwd();

  TFile *embedFile[nBeam][nParticle];
  TTree *mctrkTree[nBeam][nParticle], *trkTree[nBeam][nParticle], *eveTree[nBeam][nParticle];
  for(auto i: embedSys)for(auto j: embedPart){

      auto infile = embedFileName(i,j);
      std::cout << " Input Data " << infile << " will be opened." << std::endl;
      if( !gSystem->FindFile(".", infile )) {
	std::cout << " Input Data is not opened." << std::endl;
	continue; 
      }
      std::cout << " Input Data " << infile << " is opened." << std::endl;
      
      embedFile[i][j] = TFile::Open(embedFileName(i,j));
      embedFile[i][j]->GetObject("mctrktree",mctrkTree[i][j]);
      embedFile[i][j]->GetObject("trktree",trkTree[i][j]);
      embedFile[i][j]->GetObject("evetree",eveTree[i][j]);
    }

  gSystem->cd(WorkDir);
		
  double zet,aoq, vtxx,vtxy,vtxz;
  double mcpx,mcpy,mcpz;
  int    ggclose;
  int    ntrk,ngrtrk;
  double px,py,pz;
  double recopx,recopy,recopz;
  double vapripx,vapripy,vapripz;
  int    ndf,nclus,neclus;
  for(auto i: embedSys)for(auto j: embedPart){
      eveTree[i][j]->SetBranchAddress("ngrtrk",&ngrtrk); 
		
      mctrkTree[i][j]->SetBranchAddress("zet",&zet);
      mctrkTree[i][j]->SetBranchAddress("aoq",&aoq);
      mctrkTree[i][j]->SetBranchAddress("vtxx",&vtxx);
      mctrkTree[i][j]->SetBranchAddress("vtxy",&vtxy);
      mctrkTree[i][j]->SetBranchAddress("vtxz",&vtxz);
      mctrkTree[i][j]->SetBranchAddress("mcpx",&mcpx);
      mctrkTree[i][j]->SetBranchAddress("mcpy",&mcpy);
      mctrkTree[i][j]->SetBranchAddress("mcpz",&mcpz);
      mctrkTree[i][j]->SetBranchAddress("ggclose",&ggclose);
      mctrkTree[i][j]->SetBranchAddress("ntrk",&ntrk); 
	
      trkTree[i][j]->SetBranchAddress("zet",&zet);
      trkTree[i][j]->SetBranchAddress("aoq",&aoq);
      trkTree[i][j]->SetBranchAddress("vtxx",&vtxx);
      trkTree[i][j]->SetBranchAddress("vtxy",&vtxy);
      trkTree[i][j]->SetBranchAddress("vtxz",&vtxz);
      trkTree[i][j]->SetBranchAddress("px",&px);
      trkTree[i][j]->SetBranchAddress("py",&py);
      trkTree[i][j]->SetBranchAddress("pz",&pz);
      trkTree[i][j]->SetBranchAddress("vapripx",&vapripx);
      trkTree[i][j]->SetBranchAddress("vapripy",&vapripy);
      trkTree[i][j]->SetBranchAddress("vapripz",&vapripz);
      trkTree[i][j]->SetBranchAddress("recodist",&recodist);
      trkTree[i][j]->SetBranchAddress("recopx",&recopx);
      trkTree[i][j]->SetBranchAddress("recopy",&recopy);
      trkTree[i][j]->SetBranchAddress("recopz",&recopz);
      trkTree[i][j]->SetBranchAddress("ndf",&ndf);
      trkTree[i][j]->SetBranchAddress("nclus",&nclus);
      trkTree[i][j]->SetBranchAddress("neclus",&neclus);
      trkTree[i][j]->SetBranchAddress("mcpx",&mcpx);
      trkTree[i][j]->SetBranchAddress("mcpy",&mcpy);
      trkTree[i][j]->SetBranchAddress("mcpz",&mcpz);
      trkTree[i][j]->SetBranchAddress("ggclose",&ggclose); 
      trkTree[i][j]->SetBranchAddress("ngrtrk",&ngrtrk); 
      trkTree[i][j]->SetBranchAddress("ntrk",&ntrk); 
    }

  auto beamCutFile = TFile::Open("rootfiles/beamGate.root");
  TCutG* beamPICut[nBeam];
  for(auto i: embedSys)
    beamCutFile->GetObject(Form("sigma30_%dSn",beamA[i]),beamPICut[i]);

  auto vertexFile = TFile::Open("rootfiles/Vertex.root");
  TF1* fitVz[nBeam];
  Double_t vzPar[nBeam][4]={};
  for(auto i: embedSys){
    vertexFile->GetObject(Form("f1Vz_%dSn",beamA[i]),fitVz[i]);
    fitVz[i]->GetParameters(vzPar[i]);
  }

  // good event definition
  auto isGoodEvent = [&](int i)
    {
      bool isBeamGood = 1.; // @@@beamPICut[i]->IsInside(aoq,zet);
      bool iszVertex  = TMath::Abs(vtxz-vzPar[i][1])<=3.*vzPar[i][2];
      bool isxyVertex = TMath::Abs(vtxx)<=15.&&TMath::Abs(vtxy+205)<=20;
      return isBeamGood*iszVertex*isxyVertex*(!ggclose);
    };
	
  auto multFile = TFile::Open("rootfiles/Multiplicity.root");
  TH1 *h1B0[nBeam];
  for(auto i: embedSys)
    multFile->GetObject(Form("h1B0V_%dSn",beamA[i]),h1B0[i]);


  auto outFile   = new TFile(reducedFileName,"recreate");
  TTree *outMcTree[nBeam][nParticle];
  TTree *outRecoTree[nBeam][nParticle], *outVaTree[nBeam][nParticle], *outVapriTree[nBeam][nParticle];

  for(auto i: embedSys)for(auto j: embedPart){
      outMcTree[i][j]    = new TTree("reducedMcTree"+nameId(i,j),"");
      outRecoTree[i][j]  = new TTree("reducedRecoTree"+nameId(i,j),"");
      outVaTree[i][j]    = new TTree("reducedVaTree"+nameId(i,j),"");
      outVapriTree[i][j] = new TTree("reducedVapriTree"+nameId(i,j),"");
    }

  auto branch_reco_tree = [&](TTree* t, TString reco_type = "reco")
    {
      t->Branch(reco_type+"mom",&recomom);
      t->Branch(reco_type+"y",&recoy,reco_type+"y/D");
      t->Branch(reco_type+"yx",&recoyx,reco_type+"yx/D");
      t->Branch(reco_type+"pt",&recopt,reco_type+"pt/D");
      t->Branch(reco_type+"ut",&recout,reco_type+"ut/D");
      t->Branch(reco_type+"phi",&recophi,reco_type+"phi/D");
      t->Branch("recodist",&recodist,"recodist/D");
      t->Branch("ndf"  ,&ndf,  "ndf/I");
      t->Branch("nclus",&nclus,"nclus/I");
      t->Branch("ngrtrk",&ngrtrk,"ngrtrk/I");
      t->Branch("ntrk",&ntrk,"ntrk/I");
      t->Branch("b0",&b0,"b0/D");
      t->Branch("mcy",&mcy,"mcy/D");
      t->Branch("mcpt",&mcpt,"mcpt/D");
      t->Branch("mcut",&mcut,"mcut/D");
    };
	
  for(auto i: embedSys)for(auto j: embedPart){
      outMcTree[i][j]->Branch("mcmom",&mcmom);
      outMcTree[i][j]->Branch("mcy",&mcy,"mcy/D");
      outMcTree[i][j]->Branch("mcyx",&mcyx,"mcyx/D");
      outMcTree[i][j]->Branch("mcpt",&mcpt,"mcpt/D");
      outMcTree[i][j]->Branch("mcut",&mcut,"mcut/D");
      outMcTree[i][j]->Branch("mcphi",&mcphi,"mcphi/D");
      outMcTree[i][j]->Branch("ngrtrk",&ngrtrk,"ngrtrk/I");
      outMcTree[i][j]->Branch("ntrk",&ntrk,"ntrk/I");
      outMcTree[i][j]->Branch("b0",&b0,"b0/D");

      branch_reco_tree(outRecoTree[i][j],"reco");
      branch_reco_tree(outVaTree[i][j],"va");
      branch_reco_tree(outVapriTree[i][j],"vapri");
    }

  auto calc_yx = [](TLorentzVector v)
    {
      double beta_x = v.Beta()*TMath::Sin(v.Theta())*TMath::Cos(v.Phi());
      return 0.5*TMath::Log((1.+beta_x)/(1.-beta_x));
    };

  auto fill_reco_tree = [&](TTree* t, int i, int j, double reco_px, double reco_py, double reco_pz, double rand_rp)
    {
      TVector3 reco_mom(reco_px,reco_py,reco_pz);
      reco_mom.RotateY(invProjA[i]);

      if( j == 3 )
	reco_mom *= 2.;
      double reco_E = TMath::Sqrt(reco_mom.Mag2()+pidMass[j]*pidMass[j]);
      TLorentzVector reco_4mom(reco_mom,reco_E);
      TLorentzVector reco_4mom_rand_rp(reco_mom,reco_E); reco_4mom_rand_rp.RotateZ(rand_rp);
      if(reco_mom.Mag()>=rRange[j][0]){
	recophi = TMath::RadToDeg()*reco_mom.Phi();
	recomom = reco_4mom;
	recoy   = reco_4mom.Rapidity();
	recopt  = reco_4mom.Pt()/1000.;
	recoyx  = calc_yx(reco_4mom_rand_rp);
	recout  = reco_4mom.Pt()/pidMass[j]/u_p;
	t->Fill();
      }
    };

  for(auto i: embedSys)for(auto j: embedPart){
      std::cout<<" Analyzing "<<embedFile[i][j]->GetName()<<std::endl;
      auto entries = mctrkTree[i][j]->GetEntries();
      double randPhi;
      for(auto entry: MakeSeq(entries)){
	if(entry==0) std::cout<<" MC track Tree ana";
	if(entry%100000==0) std::cout<<" ."<<std::flush;
	if(entry==entries-1) std::cout<<" -> finish analysis!!"<<std::endl;
	mctrkTree[i][j]->GetEntry(entry);
	eveTree[i][j]->GetEntry(entry);

	if(!isGoodEvent(i))continue;

	b0 = h1B0[i]->GetBinContent(h1B0[i]->FindBin(ntrk));
			
	TVector3 mc_labmom(mcpx,mcpy,mcpz);
	mc_labmom.RotateY(invProjA[i]);
	mcphi = TMath::RadToDeg()*mc_labmom.Phi();

	Double_t E = TMath::Sqrt(mc_labmom.Mag2()+pidMass[j]*pidMass[j]);

	TLorentzVector mc_lab4mom(mc_labmom,E);
			
	mcmom = mc_lab4mom;
	mcy   = mc_lab4mom.Rapidity();
	mcpt  = mc_lab4mom.Pt()/1000.;
	mcut  = mc_lab4mom.Pt()/pidMass[j]/u_p;
		   
	randPhi = TMath::Pi()*gRandom->Uniform(-1.,1.);
	TLorentzVector mc_lab4mom_rot = mc_lab4mom;
	mc_lab4mom_rot.RotateZ(randPhi);
	mcyx = calc_yx(mc_lab4mom_rot);
		
	outMcTree[i][j]->Fill();
      }
      mctrkTree[i][j]->Delete();

      TVector3 oldMom(0.,0.,0.);
      entries = trkTree[i][j]->GetEntries();
      for(auto entry: MakeSeq(entries)){
	if(entry==0) std::cout<<" Reco. track Tree ana";
	if(entry%100000==0) std::cout<<" ."<<std::flush;
	if(entry==entries-1) std::cout<<" -> finish analysis!!"<<std::endl;
	trkTree[i][j]->GetEntry(entry);
			
	if(!isGoodEvent(i))continue;
			
	b0 = h1B0[i]->GetBinContent(h1B0[i]->FindBin(ntrk));
 			
	TVector3 initMom(mcpx,mcpy,mcpz);
	if( oldMom != initMom ){
	  randPhi = gRandom->Uniform(-1.,1.)*TMath::Pi();
	  oldMom = initMom;
	}
	initMom.RotateY(invProjA[i]);

	// if( j == 3 )
	//   initMom *= 2.;
	Double_t initE = TMath::Sqrt(initMom.Mag2()+pidMass[j]*pidMass[j]);
	TLorentzVector init4Mom(initMom,initE);
	mcy  = init4Mom.Rapidity();
	mcpt = init4Mom.Pt()/1000.;
	mcut = init4Mom.Pt()/pidMass[j]/u_p;

	//	if( !(recodist<=30&&nclus>=5&&neclus>=nclus*0.5) ) continue; //????

	TVector3 recoMom(recopx, recopy, recopz);

	fill_reco_tree(outRecoTree[i][j],i,j,recopx,recopy,recopz,recoMom.Phi());
	fill_reco_tree(outVaTree[i][j],i,j,px,py,pz,recoMom.Phi());
	fill_reco_tree(outVapriTree[i][j],i,j,vapripx,vapripy,vapripz,recoMom.Phi());
      }
      trkTree[i][j]->Delete();
    }
  outFile->Write();
}


void unfold_slim(
		 TString momName="vapri",
		 TString bgAssign="signal"
		 )
{

  auto rf = TFile::Open(reducedFileName);
  TTree *mcTree[nBeam][nParticle], *recoTree[nBeam][nParticle];
  for(auto i: embedSys)for(auto j: embedPart){
      rf->GetObject("reducedMcTree"+nameId(i,j),mcTree[i][j]);
      if(momName=="reco")       rf->GetObject("reducedRecoTree"+nameId(i,j),recoTree[i][j]);
      else if(momName=="va")    rf->GetObject("reducedVaTree"+nameId(i,j),recoTree[i][j]);
      else if(momName=="vapri") rf->GetObject("reducedVapriTree"+nameId(i,j),recoTree[i][j]);

      SetBranchAddressMCTree(mcTree[i][j]);
      SetBranchAddressRecoTree(recoTree[i][j],momName);
    }
	

  TFile *dataFile;
  TH2 *h2PtYData[nBeam][nParticle];
  TH2 *h2UtYData[nBeam][nParticle];
  TH2 *h2PhYData[nBeam][nParticle];

  dataFile = TFile::Open(inputFileName(1));
  cout << dataFile -> GetName() << endl;
  if( dataFile != NULL ) {

    for(auto i: embedSys)for(auto j: embedPart) {
	TString histname = histName(j,"_"+trklabel);
	
	cout << " histname " << histname << endl;

	dataFile->GetObject("hypt"+histname, h2PtYData[i][j]);
	dataFile->GetObject("hyut"+histname, h2UtYData[i][j]);
	dataFile->GetObject("hyph"+histname, h2PhYData[i][j]);

	if( h2PtYData[i][j] ) cout << h2PtYData[i][j]->GetName() << " is loaded. " << endl;
	else
	  cout << "h2PtYData[i][j]  is not loaded. " << endl;

	if( h2UtYData[i][j] ) cout << h2UtYData[i][j]->GetName() << " is loaded. " << endl;
	else
	  cout << "h2UtYData[i][j]  is not loaded. " << endl;

	if( h2PhYData[i][j] ) cout << h2PhYData[i][j]->GetName() << " is loaded. " << endl;
	else
	  cout << "h2PhYData[i][j]  is not loaded. " << endl;
      }
  }

	
  auto outFile = new TFile(corrFileName_slim( trklabel),"recreate");

  TH2 *h2PtYEmbed[nBeam][nParticle];
  TH2 *h2UtYEmbed[nBeam][nParticle];
  TH2 *h2PhYEmbed[nBeam][nParticle];

  // Event selection
  TCut multgate =  Form("ntrk>=%d && ntrk<=%d", trkcut[trksel][0], trkcut[trksel][1]);

  for(auto i: embedSys)for(auto j: embedPart){
      h2PtYEmbed[i][j] = new TH2D("h2PtYEmbed"+nameId(i,j),"Embed",40,-2,2,50,0.,2.5);
      h2UtYEmbed[i][j] = new TH2D("h2UtYEmbed"+nameId(i,j),"Embed",40,-2,2,50,0.,2.5);
      h2PhYEmbed[i][j] = new TH2D("h2PhYEmbed"+nameId(i,j),"Embed",40,-2,2,60,-3.15,3.15);

      mcTree[i][j]->Project(h2PtYEmbed[i][j]->GetName(),Form("mcpt:mcy/%lf-1",uNN[i]),phiCut[phicutID]+multgate);
      mcTree[i][j]->Project(h2UtYEmbed[i][j]->GetName(),Form("mcut:mcy/%lf-1",uNN[i]),phiCut[phicutID]+multgate);
      mcTree[i][j]->Project(h2PhYEmbed[i][j]->GetName(),Form("mcphi:mcy/%lf-1",uNN[i]),phiCut[phicutID]+multgate);

      h2UtYEmbed[i][j]->Print();
      h2PtYEmbed[i][j]->Print();
    }
	
  TH2 *h2PtYCorr[nBeam][nParticle][2];
  TH2 *h2UtYCorr[nBeam][nParticle][2];
  TH2 *h2PhYCorr[nBeam][nParticle][2];

  TH2 *h2PtYInit[nBeam][nParticle][2],  *h2PtYReco[nBeam][nParticle][2];
  TH2 *h2PtYEff[nBeam][nParticle][2],   *h2PtYEffSmear[nBeam][nParticle][2];
  TH2 *h2UtYInit[nBeam][nParticle][2],  *h2UtYReco[nBeam][nParticle][2];
  TH2 *h2UtYEff[nBeam][nParticle][2],   *h2UtYEffSmear[nBeam][nParticle][2];
  TH2 *h2PhYInit[nBeam][nParticle][2],  *h2PhYReco[nBeam][nParticle][2];
  TH2 *h2PhYEff[nBeam][nParticle][2],   *h2PhYEffSmear[nBeam][nParticle][2];
	
  for(auto i: embedSys)for(auto j: embedPart)for(auto iter: TSeqI(2)){
	h2PtYInit[i][j][iter] = new TH2D("h2PtYInit"+nameId(i,j,0,iter),"initial",40,-2,2,50,0.,2.5);
	h2PtYReco[i][j][iter] = new TH2D("h2PtYReco"+nameId(i,j,0,iter),"reco."  ,40,-2,2,50,0.,2.5);
	h2PtYInit[i][j][iter]->Sumw2();
	h2PtYReco[i][j][iter]->Sumw2();

	h2UtYInit[i][j][iter] = new TH2D("h2UtYInit"+nameId(i,j,0,iter),"initial",40,-2,2,50,0.,2.5);
	h2UtYReco[i][j][iter] = new TH2D("h2UtYReco"+nameId(i,j,0,iter),"reco."  ,40,-2,2,50,0.,2.5);
	h2UtYInit[i][j][iter]->Sumw2();
	h2UtYReco[i][j][iter]->Sumw2();

	h2PhYInit[i][j][iter] = new TH2D("h2PhiYInit"+nameId(i,j,0,iter),"initial",40,-2.,2.,60,-3.15,3.15);
	h2PhYReco[i][j][iter] = new TH2D("h2PhiYReco"+nameId(i,j,0,iter),"reco."  ,40,-2.,2.,60,-3.15,3.15);
	h2PhYInit[i][j][iter]->Sumw2();
	h2PhYReco[i][j][iter]->Sumw2();

      }
  
  TH2 *h2PtYWeight[nBeam][nParticle][2];
  TH2 *h2UtYWeight[nBeam][nParticle][2];
  TH2 *h2PhYWeight[nBeam][nParticle][2];

  // first: w=1.
  for(auto i: embedSys)for(auto j: embedPart)for(auto iter: TSeqI(1)){
	h2PtYWeight[i][j][iter] = (TH2D*)GaussianBlur(h2PtYData[i][j],3,1.);
	h2UtYWeight[i][j][iter] = (TH2D*)GaussianBlur(h2UtYData[i][j],3,1.);
	h2PhYWeight[i][j][iter] = (TH2D*)GaussianBlur(h2PhYData[i][j],3,1.);
      }

  for(auto i: embedSys)for(auto j: embedPart){
      for(auto iter: TSeqI(2)){
	// make efficiency function with weight

	std::cout<<" "<<embedNames[i]<<" "<<partNames[j]<<", iteration "<<iter<<std::endl;

	// mc track loop
	double weight = -1.;
	for(auto entry: TSeqL(mcTree[i][j]->GetEntries())){
	  mcTree[i][j]->GetEntry(entry);
				
	  if(ntrk < trkcut[trksel][0] || ntrk > trkcut[trksel][1]) continue;
	  
	  auto y0 = mcy/uNN[i]-1;
	  weight = getweight(h2PtYWeight[i][j][iter],y0,mcpt,h2PtYEmbed[i][j]);
	  h2PtYInit[i][j][iter]->Fill(y0,mcpt,weight);

	  weight = getweight(h2UtYWeight[i][j][iter],y0,mcut,h2UtYEmbed[i][j]);
	  h2UtYInit[i][j][iter]->Fill(y0,mcut,weight);

	  weight = getweight(h2PhYWeight[i][j][iter],y0,TMath::DegToRad()*mcphi,h2PhYEmbed[i][j]);
	  h2PhYInit[i][j][iter]->Fill(y0,TMath::DegToRad()*mcphi);				
	}

	// reco track loop
	for(auto entry: TSeqI(recoTree[i][j]->GetEntries())){
	  recoTree[i][j]->GetEntry(entry);

	  ///-----------------@@@@
	  if(ntrk < trkcut[trksel][0] || ntrk > trkcut[trksel][1]) continue;
	  if( recodist > cutdist ) continue;
	  if( ndf < cutndf ) continue;				

	  // w phicut
	  switch(phicutID) {
	  case 0:
	    if( abs(recophi) > 45 && abs(recophi) < 135) continue;
	    break;
	  case 1:
	    if( abs(recophi) > 45 ) continue;
	    break;
	  case 2:
	    if( abs(recophi) < 135 ) continue;
	    break;
	  case 3:
	    if( abs(recophi) < 45 || abs(recophi) > 135) continue;
	    break;
	  }
				
	  auto y0 = recoy/uNN[i]-1;

	  weight = getweight(h2PtYWeight[i][j][iter],y0,recopt,h2PtYEmbed[i][j]);
	  h2PtYReco[i][j][iter] ->Fill(y0,recopt,weight);

	  weight = getweight(h2UtYWeight[i][j][iter],y0,recout,h2UtYEmbed[i][j]);
	  h2UtYReco[i][j][iter] ->Fill(y0,recout,weight);

	  weight = getweight(h2PhYWeight[i][j][iter],y0,TMath::DegToRad()*recophi,h2PhYEmbed[i][j]);
	  h2PhYReco[i][j][iter]->Fill(y0,TMath::DegToRad()*recophi);
	}
			
	TString baseTitle = pidName[j]+"s from "+collName[i]+" "+bTitle[0]+";"+yTitle;

	// correction
	h2PtYEff[i][j][iter] = (TH2D*)make_hist_ratio(h2PtYReco[i][j][iter],h2PtYInit[i][j][iter],"h2PtYEff"+nameId(i,j,0,iter),"Eff. for "+baseTitle+";"+ptTitle);
	h2PtYCorr[i][j][iter] = (TH2D*)make_hist_ratio(h2PtYData[i][j],h2PtYEff[i][j][iter],"h2PtYCorr"+nameId(i,j,0,iter),baseTitle+";"+ptTitle);
			
	h2UtYEff[i][j][iter] = (TH2D*)make_hist_ratio(h2UtYReco[i][j][iter],h2UtYInit[i][j][iter],"h2UtYEff"+nameId(i,j,0,iter),"Eff. for "+baseTitle+";"+utTitle);
	h2UtYCorr[i][j][iter] = (TH2D*)make_hist_ratio(h2UtYData[i][j],h2UtYEff[i][j][iter],"h2UtYCorr"+nameId(i,j,0,iter),baseTitle+";"+utTitle);

	h2PhYEff[i][j][iter] = (TH2D*)make_hist_ratio(h2PhYReco[i][j][iter],h2PhYInit[i][j][iter],"h2PhYEff"+nameId(i,j,0,iter),"Eff. for "+baseTitle+";"+ptTitle);
	h2PhYCorr[i][j][iter] = (TH2D*)make_hist_ratio(h2PhYData[i][j],h2PhYEff[i][j][iter],"h2PhYCorr"+nameId(i,j,0,iter),baseTitle+";"+ptTitle);


	if(iter==0){
	  h2PtYWeight[i][j][iter+1] = (TH2D*)GaussianBlur(h2PtYCorr[i][j][iter],3,1.);
	  h2UtYWeight[i][j][iter+1] = (TH2D*)GaussianBlur(h2UtYCorr[i][j][iter],3,1.);
	  h2PhYWeight[i][j][iter+1] = (TH2D*)GaussianBlur(h2PhYCorr[i][j][iter],3,1.);
	}
      }
    }

  /* */
  auto latex = new TLatex();
  latex->SetTextSize(0.06);
  TCanvas *cvsUtYCorr[nBeam][nParticle];
  for(auto i: embedSys)for(auto j: embedPart) {
	TString name = "cvsUtYCorrecting"+nameId(i,j)+"."+momName+"."+outExtName;
	if(bgAssign=="signalV") name += ".slim";
	//if(isELT) name += ".corrELT";
	cvsUtYCorr[i][j] = new TCanvas(name,"",1000,600);
	cvsUtYCorr[i][j]->Divide(3,2);
	cvsUtYCorr[i][j]->Print(Form("fig/unfold/%s.pdf[",cvsUtYCorr[i][j]->GetName()),"pdf");

	for(auto iter: TSeqI(2)){
	  cout << " iter " << iter << endl;
	  cvsUtYCorr[i][j]->cd(1);
	  h2UtYInit[i][j][iter]->Draw("colz");
	  cvsUtYCorr[i][j]->cd(2);
	  h2UtYReco[i][j][iter]->Draw("colz");
	  cvsUtYCorr[i][j]->cd(3)->Clear();
	  latex->DrawLatexNDC(0.25,0.7,embedNames[i]+" "+partNames[j]);
	  latex->DrawLatexNDC(0.25,0.6,reducedFileName);
	  latex->DrawLatexNDC(0.25,0.5,Form("number of iteration: %d",iter));
	  latex->DrawLatexNDC(0.25,0.4,momName+" momentum.");

	  cvsUtYCorr[i][j]->cd(4);
	  h2UtYEff[i][j][iter]->SetMaximum(2.);
	  h2UtYEff[i][j][iter]->Draw("colz");
	  cvsUtYCorr[i][j]->cd(5);
	  //	  h2UtYData[i][j]->Draw("colz");
	  h2UtYWeight[i][j][iter]->Draw("colz");
	  cvsUtYCorr[i][j]->cd(6);

	  h2UtYCorr[i][j][iter]->SetMaximum( h2UtYCorr[i][j][iter]->GetMaximum()*0.5);
	  h2UtYCorr[i][j][iter]->Draw("colz");
			
	  cvsUtYCorr[i][j]->Print(Form("fig/unfold/%s_iter%d.pdf",cvsUtYCorr[i][j]->GetName(),iter),"pdf");

	}
    }


  //  cc->Print(Form("fig/unfold/%s.png",h2PhiYReco[i][0][0]->GetName()),"png");

	      
  /* */
  outFile->Write();
  cout << outFile->GetName() << " is created. " << endl;
}

