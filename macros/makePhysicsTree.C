#include "../analysisformat/STMassCalSimpleBB.hh"
#include  "../analysisformat/STRunToBeamA.hh"
//#include "../analysisformat/STParticle.hh"

std::map< Int_t, Int_t > beamAToSys = { {132, 0}, {108, 1}, {124, 2}, {112,3}};
const Int_t beamA[]   = {132,108,124,112};
TString lpid[] = {"1H",    "2H",      "3H"    ,"3He","4He","N"      ,"H"};  
TString pidName[]={"Proton","Deuteron","Triton","Helium3","Alpha"}; 
Double_t massu[]    = {938.2720813, 1875.612762, 2808.921112, 2808.39132, 3727.379378};

Double_t momCut[] = {2000.,4000., 4200., 4200.,4500.};
const Double_t mRange[][2] = {{0.,2100}, {700.,3000}, {1500.,4000.}, {1400.,4100.}, {2600.,4600.}, {4500.,7000.}};
const double yNN[]   = {0.3696, 0.3697, 0.3705, 0.3706};

UInt_t nmcut[]= {0,1,2};
Bool_t Msel[] = {0,0,0};//   0          1        2      3
UInt_t mtcut[4][4][2] = {{{56, 80}, {47, 56}, {0,47}, {47,80}},
			 {{55, 80}, {46, 55}, {0,46}, {46,80}},
			 {{56, 80}, {47, 56}, {0,47}, {47,80}},
			 {{55, 80}, {46, 55}, {0,46}, {46,80}}};    



const Int_t    nMultBin   = 4;
TH2* h2MassR[5][nMultBin][2];
TF1* f1MassGate[5][4][2][4][2];//!
TF1* f1Frac[5][nMultBin][100][2];
TF1* f1FracV[5][nMultBin][100][2];
TF1* f1Total[5][nMultBin][100][2];
TF1* f1TotalV[5][nMultBin][100][2]; 
Double_t vzPar[4]; 
Double_t vxbxPar[3];
Double_t vybyPar[3];
Int_t systemID = 0;

Bool_t bDEBUG = kFALSE;
//Bool_t bDEBUG = kTRUE;

Double_t Gauss(Double_t *x, Double_t *p)
{ return p[0]*TMath::Gaus(x[0],p[1],p[2]); } 
Double_t TripleGauss(Double_t *x, Double_t *p)
{ return p[0]*TMath::Gaus(x[0],p[1],p[2])+p[3]*TMath::Gaus(x[0],p[4],p[5])+p[6]*TMath::Gaus(x[0],p[7],p[8]); }
Double_t SigToTotalGauss(Double_t *x, Double_t *p)
{ return Gauss(x,p)/TripleGauss(x,p); }
Double_t TripleVoigt(Double_t *x, Double_t *p)
{ return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[2])+p[3]*TMath::Voigt(x[0]-p[4],p[5],p[5])+p[6]*TMath::Voigt(x[0]-p[7],p[8],p[8]); }
Double_t SigToTotalVoigt(Double_t *x, Double_t *p)
{ return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[2])/TripleVoigt(x,p); }


TString nameId(int beam, int pid=-1, int mbin=-1, int rbin=-1)
{
  TString      name = Form("_%dSn",beamA[beam]);
  if(pid!=-1)  name = name+"_"+pidName[pid];
  if(mbin!=-1) name = Form(name+"_mbin%d",mbin);
  if(rbin!=-1) name = Form(name+"_rbin%d",rbin);
  return name;
}

void SetupBeamVertex(UInt_t beam)
{
  auto vertexFile = TFile::Open("data/rootfiles/Vertex.root");
  TF1 *f1Vz;
  TF1 *f1VxBx, *f1VyBy;
  vertexFile->GetObject(Form("f1Vz_%dSn",  beam),f1Vz);
  vertexFile->GetObject(Form("f1VxBx_%dSn",beam),f1VxBx);
  vertexFile->GetObject(Form("f1VyBy_%dSn",beam),f1VyBy);
  f1Vz->GetParameters(vzPar);
  f1VxBx->GetParameters(vxbxPar);
  f1VyBy->GetParameters(vybyPar);
  vertexFile->Close();
}

Bool_t SetupPIDFit(UInt_t beam)
{
  std::vector< TString > label_rl = {"left","right"};

  UInt_t k = 0;
  for( auto il : label_rl ) {
    TString massFitName = Form("db/MassFit_%dSn.%s.root",beam,il.Data());
    auto massGateFile = TFile::Open(massFitName);
    if( massGateFile == NULL ) {
      LOG(ERROR) << massFitName << " is not found. " << FairLogger::endl;
      return kFALSE;
    }    
    LOG(INFO) << massGateFile->GetName() << " is opened. " << FairLogger::endl;
    

    for(auto pid: ROOT::MakeSeq(5)) 
      for(auto mbin: ROOT::MakeSeq(4)) 
	for(auto i: ROOT::TSeqI(2)) 
	  for(auto j: ROOT::TSeqI(4)) { 
	    TString f1name = Form("f1MassGate_%dSn_"+pidName[pid]+"_mbin%d_%d_%d",beam,mbin,i,j);
	    massGateFile->GetObject(f1name,  f1MassGate[pid][mbin][i][j][k]);

	    if( f1MassGate[pid][mbin][i][j][k] == NULL ) {
	      LOG(ERROR) << f1name << " is not found. " << FairLogger::endl;
	      return kFALSE;
	    }
	    f1name = f1MassGate[pid][mbin][i][j][k]->GetName()+il;
	    f1MassGate[pid][mbin][i][j][k]->SetName(f1name);
	    //	    cout << f1MassGate[pid][mbin][i][j][k]->GetName() << " is loaded. " << endl;
	  }
  

    int nRigiBins[5][2];
    for(auto pid: ROOT::MakeSeq(5))for(auto mbin: ROOT::MakeSeq(nMultBin)){
	massGateFile->GetObject("h2MassCalibR"+nameId(systemID,pid,mbin),h2MassR[pid][mbin][k]);
	TString h2name = h2MassR[pid][mbin][k]->GetName()+il;
	h2MassR[pid][mbin][k]->SetName(h2name);
	  
	int nRigiBin = h2MassR[pid][mbin][k]->GetNbinsX();

	for(auto rbin: ROOT::MakeSeq( nRigiBin )){
	  TString fname = nameId(systemID,pid,mbin,rbin);
	  massGateFile->GetObject("f1Total"+fname,f1Total[pid][mbin][rbin][k]);
	  f1Frac[pid][mbin][rbin][k] = new TF1("f1Frac"+fname,SigToTotalGauss,mRange[pid][0],mRange[pid][1],9);
	  f1Frac[pid][mbin][rbin][k]->SetNpx(2000);
	  for(auto par: ROOT::MakeSeq(9))
	    f1Frac[pid][mbin][rbin][k]->SetParameter(par,f1Total[pid][mbin][rbin][k]->GetParameter(par));

	  massGateFile->GetObject("f1TotalV"+fname,f1TotalV[pid][mbin][rbin][k]);
	  f1FracV[pid][mbin][rbin][k] = new TF1("f1FracV"+fname,SigToTotalVoigt,mRange[pid][0],mRange[pid][1],9);
	  f1FracV[pid][mbin][rbin][k]->SetNpx(2000);
	  for(auto par: ROOT::MakeSeq(9))
	    f1FracV[pid][mbin][rbin][k]->SetParameter(par,f1TotalV[pid][mbin][rbin][k]->GetParameter(par));
	}
      }
    
    massGateFile->Close();
    k++;
  }
  
  return kTRUE;
}

Int_t GetPIDFit(Int_t z, Double_t mass, TVector3 vMom, Int_t mult)
{
  // coded by Kaneko 20200418
  // pid: p=0, d=1, t=2, 3he=3, 4he=4. 
  // mass[0]: calibrated mass of track-> p,d,t: mass[1]:helium mass
  // fMom: rigidity magnitude of track
  // mbin: centrality bin, mbin=0: M>=56, mbin=1: 50<=M<56, mbin=2: 40<=M<50, mbin=3: no selection

  Double_t MassRange_Fit[6][2] = { {700.,1200},      // Kanekokun's pid range
				   {1500.,2300.},
				   {2400.,3400.},
				   {2400.,3200.},
				   {3250.,4200.},
				   {5200.,6000.} };   //!

  Int_t fpid = 0;
  if( mass <= 0 ) return 0;
  Bool_t fitCut = 1;  

  Double_t fMom = vMom.Mag();
  Double_t phi  = vMom.Phi()*TMath::RadToDeg();
  //-30<phi<20 left -150<phi or phi > 160 right
 
  UInt_t phibin = 1;
  if( phi <= 20 && phi >= -30 ) phibin = 0;

  UInt_t mbin = 0;
  if( mult >= 56 ) mbin = 0;
  else if( mult >= 50 && mult < 56 ) mbin = 1;
  else if( mult >= 40 && mult < 50 ) mbin = 2;
  else mbin = 3;

  Bool_t bfind = kFALSE;

  if( z == 1 ) {
    for(auto i : {0,1,2} ) {

      Bool_t roughCut =  mass >= MassRange_Fit[i][0] && mass <= MassRange_Fit[i][1]; 

      if( roughCut ) {
    
	Int_t fitSigmaID =  3 ;
	fitCut = mass >= f1MassGate[i][mbin][0][fitSigmaID][phibin]->Eval(fMom) && mass <= f1MassGate[i][mbin][1][fitSigmaID][phibin]->Eval(fMom);
	
	if( fitCut ) {
	  fpid = i;
	  //	 	  cout << " fpid " << fpid << " mass " << mass << endl;
	  bfind = kTRUE;
	  Break;
	}
	else if ( 0 ) {
	  cout << " mass " << mass 
	       << " p " << fMom
	       << " fitSigmaID " << fitSigmaID
	       << " mbin " << mbin
	       << " mult " << mult
	       << " " << f1MassGate[i][mbin][0][fitSigmaID][phibin]->Eval(fMom)
	       << " " << f1MassGate[i][mbin][1][fitSigmaID][phibin]->Eval(fMom)
	       << " phi " << phi
	       << " phibin " << phibin
	       << endl;
	}

      }
    }
  }
  else if( z == 2 ) {
    fMom *= 2.;
    for(auto i : {3,4} ) {
      Bool_t roughCut =  mass >= MassRange_Fit[i][0] && mass <= MassRange_Fit[i][1];
      if( roughCut ) {
    
	Int_t fitSigmaID =  2 ;
	if( i == 3 )
	  fitCut = mass >= f1MassGate[i][mbin][0][fitSigmaID-2][phibin]->Eval(fMom) && mass <= f1MassGate[i][mbin][1][fitSigmaID][phibin]->Eval(fMom);
    	else if( i == 4 )
	  fitCut = mass >= f1MassGate[i][mbin][0][fitSigmaID-1][phibin]->Eval(fMom) && mass <= f1MassGate[i][mbin][1][fitSigmaID+1][phibin]->Eval(fMom);
    
	if( fitCut ) {
	  fpid = i;
	  bfind = kTRUE;
	  Break;
	}
      }
    }
  }

  if( !bfind ) return 0;
  if( fMom < 100) return 0;

  STPID::PID pid = static_cast<STPID::PID>(fpid+1);
  Int_t rpid = STPID::GetPDG(pid) * fitCut;
    
  return rpid;
}

void makePhysicsTree()
{
  auto getSignalToTotal = [&](double m, TF1* f)
    {
      double frac = f->Eval(m);
      // when there are no background fit component, signal/total becomes a little bit larger than 1.                                               
      if(TMath::IsNaN(frac)||frac<0.) frac=0.;
      else if(frac>1.) frac=1.;
      return frac;
    };
  
  TString sRun = gSystem -> Getenv("RUN");
  UInt_t iRun = atoi(sRun);
  UInt_t bmA =  STRunToBeamA::GetBeamA(iRun);
  systemID = beamAToSys[bmA];
  
  cout << " beam is " << bmA << " and systemID is " << systemID <<  endl;
  
  auto massCalH  = new STMassCalSimpleBB("EmpiricalBB");
  auto massCalHe = new STMassCalSimpleBB("EmpiricalBB");

  TFile* calFile = TFile::Open(Form("db/PIDCalib_%dSn.root",bmA));
  if( bmA == 124 )
    calFile = TFile::Open(Form("db/PIDCalib_%dSn.root",132));

  TH2D *h2ParamH[7], *h2ParamHe[7];
  if( calFile  ) {

    for(auto i: ROOT::TSeqI(7)){
      calFile->GetObject(Form("h2InterpolateNM_%dSn_Par%d"  ,bmA,i),h2ParamH[i]);
      calFile->GetObject(Form("h2InterpolateHeNM_%dSn_Par%d",bmA,i),h2ParamHe[i]);
      if( bmA == 124 ) {
	calFile->GetObject(Form("h2InterpolateNM_%dSn_Par%d"  ,132,i),h2ParamH[i]);
	calFile->GetObject(Form("h2InterpolateHeNM_%dSn_Par%d",132,i),h2ParamHe[i]);
      }
    }
    massCalH  -> AddParameters(h2ParamH);
    massCalHe -> AddParameters(h2ParamHe);
  }
  //  calFile->Close();

  SetupBeamVertex(bmA);

  if( !SetupPIDFit(bmA)){
    cout << " Mass fit  data is not loaded. " << endl;
    return;
  }

  auto gcutFile = new TFile(Form("data/gcut%dSn.ROOT",bmA));
  auto gBeamCut =(TCutG*)gcutFile->Get("sigma20");
  gcutFile->Close();


  std::vector< UInt_t > sqPart = {0,1,2,3,4};
  //  std::vector< UInt_t > sqPart = {0};//,1,2,3,4};


  //@@@@@booking++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TH2D* hypt[5][3];
  TH1I* hmultsel[3];

  TFile* outfile[6];
  TH2D* hmassHp[5];
  TH2D* hmassHep[5];

  for( auto ipart : sqPart ) {
    outfile[ipart] = TFile::Open(Form("data/run%04d_ph_%s.v1.root",iRun,pidName[ipart].Data()),"recreate");
    hmassHp[ipart]    = new TH2D("hmassHp" ,lpid[ipart]+";p;mass",100, 0.,3000., 100,    0., 4000.);
    hmassHep[ipart]   = new TH2D("hmassHep",lpid[ipart]+";p;mass",100, 0.,5000., 100, 2500., 5500.);

    for( auto im : nmcut)
    hypt[ipart][im]    = new TH2D(Form("hypt_M%dto%d",mtcut[systemID][im][0],mtcut[systemID][im][1]), lpid[ipart]+";y;pt", 40,-2.,2.,50,0.,2.0);
  }
  for( auto im : nmcut )
    hmultsel[im] = new TH1I(Form("hmultsel_M%dto%d",mtcut[systemID][im][0],mtcut[systemID][im][1]), "multiplicity", 100,0,100);

  
  outfile[5] = TFile::Open(Form("data/run%04d_ph.v1.root",iRun),"recreate");
  auto hmult   = new TH1I("hmult"  ,"hmult" ,100,0,100);
  auto outtree = new TTree("cbmsim","Compressed");

  //---
  auto fChain = new TChain("cbmsim");

  TString rootDir = Form("/home/recoData/20200529/data/Sn%d/",bmA);

  //@@@@@Input++++++++++++++++++++++++++++++++++++++++++++++++++
  UInt_t i = 0;
  while(  1 ) {
    TString recoFile = Form("run%04d_s%02d.reco.develop.1988.bf2b00e.root",iRun,i);
    //    cout << " recoFile ... " << rootDir << recoFile << endl;
    if( i < 10 ) 
      recoFile = Form("run%04d_s%d.reco.develop.1988.bf2b00e.root",iRun,i);
      
    if(gSystem->FindFile(rootDir,recoFile)) {
      fChain -> Add(recoFile);
      cout << " LOAD " << recoFile << endl;
    }
    else
      break;
      
    i++;
  }
  cout << " Number of Entries " << fChain->GetEntries() << endl;
  if( fChain->GetEntries() == 0 ) return;

  UInt_t BeamIndex = 1;
  //   root [1] cbmsim->Print("toponly")
  // ******************************************************************************
  //     *Tree    :cbmsim    : /cbmout                                                *
  //     *Entries :    18360 : Total =     11939005302 bytes  File  Size = 6001466647 *
  //     *        :          : Tree compression factor =   1.99                       *
  // ******************************************************************************
  //     branch: EventHeader.            109462
  //     branch: STEventHeader           555190  fEventID
  //     branch: STCandList           1553875958
  //     branch: STRecoTrack          1881985287
  //     branch: STVertex              82141284
  //     branch: VATracks             1339163171
  //     branch: VACandList           1085440015
  //     branch: VAVertex              55616258
  //     branch: STBeamInfo             1216700
  //     branch: BDCVertex              1028291

  STEventHeader* eventHeader  = NULL;
  TClonesArray* trackArray    = NULL;
  TClonesArray* trackVAArray  = NULL;
  TClonesArray* vertexArray = NULL;
  TClonesArray* vertexVAArray = NULL;
  TClonesArray* vertexBDCArray= NULL;
  STBeamInfo*   beamInfo      = NULL;

  fChain -> SetBranchAddress("STEventHeader", &eventHeader);
  fChain -> SetBranchAddress("STRecoTrack",   &trackArray);
  fChain -> SetBranchAddress("VATracks",      &trackVAArray);

  fChain -> SetBranchAddress("VAVertex"   ,   &vertexVAArray);
  fChain -> SetBranchAddress("BDCVertex"  ,   &vertexBDCArray);
  fChain -> SetBranchAddress("STVertex"   ,   &vertexArray);
  fChain -> SetBranchAddress("STBeamInfo" ,   &beamInfo);



  UInt_t beamPID;
  UInt_t run;
  UInt_t event;
  UInt_t mtrack0, mtrack1, mtrack2;
  auto ParticleArray = new TClonesArray("STParticle",80);

  outfile[5]->cd();
  outtree->Branch("beamPID",&beamPID,"beamPID/I");
  outtree->Branch("run",&run,"run/I");
  outtree->Branch("event",&event,"event/I");
  outtree->Branch("mtrack0",&mtrack0,"mtrack0/I");
  outtree->Branch("mtrack1",&mtrack1,"mtrack1/I");
  outtree->Branch("mtrack2",&mtrack2,"mtrack2/I");
  outtree->Branch("Particle", &ParticleArray);


  auto apart = new STParticle();

  run = iRun;
  beamPID = bmA;

  for( auto ievt : ROOT::TSeqL( fChain->GetEntries()) ) {

    ParticleArray->Clear("C");
    
    //@@@@maxevent---
    //    if( ievt != 47 ) continue;
    //    if( ievt > 100 ) break;

    if(ievt%1000 == 0) 
      cout << " Process " << setw(4) << ievt << " / " << fChain->GetEntries() << endl;

    fChain->GetEntry(ievt);
    
    Bool_t sigma20 = gBeamCut->IsInside(beamInfo->fBeamAoQ, beamInfo->fBeamZ) ;

    event = ievt;

    STVertex *vertex = (STVertex*)vertexArray->At(0);
    if( !vertex ) continue;
    TVector3 vert  = vertex->GetPos();    

    STVertex *vertexVA = (STVertex*)vertexVAArray->At(0);
    if( !vertexVA ) continue;
    TVector3 vertVA  = vertexVA->GetPos();    

    STVertex* bdcvertex = (STVertex*)vertexBDCArray->At(0);
    if( !bdcvertex ) continue;
    TVector3 bdcVert = bdcvertex->GetPos();

    Bool_t vtxZ3Sigma  = TMath::Abs(vert.Z() - vzPar[1]) <= vzPar[2]*3.;
    Bool_t vtxbdcX     = TMath::Abs(vert.X() - bdcVert.X() - vxbxPar[1]) <= vxbxPar[2]*3.;
    Bool_t vtxbdcY     = TMath::Abs(vert.Y() - bdcVert.Y() - vybyPar[1]) <= vybyPar[2]*3.;
    Bool_t isGoodEvent = sigma20*vtxZ3Sigma*vtxbdcX*vtxbdcY;

    if( !isGoodEvent ) continue;

    TIter nextReco(trackArray);
    TIter nextVA(trackVAArray);

    STRecoTrack* atrackReco= NULL;
    STRecoTrack* atrackVA = NULL;


    mtrack0 = trackArray->GetEntries();
    mtrack1 = 0;

    while( (atrackReco = (STRecoTrack*)nextReco() ) ) {

      TVector3 vdist = atrackReco->GetPOCAVertex() - vert;
      auto dist = vdist.Mag();
      if( dist <= 20 ) 
	mtrack1++;
    }

    for( auto im : nmcut ) {
      if( mtrack1 >= mtcut[systemID][im][0] && mtrack1 < mtcut[systemID][im][1] )  {
	hmultsel[im] -> Fill( mtrack1 );
	Msel[im] = 1;
      }
      else 
	Msel[im] = 0;
    }


    for( auto it : ROOT::TSeqL( mtrack0 ) ) {

      atrackReco = (STRecoTrack*)trackArray->At(it);
      auto helixID = atrackReco->GetHelixID();

      if( !(atrackVA   = (STRecoTrack*)trackVAArray->At(it) ) ) continue;
      auto helixIDVA = atrackVA->GetHelixID();
      if( helixIDVA == -1 ) continue;
      
      Int_t diffID = helixIDVA - helixID;
      if( diffID > 0  ) {
	while( 1 ) {
	  it++;
	  if( !(atrackReco = (STRecoTrack*)trackArray->At(it) ) ) break;
	  helixID = atrackReco->GetHelixID(); 

	  diffID = helixIDVA - helixID;
	  if( diffID == 0 ) break;
	}
      }

      if( bDEBUG ) {
	cout << " ievet " << event
	     << " reco " << atrackReco->GetHelixID()
	     << " VA "   << atrackVA  ->GetHelixID()
	     << " dedx " << atrackVA->GetdEdxWithCut(0,0.7)
	     << endl;
      }

      apart->SetRecoTrack(atrackReco, atrackVA);
      apart->SetVertex(vert);

      //      TVector3 POCA = atrackReco->GetPOCAVertex();
      TVector3 POCA = apart->GetPOCAVertex();

      Bool_t bgtrack = 
	(abs(POCA.Z() + 14.85 ) <= 3.* 1.33) &&
	(abs(POCA.X()) <= 15.) &&
	(abs(POCA.Y() + 205.) <= 20. );
      
      if( !bgtrack ) apart->SetVertexAtTargetFlag(0);
      
	
      auto dEdx   = apart->GetdEdx();
      auto Vmom   = apart->GetMomentumAtTarget();
      auto P      = Vmom.Mag();
      auto phi    = Vmom.Phi()*TMath::RadToDeg();
      auto massH  = massCalH ->CalcMass(1., Vmom, dEdx, kTRUE );
      auto massHe = massCalHe->CalcMass(2., Vmom, dEdx, kTRUE );
      // auto vaCharge = atrackVA -> GetVAGFCharge();
      // auto vaR  = Vmom.Mag()/vaCharge;

      //      cout << " massH " << massH << " massHe " << massHe << endl;

      apart->SetBBMass(massH);
      apart->SetBBMassHe(massHe);
      auto dist = apart->GetDistanceAtVertex();
      //      dist = (POCA - vert).Mag();
      auto distv= (POCA - vert).Mag();
      auto ncl  = apart->GetNumCluster();
      auto nclv = (atrackVA   -> GetClusterIDArray())->size() - 1;

      if( ncl != nclv )
	cout  << " ERROR::  events " << ievt << " partID " << mtrack2 << " ncl " << ncl << " nclv " << nclv 
	      << " helixID " << helixID << " helixIDVA " << helixIDVA
	      << endl;
	
      if( dist > 20 )
	apart->SetDistanceAtVertexFlag(0);
      if( nclv < 15 )
	apart->SetNumClusterFlag(0);
      

      UInt_t nPID[2] = {0,0};
      nPID[0] =  GetPIDFit(1, massH,  Vmom, mtrack1);  
      nPID[1] =  GetPIDFit(2, massHe, Vmom, mtrack1);


      if( 0 ) {
	cout << " next---> " << event
	     << " reco " << atrackReco->GetHelixID()
	     << " dedx " << atrackVA->GetdEdxWithCut(0,0.7)
	     << " p " << P
	     << " massH " << massH
	     << " massHe " << massHe
	     << " npid[0] " << nPID[0]
	     << " npid[1] " << nPID[1]
	     << endl;
      }


      apart->RotateAlongBeamDirection(beamInfo->fRotationAngleATargetPlane/1000.,
				      beamInfo->fRotationAngleBTargetPlane/1000.);


      Int_t pid = 0;
      if( nPID[1] != 0 ) 
	pid = nPID[1];
      else if( nPID[0] != 0 ) 
	pid = nPID[0];

      apart->SetPID(pid);
      apart->SetBBMass(massH);
      apart->SetBBMassHe(massHe);

      Int_t ipart = apart->GetPID_seq();
      if( ipart >= 2 && ipart < 8) {
	ipart -= 2;
	if( apart->GetP() > momCut[ipart] ) apart->SetMomentumFlag(0);

	auto y  = apart->GetRapidity()/yNN[systemID] - 1.;
	auto pt = apart->GetRotatedMomentum().Perp();
	
	for( auto im : nmcut ) {
	  if( Msel[im] )
	    hypt[ipart][im] -> Fill(y, pt/1000.);
	}
	if( ipart <= 2 )
	  hmassHp[ipart]  -> Fill( P, massH  );
	else if( ipart <= 4 )
	  hmassHep[ipart] -> Fill( P, massHe );
      }

      //@@@@@fill-----
      mtrack2 = ParticleArray->GetEntries();
      apart->SetTrackID( mtrack2 );
      if( apart->GetGoodTrackFlag() > 1000)
	new( (*ParticleArray)[mtrack2] ) STParticle(*apart);
    }
    
    outtree->Fill();

    if( mtrack1 > 0 )
      hmult -> Fill( mtrack1 );
    


    //    cout << endl;
  }
    
  for( auto ipart : sqPart ) {
    outfile[ipart]->cd();
    hmult->SetDirectory(outfile[ipart]);
    hmult->Write();
    for( auto im : nmcut ) {
      hmultsel[im]->SetDirectory(outfile[ipart]);
      hmultsel[im]->Write();
      hypt[ipart][im]->Write();
    }
  }

  outfile[5]->cd();
  hmult->SetDirectory(outfile[5]);    
  for( auto im : nmcut ) {
    hmultsel[im]->SetDirectory(outfile[5]);
    hmultsel[im]->Write();
  }
  outtree->Write();
      

  //@@@@@Draw ++++++++++++++++++++++++++++++++++++++++++++++++++
  TCanvas *ccv;
  TPad *ppad;
  Int_t iid = 0;
  Int_t id;

  if( 1 ) {
    
    ccv = new TCanvas(Form("ccv%d",iid),Form("ccv%d",iid),1200,700); iid++;
    ccv->Divide(5,2);  
    id = 1;
    for( auto ii : {0,1,2} ) {
      ppad = (TPad*)ccv->cd(id); id++;
      ppad->SetLogz();
      hmassHp[ii]->Draw("colz");
    }
    for( auto ii : {3,4} ) {
      ppad = (TPad*)ccv->cd(id); id++;
      ppad->SetLogz();
      hmassHep[ii]->Draw("colz");
    }

    for( auto ii : ROOT::TSeqI(5) ) {
      ppad= (TPad*)ccv->cd(id); id++; 
      ppad->SetLogz();
      hypt[ii][0]->Draw("colz");
    }
  }

  
  cout << " Saved in " << outfile[5]->GetName() << endl;
}

