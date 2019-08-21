#include "myAnalysis.h"
////--------------------------------------------------&&
void makeLCPTree();
void fillLCPSpectra();
void drawLCPSpectra();
void drawLCPRatio();
void makeLCPSpectra()
{
	//beamID.clear(); beamID.push_back(2);beamID.push_back(3);
	//makeLCPTree();  
   fillLCPSpectra();
	drawLCPSpectra();
	drawLCPRatio();
}

const bool drawAMD=true;
const bool drawJQMD=false;

double mLimit[][2] = {{700.,1200.}, {1400.,2300}, {2400.,3200.}, {2300.,3000.}, {3000.,4200.}, {5000.,6000.}};

//void makeLCPTree()
void fillLCPSpectra()
{
	auto massTreeFile = TFile::Open("rootfiles/MassTree.root");
	TTree* massTree[nBeam];
//	TTree* eveTree[nBeam];
	TH1* h1Mult[nBeam][3];
	TVector3 *P=nullptr, *rotP=nullptr;
	Double_t R,dEdx,dist;
	Int_t ncl;
	Double_t pitch,yaw;
	Double_t theta,phi;
	Double_t massP[2];
	Double_t massD[2];
	Double_t calibdEdx[2];
	for(auto beam: beamID){
		massTreeFile->GetObject(Form("massTree_%dSn",beamA[beam]),massTree[beam]);
		massTree[beam]->SetBranchStatus("*",0);
		massTree[beam]->SetBranchStatus("run",1);     massTree[beam]->SetBranchAddress("run",&run);
		massTree[beam]->SetBranchStatus("event",1);     massTree[beam]->SetBranchAddress("event",&event);
		massTree[beam]->SetBranchStatus("P",1);     massTree[beam]->SetBranchAddress("P",&P);
		massTree[beam]->SetBranchStatus("rotP",1);  massTree[beam]->SetBranchAddress("rotP",&rotP);
		massTree[beam]->SetBranchStatus("R",1);     massTree[beam]->SetBranchAddress("R",&R);
		massTree[beam]->SetBranchStatus("dEdx",1);  massTree[beam]->SetBranchAddress("dEdx",&dEdx);
		massTree[beam]->SetBranchStatus("dist",1);  massTree[beam]->SetBranchAddress("dist",&dist);
		massTree[beam]->SetBranchStatus("ncl",1);   massTree[beam]->SetBranchAddress("ncl",&ncl);
		massTree[beam]->SetBranchStatus("pitch",1);   massTree[beam]->SetBranchAddress("pitch",&pitch);
		massTree[beam]->SetBranchStatus("yaw",1);   massTree[beam]->SetBranchAddress("yaw",&yaw);
		massTree[beam]->SetBranchStatus("phi",1);   massTree[beam]->SetBranchAddress("phi",&phi);
		massTree[beam]->SetBranchStatus("massP",1); massTree[beam]->SetBranchAddress("massP",massP);
		massTree[beam]->SetBranchStatus("massD",1); massTree[beam]->SetBranchAddress("massD",massD);
		massTree[beam]->SetBranchStatus("calibdEdx",1); massTree[beam]->SetBranchAddress("calibdEdx",calibdEdx);
		massTree[beam]->SetBranchStatus("beamY0",1); massTree[beam]->SetBranchAddress("beamY0",&beamY0);
		massTree[beam]->SetBranchStatus("boostBeta",1); massTree[beam]->SetBranchAddress("boostBeta",&boostBeta);
		massTree[beam]->SetBranchStatus("multTPC",1); massTree[beam]->SetBranchAddress("multTPC",multTPC);
		
		//massTreeFile->GetObject(Form("EventTree_%dSn",beamA[beam]),eveTree[beam]);
		for(auto mul: ROOT::TSeqI(3))
			massTreeFile->GetObject(Form("h1Mult_%dSn_%d",beamA[beam],mul),h1Mult[beam][mul]);
		
	}


	auto massGateFile = TFile::Open("rootfiles/MassGate.root");
	TF1* f1TotalBase[nBeam][6][nRBinBase];
	TF1* f1Total[nBeam][6][nRBin];
	for(auto beam: beamID)for(auto pid: ROOT::TSeqI(6))for(auto rbin: ROOT::TSeqI(nRBin)){
		Int_t pd = pid==0?0:1;
		TString name = Form("_%dSn_"+pidName[pid]+"_%s_%d",beamA[beam],pd==0?"P":"D",rbin);
		massGateFile->GetObject("f1Total"+name,f1Total[beam][pid][rbin]);
		if(rbin<nRBinBase) massGateFile->GetObject("f1TotalBase"+name,f1TotalBase[beam][pid][rbin]);
	}
	
	double fac[6]={3.,3.,2.5,2.,2.,1.5};
	
	auto getXPoint = [&](TF1* ftot, int pid, int i)
	{
		TF1 fsig("fsig",PsdVoigt,mRange[pid][0],mRange[pid][1],4);
		TF1 fbg1("fbg1",PsdVoigt,mRange[pid][0],mRange[pid][1],4);
		TF1 fbg2("fbg2",PsdVoigt,mRange[pid][0],mRange[pid][1],4);
		for(auto par: ROOT::TSeqI(4)){
			fsig.SetParameter(par,ftot->GetParameter(par));
			fbg1.SetParameter(par,ftot->GetParameter(par+4));
			fbg2.SetParameter(par,ftot->GetParameter(par+8));
		}
		double dx = i==0?-0.1:0.1;
		double start = fsig.GetParameter(1);
		double df=1.e10;
		double max = TMath::Max(fsig.GetParameter(1)-(fac[pid]+0.5)*fsig.GetParameter(2),mRange[pid][0]);
		double min = TMath::Min(fsig.GetParameter(1)+(fac[pid]+0.5)*fsig.GetParameter(2),mRange[pid][1]);
		while(df>0){
			if(i==0&&start<max)break;
			if(i==1&&start>min)break;
			df = i==0 ? fsig.Eval(start)-fbg1.Eval(start) : fsig.Eval(start)-fbg2.Eval(start);
			start += dx;
		}
		return start;
	};

	double baseXPoints[nBeam][6][nRBinBase][2]={};
	double xPoints[nBeam][6][nRBin][2]={};
	double normalWidth[nBeam][6][nRBin][2]={};
	double tightWidth[nBeam][6][nRBin][2]={};
	for(auto beam: beamID)for(auto pid: ROOT::TSeqI(6))for(auto rbin: ROOT::TSeqI(nRBin))for(auto i: ROOT::TSeqI(2)){
		if(rbin<nRBinBase)baseXPoints[beam][pid][rbin][i] = getXPoint(f1TotalBase[beam][pid][rbin],pid,i);
		xPoints[beam][pid][rbin][i] = getXPoint(f1Total[beam][pid][rbin],pid,i);
		double meanM = f1Total[beam][pid][rbin]->GetParameter(1);
		double sigM  = f1Total[beam][pid][rbin]->GetParameter(2);
		normalWidth[beam][pid][rbin][i] = i==0?meanM-fac[pid]*sigM:meanM+fac[pid]*sigM;
		tightWidth[beam][pid][rbin][i]  = i==0?meanM-(fac[pid]-1.0)*sigM:meanM+(fac[pid]+1.0)*sigM;
	};
	
	auto checkLoose = [&](double m, int beam, int pid, int rbin)
	{
		bool isLoose=false;
		double min = TMath::Max(mLimit[pid][0],baseXPoints[beam][pid][int(rbin/(nRBin/nRBinBase))][0]);
		double max = TMath::Min(mLimit[pid][1],baseXPoints[beam][pid][int(rbin/(nRBin/nRBinBase))][1]);
		if(m>=min&&m<=max) isLoose=true;
		return isLoose;
	};
	auto checkNormal = [&](double m, int beam, int pid, int rbin)
	{
		bool isNormal = false;
		if(m>=normalWidth[beam][pid][rbin][0]&&m<=normalWidth[beam][pid][rbin][1]) isNormal=true;
		return isNormal;
	};
	auto checkTight  = [&](double m, int beam, int pid, int rbin)
	{
		bool isTight = false;
		if(m>=tightWidth[beam][pid][rbin][0]&&m<=tightWidth[beam][pid][rbin][1]) isTight=true;
		return isTight;
	};
	
	TF1* fsig[6];
	for(auto pid: ROOT::TSeqI(6)){
		fsig[pid] = new TF1(Form("fsig%d",pid),PsdVoigt,mRange[pid][0],mRange[pid][1],4);
		fsig[pid]->SetNpx(2000);
	}
	auto getSignalToTotalFactor = [&](double m, int beam, int pid, int rbin)
	{
		for(auto par: ROOT::TSeqI(4))
			fsig[pid]->SetParameter(par,f1Total[beam][pid][rbin]->GetParameter(par));
		double frac = fsig[pid]->Eval(m)/f1Total[beam][pid][rbin]->Eval(m);
		// when there are no background fit component, signal/total becomes a little bit larger than 1. like 1.0001 .
		if(frac<0.) frac=0.;			
		else if(frac>1.) frac=1.;
		return frac;
	};
	auto getSignalToTotalFactorBase = [&](double m, int beam, int pid, int rbin)
	{
		for(auto par: ROOT::TSeqI(4))
			fsig[pid]->SetParameter(par,f1TotalBase[beam][pid][int(rbin/(nRBin/nRBinBase))]->GetParameter(par));
		double frac = fsig[pid]->Eval(m)/f1TotalBase[beam][pid][int(rbin/(nRBin/nRBinBase))]->Eval(m);
		// when there are no background fit component, signal/total becomes a little bit larger than 1. like 1.0001 .
		if(frac<0.) frac=0.;			
		else if(frac>1.) frac=1.;
		return frac;
	};

	
	
	auto outFile = TFile::Open("rootfiles/LCPSpectra.root","recreate");
	TTree *testTree[nBeam];
	TTree *outTree[nBeam];

	TVector3 oP, orotP;
	Double_t centrality;
	Bool_t   isInsideLoose[6];	 // x-points
	Bool_t   isInsideNormal[6]; // mean pm sigma
	Bool_t   isInsideTight[6];  // mean pm small sigma
	Double_t fracSigToTot;     // for estimation of background by other particle.
	Double_t fracSigToTotBase;     // for estimation of background by other particle.
	Double_t z,a;
	Double_t y0,pt,ke;
	
	auto setTTreeBranch = [&](TTree* t){
		t->Branch("run",&run,"run/I");
		t->Branch("event",&event,"event/I");
		t->Branch("P",&oP);
		t->Branch("rotP",&orotP);
		t->Branch("R",&R,"R/D");
		t->Branch("dEdx",&dEdx,"dEdx/D");
		t->Branch("dist",&dist,"dist/D");
		t->Branch("ncl",&ncl,"ncl/I");
		t->Branch("theta",&theta,"theta/D");
		t->Branch("phi",&phi,"phi/D");
		t->Branch("pitch",&pitch,"pitch/D");
		t->Branch("yaw",&yaw,"yaw/D");
		t->Branch("massP",massP,"massP[2]/D");
		t->Branch("massD",massD,"massD[2]/D");
		t->Branch("calibdEdx",calibdEdx,"calibdEdx[2]/D");
		t->Branch("beamY0",&beamY0,"beamY0/D");
		t->Branch("boostBeta",&boostBeta,"boostBeta/D");
		t->Branch("multTPC",multTPC,"multTPC[3]/I");
		t->Branch("isInsideLoose",isInsideLoose,"isInsideLoose[6]/O");
		t->Branch("isInsideNormal",isInsideNormal,"isInsideNormal[6]/O");
		t->Branch("isInsideTight",isInsideTight,"isInsideTight[6]/O");
		t->Branch("fracSigToTot",&fracSigToTot,"fracSigToTot/D");
		t->Branch("fracSigToTotBase",&fracSigToTotBase,"fracSigToTotBase/D");
		t->Branch("centrality",&centrality,"centrality/D");
		t->Branch("z",&z,"z/D");
		t->Branch("a",&a,"a/D");
		t->Branch("y0",&y0,"y0/D");
		t->Branch("pt",&pt,"pt/D");
		t->Branch("ke",&ke,"ke/D");
	};
	for(auto beam: beamID){
		testTree[beam] = new TTree(Form("LCPTreeTest_%dSn",beamA[beam]),"");
		setTTreeBranch(testTree[beam]);
		outTree[beam] = new TTree(Form("LCPTree_%dSn",beamA[beam]),"");
		setTTreeBranch(outTree[beam]);
	}
	
	TH2* h2PtY[nBeam][6][3];
	TH1* h1Y[nBeam][6][3];
	TH1* h1Pt[nBeam][6][5][3];
	TH1* h1Ecm[nBeam][6][9][3];
	
	auto setupTH1 = [&](TH1* h){ 
		h->SetMarkerStyle(21); h->SetMarkerColor(kBlue+3); h->SetMarkerSize(1.2);
		h->SetLineColor(kBlue-5); 
		h->Sumw2();
	};

	for(auto beam: beamID)for(auto pid: ROOT::TSeqI(6))for(auto type: ROOT::TSeqI(3)){
		TString name      = Form("_%dSn_%s",beamA[beam],pidName[pid].Data());
		TString typeTitle = type==0?"Total":type==1?"Signal":"Background";
		TString title     = Form("%dSn %s ",beamA[beam],pidName[pid].Data());
		h2PtY[beam][pid][type] = DefineTH2D(Form("h2PtY"+name+"_%d",type),title+"pT-y "+typeTitle+"; y^{0} = y^{c.m.}/y^{beam};pT (GeV/c)",80,-1.5,1.5,80,0,ptRange[pid]);
		h1Y[beam][pid][type]   = DefineTH1D(Form("h1Y"+name+"_%d",type),title+"rapidity "+typeTitle+";y^{0} = y^{c.m.}/y^{beam};dN/dy^{0}",20,-1.5,1.5);
		setupTH1(h1Y[beam][pid][type]);
		for(auto ybin: ROOT::TSeqI(5)){
			auto ptTitle = Form(title+"pT "+typeTitle+", %.1f #leqy^{0}< %.1f;pT (GeV/c); dN/dpT",yBin[ybin],yBin[ybin+1]);
			h1Pt[beam][pid][ybin][type] = DefineTH1D(Form("h1Pt"+name+"_%d_%d",ybin,type),ptTitle,20,0,ptRange[pid]);
			setupTH1(h1Pt[beam][pid][ybin][type]);
		}
		for(auto tbin: ROOT::TSeqI(9)){
			auto ecmTitle = Form(title+"Ecm "+typeTitle+", %d #leq#theta< %d;E_{c.m.}/A (MeV); dN/dE_{c.m.}",tbin*20,(tbin+1)*20);
			h1Ecm[beam][pid][tbin][type] = DefineTH1D(Form("h1Ecm"+name+"_%d_%d",tbin,type),ecmTitle,24,0,240);
			setupTH1(h1Ecm[beam][pid][tbin][type]);
		}
	}

	
	Double_t centralities[nBeam][100]={};
	for(auto beam: beamID)for(auto mbin: ROOT::TSeqI(100))
		centralities[beam][mbin] = trigEff[beam]*calcCentrality(h1Mult[beam][1],mbin);
	// temporal
	int binCentralCut[nBeam]={};
	for(auto beam: beamID)for(auto ibin: ROOT::TSeqI(1,h1Mult[beam][1]->GetNbinsX(),1))
		if(trigEff[beam]*h1Mult[beam][1]->Integral(ibin,h1Mult[beam][1]->GetNbinsX())/h1Mult[beam][1]->Integral()<centralCut){
			binCentralCut[beam]=ibin;
			std::cout<<"beam:"<<beamA[beam]<<", cut mult:"<<ibin<<std::endl;
			break;
		}

	auto isInPIDMassRange = [&](double m, int pid)
	{ return m>=pidMRange[pid][0] && m<=pidMRange[pid][1]; };

	
	for(auto beam: beamID){
		auto entries = massTree[beam]->GetEntries();
		Int_t nTestTrack=0;
		Int_t oldRun=-1;
		for(auto i: ROOT::MakeSeq(entries)){
			if(i%1000000==0) std::cout<<"Track entry: "<<i<<"/"<<entries<<std::endl;
			massTree[beam]->GetEntry(i);

			centrality = centralities[beam][multTPC[1]];

			oP = *P;
			orotP = *rotP;
				
			for(auto pid: ROOT::TSeqI(6)){
				isInsideLoose[pid]=false;
				isInsideNormal[pid]=false;
				isInsideTight[pid]=false;
			}
			fracSigToTot = -1.;
			fracSigToTotBase = -1.;
			y0=-99.;
			pt=-1.;
			ke=-1.;
			z=-1; a=-1;

			int pid=-1;
			if(isInPIDMassRange(massP[0],0)) pid=0;
			else if(isInPIDMassRange(massD[0],1)) pid=1;
			else if(isInPIDMassRange(massD[0],2)&&massD[0]<=0.5*R+2800.) pid=2;
			else if(isInPIDMassRange(massD[1],3)) pid=3;
			else if(isInPIDMassRange(massD[1],4)) pid=4;
			else if(isInPIDMassRange(massD[1],5)) pid=5;

			if(pid!=-1){
				double m = pid==0?massP[0]: pid<=2?massD[0]:massD[1];
				double wRBin = (rRange[pid][1]-rRange[pid][0])/nRBin;
				int rbin = (R-rRange[pid][0])/wRBin;
				if(rbin>nRBin-1) rbin=nRBin-1;
				if(rbin>=0&&rbin<=nRBin-1){
					isInsideLoose[pid]=true;
					isInsideNormal[pid] = checkNormal(m,beam,pid,rbin);
					isInsideTight[pid]  = checkTight(m,beam,pid,rbin);
					fracSigToTot        = getSignalToTotalFactor(m,beam,pid,rbin);
					fracSigToTotBase    = getSignalToTotalFactorBase(m,beam,pid,rbin);
				}

				// fill spectra
				if(isInsideLoose[pid]&&dEdx<=1000.&&isPhiFiducial(phi)&&centrality<=centralCut&&isPYFiducial(pitch,yaw)){
					z = pid<=2 ? 1: 2;
					a = pid<=2 ? pid+1: pid!=5? pid: 6;
					TVector3 mom3(rotP->X()*z, rotP->Y()*z, rotP->Z()*z);
					TLorentzVector mom(mom3,TMath::Sqrt(mom3.Mag2()+pidMass[pid]*pidMass[pid]));
					mom.Boost(TVector3(0.,0.,boostBeta));

					y0 = mom.Rapidity()/beamY0;
					pt = mom.Pt()/1000.;
					ke = mom.E()-mom.M();

					for(auto type: ROOT::TSeqI(3)){
						double frac = type==0?1.:type==1?fracSigToTot:1.-fracSigToTot;

						h2PtY[beam][pid][type]->Fill(y0,pt,frac);
						h1Y[beam][pid][type]->Fill(y0,frac);
						Int_t ybin = -1;
						if(y0<-1.2)ybin=0;
						else if(y0>1.2)ybin=4;
						else for(auto yb: ROOT::TSeqI(5))if(y0>=yBin[yb]&&y0<yBin[yb+1])ybin=yb;
						if(ybin!=-1) h1Pt[beam][pid][ybin][type]->Fill(pt,frac);
						Double_t theta = mom.Theta()*TMath::RadToDeg();
						Int_t tbin = -1;
						for(auto tb: ROOT::TSeqI(9))if(theta>=tb*20&&theta<(tb+1)*20)tbin=tb;
						if(tbin!=-1) h1Ecm[beam][pid][tbin][type]->Fill(ke/a,frac);
					}
				}
			}
			//testTree[beam]->Fill(); nTestTrack++;
			if(nTestTrack<10000){ testTree[beam]->Fill(); nTestTrack++; }
			if(oldRun!=run) nTestTrack=0; 
			oldRun = run;
		//	outTree[beam]->Fill();
			
		}
	}
	
	double phiScale = 360./phiCut/2.;
	auto dndx = [&](TH1* h, long e)
	{
		auto dx = (h->GetXaxis()->GetXmax()-h->GetXaxis()->GetXmin())/h->GetNbinsX();
		h->Scale(1./dx/e*phiScale);
	};
	for(auto beam: beamID)for(auto pid: ROOT::TSeqI(6))for(auto type: ROOT::TSeqI(3)){
		dndx(h1Y[beam][pid][type],h1Mult[beam][1]->Integral(binCentralCut[beam],h1Mult[beam][1]->GetNbinsX()));
		for(auto ybin: ROOT::TSeqI(5)){
			dndx(h1Pt[beam][pid][ybin][type],h1Mult[beam][1]->Integral(binCentralCut[beam],h1Mult[beam][1]->GetNbinsX()));
			h1Pt[beam][pid][ybin][type]->Scale(1./abs(yBin[ybin+1]-yBin[ybin]));
		}
		for(auto tbin: ROOT::TSeqI(9))
			dndx(h1Ecm[beam][pid][tbin][type],h1Mult[beam][1]->Integral(binCentralCut[beam],h1Mult[beam][1]->GetNbinsX()));
	}



	outFile->cd();
	for(auto beam: beamID){
		testTree[beam]->Write();
//		outTree[beam]->Write();

	}
	for(auto beam: beamID)for(auto pid: ROOT::TSeqI(6))for(auto type: ROOT::TSeqI(3)){
		h2PtY[beam][pid][type]->Write();
		h1Y[beam][pid][type]->Write();
		for(auto ybin: ROOT::TSeqI(5)) h1Pt[beam][pid][ybin][type]->Write();
		for(auto tbin: ROOT::TSeqI(9)) h1Ecm[beam][pid][tbin][type]->Write();
	}
	outFile->Close();

}
	



Double_t pidRange[][2] = {{0.45,1.2},{0.7,1.4},{1.05,1.8},{0.4,1.3},{0.6,1.4},{1.25,2.45}};
TString types[3]       = {"Total","Signal","Background"};

///--------------------------------------------------%%

void drawLCPSpectra()
{

	//beamID.clear();
	//beamID.push_back(0);
	//beamID.push_back(3);
	auto f = TFile::Open("rootfiles/LCPSpectra.root");
	TH2* h2PtY[nBeam][6][3];
	TH1* h1Y[nBeam][6][3];
	TH1* h1Pt[nBeam][6][5][3];
	TH1* h1Ecm[nBeam][6][9][3];
	for(auto beam: beamID)for(auto pid: ROOT::TSeqI(6))for(auto type: ROOT::TSeqI(3)){
		f->GetObject(Form("h2PtY_%dSn_%s_%d",beamA[beam],pidName[pid].Data(),type),h2PtY[beam][pid][type]);
		f->GetObject(Form("h1Y_%dSn_%s_%d",beamA[beam],pidName[pid].Data(),type),h1Y[beam][pid][type]);
		for(auto ybin: ROOT::TSeqI(5))
			f->GetObject(Form("h1Pt_%dSn_%s_%d_%d",beamA[beam],pidName[pid].Data(),ybin,type),h1Pt[beam][pid][ybin][type]);
		for(auto tbin: ROOT::TSeqI(9))
			f->GetObject(Form("h1Ecm_%dSn_%s_%d_%d",beamA[beam],pidName[pid].Data(),tbin,type),h1Ecm[beam][pid][tbin][type]);
	}


	TCanvas *cvsPtY[nBeam][3];
	for(auto beam: beamID)for(auto type: ROOT::TSeqI(3)){
		cvsPtY[beam][type] = new TCanvas(Form("cvsPtY_%dSn_%d",beamA[beam],type),"",1600,1000);
		cvsPtY[beam][type]->Divide(3,2,0.006,0.006);
		for(auto pid: ROOT::TSeqI(6)){
			//cvsPtY[beam][type]->cd(pid+1)->SetLogz();
			Polish(cvsPtY[beam][type]->cd(pid+1));
			h2PtY[beam][pid][type]->Draw("colz");
		}
		cvsPtY[beam][type]->SaveAs(Form("fig/PtY_%dSn_%s.png",beamA[beam],types[type].Data()));
	}
	
	
	THStack *hsY[6][3], *hsPt[6][5][3], *hsEcm[6][9][3];
	TString systems[4]={}; int csystems=0;
	for(auto beam: beamID){
		systems[csystems]=Form("%dSn%dSn",beamA[beam],targetA[beam]);
		csystems++;
	}
	auto setupTHStack = [&](THStack* hs){
		for(auto id: ROOT::TSeqI(hs->GetNhists())){
		  ((TH1D*)hs->GetHists()->At(id))->SetMarkerColor(GetColor(id+1));
		  ((TH1D*)hs->GetHists()->At(id))->SetTitle(systems[id]);
		}
	};
	for(auto pid: ROOT::TSeqI(6))for(auto type: ROOT::TSeqI(3)){
		TString name   = pidName[pid]+"_"+types[type];
		TString yTitle = pidName[pid]+" rapidity ("+types[type]+");y^{0} = y^{cm}/y^{beam};dN/dy^{0}";
		hsY[pid][type] = new THStack("hsY_"+name,yTitle);
		for(auto beam: beamID) hsY[pid][type]->Add(h1Y[beam][pid][type]);
		setupTHStack(hsY[pid][type]);

		for(auto ybin: ROOT::TSeqI(5)){
		   auto ptTitle = Form(pidName[pid]+" pT ("+types[type]+"), %.1f #leq y^{0}< %.1f;p_{T} (GeV/c); dN/dp_{T}",yBin[ybin],yBin[ybin+1]);
			hsPt[pid][ybin][type] = new THStack(Form("hsPt_"+name+"_%d",ybin),ptTitle);
			for(auto beam: beamID) hsPt[pid][ybin][type]->Add(h1Pt[beam][pid][ybin][type]);
			setupTHStack(hsPt[pid][ybin][type]);
		}
		for(auto tbin: ROOT::TSeqI(9)){
			auto ecmTitle = Form(pidName[pid]+" Ecm ("+types[type]+"), %d #leq#theta< %d;E_{c.m.}/A (MeV); dN/dEcm",tbin*20,(tbin+1)*20);
			hsEcm[pid][tbin][type] = new THStack(Form("hsEcm_"+name+"_%d",tbin),ecmTitle);
			for(auto beam: beamID) hsEcm[pid][tbin][type]->Add(h1Ecm[beam][pid][tbin][type]);
			setupTHStack(hsEcm[pid][tbin][type]);
		}
	}

	
	auto setupPad    = [&](TVirtualPad* c){ c->SetLogy(); c->SetGrid(); Polish(c); };
	auto setupLegend = [&](TLegend* l){ l->SetBorderSize(0); l->SetFillStyle(3001); };

//	size_t nTypes = sizeof(types)/sizeof(types[0]);
	size_t nTypes = 2;
	TCanvas *cvsY[3];
	for(auto type: ROOT::TSeqI(nTypes)){
		cvsY[type] = new TCanvas(Form("cvsY_%d",type),"",1600,1000);
		cvsY[type]->Divide(3,2);
		for(auto pid: ROOT::TSeqI(6)){
			Polish(cvsY[type]->cd(pid+1));
			cvsY[type]->cd(pid+1)->SetGrid();
			hsY[pid][type]->Draw("nostack ep");
			Polish(hsY[pid][type]);
			setupLegend(cvsY[type]->cd(pid+1)->BuildLegend(0.15,0.68,0.55,0.88));
		}
		cvsY[type]->SaveAs("fig/Y_"+types[type]+".png");
	}
	
	TCanvas *cvsPt[6][3];
	for(auto pid: ROOT::TSeqI(6))for(auto type: ROOT::TSeqI(nTypes)){
		cvsPt[pid][type] = new TCanvas(Form("cvsPt_%d_%d",pid,type),"",1600,1000);
		cvsPt[pid][type]->Divide(3,2);
		for(auto ybin: ROOT::TSeqI(5)){
			setupPad(cvsPt[pid][type]->cd(ybin+1));
			hsPt[pid][ybin][type]->Draw("nostack ep");
			Polish(hsPt[pid][ybin][type]);
			setupLegend(cvsPt[pid][type]->cd(ybin+1)->BuildLegend(0.15,0.15,0.55,0.35));
		}
		cvsPt[pid][type]->SaveAs("fig/Pt_"+pidName[pid]+"_"+types[type]+".png");
	}
	
	TCanvas *cvsEcm[6][3];
	for(auto pid: ROOT::TSeqI(6))for(auto type: ROOT::TSeqI(nTypes)){
		cvsEcm[pid][type] = new TCanvas(Form("cvsEcm_%d_%d",pid,type),"",1600,1200);
		cvsEcm[pid][type]->Divide(3,3);
		for(auto tbin: ROOT::TSeqI(9)){
			setupPad(cvsEcm[pid][type]->cd(tbin+1));
			hsEcm[pid][tbin][type]->Draw("nostack ep");
			Polish(hsEcm[pid][tbin][type]);
			setupLegend(cvsEcm[pid][type]->cd(tbin+1)->BuildLegend(0.15,0.15,0.55,0.35));
		}
		cvsEcm[pid][type]->SaveAs("fig/Ecm_"+pidName[pid]+"_"+types[type]+".png");
	}


	f->Close();
}

void drawLCPRatio()
{
	//auto f = TFile::Open("rootfiles/LCPSpectra_phi60_central20_rbin30.root");
	auto f = TFile::Open("rootfiles/LCPSpectra.root");
	//auto f = TFile::Open("rootfiles/LCPSpectra_pitch-24-24_yaw0-40_central20_rbin30.root");
	TH1* h1Y[nBeam][6][3];
	TH1* h1Pt[nBeam][6][5][3];
	TH1* h1Ecm[nBeam][6][9][3];
	for(auto beam: beamID)for(auto pid: ROOT::TSeqI(6))for(auto type: ROOT::TSeqI(3)){
		f->GetObject(Form("h1Y_%dSn_%s_%d",beamA[beam],pidName[pid].Data(),type),h1Y[beam][pid][type]);
		for(auto ybin: ROOT::TSeqI(5))
			f->GetObject(Form("h1Pt_%dSn_%s_%d_%d",beamA[beam],pidName[pid].Data(),ybin,type),h1Pt[beam][pid][ybin][type]);
		for(auto tbin: ROOT::TSeqI(9))
			f->GetObject(Form("h1Ecm_%dSn_%s_%d_%d",beamA[beam],pidName[pid].Data(),tbin,type),h1Ecm[beam][pid][tbin][type]);
	}

	TH1* h1YR[6][3];
	TH1* h1PtR[6][5][3];
	TH1* h1EcmR[6][9][3];
	
	for(auto pid: ROOT::TSeqI(6))for(auto type: ROOT::TSeqI(2)){
		h1YR[pid][type]=(TH1D*)h1Y[0][pid][type]->Clone(); 
		h1YR[pid][type]->Divide(h1Y[3][pid][type]);
		h1YR[pid][type]->SetNameTitle(Form("h1YR_%d_%d",pid,type),"Data 132Sn/108Sn");
		for(auto ybin: ROOT::TSeqI(5)){
			h1PtR[pid][ybin][type]=(TH1D*)h1Pt[0][pid][ybin][type]->Clone(); 
			h1PtR[pid][ybin][type]->Divide(h1Pt[3][pid][ybin][type]);
		   h1PtR[pid][ybin][type]->SetNameTitle(Form("h1PtR_%d_%d_%d",pid,ybin,type),"Data 132Sn/108Sn");
		}
		for(auto tbin: ROOT::TSeqI(9)){
			h1EcmR[pid][tbin][type]=(TH1D*)h1Ecm[0][pid][tbin][type]->Clone();
			h1EcmR[pid][tbin][type]->Divide(h1Ecm[3][pid][tbin][type]);
		   h1EcmR[pid][tbin][type]->SetNameTitle(Form("h1EcmR_%d_%d_%d",pid,tbin,type),"Data 132Sn/108Sn");
		}
	}

	TString amdName;
	if(applyPYFiducial) amdName = "AMDSpectrum_pitch-24-24_yaw0-40.root";
	else                amdName = "AMDSpectrum.root";
	auto amdFile = TFile::Open("rootfiles/"+amdName);
	TH1 *h1amdY[nBeam][6][2][2], *h1amdYR[6][2][2];
	TH1 *h1amdPt[nBeam][6][5][2][2], *h1amdPtR[6][5][2][2];
	TH1 *h1amdEcm[nBeam][6][9][2][2], *h1amdEcmR[6][9][2][2];
	std::vector<TString> partName  = {"Neutron","Proton","Deuteron","Triton","3He","4He","6He","6Li","7Li"};
	std::vector<TString> stiffness = {"Soft","Hard"};
	std::vector<TString> cluster   = {"WithCl","NoCl"};
	std::vector<TString> sys       = {"132Sn124Sn","124Sn112Sn","112Sn124Sn","108Sn112Sn"};

	for(auto beam: {0,3})for(auto pid: ROOT::TSeqI(6))for(auto cl: ROOT::TSeqI(2))for(auto eos: ROOT::TSeqI(2)){
		TString name = sys[beam]+"_"+cluster[cl]+"_"+stiffness[eos]+"_"+partName[pid+1];
		amdFile->GetObject("h1Y_"+name, h1amdY[beam][pid][cl][eos]);
		for(auto ybin: ROOT::TSeqI(5)) amdFile->GetObject(Form("h1Pt_"+name+"_%d",ybin), h1amdPt[beam][pid][ybin][cl][eos]);
		for(auto tbin: ROOT::TSeqI(9)) amdFile->GetObject(Form("h1KE_"+name+"_%d",tbin), h1amdEcm[beam][pid][tbin][cl][eos]);
	}
	auto setupAMDTH1 = [&](TH1* h, int cl, int eos){ 
		h->SetMarkerSize(1.2);
		Int_t mark=-1, line=-1;
		Color_t mcol,lcol;
		if(cl==0&&eos==0)      { mark=34; line=1; mcol=kOrange+2;  lcol=kOrange-5; }
		else if(cl==0&&eos==1) { mark=47; line=1; mcol=kTeal+2;    lcol=kTeal-5; }
		else if(cl==1&&eos==0) { mark=28; line=2; mcol=kMagenta+2; lcol=kMagenta-5; }
		else if(cl==1&&eos==1) { mark=46; line=2; mcol=kAzure+2;   lcol=kAzure-5; }
		h->SetMarkerStyle(mark); h->SetMarkerColor(mcol);
		h->SetLineStyle(line);   h->SetLineColor(lcol);
	};
   for(auto pid: ROOT::TSeqI(6))for(auto cl: ROOT::TSeqI(2))for(auto eos: ROOT::TSeqI(2)){
		TString name  = cluster[cl]+"_"+stiffness[eos]+"_"+partName[pid+1];
		TString title = "AMD "+cluster[cl]+","+stiffness[eos];
		//TString title = "132/108,"+cluster[cl]+","+stiffness[eos];
		h1amdYR[pid][cl][eos]=(TH1D*)h1amdY[0][pid][cl][eos]->Clone();
		h1amdYR[pid][cl][eos]->Divide(h1amdY[3][pid][cl][eos]);
		h1amdYR[pid][cl][eos]->SetNameTitle("h1amdYR_"+name,title);
		setupAMDTH1(h1amdYR[pid][cl][eos],cl,eos);
		for(auto ybin: ROOT::TSeqI(5)){
			h1amdPtR[pid][ybin][cl][eos]=(TH1D*)h1amdPt[0][pid][ybin][cl][eos]->Clone();
			h1amdPtR[pid][ybin][cl][eos]->Divide(h1amdPt[3][pid][ybin][cl][eos]);
			h1amdPtR[pid][ybin][cl][eos]->SetNameTitle(Form("h1amdPtR_"+name+"_%d",ybin),title);
			setupAMDTH1(h1amdPtR[pid][ybin][cl][eos],cl,eos);
		}
		for(auto tbin: ROOT::TSeqI(9)){
			h1amdEcmR[pid][tbin][cl][eos]=(TH1D*)h1amdEcm[0][pid][tbin][cl][eos]->Clone();
			h1amdEcmR[pid][tbin][cl][eos]->Divide(h1amdEcm[3][pid][tbin][cl][eos]);
			h1amdEcmR[pid][tbin][cl][eos]->SetNameTitle(Form("h1amdEcmR_"+name+"_%d",tbin),title);
			setupAMDTH1(h1amdEcmR[pid][tbin][cl][eos],cl,eos);
		}
	}
	
	auto jqmdFile = TFile::Open("rootfiles/JQMDSpectrum.root");
	TH1 *h1jqmdY[nBeam][6], *h1jqmdYR[6];
	TH1 *h1jqmdPt[nBeam][6][5], *h1jqmdPtR[6][5];
	TH1 *h1jqmdEcm[nBeam][6][9], *h1jqmdEcmR[6][9];
	std::vector<TString> jqmdPartName  = {"n","1H","2H","3H","3He","4He","6He"};

	for(auto beam: {0,3})for(auto pid: ROOT::TSeqI(6)){
		TString name = sys[beam]+"_"+jqmdPartName[pid+1];
		jqmdFile->GetObject("h1Y_"+name, h1jqmdY[beam][pid]);
		for(auto ybin: ROOT::TSeqI(5)) jqmdFile->GetObject(Form("h1Pt_"+name+"_%d",ybin), h1jqmdPt[beam][pid][ybin]);
		for(auto tbin: ROOT::TSeqI(9)) jqmdFile->GetObject(Form("h1Ecm_"+name+"_%d",tbin), h1jqmdEcm[beam][pid][tbin]);
	}
	auto setupJQMDTH1 = [&](TH1* h){ 
		h->SetMarkerSize(1.2);
		Int_t mark=24, line=1;
		Color_t mcol=kCyan+2,lcol=kCyan-5;
		h->SetMarkerStyle(mark); h->SetMarkerColor(mcol);
		h->SetLineStyle(line);   h->SetLineColor(lcol);
	};
   for(auto pid: ROOT::TSeqI(6)){
		TString name  = partName[pid+1];
		TString title = "JQMD";
		h1jqmdYR[pid]=(TH1D*)h1jqmdY[0][pid]->Clone();
		h1jqmdYR[pid]->Divide(h1jqmdY[3][pid]);
		h1jqmdYR[pid]->SetNameTitle("h1jqmdYR_"+name,title);
		setupJQMDTH1(h1jqmdYR[pid]);
		for(auto ybin: ROOT::TSeqI(5)){
			h1jqmdPtR[pid][ybin]=(TH1D*)h1jqmdPt[0][pid][ybin]->Clone();
			h1jqmdPtR[pid][ybin]->Divide(h1jqmdPt[3][pid][ybin]);
			h1jqmdPtR[pid][ybin]->SetNameTitle(Form("h1jqmdPtR_"+name+"_%d",ybin),title);
			setupJQMDTH1(h1jqmdPtR[pid][ybin]);
		}
		for(auto tbin: ROOT::TSeqI(9)){
			h1jqmdEcmR[pid][tbin]=(TH1D*)h1jqmdEcm[0][pid][tbin]->Clone();
			h1jqmdEcmR[pid][tbin]->Divide(h1jqmdEcm[3][pid][tbin]);
			h1jqmdEcmR[pid][tbin]->SetNameTitle(Form("h1jqmdEcmR_"+name+"_%d",tbin),title);
			setupJQMDTH1(h1jqmdEcmR[pid][tbin]);
		}
	}

	THStack *hsYR[6], *hsPtR[6][5], *hsEcmR[6][9];
	for(auto pid: ROOT::TSeqI(6)){
		hsYR[pid] = new THStack(Form("hsYR_%d",pid),pidName[pid]+" rapidity ratio;y^{0} = y^{cm}/y^{beam};Ratio");
		hsYR[pid]->Add(h1YR[pid][1]);
		if(drawAMD)for(auto cl: ROOT::TSeqI(2))for(auto eos: ROOT::TSeqI(2))hsYR[pid]->Add(h1amdYR[pid][cl][eos]);
		if(drawJQMD)hsYR[pid]->Add(h1jqmdYR[pid]);
		for(auto ybin: ROOT::TSeqI(5)){
			auto ptTitle = Form(pidName[pid]+" p_{T} ratio, %.1f #leq y^{0}< %.1f;p_{T} (GeV/c); Ratio",yBin[ybin],yBin[ybin+1]);
			hsPtR[pid][ybin] = new THStack(Form("hsPtR_%d_%d",pid,ybin),ptTitle);
		   hsPtR[pid][ybin]->Add(h1PtR[pid][ybin][1]);
		   if(drawAMD)for(auto cl: ROOT::TSeqI(2))for(auto eos: ROOT::TSeqI(2))hsPtR[pid][ybin]->Add(h1amdPtR[pid][ybin][cl][eos]);
		   if(drawJQMD)hsPtR[pid][ybin]->Add(h1jqmdPtR[pid][ybin]);
		}
		for(auto tbin: ROOT::TSeqI(9)){
			auto ecmTitle = Form(pidName[pid]+" Ecm ratio, %d #leq#theta< %d;Ecm/A (MeV); Ratio",tbin*20,(tbin+1)*20);
			hsEcmR[pid][tbin] = new THStack(Form("hsEcmR_%d_%d",pid,tbin),ecmTitle);
		   hsEcmR[pid][tbin]->Add(h1EcmR[pid][tbin][1]);
		   if(drawAMD)for(auto cl: ROOT::TSeqI(2))for(auto eos: ROOT::TSeqI(2))hsEcmR[pid][tbin]->Add(h1amdEcmR[pid][tbin][cl][eos]);
		   if(drawJQMD)hsEcmR[pid][tbin]->Add(h1jqmdEcmR[pid][tbin]);
		}
	}

	auto setupTHStackRange = [&](THStack* hs, TH1* h, int pid)
	{
		Double_t ymean=0., counter=0;
		for(auto ibin: ROOT::TSeqI(1,h->GetNbinsX(),1))
			if(h->GetBinContent(ibin)>=0.3&&h->GetBinContent(ibin)<=3.){
				ymean += h->GetBinContent(ibin);
				counter++;
			}
		ymean /= counter;
		Double_t hRange = (pidRange[pid][1]-pidRange[pid][0])/2.;
		hs->SetMinimum(ymean-hRange);
		hs->SetMaximum(ymean+hRange);
	};
	
	auto setupLegend = [&](TLegend* l){ l->SetBorderSize(0); l->SetFillStyle(3001); };
	
	auto cvsYR = new TCanvas("cvsYR","",1600,1000);
	cvsYR->Divide(3,2);
	for(auto pid: ROOT::TSeqI(6)){
		Polish(cvsYR->cd(pid+1));
		cvsYR->cd(pid+1)->SetGrid();
		hsYR[pid]->Draw("ep nostack");
		setupTHStackRange(hsYR[pid],h1YR[pid][1],pid);
		Polish(hsYR[pid]);
		setupLegend(cvsYR->cd(pid+1)->BuildLegend(0.15,0.15,0.53,0.35));
	}
	cvsYR->SaveAs("fig/YRatio.png");
	
	TCanvas *cvsPtR[6];
	for(auto pid: ROOT::TSeqI(6)){
		cvsPtR[pid] = new TCanvas(Form("cvsPtR_%d",pid),"",1800,700);
		cvsPtR[pid]->Divide(4,1);
		for(auto ybin: ROOT::TSeqI(1,5,1)){
			Polish(cvsPtR[pid]->cd(ybin));
			cvsPtR[pid]->cd(ybin)->SetGrid();
			hsPtR[pid][ybin]->Draw("ep nostack");
		   setupTHStackRange(hsPtR[pid][ybin],h1PtR[pid][ybin][1],pid);
			Polish(hsPtR[pid][ybin]);
			setupLegend(cvsPtR[pid]->cd(ybin)->BuildLegend(0.47,0.15,0.85,0.35));
		}
		cvsPtR[pid]->SaveAs("fig/PtRatio_"+pidName[pid]+".png");
	}
	
	TCanvas *cvsEcmR[6];
	for(auto pid: ROOT::TSeqI(6)){
		cvsEcmR[pid] = new TCanvas(Form("cvsEcmR_%d",pid),"",1600,1000);
		cvsEcmR[pid]->Divide(3,2);
		for(auto tbin: ROOT::TSeqI(6)){
			Polish(cvsEcmR[pid]->cd(tbin+1));
			cvsEcmR[pid]->cd(tbin+1)->SetGrid();
			hsEcmR[pid][tbin]->Draw("ep nostack");
		   setupTHStackRange(hsEcmR[pid][tbin],h1EcmR[pid][tbin][1],pid);
			Polish(hsEcmR[pid][tbin]);
			setupLegend(cvsEcmR[pid]->cd(tbin+1)->BuildLegend(0.15,0.15,0.53,0.35));
		}
		cvsEcmR[pid]->SaveAs("fig/EcmRatio_"+pidName[pid]+".png");
	}

	f->Close();
}
