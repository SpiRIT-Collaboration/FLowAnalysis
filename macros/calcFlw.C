#include "FlowFunctions.h"
#include "openFlw.C"
#include "GetLorentzBoost.C"

auto *fcos1 = new TF1("fcos1","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());

Int_t  seltrackID = 4;
UInt_t selReactionPlanef = 10;
Int_t  seltrack;

TCanvas *cc[12];

Double_t rapid_max = 0.4;
Double_t rapid_min = 0.2;
//const Double_t ycm      = 0.388568; // 132Sn + 124Sn before Nov.15 201 //278MeV/u
//                         0: 132Sn + 124Sn,  1: 108Sn + 112Sn  //270MeV/u
//                           132+124,   108+112,  124+112,  112+124,  p + p
const Double_t ycm[]      = {0.382453,  0.364873, 0.390302, 0.354066, 0.371326};
const Double_t ybeam_cm[] = {0.360199,  0.377779, 0.354065, 0.390301, 0.371326};
TLorentzVector* bmVec;
TLorentzVector* tgVec;



const UInt_t nspec = 5;
const UInt_t  nbin = 9;

//UInt_y   y_nbin  = 100;
Double_t y_min[] = {0., 0., 0., 0., 0.};
Double_t y_max[] = {1.2, 1.2, 1.2, 1.2, 1.2};
Double_t y_bin[] = {0.05, 0.05, 0.01, 0.01, 0.01};
Double_t y_binx  =  0.1;

Double_t pt_prmax  =  800.;
Double_t pt_dtmax  = 1000.;
Double_t pt_trmax  = 1400.;
Double_t pt_pimax  = 300.;

Double_t pt_prmin  =  0.;
Double_t pt_dtmin  =  0.;
Double_t pt_trmin  =  0.;
Double_t pt_pimin  =  0.;

UInt_t   pt_nbin   = 100;

Double_t pt_min   = 0;
Double_t pt_max   = 1100;
Double_t pt_dbin  = (pt_max - pt_min)/(Double_t)(pt_nbin-1); 


Double_t yrange1[] = { -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45};
const UInt_t ybin1 = sizeof(yrange1)/sizeof(Double_t);

Double_t yrange2[] = {-0.25, -0.1,  0.1, 0.25, 0.45};
const UInt_t ybin2 = sizeof(yrange2)/sizeof(Double_t);


TString  partname[] = {"pi-","pi+","proton","deuteron","triton"};
UInt_t   partid[]   = {211, 211, 2212, 1000010020, 1000010030};
UInt_t   icol[]     = { 2, 58, 78, 68, 53, 96, 156, 52, 100, 226, 229, 108};
Double_t cutbmass[] = {191.1, 191.1, 1165.9, 2249.9, 3,475.2};

const UInt_t nsys = 4;
const UInt_t nprt = 5;
Color_t icol2[nprt][nsys]= { {kRed,          kBlue,  kOrange-3,   kGreen+1}, 
			    {kBlue+2,   kOrange+7,  kGreen-3,     kPink+9},
			    {kGreen-3,    kPink+7,  kCyan-1,    kYellow-2}, 
			    {kRed-2,      kBlue-2,  kOrange-2,   kGreen-1},
			    {kBlue-3,   kOrange+5,  kGreen+3,     kPink-9} };


TString  iopt[]     = {"","same","same","same","same"};
UInt_t   imark[]    = {20, 21, 22, 23};
  
UInt_t ic = -1;
const Int_t nbinx = 30;

// Retrivew tree
Int_t ntrack[7];
TClonesArray *aArray; // = new TClonesArray("STParticle",100);
TClonesArray *aNLClusterArray;// = new TClonesArray("STNueLANDCluster",100);
Double_t aoq;
Double_t z;
Int_t    snbm;
Double_t ProjA;
Double_t ProjB;
Int_t    mtrack;
Int_t    mtrack_1;
Int_t    mtrack_2;
TVector3 *unitP_fc    = NULL;
TVector3 *unitP_rc    = NULL;
TVector2 *unitP_1     = NULL;
TVector2 *unitP_2     = NULL;
TVector2 *unitP_lang  = NULL;
TBranch  *bunitP_fc;
TBranch  *bunitP_rc;
TBranch  *bunitP_1;
TBranch  *bunitP_2;
TBranch  *brpphi   = 0;
TBranch  *biphi    = 0;
TBranch  *bdeltphi = 0;
TBranch  *bunitP_lang = 0;

// histogram
TH1D* hptpm[4][nbin];
TH1D* hptpp[4][nbin];
TH1D* hptpr[4][nbin];
TH1D* hptdt[4][nbin];
TH1D* hpttr[4][nbin];

// correction
const UInt_t nphi = 30;
std::vector< std::vector< Double_t > >   labphi;
std::vector< std::vector< Double_t > >  subdphi;
std::vector< std::vector< Double_t > >  subcos1;
std::vector< std::vector< Double_t > >  subcos2;
std::vector<Double_t>                allsubcos1;
std::vector<Double_t>                allsubcos2;


ROOT::Math::Interpolator *itplPhi;
Double_t *itplx;
Double_t *itplxe;
Double_t *itpl1;
Double_t *itpl1e;
Double_t atplx;
Double_t atplxe;
Double_t atpl1;
Double_t atpl1e;
Double_t atpl2;
Double_t atpl2e;

// all events
// Double_t mcos1[4]  = {0.512829, 1., 1., 1.} ;
// Double_t mcos1e[4] = {0.0022939,  0., 0., 0.};
// Double_t mcos2[4]  = {0.167427, 1., 1., 1.};
// Double_t mcos2e[4] = {0.0015012, 0., 0., 0.};

// ntrack[4]>=16&&ntrack[4]<32" @ v4.7 
// Double_t mcos1[4]  = {0.5345   ,  1., 1., 1.} ;
// Double_t mcos1e[4] = {0.0027   ,  0., 0., 0.};
// Double_t mcos2[4]  = {0.181853 ,  1., 1., 1.};
// Double_t mcos2e[4] = {0.003047 ,  0., 0., 0.};

// v5.0.0 ntrack[4]>=16&&ntrack[4]<32"
// Double_t mcos1[4]  = {0.53393, 0.53568, 0.53278, 0.55748} ;
// Double_t mcos1e[4] = {0.00289, 0.00366, 0.00149, 0.00339};
// Double_t mcos2[4]  = {0.18149, 0.18268, 0.18071, 0.19785};
// Double_t mcos2e[4] = {0.00197, 0.00251, 0.00788, 0.00241};

// v5.0.0 ntrack[4]>=16"
// Double_t mcos1[4]  = {0.53374, 0.53512, 0.53195, 0.55687} ;
// Double_t mcos1e[4] = {0.00287, 0.00365, 0.00114, 0.00340};
// Double_t mcos2[4]  = {0.18136, 0.18230, 0.18015, 0.19742};
// Double_t mcos2e[4] = {0.00196, 0.00249, 0.00783, 0.00241};

// v7.0.0 // updated on 21 Aug 2018
// Double_t mcos1[]  = {0.506657, 0.516276, 0.496168, 0.524613};
// Double_t mcos1e[] = {0.002763, 0.003706, 0.011740, 0.003536};
// Double_t mcos2[]  = {0.163421, 0.169685, 0.156725, 0.17521};
// Double_t mcos2e[] = {0.001787, 0.002445, 0.007504, 0.002370};

// v10.0.0 //
Double_t mcos1[]  = {0.533779, 0.555799, 0.496168, 0.524613};
Double_t mcos1e[] = {0.001789, 0.002174, 0.011740, 0.003536};
Double_t mcos2[]  = {0.181385, 0.19666 , 0.156725, 0.17521};
Double_t mcos2e[] = {0.001218, 0.001541, 0.007504, 0.002370};




Int_t   iVer;
TString sVer;

// functions
void GetFittingParameters(TH1D &h1, Double_t pp[6]);
void GetFittingParameters(TH1D &h1, Double_t pp[6], Double_t corr[2]);
Double_t GetRapidity_cm(TVector3 p, Double_t mass, TVector3 bvec);
void combineHist(UInt_t sys, TH2D *hyptphi1, TH2D *hyptphi2);

UInt_t pxbooking();
void Plotv1v2(UInt_t selid=2);
void PtDependece(UInt_t hrm=1);
void YDependece(UInt_t hrm=1);
void PhiYbin();
void PhiYSbinf();
void PxDistribution(UInt_t nplot=0);
void hpt_plot();
void meanPx();
void dndy(UInt_t iout=0);
void StoreSubEeventRP(Double_t Phi, Double_t phi_sub);
void ShiftingCorrection(STParticle *apar);
void SaveRPResolution(UInt_t m);
UInt_t   LoadRPResolution(UInt_t m=0);
void     GetRPResolution( UInt_t m=0);  
Double_t GetRPInterpolator(UInt_t m=0, Double_t x=0);
UInt_t   SetBranch(UInt_t m);
void GetRPResolutionwChi(UInt_t m=0);
void PlotNeuLANDv1v2();
void PlotNeuLANDProperty(UInt_t iout=0);


//--------- main ----------//
void calcFlw() 
{
  gROOT->Reset();
  openFlw();

  for(UInt_t i = 0; i < 4; i++){
    if(rChain[i] != NULL) {    
      std::cout << " System " << i << " "  << isys[i] << "  -> " << sysName[isys[i]] << std::endl; 
    }
  }

  if(rChain[0] == NULL)
    exit(0);

  std::cout  << " Total  = " << m_end  << std::endl;
  
  gROOT->ProcessLine(".! grep -i void calcFlw.C | grep '//%%'");

  //  PlotNeuLANDProperty(1);
}

void GetRPResolutionwChi(UInt_t m)            //%% Executable : 
{
  auto hphi0_180  = new TH1D("hphi0_180" ,"#Phi from  0 to 180; #Phi",100,0.,3.2);
  auto hphi90_180 = new TH1D("hphi90_180","#Phi from 90 to 180; #Phi",100,0.,3.2);

  //  TCut mcut = "mtrack_1>0&&mtrack_2>0&&ntrack[4]>=16&&ntrack[4]<32"; // mid-central
  //  TCut mcut = "mtrack_1>0&&mtrack_2>0&&ntrack[4]>=16"; // mid-central
  //  TCut mcut = "mtrack_1>0&&mtrack_2>0&&ntrack[4]>=32"; // most-central
  TCut mcut = "mtrack_1>0&&mtrack_2>0"; //all

  rChain[m]->SetAlias("dPhi","abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi()))");

  rChain[m]->Project("hphi0_180" ,"dPhi",mcut);
  rChain[m]->Project("hphi90_180","dPhi",mcut&&"dPhi>TMath::Pi()/2.");

  Double_t m0 = hphi0_180->GetEntries();
  Double_t m1 = hphi90_180->GetEntries();

  Double_t chi = sqrt(-2.* log(2.* m1/m0));
  
  Double_t m01e = m1/m0*sqrt(1./m1+1./m0);
  Double_t chie = chi - sqrt(-2.* log(2. * (m1/m0+m01e)));

  mcos1[m]  = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*chi;
  mcos1e[m] = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*(chi+chie) - mcos1[m];
  mcos2[m]  = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi,2);
  mcos2e[m] = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi+chie,2) - mcos2[m];

  std::cout << m1 << " / " << m0 << " = " << m1/m0 << " -> Chi " << chi << " +- " << chie
	    << std::endl;
  std::cout << " <cos(Phi)> = " << mcos1[m] << " +- " << mcos1e[m] 
	    << " <cos(2Phi)> = "<< mcos2[m] << " +- " << mcos2e[m] 
	    << std::endl;

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hphi0_180->Draw();
  hphi90_180->SetLineColor(2);
  hphi90_180->Draw("same");

  Double_t para[6];
  GetFittingParameters( *hphi0_180, para);

  cout << std::setw(14) << para[1] << " +- " << std::setw(10) << para[3]
       << std::endl;

}


void GetRPResolution(UInt_t m)                //%% Executable : Plot Phi and subevent, Phi_A and Phi_B correlation
{

  // Booking
  Double_t phi_min = -3.2;
  Double_t phi_max =  3.2;

  TString hname = Form("hrpphi%d",m);
  auto hrpphi = new TH1D(hname,sysName[isys[m]]+";#Phi ",nphi, phi_min, phi_max);
  hrpphi -> SetLineColor(icol[isys[m]]);
  
  hname = Form("hdltphi%d",m);
  auto hdltphi = new TH2D(hname,sysName[isys[m]]+";#Phi ; #Phi_A - #Phi_B",nphi,phi_min, phi_max, nphi,phi_min, phi_max);
  hdltphi -> SetLineColor(icol[isys[m]]);

  hname = Form("hsubphi%d",m);
  auto hsubphi = new TH2D(hname,sysName[isys[m]]+";#Phi_A ; #Phi_B",nphi, phi_min, phi_max, nphi,phi_min, phi_max);

    
  // Retreview

  labphi.resize(nphi);
  subdphi.resize(nphi);
  subcos1.resize(nphi);

  Int_t nEntry = SetBranch(m);
  
  for(Int_t i = 0; i < nEntry; i++){
    rChain[m]->GetEntry(i);

    if(i%10000 == 0) 
      std::cout << "GetRPResolution Processing ... " << i << " : " << unitP_fc->Phi() << std::endl; 

    if(mtrack_1 > 0 && mtrack_2 > 0){
      hrpphi  -> Fill(TVector2::Phi_mpi_pi(unitP_fc->Phi()));
      hsubphi -> Fill(TVector2::Phi_mpi_pi(unitP_1->Phi()),  TVector2::Phi_mpi_pi(unitP_2->Phi()) );
      hdltphi -> Fill(TVector2::Phi_mpi_pi(unitP_fc->Phi()), TVector2::Phi_mpi_pi(unitP_1->Phi() - unitP_2->Phi()));

      StoreSubEeventRP(unitP_fc->Phi(), TVector2::Phi_mpi_pi(unitP_1->Phi() - unitP_2->Phi()));
    }
  }    


  SaveRPResolution(m);
  auto gscos = new TGraphErrors(nphi, itplx, itpl1, itplxe, itpl1e);
  gscos->SetName("gscos");
  gscos->SetTitle(sysName[isys[m]]+";#Phi; <cos(#phi_A-#phi_B)>");

  // Draw
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1500, 400*m_end);
  cc[ic]->Divide(4,m_end);
  
  UInt_t id = 1;
  
  cc[ic]->cd(id); id++;
  hrpphi -> SetLineColor(4);
  hrpphi -> Draw();

  cc[ic]->cd(id); id++;
  hsubphi-> Draw("colz");
  
  cc[ic]->cd(id); id++;
  hdltphi-> Draw("colz");

  cc[ic]->cd(id); id++;
  gscos->Draw();
}


void StoreSubEeventRP(Double_t Phi, Double_t phi_sub)
{
  Double_t minphi =  -1.*TMath::Pi();
  Double_t dltphi =  -2.*minphi/nphi;

  Phi     = TVector2::Phi_mpi_pi( Phi );

  for(UInt_t iphi = 0; iphi < (UInt_t)nphi; iphi++){
    
    if( Phi < minphi + dltphi*(iphi+1)) {
      labphi[iphi] .push_back(Phi);
      subdphi[iphi].push_back(phi_sub);
      break;
    }
  }
}


void SaveRPResolution(UInt_t m)
{
  Bool_t bPlot = kFALSE;

  Double_t phi_min = -3.2;
  Double_t phi_max =  3.2;


  itplx = new Double_t[nphi];
  itplxe= new Double_t[nphi];
  itpl1 = new Double_t[nphi];
  itpl1e= new Double_t[nphi];

  Double_t *getx = new Double_t[2];
  
  TH1D *hphi1[nphi];

  auto hsubphia1 = new TH1D("hsubphia1","all dphi1 "  , 100, phi_min,    phi_max);
  auto hsubphi1  = new TH1D("hsubphi1" ,"sub dphi1 "  , 100, phi_min,    phi_max);

  if(bPlot){
    ic++; 
    cc[ic]   = new TCanvas("subphi1", "subphi1"  ,2000, 2000);
    cc[ic]->Divide(3,nphi/3);
  }

  for(UInt_t jn = 0; jn < nphi; jn++){
    cout << " jn " << jn << " :size " << subdphi[jn].size() << endl;

    getx  = vMean( labphi[jn] );
    itplx[jn]  = getx[0];
    itplxe[jn] = getx[1];

    hsubphi1->Reset();
    for(UInt_t ii = 0; ii < (UInt_t)subdphi[jn].size(); ii++){
      hsubphia1->Fill(    subdphi[jn].at(ii) );
      hsubphi1 ->Fill(    subdphi[jn].at(ii) );
    }

    // Plot fitting  cos(delta Phi) depending on the Phi_RP
    if(bPlot){
      TString hname = Form((TString)hsubphi1->GetName()+"_%d",jn);
      hphi1[jn] = (TH1D*)hsubphi1->Clone();
      hphi1[jn]->SetName(hname);
      cout << " entries " << hphi1[jn]->GetEntries() << endl;
      cc[ic]->cd(jn+1);
      hphi1[jn]->Fit("fcos1","","",-3.,3.);
    }


    // fit cos(delta Phi)
    Double_t para[6];
    GetFittingParameters( *hsubphi1, para);

    itpl1[jn]  = para[1];
    itpl1e[jn] = para[3];;

    std::cout << setw(5) << jn 
	      <<  " lab phi " << setw(10) << itplx[jn] 
	      <<  " sub 1 " << setw(14) << para[1]  << " +- " << setw(10) << para[3]
	      << std::endl;
    std::cout << " ====================-====================" << std::endl;
  }


  ic++;
  cc[ic]   = new TCanvas(Form("cc%d",ic)  ,Form("cc%d",ic));

  Double_t para[6];
  GetFittingParameters(*hsubphia1, para);
  atpl1  = para[1];
  atpl1e = para[3];

  hsubphia1->Draw();


  TString fName = "RP" +  sysName[isys[m]] + sDB[m] + ".data";
  gSystem->cd("db");

  std::fstream fout;
  fout.open(fName, std::fstream::out);

  fout << "     ,    PHI     +- PHI_error ,    v1      +-   v1_error "  << std::endl; 

  fout << ">      ALL                     " 
       << std::setw(14) << atpl1 << " +- " << std::setw(10) << atpl1e
       << std::endl;
  fout << "##---- " << std::endl;

  for(UInt_t i = 0; i < nphi; i++) 
    fout << "> "
	 << std::setw(5) << i 
	 << std::setw(12) << itplx[i] << " +- " << std::setw(10) << itplxe[i]
	 << std::setw(14) << itpl1[i] << " +- " << std::setw(10) << itpl1e[i] 
	 << std::endl;
  
  gSystem->cd("..");

  std::cout << "db/" << fName << " is saved. " << std::endl;
  fout.close();


}

UInt_t LoadRPResolution(UInt_t m)
{
  std::vector<UInt_t>   aksize;
  std::vector<Float_t>  arpphi;
  std::vector<Float_t>  arpphie;
  std::vector<Float_t>  arpres1;
  std::vector<Float_t>  arpres1e;

  TString fName = "RP" +  sysName[isys[m]] + sDB[m] + ".data";
  gSystem->cd("db");

  std::fstream fin;
  fin.open(fName, std::fstream::in); 
  if(fin)
    std::cout << fName << " is opened. " << std::endl;
  else {
    std::cout << fName <<" is not opened. " << std::endl;
    return 0;
  }

  TString sget;
  while(!fin.eof()){
    fin >> sget;
    
    if(sget != ">") continue;

    fin >> sget;
    if(sget == "ALL") {
      //v1
      fin >> sget;
      Float_t fget = atof(sget);
      atpl1 = fget;
      
      fin >> sget;
      if(sget != "+-") {
	std::cout << "Error loading RP resolution 1" << sget << std::endl;
        break;
      }

      fin >> sget;
      fget = atof(sget);
      atpl1e = atof(sget);
      
    }
    else {

      UInt_t iget = atoi(sget);
      aksize.push_back(iget);
    
      if(fin.eof()) break;

      fin >> sget;
      Float_t fget = atof(sget);
      arpphi.push_back(fget);
    

      if(fin.eof()) break;

      fin >> sget;
      if(sget != "+-") {
	std::cout << "Error loading RP resolution 1 " << sget << std::endl;
	break;
      }

      // read <cos(dphi)>
      if(fin.eof()) break;
      fin >> sget;

      fget = atof(sget);
      arpphie.push_back(fget);

      if(fin.eof()) break;
      fin >> sget;
      fget = atof(sget);
      arpres1.push_back(fget);
    
      if(fin.eof()) break;
      fin >> sget;

      if(sget != "+-") {
	std::cout << "Error loading RP resolution 2 " <<  sget << std::endl;
	break;
      }
  
      if(fin.eof()) break;
      fin >> sget;
    
      fget = atof(sget);
      arpres1e.push_back(fget);


    }
  }

  fin.close();
  gSystem->cd("..");

  std::cout << " ALL "
	    << std::setw(14) << atpl1 << " +- " << std::setw(10) << atpl1e
	    << std::endl;

  const UInt_t vsize = arpres1.size();
  itplx  = new Double_t[vsize];
  itplxe = new Double_t[vsize];
  itpl1  = new Double_t[vsize];
  itpl1e = new Double_t[vsize];
  

  for(UInt_t i = 0; i < vsize; i++){
    itplx[i]  = arpphi.at(i);
    itplxe[i] = arpphi.at(i);
    itpl1[i]  = arpres1.at(i);
    itpl1e[i] = arpres1e.at(i);

    std::cout << std::setw(4) << i << " : "
  	      << std::setw(12) << itplx[i] << " +- " << std::setw(10) << itplxe[i]
  	      << std::setw(14) << itpl1[i] << " +- " << std::setw(10) << itpl1e[i] 
  	      << std::endl;
  }

  return vsize;
}

Double_t GetRPInterpolator(UInt_t m, Double_t x)
{
  static Bool_t bfirst = kTRUE;

  UInt_t vsize;
  if(bfirst){
    vsize = LoadRPResolution(m);
    if(vsize == 0) {
      GetRPResolution(m);
      vsize = nphi;
    }
    bfirst = kFALSE;

    itplPhi = new ROOT::Math::Interpolator(vsize,ROOT::Math::Interpolation::kAKIMA);
    itplPhi->SetData(vsize, itplx, itpl1);
  }

  Double_t value = atpl1;

  if(x < 3.036 && x > -3.037) {
    
    value = itplPhi->Eval(x);
    if( std::isnan(value)) 
      value = atpl1;
  }

  //  cout << x << " value " << value << endl;

  return value;
}

void PlotAcceptance(UInt_t m = 0, UInt_t selid = 2)         //
{
  gStyle->SetOptStat(0);

  TVector3 boostVec = LorentzBoost(4);

  auto haccp    = new TH2D("haccp"   , partname[selid]+"; Rapidity ; Pt [MeV/c]",100, -0.4, 0.6, 100,  0.,600.);

  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  Int_t RPflag = 0;
  if(selid >= 2) RPflag = 10;

  Int_t nevt = SetBranch(m);
  for(Int_t i = 0; i < nevt; i++){
    rChain[m]->GetEntry(i);

    TIter next(aArray);
    STParticle *aPart = NULL;
  
    while( (aPart = (STParticle*)next()) ) {

      auto rpf = aPart->GetReactionPlaneFlag();
      auto pid = aPart->GetPID();

      if(pid == partid[selid] && aPart->GetCharge() == pcharge ){

	if( pid > 2000 && (rpf == 110 || rpf == 210 ) ) continue;

	//	auto rapid = aPart-> GetRapidity();
	//	rapid = (rapid - ycm[isys[m]])/ybeam_cm[isys[m]];

	auto rapid = GetRapidity_cm(aPart->GetRotatedMomentum(), aPart->GetMass(), -boostVec);
	auto pt    = aPart-> GetRotatedMomentum().Pt();

	haccp -> Fill(rapid, pt);
      }
    }
  }

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  haccp->Draw("colz");
}


void Plotv1v2(UInt_t selid=2)                 //%% Executable : v1  as a function of rapidity
{
  if(selid > 4 ) return;
      
  gStyle->SetGridColor(7);
  gStyle->SetGridStyle(1);

  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  Int_t RPflag = 0;
  if(selid >= 2) RPflag = 10;

  cout << " Particle " << partname[selid] << endl;

  UInt_t nybin  = 15;

  UInt_t nbin1 = 100;
  UInt_t nbin2 = 100;

  auto hdphi1  = new TH2D("hdphi1" ,";#Delta(#phi-#Phi)  ; Rapidity"      , nbin1, -3.2, 3.2, nybin, -0.4, 0.4);
  auto hdphi2  = new TH2D("hdphi2" ,";2*#Delta(#phi-#Phi); Rapidity"      , nbin2, -3.2, 3.2, nybin, -0.4, 0.4);
  auto hcos2   = new TH2D("hcos2"  ,";cos(2*#Delta(#phi-#Phi)); Rapidity" , nbin1, -1. , 1. , nybin, -0.4, 0.4);

  for(Int_t m = m_bgn; m < m_end; m++){
    
    hdphi1->Reset();
    hdphi2->Reset();
    hcos2->Reset();

    TVector3 boostVec = LorentzBoost(4);

    Int_t nevt = SetBranch(m);
    cout << " Number of events " << nevt << endl;


    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      if(ntrack[2] == 0) continue;

      //      Double_t fcphi = unitP_fc->Phi();
      //      Double_t RPres = GetRPInterpolator(m, fcphi);

      //      if(ntrack[4] < 16 ||  ntrack[4] <= 32) continue;
      if(ntrack[4] < 16 ) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto rpf = aPart->GetReactionPlaneFlag();
	auto pid = aPart->GetPID();
	auto bmass = aPart->GetBBMass();

	if(pid == partid[selid] && aPart->GetCharge() == pcharge && bmass <= cutbmass[selid] ){

	  if( pid > 2000 && rpf < 2000 ) continue;

	  auto rapid = GetRapidity_cm(aPart->GetMomentum(), aPart->GetMass(), -boostVec);
	  //	  auto rapid = aPart-> GetRapidity();

	  //	  cout << "rapid " << rapid << " : " << aPart-> GetRapidity() << endl;

	  if(aPart->GetIndividualRPAngle() > -9 ) {

	    hdphi1->Fill(                        aPart->GetAzmAngle_wrt_RP() , rapid);
	    hdphi2->Fill(TVector2::Phi_mpi_pi(2.*aPart->GetAzmAngle_wrt_RP()), rapid);
	    hcos2 ->Fill(                 cos(2.*aPart->GetAzmAngle_wrt_RP()), rapid);
	  }
	}
      } //while( (aPart = (STParticle*)next()) ) {
    }

    TH1D *hydphi1[nybin];
    TH1D *hydphi2[nybin];
    TH1D *hycos2[nybin];

    TString fName = "VN" + sysName[isys[m]] + "_" + partname[selid]+".root";
    gSystem->cd("data");
    auto GraphSave = new TFile(fName,"recreate");

    auto gv_v1 = new TGraphErrors();
    gv_v1->SetName("gv_v1");
    gv_v1->SetTitle(partname[selid]+"; Rapidity ; v1 ");

    auto gv_v1c = new TGraphErrors();
    gv_v1c->SetName("gv_v1c");
    gv_v1c->SetTitle(partname[selid]+"; Rapidity ; v1 (w/o corr.) ");

    auto gv_v2 = new TGraphErrors();
    gv_v2->SetName("gv_v2");
    gv_v2->SetTitle(partname[selid]+"; Rapidity ; v2 ");

    auto gv_v2c = new TGraphErrors();
    gv_v2c->SetName("gv_v2c");
    gv_v2c->SetTitle(partname[selid]+"; Rapidity ; v2 (w/o corr.)");

    auto gv_v2d = new TGraphErrors();
    gv_v2d->SetName("gv_v2d");
    gv_v2d->SetTitle(partname[selid]+"; Rapidity ; <cos 2#phi >");
     

    std::cout << " ---- Resutls ---------------------" << std::endl;
    std::cout << " < cos phi > " << mcos1[isys[m]] << std::endl;

    ic++; id = 1; UInt_t idd = 1;
    cc[ic]   = new TCanvas("dphi1","dphi1",1500,1000);
    cc[ic]->Divide(3, nybin/3);
    cc[ic+1] = new TCanvas("dphi2","dphi2",1500,1000);
    cc[ic+1]->Divide(3, nybin/3);

    UInt_t ip1 = 0;
    UInt_t ip2 = 0;
    for(UInt_t jn = 0; jn < nybin; jn++){

      Double_t rpd  = hdphi1->GetYaxis()->GetBinCenter(jn+1);
      Double_t rpde = hdphi1->GetYaxis()->GetBinWidth(jn)/sqrt(12.);
      
      cc[ic]->cd(id); id++;
      hydphi1[jn] = (TH1D*)hdphi1->ProjectionX((TString)Form("hydphi1_%d",jn+1),jn+1, jn+1,"eo");
      
      if(hydphi1[jn]->GetEntries() > 50 && rpd <= 1.){
	
	Double_t para[6];
	Double_t corr[2] = {mcos1[isys[m]], mcos1e[isys[m]]};
	GetFittingParameters(*hydphi1[jn], para, corr);

	Double_t v1c  = para[1];
	Double_t v1ce = para[3];
	
	gv_v1c->SetPoint(ip1, rpd, v1c);
	gv_v1c->SetPointError(ip1, rpde, v1ce);

	Double_t v1   = para[4];
	Double_t v1e  = para[5];
	
	gv_v1->SetPoint(ip1, rpd, v1);
	gv_v1->SetPointError(ip1, rpde, v1e);
	ip1++;
	
	std::cout << setw(5) << jn << " w  c : " << setw(12)
		  << rpd  << " +- " << rpde 
		  <<  " v1 " << setw(12) << v1c << " +- " << setw(10) << v1ce << std::endl;
	std::cout << " ----------------------------------" << std::endl;
      }


      cc[ic+1]->cd(idd); idd++;
      hydphi2[jn] = (TH1D*)hdphi2->ProjectionX((TString)Form("hydphi2_%d",jn+1),jn+1, jn+1,"eo");
    
      if(hydphi2[jn]->GetEntries() > 50 && rpd <= 1.){

	Double_t para[6];
	Double_t corr[2]={mcos2[isys[m]], mcos2e[isys[m]]};
	GetFittingParameters( *hydphi2[jn], para, corr);

	Double_t v2c  = para[1];
	Double_t v2ce = para[3];
	gv_v2c->SetPoint(ip2, rpd, v2c);
	gv_v2c->SetPointError(ip2, rpde, v2ce);

	Double_t v2   = para[4];
	Double_t v2e  = para[5];
	
	gv_v2->SetPoint(ip2, rpd, v2);
	gv_v2->SetPointError(ip2, rpde, v2e);
	ip2++;

     	std::cout << setw(5) << jn << " w  c : " << setw(12)
		  << rpd  << " +- " << rpde 
		  <<  " v2 " << setw(12) << v2c << " +- " << setw(10) << v2ce << std::endl;
	std::cout << " ==================================" << std::endl;

	
      }
    }


    gv_v1->Write();
    gv_v2->Write();


    // plotting
    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

    auto mv1 = new TMultiGraph("mv1","; Rapidity; v1");
    mv1->Add(gv_v1,"lp");
    mv1->Add(gv_v1c,"lp");

    gv_v1->SetLineColor(4);
    gv_v1->SetMarkerStyle(20);
    gv_v1->SetMarkerColor(4);

    gv_v1->Print();

    mv1->Draw("alp");

    auto Ymin = mv1->GetYaxis()->GetXmin();
    auto Ymax = mv1->GetYaxis()->GetXmax();
    auto Xmin = mv1->GetXaxis()->GetXmin();
    auto Xmax = mv1->GetXaxis()->GetXmax();

    auto aLineX1 = new TLine(Xmin, 0., Xmax, 0.);
    aLineX1->SetLineColor(1);
    aLineX1->SetLineStyle(3);
    aLineX1->Draw();


    auto aLineY1 = new TLine(0., Ymin, 0., Ymax);
    //    auto aLineY1 = new TLine(ybeam_cm[isys[m]], Ymin, ybeam_cm[isys[m]], Ymax);
    aLineY1->SetLineColor(1);
    aLineY1->SetLineStyle(3);
    aLineY1->Draw();


    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    
    auto mv2 = new TMultiGraph("mv2","; Rapidity; v2");
    mv2->Add(gv_v2,"lp");
    mv2->Add(gv_v2c,"lp");

    gv_v2->SetLineColor(2);
    gv_v2->SetMarkerStyle(20);
    gv_v2->SetMarkerColor(2);

    //     gv_v2->Draw("alp");
    mv2->Draw("apl");

    Ymin = mv2->GetYaxis()->GetXmin();
    Ymax = mv2->GetYaxis()->GetXmax();
    Xmin = mv2->GetXaxis()->GetXmin();
    Xmax = mv2->GetXaxis()->GetXmax();

    auto aLineX2 = new TLine(0., Ymin, 0., Ymax);
    //    auto aLineX2 = new TLine(ybeam_cm[isys[m]], Ymin, ybeam_cm[isys[m]], Ymax);
    aLineX2->SetLineColor(1);
    aLineX2->SetLineStyle(3);
    aLineX2->Draw();
  
  }

  gSystem->cd("..");
}

void dndy(UInt_t iout)                        //%% Executable : Make plots of dNdy for p, d, t, pi+-
{

  //----- booking
  TFile* hout;
  TH1D* hrap[4][5];
  TH1D* hnpart[4][5];
  TH1D* hmtrack[4];

  auto aLeg0 = new TLegend(0.7,0.7,0.9 ,0.9,"");
  auto aLeg1 = new TLegend(0.7,0.7,0.9 ,0.9,"");
  

  for(Int_t m = m_bgn; m < m_end; m++){


    for(UInt_t i = 0; i < 5; i++){

      TString hname = Form("hrap%d_%d",m,i);
      TString htitle= partname[i] + "; Rapidity; dN/dy";
      UInt_t  y_nbin = UInt_t((y_max[i] - y_min[i])/y_bin[i]);

      hrap[m][i] = new TH1D(hname, htitle, y_nbin, y_min[i], y_max[i]);
      hrap[m][i] ->SetLineColor(icol2[i][isys[m]]);

      hname  = Form("hnpart%d_%d",m,i);
      htitle = partname[i]+" ; Multiplicity";
      hnpart[m][i] = new TH1D(hname, htitle, 25,0,25);
      hnpart[m][i]->SetLineColor(icol2[i][m]);
    }

    hmtrack[m] = new TH1D(Form("hmtrack%d",m),"Number of good tracks; Multiplicity",80, 0, 80);
    hmtrack[m] -> SetLineColor(icol[isys[m]]);
    
    aLeg0->AddEntry(hrap[m][4],sysName[isys[m]],"lp");
    aLeg1->AddEntry(hnpart[m][4],sysName[isys[m]],"lp");
  }

  //------------------------

  for(Int_t m = m_bgn; m < m_end; m++){

    Int_t nEntry =  SetBranch(m);

    for(Int_t i = 0; i < nEntry; i++){
      rChain[m]->GetEntry(i);

      TIter next(aArray);
      STParticle *aPart = NULL;

      hmtrack[m]->Fill(ntrack[3]);
      UInt_t npart[5] = {0,0,0,0,0};

      while( (aPart = (STParticle*)next()) ) {

	auto flag  = aPart->GetReactionPlaneFlag();
	auto pid   = aPart->GetPID();
	auto charg = aPart->GetCharge();
	auto p     = aPart->GetRotatedMomentum();
	auto rapid = aPart->GetRapidity();
	auto theta = aPart->GetRotatedMomentum().Theta();

	if(pid == partid[0] ){
	  if( charg < 0){
	    hrap[m][0]->Fill(rapid);
	    npart[0]++;
	  }
	  else {
	    hrap[m][1]->Fill(rapid);
	    npart[1]++;
	  }
	}
	
	if(flag > 100){
	  for(UInt_t i = 2; i < 5; i++){
	    if(pid == partid[i]){
	      hrap[m][i]->Fill(rapid);
	      npart[i]++;
	    }
	  }
	}
      }
      for(UInt_t i = 0; i < 5; i++)
	hnpart[m][i]->Fill(npart[i]);
    }

    //----- Output file   
    if( iout == 1 ){
      TString fName = "dNdySn" + sysName[isys[m]] + ".root";
      gSystem->cd("data");
      hout = new TFile(fName,"recreate");
      for(UInt_t i = 0; i < 5; i++) {
	hnpart[m][i]->Write();
	hrap[m][i]  ->Write();
      }
      hmtrack[m]->Write();
      hout->Close();
      
      gSystem->cd("..");
    }
   
    else { 

      std::cout << " -------------------- " << std::endl;
      for(UInt_t i = 0; i < 5; i++) {

	//----- Normalizing
	auto scl = 1./nEntry * (Double_t)hrap[m][i]->GetNbinsX() / (y_max[i] -y_min[i]);
	//hrap[m][i]->Scale(scl);

	scl = 1./nEntry;
	//hnpart[m][i]->Scale(scl);
	auto mean   = hnpart[m][i]->GetMean();
	auto meaner = hnpart[m][i]->GetMeanError();
	std::cout << setw(12) << partname[i] << " : " << mean << " +- " << meaner << std::endl;

      }
      hmtrack[m]->Scale(1./nEntry);
      std::cout << " -------------------- " << std::endl;
    }
  }

  if( iout == 1 ) return;

  //----- Drawing
  //----- canvas
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1400,500);
  cc[ic]->Divide(3,1);

  const UInt_t nsys = m_end - 1;


  Double_t hmax[nsys];


  // plot p, d, and t
  for(UInt_t ip = 2; ip < 5; ip++){
    cc[ic]->cd(ip-1); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hrap[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    cout << "ip : " << ip  << " gmax " << gmax << endl;
    
    hrap[0][ip]->SetMaximum(gmax*1.1);

    for(UInt_t m = m_bgn; m < m_end; m++) 
      hrap[m][ip]->Draw(iopt[m]); 
     
    if(ip == 4)
      aLeg0->Draw();
  }


  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
  cc[ic]->Divide(2,1);

 
  for(UInt_t ip = 0; ip < 2; ip++){
    cc[ic]->cd(ip+1); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hrap[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    hrap[0][ip]->SetMaximum(gmax*1.1);
    for(UInt_t m = m_bgn; m < m_end; m++) 
      hrap[m][ip]->Draw(iopt[m]); 
  

    if(ip == 1)
      aLeg1->Draw();
  }


  // //----cc1
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,800);
  cc[ic]->Divide(3,2);
  
  for(UInt_t ip = 2; ip < 5; ip++){
    cc[ic]->cd(ip-1); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hnpart[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    hnpart[0][ip]->SetMaximum(gmax*1.1);
    for(UInt_t m = m_bgn; m < m_end; m++) 
      hnpart[m][ip]->Draw(iopt[m]);

    if(ip == 4)
      aLeg1->Draw();
  }

  for(UInt_t ip = 0; ip < 2; ip++){
    cc[ic]->cd(ip+4); 

    for(UInt_t m = m_bgn; m < m_end; m++)
      hmax[m] = hnpart[m][ip]->GetMaximum();

    Double_t gmax = TMath::MaxElement(nsys, hmax);
    hnpart[0][ip]->SetMaximum(gmax*1.1);
    for(UInt_t m = m_bgn; m < m_end; m++)
      hnpart[m][ip]->Draw(iopt[m]);
    

    if(ip == 1)
      aLeg1->Draw();
  }
  
}


void meanPx()                                 //%% Executable : Make  plots of <px> vs rapidity 
{
  PxDistribution(1);

  TGraphErrors  *gpr[4];
  TGraphErrors  *gdt[4];
  TGraphErrors  *gtr[4];
  TGraphErrors  *gpm[4];
  TGraphErrors  *gpp[4];

  auto mgpdt = new TMultiGraph();
  mgpdt->SetTitle("p,d and t; Rapidity_lab; <px>/A [MeV/c]");
  auto aLeg0 = new TLegend(0.1,0.7,0.35,0.9,"");

  auto mgpi  = new TMultiGraph();
  mgpi->SetTitle("#pi^{+-}; Rapidity_lab; <px> [MeV/c]");
  auto aLeg1 = new TLegend(0.1,0.7,0.35,0.9,"");

  for(Int_t m = m_bgn; m < m_end; m++){
   
    Double_t rap[nbin];
    Double_t rape[nbin];
 
    Double_t mptpm[nbin];
    Double_t mptpp[nbin];
    Double_t mptpr[nbin];
    Double_t mptdt[nbin];
    Double_t mpttr[nbin];

    Double_t mptpme[nbin];
    Double_t mptppe[nbin];
    Double_t mptpre[nbin];
    Double_t mptdte[nbin];
    Double_t mpttre[nbin];

    UInt_t   npr = 0;
    UInt_t   ndt = 0;
    UInt_t   ntr = 0;
    UInt_t   npm = 0;
    UInt_t   npp = 0;


    // get mean 
    for(UInt_t k = 0; k < nbin; k++){
      rap[k]   = y_min[0] + y_binx*(k+0.5);
      rape[k]  = y_binx*0.5;

      if(hptpm[m][k]->GetEntries() > 0){
	mptpm[npm]  = hptpm[m][k]->GetMean(); 
	mptpme[npm] = hptpm[m][k]->GetMeanError(); 
	npm++;
      }

      if(hptpp[m][k]->GetEntries() > 0){
	mptpp[npp]  = hptpp[m][k]->GetMean();
	mptppe[npp] = hptpp[m][k]->GetMeanError();
	npp++;
      }

      if(hptpr[m][k]->GetEntries() > 0){
	mptpr[npr]  = hptpr[m][k]->GetMean();
	mptpre[npr] = hptpr[m][k]->GetMeanError();
	npr++;
      }

      if(hptdt[m][k]->GetEntries() > 0){
	mptdt[ndt]  = hptdt[m][k]->GetMean()/2.;
	mptdte[ndt] = hptdt[m][k]->GetMeanError();
	ndt++;
      }

      if(hpttr[m][k]->GetEntries() > 0){
	mpttr[ntr] = hpttr[m][k]->GetMean()/3.;
	mpttre[ntr] = hpttr[m][k]->GetMeanError();
	ntr++;
      }
    }


    gpr[m] = new TGraphErrors(npr, rap, mptpr, rape, mptpre);
    gdt[m] = new TGraphErrors(ndt, rap, mptdt, rape, mptdte);
    gtr[m] = new TGraphErrors(ntr, rap, mpttr, rape, mpttre);
    gpm[m] = new TGraphErrors(npm, rap, mptpm, rape, mptpme);
    gpp[m] = new TGraphErrors(npp, rap, mptpp, rape, mptppe);
   
  
    TString atitle = Form("gpr%d",m);
    gpr[m]->SetName(atitle);
    gpr[m]->SetLineColor(icol[isys[m]]);
    gpr[m]->SetMarkerStyle(imark[isys[m]]);
    gpr[m]->SetMarkerColor(icol[isys[m]]);
    gpr[m]->SetTitle("Proton; y_lab; <Px> (MeV/c)");
    mgpdt->Add(gpr[m],"lp");
    aLeg1->AddEntry(gpr[m],"proton   "+sysName[isys[m]],"lp");

    atitle = Form("gdt%d",m);
    gdt[m]->SetName(atitle);
    gdt[m]->SetLineColor(icol[isys[m]]);
    gdt[m]->SetMarkerStyle(21);
    gdt[m]->SetMarkerColor(icol[isys[m]]);
    gdt[m]->SetTitle("Deuteron; y_lab; <Px> (MeV/c)");
    mgpdt->Add(gdt[m],"lp");
    aLeg1->AddEntry(gdt[m],"deuteron "+sysName[isys[m]],"lp");

    atitle = Form("gtr%d",m);
    gtr[m]->SetName(atitle);
    gtr[m]->SetLineColor(icol[isys[m]]);
    gtr[m]->SetMarkerStyle(22);
    gtr[m]->SetMarkerColor(icol[isys[m]]);
    gtr[m]->SetTitle("Triton; y_lab; <Px> (MeV/c)");
    mgpdt->Add(gtr[m],"lp");
    aLeg1->AddEntry(gtr[m],"trition  "+sysName[isys[m]],"lp");

    atitle = Form("gpm%d",m);
    gpm[m]->SetName(atitle);
    gpm[m]->SetLineColor(icol[isys[m]]);
    gpm[m]->SetMarkerStyle(23);
    gpm[m]->SetMarkerColor(icol[isys[m]]);
    gpm[m]->SetTitle("pi-; y_lab; <Px> (MeV/c)");
    mgpi->Add(gpm[m],"ip");
    aLeg0->AddEntry(gpm[m],"#pi^{-}   "+sysName[isys[m]],"lp");

    atitle = Form("gpp%d",m);
    gpp[m]->SetName(atitle);
    gpp[m]->SetLineColor(icol[isys[m]]);
    gpp[m]->SetMarkerStyle(33);
    gpp[m]->SetMarkerColor(icol[isys[m]]);
    gpp[m]->SetTitle("pi+; y_lab; <Px> (MeV/c)");
    mgpi->Add(gpp[m],"ip");
    mgpi->Add(gpp[m],"ip");
    aLeg0->AddEntry(gpp[m],"#pi^{+}   "+sysName[isys[m]],"lp");
  }

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mgpdt -> Draw("a");
  aLeg1->Draw();


  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mgpi  ->Draw("a");
  aLeg0->Draw();
}

UInt_t pxbooking()  // used by meanPx()
{
  //----- booking

  for(Int_t m = m_bgn; m < m_end; m++){

    for(UInt_t i = 0; i < nbin; i++){

      Double_t yL = y_min[0] + y_binx*i;
      Double_t yU = yL + y_binx;

      TString hname = Form("hptpr%d%d",m,i);
      TString htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpr[m][i] = new TH1D(hname,"Proton :  "+htitle,pt_nbin, pt_prmin, pt_prmax);
      hptpr[m][i] ->SetLineColor(icol[isys[m]]);

      hname = Form("hptdt%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptdt[m][i] = new TH1D(hname,"Deuteron :"+htitle,pt_nbin, pt_dtmin, pt_dtmax);
      hptdt[m][i] ->SetLineColor(icol[isys[m]]);

      hname = Form("hpttr%d%d",m,i);
     htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hpttr[m][i] = new TH1D(hname,"Triton :  "+htitle,pt_nbin, pt_trmin, pt_trmax);
      hpttr[m][i] ->SetLineColor(icol[isys[m]]);

      hname = Form("hptpm%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpm[m][i] = new TH1D(hname,"Pi- :     "+htitle,pt_nbin, pt_pimin, pt_pimax);
      hptpm[m][i] ->SetLineColor(icol[isys[m]]);

      hname = Form("hptpp%d%d",m,i);
      htitle= Form("fRapidity>=%f&&fRapidity<%f; Pt(MeV/c); dN/dPtdy",yL,yU);
      hptpp[m][i] = new TH1D(hname,"Pi+ :     "+htitle,pt_nbin, pt_pimin, pt_pimax);
      hptpp[m][i] ->SetLineColor(icol[isys[m]]);

    }
  }
  return 1;
}



void PxDistribution(UInt_t nplot)  // used by meanPx()
{
  pt_prmin  =  -pt_prmax;
  pt_dtmin  =  -pt_dtmax;
  pt_trmin  =  -pt_trmax;
  pt_pimin  =  -pt_pimax;

  pxbooking();
  
  for(Int_t m = m_bgn; m < m_end; m++){

    Int_t nEntry = SetBranch(m);

    cout << " Number of events " << nEntry << endl;


    for(Int_t i = 0; i < nEntry; i++){
      aArray->Clear();

      rChain[m]->GetEntry(i);

      if(ntrack[2] == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto rpf   = aPart->GetReactionPlaneFlag();
	auto pid   = aPart->GetPID();

	if( rpf == 0)  continue;
	if( pid > 2000 && (rpf == 110 || rpf == 210 ) ) continue;

	auto charg = aPart->GetCharge();
	auto rapid = aPart->GetRapidity();
	auto vp    = aPart->GetRotatedMomentum();
	auto dltphi= aPart->GetAzmAngle_wrt_RP();;
	vp.SetPhi(dltphi);

	auto px    = vp.X();
	//	auto px    = vp.Pt();

	//rapid = (rapid - ycm[isys[m]]) / ybeam_cm[isys[m]];

	for(UInt_t k = 0; k < nbin; k++){
	  
	  if(rapid < y_min[0] + y_binx * (Double_t)(k+1)) {

	    if(pid == partid[0] && charg < 0)
	      hptpm[m][k]->Fill(px);
	    else if(pid == partid[1] && charg > 0)
	      hptpp[m][k] ->Fill(px);
	    else if(pid == partid[2])
	      hptpr[m][k] ->Fill(px);
	    else if(pid == partid[3])
	      hptdt[m][k] ->Fill(px);
	    else if(pid == partid[4])
	      hpttr[m][k] ->Fill(px);

	    break;
	    
	  }
	}

      }
    }

  }
  
  if(nplot == 0) hpt_plot(); 
  
}
 

void hpt_plot()
{

  for(Int_t m = m_bgn; m < m_end; m++){

    Int_t nEntry = rChain[m]->GetEntries();

    auto sclp  = 1./nEntry * (Double_t)pt_nbin / (pt_prmax - pt_prmin) * y_binx;
    auto scld  = 1./nEntry * (Double_t)pt_nbin / (pt_dtmax - pt_dtmin) * y_binx;
    auto sclt  = 1./nEntry * (Double_t)pt_nbin / (pt_trmax - pt_trmin) * y_binx;
    auto sclpi = 1./nEntry * (Double_t)pt_nbin / (pt_pimax - pt_pimin) * y_binx;

    for(UInt_t j = 0; j < nbin; j++){
      hptpr[m][j]->Scale(sclp);
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
	hptpr[1][j]->Draw(iopt[0]);
	hptpr[0][j]->Draw(iopt[1]);

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
      hptdt[m][j]->Scale(scld);
      hptdt[m][j]->Draw(iopt[m]);
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
      hpttr[m][j]->Scale(sclt);
      hpttr[m][j]->Draw(iopt[m]);
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
      hptpp[m][j]->Draw(iopt[m]);
    }
  }
}


void JPS()
{
  UInt_t m = 0;

  auto hphir0 = new TH1D("hphir0","fRapidity<0.1; #Delta#phi; dN/d#phi",100,-180.,180);
  rChain[m]->Project("hphi0","fdeltphi","fRapidity<0.1");

  auto hphir1 = new TH1D("hphir1","Rapidity(0.36 ~ 0.46); #Delta#phi; dN/d#phi",100,-180.,180);
  rChain[m]->Project("hphi1","fdeltphi","fRapidity>0.36&&fRapidity<0.46");

  auto hphir2 = new TH1D("hphir2","fRapidity>0.7; #Delta#phi; dN/d#phi",100,-180.,180);
  rChain[m]->Project("hphi2","fdeltphi","fRapidity>0.7");



}


void PlotPtDependence(UInt_t selid = 2)       //%% Executable :
{

  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  Int_t RPflag = 0;
  if(selid >= 2) RPflag = 10;

  cout << " Particle " << partname[selid] << endl;
  
  // Rapidity binning
  TArrayD arr_range1(ybin1, yrange1);
  TArrayD arr_range2(ybin2, yrange2);

  // PT binning
  Double_t pt_max = 800.;
  UInt_t nbin1 = 10; //16
  UInt_t nbin2 = 8; //10
  Double_t dpt1 = pt_max/(Double_t)nbin1;
  Double_t dpt2 = pt_max/(Double_t)nbin2;


  std::cout << " Rapidity binning " << ybin1 << std::endl;
  TString rangeLabel[ybin1];
  for(UInt_t i = 0; i < ybin1; i++ ){
    if( i == 0 )
      rangeLabel[0] = Form(" y < %f ",yrange1[0]);
    else if ( i == ybin1 -1 )
      rangeLabel[i] = Form("%f <= y "    ,yrange1[ybin1-1]);
    else 
      rangeLabel[i] = Form("%f <= y < %f",yrange1[i-1],yrange1[i]);
  }

  
  auto hmass = new TH2D("hmass",";P/Q; Mass [MeV/c]",200,0.,2500.,200,0.,7000);

  for(Int_t m = m_bgn; m < m_end; m++){

    hmass->SetName(Form("hmass_%d",m));
    hmass->Reset();

    TH2D *hypt[ybin1];
    TH2D *hypt2[ybin1];
    TH2D *hyphi1[ybin1];
    TH2D *hyphi2[ybin2];
    TH2D *hyptphi1[ybin1][nbin1];
    TH2D *hyptphi2[ybin2][nbin2];

    // Lorentz Transform    
    TVector3 boostVec = LorentzBoost(4);
  
    for(UInt_t kn = 0; kn < ybin1; kn++){ 

      TString sname = rangeLabel[kn];
      hypt[kn]    = new TH2D((TString)Form("hypt_%d",kn),    sname+"; Rapidity; Pt [MeV/c]"  , 200, -0.4, 0.45, 200, 0., 800);
      hyphi1[kn]  = new TH2D((TString)Form("hyphi1_%d",kn),  sname+"; #Phi"   , 100, -3.2,  3.2, 200, 0., 800);

      for(UInt_t pn = 0; pn < nbin1; pn++)
	hyptphi1[kn][pn] = new TH2D((TString)Form("hyptphi1_%d%d",kn,pn),"pt; #Delta #Phi" , 100, -3.2, 3.2, 200, 0., 1000.); 
    }


    for(UInt_t kn = 0; kn < ybin2; kn++){ 
      TString sname = rangeLabel[kn];
      hypt2[kn]   = new TH2D((TString)Form("hypt2_%d",kn),   sname+"; Rapidity; Pt [MeV/c]"  , 200, -0.4, 0.45, 200, 0., 800);
      hyphi2[kn]  = new TH2D((TString)Form("hyphi2_%d",kn),  sname+";2*#Phi"  ,  50,  0.,   3.2, 200, 0., 800);

      for(UInt_t pn = 0; pn < nbin2; pn++)
	hyptphi2[kn][pn] = new TH2D((TString)Form("hyptphi2_%d%d",kn,pn),"pt; #Delta 2*#Phi", 50, 0., 3.2,   200, 0., 1000.); 
    }


    Int_t nevt = SetBranch(m);
    cout << " Number of events " << nevt << endl;

    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      if(ntrack[2] == 0) continue;

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto rpf   = aPart->GetReactionPlaneFlag();
	auto pid   = aPart->GetPID();
	auto bmass = aPart->GetBBMass();

	if( pid == partid[selid] && aPart->GetCharge() == pcharge && aPart->GetIndividualRPAngle() > -9 && bmass <= cutbmass[selid] ){

	  if( pid > 2000 && (rpf == 110 || rpf == 210 ) ) continue;

	  auto pt    = aPart->GetRotatedMomentum().Pt();
	  auto dphi  = aPart->GetAzmAngle_wrt_RP();
	  auto rapid = GetRapidity_cm(aPart->GetMomentum(), aPart->GetMass(), -boostVec);


	  hmass->Fill( aPart->GetRotatedMomentum().Mag(), bmass);


	  UInt_t irapid = ybin1 - 1;
	  for( UInt_t i = 0; i < ybin1; i++){
	    if(rapid < yrange1[i]){
	      irapid = i;
	      break;
	    }
	  }

	  hypt[irapid]->Fill(rapid, pt);
	  hyphi1[irapid]->Fill(dphi, pt);
	  UInt_t ipt = nbin1 - 1;
	  for(UInt_t i = 0; i < nbin1; i++){
	    if( pt < dpt1*(i+1)) {
	      ipt = i;
	      break;
	    }
	  }
	  hyptphi1[irapid][ipt]->Fill(dphi, pt);
	

	  
	  irapid = ybin2 - 1;
	  for( UInt_t i = 0; i < ybin2; i++){
	    if(rapid < yrange2[i]){
	      irapid = i;
	      break;
	    }
	  }

	  hypt2[irapid]->Fill(rapid, pt);
	  hyphi2[irapid]->Fill(abs( TVector2::Phi_mpi_pi(2.*dphi)), pt );

	  ipt = nbin2 - 1;
	  for(UInt_t i = 0; i < nbin2; i++){
	    if( pt < dpt2*(i+1)) {
	      ipt = i;
	      break;
	    }
	  }
	  hyptphi2[irapid][ipt]->Fill(abs( TVector2::Phi_mpi_pi(2.*dphi) ), pt);
	}
      }
    }
    
  
    TString fName = "YPT" + sysName[isys[m]] + "_" + partname[selid]+".root";
    gSystem->cd("data");
    auto GraphSave = new TFile(fName,"recreate");
    
    TGraphErrors *gv_v1 = new TGraphErrors();
    gv_v1->SetName("gv_v1");
    TGraphErrors *gv_v2 = new TGraphErrors();
    gv_v2->SetName("gv_v2");

    TGraphErrors *gPt_v1[ybin1];
    TGraphErrors *gPt_v2[ybin1];

    for(UInt_t kn = 0; kn < ybin1 ; kn++){      
      gPt_v1[kn] = new TGraphErrors();
      gPt_v1[kn]->SetName((TString)Form("gPt_v1%d",kn));
      TString sname = partname[selid]+"; Pt [MeV/c]; v1";
      gPt_v1[kn]->SetTitle(sname);
    }

    for(UInt_t kn = 0; kn < ybin2 ; kn++){      
      gPt_v2[kn] = new TGraphErrors();
      gPt_v2[kn]->SetName((TString)Form("gPt_v2%d",kn));
      TString sname = partname[selid]+"; Pt [MeV/c]; v2";
      gPt_v2[kn]->SetTitle(sname);
    }

    ic++; 
    cc[ic]   = new TCanvas("dphi1","dphi1",2000,1200);
    cc[ic]->Divide(nbin1, ybin1);

    cc[ic+1] = new TCanvas("dphi2","dphi2",1400,1200);
    cc[ic+1]->Divide(nbin2, ybin2);

    cc[ic+2] = new TCanvas(Form("cc%d",ic+2),Form("cc%d",ic+2), 1000, 1200);
    cc[ic+2]->Divide(2, ybin1);

    //************************************************** 
    std::cout << " ---- Resutls ---------------------" << std::endl;

    UInt_t kl = 0;
    UInt_t id1 = 0;
    UInt_t id2 = 0;
    Double_t para[6];
    
    for(UInt_t kn = 0; kn < ybin1; kn++){
      
      Double_t rapm  = hypt[kn]->GetMean(1);
      Double_t rape  = hypt[kn]->GetStdDev(1);

      // v1 vs rapidity
      cc[ic+2]->cd(2*(kn+1)-1);
      auto hyphi = (TH1D*)hyphi1[kn]->ProjectionX("",0,-1,"eo");

      Double_t corr[2]={mcos1[isys[m]], mcos1e[isys[m]]};
      GetFittingParameters(*hyphi, para, corr);
      
      gv_v1->SetPoint( kl,     rapm, para[4]);
      gv_v1->SetPointError( kl, rape, para[5]);
      kl++;

      // pt dependence 
      UInt_t il = 0; 
      for(UInt_t jn = 0; jn < nbin1; jn++){

	id1++; cc[ic]->cd(id1); 

	if( hyptphi1[kn][jn]->GetEntries() > 0 ) {	

	  Double_t ptc  = hyptphi1[kn][jn]->GetMean(2);
	  Double_t ptce = hyptphi1[kn][jn]->GetStdDev(2);
	    
	  auto hypt = (TH1D*)hyptphi1[kn][jn]->ProjectionX();
	  GetFittingParameters(*hypt, para, corr);
	    

	  gPt_v1[kn]->SetPoint(il, ptc, para[4]);
	  gPt_v1[kn]->SetPointError(il, ptce, para[5]);
	    
	  il++;
	}
      }
      gPt_v1[kn]->Write();
    }

    kl = 0;
    for(UInt_t kn = 0; kn < ybin2; kn++){

      Double_t rapm  = hypt2[kn]->GetMean(1);
      Double_t rape  = hypt2[kn]->GetStdDev(1);

      // v2 vs rapidity
      cc[ic+2]->cd(2*(kn+1));

      auto hyphi = (TH1D*)hyphi2[kn]->ProjectionX("",0,-1,"eo");
      
      Double_t corr[2] = {mcos2[isys[m]], mcos2e[isys[m]]};
      GetFittingParameters(*hyphi, para, corr);

      gv_v2->SetPoint( kl,    rapm, para[4]);
      gv_v2->SetPointError( kl, rape, para[5]);
      
      kl++;
    

      UInt_t il = 0; 
      for(UInt_t jn = 0; jn < nbin2; jn++){
	
	id2++; cc[ic+1]->cd(id2); 
	if( hyptphi2[kn][jn]->GetEntries() > 0 ){

	  Double_t ptc  = hyptphi2[kn][jn]->GetMean(2);
	  Double_t ptce = hyptphi2[kn][jn]->GetStdDev(2);

	  auto hypt = (TH1D*)hyptphi2[kn][jn]->ProjectionX();

	  Double_t corr[2] = {mcos2[isys[m]], mcos2e[isys[m]]};
	  GetFittingParameters(*hypt, para, corr);

	  gPt_v2[kn]->SetPoint(il, ptc, para[4]);
	  gPt_v2[kn]->SetPointError(il, ptce, para[5]);
	  
	  il++;
	}
      }
      gPt_v2[kn]->Write();
    }

    gv_v1->Write();
    gv_v2->Write();

    ic+=3;


    //plotting
    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    hmass->Draw("colz");

    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
    cc[ic]->Divide(ybin1,1);
    for(UInt_t kn = 0; kn < ybin1; kn++){
      cc[ic]->cd(id); id++;
      gPt_v1[kn]->Draw("ALP");
    }

    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
    cc[ic]->Divide(ybin2,1);
    for(UInt_t kn = 0; kn < ybin2; kn++){
      cc[ic]->cd(id); id++;
      gPt_v2[kn]->Draw("ALP");
    }

    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,400);
    cc[ic]->Divide(ybin2,1);
    for(UInt_t kn = 0; kn < ybin2; kn++){
      cc[ic]->cd(id); id++;
      hypt2[kn]->Draw("colz");
    }


    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
    gv_v1->Draw("ALP");

    auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
    auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
    LineV->SetLineStyle(3);
    LineH->SetLineStyle(3);
    LineV->Draw();
    LineH->Draw();


    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
    gv_v2->Draw("ALP");


    GraphSave->WriteObject(&arr_range1, "v1_yrange");
    GraphSave->WriteObject(&arr_range2, "v2_yrange");
  }

  
  gSystem->cd("..");
}

void YDependece(UInt_t hrm=1)
{
  cout << " YDependece v"<< hrm  << endl;

  UInt_t   selid[3] = {2,3,4};

  for(Int_t m = m_bgn; m < m_end; m++){

    vector< vector< vector<Double_t> > > bphi;
    vector< vector< vector<Double_t> > > xbin;
    bphi.resize(3);
    xbin.resize(3);

    for(UInt_t im = 0; im < 3; im++){
      bphi[im].resize(nbin);
      xbin[im].resize(nbin);
    }


    Int_t nevt = SetBranch(m);
    cout << " Number of events " << nevt << endl;

    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	auto pid   = aPart->GetPID();
	auto rapid = aPart-> GetRapidity();
	rapid = (rapid - ycm[isys[m]])/ybeam_cm[isys[m]];

	if(  aPart->GetBestTrackFlag() > 0  ){

	  UInt_t k = 0; 
	  while( k < nbin ){
	    if( rapid < y_min[0] + y_binx*(Double_t)(k+1) )
	      break;
	    
	    k++;
	  }

	  if(k < nbin){
	    for(UInt_t j = 0; j < 3; j++){
	      if( pid == partid[selid[j]] ){
		bphi[j][k].push_back(aPart->GetAzmAngle_wrt_RP());
		xbin[j][k].push_back(rapid);
		
	      }
	    }
	  }
	}
      }
    }

    Double_t xval [4][nbin];
    Double_t xvale[4][nbin];
    Double_t yval [4][nbin];
    Double_t yvale[4][nbin];

    for(UInt_t j = 0; j < 3; j++) {
      for(UInt_t jn = 0; jn < nbin; jn++) {

	auto getx = vMean(xbin[j][jn]);
	auto gety = vn(hrm, bphi[j][jn]);

	xval [j][jn] = getx[0];
	xvale[j][jn] = getx[1];

	yval [j][jn] = gety[0];
	yvale[j][jn] = gety[1];

      }
    }

    Double_t yratio [2][nbin];
    Double_t yratioe[2][nbin];

    // reversed
    for(UInt_t jn = 0; jn < nbin; jn++) {
      xval [3][jn] = -xval[0][jn];
      xvale[3][jn] = xvale[0][jn];

      yval [3][jn] =  yval[0][jn] * pow(-1, hrm);
      yvale[3][jn] = yvale[0][jn];

      if(yval[0][jn] != 0){
	yratio [0][jn] = abs(yval[1][jn] / yval[0][jn]);
	yratioe[0][jn] = pow(yvale[1][jn],2)/pow(yval[1][jn],2) + pow(yvale[0][jn],2)/pow(yval[0][jn],2) ;
	yratioe[0][jn] = abs(yratio [0][jn]) * sqrt(yratioe[0][jn]);

	// cout << " yratio 0 " << jn << " " << setw(10)
	//      << yval[1][jn] << " / " << yval[0][jn] << " = " << yratio [0][jn]
	//      << " +- " << yratioe[0][jn] << " " << xval[0][jn] <<endl;
      }

      if(yval[0][jn] != 0){
	yratio [1][jn] = abs(yval[2][jn] / yval[0][jn]);
	yratioe[1][jn] = pow(yvale[2][jn],2)/pow(yval[2][jn],2) + pow(yvale[0][jn],2)/pow(yval[0][jn],2) ;
	yratioe[1][jn] = abs(yratio [1][jn]) * sqrt(yratioe[1][jn]);

	// cout << xval[0][jn] << " " <<  xval[2][jn] << endl; 
	// cout << " yratio 1 " << jn << " : " << yratio[1][jn] << " +- " << yratioe[1][jn] << endl;
      }
    }
    

    const UInt_t nsl = 10;
    UInt_t iisl = 2;
    auto gr0 = new TGraphErrors(nsl, &xval[0][iisl], &yval[0][iisl], &xvale[0][iisl], &yvale[0][iisl]);
    gr0->SetName("gr0");
    gr0->SetLineColor(2);
    gr0->SetMarkerStyle(20);
    gr0->SetMarkerColor(2);

    auto gr1 = new TGraphErrors(nsl, &xval[1][iisl], &yval[1][iisl], &xvale[1][iisl], &yvale[1][iisl]);
    gr1->SetName("gr1");
    gr1->SetLineColor(4);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerColor(4);

    auto gr2 = new TGraphErrors(nsl, &xval[2][iisl], &yval[2][iisl], &xvale[2][iisl], &yvale[2][iisl]);
    gr2->SetName("gr2");
    gr2->SetLineColor(8);
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerColor(8);

    auto *rv_gr0 = new TGraphErrors(nsl, &xval[3][iisl], &yval[3][iisl], &xvale[3][iisl], &yvale[3][iisl]);
    rv_gr0->SetName("rv_gv0");
    rv_gr0->SetLineColor(2);
    rv_gr0->SetMarkerStyle(24);
    rv_gr0->SetMarkerColor(2);
    rv_gr0->SetLineStyle(3);


    auto mg = new TMultiGraph();
    TString aTitle = Form("; Ycm/Ycm_beam; v%d(a.u.)",hrm);
    mg->SetTitle(aTitle);
    mg->SetName("mg");
    mg->Add(gr0,"lp");
    mg->Add(gr1,"lp");
    mg->Add(gr2,"lp");
    mg->Add(rv_gr0,"lp");

    if(hrm == 2)
      mg->SetMaximum(0.001);

    // Legend                                                                                                                               

    TString sconf = ""; //Form("%4.2f < y < %4.2f",rapid_min,rapid_max);
    auto aLeg = new TLegend(0.1,0.7,0.35,0.9,sconf);
    aLeg->AddEntry(gr0,partname[selid[0]],"lp");
    aLeg->AddEntry(gr1,partname[selid[1]],"lp");
    aLeg->AddEntry(gr2,partname[selid[2]],"lp");
    aLeg->AddEntry(rv_gr0,partname[selid[0]]+" reversed","lp");

    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    //    gr0->Draw();

    mg  ->Draw("a");
    aLeg->Draw();


    iisl = 2;
    auto gt0 = new TGraphErrors(nsl, &xval[0][iisl], &yratio[0][iisl], &xvale[0][iisl], &yratioe[0][iisl]);
    gt0->SetName("gt0");
    gt0->SetLineColor(4);
    gt0->SetMarkerStyle(21);
    gt0->SetMarkerColor(4);

    auto gt1 = new TGraphErrors(nsl, &xval[0][iisl], &yratio[1][iisl], &xvale[0][iisl], &yratioe[1][iisl]);
    gt1->SetName("gt0");
    gt1->SetLineColor(8);
    gt1->SetMarkerStyle(22);
    gt1->SetMarkerColor(8);

    auto mgr = new TMultiGraph();
    aTitle = Form("; Ycm/Ycm_beam; v%d /v%d(proton)",hrm, hrm);
    mgr->SetTitle(aTitle);
    mgr->SetName("mgr");
    mgr->Add(gt0,"lp");
    mgr->Add(gt1,"lp");

    auto aLegr = new TLegend(0.1,0.7,0.35,0.9,"");
    aLegr->AddEntry(gt0,partname[selid[1]],"lp");
    aLegr->AddEntry(gt1,partname[selid[2]],"lp");


    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    mgr->Draw("a");
    aLegr->Draw();

  }
}




void FlatteningCheck()                        //%% Executable : 
{
  //----- Parametres
  Int_t ntrack[7];

  //----- Booking
  TH2D *hphitheta[4];
  TH2D *hphimtrck[4];

  for(Int_t m = m_bgn; m < m_end; m++){
    TString hname = Form("hphitheta%d",m);
    hphitheta[m] = new TH2D(hname, sysName[isys[m]]+"; #Theta ; #Phi",100,0,0.8, 100,-3.2, 3.2);

    hname = Form("hphimtrck%d",m);
    hphimtrck[m] = new TH2D(hname, sysName[isys[m]]+"; Number of Track ; #Phi",35,0,35, 100,-3.2, 3.2);
  }


  //----- Filling
  for(Int_t m = m_bgn; m < m_end; m++){


    Int_t nEntry = SetBranch(m);

    for(Int_t i = 0; i < nEntry; i++){
      aArray->Clear();

      rChain[m]->GetEntry(i);

      TIter next(aArray);
      STParticle *aPart = NULL;

      while( (aPart = (STParticle*)next()) ) {

	// auto pid   = aPart->GetPID();
	// auto charg = aPart->GetCharge();
	// auto rapid = aPart->GetRapidity();
	// auto vp    = aPart->GetFlattenMomentum();
	// auto dltphi= aPart->GetAzmAngle_wrt_RP();;
	auto phi   = aPart->GetFlattenMomentum().Phi();
	auto theta = aPart->GetFlattenMomentum().Theta();
	auto flag  = aPart->GetReactionPlaneFlag();

	//	if(flag > 110 ){
	if(flag >= selReactionPlanef ){
	  hphitheta[m]->Fill( theta, phi );
	  hphimtrck[m]->Fill( ntrack[5], phi ); 
	}
      }
    }
  }


  //----- Drawing 

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  cc[ic]->Divide(2,m_end);

  UInt_t id = 1;

  for(Int_t m = m_bgn; m < m_end; m++){
    cc[ic]->cd(id); id++;
    hphitheta[m]->Draw("colz");

    cc[ic]->cd(id); id++;
    hphimtrck[m]->Draw("colz");
  }

}


void CorrectedPsi()                           //%%
{
  auto hrt  = new TH1D("hrt" ,"Original; #Psi",60,-3.15,3.15);
  auto hrc  = new TH1D("hrc" ,"ReCentering; #Psi",60,-3.15,3.15);
  auto hfc  = new TH1D("hfc" ,"RC + Shifting; #Psi",60,-3.15,3.15);


  rChain[0]->Project("hrt" , "TVector2::Phi_mpi_pi(unitP2_rot.Phi())");
  rChain[0]->Project("hrc" , "unitP_rc.Phi()");
  rChain[0]->Project("hfc" , "unitP_fc.Phi()");


  //----- Drawing                                                                                                                          
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  hrt->SetLineColor(2);
  hrc ->SetLineColor(8);
  hfc ->SetLineColor(4);

  hrt->Draw("e");
  hrc->Draw("samee");
  hfc->Draw("samee");


  auto aLeg = new TLegend(0.15,0.7,0.45,0.9,"");
  aLeg->AddEntry(hrt ,"No Collection ","lp");
  aLeg->AddEntry(hrc ,"ReCentering","lp");
  aLeg->AddEntry(hfc ,"ReCentering & Shifting","lp");

  aLeg->Draw();
}


void PlotSubEvent(Double_t ml, Double_t mu)   //%%
{

  Double_t mlt[] = {ml, mu};
  //    Double_t mlt[2] = {0., 8.};
  //Double_t mlt[2] = {8., 16.};
  // Double_t mlt[2] = {16., 24.};
  // Double_t mlt[2] = {24., 32.};
  // Double_t mlt[2] = {32., 40.};
  // Double_t mlt[2] = {40., 100.};

  TCut mcrot = Form("ntrack[4]>%f&&ntrack[4]<%f",mlt[0]*2.,mlt[1]*2.);
  TCut mc1r  = Form("mtrack_1>%f&&mtrack_1<%f"  ,mlt[0],mlt[1]);
  TCut mc2r  = Form("mtrack_2>%f&&mtrack_2<%f"  ,mlt[0],mlt[1]);

  TString sname;
  
  sname = mcrot.GetTitle();
  auto *hrotx = new TH1D("hrotx","All   "+sname+ ";Qx" ,100,-12.,12.);
  auto *h1rx  = new TH1D("h1rx", "sub_1 "+sname+ ";Qx" ,100,-12.,12.);
  auto *h2rx  = new TH1D("h2rx", "sub_2 "+sname+ ";Qx" ,100,-12.,12.);
  auto *hroty = new TH1D("hroty","All   "+sname+ ";Qy" ,100,-12.,12.);
  auto *h1ry  = new TH1D("h1ry", "sub_1 "+sname+ ";Qy" ,100,-12.,12.);
  auto *h2ry  = new TH1D("h2ry", "sub_2 "+sname+ ";Qy" ,100,-12.,12.);
 

  rChain[0]->Project("hrotx","unitP2_rot.X()",mcrot);
  rChain[0]->Project("h1rx" ,"unitP_1r.X()"  ,mc1r);
  rChain[0]->Project("h2rx" ,"unitP_2r.X()"  ,mc2r);
					      
  rChain[0]->Project("hroty","unitP2_rot.Y()",mcrot);
  rChain[0]->Project("h1ry" ,"unitP_1r.Y()"  ,mc1r);
  rChain[0]->Project("h2ry" ,"unitP_2r.Y()"  ,mc2r);


  //----- Drawing                                                                                                                          
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  hrotx->SetLineColor(2);
  hrotx->SetNormFactor(1);
  h1rx ->SetLineColor(4);
  h1rx ->SetNormFactor(1);
  h2rx ->SetLineColor(6);
  h2rx ->SetNormFactor(1);

  h1rx->Draw();
  h2rx->Draw("same");
  hrotx->Draw("same");

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  hroty->SetLineColor(2);
  hroty->SetNormFactor(1);
  h1ry ->SetLineColor(4);
  h1ry ->SetNormFactor(1);
  h2ry ->SetLineColor(6);
  h2ry ->SetNormFactor(1);

  h1ry->Draw();
  h2ry->Draw("same");
  hroty->Draw("same");
}



void PlotNeuLANDPsi()                         //%%
{
  TH2D *hACCp[4];
  TH2D *hnlACCp[4];
  TH2D *hnlACCn[4];
  TH1D *hfcrn[4][4];
  TH1D *hfcrp[4][4];
  TH1D *hfcrnd[4][4];

  TCut ncCut[4];
  ncCut[0]="";
  ncCut[1]="ncPID==2112";
  ncCut[2]=ncCut[1]&&"ncRapidity<=0.37";
  ncCut[3]=ncCut[1]&&"ncRapidity>0.37";

  TCut pcCut[4];
  pcCut[0]="";
  pcCut[1]="ncPID==2212";
  pcCut[2]=pcCut[1]&&"ncRapidity<=0.37";
  pcCut[3]=pcCut[1]&&"ncRapidity>0.37";

  //----- Booking            
  TString hname;
  for(Int_t m = m_bgn; m < m_end; m++){
    hname = Form("hnlACCn%d",m);
    hnlACCn[m] = new TH2D(hname, hname,  200, 0., 1., 200., 0., 500.);
    hnlACCn[m]->SetMarkerColor(4);
    rChain[m]->Project(hname,"ncP.Pt():ncRapidity","ncPID==2112");

    hname = Form("hnlACCp%d",m);
    hnlACCp[m] = new TH2D(hname, hname,  200, 0., 1., 200., 0., 500.);
    hnlACCp[m]->SetMarkerColor(2);
    rChain[m]->Project(hname,"ncP.Pt():ncRapidity","ncPID==2212");

    hname = Form("hACCp%d",m);
    hACCp[m] = new TH2D(hname, hname,  200, 0., 1., 200., 0., 500.);
    hACCp[m]->SetTitle("; Rapidity; Pt [GeV/c]");
    rChain[m]->Project(hname,"fRotatedP3.Pt():fRapidity","fPID==2212&&fReactionPlanef>100");

    for(UInt_t k = 0; k < 4; k++){
      hname = Form("hfcrn%d%d",m,k);
      hfcrn[m][k] = new TH1D(hname, hname+ncCut[k].GetTitle(), 60, -3.15, 3.15);
      rChain[m]->Project(hname,"TVector2::Phi_mpi_pi(unitP_fc.Phi()-ncP.Phi())",ncCut[k]);

      hname = Form("hfcrnd%d%d",m,k);
      hfcrnd[m][k] = new TH1D(hname, hname+ncCut[k].GetTitle(), 60, -180., 180.);
      rChain[m]->Project(hname,"TVector2::Phi_mpi_pi(unitP_fc.Phi()-ncP.Phi())*180./TMath::Pi()",ncCut[k]);

      hname = Form("hfcrp%d%d",m,k);
      hfcrp[m][k] = new TH1D(hname, hname+pcCut[k].GetTitle(), 60,  -180., 180.);
      rChain[m]->Project(hname,"TVector2::Phi_mpi_pi(unitP_fc.Phi()-ncP.Phi())*180./TMath::Pi()",pcCut[k]);
    }

    if(kFALSE){
      ic++; id = 1;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

      hACCp[m]->Draw("colz");
      hnlACCn[m]->Draw("same");

      auto PRLabel = new TLatex(0.7,50,"TPC proton");
      auto NLLabel = new TLatex(0.2,100,"Neutron");
      NLLabel->SetTextColor(0);
      PRLabel->SetTextColor(0);
      NLLabel->Draw();
      PRLabel->Draw();


      ic++; id = 1;
      cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
      cc[ic]->Divide(4,2);

      for(UInt_t k = 0; k < 4; k++){
	cc[ic]->cd(id); id++;
	hfcrn[m][k]->SetLineColor(4);
	hfcrn[m][k]->Draw("e");
      }
      for(UInt_t k = 0; k < 4; k++){
	cc[ic]->cd(id); id++;
	hfcrp[m][k]->SetLineColor(2);
	hfcrp[m][k]->Draw("e");
      }
    }

    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    hfcrnd[m][2]->SetLineColor(2);
    hfcrnd[m][2]->SetMarkerColor(2);
    hfcrnd[m][2]->SetNormFactor(60);
    hfcrnd[m][2]->SetMarkerStyle(20);
    hfcrnd[m][2]->SetTitle(";#Psi;1/NdN/d#Psi");
    hfcrnd[m][3]->SetLineColor(4);
    hfcrnd[m][3]->SetMarkerColor(4);
    hfcrnd[m][3]->SetNormFactor(60);
    hfcrnd[m][3]->SetMarkerStyle(20);
    hfcrnd[m][3]->SetTitle(";#Psi;1/NdN/d#Psi");
    hfcrp[m][3]->SetLineColor(8);
    hfcrp[m][3]->SetNormFactor(60);
    hfcrp[m][3]->SetMarkerColor(8);
    hfcrp[m][3]->SetMarkerStyle(20);;

    hfcrp[m][3]->Draw("e");
    hfcrnd[m][3]->Draw("samee");
    hfcrnd[m][2]->Draw("samee");



    auto aLeg = new TLegend(0.65,0.7,0.9,0.9,"");
    aLeg->AddEntry(hfcrnd[m][2] ,"Neutron hit y < y_cm","lp");
    aLeg->AddEntry(hfcrnd[m][3] ,"Neutron hit y > y_cm","lp");
    aLeg->AddEntry(hfcrp[m][3] ,"Proton on NL y > y_cm","lp");

    aLeg->Draw();


  }
}


//--------------------------------------------//%% Executable : 
void PlotNeuLANDProperty(UInt_t iout)           //%% Executable : 
{
  std::cout << " Executing : PlotNeuLANDProperty(" << iout << ")" << std::endl;

  //----- Parametres                                                                                                                       
  UInt_t nybin  = 100;
  TFile* hout;
  UInt_t total_neut = 0;
 
  //----- Canvas
  if( iout == 1 ){
    ic++;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  }

  //----- Booking
  TH1D* hmult = new TH1D("hmult0",";Multiplicity ", 12,0.,12.);
  TH2D* haccp = new TH2D("haccp0","; Rapidity ; Pt [MeV/c]",nybin, -0.4, 0.5, 100,   0., 800);


  //----- Event loop
  for(Int_t m = m_bgn; m < m_end; m++){

    total_neut = 0;
    TVector3 boostVec = LorentzBoost(4);

  //----- Output file
    if( iout == 1 ) {
      TString fName = "NLdbProp.cm" + sysName[isys[m]] + ".root";
      gSystem->cd("data");
      hout = new TFile(fName,"recreate");
    }
    
    hmult->Reset();
    haccp->Reset();

    Int_t nevt = SetBranch(m);
    cout << " Number of events " << nevt << endl;
    
    UInt_t vID = 0;
    if( snbm == 108 )
      vID = 1;

    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      UInt_t nneut = 0;
      TIter next(aNLClusterArray);
      STNeuLANDCluster* aNLClust = NULL;

      while( (aNLClust = (STNeuLANDCluster*)next()) ){

	aNLClust->SetBeamAngle(ProjA/1000., ProjB/1000.);
        auto pid      = aNLClust->GetPID();
        auto pt       = aNLClust->GetP().Pt();
        auto veto_all = aNLClust->GetVetoHitAll(vID);
        auto veto_bar = aNLClust->GetVetoHitOne(vID);
	
	if( pid == 2112 && aNLClust->GetTOF()>40 ){
	  auto rapidity = GetRapidity_cm( aNLClust->GetP(), 939.5731, -boostVec);
	  haccp->Fill(rapidity, pt);
	  nneut++;
	  total_neut++;
	}
      }
      
      hmult->Fill(nneut);
    }
  
    
    //    hmult->Draw();

    if(iout == 1){
      hmult->Write();
      haccp->Write();

      hout->Close();
      gSystem->cd("..");
    }

    gStyle->SetOptStat(0);
    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    haccp->Draw("colz");


    ic++;
    auto cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    cc->SetLogy();
    hmult->Draw();


    std::cout << " Number of neutron " << total_neut << " / " << hmult->GetEntries() 
	      << " = " << Float_t(total_neut/hmult->GetEntries())
	      << " : " << nevt << std::endl;
  }
}

//--------------------------------------------//%% Executable : 
void PlotNeuLANDv1v2()                        //%% Executable : 
{
  //----- Parametres                                                                                                                       
  TFile *hout;

  // Rapidity binning
  TArrayD arr_range1(ybin1, yrange1);
  TArrayD arr_range2(ybin2, yrange2);

  // PT binning                                                                                                                            
  Double_t pt_max = 800.;
  UInt_t nbin1 = 10; //16                                                                                                                  
  UInt_t nbin2 = 8; //10    
  Double_t dpt1 = pt_max/(Double_t)nbin1;
  Double_t dpt2 = pt_max/(Double_t)nbin2;

  std::cout << " Rapidity binning " << ybin1 << std::endl;
  TString rangeLabel[ybin1];
  for(UInt_t i = 0; i < ybin1; i++ ){
    if( i == 0 )
      rangeLabel[0] = Form(" y < %f ",yrange1[0]);
    else if ( i == ybin1 -1 )
      rangeLabel[i] = Form("%f <= y "    ,yrange1[ybin1-1]);
    else
      rangeLabel[i] = Form("%f <= y < %f",yrange1[i-1],yrange1[i]);
  }


  //----- System Loop
  for(Int_t m = m_bgn; m < m_end; m++){

    TH2D *hypt[ybin1];
    TH2D *hypt2[ybin1];
    TH2D *hyphi1[ybin1];
    TH2D *hyphi2[ybin2];
    TH2D *hyptphi1[ybin1][nbin1];
    TH2D *hyptphi2[ybin2][nbin2];

    // Lorentz Transform                                                                                                                    
    TVector3 boostVec = LorentzBoost(4);

    for(UInt_t kn = 0; kn < ybin1; kn++){

      TString sname = rangeLabel[kn];
      hypt[kn]    = new TH2D((TString)Form("hypt_%d",kn),    sname+"; Rapidity; Pt [MeV/c]"  , 200, -0.4, 0.45, 200, 0., 800);
      hyphi1[kn]  = new TH2D((TString)Form("hyphi1_%d",kn),  sname+"; #Phi"   , 100, -3.2,  3.2, 200, 0., 800);

      for(UInt_t pn = 0; pn < nbin1; pn++)
        hyptphi1[kn][pn] = new TH2D((TString)Form("hyptphi1_%d%d",kn,pn),"pt; #Delta #Phi" , 100, -3.2, 3.2, 200, 0., 1000.);
    }

    for(UInt_t kn = 0; kn < ybin2; kn++){
      TString sname = rangeLabel[kn];
      hypt2[kn]   = new TH2D((TString)Form("hypt2_%d",kn),   sname+"; Rapidity; Pt [MeV/c]"  , 200, -0.4, 0.45, 200, 0., 800);
      hyphi2[kn]  = new TH2D((TString)Form("hyphi2_%d",kn),  sname+";2*#Phi"  ,  50,  0.,   3.2, 200, 0., 800);

      for(UInt_t pn = 0; pn < nbin2; pn++)
        hyptphi2[kn][pn] = new TH2D((TString)Form("hyptphi2_%d%d",kn,pn),"pt; #Delta 2*#Phi", 50, 0., 3.2,   200, 0., 1000.);
    }


    
    UInt_t vid = 0;
    if( snbm == 108 )
      vid = 1;


    Int_t nevt = SetBranch(m);

    cout << " Number of events " << nevt << endl;

    for(Int_t i = 0; i < nevt; i++){
      rChain[m]->GetEntry(i);

      TIter next(aNLClusterArray);
      STNeuLANDCluster* aNLClust = NULL;
    
      while( (aNLClust = (STNeuLANDCluster*)next()) ){

	aNLClust->SetBeamAngle(ProjA/1000., ProjB/1000.);
  	
   	auto pid      = aNLClust->GetPID();
   	auto phi      = aNLClust->GetP().Phi(); 
	auto pt       = aNLClust->GetP().Pt();

	auto rapid    = GetRapidity_cm( aNLClust->GetP(), 939.5731, -boostVec);
	auto veto_all = aNLClust->GetVetoHitAll(vid);
	auto veto_bar = aNLClust->GetVetoHitOne(vid);
	auto veto_mid = aNLClust->GetVetoHitMid(vid);


   	auto dphi     = TVector2::Phi_mpi_pi(unitP_fc->Phi() - phi);
  	
	if( pid == 2112 && veto_all == 0){ // neutron
	//   	if( pid == 2112 ){ // neutron

          UInt_t irapid = ybin1 - 1;
          for( UInt_t i = 0; i < ybin1; i++){
            if(rapid < yrange1[i]){
              irapid = i;
              break;
            }
          }

	  hypt[irapid]->Fill(rapid, pt);
          hyphi1[irapid]->Fill(dphi, pt);
          UInt_t ipt = nbin1 - 1;
          for(UInt_t i = 0; i < nbin1; i++){
            if( pt < dpt1*(i+1)) {
              ipt = i;
              break;
            }
          }
          hyptphi1[irapid][ipt]->Fill(dphi, pt);

          irapid = ybin2 - 1;
          for( UInt_t i = 0; i < ybin2; i++){
            if(rapid < yrange2[i]){
              irapid = i;
              break;
            }
          }
          
	  ipt = nbin2 - 1;
          for(UInt_t i = 0; i < nbin2; i++){
            if( pt < dpt2*(i+1)) {
              ipt = i;
              break;
            }
          }
	  hypt2[irapid]->Fill(rapid, pt);
	  hyphi2[irapid]->Fill(abs( TVector2::Phi_mpi_pi(2.*dphi) ), pt);
	  hyptphi2[irapid][ipt]->Fill(abs( TVector2::Phi_mpi_pi(2.*dphi) ), pt);
        }
      }
    }


    // // 132Sn
    //    combineHist(132, hyptphi1, hyptphi2);
    // Sumup edge of acceptance
    // UInt_t highstat = 0;
    // for(UInt_t irapid = 0; i < nbin1l; irapid++){
    //   for(UInt_t ipt = 0; ipt < nbin1; ipt++) {

    // 	if( hyptphi1[irapid][ipt] -> GetEntries() > 100 ) {
    // 	  highstat = ipt;
    // 	}
    // 	else if( highstat < ipt ){
    // 	*hyptphi1[0][0] = *hyptphi1[0][0] + *hyptphi1[0][1];	

    //   }
    // }


    *hyptphi1[0][0] = *hyptphi1[0][0] + *hyptphi1[0][1];
    hyptphi1[0][1]->Reset();

    *hyptphi1[1][1] = *hyptphi1[1][1] + *hyptphi1[1][2];
    hyptphi1[1][2]->Reset();

    *hyptphi1[2][1] = *hyptphi1[2][1] + *hyptphi1[2][3];
    hyptphi1[2][3]->Reset();

    *hyptphi1[3][1] = *hyptphi1[3][1] + *hyptphi1[3][2];
    hyptphi1[3][2]->Reset();

    *hyptphi1[3][3] = *hyptphi1[3][3] + *hyptphi1[3][4];
    hyptphi1[3][4]->Reset();
    
    *hyptphi1[4][5] = *hyptphi1[4][5] + *hyptphi1[4][6];
    hyptphi1[4][6]->Reset();

    *hyptphi1[5][2] = *hyptphi1[5][2] + *hyptphi1[5][3];
    hyptphi1[5][3]->Reset();

    *hyptphi1[5][6] = *hyptphi1[5][6] + *hyptphi1[5][7];
    hyptphi1[5][7]->Reset();
    *hyptphi1[5][6] = *hyptphi1[5][6] + *hyptphi1[5][8];
    hyptphi1[5][8]->Reset();

    *hyptphi1[6][3] = *hyptphi1[6][3] + *hyptphi1[6][4];
    hyptphi1[6][4]->Reset();

    *hyptphi1[6][7] = *hyptphi1[6][7] + *hyptphi1[6][8];
    hyptphi1[6][8]->Reset();
    *hyptphi1[6][7] = *hyptphi1[6][7] + *hyptphi1[6][9];
    hyptphi1[6][9]->Reset();

    *hyptphi1[7][4] = *hyptphi1[7][4] + *hyptphi1[7][5];
    hyptphi1[7][5]->Reset();



    //v2 pt
    *hyptphi2[1][1] = *hyptphi2[1][1] + *hyptphi2[1][2];
    hyptphi2[1][2]->Reset();

    *hyptphi2[2][3] = *hyptphi2[2][3] + *hyptphi2[2][4];
    hyptphi2[2][4]->Reset();

    *hyptphi2[3][5] = *hyptphi2[3][5] + *hyptphi2[3][6];
    hyptphi2[3][6]->Reset();

    *hyptphi2[4][2] = *hyptphi2[4][2] + *hyptphi2[4][3];
    hyptphi2[4][3]->Reset();

    *hyptphi2[4][6] = *hyptphi2[4][6] + *hyptphi2[4][7];
    hyptphi2[4][7]->Reset();




    // output    
    //************************************************** 
    std::cout << " ---- Resutls ---------------------" << std::endl;

    TString fName = "vaNL" + sysName[isys[m]] + "_neutron.root";
    gSystem->cd("data");
    auto GraphSave = new TFile(fName,"recreate");
    
    TGraphErrors *gv_v1 = new TGraphErrors();
    gv_v1->SetName("gv_v1");
    TGraphErrors *gv_v2 = new TGraphErrors();
    gv_v2->SetName("gv_v2");

    TGraphErrors *gPt_v1[ybin1];
    TGraphErrors *gPt_v2[ybin1];

    for(UInt_t kn = 0; kn < ybin1 ; kn++){      
      gPt_v1[kn] = new TGraphErrors();
      gPt_v1[kn]->SetName((TString)Form("gPt_v1%d",kn));
      TString sname = "Neutron; Pt [MeV/c]; v1";
      gPt_v1[kn]->SetTitle(sname);
    }

    for(UInt_t kn = 0; kn < ybin2 ; kn++){      
      gPt_v2[kn] = new TGraphErrors();
      gPt_v2[kn]->SetName((TString)Form("gPt_v2%d",kn));
      TString sname = "Neutron; Pt [MeV/c]; v2";
      gPt_v2[kn]->SetTitle(sname);
    }

    ic++; 
    cc[ic]   = new TCanvas("dphi1","dphi1",2000,1200);
    cc[ic]->Divide(nbin1, ybin1);

    cc[ic+1] = new TCanvas("dphi2","dphi2",1400,1200);
    cc[ic+1]->Divide(nbin2, ybin2);

    cc[ic+2] = new TCanvas(Form("cc%d",ic+2),Form("cc%d",ic+2), 1000, 1200);
    cc[ic+2]->Divide(2, ybin1);


    UInt_t kl = 0;
    UInt_t id1 = 0;
    UInt_t id2 = 0;
    Double_t para[6];
    
    for(UInt_t kn = 0; kn < ybin1; kn++){
      // v1 vs rapidity
      cc[ic+2]->cd(2*(kn+1)-1);

      if(  hypt[kn]->GetEntries() < 0 ) continue;
      
      Double_t rapm  = hypt[kn]->GetMean(1);
      Double_t rape  = hypt[kn]->GetStdDev(1);


      if( hyphi1[kn]->GetEntries() < 0 ) continue;
      auto hyphi = (TH1D*)hyphi1[kn]->ProjectionX("",0,-1,"eo");

      Double_t corr[2]={mcos1[isys[m]], mcos1e[isys[m]]};
      GetFittingParameters(*hyphi, para, corr);
      
      if( !std::isnan(para[4]) && !std::isinf(para[4] )) {
	gv_v1->SetPoint( kl,     rapm, para[4]);
	gv_v1->SetPointError( kl, rape, para[5]);
	kl++;
      }

      // pt dependence 
      UInt_t il = 0; 
      for(UInt_t jn = 0; jn < nbin1; jn++){

	id1++; cc[ic]->cd(id1); 

	if( hyptphi1[kn][jn]->GetEntries() > 0 ) {	

	  Double_t ptc  = hyptphi1[kn][jn]->GetMean(2);
	  Double_t ptce = hyptphi1[kn][jn]->GetStdDev(2);
	    
	  auto hypt = (TH1D*)hyptphi1[kn][jn]->ProjectionX();
	  GetFittingParameters(*hypt, para, corr);
	    
	  if( !std::isnan(para[4]) && !std::isinf(para[4] )) {
	    gPt_v1[kn]->SetPoint(il, ptc, para[4]);
	    gPt_v1[kn]->SetPointError(il, ptce, para[5]);
	    
	    il++;
	  }
	}
      }
      gPt_v1[kn]->Write();
    }

    kl = 0;
    for(UInt_t kn = 0; kn < ybin2; kn++){

      // v2 vs rapidity
      cc[ic+2]->cd(2*(kn+1));

      if(  hypt2[kn]->GetEntries() < 0 ) continue;

      Double_t rapm  = hypt2[kn]->GetMean(1);
      Double_t rape  = hypt2[kn]->GetStdDev(1);

      if( hyphi2[kn]->GetEntries() < 0 ) continue;
      auto hyphi = (TH1D*)hyphi2[kn]->ProjectionX("",0,-1,"eo");
      
      Double_t corr[2] = {mcos2[isys[m]], mcos2e[isys[m]]};
      GetFittingParameters(*hyphi, para, corr);


      if( !std::isnan(para[4]) && !std::isinf(para[4] )) {
	  gv_v2->SetPoint( kl,    rapm, para[4]);
	  gv_v2->SetPointError( kl, rape, para[5]);
	  kl++;
      }
    

      UInt_t il = 0; 
      for(UInt_t jn = 0; jn < nbin2; jn++){
	
	id2++; cc[ic+1]->cd(id2); 
	if( hyptphi2[kn][jn]->GetEntries() > 0 ){

	  Double_t ptc  = hyptphi2[kn][jn]->GetMean(2);
	  Double_t ptce = hyptphi2[kn][jn]->GetStdDev(2);

	  auto hypt = (TH1D*)hyptphi2[kn][jn]->ProjectionX();

	  Double_t corr[2] = {mcos2[isys[m]], mcos2e[isys[m]]};
	  GetFittingParameters(*hypt, para, corr);

	  if( !std::isnan(para[4]) && !std::isinf(para[4] )) {
	    gPt_v2[kn]->SetPoint(il, ptc, para[4]);
	    gPt_v2[kn]->SetPointError(il, ptce, para[5]);
	    
	    il++;
	  }
	}
      }
      gPt_v2[kn]->Write();
    }


    ic+=3;


    //plotting
    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
    cc[ic]->Divide(ybin1,1);
    for(UInt_t kn = 0; kn < ybin1; kn++){
      cc[ic]->cd(id); id++;
      gPt_v1[kn]->Draw("ALP");
    }

    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1000,500);
    cc[ic]->Divide(ybin2,1);
    for(UInt_t kn = 0; kn < ybin2; kn++){
      cc[ic]->cd(id); id++;
      gPt_v2[kn]->Draw("ALP");
    }

    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,400);
    cc[ic]->Divide(ybin2,1);
    for(UInt_t kn = 0; kn < ybin2; kn++){
      cc[ic]->cd(id); id++;
      hypt2[kn]->Draw("colz");
    }

    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,400);
    cc[ic]->Divide(ybin1,1);
    for(UInt_t kn = 0; kn < ybin1; kn++){
      cc[ic]->cd(id); id++;
      hypt[kn]->Draw("colz");
    }


    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
    gv_v1->Draw("ALP");

    auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
    auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
    LineV->SetLineStyle(3);
    LineH->SetLineStyle(3);
    LineV->Draw();
    LineH->Draw();


    ic++; id = 1;
    cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
    gv_v2->Draw("ALP");


    gv_v1->Write();
    gv_v2->Write();

    GraphSave->WriteObject(&arr_range1, "v1_yrange");
    GraphSave->WriteObject(&arr_range2, "v2_yrange");
  }

  
  gSystem->cd("..");

}

void combineHist(UInt_t sys, TH2D *hyptphi1, TH2D *hyptphi2)
{

}

void GetFittingParameters(TH1D &h1, Double_t pp[6])
{
  for(UInt_t i = 0; i < 6; i++)
    pp[i] = 0.;


  Double_t nbin = h1.GetNbinsX();
  Double_t scf  = h1.GetEntries();
  h1.Scale(nbin/scf);
  // h1.SetMaximum(1.2);
  // h1.SetMinimum(0.8);

  h1.Fit("fcos1","Q","",-3.1,3.1);

  pp[0] = fcos1->GetParameter(0);
  pp[1] = fcos1->GetParameter(1);
  pp[2] = fcos1->GetParError(0);
  pp[3] = fcos1->GetParError(1);
  

}
void GetFittingParameters(TH1D &h1, Double_t pp[6], Double_t corr[2])
{
  GetFittingParameters(h1,pp);
  
  if( corr[0] != 0 && corr[1] != 0 ) {
    pp[4] = pp[1]/corr[0];
    pp[5] = sqrt( pow(pp[4],2)*( pow(pp[3]/pp[1],2) + pow(corr[1]/corr[0], 2) ));
  }
  else {
    pp[4] = 0.;
    pp[5] = 0.;
  }
}


//--------------------------------------------//%% Executable : 
void Template()                   
{
  //----- Parametres

  //----- Booking
  for(Int_t m = m_bgn; m < m_end; m++){

  }

  //----- Filling
  for(Int_t m = m_bgn; m < m_end; m++){
    Int_t nEntry = rChain[m]->GetEntries();

    for(Int_t i = 0; i < nEntry; i++){
      rChain[m]->GetEntry(i);
    }
  }
  //----- Drawing 
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,1000);
  cc[ic]->Divide(3,m_end);

  UInt_t id = 1;
  for(Int_t m = m_bgn; m < m_end; m++){
  }
}


UInt_t SetBranch(UInt_t m)
{

  if(aArray != NULL)
    aArray->Clear();


  if(rChain[m] == NULL) {
    std::cout << " no file is loaded " << std::endl;
    return 0;
  }

  rChain[m]->SetBranchAddress("STParticle",&aArray);
  rChain[m]->SetBranchAddress("ntrack",   ntrack);
  rChain[m]->SetBranchAddress("aoq",&aoq);
  rChain[m]->SetBranchAddress("z",&z);
  rChain[m]->SetBranchAddress("snbm",&snbm);
  rChain[m]->SetBranchAddress("ProjA",&ProjA);
  rChain[m]->SetBranchAddress("ProjB",&ProjB);
  rChain[m]->SetBranchAddress("unitP_fc"  ,&unitP_fc,&bunitP_fc);
  rChain[m]->SetBranchAddress("unitP_rc"  ,&unitP_rc,&bunitP_rc);
  rChain[m]->SetBranchAddress("unitP_1"   ,&unitP_1,&bunitP_1);
  rChain[m]->SetBranchAddress("unitP_2"   ,&unitP_2,&bunitP_2);
  rChain[m]->SetBranchAddress("mtrack"    ,&mtrack);
  rChain[m]->SetBranchAddress("mtrack_1"  ,&mtrack_1);    rChain[m]->SetBranchAddress("mtrack_2"  ,&mtrack_2);
  rChain[m]->SetBranchAddress("unitP_lang",&unitP_lang,&bunitP_lang);
  rChain[m]->SetBranchAddress("STNeuLANDCluster", &aNLClusterArray);

  return rChain[m]->GetEntries();

}
