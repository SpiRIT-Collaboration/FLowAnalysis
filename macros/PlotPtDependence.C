#include "DoFlow.h"
#include "SetStyle.C"
#include "FlowFunction.C"
#include "SetColor.C"
#include "calculateImpactParameter.C"
#include "PlotTransport.C"

FairLogger *logger = FairLogger::GetLogger();

struct gplot{
  TString Version;
  TString fileHeader;
  TString comment;
};


TString  bName[]   = {"132Sn_","108Sn_","124Sn_","112Sn_","100Sn_"};
Double_t sysdlt[]  = {0.22,    0.09,      0.15,   0.15   , 0.22};
Double_t sysA[]    = {256.,    220.,      236.,   236.   ,  256.};

//==================================================
//-- plot configuration
//--------------------------------------------------
Bool_t bplot[] = 
  { 1, // 0 data   It should be set to 1 in the code.
    0, // 1 model  It should be set to 1 in the code.
    1, // 2 v1 and v2 rapidity 
    0, // 3 v1 and v2 on pt in one window
    0, // 4 v1 and v2 in individual windows
    0, // 5 Acceptance ypt => never be plotted 
    0, // 6 <px>/A
    0, // 7 v1 vs part system dependence
    1, // 8 v2_min system dependence 
    0, // 9 v1 slope(V11) and v2 max dependence on m
    0, //10 v1 slope(v11) and v2 max dependence on impact parameter
    0, //11 v1 slop and v2 max systematic error check
    0, //12 v11 and v20 plot
  };

Bool_t bstyle[] =
  { 0,  // 0 132Sn p,d,t, v1v2-y (for NUSYM2019)
    0,  // 1 Model comparison 
    0,  //   RPSim fv1y and fv2y is plotted.
  };
//==================================================


// --> Plotting selection
//--- Data
Bool_t bsys[]  = { 1, 1, 0, 1, 0};    //132Sn, 108Sn, 124Sn, 112Sn, Sim
Bool_t bpid[]  = { 1, 1, 1, 1, 0, 0, 0}; //0:p, 1:d, 2:t, 3:3He, 4:4He, 5:n 6:H
Bool_t bcnt[]  = { 1, 0, 0}; 
Bool_t bmdl[]  = { 0,      0,         0,       0,      0,       0,       0      }; 
TString nmTrs[]= {"AMD",   "ChiBUU",  "IBUU",  "IQMD", "SMF",   "UrQMD", "dcQMD"};
UInt_t cntw = 1;
UInt_t iv2at = 4;
//-----------

gplot gnames[] = {  
  //  {".v52.14.10" ,"finYPt_","ut>0 w/ocrr.b3fm"},
  {".v52.14.12" ,"finYPt_","ut>0 wcrr.b3fm"},
  //  {".v52.11.9" ,"finYPt_","ut>0.crr.b3fm"},
  //  {".v52.11.9" ,"finYPt_","w/ocorr.b3fm"},
  // {".v52.11.3" ,"finYPt_","ut>0.w/ocrr.b3fm"},
  // {".v52.11.2" ,"finYPt_","ut>0.crr.b3fm"},
  // {".v52.11.1" ,"finYPt_","crr.b3fm"},
  //  {".v52.10.28" ,"finYPt_","b3fm"},
  //  {".v52.10.26" ,"finYPt_",""},

  //  {".v52.10.27" ,"finYPt_",""},
  //  {".v22.0.0" ,"finYPt_","SIM full"},
  //  {".v50.8.0" ,"finYPt_","SIM TPC"},
  //  {".v52.10.30" ,"finYPt_","M65-80"},
  // {".v52.10.31" ,"finYPt_","M60-65"},
  // {".v52.10.32" ,"finYPt_","M55-60"},
  // {".v52.10.33" ,"finYPt_","M50-55"},
  // {".v52.10.34" ,"finYPt_","M45-50"},
  // {".v52.10.35" ,"finYPt_","M40-45"},
  // {".v52.10.36" ,"finYPt_","M30-40"},
  // {".v52.10.37" ,"finYPt_","M0-30"},
  // {".v52.0.13" ,"advYPt_","y_nn"},
};
//const UInt_t nmax = (UInt_t)sizeof(bver)/sizeof(UInt_t);
const UInt_t nmax = (UInt_t)sizeof(gnames)/sizeof(gplot);

//APR2019
//  {".v50.7.1","advYPt_","(A)TPC"}, //APR2019
//  {".v22.0.7","advYPt_","(B)Full"},//APR2019
//  {".v50.0"  ,"advYPt_","(A)TPC"},
//  {".v22.0.8","advYPt_","(B)Full"},

TString *sVer  = new TString[nmax];
TString *sName = new TString[nmax];
TString *cmnt  = new TString[nmax];


//==================================================
UInt_t   ccvid = 0;

// --> configuration
Size_t  imsz[] = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
//Color_t icolnn[]={ kPink, kBlue+1, kDeepSea+2, kViolet-1};


TMultiGraph  *mv1[ybin1];
TMultiGraph  *mv2[ybin2];
TLegend      *lg1[ybin1];
TLegend      *lg2[ybin2];

TGraphErrors *g_v11;
TGraphErrors *g_v20;
TGraphErrors *g_v2n;
std::vector< TString > lx_v11;
std::vector< TString > lx_v20;
std::vector< Double_t > d_v11;
std::vector< Double_t > d_v11e;
std::vector< Double_t > d_v20;
std::vector< Double_t > d_v20e;
std::vector< Double_t > d_v2n;
std::vector< Double_t > d_v2ne;


TString  hmultX;

void RapidityShift(TGraphErrors *gv);
void ShiftX(TGraphErrors *gv, Double_t off);
void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer);
void GetAveragev2(TGraphErrors *gr);
void FindFirstBin( TH1I *h, Int_t &min, Int_t &max );
void DividebyTwo(TGraphErrors *yv );


Double_t FittingAndIntegral(TGraphErrors *gr)
{
  TF1 *p1 = new TF1("p1","[0]+[1]*x",0.,800.);
  
  gr->Fit("p1","","",200.,700.);
  
  return p1->Integral(0., 800) / 800.;
}



void PlotPtDependence()
{
  
  gStyle->SetOptStat(0);
  SetStyle();
  SetColor();

  if( bplot[10] )
    calculateImpactParameter();


  UInt_t ndata = 0;
  for(UInt_t i = 0; i < nmax; i++){
    ndata++;
    sVer[i]  = gnames[i].Version;
    sName[i] = gnames[i].fileHeader;
    cmnt[i]  = gnames[i].comment;

    if( gnames[i].Version == "" ) {
      std::cout << " Version? v##.# ";
      TString svv; 
      std::cin  >> svv;
      sVer[i] = ".v" + svv;
    }
  }

  TString ltitle;


  UInt_t ngr = 0;
  std::vector< TString > fname;
  std::vector< std::vector< UInt_t >> index(6);


  TGraphErrors *g_mpxslp = new TGraphErrors();
  g_mpxslp->SetName("g_mpxslp");
  g_mpxslp->SetTitle("; (N-P)/A; F[MeV/c]");

  TGraphErrors *g_v1slp = new TGraphErrors();
  g_v1slp->SetName("g_v1slp");
  g_v1slp->SetTitle("; multiplicity; v_{11}");


  TGraphErrors *g_v2max = new TGraphErrors();
  g_v2max->SetName("g_v2max");
  g_v2max->SetTitle("; multiplicity; -v2_max");

  TGraphErrors *g_v2sysA = new TGraphErrors();
  g_v2sysA->SetName("g_v2sysA");
  g_v2sysA->SetTitle(";A_{beam} + A_{target}; v2");

  TGraphErrors *g_v2sysD = new TGraphErrors();
  g_v2sysD->SetName("g_v2sysD");
  g_v2sysD->SetTitle(";(N-P)/A; v2");

  TGraphErrors *g_v1slpsys[20];
  TGraphErrors *g_v2maxsys[20];
  TGraphErrors *g_v1slpm[20];
  TGraphErrors *g_v2maxm[20];

  for(UInt_t jj = 0; jj < 20; jj++){
    g_v1slpsys[jj] = new TGraphErrors();
    g_v2maxsys[jj] = new TGraphErrors();    
    g_v1slpm[jj]   = new TGraphErrors();
    g_v2maxm[jj]   = new TGraphErrors();
    g_v1slpm[jj] -> SetName(Form("g_v1slpm%d",jj));
    g_v1slpm[jj] -> SetTitle("; multiplicity; v_{11}");
    g_v2maxm[jj] -> SetName(Form("g_v2maxm%d",jj));
    g_v2maxm[jj] -> SetTitle("; multiplicity; -v2 max");
  }

  // Rapidity dependence 
  TH2D *hyptacp[20];

  auto mrv1  = new TMultiGraph("mrv1"  ,";y/y_{nn}-1; v1");
  auto mrv2  = new TMultiGraph("mrv2"  ,";y/y_{nn}-1;v2");//";y_{cm}/y_{beam}; v2");
  auto mv1sl = new TMultiGraph("mv1sl" ,";Centrality; v1 slope");
  auto mv1slp= new TMultiGraph("mv1slp","; Particle ; v1 slope");
  auto mmpx  = new TMultiGraph("mmpx","; y/y_{nn}-1 ; <px>/A");
  //  mv1slp->GetXaxis()->SetAlphanumeric(kTRUE);

  auto lgr1 = new TLegend(0.3 , 0.65, 0.60, 0.90,""); 
  auto lgr2 = new TLegend(0.5 , 0.60, 0.70, 0.9, "");
  //  auto lgr2 = new TLegend(0.5 , 0.20, 0.70, 0.45, "");
  auto lgr3 = new TLegend(0.35, 0.13, 0.7, 0.33, "");
  auto lgr4 = new TLegend(0.15, 0.63, 0.5, 0.85, "");
  auto lgr5 = new TLegend(0.15, 0.63, 0.5, 0.85, "");
  auto lgr6 = new TLegend(0.72, 0.72, 0.94, 0.88,"");
  auto lgr7 = new TLegend(0.16, 0.70, 0.46, 0.85,"");
  auto lgr8 = new TLegend(0.16, 0.70, 0.46, 0.85,"");
  auto lgr9 = new TLegend(0.16, 0.70, 0.46, 0.85,"");
  auto lgr0 = new TLegend(0.16, 0.70, 0.46, 0.85,"");
  auto lgr10= new TLegend(0.16, 0.70, 0.46, 0.85,"");

  g_v11 = new TGraphErrors();
  g_v11 -> SetTitle("v_{11}");
  g_v20 = new TGraphErrors();
  g_v20 -> SetTitle("v_{20}");
  g_v2n = new TGraphErrors();
  g_v2n -> SetTitle("v_{2n}");

  
  TString fOutName = "";

  UInt_t inx = 0;
  UInt_t nextp = 0;
  for(UInt_t is = 0; is < 5; is++){

    for(UInt_t ip = 0; ip < (UInt_t)sizeof(bpid)/sizeof(Bool_t); ip++){

      if( !bsys[is] || !bpid[ip] ) continue;

      for(UInt_t it = 0; it < (UInt_t)sizeof(bcnt)/sizeof(Bool_t); it++){
	
	if( !bcnt[it] ) continue;

	for(UInt_t iz = 0; iz < nmax; iz++){
	  
	  fname.push_back( "data/"+ sName[iz] + bName[is] + fpid[ip] + sVer[iz] + ".root" );

	  if( fOutName == "" )
	    fOutName = "data/" + bName[is] + fpid[ip] + sVer[iz] + "_plt.root";

	  std::cout << fname.at(ngr) << std::endl;
	  ltitle = "";
	  //	  ltitle = Form("mult %d ~ %d",cent[it],cent[it+cntw]);
	  //	  std::cout << " label title " << ltitle << std::endl;

	  index[0].push_back(is); // system
	  index[1].push_back(ip); // pid
	  index[2].push_back(it); // centrality
	  index[3].push_back(it+cntw);
	  index[4].push_back(iz); // version

	  if( is == 4 ) {
	    index[5].push_back(4);
	    FlowFunction();
	  }
	  else
	    index[5].push_back(1);

	  if(ngr == 0)
	    nextp = 10*is + ip;

	  ngr++;

	  bplot[0] = 1;
	}
      }


      for(UInt_t imdl = 0; imdl < (UInt_t)sizeof(bmdl)/sizeof(Bool_t); imdl++){
	if( !bmdl[imdl] ) continue;
	
	for( UInt_t iypt = 0; iypt < 2; iypt++ ) {
	  if( iypt == 0 && !bplot[2] ) continue; 
	  if( iypt == 1 && ( !bplot[3] && !bplot[4] )) continue;

	  auto foundfile = FindTransportData(imdl, iypt, ip);
	  if( foundfile != "" ) {
	    std::cout << foundfile << std::endl;
	    fname.push_back( foundfile );

	    index[0].push_back(is); // system
	    index[1].push_back(ip); // pid
	    index[2].push_back(0); // centrality
	    index[3].push_back(0);
	    index[4].push_back(0); // version
	    index[5].push_back(5+imdl);  // transport
	    
	  }
	}
      }
    }
  }




  cout << " Number of data files --> " << fname.size() << endl;


  TFile *fOutput = new TFile(fOutName,"recreate");
  TFile *fOpen;


  // Pt dependence setup
  for(UInt_t k = 0; k < ybin1; k++){
    mv1[k]  = new TMultiGraph((TString)Form("mv1%d",k),  ";Pt [MeV/c]; v_1");
    lg1[k] = new TLegend(0.28, 0.18, 0.65, 0.31); 
    lg1[k]->SetTextSize(0.05);
  }

  for(UInt_t k = 0; k < ybin2; k++) {
    mv2[k]  = new TMultiGraph((TString)Form("mv2%d",k),  ";Pt [MeV/c]; v_2");
    lg2[k] = new TLegend(0.28, 0.18, 0.65, 0.31); 
    lg2[k]->SetTextSize(0.05);
  }

  UInt_t islp = 0;
  UInt_t isp = 0;
  UInt_t isk = 0;
  UInt_t iss = 9;
  UInt_t ipp = 0;
  UInt_t ips1[4] = {0,0,0,0};
  UInt_t ips2[4] = {0,0,0,0};
  Color_t icolor = 2; 
  UInt_t iv1f = 0;
  UInt_t isp1[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  UInt_t isp2[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  for(UInt_t igr = 0; igr < (UInt_t)fname.size(); igr++ ) {

    UInt_t is = index[0].at(igr);
    UInt_t ip = index[1].at(igr);
    UInt_t it = index[2].at(igr);
    UInt_t ik = index[3].at(igr);
    UInt_t iz = index[4].at(igr);
    UInt_t ia = index[5].at(igr);
    
    if( is != iss ) {
      ipp = 0;
      iss = is;
    }

    TString ohtitle = lsys[is];//+" "+lpid[ip]+" ";
    TString otitle  = ohtitle;

    // ia: -------------
    // 1 : data
    //-------------------
    if( ia == 1 ) { //data
      fOpen = TFile::Open(fname.at(igr)); 

      if( fOpen == NULL ) continue;
      else
	std::cout << fname.at(igr) << " is opened. " << std::endl;

      std::cout << "mnt " << cmnt[iz] << " " << iz << endl;
      otitle += "("+ lpid[ip]+")"+cmnt[iz];
      //    }

    //acceptance
      if( bplot[0] && bplot[5] && kFALSE ) {
	hyptacp[igr] = (TH2D*)fOpen->Get("hyptacp");
	if( hyptacp[igr] == NULL) {
	  cout <<" NULL 222" << endl;
	} 
	else {
	  hyptacp[igr]->Print();
	  hyptacp[igr]->SetName(Form("hyptacp_%d",igr));
	  cout << " igr " << igr << " " <<hyptacp[igr]->GetName() << endl;
	  cout << hyptacp[igr] << endl;
	  hyptacp[igr]->Draw("colz");
	}
	//    }
	
	//    if( igr < ngr ) {
	TH1D *fevt = (TH1D*)fOpen->Get("hphi0_180");
	
	if( fevt != NULL ) {
	  UInt_t tevt = fevt->GetEntries();
	  
	  auto ff = (TH2D*)fOpen->Get("hyptacp");
	  ff->SetTitle(otitle);
	  ff->SetDirectory(gROOT);
	  
	  gROOT->cd();
	  ff->ProjectionX("_px");
	  TH1D* ff1 = (TH1D*)gROOT->Get((TString)ff->GetName()+"_px");
	  
	  Double_t xwd = (Double_t)ff1->GetXaxis()->GetNbins()/(ff1->GetXaxis()->GetXmax() - ff1->GetXaxis()->GetXmin());
	  
	  if( fpid[ip] == "neutron" )
	    ff1->Scale( xwd / tevt /0.26);
	  else
	    ff1->Scale(xwd / tevt);
	  
	  ff1->SetTitle("; Y_{cm}; dN/dy");
	  ff1->SetLineColor(icolor);
	  lgr3->AddEntry(ff1, otitle );
	}
      }
      fOpen->cd();
    }
    //------------------------------    


    // rapidity dependence
    //--- v1 vs rapidity ---
    if( ( bplot[0] || bplot[1]) && bplot[2] ) {

      TGraphErrors *yv1 = (TGraphErrors*)fOpen->Get("gu_v1");
      //      TGraphErrors *yv1 = (TGraphErrors*)fOpen->Get("gy_v1");


      if( ia >= 5 ) { // transport model
	yv1 = ReadTransportData( fname.at(igr) , 1);

	otitle += nmTrs[ ia - 5 ];
      }

      if( yv1 != NULL ) {

	LOG(INFO) << " yv1 : " << yv1->GetName() << FairLogger::endl;
	//	yv1->Print();


	auto yv1_rev = new TGraphErrors();
	for( UInt_t i = 0; i < yv1->GetN(); i++ ) {
	  Double_t x = 0, y = 0;
	  yv1->GetPoint( yv1->GetN()-1-i, x, y );
	  yv1_rev->SetPoint( i, -x, -y );
	}


	yv1->SetMarkerColor(icolor);
	yv1->SetMarkerStyle(imark[is]);
	yv1_rev->SetMarkerColor(icolor);
	yv1_rev->SetMarkerStyle(imark[is]);

	// --- check it
	if( kFALSE ) { //|| is != 5 ) {
	  yv1->RemovePoint(10);
	  yv1->RemovePoint(9);
	  yv1->RemovePoint(8);
	  yv1->RemovePoint(7);
	  yv1->RemovePoint(0);
	  yv1->RemovePoint(0);
	}

	if( ip == 5 )
	  yv1->SetMarkerStyle(imark[is]+4);

	yv1->SetMarkerSize(imsz[is]);
	yv1->SetLineColor(icolor);
	yv1_rev->SetMarkerStyle(imark[is]+4);
	yv1_rev->SetMarkerSize(imsz[is]);
	yv1_rev->SetLineColor(icolor);

	fv1fit->SetLineColor(icolor);

	Double_t v1para1[4]={0.4, 0.4, 0.5, 0.6};

	fv1fit->SetParameter(0,0.);
	fv1fit->SetParameter(1,v1para1[ip]);
	fv1fit->SetParameter(2,-0.2);
	fv1fit->SetParameter(3,-0.01);
	  
	if( kFALSE ) {//is == 2 ) {
	  mrv1->Add(yv1_rev,"p");
	  lgr1->AddEntry(yv1_rev,  otitle+"_rev" ,"lp");
	  yv1_rev->Fit("fv1fit","Q0","",-1.1,1.1); //"Q0","");
	}
	else {
	  mrv1->Add(yv1,"p");
	  lgr1->AddEntry(yv1,  otitle ,"lp");
	  yv1->Fit("fv1fit","","",-0.6, 0.8); //"Q0","");
	}

	if( bplot[12] ) {
	  d_v11.push_back(fv1fit->GetParameter(1));
	  d_v11e.push_back(fv1fit->GetParError(1));
	  
	  g_v11->SetPoint((UInt_t)lx_v11.size(), (UInt_t)lx_v11.size()+0.5, fv1fit->GetParameter(1) );
	  g_v11->SetPointError((UInt_t)lx_v11.size(),                 0., fv1fit->GetParError(1) );
	  lx_v11.push_back(otitle);
	}
      }

      //--- y vs v2 ---
      TGraphErrors *yv2 = (TGraphErrors*)fOpen->Get("gu_v2");
      //       TGraphErrors *yv2 = (TGraphErrors*)fOpen->Get("gy_v2");

      if( ia >= 5 ) { // transport model
	yv2 = ReadTransportData( fname.at(igr) , 2);
      }

      
      if( yv2 != NULL ) {
	for( Int_t iip = (Int_t)yv2->GetN()-1; iip >= 0; iip-- ){
	  Double_t xpnt, ypnt;
	  //	  yv2->GetPoint(iip, xpnt, ypnt);
	  //	  if( (xpnt > 1 || xpnt < -0.8) && is != 5 )
	  //	    yv2->RemovePoint(iip);
	  ypnt = yv2->GetErrorY(iip);
	  if( ypnt > 0.05 )
	    yv2->RemovePoint(iip);
	}
	  
	
	yv2->SetMarkerColor(icolor);
	yv2->SetMarkerStyle(imark[is]);
	if( ip == 5 || ia != 1)
	  yv2->SetMarkerStyle(imark[is]+4);

	yv2->SetMarkerSize(imsz[is]);
	yv2->SetLineColor(icolor);
	fv2fit->SetLineColor(icolor);
	

	Double_t v2para0[4]={-0.03, -0.04, -0.05, -0.05};
	if( ip == 3 )
	  yv2->Fit("fv2fit","","",-0.48,1.);
	else {
	  fv2fit->SetParameter(0,v2para0[ip]);
	  fv2fit->SetParameter(1,0.05);
	  fv2fit->SetParameter(2,0.02);
	  yv2->Fit("fv2fit","","",-0.5,0.5);
	}

	d_v20.push_back(fv2fit->GetParameter(0));
	d_v20e.push_back(fv2fit->GetParError(0));
	
	g_v20->SetPoint((UInt_t)lx_v20.size(), (UInt_t)lx_v20.size()+0.5, fv2fit->GetParameter(0) );
	g_v20->SetPointError((UInt_t)lx_v20.size(),                   0., fv2fit->GetParError(0) );
	lx_v20.push_back(otitle);

	mrv2->Add(yv2,"p");
	lgr2->AddEntry(yv2,  otitle ,"lp");
	
      }


      if( (bplot[10] || bplot[7] || bplot[8]) &&  ia == 1 ) {

	TH1I *hmult  = (TH1I*)fOpen->Get("hmult");
	if( hmult == NULL ) 
	  continue;
	
	// hmultX = " Impact Parameter [fm] ";
	// Int_t mmin = hmult->GetXaxis()->GetXmin();;
	// Int_t mmax = hmult->GetXaxis()->GetXmax();
	// FindFirstBin(hmult, mmin, mmax);
	  
	// mstd  = GetMinB( Asys[is], mmean-mstd, mmean+mstd );
	// mmean = GetMinB( Asys[is], mmin, mmax);
	// cout << "getminB  mean " << mmean << " +-" << mstd << endl;
	
	hmultX = "Multiplicity";
	Double_t mmean = hmult->GetMean();
	//Double_t mstd  = hmult->GetMeanError();
	Double_t mstd  = hmult->GetStdDev()/sqrt(hmult->GetEntries());
	cout << "get mean " << mmean << " +- " << mstd<< endl;

	//v1
	Double_t constlslope = fv1fit->GetParameter(0);
	Double_t slope = fv1fit->GetParameter(1);
	Double_t slopee= fv1fit->GetParError(1);
	
	std::cout << "************" << std::endl;
	std::cout << otitle << " " << slope << " +- " << slopee << std::endl;      

	g_v1slp -> SetPoint(iv1f, mmean, slope);
	g_v1slp -> SetPointError(iv1f, mstd, slopee);	
	iv1f++;

	if( nextp != 10*is + ip ) {
	  inx++;
	  nextp = 10*is + ip;
	}
	cout << " nextp " << nextp << " " << inx << endl;

	g_v1slpm[inx]->SetTitle(fsys[is]);
	g_v1slpm[inx]->SetPoint     (isp1[inx], mmean, slope);
	g_v1slpm[inx]->SetPointError(isp1[inx], mstd, slopee);
	isp1[inx]++;


	g_v1slpsys[ip]->SetPoint     (ips1[ip], sysdlt[is], slope);
	g_v1slpsys[ip]->SetPointError(ips1[ip], 0., slopee);
	ips1[ip]++;
	
	//v2
	Double_t v2x, v2y, v2ye;
	//	GetMinimumv2(yv2, v2y, v2ye);

	v2y  = fv2fit->GetParameter(0);
	v2ye = fv2fit->GetParError(0);

	cout << " v2y " << v2y << " +- " << v2ye << endl;

	g_v2max->SetPoint(isp , mmean, -v2y);
	g_v2max->SetPointError(isp, mstd, v2ye);
	isp++;
	

	g_v2maxm[inx]->SetTitle(fsys[is]);
	g_v2maxm[inx]->SetPoint     (isp2[inx] , mmean, -v2y);
	g_v2maxm[inx]->SetPointError(isp2[inx] , mstd, v2ye);
	isp2[inx]++;
	

	g_v2maxsys[ip]->SetPoint(ips2[ip], sysdlt[is], -v2y);
	g_v2maxsys[ip]->SetPointError(ips2[ip], 0., v2ye);
	ips2[ip]++;
	ipp++;
	

	if( bsys[ik] ) {
	  g_v2sysA->SetPoint(isk, sysA[is], v2y);
	  g_v2sysA->SetPointError(isk, 0., v2ye);
	  
	  g_v2sysD->SetPoint(isk, sysdlt[is], v2y);
	  g_v2sysD->SetPointError(isk,  0.,  v2ye);
	  isk++;
	}
      }
    }
      

    //--- <px>/A --- 
    if( ( bplot[0] || bplot[1]) && bplot[6] ) {

      TGraphErrors *ympx= (TGraphErrors*)fOpen->Get("gmpx");    
      if( ympx != NULL) {

	ympx->SetMarkerColor(icolor);
	ympx->SetMarkerStyle(imark[is]);
	
	ympx->SetMarkerSize(imsz[is]);
	ympx->SetLineColor(icolor);
	  
	mmpx->Add(ympx,"lp");
	lgr9->AddEntry(ympx, otitle, "lp");
	
	ympx->Fit("lslope","Q0","",-0.4,0.4);
	Double_t constlslope = lslope->GetParameter(0);
	Double_t slope       = lslope->GetParameter(1)*1.4;
	Double_t slopee      = lslope->GetParError(1);
	std::cout << " <px> slop " << slope << " +- " << slopee << std::endl;
	
	//slope
	g_mpxslp-> SetPoint(islp, sysdlt[is], slope);
	g_mpxslp-> SetPointError(islp, 0, slopee);
	islp++;
      }
    }


    //--- Pt Ut dependence 
    if( ( bplot[0] || bplot[1]) && (bplot[3] || bplot[4]) ) {

      Double_t v1mx = -99.;
      Double_t v1mn =  99.;
      Double_t v2mx = -99.;
      Double_t v2mn =  99.;

      for(UInt_t k = 0; k < ybin1; k++){

	TGraphErrors *gr_v1 = (TGraphErrors*)fOpen->Get((TString)Form("gUt_v1%d",k));
	if( gr_v1 == NULL )
	  gr_v1 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v1%d",k));
	
	if( gr_v1 == NULL ) continue;
	if( igr == 0 ) {
	  mv1[k]->SetTitle(gr_v1->GetTitle());
	  mv1[k]->GetXaxis()->SetTitle(gr_v1->GetXaxis()->GetTitle());
	  mv1[k]->GetYaxis()->SetTitle("v1");
	}

	TGraphErrors *gr_v1A = new TGraphErrors();
	for(UInt_t j = 0; j < (UInt_t)gr_v1->GetN(); j++ ){
	  Double_t x, y;
	  gr_v1->GetPoint(j, x, y);
	  //x /= partA[ip+2];
	  //	y /= partA[ip+2];
	  auto xe = gr_v1->GetErrorX(j);
	  auto ye = gr_v1->GetErrorY(j);

	  if( x < 550. ) {
	    gr_v1A->SetPoint(j, x, y);	
	    gr_v1A->SetPointError(j,xe,ye);
	  }
	}

      
	if( gr_v1A == NULL ) continue;
	gr_v1A->SetMarkerStyle(imark[is]);
	

	if( ip == 5 )
	  gr_v1A->SetMarkerStyle(imark[is]+4);
	gr_v1A->SetMarkerColor(icolor);
	gr_v1A->SetMarkerSize(imsz[is]);
	gr_v1A->SetLineColor(icolor);
	gr_v1A->GetXaxis()->SetRangeUser(0.,550.);

	mv1[k]->Add(gr_v1A,"lp");
	lg1[k]->AddEntry(gr_v1A,  otitle,"lp");
	
      }
    
      cout << " ==== " << endl;

      for(UInt_t k = 0; k < ybin2; k++){
	
	TGraphErrors *gr_v2 = (TGraphErrors*)fOpen->Get((TString)Form("gUt_v2%d",k));
	
	if( gr_v2 == NULL )
	  gr_v2 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v2%d",k));
	
	if( gr_v2 == NULL ) continue;
	
	if( igr == 0 ) {
	  mv2[k]->SetTitle(gr_v2->GetTitle());
	  mv2[k]->GetXaxis()->SetTitle(gr_v2->GetXaxis()->GetTitle());
	  mv2[k]->GetYaxis()->SetTitle("v2");
	}

	TGraphErrors *gr_v2A = new TGraphErrors();
	for(UInt_t j = 0; j < (UInt_t)gr_v2->GetN(); j++ ){
	  Double_t x, y;
	  gr_v2->GetPoint(j, x, y);
	  //x /= partA[ip+2];
	  //y /= partA[ip+2];
	  auto xe = gr_v2->GetErrorX(j);
	  auto ye = gr_v2->GetErrorY(j);

	  if( x < 550. ) {
	    gr_v2A->SetPoint(j, x, y);	
	    gr_v2A->SetPointError(j,xe,ye);
	  }
	}
	
	//      RemovePtPoint(ip,gr_v2);
	if( v2mx > gr_v2->GetYaxis()->GetXmax() )
	  v2mx = gr_v2->GetYaxis()->GetXmax();
	
	if( v2mn < gr_v2->GetYaxis()->GetXmin() )
	  v2mx = gr_v2->GetYaxis()->GetXmin();
	
	gr_v2A->SetMarkerStyle(imark[is]);
	if( ip == 5 )
	  gr_v2->SetMarkerStyle(imark[is]+4);
	gr_v2A->SetMarkerColor(icolor);
	gr_v2A->SetMarkerSize(imsz[is]);
	gr_v2A->SetLineColor(icolor);
	gr_v2A->GetXaxis()->SetRangeUser(0.,550.);

	mv2[k]->Add(gr_v2A,"lp");
	lg2[k]->AddEntry(gr_v2A, otitle ,"lp");

	//      mv2[k]->GetYaxis()->SetRangeUser(v2mn, v2mx);
      }
    }	

    icolor++;
    if( icolor == 10 )
      icolor++;

    fOpen->Close();
  }

  //-------------->>>> Plotting >>>>>---------------------

  //--- v1 and v2 vs rapidity ---
  if( ( bplot[0] || bplot[1]) && bplot[2] ) {

    Double_t Ymin, Ymax, Xmin, Xmax;
    //--- v1 vs rapidity ---
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

    mrv1->GetXaxis()->SetRangeUser(-1.2,1.2);
    if( bstyle[0] ) {
      mrv1->GetXaxis()->SetRangeUser(-0.65,1.2);
      mrv1->GetYaxis()->SetRangeUser(-0.4, 0.6);
    }

    mrv1->Draw("ALP");
    lgr1->Draw();

    Ymin = mrv1->GetYaxis()->GetXmin();
    Ymax = mrv1->GetYaxis()->GetXmax();
    Xmin = mrv1->GetXaxis()->GetXmin();
    Xmax = mrv1->GetXaxis()->GetXmax();

    auto aLineX1 = new TLine(Xmin, 0., Xmax, 0.);
    aLineX1->SetLineColor(1);
    aLineX1->SetLineStyle(3);
    aLineX1->Draw();

    auto aLineY1 = new TLine(0., Ymin, 0., Ymax);
    aLineY1->SetLineColor(1);
    aLineY1->SetLineStyle(3);
    aLineY1->Draw();

    //--- v2 vs rapidity ---
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

    mrv2->GetXaxis()->SetRangeUser(-1.2,1.2);

    if( bstyle[0] )
      mrv2->GetXaxis()->SetRangeUser(-0.55,1.2);
    
    mrv2->Draw("ALP");
    lgr2->Draw();

    if( bstyle[2] ) {
      fv2y->Draw("same");
    }

    Ymin = mrv2->GetYaxis()->GetXmin();
    Ymax = mrv2->GetYaxis()->GetXmax();
    Xmin = mrv2->GetXaxis()->GetXmin();
    Xmax = mrv2->GetXaxis()->GetXmax();

    auto aLineX2 = new TLine(0., Ymin, 0., Ymax);
    aLineX2->SetLineColor(1);
    aLineX2->SetLineStyle(3);
    aLineX2->Draw();


    if( bstyle[2] ) {
      ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
      fv2y->SetTitle(";y_{nrm}; v2");
      fv2y->SetLineColor(3);
      fv2y->Draw("");
      mrv2->Draw("");
      lgr2->AddEntry("fv2y","initial","lp");
      lgr2->Draw();
    }
  }

  if( bplot[0] && bplot[5] ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),500,700);
    cout << "ngr " << ngr << endl;
    cc->Divide(1,ngr);
    for(UInt_t i = 0; i < ngr; i++){
      if( hyptacp[i] == NULL) {
	cout << " NULL " << endl;
	continue;
      }
      else {
       	cc->cd(i+1);
	cout << hyptacp[i] << endl;
	hyptacp[i]->Draw("colz");
      }
    }
  }

  if( (bplot[0] || bplot[1]) && bplot[3] ) { // gethered plots
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1500,1200);//2000,500);
    cc->Divide((ybin1-1)/3,3);
    for(UInt_t k = 0; k < ybin1-1; k++){
      cc->cd(k+1);
      mv1[k]->Draw("ALP");
      if( k == ybin1-2 )
	lg1[k]->Draw();
      fOutput->cd();
      mv1[k]->Write();
    }

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1500,1200);//1400,500);
    cc->Divide((ybin2-1)/3,3);
    for(UInt_t k = 0; k < ybin2-1; k++){
      cc->cd(k+1);
      mv2[k]->Draw("ALP");
      if( k == ybin2-2 ) lg2[k]->Draw();
      fOutput->cd();
      mv2[k]->Write();
    }
  }


  if( bplot[0] && bplot[4] ) {  // individual pt plots  

    UInt_t ichoise[]={0,1,4};
    for(UInt_t i : ichoise) {
      cc = new TCanvas(Form("cv%d",i),Form("cv%d",i),500,550);
      cc->SetRightMargin(0.02);
      cc->SetLeftMargin(0.15);
      cc->SetTopMargin(0.08);
      cc->SetBottomMargin(0.15);

      mv1[i]->GetXaxis()->SetNdivisions(505);
      mv1[i]->GetXaxis()->SetLabelSize(0.04);
      mv1[i]->GetYaxis()->SetLabelSize(0.04);
      mv1[i]->GetXaxis()->SetTitleSize(0.05);
      mv1[i]->GetYaxis()->SetTitleSize(0.04);
      mv1[i]->GetYaxis()->SetTitleOffset(2.5);

      //      mv1[i]->GetXaxis()->SetRangeUser(0.,1.5);
      mv1[i]->Draw("ALP");
      lg1[i]->Draw();

      fOutput->cd();
      mv1[i]->Write();
    }

    UInt_t ichoise2[] = {0,1,2,3,4,5};
    for(UInt_t i : ichoise2) {
      cc = new TCanvas(Form("cv%d",i+10),Form("cv%d",i+10),500,550);
      cc->SetRightMargin(0.02);
      cc->SetLeftMargin(0.15);
      cc->SetTopMargin(0.08);
      cc->SetBottomMargin(0.15);

      mv1[i]->GetXaxis()->SetNdivisions(505);
      mv2[i]->GetXaxis()->SetLabelSize(0.04);
      mv2[i]->GetYaxis()->SetLabelSize(0.04);
      mv2[i]->GetXaxis()->SetTitleSize(0.05);
      mv2[i]->GetYaxis()->SetTitleSize(0.04);
      mv2[i]->GetYaxis()->SetTitleOffset(2.5);

      //      mv2[i]->GetXaxis()->SetRangeUser(0.,1.5);
      //      mv2[i]->Fit("fv2fit","","",0.,1.1);
      mv2[i]->Draw("ALP");
      lg2[i]->Draw();

      fOutput->cd();
      mv2[i]->Write();
    }
  }

  if( kFALSE ) { //bplot[0] && bplot[7] ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    g_v1slp->SetMarkerStyle(20);
    g_v1slp->SetMarkerColor(2);
    g_v1slp->Draw("AP");
    //    g_v1slp->Print();
  }

  if( bplot[0] && bplot[7]) {
    auto m_v1sys = new TMultiGraph();
    auto lgv1sys = new TLegend(0.2,0.7,0.5,0.9,"");
    m_v1sys->SetTitle(";(n-p)/A; v_{11}");

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    for(UInt_t jj = 0; jj < 6; jj++){
      cout << " g_v1slpsys[jj]->GetN() " << g_v1slpsys[jj]->GetN() << endl;
      if(g_v1slpsys[jj]->GetN() > 0) {
	g_v1slpsys[jj]->SetMarkerStyle(20);
	g_v1slpsys[jj]->SetMarkerSize(1.5);
	g_v1slpsys[jj]->SetMarkerColor(icol[jj]);
	g_v1slpsys[jj]->SetLineColor(icol[jj]);
	
	m_v1sys->Add(g_v1slpsys[jj]);
	lgv1sys->AddEntry(g_v1slpsys[jj],lpid[jj],"lp");	
      }
    }
    m_v1sys->Draw("AP");
    lgv1sys->Draw();
  }


  if( bplot[0] && bplot[9]) {
    auto m_v1m = new TMultiGraph();
    auto lgv1slpm =  new TLegend(0.2,0.2,0.5,0.4,"");


    m_v1m->SetTitle(";"+ hmultX +"; v_{11}");

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    for(UInt_t jj = 0; jj < 20; jj++){
      if(g_v1slpm[jj]->GetN() > 0) {
	g_v1slpm[jj]->SetMarkerStyle(20);
	g_v1slpm[jj]->SetMarkerSize(0.9);
	g_v1slpm[jj]->SetMarkerColor(icol[jj]);
	g_v1slpm[jj]->SetLineColor(icol[jj]);
	
	m_v1m->Add(g_v1slpm[jj]);
	lgv1slpm->AddEntry(g_v1slpm[jj],g_v1slpm[jj]->GetTitle(),"lp");		
      }
    }
    m_v1m->Draw("AP");
    lgv1slpm->Draw();
  }

  if( kFALSE ) { //bplot[0] &&  bplot[9] && g_v2max->GetN()>0) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    g_v2max->SetMarkerStyle(20);
    g_v2max->SetMarkerColor(2);
    g_v2max->Draw("AP");
    //    g_v2max->Print();
  }
    
  if( bplot[0] && bplot[8] ){

    auto m_v2sys  = new TMultiGraph();
    auto lgv2sys  =  new TLegend(0.2,0.7,0.5,0.9,"");
    m_v2sys->SetTitle(";(n-p)/A; -v2 max");

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

    for(UInt_t jj = 0; jj < 6; jj++){
      if(g_v2maxsys[jj]->GetN() > 0) {
	g_v2maxsys[jj]->SetMarkerStyle(20);
	g_v2maxsys[jj]->SetMarkerSize(1.);
	g_v2maxsys[jj]->SetMarkerColor(icol[jj]);
	g_v2maxsys[jj]->SetLineColor(icol[jj]);

	lgv2sys->AddEntry(g_v2maxsys[jj],lpid[jj],"lp");	
	m_v2sys->Add(g_v2maxsys[jj]);

      }
    }
    m_v2sys->Draw("AP");
    lgv2sys->Draw();
  }

  if( bplot[0] && bplot[9] ){

    auto m_v2m  = new TMultiGraph();
    auto lgv2m  =  new TLegend(0.5,0.7,0.8,0.9,"");
    m_v2m->SetTitle(";"+hmultX+"; -v2 max");

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

    for(UInt_t jj = 0; jj < 20; jj++){
      if(g_v2maxm[jj]->GetN() > 0) {
	cout << "m_v2m " << jj << " " << icol[jj] << endl;
    
	//	if( jj == 0 )
	g_v2maxm[jj]->SetMarkerSize(0.9);
	g_v2maxm[jj]->SetMarkerStyle(20);
	g_v2maxm[jj]->SetMarkerColor(icol[jj]);
	g_v2maxm[jj]->SetLineColor(icol[jj]);

	lgv2m->AddEntry(g_v2maxm[jj], g_v2maxm[jj]->GetTitle(),"lp");	
	m_v2m->Add(g_v2maxm[jj]);
      }
    }
    m_v2m->Draw("AP");
    lgv2m->Draw();
  }
  

  if( bplot[0] && bplot[6] ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    auto msysD = new TMultiGraph("msysD","; (N-P)/A; F[MeV/c]");
    g_mpxslp->SetMarkerStyle(20);
    g_mpxslp->SetMarkerColor(2);
    msysD->Add(g_mpxslp, "LP");
    lgr0->AddEntry(g_mpxslp,"DATA x 1.4");

    msysD->Draw("ALP");
    lgr0->Draw();

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    mmpx->Draw("ALP");
    lgr9->Draw();

    auto Ymin = mmpx->GetYaxis()->GetXmin();
    auto Ymax = mmpx->GetYaxis()->GetXmax();
    auto Xmin = mmpx->GetXaxis()->GetXmin();
    auto Xmax = mmpx->GetXaxis()->GetXmax();

    auto aLineX3 = new TLine(Xmin, 0., Xmax, 0.);
    aLineX3->SetLineColor(1);
    aLineX3->SetLineStyle(3);
    aLineX3->Draw();

    auto aLineY3 = new TLine(0., Ymin, 0., Ymax);
    aLineY3->SetLineColor(1);
    aLineY3->SetLineStyle(3);
    aLineY3->Draw();


    if( g_v2sysD->GetN() > 1 ) {
      auto mv2sysD = new TMultiGraph("mv2sysD","; (N-P)/A; v2");
      g_v2sysD->SetMarkerStyle(20);
      g_v2sysD->SetMarkerColor(2);
      g_v2sysD->SetLineColor(4);
      mv2sysD->Add(g_v2sysD,"LP");
      lgr10->AddEntry(g_mpxslp,"DATA");

      ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
      mv2sysD->Draw("AP");
      lgr10->Draw();
    }
  }


  for( UInt_t j = 0; j < (UInt_t)lx_v11.size(); j++ ) {
    auto fbin = g_v11 -> GetXaxis() -> FindBin( j+0.5 );
    if( fbin != 101 )
      g_v11 -> GetXaxis() -> SetBinLabel( fbin, lx_v11.at(j) );
  }

  if( bplot[12] && g_v11 != NULL) {

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    cc  -> GetPad(0)->SetBottomMargin(0.30);

    g_v11->SetMarkerColor(3);
    g_v11->SetMarkerStyle(20);
    g_v11->Draw("AP");

    cout << "d_v11.size() = " << d_v11.size() << endl;
    if( d_v11.size() > 0 ) {
      auto Xmin = g_v11->GetXaxis()->GetXmin();
      auto Xmax = g_v11->GetXaxis()->GetXmax();
      auto aLinev11 = new TLine(Xmin, d_v11.at(0), Xmax, d_v11.at(0));
      aLinev11->SetLineColor(2);
      aLinev11->SetLineStyle(3);
      auto aBoxv11 = new TBox(Xmin, d_v11.at(0)+d_v11e.at(0), Xmax, d_v11.at(0)-d_v11e.at(0) ); 
      aBoxv11->SetFillStyle(3008);
      aBoxv11->SetFillColor(4);
      aLinev11->Draw();
      aBoxv11->Draw();
    }


    for( UInt_t j = 0; j < (UInt_t)lx_v20.size(); j++ ) {
      auto fbin = g_v20 -> GetXaxis() -> FindBin( j+0.5 );
      if( fbin != 101 )
	g_v20 -> GetXaxis() -> SetBinLabel( fbin, lx_v20.at(j) );
    }

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    cc  -> GetPad(0)->SetBottomMargin(0.30);

    g_v20->SetMarkerColor(3);
    g_v20->SetMarkerStyle(20);
    g_v20->Draw("AP");

    if( d_v20.size() > 0 ) {
      auto Xmin = g_v20->GetXaxis()->GetXmin();
      auto Xmax = g_v20->GetXaxis()->GetXmax();
      auto aLinev20 = new TLine(Xmin, d_v20.at(0), Xmax, d_v20.at(0));
      aLinev20->SetLineColor(2);
      aLinev20->SetLineStyle(3);
      auto aBoxv20 = new TBox(Xmin, d_v20.at(0)+d_v20e.at(0), Xmax, d_v20.at(0)-d_v20e.at(0) ); 
      aBoxv20->SetFillStyle(3008);
      aBoxv20->SetFillColor(4);
      aLinev20->Draw();
      aBoxv20->Draw();
    }


  }

}


void RapidityShift(TGraphErrors *gv) 
{
  // shifting artificially
  std::cout << "===> SHIFTING " << 0.018 << std::endl;
  Double_t xx, yy;
  for(UInt_t k = 0; k < (UInt_t)gv->GetN(); k++) {
    gv->GetPoint(k, xx, yy);
    xx += 0.018;
    gv->SetPoint(k, xx, yy);
  }
}

void ShiftX(TGraphErrors *gv, Double_t off)
{
  Double_t xx, yy;
  for(Int_t k = (Int_t)gv->GetN()-1; k > -1; k--) {

    gv->GetPoint(k, xx, yy);
    xx += off;
    gv->SetPoint(k, xx, yy);
  }
}

void integ()
{
  UInt_t idx = 0;
  while(1) {
    TH1D* ff = (TH1D*)gROOT->Get(Form("hyptacp_%d_px",idx));
    if( ff != NULL ) {
      ff->Print();
      cout << " -> " <<  ff->Integral()/200. << endl;;
      //      cout << " -> " <<  ff->Integral()/(Double_t)ff->GetXaxis()->GetNbins() << endl;;
      idx++;
    }
    else break;

    //    hyptacp_0_px->Integral()/200+hyptacp_1_px->Integral()/200+hyptacp_2_px->Integral()/200;
    //    hyptacp_3_px->Integral()/200+hyptacp_4_px->Integral()/200+hyptacp_5_px->Integral()/200;
  }
}


void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer)
{
  std::cout << "GetMinimum v2 " << gr->GetN() << std::endl;
  gr->Print();

  Double_t x, y;
  min = 9.;
  UInt_t mid = 0;
  for(UInt_t i = 0; i < gr->GetN(); i++) {
    gr->GetPoint(i, x, y);
    if( abs(y) > 0.2 ) {
      //      gr->RemovePoint(i);
      //      i = 0;
      continue;
    }

    if( y < min && y < 0) {
      min = y;
      mid = i;

      std::cout << " min v2 = " << min << " at " << mid << std::endl;
    }
  }


  mer = gr->GetErrorY(mid);
}

void GetAveragev2(TGraphErrors *gr)
{
  Double_t x, y;
  Double_t sum = 0.;
  for(UInt_t i = 1; i < gr->GetN(); i++) {
    gr->GetPoint(i, x, y);

    sum += y;
    
  }

  Double_t ave = sum/(Double_t)(gr->GetN()-1);

  //  double v2in;
  //  std::cout << " What is initial v2? " ;
  //  cin  >> v2in;

  std::cout << " Average v2 " << ave 
	    << setprecision(4) 
    //	    << " Deviation " << (ave - v2in)/v2in *100. << " % "
    //	    << " Original v2 " << v2in
	    << std::endl;

  //  fv2y->SetParameter(0,  v2in);
  //  fv2y->SetParameter(1,  0.);
  //  fv2y->SetParameter(2,  0.);


}

void FindFirstBin( TH1I *h, Int_t &min, Int_t &max )
{
  for( UInt_t i = 0; i < h->GetXaxis()->GetNbins(); i++ ) {
    if( h->GetBinContent(i) > 0 ){
      min = i;
      cout << "min " << min << endl;
      break;
    }
  }
  for( UInt_t i = h->GetXaxis()->GetNbins(); i > 0;  i-- ) {
    if( h->GetBinContent(i) > 0 ) {
      max = i;
      cout << "max " << max << endl;
      break;
    }
  }
}

void DividebyTwo(TGraphErrors *yv ) 
{
  for( UInt_t i = 0; i < yv->GetN(); i++ ) {
    Double_t xx, yy;
    yv->GetPoint(i, xx, yy);
    yv->SetPoint(i, xx, yy/2.);

    xx = yv->GetErrorX(i);
    yy = yv->GetErrorY(i);
    yv->SetPointError(i, xx, yy/2.);
  }
}
