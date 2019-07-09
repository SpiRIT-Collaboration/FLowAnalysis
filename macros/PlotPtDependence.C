#include "DoFlow.h"
#include "SetStyle.C"

struct gplot{
  TString Version;
  TString fileHeader;
  TString comment;
};


TString  bName[]   = {"132Sn_","108Sn_","124Sn_","112Sn_"};
Double_t sysdlt[]  = {0.22,    0.09,      0.15,   0.15};
Double_t sysA[]    = {256.,    220.,      236.,   236};

//==================================================
//-- plot configuration
//--------------------------------------------------
  // --> Plotting selection
//--- Data
Bool_t bsys[]  = { 1, 0, 0, 0};
Bool_t bpid[]  = { 1, 0, 0, 0, 0, 0, 0}; //0:p, 1:d, 2:t, 3:3He, 4:4He, 5:n 6:H
Bool_t bcnt[]  = { 1, 0, 0}; 
UInt_t cntw = 1;
UInt_t iv2at = 4;
//-----------

UInt_t  bver[]  = {0, 0, 0, 0, 0, 0, 0, 0};
const UInt_t nmax = (UInt_t)sizeof(bver)/sizeof(UInt_t);
gplot gnames[] = { 
  {".v29.1.24" ,"advYPt_","m5-55"} ,
  {".v29.1.22" ,"advYPt_","m5-55"} ,
  {".v29.1.21" ,"advYPt_","yaw<0"} , //"m5-80"} ,
  {".v29.1.3"  ,"advYPt_","yaw<0"} , //"m5-80"} ,
  {".v29.1.15" ,"advYPt_","m<2790"} , 
  {".v29.1.20" ,"advYPt_","p>1000 yaw>0"} , 
  {".v29.1.19" ,"advYPt_","p>1000 yaw<0"} , 
  {".v29.1.16" ,"advYPt_","m>2790"} , 
  {".v29.1.0"  ,"advYPt_","yaw>0"} , 
  {".v29.1.12" ,"advYPt_","m5-80"} , 
  {".v29.1.9"  ,"advYPt_","m5-80"} , 
  {".v29.1.1"  ,"advYPt_","m5-55"} ,
  {".v29.1.2"  ,"advYPt_","m5-40"} ,
  {".v29.1.4"  ,"advYPt_","yaw>0&pitch>0"},//"m5-40"} ,
  {".v29.1.5"  ,"advYPt_","yaw>0&pitch<0"},//"m5-40"} ,
  {".v29.1.6"  ,"advYPt_","yaw>0&pitch<0"}, //"m5-80"} ,
  {".v29.1.7"  ,"advYPt_","yaw>0&pitch>0"} //"m5-80"} ,
};

TString sVer[nmax];
TString sName[nmax];
TString cmnt[nmax];

//-- pBUU
Bool_t  bpBUUConfig[]  = {1, 0, 0, 0};
//-----------
gplot pBUUConfig[] = {
  {"_energy270_gamma0.50_b5_withCluster.","pBUU2","pBUU b5g0.5"},
  {"_energy270_gamma1.75_b5_withCluster.","pBUU2","pBUU b5g1.75"},
  {".b2.g05"   ,"pBUU_"  ,"pBUU b2g0.5"},
  {".b2.g175"  ,"pBUU_"  ,"pBUU b2g1.75"},
  {".b4.g05"   ,"pBUU_"  ,"pBUU b4g0.5"},
  {".b4.g175"  ,"pBUU_"  ,"pBUU b4g1.75"}
};
Double_t pBUU_F[4][4] = { //F b2 g0.5, b2 g1.75, b4 g0.5, b4 g1.75
  { 92.1902, 93.199, 120.6029, 120.467},   // 132
  { 88.2085, 87.8252,113.6998, 113.6738},  // 108
  { 88.2085, 87.8252,113.6998, 113.6738},  // 108 (124)
  { 88.2085, 87.8252,113.6998, 113.6738}   // 108 (112)
};

Double_t pBUU_v2[4][4] = { //v2 b2 g0.5, b2 g1.75, b4 g0.5, b4 g1.75
  { -0.0071928, -0.0058827, -0.0372598, -0.0338454}, //132
  { -0.0066653, -0.0059764, -0.0372598, -0.0343415}, //108
  { -0.0066653, -0.0059764, -0.0372598, -0.0343415}, //(124)
  { -0.0066653, -0.0059764, -0.0372598, -0.0343415}  //(112)
};

//-- AMD
Bool_t  amdEOS[]= {1, 0};
//-----------
TString amdName[] = {"SLy4",
		     "SLy4-L108"};

TString amdHeader[] = {"amd_132Sn124Sn270AMeV_cluster_",
		       "amd_108Sn112Sn270AMeV_cluster_"};

//==================================================
Bool_t bplot[] = 
  { 0, // 0 data   It should be set to 1 in the code.
    0, // 1 model  It should be set to 1 in the code.
    0, // 2 v1 and v2 rapidity 
    1, // 3 v1 and v2 on pt in one window
    0, // 4 v1 and v2 in individual windows
    0, // 5 Acceptance ypt
    0, // 6 v1 on pt special
    0, // 7 <px>/A
  };
//==================================================

UInt_t   ccvid = 0;
TF1 *lslope = new TF1("lslope","[0]+[1]*x",-1., 1.);

// --> configuration
Size_t  imsz[] = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
Color_t icol[] = {  kRed, kBlue, kSpring, kMagenta, kOrange, kViolet};
Color_t icolnn[]={ kPink, kBlue+1, kGreen+2, kViolet-1};


TMultiGraph  *mv1[ybin1];
TMultiGraph  *mv2[ybin2];
TLegend      *lg1[ybin1];
TLegend      *lg2[ybin2];

void RapidityShift(TGraphErrors *gv);
void ShiftX(TGraphErrors *gv, Double_t off);


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

  for(UInt_t i = 0; i < nmax; i++){
    if( bver[i] ) {
      sVer[i]  = gnames[i].Version;
      sName[i] = gnames[i].fileHeader;
      cmnt[i]  = gnames[i].comment;
    }
  }

  TString ltitle;


  UInt_t ngr = 0;
  std::vector< TString > fname;
  std::vector< std::vector< UInt_t >> index(6);

  TGraphErrors *g_pBUUv2g5 = new TGraphErrors();
  g_pBUUv2g5->SetName("g_pBUUv2g5");
  g_pBUUv2g5->SetTitle("; (N-P)/A; v2");

  TGraphErrors *g_pBUUv2g17 = new TGraphErrors();
  g_pBUUv2g17->SetName("g_pBUUv2g17");
  g_pBUUv2g17->SetTitle("; (N-P)/A; v2");


  TGraphErrors *g_pBUUFg5 = new TGraphErrors();
  g_pBUUFg5->SetName("g_pBUUFg5");
  g_pBUUFg5->SetTitle("; (N-P)/A; F{MeV/c]");

  TGraphErrors *g_pBUUFg17 = new TGraphErrors();
  g_pBUUFg17->SetName("g_pBUUFg17");
  g_pBUUFg17->SetTitle("; (N-P)/A; F{MeV/c]");


  TGraphErrors *g_mpxslp = new TGraphErrors();
  g_mpxslp->SetName("g_mpxslp");
  g_mpxslp->SetTitle("; (N-P)/A; F[MeV/c]");

  TGraphErrors *g_v1slp = new TGraphErrors();
  g_v1slp->SetName("g_v1slp");
  g_v1slp->SetTitle("; Multiplicity; slope of v1");

  TGraphErrors *g_v2max = new TGraphErrors();
  g_v2max->SetName("g_v2max");
  g_v2max->SetTitle(";multiplicity; v2_max");

  TGraphErrors *g_v2sysA = new TGraphErrors();
  g_v2sysA->SetName("g_v2sysA");
  g_v2sysA->SetTitle(";A_{beam} + A_{target}; v2");

  TGraphErrors *g_v2sysD = new TGraphErrors();
  g_v2sysD->SetName("g_v2sysD");
  g_v2sysD->SetTitle(";(N-P)/A; v2");

  // Rapidity dependence 
  TH2D *hyptacp[10];

  auto mrv1  = new TMultiGraph("mrv1"  ,";y_{cm}/y_{beam}; v1");
  auto mrv2  = new TMultiGraph("mrv2"  ,";y_{cm}/y_{beam}; v2");
  auto mv1sl = new TMultiGraph("mv1sl" ,";Centrality; v1 slope");
  auto mv1slp= new TMultiGraph("mv1slp","; Particle ; v1 slope");
  auto mmpx  = new TMultiGraph("mmpx","; y/y_{cm} ; <px>/A");
  //  mv1slp->GetXaxis()->SetAlphanumeric(kTRUE);

  auto lgr1 = new TLegend(0.54, 0.13, 0.87, 0.4, ""); 
  auto lgr2 = new TLegend(0.38, 0.68, 0.75, 0.9, "");
  auto lgr3 = new TLegend(0.35, 0.13, 0.7, 0.33, "");
  auto lgr4 = new TLegend(0.15, 0.63, 0.5, 0.85, "");
  auto lgr5 = new TLegend(0.15, 0.63, 0.5, 0.85, "");
  auto lgr6 = new TLegend(0.72, 0.72, 0.94, 0.88,"");
  auto lgr7 = new TLegend(0.16, 0.70, 0.46, 0.85,"");
  auto lgr8 = new TLegend(0.16, 0.70, 0.46, 0.85,"");
  auto lgr9 = new TLegend(0.16, 0.70, 0.46, 0.85,"");
  auto lgr0 = new TLegend(0.16, 0.70, 0.46, 0.85,"");
  auto lgr10= new TLegend(0.16, 0.70, 0.46, 0.85,"");
  

  for(UInt_t is = 0; is < 4; is++){

    for(UInt_t ip = 0; ip < (UInt_t)sizeof(bpid)/sizeof(Bool_t); ip++){

      if( !bsys[is] || !bpid[ip] ) continue;

      for(UInt_t it = 0; it < (UInt_t)sizeof(bcnt)/sizeof(Bool_t); it++){
	
	if( !bcnt[it] ) continue;

	for(UInt_t iz = 0; iz < (UInt_t)sizeof(bver)/sizeof(UInt_t); iz++){
	  
	  if( bver[iz] == 0 ) continue;

	  fname.push_back( "data/"+ sName[iz] + bName[is] + fpid[ip] + sVer[iz] + ".root" );

	  std::cout << fname.at(ngr) << std::endl;
	  ltitle = "";
	  //	  ltitle = Form("mult %d ~ %d",cent[it],cent[it+cntw]);
	  //	  std::cout << " label title " << ltitle << std::endl;

	  index[0].push_back(is); // system
	  index[1].push_back(ip); // pid
	  index[2].push_back(it); // centrality
	  index[3].push_back(it+cntw);
	  index[4].push_back(iz);
	  index[5].push_back(1);
	  ngr++;

	  bplot[0] = 1;
	}
      }
    }
  }


  cout << fname.size() << endl;

  //---- amd
  UInt_t kgr = 0;
  for(UInt_t ls = 0; ls < 2; ls++){
    if( !bsys[ls] ) continue;

    for(UInt_t lp = 0; lp < (UInt_t)sizeof(bpid)/sizeof(Bool_t); lp++){
      if( !bpid[lp] ) continue;

      for(UInt_t lo = 0; lo < (UInt_t)sizeof(amdEOS)/sizeof(Bool_t); lo++) { 
	if( !amdEOS[lo] ) continue;

	fname.push_back( "data/"+ amdHeader[ls]+amdName[lo]+"_"+amdpartname[lp]+"_cm.root" );
	std::cout << fname.at(ngr+kgr) << std::endl;

	index[0].push_back(ls);
	index[1].push_back(lp);
	index[2].push_back(lo);
	index[3].push_back(3);
	index[4].push_back((UInt_t)sizeof(bver)/sizeof(UInt_t)+1);
	index[5].push_back(2);
	kgr++;

      }
    }
  }
  
  //---- pBUU
  UInt_t pgr = 0;
  for(UInt_t ms = 0; ms < 4; ms++ ){
    if( !bsys[ms] ) continue;
    
    for(UInt_t ip = 0; ip < (UInt_t)sizeof(bpid)/sizeof(Bool_t); ip++){
      cout << ip << " mp " << bpid[ip] << endl;

      if( !bpid[ip] ) continue;
      //      if( (!bpid[ip] && ip < 5) || (!bpid[6] && ip == 6)  ) continue;

      for(UInt_t mt = 0; mt < (UInt_t)sizeof(bpBUUConfig)/sizeof(Bool_t); mt++) {
	if( !bpBUUConfig[mt] ) continue;
	
	//@@@
	fname.push_back("../../pBUU/data/"+pBUUConfig[mt].fileHeader+"sn" + rsys[ms] + "_sn" 
			+ tsys[ms] + pBUUConfig[mt].Version + fpid[ip] + ".root");
	
	std::cout << fname.at(ngr+kgr+pgr) << std::endl;
	index[0].push_back(ms);
	index[1].push_back(ip);
	index[2].push_back(0);
	index[3].push_back(mt);
	index[4].push_back(1);
	index[5].push_back(3);
	pgr++;

	bplot[1] = 1;
      }
    }
  }

  //---------------------------------------------------

  // pBUU by Genie
  if( bplot[1] ) {
    g_pBUUFg5->SetPoint(0, sysdlt[0], pBUU_F[0][2]);
    g_pBUUFg5->SetPoint(1, sysdlt[1], pBUU_F[1][2]);
    
    g_pBUUFg17->SetPoint(0, sysdlt[0], pBUU_F[0][3]);
    g_pBUUFg17->SetPoint(1, sysdlt[1], pBUU_F[1][3]);

    g_pBUUv2g5->SetPoint(0, sysdlt[0], pBUU_v2[0][2]);
    g_pBUUv2g5->SetPoint(1, sysdlt[1], pBUU_v2[1][2]);
    
    g_pBUUv2g17->SetPoint(0, sysdlt[0], pBUU_v2[0][3]);
    g_pBUUv2g17->SetPoint(1, sysdlt[1], pBUU_v2[1][3]);
  }  


  // TString rapRange[ybin1];
  // for(UInt_t i = 1; i < ybin1 - 1; i++)
  //   rapRange[i] = Form(" %f <= Y_cm < %f", yrange1[i-1],yrange1[i]);
  // rapRange[0] = Form(" Y_cm < %f", yrange1[0]);
  // rapRange[ybin1-1] = Form(" %f <= Y_cm", yrange1[ybin1-2]);

  // TString rapRange2[ybin2];
  // for(UInt_t i = 1; i < ybin2 - 1; i++)
  //   rapRange2[i] = Form(" %f <= Y_cm < %f", yrange2[i-1],yrange2[i]);
  // rapRange2[0] = Form(" Y_cm < %f", yrange2[0]);
  // rapRange2[ybin2-1] = Form(" %f <= Y_cm", yrange2[ybin2-2]);


  TFile *fOpen;


  // Pt dependence setup
  for(UInt_t k = 0; k < ybin1; k++){
    mv1[k]  = new TMultiGraph((TString)Form("mv1%d",k),  ";Pt [MeV/c]; v_1");
    lg1[k] = new TLegend(0.18, 0.2, 0.6, 0.4); 
    lg1[k]->SetTextSize(0.05);
  }

  for(UInt_t k = 0; k < ybin2; k++) {
    mv2[k]  = new TMultiGraph((TString)Form("mv2%d",k),  ";Pt [MeV/c]; v_2");
    lg2[k] = new TLegend(0.2, 0.15, 0.65, 0.43); 
    lg2[k]->SetTextSize(0.05);
  }

  UInt_t islp = 0;
  UInt_t isp = 0;
  UInt_t isk = 0;
  UInt_t iss = 9;
  UInt_t ipp = 0;
  Color_t icolor = 100; 
  for(UInt_t igr = 0; igr < ngr+kgr+pgr; igr++ ) {

    fOpen = TFile::Open(fname.at(igr)); 

    if( fOpen == NULL ) continue;
    else
      std::cout << fname.at(igr) << " is opened. " << std::endl;

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

    if( igr < 7)
      icolor = icol[igr];
    else
      icolor -= 2;

    cout << igr  << " " << ip << " color " << icolor << endl;
    TString ohtitle = lsys[is]+" "+fpid[ip];
    TString otitle  = ohtitle;

    if( ia == 1 ) { //data
      otitle += sName[iz]+sVer[iz]+";"+cmnt[iz];
    }
    else if( ia == 2 ) //amd
      otitle += amdHeader[is](4,5) + amdName[it];
    else if( ia == 3 ) {// pBUU
      otitle += pBUUConfig[ik].comment; 
    }


    //acceptance
    if( igr < ngr ) {
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

    //------------------------------    
    // rapidity dependence
    TH1I *hmult  = (TH1I*)fOpen->Get("hmult");

    TGraphErrors *yv1 = (TGraphErrors*)fOpen->Get("gv_v1");
    if( yv1 != NULL ) {
      yv1->SetMarkerColor(icolor);

      cout << "imark[" << is << "] " << imark[is] << endl;
      yv1->SetMarkerStyle(imark[is]);

      yv1->RemovePoint(10);
      yv1->RemovePoint(9);
      yv1->RemovePoint(8);
      yv1->RemovePoint(7);

      if( ip == 5 )
	yv1->SetMarkerStyle(imark[is]+4);
      yv1->SetMarkerSize(imsz[is]);
      yv1->SetLineColor(icolor);
    
      mrv1->Add(yv1,"lp");
      lgr1->AddEntry(yv1,  ohtitle ,"lp");

      yv1->Fit("lslope","Q0","",-0.4,0.4);
      Double_t constlslope = lslope->GetParameter(0);
      Double_t slope = lslope->GetParameter(1);
      Double_t slopee= lslope->GetParError(1);
      
      if( hmult != NULL ) {
	Double_t mmean = hmult->GetMean();
	Double_t mstd  = hmult->GetStdDev() / sqrt( (Double_t)hmult->GetEntries() );
	
	g_v1slp -> SetPoint(islp, mmean, slope);
	g_v1slp -> SetPointError(islp, mstd, slopee);	
      
      }
    }

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

    //--- y vs v2 ---
    TGraphErrors *yv2 = (TGraphErrors*)fOpen->Get("gv_v2");

    if( yv2 != NULL ) {
      //    ShiftX(yv2, 0.01*iz);
      //    RemoveYPoint(yv2);

      yv2->RemovePoint(7);
      yv2->RemovePoint(6);
      yv2->RemovePoint(5);

      yv2->SetMarkerColor(icolor);
      yv2->SetMarkerStyle(imark[is]);
      if( ip == 5 )
	yv2->SetMarkerStyle(imark[is]+4);
      yv2->SetMarkerSize(imsz[is]);
      yv2->SetLineColor(icolor);

      mrv2->Add(yv2,"lp");
      lgr2->AddEntry(yv2,  ohtitle ,"lp");
      // --end of rapidity dependence

      Double_t v2x, v2y, v2xe;
      yv2->GetPoint(iv2at, v2x, v2y);
      Double_t v2ye = yv2->GetErrorY(iv2at);

      if( hmult != NULL ) {
	Double_t mmean = hmult->GetMean();
	Double_t mstd  = hmult->GetStdDev();

	g_v2max->SetPoint(isp , mmean, v2y);
	g_v2max->SetPointError(isp, mstd, v2ye);
	isp++;
      }

      if( bsys[ik] ) {
	g_v2sysA->SetPoint(isk, sysA[is], v2y);
	g_v2sysA->SetPointError(isk, 0., v2ye);
	
	g_v2sysD->SetPoint(isk, sysdlt[is], v2y);
	g_v2sysD->SetPointError(isk,  0.,  v2ye);
	isk++;
      }
    }

    //--------------------------------------------------
    // plotting
    //--------------------------------------------------
    // Pt dependence 
    Double_t v1mx = -99.;
    Double_t v1mn =  99.;
    Double_t v2mx = -99.;
    Double_t v2mn =  99.;

    for(UInt_t k = 0; k < ybin1; k++){

      TGraphErrors *gr_v1 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v1%d",k));
      if( gr_v1 == NULL ) continue;
      if( igr == 0 ) {
	mv1[k]->SetTitle(gr_v1->GetTitle());
	mv1[k]->GetXaxis()->SetTitle("Pt/A {MeV/c]");
	mv1[k]->GetYaxis()->SetTitle("v1/A");
      }

      TGraphErrors *gr_v1A = new TGraphErrors();
      for(UInt_t j = 0; j < (UInt_t)gr_v1->GetN(); j++ ){
	Double_t x, y;
	gr_v1->GetPoint(j, x, y);
	x /= partA[ip+2];
	y /= partA[ip+2];
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
      lg1[k]->AddEntry(gr_v1A,  ohtitle,"lp");

    }
    
    cout << " ==== " << endl;

    for(UInt_t k = 0; k < ybin2; k++){

      TGraphErrors *gr_v2 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v2%d",k));
      if( gr_v2 == NULL ) continue;

      if( igr == 0 ) {
	mv2[k]->SetTitle(gr_v2->GetTitle());
	mv2[k]->GetXaxis()->SetTitle("Pt/A {MeV/c]");
	mv2[k]->GetYaxis()->SetTitle("v2/A");
      }

      TGraphErrors *gr_v2A = new TGraphErrors();
      for(UInt_t j = 0; j < (UInt_t)gr_v2->GetN(); j++ ){
	Double_t x, y;
	gr_v2->GetPoint(j, x, y);
	x /= partA[ip+2];
	y /= partA[ip+2];
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
      lg2[k]->AddEntry(gr_v2A, ohtitle ,"lp");

      //      mv2[k]->GetYaxis()->SetRangeUser(v2mn, v2mx);
    }
	
    fOpen->Close();
  }


  if( bplot[0] && bplot[5] ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),500,200);
    cc->Divide(1,ngr);
    for(UInt_t i = 0; i < ngr; i++){
      if( hyptacp[i] == NULL) continue;
      cc->cd(i+1);
      hyptacp[ngr]->Draw("colz");
    }
  }

  if( (bplot[0] || bplot[1]) && bplot[3] ) { // gethered plots
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1400,900);
    cc->Divide(ybin1/2,2);
    for(UInt_t k = 0; k < ybin1-1; k++){
      cc->cd(k+1);
      mv1[k]->Draw("ALP");
      if( k == ybin1-2 )
	lg1[k]->Draw();
    }

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1400,500);
    cc->Divide(ybin2-1,1);
    for(UInt_t k = 0; k < ybin2-1; k++){
      cc->cd(k+1);
      mv2[k]->Draw("ALP");
      if( k == ybin2-2 ) lg2[k]->Draw();
    }
  }


  if( bplot[0] && bplot[6] ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    mv1[2]->Draw("ALP");
    lg1[2]->Draw();
    
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    mv1[3]->Draw("ALP");
    lg1[3]->Draw();

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    mv1[4]->Draw("ALP");
    lg1[4]->Draw();
  }

  //--------------------
  //--- 
  if( bplot[0] || bplot[1] ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    mrv1->Draw("ALP");
    lgr1->Draw();

    auto Ymin = mrv1->GetYaxis()->GetXmin();
    auto Ymax = mrv1->GetYaxis()->GetXmax();
    auto Xmin = mrv1->GetXaxis()->GetXmin();
    auto Xmax = mrv1->GetXaxis()->GetXmax();

    auto aLineX1 = new TLine(Xmin, 0., Xmax, 0.);
    aLineX1->SetLineColor(1);
    aLineX1->SetLineStyle(3);
    aLineX1->Draw();

    auto aLineY1 = new TLine(0., Ymin, 0., Ymax);
    aLineY1->SetLineColor(1);
    aLineY1->SetLineStyle(3);
    aLineY1->Draw();
  
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    mrv2->Draw("ALP");
    lgr2->Draw();

    Ymin = mrv2->GetYaxis()->GetXmin();
    Ymax = mrv2->GetYaxis()->GetXmax();
    Xmin = mrv2->GetXaxis()->GetXmin();
    Xmax = mrv2->GetXaxis()->GetXmax();

    auto aLineX2 = new TLine(0., Ymin, 0., Ymax);
    aLineX2->SetLineColor(1);
    aLineX2->SetLineStyle(3);
    aLineX2->Draw();
  }



  if( bplot[0] && bplot[4] ) {  // individual pt plots  

    UInt_t ichoise[]={1,3,6};
    for(UInt_t i : ichoise) {
      cc = new TCanvas(Form("cv%d",i),Form("cv%d",i),500,550);
      cc->SetRightMargin(0.02);
      cc->SetLeftMargin(0.2);
      cc->SetTopMargin(0.08);
      cc->SetBottomMargin(0.15);

      mv1[i]->GetXaxis()->SetNdivisions(505);
      mv1[i]->GetXaxis()->SetLabelSize(0.04);
      mv1[i]->GetYaxis()->SetLabelSize(0.04);
      mv1[i]->GetXaxis()->SetTitleSize(0.05);
      mv1[i]->GetYaxis()->SetTitleSize(0.04);
      mv1[i]->GetYaxis()->SetTitleOffset(2.5);

      mv1[i]->GetXaxis()->SetRangeUser(0.,550.);
      mv1[i]->Draw("ALP");

      if( ngr > 1 )
	lg1[i]->Draw();

    }

    UInt_t ichoise2[] = {2};
    for(UInt_t i : ichoise2) {
      cc = new TCanvas(Form("cv%d",i+10),Form("cv%d",i+10),500,550);
      cc->SetRightMargin(0.02);
      cc->SetLeftMargin(0.2);
      cc->SetTopMargin(0.08);
      cc->SetBottomMargin(0.15);

      mv1[i]->GetXaxis()->SetNdivisions(505);
      mv2[i]->GetXaxis()->SetLabelSize(0.04);
      mv2[i]->GetYaxis()->SetLabelSize(0.04);
      mv2[i]->GetXaxis()->SetTitleSize(0.05);
      mv2[i]->GetYaxis()->SetTitleSize(0.04);
      mv2[i]->GetYaxis()->SetTitleOffset(2.5);

      mv2[i]->GetXaxis()->SetRangeUser(0.,550.);
      mv2[i]->Draw("ALP");
      if( ngr > 1 )
	lg2[i]->Draw();

    }
  }

  if( bplot[0] && bplot[7] ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    g_v2max->SetMarkerStyle(20);
    g_v2max->SetMarkerColor(2);
    g_v2max->Draw("AP");
    g_v2max->Print();

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    g_v1slp->SetMarkerStyle(20);
    g_v1slp->SetMarkerColor(2);
    g_v1slp->Draw("AP");
    g_v1slp->Print();

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    auto msysD = new TMultiGraph("msysD","; (N-P)/A; F[MeV/c]");
    g_mpxslp->SetMarkerStyle(20);
    g_mpxslp->SetMarkerColor(2);
    msysD->Add(g_mpxslp, "LP");
    lgr0->AddEntry(g_mpxslp,"DATA x 1.4");

    if( bplot[1] ) {
      g_pBUUFg5 ->SetMarkerStyle(21);
      g_pBUUFg5 ->SetMarkerColor(6);
      g_pBUUFg5 ->SetLineColor(6);
      g_pBUUFg17->SetMarkerStyle(21);
      g_pBUUFg17->SetMarkerColor(8);
      g_pBUUFg17->SetLineColor(8);
	
      msysD->Add(g_pBUUFg5,"LP");
      msysD->Add(g_pBUUFg17,"LP");
      lgr0->AddEntry(g_pBUUFg5,"pBUU b4 g0.5");
      lgr0->AddEntry(g_pBUUFg17,"pBUU b4 g1.75");
    }
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

      if( bplot[1] ) {
	g_pBUUv2g5 ->SetMarkerStyle(21);
	g_pBUUv2g5 ->SetMarkerColor(6);
	g_pBUUv2g5 ->SetLineColor(6);
	g_pBUUv2g17->SetMarkerStyle(21);
	g_pBUUv2g17->SetMarkerColor(8);
	g_pBUUv2g17->SetLineColor(8);
	
	mv2sysD->Add(g_pBUUv2g5,"LP");
	mv2sysD->Add(g_pBUUv2g17,"LP");
	lgr10->AddEntry(g_pBUUv2g5,"pBUU b4 g0.5");
	lgr10->AddEntry(g_pBUUv2g17,"pBUU b4 g1.75");
      }

      ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
      mv2sysD->Draw("AP");
      lgr10->Draw();
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

