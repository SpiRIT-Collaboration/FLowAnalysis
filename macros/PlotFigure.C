#include "DoFlow.h"
#include "SetStyle.C"
#include "FlowFunction.C"
#include "SetColor.C"

void CanvasPartitionY(TCanvas *C, const Int_t Ny = 3,
		    Float_t lMargin = 0.15, Float_t rMargin = 0.05,
		    Float_t bMargin = 0.15, Float_t tMargin = 0.05);

void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
		    Float_t lMargin = 0.15, Float_t rMargin = 0.05,
		    Float_t bMargin = 0.15, Float_t tMargin = 0.05);

void CanvasPartitionTwoColumn(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
			      Float_t lMargin = 0.15, Float_t rMargin = 0.05,
			      Float_t bMargin = 0.15, Float_t tMargin = 0.05,
			      Float_t midMargin=0.05);

struct gplot{
  TString Version;
  TString fileHeader;
  TString comment;
};

TString  bName[]   = {"132Sn_","108Sn_","124Sn_","112Sn_","pp",    "100Sn_"};
Double_t sysD[]    = {0.22,    0.09,      0.15,   0.15   ,   0,      0.22};
Double_t sysA[]    = {256.,    220.,      236.,   236.   ,   2,    256.};
Double_t sysBN[]   = {82./50., 58./50.,   74./50.,62./50.,   0,     50./50.};
Double_t sysN[]    = {156.,    110.,      136.,   136.   ,   0,    156.};
Double_t sysNA[]   = {156./100.,110./100.,136./100.,136./100., 0,  156./100.};
Size_t   imsz[]    = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
Color_t  pcolor[]  = {2, 3, 4, 5, 6}; // p,d,t,3He

TString FOPI_data_sys[] = {"Au","Ru","Ca"};



std::vector<UInt_t> sqv1sel = {10, 8, 2, 1};
std::vector<UInt_t> sqv2sel = {4,3,2,1};
std::vector<UInt_t> sqPart  = {0,1,2,3,4};
const UInt_t npart = 5;

gplot gnames[] = {
  {".v52.15.39","finYPt_", "45>|#phi|"}
};

Double_t FOPI_AuAu_v10[][2]={
  {1.0, 0.384},  //p
  {2.0, 0.641},  //d
  {3.0, 0.800}}; //A=3
Double_t FOPI_AuAu_v20[][2]={
  {1.0, -0.048},
  {2.0, -0.105},
  {3.0, -0.170}};


TCanvas *ccv; UInt_t iccv = 0;
TString xlabel;
TLatex  plabel;


Double_t *syslabel;
TGraphErrors* g_v1slp[npart];
TGraphErrors* g_v2max[npart];
TGraphErrors* g_v2n[npart];
TGraphErrors* g_v1pslp[npart];
TMultiGraph* mrv1[npart];
TMultiGraph* mrv2[npart];
TLegend*     lrv1[npart];
TLegend*     lrv2[npart];
TMultiGraph* g_v1sysD;
TLegend*     l_v1sysD;
TMultiGraph* g_v2sysD;
TLegend*     l_v2sysD;
TMultiGraph* g_v1psysD;
TLegend*     l_v1psysD;
TMultiGraph* g_v2nsysD;
TLegend*     l_v2nsysD;
TGraphErrors *v1Ut;;
TGraphErrors *v2Ut;



FairLogger *logger = FairLogger::GetLogger();
void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer);
void Draw_v1Ratio();
void Draw_v1v2SystemD();
void Draw_v1v2SystemD_fit();
void Draw_Indiv_v1SystemD();
void Draw_Indiv_v2SystemD();
void Draw_Indiv_v2SystemD();
void Draw_Indiv_v1SystemD_Three();
void Draw_v2SystemD();
void Draw_v2nSystemD();
void Draw_v2SystemD_one();
void Draw_v1Ut_SystemD(UInt_t ipart);
void Draw_v2Ut_SystemD(UInt_t ipart);
void Draw_Ut_SystemD();
void Draw_Ut_Comparison(UInt_t isys);
void Draw_Ut_Comparison_FOPI(UInt_t isys);
void Draw_Ut_Particle();
void Draw_v_y(UInt_t vn);
void LoadData();
TGraphErrors* LoadAMD(UInt_t isys, UInt_t ipart, UInt_t iflow);
TGraphErrors* LoadAMD(UInt_t isys, UInt_t ipart, TString grname);
TGraphErrors* LoadFOPI(UInt_t ipart, TString fdir, TString sysname );
TGraphErrors* LoadTextGraph(TString fname);

void PlotFigure()
{

  gStyle->SetOptStat(0);

  SetStyle();
  SetColor();

  //  LoadData();


  // Draw function

  if( 0 ) { //v11 v20 system dependence
    Draw_Indiv_v1SystemD();
    Draw_Indiv_v2SystemD();
  }

  if( 0 )
    Draw_v2SystemD_one();
  
  if( 0 )
    Draw_v1Ratio();

  if( 0 )   // fitting results for v1 and v2 vs y
    Draw_v1v2SystemD_fit();

  if( 0 )
    Draw_v1v2SystemD();

  if( 0 )
    Draw_v2SystemD();
  
  //edit
  if( 0 ) {
    Draw_v1Ut_SystemD(1);
    Draw_v2Ut_SystemD(1);
  }

  if( 0 ) //OK
    Draw_Ut_SystemD(); // Final

  if( 1 ) 
    Draw_Ut_Comparison(3); // Comparison with AMD (0:system) 

  if( 0 ) 
    Draw_Ut_Comparison_FOPI(0); // Comparison with FOPI

  if( 0 )
    Draw_Ut_Particle();

  if( 0 )
    Draw_v_y(2);

}
//------------------------------------------------
//------------------------------------------------
//------------------------------------------------
TGraphErrors* LoadFOPI(UInt_t ipart=0, TString fdir="PRC89/Fig8_v1Ut_0.25", TString sysname="" )
{
  TString FOPI_dir =  "data/FOPI/";

  FOPI_dir += fdir + "_" + fpid[ipart] + "_" + sysname + ".txt";

  TGraphErrors* grph = LoadTextGraph(FOPI_dir); 
  grph -> SetName("g_utv1_FOPI_" +  lpid[ipart] + sysname );

  if( grph->GetN() > 0 ) 
    LOG(INFO) << grph->GetName() << " is registered." << FairLogger::endl;
  else {
    LOG(INFO) << grph->GetName() << " is not found." << FairLogger::endl;
    return NULL;
  }

  grph->Print();

  return grph;
}

TGraphErrors* LoadAMD(UInt_t isys, UInt_t ipart, TString grname="v1_ut" )
{
  TGraphErrors *grv = NULL;

  // TString filename[] = {"Sn132/SLy4_L108/flow_proton.root",
  // 			"Sn132/SLy4_L108/flow_deuteron.root",
  // 			"Sn132/SLy4_L108/flxoxw_triton.root"};

  TString filename[] = {"Sn132/SLy4/flow_proton.root",
			"Sn132/SLy4/flow_deuteron.root",
			"Sn132/SLy4/flow_triton.root",
			"Sn132/SLy4/flow_alpha.root"};

  auto ifile = new TFile("../../TransportModel/AMD/"+filename[ipart],"READ");

  if( !ifile ) return NULL;

  grv = (TGraphErrors*)ifile->Get(grname);
  if( grv == NULL ) return NULL;

  TString gname = grname + "_" + lpid[ipart];
  grv -> SetName(gname);

  delete ifile;

  return grv;
}


TGraphErrors* LoadData(UInt_t isys, UInt_t ipart, TString gname)
{
  TGraphErrors* grv = NULL;

  TFile *fOpen;
  TString fname = gnames[0].fileHeader + bName[isys] + fpid[ipart] + gnames[0].Version + ".root";

  if( !gSystem->FindFile("data", fname) ) {
    LOG(ERROR) << fname << " is not found " << FairLogger::endl;
    return NULL;
  }
  else {
    fOpen = TFile::Open( fname );
    LOG(INFO) << fname << " is opened. " << FairLogger::endl;
  }

  grv =  (TGraphErrors*)fOpen->Get(gname);
  if( grv == NULL ) return NULL;


  gname += "_" + rsys[isys] + "_" + lpid[ipart];
  grv -> SetName(gname);

  fOpen->Close();

  return grv;
}


TGraphErrors* LoadTextGraph(TString fname)
{
  std::fstream fread;
  fread.open(fname, std::fstream::in);

  LOG(INFO) << fname << " is opened. " << FairLogger::endl;
  auto vut = new TGraphErrors();
  Double_t x, y;
  TString sget;
  UInt_t  in = 0;
  while( !fread.eof() ) {
    fread >> sget;
    x = (Double_t)atof(sget);
    fread >> sget;
    y = (Double_t)atof(sget);

    if( !std::isnan(x) ) {
      vut -> SetPoint(in, x, y);
      in++;
    }
  }
  
  if( in == 0 ) return NULL;

  vut -> RemovePoint(in-1);
  return vut;
}



void LoadData()
{
  Color_t fcolor = 2;
  UInt_t ix = 0;
  switch(ix) {
  case 0:
    xlabel = "(N-P)/A";
    syslabel = sysD;
    break;
  case 1:
    xlabel = "N";
    syslabel = sysN;
    break;
  case 2:
    xlabel = "N/Z(Beam)";
    syslabel = sysBN;
    break;
  case 3:
    xlabel = "N/Z";
    syslabel = sysNA;
    break;
  }

  TGraphErrors *yv1;
  TGraphErrors *yv2;

  Double_t v1_slp [4][npart];
  Double_t v1_slpe[4][npart];
  Double_t v2_max [4][npart];
  Double_t v2_maxe[4][npart];
  Double_t  v20[4][npart];
  Double_t v20e[4][npart];
  Double_t  v21[4][npart];
  Double_t v21e[4][npart];
  Double_t  v2n[4][npart];
  Double_t v2ne[4][npart];
  

  for(UInt_t i = 0; i < npart; i++ ){
    mrv1[i] = new TMultiGraph(Form("mrv1_%d",i)  ,";y/y_{nn}-1; v1");
    mrv2[i] = new TMultiGraph(Form("mrv2_%d",i)  ,";y/y_{nn}-1; v2");
    lrv1[i] = new TLegend(0.6 , 0.15, 0.90, 0.50, lsys[i]);
    lrv2[i] = new TLegend(0.75 , 0.15, 0.90, 0.50, lsys[i]); 
  }

  auto g_v1off = new TGraphErrors();
  std::vector< TString > g_v1label;
  UInt_t iv1off = 0;

  TFile *fOpen;

  for( UInt_t is = 0; is < 4; is++ ) {

    if( is == 2 ) continue;

    for( UInt_t ip = 0; ip < npart; ip++ ) {

      TString fname = gnames[0].fileHeader + bName[is] + fpid[ip] + gnames[0].Version + ".root";

      if( !gSystem->FindFile("data", fname) ) {
	LOG(ERROR) << fname << " is not found " << FairLogger::endl;
	continue;
      }
      else {
	fOpen = TFile::Open( fname );
	LOG(INFO) << fname << " is opened. " << FairLogger::endl;
      }


      //yv1 = (TGraphErrors*)fOpen->Get("gy_v1");
      //yv2 = (TGraphErrors*)fOpen->Get("gy_v2");
      yv1 = (TGraphErrors*)fOpen->Get("gu_v1");
      yv2 = (TGraphErrors*)fOpen->Get("gu_v2");

      if( yv1 && yv2 ) {
	fcolor = pcolor[ip];
	yv1->SetMarkerColor(fcolor);
	yv1->SetMarkerStyle(imark[0]);
	yv1->SetLineColor(fcolor);

	yv2->SetMarkerColor(fcolor);
        yv2->SetMarkerStyle(imark[0]);
        yv2->SetMarkerSize(imsz[0]);
        yv2->SetLineColor(fcolor);

	for( Int_t j = (Int_t)yv2->GetN()-1; j >= 0; j-- ) {
	  Double_t xx, yy;
	  yv2->GetPoint(j, xx, yy);
	  if( xx > 1. )
	    yv2->RemovePoint(j);
	}
      }	


      mrv1[is] -> Add(yv1,"p");
      lrv1[is] -> AddEntry(yv1,  lpid[ip] ,"lp");
      // mrv1[ip] -> Add(yv1,"p");
      // lrv1[ip] -> AddEntry(yv1,  lsys[is] ,"lp");

      fv1fit->SetLineColor(fcolor);
      fv1fit->SetParameter(1,0.2);
      yv1->Fit("fv1fit","","",-0.6,0.6); //"Q0","");    

      v1_slp [is][ip]   = fv1fit->GetParameter(1);
      v1_slpe[is][ip]   = fv1fit->GetParError(1);

      if( ip != 3 ) {
	g_v1label.push_back( bName[is] + lpid[ip] ) ;
	g_v1off -> SetPoint(iv1off, iv1off, fv1fit->GetParameter(3) );
	g_v1off -> SetPointError(iv1off, 0, fv1fit->GetParError(3) );
	iv1off++;
      }

      if( yv2 != NULL ) {
        for( Int_t iip = (Int_t)yv2->GetN()-1; iip >= 0; iip-- ){

          Double_t xpnt, ypnt;
          ypnt = yv2->GetErrorY(iip);
          if( ypnt > 0.05 )
            yv2->RemovePoint(iip);
        }

	Double_t v2x, v2y, v2ye;
        GetMinimumv2(yv2, v2_max[is][ip], v2_maxe[is][ip]);	
	yv2   ->SetLineColor(fcolor);
	fv2fit->SetLineColor(fcolor);

	fv2fit->SetParameter(0, -0.01);
	fv2fit->SetParameter(1,  0.15);
	fv2fit->SetParameter(2,  0.);

	// if( ip == 3 )
	//   yv2->Fit("fv2fit","","",-0.48,1.);
	// else
	  yv2->Fit("fv2fit","","",-0.4,0.45);

	mrv2[is]->Add(yv2,"p");
	lrv2[is]->AddEntry(yv2, lpid[ip], "p");
	// mrv2[ip]->Add(yv2,"p");
	// lrv2[ip]->AddEntry(yv2, lsys[is], "p");

	v20[is][ip]  = fv2fit->GetParameter(0);
	v20e[is][ip] = fv2fit->GetParError(0);

	v21[is][ip]  = fv2fit->GetParameter(1);
	v21e[is][ip] = fv2fit->GetParError(1);
	
	v2n[is][ip]  = abs(v21[is][ip]) + abs(v20[is][ip]);
	v2ne[is][ip] = 0.; sqrt( pow(v21[is][ip],2) + pow(v20[is][ip],2) );

      }
     
      // Getting Ut plots
      TString gname = "g_utv1";
      v1Ut =  (TGraphErrors*)fOpen->Get(gname);
      cout << gname << endl;
      gname += "_" + rsys[is] + "_" + lpid[ip];
      v1Ut -> SetName(gname);
      cout << gname << endl;
      gROOT->Add(v1Ut);

      gname = "g_utv2";
      v2Ut =  (TGraphErrors*)fOpen->Get(gname);
      gname += "_" + rsys[is] + "_" + lpid[ip];
      v2Ut -> SetName(gname);
      gROOT->Add(v2Ut);


      for(auto k: sqv1sel) {
	gname = (TString)Form("gUt_v1%d",k);
	v1Ut = (TGraphErrors*)fOpen->Get(gname);
	gname += "_" + rsys[is] + "_" + lpid[ip];
	v1Ut -> SetName(gname);
	gROOT->Add(v1Ut);
      }

      for(auto k: sqv2sel) {
	gname = (TString)Form("gUt_v2%d",k);
	v2Ut = (TGraphErrors*)fOpen->Get(gname);
	gname += "_" + rsys[is] + "_" + lpid[ip];
	v2Ut -> SetName(gname);
	gROOT-> Add(v2Ut);
      }
 
      fOpen->Close();
    }
    fcolor++;
    if( fcolor == 10 ) fcolor++;
    
  }

  //----------
  for( UInt_t ip = 0; ip < npart; ip++ ){
    g_v1slp[ip]  = new TGraphErrors();
    g_v1slp[ip]  -> SetName(Form("g_v1slp_%d",ip));
    g_v1slp[ip]  -> SetTitle(";"+xlabel+";");
    //    if( ip == 0 )
      g_v1slp[ip]  -> SetTitle(";"+xlabel+";v_{11}");

    g_v2max[ip]  = new TGraphErrors();
    g_v2max[ip]  -> SetName(Form("g_v2max_%d",ip));
    g_v2max[ip]  -> SetTitle(";"+xlabel+"; -v_{20}");
    g_v2n[ip]    = new TGraphErrors();
    g_v2n[ip]    -> SetName(Form("g_v2n_%d",ip));
    g_v2n[ip]    -> SetTitle(";"+xlabel+"; -v_{2n}");

    g_v1pslp[ip] = new TGraphErrors();
  }

  
  g_v1sysD = new TMultiGraph("g_v1slpD", ";"+xlabel+"; v_{11}");
  l_v1sysD = new TLegend(0.2, 0.7, 0.5, 0.9, "");
  g_v2sysD = new TMultiGraph("g_v2sysD", ";"+xlabel+"; -v20");
  l_v2sysD = new TLegend(0.35, 0.13, 0.7, 0.33, "");

  g_v1psysD = new TMultiGraph("g_v1pslpD", ";"+xlabel+"; v_{11}/ v_{11}(proton)");
  l_v1psysD = new TLegend(0.2, 0.7, 0.5, 0.9, "");

  g_v2nsysD = new TMultiGraph("g_v2bsysD", ";"+xlabel+"; -v2n");
  l_v2nsysD = new TLegend(0.35, 0.13, 0.7, 0.33, "");

  for( UInt_t ip = 0; ip < npart; ip++ ) {

    UInt_t iss = 0;
    for( UInt_t is = 0; is < 4; is++ ) {
      
      if( is == 2 ) continue;

      if( kTRUE ) {
	LOG(INFO) << " v1slp " << v1_slp[is][ip]
		  << " v2max " << v2_max[is][ip]
		  << FairLogger::endl;

	g_v1slp[ip] -> SetPoint     ( iss,  *(syslabel+is), v1_slp[is][ip] ); 
	g_v1slp[ip] -> SetPointError( iss, 0,   v1_slpe[is][ip] );

	// g_v2max[ip] -> SetPoint     ( iss,  *(syslabel+is), -v2_max[is][ip] ); 
	// g_v2max[ip] -> SetPointError( iss,  0., v2_maxe[is][ip] ); 
	g_v2max[ip] -> SetPoint     ( iss,  *(syslabel+is), -v20[is][ip] ); 
	g_v2max[ip] -> SetPointError( iss,  0., v20e[is][ip] ); 

	g_v2n[ip]   -> SetPoint     ( iss,  *(syslabel+is), v2n[is][ip] );
	g_v2n[ip]   -> SetPointError( iss,  0., v2ne[is][ip] );
	

	if( ip > 0 ) {
	  auto ratio = v1_slp[is][ip]/v1_slp[is][0];
	  auto ratie = ratio * sqrt( pow(v1_slpe[is][ip]/v1_slp[is][ip],2) + pow(v1_slpe[is][0]/v1_slp[is][0],2) );
	  g_v1pslp[ip] -> SetPoint     ( iss, *(syslabel+is), ratio );
	  g_v1pslp[ip] -> SetPointError( iss,        0., ratie );
	}

	iss++;
      }
    }

    g_v1slp[ip]->SetMarkerStyle(20);
    g_v1slp[ip]->SetMarkerSize(1.5);
    g_v1slp[ip]->SetMarkerColor(icol[ip]);
    g_v1slp[ip]->SetLineColor(icol[ip]);

    g_v2max[ip]->SetMarkerStyle(20);
    g_v2max[ip]->SetMarkerSize(1.5);
    g_v2max[ip]->SetMarkerColor(icol[ip]);
    g_v2max[ip]->SetLineColor(icol[ip]);

    g_v2n[ip]->SetMarkerStyle(20);
    g_v2n[ip]->SetMarkerSize(1.5);
    g_v2n[ip]->SetMarkerColor(icol[ip]);
    g_v2n[ip]->SetLineColor(icol[ip]);

    
    g_v1sysD -> Add( g_v1slp[ip], "p" );
    l_v1sysD -> AddEntry( g_v1slp[ip], lpid[ip], "P" );

    g_v2sysD -> Add( g_v2max[ip],"p");
    l_v2sysD -> AddEntry( g_v2max[ip], lpid[ip] ,"P");

    g_v2nsysD -> Add( g_v2n[ip],"p");
    l_v2nsysD -> AddEntry( g_v2n[ip], lpid[ip] ,"P");

    g_v1pslp[ip]->SetMarkerStyle(20);
    g_v1pslp[ip]->SetMarkerSize(1.5);
    g_v1pslp[ip]->SetMarkerColor(icol[ip]);
    g_v1pslp[ip]->SetLineColor(icol[ip]);

    
    if( ip > 1 ) {
      g_v1psysD -> Add( g_v1pslp[ip], "p" );
      l_v1psysD -> AddEntry( g_v1slp[ip], lpid[ip], "P" );
    }

  }

}



void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer)
{
  Double_t x, y;
  min = 9.;
  UInt_t mid = 0;
  for(UInt_t i = 0; i < gr->GetN(); i++) {
    gr->GetPoint(i, x, y);
    if( abs(y) > 0.2 ) {
      gr->RemovePoint(i); 
      i = 0;
      continue;
    }

    if( y < min && y < 0) {
      min = y;
      mid = i;
    }
  }


  mer = gr->GetErrorY(mid);
}

//----------------------------
//------- Draw ---------------
//----------------------------
void Draw_Ut_Comparison(UInt_t isys)
{
  LOG(INFO) << "Draw_Ut_Comparison " << FairLogger::endl;


  auto mv1 = new TMultiGraph(Form("mv1_%d",isys),";U_{t0}; v1");
  auto mv2 = new TMultiGraph(Form("mv2_%d",isys),";U_{t0}; v2");
  auto lg1 = new TLegend(0.64,0.15,0.8,0.5,"");
  auto lg2 = new TLegend(0.18 ,0.15,0.4,0.5,"");


  // --------- DATA _______
  if( 1 ) {
    std::vector<UInt_t> sqPart  = {0,1,2,3,4};
    for(auto ipart: sqPart ) {
      TString gname = "g_utv1_1" ;
      v1Ut = (TGraphErrors*)LoadData(isys, ipart, gname);
      v1Ut -> SetName(gname);

      if( v1Ut != NULL ) {
	LOG(INFO) << " v1ut " << gname << " is registered." << FairLogger::endl;
	v1Ut -> SetMarkerStyle(20);
	v1Ut -> SetMarkerSize(1.);
	v1Ut -> SetLineColor(icol[ipart]);
	v1Ut -> SetMarkerColor(icol[ipart]);

	mv1 -> Add( v1Ut, "pl");
	lg1 -> AddEntry(v1Ut,lpid[ipart]);
      }
      else
	LOG(ERROR) << " v1ut " << gname << " is not found." << FairLogger::endl;

      gname = "g_utv2";
      v2Ut = (TGraphErrors*)LoadData(isys, ipart, gname);
      v2Ut -> SetName(gname);

      if( v2Ut != NULL ) {
	LOG(INFO) << " v2ut " << gname << " is registered." << FairLogger::endl;

	v2Ut -> SetMarkerStyle(20);
	v2Ut -> SetMarkerSize(1.);
	v2Ut -> SetLineColor(icol[ipart]);
	v2Ut -> SetMarkerColor(icol[ipart]);

	mv2 -> Add( v2Ut, "pl" );

	lg2 -> AddEntry(v2Ut,lpid[ipart]);

      }
    }
  }
  //------ FOPI ____
  if( 1 ) {
    TString FOPI_v1Ut[][3] = { {"PRC89/Fig8_v1Ut_" , "0.25" ,"Au"}, //[0]
			       {"PRC89/Fig8_v1Ut_" , "0.4"  ,"Au"}, //[1]
			       {"NPA876/Fig8_v1Ut_", "0.25","Au"},  //[2]
			       {"NPA876/Fig9_v1Ut_", "0.4","Au"},   //[3]
			       {"NPA876/Fig11_v1Ut_","0.4","Au"},   //[4]
			       {"NPA876/Fig11_v1Ut_","0.4","Ru"},   //[5]
			       {"NPA876/Fig11_v1Ut_","0.4","Ca"}};  //[6]
			    
    TString FOPI_v2Ut[][3] = { {"PRC89/Fig8_v2Ut_","0.25","Au"},
			       {"PRC89/Fig8_v2Ut_","0.4" ,"Au"}};


    //			       {"NPA876/Fig7_v1y_", "0.4" ,"Au"}, //[2]

    // system comparison
    UInt_t ipart = 1;
    std::vector< UInt_t > sqFOPIfig1  = {3};
    std::vector< UInt_t > sqFOPIpart = {0,1,2,4};
    for(auto ifig: sqFOPIfig1)for(auto ipart: sqFOPIpart) {
	v1Ut = (TGraphErrors*)LoadFOPI(ipart, FOPI_v1Ut[ifig][0]+FOPI_v1Ut[ifig][1], FOPI_v1Ut[ifig][2]);
	if( v1Ut != NULL ) {
	  v1Ut -> SetMarkerStyle(25+ifig);
	  v1Ut -> SetMarkerSize(1.);
	  v1Ut -> SetLineColor(icol[ipart]);
	  v1Ut -> SetMarkerColor(icol[ipart]);

	  mv1 -> Add( v1Ut, "p");
	  TString partname = lpid[ipart];
	  if( ipart == 2) partname = "A3";
	  lg1 -> AddEntry(v1Ut,partname+"(FOPI)"+FOPI_v1Ut[ifig][2]+FOPI_v1Ut[ifig][1]);
	}
      }
  
    std::vector< UInt_t > sqFOPIfig2  = {1};
    for(auto ifig: sqFOPIfig2)for(auto ipart: sqFOPIpart) {
	v2Ut = (TGraphErrors*)LoadFOPI(ipart, FOPI_v2Ut[ifig][0]+FOPI_v2Ut[ifig][1],FOPI_v2Ut[ifig][2]);
	if( v2Ut != NULL ) {
	  v2Ut -> SetMarkerStyle(25+ifig);
	  v2Ut -> SetMarkerSize(1.);
	  v2Ut -> SetLineColor(icol[ipart]);
	  v2Ut -> SetMarkerColor(icol[ipart]);

	  mv2 -> Add( v2Ut, "p" );
	  TString partname = lpid[ipart];
	  if( ipart == 2) partname = "A3";
	  lg2 -> AddEntry(v2Ut,partname+"(FOPI)"+FOPI_v2Ut[ifig][2]+FOPI_v2Ut[ifig][1]);
	}
      }
  }
  

  //-------- AMD ___
  if( 0 ) {
    for(auto ipart: sqPart ) {
      v1Ut = (TGraphErrors*)LoadAMD(isys, ipart, "v1_ut");
      if( v1Ut != NULL ) {
	v1Ut -> SetMarkerStyle(25);
	v1Ut -> SetMarkerSize(1.);
	v1Ut -> SetLineColor(icol[ipart]);
	v1Ut -> SetMarkerColor(icol[ipart]);


	mv1 -> Add( v1Ut, "p");
	lg1 -> AddEntry(v1Ut,lpid[ipart]+"(AMD)");
      }

      v2Ut = (TGraphErrors*)LoadAMD(isys, ipart, "v2_ut");    
      if( v2Ut != NULL ) {
	v2Ut -> SetMarkerStyle(25);
	v2Ut -> SetMarkerSize(1.);
	v2Ut -> SetLineColor(icol[ipart]);
	v2Ut -> SetMarkerColor(icol[ipart]);


	mv2 -> Add( v2Ut, "p" );
	lg2 -> AddEntry(v2Ut,lpid[ipart]+"(AMD)");
      }
    }
  }
  //----

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 800, 800); iccv++;
  ccv -> Divide(1, 2);

  ccv -> cd(1);

  mv1->GetXaxis()->SetRangeUser( 0.,  2.2);
  mv1->GetYaxis()->SetRangeUser(-0.05,0.62);
  mv1->GetXaxis()->SetNdivisions(505);
  mv1 -> Draw("AP");
  lg1 -> Draw();
  TLatex *v1label = new TLatex(0.3,0.62," 0.25 < b_{0} < 0.45 0.4 < y_{norm} < 0.8");
  v1label -> Draw();
  auto labelv1 = new TLatex(1.3,  0.5, fsys[isys]);
  labelv1 -> Draw();


  ccv -> cd(2);
  mv2->GetXaxis()->SetRangeUser( 0.,  2.2);
  mv2->GetYaxis()->SetRangeUser(-0.32,0.01);
  mv2->GetXaxis()->SetNdivisions(505);
  mv2 -> Draw("AP");
  lg2 -> Draw();
  TLatex *v2label = new TLatex(0.3,0.02," 0.25 < b_{0} < 0.45  -0.4 < y_{norm} < 0.4");
  v2label -> Draw();
  auto labelv2 = new TLatex(1.3, -0.05, fsys[isys]);
  labelv2 -> Draw();
}
void Draw_Ut_Comparison_FOPI(UInt_t isys)
{
  LOG(INFO) << "Draw_Ut_Comparison_FOPI " << FairLogger::endl;

  auto labelv1 = new TLatex(1.3,  0.45, fsys[isys]);
  auto labelv2 = new TLatex(1.3, -0.008, fsys[isys]);

  TMultiGraph *mv1[2];
  mv1[0] = new TMultiGraph(Form("mv1_1H_%d",isys),";U_{t0}; v1");
  mv1[1] = new TMultiGraph(Form("mv1_2H_%d",isys),";U_{t0}; v1");

  TLegend* lg1[2];
  lg1[0] = new TLegend(0.64,0.2,0.8,0.5,"");
  lg1[1] = new TLegend(0.64,0.2,0.8,0.5,"");


  std::vector<UInt_t> sqPart  = {0,1};
  for(auto ipart: sqPart ) {
    TString gname = "g_utv1_1" ;
    v1Ut = (TGraphErrors*)LoadData(isys, ipart, gname);
    v1Ut -> SetName(gname);

    if( v1Ut != NULL ) {
      LOG(INFO) << " v1ut " << gname << " is registered." << FairLogger::endl;
      v1Ut -> SetMarkerStyle(20);
      v1Ut -> SetMarkerSize(1.);
      v1Ut -> SetLineColor(icol[ipart]);
      v1Ut -> SetMarkerColor(icol[ipart]);

      mv1[ipart] -> Add( v1Ut, "pl");
      lg1[ipart] -> AddEntry(v1Ut,lpid[ipart]);
    }
    else
      LOG(ERROR) << " v1ut " << gname << " is not found." << FairLogger::endl;

  
    // system comparison
    std::vector< UInt_t> sqFOPI = {0,1,2};
    for(auto isys: sqFOPI ) {
      v1Ut = (TGraphErrors*)LoadFOPI(ipart,"NPA876/Fig7_v1Ut_0.4", FOPI_data_sys[isys]);
      if( v1Ut != NULL ) {
	v1Ut -> SetMarkerStyle(25);
	v1Ut -> SetMarkerSize(1.);
	v1Ut -> SetLineColor(icol[isys]);
	v1Ut -> SetMarkerColor(icol[isys]);

	mv1[ipart] -> Add( v1Ut, "p");
	lg1[ipart] -> AddEntry(v1Ut,lpid[ipart]+"(FOPI)"+FOPI_data_sys[isys]);
      }
      else
	LOG(ERROR) << " v1ut " << gname << " is not found." << FairLogger::endl;

    }
  }


  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 800, 800); iccv++;
  ccv -> Divide(1,2);

  for(auto ipart: sqPart ) {
    ccv -> cd(ipart+1);
    mv1[ipart]->GetXaxis()->SetRangeUser( 0.,  2.2);
    //    mv1[ipart]->GetYaxis()->SetRangeUser(-0.05,0.62);
    mv1[ipart]->GetXaxis()->SetNdivisions(505);
    mv1[ipart] -> Draw("AP");
    lg1[ipart] -> Draw();
    TLatex *v1label = new TLatex(0.3,0.62," 0.25 < b_{0} < 0.45 0.4 < y_{norm} < 0.8");
    v1label -> Draw();
    labelv1 -> Draw();
  }
}

void Draw_Ut_Particle()
{
  std::vector<UInt_t> sqSys   = {1,3,0};
  std::vector<UInt_t> sqPart  = {0,1,2,3,4};

  LOG(INFO) << "Draw_Ut_Particle " << FairLogger::endl;
  
  TLegend *lg1 = new TLegend(0.75,0.4,0.9,0.65,"");
  TLegend *lg2 = new TLegend(0.25,0.4,0.4,0.65,"");
  TLatex *labelv1[npart];
  TLatex *labelv2[npart];
       
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 1000, 800); iccv++;
  gStyle->SetOptTitle(0);

  const Int_t Nx = 2;
  const Int_t Ny = 5;

  Float_t lMargin = 0.08;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.10;
  Float_t mMargin = 0.08;

  CanvasPartitionTwoColumn(ccv,Nx,Ny,lMargin,rMargin,bMargin,tMargin,mMargin);

  TPad *pad[Nx][Ny];

  for( auto i : ROOT::TSeqI(Nx) )for( auto j : ROOT::TSeqI(Ny) ) {
      ccv->cd(0);

      TString pname = Form("pad_%i_%i",i,j);
      pad[i][j] = (TPad*)gROOT->FindObject(pname);
      cout << " pad " << pname << endl;
      if( pad[i][j] == NULL ) {
	cout << " pad is not found " << pname << endl;
	continue;
      }
      pad[i][j] -> Draw();
      pad[i][j] -> SetFillStyle(4000);
      pad[i][j] -> SetFrameFillStyle(4000);
      pad[i][j] -> cd();

      Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
      Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();

      if( i == 0 ) {
	TMultiGraph *mv1 = new TMultiGraph();
	mv1 -> SetName(Form("mv1_%d",j));
	mv1 -> SetTitle(";U_{t0}; v1");
	labelv1[j] = new TLatex(0.4,  0.4, lpid[sqPart[j]]);
       
	for(auto isys: sqSys) {
	  v1Ut = (TGraphErrors*)LoadData(isys, sqPart[j], "g_utv1");
	  if( v1Ut != NULL ) {
	    LOG(INFO) << " v1Ut :: " << v1Ut->GetName() << " is registered." << FairLogger::endl;

	    v1Ut -> SetMarkerStyle(20);
	    v1Ut -> SetMarkerSize(1.);
	    v1Ut -> SetLineColor(icol[isys]);
	    v1Ut -> SetMarkerColor(icol[isys]);
	    v1Ut -> GetXaxis()->SetRangeUser(0.0,2.2);
	    
	    mv1 -> Add( v1Ut, "pl" );
	    
	    if( j == 0 )
	      lg1 -> AddEntry( v1Ut, fsys[isys] );
	  }
	}

	mv1->GetXaxis()->SetLabelFont(43);
	mv1->GetXaxis()->SetLabelSize(16);
	mv1->GetXaxis()->SetLabelOffset(0.02);
	mv1->GetXaxis()->SetTitleFont(43);
	mv1->GetXaxis()->SetTitleSize(18);
	mv1->GetXaxis()->SetTitleOffset(4);
	mv1->GetXaxis()->CenterTitle();
	//	mv1->GetXaxis()->SetNdivisions(505);
	mv1->GetXaxis()->SetRangeUser(0.0, 2.2);
	cout << " x min " << mv1->GetXaxis()->GetXmin() << " " << mv1->GetXaxis()->GetXmax() << endl;
	mv1->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);

	mv1->GetYaxis()->SetLabelFont(43);
	mv1->GetYaxis()->SetLabelSize(16);
	mv1->GetYaxis()->SetLabelOffset(0.02);
	mv1->GetYaxis()->SetTitleFont(43);
	mv1->GetYaxis()->SetTitleSize(18);
	mv1->GetYaxis()->SetTitleOffset(3.5);
	mv1->GetYaxis()->CenterTitle();
	mv1->GetYaxis()->SetNdivisions(505);
	mv1->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
	//mv1->GetYaxis()->SetRangeUser(-0.05,0.6);

	mv1 -> Draw("AP");
	if( j == 0 ) {
	  lg1 -> Draw();
	}

	labelv1[j] -> SetTextSize(0.07*yFactor);
	labelv1[j] -> Draw();

	if( j == Ny - 1 ){
	  TLatex *v1label = new TLatex(0.4,0.63,"0.25 < b_{0} < 0.45; 0.4 < y_{norm} < 0.8");
	  v1label -> SetTextSize(0.08);
	  v1label -> Draw();; 
	}
      }
      else {

	TMultiGraph *mv2 = new TMultiGraph();
	mv2 -> SetName(Form("mv2_%d",j));
	mv2 -> SetTitle(";U_{t0}; v2");
	labelv2[j] = new TLatex(0.5,  -0.1, lpid[sqPart[j]]);

	for(auto isys: sqSys) {
	  v2Ut = (TGraphErrors*)LoadData(isys, sqPart[j], "g_utv2");

	  if( v2Ut != NULL ) {
	    LOG(INFO) << " v2Ut " << v2Ut->GetName() << " is registered." << FairLogger::endl;

	    v2Ut -> SetMarkerStyle(20);
	    v2Ut -> SetMarkerSize(1.);
	    v2Ut -> SetLineColor(icol[isys]);
	    v2Ut -> SetMarkerColor(icol[isys]);

	    if( mv2 )
	      mv2 -> Add( v2Ut, "pl" );

	    if( j == 0 )
	      lg2 -> AddEntry(v2Ut, fsys[isys]);

	  }
	  else 
	    LOG(INFO) << " v2ut is not found." << FairLogger::endl;
	}      
	
	mv2->GetYaxis()->SetLabelFont(43);
	mv2->GetYaxis()->SetLabelSize(16);
	mv2->GetYaxis()->SetLabelOffset(0.02);
	mv2->GetYaxis()->SetTitleFont(43);
	mv2->GetYaxis()->SetTitleSize(16);
	mv2->GetYaxis()->SetTitleOffset(3.5);
	mv2->GetYaxis()->CenterTitle();
	mv2->GetYaxis()->SetNdivisions(505);
	//	mv2->GetYaxis()->SetRangeUser(-0.25,0.01);

	// TICKS Y Axis                                                                                                                             
	mv2->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);

	// Format for x axis                                                                                                                        
	mv2->GetXaxis()->SetLabelFont(43);
	mv2->GetXaxis()->SetLabelSize(16);
	mv2->GetXaxis()->SetLabelOffset(0.02);
	mv2->GetXaxis()->SetTitleFont(43);
	mv2->GetXaxis()->SetTitleSize(16);
	mv2->GetXaxis()->SetTitleOffset(4);
	mv2->GetXaxis()->CenterTitle();
	mv2->GetXaxis()->SetNdivisions(505);
	mv2->GetXaxis()->SetRangeUser(0.0,2.2);

	// TICKS X Axis                                                                                                                             
	mv2->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);

	mv2-> Draw("AP");
	//	mv2-> Print();
	if( j == 0 ) {
          lg2 -> Draw();
        }

	labelv2[j] -> SetTextSize(0.07*yFactor);
	labelv2[j] -> Draw();

	if( i == 0 )
	  lg2 -> Draw();
	//mv2 -> Print();

	if( j == Ny - 1 ){
	  TLatex *v2label = new TLatex(0.4,0.02,"0.25 < b_{0} < 0.45; -0.4 < y_{norm} < 0.4");
	  v2label -> SetTextSize(0.08);
	  v2label -> Draw();; 
	}
      }
    }
}


void Draw_Ut_SystemD()
{
  LOG(INFO) << "Draw_Ut_SystemD " << FairLogger::endl;

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 1000, 800); iccv++;
  gStyle->SetOptTitle(0);

  const Int_t Nx = 2;
  const Int_t Ny = 3;

  Float_t lMargin = 0.08;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.10;
  Float_t mMargin = 0.08;

  CanvasPartitionTwoColumn(ccv,Nx,Ny,lMargin,rMargin,bMargin,tMargin,mMargin);

  TPad *pad[Nx][Ny];

  std::vector<UInt_t> sqSys   = {1,3,0};
  for( auto i : ROOT::TSeqI(Nx) ) for( auto j : ROOT::TSeqI(Ny) ) {
      //  std::vector<UInt_t> vv={3};
      //  for( auto i : vv ) for( auto j : vv ) {
      ccv->cd(0);

      TString pname = Form("pad_%i_%i",i,j);
      pad[i][j] = (TPad*)gROOT->FindObject(pname);
      cout << " pad " << pname << endl;
      if( pad[i][j] == NULL ) {
       	cout << " pad is not found " << pname << endl;
       	continue;
      }
      pad[i][j] -> Draw();
      pad[i][j] -> SetFillStyle(4000);
      pad[i][j] -> SetFrameFillStyle(4000);
      pad[i][j] -> cd();

      Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
      Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();

      if( i == 0 ) {

	auto mv1 = new TMultiGraph(Form("mv1_%d",sqSys[j]),";U_{t0}; v1");
	auto lg1 = new TLegend(0.75,0.3,0.9,0.55,"");
	auto labelv1 = new TLatex(1.6,  0.50, fsys[sqSys[j]]);

	for(auto ipart: sqPart ) {

	  auto v1Ut = (TGraphErrors*)LoadData(sqSys[j], ipart, "g_utv1");
	  if( v1Ut != NULL ) {

	    v1Ut -> SetMarkerStyle(20);
	    v1Ut -> SetMarkerSize(1.);
	    v1Ut -> SetLineColor(icol[ipart]);
	    v1Ut -> SetMarkerColor(icol[ipart]);
	    if( mv1 ) {
	      mv1 -> Add( v1Ut, "pl" );
	      LOG(INFO) << " v1Ut :: " << v1Ut->GetName() << " is registered." << FairLogger::endl;
	    }

	    if( i == 0 )
	      lg1 -> AddEntry(v1Ut,lpid[ipart]);
	
	  }
	  else
	    LOG(ERROR) << " v1ut  is not found." << FairLogger::endl;
	}

	mv1->GetYaxis()->SetLabelFont(43);
	mv1->GetYaxis()->SetLabelSize(16);
	mv1->GetYaxis()->SetLabelOffset(0.02);
	mv1->GetYaxis()->SetTitleFont(43);
	mv1->GetYaxis()->SetTitleSize(18);
	mv1->GetYaxis()->SetTitleOffset(3.5);
	mv1->GetYaxis()->CenterTitle();
	mv1->GetYaxis()->SetNdivisions(505);
	mv1->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
	mv1->GetYaxis()->SetRangeUser(-0.05,0.62);

	mv1->GetXaxis()->SetLabelFont(43);
	mv1->GetXaxis()->SetLabelSize(16);
	mv1->GetXaxis()->SetLabelOffset(0.02);
	mv1->GetXaxis()->SetTitleFont(43);
	mv1->GetXaxis()->SetTitleSize(18);
	mv1->GetXaxis()->SetTitleOffset(4);
	mv1->GetXaxis()->CenterTitle();
	mv1->GetXaxis()->SetNdivisions(505);
	mv1->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);

	mv1-> Draw("AP");
	if( j == 0 ) {
	  lg1 -> Draw();
	}

	labelv1 -> SetTextSize(0.07*yFactor);
	labelv1 -> Draw();

	if( j == Ny - 1 ){
	  TLatex *v1label = new TLatex(0.4,0.65,"0.25 < b_{0} < 0.45; 0.4 < y_{norm} < 0.8");
	  v1label -> SetTextSize(0.08);
	  v1label -> Draw();; 
	}
      }
      else {  // Right side panel

	auto mv2 = new TMultiGraph(Form("mv2_%d",sqSys[j]),";U_{t0}; v2");
	auto lg2 = new TLegend(0.25,0.3,0.4,0.55,"");
	auto labelv2 = new TLatex(1.6,  0.50, fsys[sqSys[j]]);

	for(auto ipart: sqPart ) {
	  auto v2Ut = (TGraphErrors*)LoadData(sqSys[j], ipart, "g_utv2");

	  if( v2Ut != NULL ) {

	    v2Ut -> SetMarkerStyle(20);
	    v2Ut -> SetMarkerSize(1.);
	    v2Ut -> SetLineColor(icol[ipart]);
	    v2Ut -> SetMarkerColor(icol[ipart]);

	    if( mv2 ) {
	      mv2 -> Add( v2Ut, "pl" );
	      LOG(INFO) << " v2Ut " << v2Ut->GetName() << " is registered." << FairLogger::endl;
	    }

	    if( i == 0 )
	      lg2 -> AddEntry(v2Ut,lpid[ipart]);
	  }
	  else 
	    LOG(INFO) << " v2ut is not found." << FairLogger::endl;
	}

	mv2->GetYaxis()->SetLabelFont(43);
	mv2->GetYaxis()->SetLabelSize(16);
	mv2->GetYaxis()->SetLabelOffset(0.02);
	mv2->GetYaxis()->SetTitleFont(43);
	mv2->GetYaxis()->SetTitleSize(16);
	mv2->GetYaxis()->SetTitleOffset(3.5);
	mv2->GetYaxis()->CenterTitle();
	mv2->GetYaxis()->SetNdivisions(505);
	mv2->GetYaxis()->SetRangeUser(-0.25,0.01);

	// TICKS Y Axis                                                                                                                             
	mv2->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
	  
	// Format for x axis                                                                                                                        
	mv2->GetXaxis()->SetLabelFont(43);
	mv2->GetXaxis()->SetLabelSize(16);
	mv2->GetXaxis()->SetLabelOffset(0.02);
	mv2->GetXaxis()->SetTitleFont(43);
	mv2->GetXaxis()->SetTitleSize(16);
	mv2->GetXaxis()->SetTitleOffset(4);
	mv2->GetXaxis()->CenterTitle();
	mv2->GetXaxis()->SetNdivisions(505);

	// TICKS X Axis                                                                                                                             
	mv2->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
	
	mv2-> Draw("AP");
      
	labelv2 -> SetTextSize(0.07*yFactor);
	labelv2 -> Draw();
	  
	if( j == 0 )
	  lg2 -> Draw();
	//mv2 -> Print();
	
	if( j == Ny - 1 ){
	  TLatex *v2label = new TLatex(0.4,0.02,"0.25 < b_{0} < 0.45; -0.4 < y_{norm} < 0.4");
	  v2label -> SetTextSize(0.08);
	  v2label -> Draw();; 
	}
	
      }
    }
}


void Draw_v1Ut_SystemD(UInt_t ipart)
{
  LoadData();
  std::vector<UInt_t> sqSys   = {1,3,0};

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 800); iccv++;
  gStyle->SetOptTitle(0);

  const Int_t Ny = 4;

  Float_t lMargin = 0.02;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.05;
  CanvasPartitionY(ccv,Ny,lMargin,rMargin,bMargin,tMargin);

  TPad *pad[Ny];

  ccv->Divide((Int_t)sqv1sel.size(), 1);
  //  ccv->Divide(6, 1);

  TLegend     *lgv1 = new TLegend(0.5, 0.18, 0.85, 0.31);

  UInt_t ii = 0;
  for(auto i: sqv1sel ) {
    auto mv1Ut = new TMultiGraph(Form("mv1Ut_%d",i) ,";Ut ; v1");

    // gROOT->ls();
      
    for(auto k:sqSys) {
      TString gname = (TString)Form("gUt_v1%d",i) +  "_" + rsys[k] + "_" + lpid[ipart];
      v1Ut = (TGraphErrors*)gROOT->Get(gname);
      //      v1Ut -> SetName(gname+"cp");

      if( v1Ut != NULL ) {
	cout << " v1ut " << gname << " is registered." << endl;

	v1Ut -> SetMarkerStyle(20);
	v1Ut -> SetMarkerSize(1.);
	v1Ut -> SetLineColor(icol[k]);
	v1Ut -> SetMarkerColor(icol[k]);

	mv1Ut -> Add( v1Ut, "pl" );
	if( ii == 0)
	  lgv1  -> AddEntry( v1Ut,  rsys[k]+"-"+lpid[ipart] , "lp");
      }
      else 
	cout << " v1Ut " << gname << " not be found. " << endl;
    }
    ccv->cd(0);
    TString pname = Form("pad_%i",ii);
    pad[ii] = (TPad*)gROOT->FindObject(pname);
    pad[ii] -> Draw();
    pad[ii] -> SetFillStyle(4000);
    pad[ii] -> SetFrameFillStyle(4000);
    pad[ii] -> cd();

    mv1Ut -> SetTitle(v1Ut->GetTitle());
    mv1Ut -> Draw("AP");
    if( ii == 0)
      lgv1  -> Draw();

    ii++;
  }
}

void Draw_v2Ut_SystemD(UInt_t ipart)
{
  std::vector<UInt_t> sqSys   = {1,3,0};

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 800); iccv++;
  gStyle->SetOptTitle(0);

  const Int_t Ny = 4;

  Float_t lMargin = 0.02;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.05;
  CanvasPartitionY(ccv,Ny,lMargin,rMargin,bMargin,tMargin);

  TPad *pad[Ny];

  ccv->Divide((Int_t)sqv2sel.size(), 1);
  //  ccv->Divide(6, 1);

  TLegend  *lgv2 = new TLegend(0.8, 0.75, 0.9, 0.98);

  UInt_t ii = 0;
  for(auto i: sqv2sel ) {
    auto mv2Ut = new TMultiGraph(Form("mv2Ut_%d",i) ,";Ut ; v2");

    // gROOT->ls();
      
    for(auto k:sqSys) {
      TString gname = (TString)Form("gUt_v2%d",i) +  "_" + rsys[k] + "_" + lpid[ipart];
      v2Ut = (TGraphErrors*)gROOT->Get(gname);
      //      v2Ut -> SetName(gname+"cp");

      if( v2Ut != NULL ) {
	cout << " v2ut " << gname << " is registered." << endl;

	v2Ut -> SetMarkerStyle(20);
	v2Ut -> SetMarkerSize(1.);
	v2Ut -> SetLineColor(icol[k]);
	v2Ut -> SetMarkerColor(icol[k]);

	mv2Ut -> Add( v2Ut, "pl" );
	if( ii == 0)
	  lgv2  -> AddEntry( v2Ut,  rsys[k]+"-"+lpid[ipart] , "lp");
      }
      else 
	cout << " v2Ut " << gname << " not be found. " << endl;
    }
    ccv->cd(0);
    TString pname = Form("pad_%i",ii);
    pad[ii] = (TPad*)gROOT->FindObject(pname);
    pad[ii] -> Draw();
    pad[ii] -> SetFillStyle(4000);
    pad[ii] -> SetFrameFillStyle(4000);
    pad[ii] -> cd();

    mv2Ut -> SetTitle(v2Ut->GetTitle());
    mv2Ut -> Draw("AP");
    TLatex *Lv1 = new TLatex(-0.25,0.2,(TString)v2Ut->GetTitle() );
    cout << v2Ut->GetTitle() << endl;
    Lv1->Draw();

    if( ii == 0)
      lgv2  -> Draw();

    ii++;
  }

}



void Draw_v2SystemD_one()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 700, 500); iccv++;
  auto gmv2n = new TMultiGraph("gmv2n", ";"+xlabel+"; v2n");
  auto lmv2n = new TLegend(0.35,0.13,0.7,0.33,"");

  gmv2n -> Add( g_v2n[0],  "p");
  gmv2n -> Add( g_v2n[1],  "p");
  gmv2n -> Add( g_v2n[2],  "p");
  gmv2n -> Add( g_v2n[3],  "p");

  lmv2n -> AddEntry( g_v2n[0], lpid[0], "p");
  lmv2n -> AddEntry( g_v2n[1], lpid[1], "p");
  lmv2n -> AddEntry( g_v2n[2], lpid[2], "p");
  lmv2n -> AddEntry( g_v2n[3], lpid[3], "p");

  gmv2n -> Draw("ALP");
  lmv2n -> Draw();
  
}
void Draw_v2SystemD()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;

  Double_t BYs = 0.1;            // Bottom Y space
  Double_t BYm = 0.1;            // Bottom Y mergin
  Double_t TYm = 0.01;            // Top Y mergin

  Double_t Ny  = 4;               // Number of pads along Y  
  Double_t H   = (1.0 - (TYm+BYm+BYs))/Ny; // pad height
  cout << " H= " << H << endl;

  Double_t Yl = BYs;  Double_t Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p1 = new TPad("p1", "p1", 0., Yl, 1.,  Yu, 0, 0, 0);
  p1->SetTopMargin(0);
  p1->SetBottomMargin(BYm);
  p1->Draw();


  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p2 = new TPad("p2", "p2", 0., Yl, 1., Yu, 0, 0, 0);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0);
  p2->Draw();

  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p3 = new TPad("p3", "p3", 0., Yl, 1., Yu, 0, 0, 0);  
  p3->SetTopMargin(0);
  p3->SetBottomMargin(0);
  p3->Draw();

  Yl = Yu;  Yu = Yl +  H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p4 = new TPad("p4", "p4", 0., Yl, 1., Yu, 0, 0, 0);  
  p4->SetTopMargin(TYm);
  p4->SetBottomMargin(0);
  p4->Draw();


  //  gStyle->SetLabelSize(0.08);
  g_v2max[0]->SetTitle(";"+xlabel+";"); 
  g_v2max[3]->SetTitle(";; v_{20} ");

  p1->cd(); 
  g_v2max[0] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[0] -> GetYaxis() -> SetLabelSize(0.08);
  g_v2max[0] -> GetXaxis() -> SetLabelSize(0.1);
  g_v2max[0] -> GetXaxis() -> SetTitleOffset(1.);
  g_v2max[0] -> GetXaxis() -> SetTitleSize(1);
  //  g_v2max[0] -> GetXaxis() -> SetTitleOffSet(0.5);
  g_v2max[0] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2, 0.9, fpid[0]);


  p2->cd();
  g_v2max[1] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[1] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2max[1] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[1]);

  p3->cd(); 
  g_v2max[2] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[2] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2max[2] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[2]);

  p4->cd(); 
  g_v2max[3] -> GetYaxis() -> SetTitleOffset(1.);
  g_v2max[3] -> GetYaxis() -> SetTitleSize(0.1);
  g_v2max[3] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[3] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2max[3] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[3]);

}

void Draw_v2nSystemD()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;

  Double_t BYs = 0.1;            // Bottom Y space
  Double_t BYm = 0.1;            // Bottom Y mergin
  Double_t TYm = 0.01;            // Top Y mergin

  Double_t Ny  = 4;               // Number of pads along Y  
  Double_t H   = (1.0 - (TYm+BYm+BYs))/Ny; // pad height
  cout << " H= " << H << endl;

  Double_t Yl = BYs;  Double_t Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p1 = new TPad("p1", "p1", 0., Yl, 1.,  Yu, 0, 0, 0);
  p1->SetTopMargin(0);
  p1->SetBottomMargin(BYm);
  p1->Draw();


  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p2 = new TPad("p2", "p2", 0., Yl, 1., Yu, 0, 0, 0);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0);
  p2->Draw();

  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p3 = new TPad("p3", "p3", 0., Yl, 1., Yu, 0, 0, 0);  
  p3->SetTopMargin(0);
  p3->SetBottomMargin(0);
  p3->Draw();

  Yl = Yu;  Yu = Yl +  H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p4 = new TPad("p4", "p4", 0., Yl, 1., Yu, 0, 0, 0);  
  p4->SetTopMargin(TYm);
  p4->SetBottomMargin(0);
  p4->Draw();


  //  gStyle->SetLabelSize(0.08);
  g_v2n[0]->SetTitle(";"+xlabel+";"); 
  g_v2n[3]->SetTitle(";; v_{21} ");

  p1->cd(); 
  g_v2n[0] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2n[0] -> GetYaxis() -> SetLabelSize(0.08);
  g_v2n[0] -> GetXaxis() -> SetLabelSize(0.1);
  g_v2n[0] -> GetXaxis() -> SetTitleOffset(1.);
  g_v2n[0] -> GetXaxis() -> SetTitleSize(1);
  //  g_v2n[0] -> GetXaxis() -> SetTitleOffSet(0.5);
  g_v2n[0] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2, 0.9, fpid[0]);


  p2->cd();
  g_v2n[1] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2n[1] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2n[1] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[1]);

  p3->cd(); 
  g_v2n[2] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2n[2] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2n[2] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[2]);

  p4->cd(); 
  g_v2n[3] -> GetYaxis() -> SetTitleOffset(1.);
  g_v2n[3] -> GetYaxis() -> SetTitleSize(0.1);
  g_v2n[3] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2n[3] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2n[3] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[3]);

}


void Draw_v1Ratio()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;

  Double_t BYs = 0.1;            // Bottom Y space
  Double_t TYm = 0.05;            // Top Y mergin
  Double_t BYm = 0.1;            // Bottom Y mergin
  Double_t Ny  = 3;               // Number of pads along Y  
  Double_t H   = (1.0 - (TYm+BYm+BYs))/Ny; // pad height

  Double_t Yl = BYs;   Double_t Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *pp1 = new TPad("pp1", "pp1", 0., Yl, 1.,  Yu, 0, 0, 0);
  pp1->SetTopMargin(0);
  pp1->Draw();

  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *pp2 = new TPad("pp2", "pp2", 0., Yl, 1., Yu, 0, 0, 0);
  pp2->SetTopMargin(0);
  pp2->SetBottomMargin(0);
  pp2->Draw();

  Yl = Yu;  Yu = Yl + H ;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *pp3 = new TPad("pp3", "pp3", 0., Yl, 1., Yu, 0, 0, 0);  
  pp3->SetBottomMargin(0);
  pp3->Draw();

  plabel.SetTextSize(0.1);

  g_v1pslp[1]->SetTitle(";"+xlabel+";"); 
  g_v1pslp[3]->SetTitle(";; v11/ v11_{p} ");

  pp1->cd(); 
  plabel.DrawLatexNDC(0.2,0.9,fpid[1]);
  g_v1pslp[1] -> GetYaxis() -> SetLabelSize(0.05);
  g_v1pslp[1] -> GetXaxis() -> SetLabelSize(0.05);
  g_v1pslp[1] -> Draw("AP");  

  pp2->cd();
  plabel.DrawLatexNDC(0.2,0.9,fpid[3]);
  g_v1pslp[2] -> GetYaxis() -> SetLabelSize(0.05);
  g_v1pslp[2] -> Draw("AP");  

  pp3->cd(); 
  plabel.DrawLatexNDC(0.2,0.9,fpid[4]);
  g_v1pslp[3] -> GetYaxis() -> SetLabelSize(0.05);
  g_v1pslp[3] -> Draw("AP");  
  //---
} 


void Draw_v1v2SystemD_fit()
{
  LoadData();

  plabel.SetTextAlign(13);
  plabel.SetTextSize(0.08);

  for(UInt_t i = 0; i < 4; i++ ){

    if( mrv1[i]->GetListOfGraphs() != NULL ) {
      ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;

      mrv1[i]->SetTitle();

      mrv1[i]->GetXaxis()->SetRangeUser(-0.8, 1.2);
      mrv1[i]->GetYaxis()->SetRangeUser(-0.5, 0.6);
      mrv1[i]->Draw("ALP");
      lrv1[i]->Draw();
    }

    auto aLineX1 = new TLine( mrv1[i]->GetXaxis()->GetXmin(), 0.,
			      mrv1[i]->GetXaxis()->GetXmax(), 0.);
    auto aLineY1 = new TLine( 0., mrv1[i]->GetYaxis()->GetXmin(), 
			      0., mrv1[i]->GetYaxis()->GetXmax());

    aLineX1->SetLineColor(1);
    aLineX1->SetLineStyle(3);
    aLineX1->Draw();
    
    aLineY1->SetLineColor(1);
    aLineY1->SetLineStyle(3);
    aLineY1->Draw();
  }



  for(UInt_t i = 0; i <4; i++ ){

    if( mrv2[i]->GetListOfGraphs() != NULL ) {
      ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;

      mrv2[i]->GetYaxis()->SetRangeUser(-0.13, 0.02);
      mrv2[i]->Draw("ALP");
      lrv2[i]->Draw();
    }
  }
}

void Draw_v1v2SystemD()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  g_v2sysD -> Draw("ALP");
  l_v2sysD -> Draw();


  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  g_v1psysD -> Draw("ALP");
  l_v1psysD -> Draw();


  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  g_v1sysD -> Draw("ALP");
  l_v1sysD -> Draw();
  //<<<---
}

void Draw_Indiv_v1SystemD_Three() 
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 640); iccv++;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const Int_t Ny = 4;

  Float_t lMargin = 0.02;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.05;

  CanvasPartitionY(ccv,Ny,lMargin,rMargin,bMargin,tMargin);

  TPad *pad[Ny];
  
  for(Int_t j = 0; j < Ny; j++ ) {
      
    ccv->cd(0);
    char pname[16];
    sprintf(pname,"pad_%i",j);
    pad[j] = (TPad*)gROOT->FindObject(pname);
    pad[j] -> Draw();
    pad[j] -> SetFillStyle(4000);
    pad[j] -> SetFrameFillStyle(4000);
    pad[j] -> cd();

    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[j]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[j]->GetAbsHNDC();


    // if( j == 2 ) {
    // 	auto mg_T3He = new TMultiGraph("mg_T3He",";;v_{11}");
    // 	auto l_T3He  = new TLegend(0.2,0.7,0.4,0.9,"");
    // 	mg_T3He->Add(g_v1slp[2], "AP");
    // 	mg_T3He->Add(g_v1slp[3], "AP");
    // 	l_T3He->AddEntry(g_v1slp[2],"Triton");
    // 	l_T3He->AddEntry(g_v1slp[3],"^{3}He");

    // 	mg_T3He -> GetXaxis() -> SetLimits(g_v1slp[0]->GetXaxis()->GetXmin(), g_v1slp[0]->GetXaxis()->GetXmax());
    // 	mg_T3He -> SetMaximum(0.794);
    // 	// mg_T3He -> GetYaxis() -> SetTitleOffset(0.75);
    // 	// mg_T3He -> GetYaxis() -> SetTitleSize(0.1);
    // 	// mg_T3He -> GetYaxis() -> SetLabelOffset(0.01);
    // 	// mg_T3He -> GetYaxis() -> SetLabelSize(0.1);
    // 	mg_T3He -> Draw("AP");  
    // 	l_T3He  -> Draw();
    // }

    // else {
      
    //@@@
    // Format for y axis
    g_v1slp[j] -> GetYaxis()->SetLabelFont(43);
    g_v1slp[j] -> GetYaxis()->SetLabelSize(20);
    g_v1slp[j] -> GetYaxis()->SetLabelOffset(0.02);
    g_v1slp[j] -> GetYaxis()->SetTitleFont(43);
    g_v1slp[j] -> GetYaxis()->SetTitleSize(20);
    g_v1slp[j] -> GetYaxis()->SetTitleOffset(2.5);
    g_v1slp[j] -> GetYaxis()->CenterTitle();
    g_v1slp[j] -> GetYaxis()->SetNdivisions(504);
    g_v1slp[j] -> GetYaxis()->SetTickLength(xFactor*0.04/yFactor);

    cout << j << " min " << g_v1slp[j]->GetEXlow() << " max " << g_v1slp[j]->GetEX() << endl; 
    Float_t Ymin = g_v1slp[j]->GetYaxis()->GetXmin();
    Ymin = Int_t(Ymin/0.01);
    Float_t Ymax = g_v1slp[j]->GetYaxis()->GetXmax();
    Ymax = Int_t(Ymax/0.01);
    g_v1slp[j] -> GetYaxis() -> SetRangeUser((Ymin*0.01)-0.003, (Ymax*0.01+0.003));

    // Format for x axis
    g_v1slp[j] -> GetXaxis()->SetLabelFont(43);
    g_v1slp[j] -> GetXaxis()->SetLabelSize(20);
    g_v1slp[j] -> GetXaxis()->SetLabelOffset(0.02);
    g_v1slp[j] -> GetXaxis()->SetTitleFont(43);
    g_v1slp[j] -> GetXaxis()->SetTitleSize(16);
    g_v1slp[j] -> GetXaxis()->SetTitleOffset(5);
    g_v1slp[j] -> GetXaxis()->SetNdivisions(505);
 
    // TICKS X Axis
    g_v1slp[j] -> GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
      
    g_v1slp[j] -> Draw("AP");

    cout << yFactor << endl;
    plabel.SetTextAlign(13);
    plabel.SetTextSize(0.09*yFactor);

    if( j == 0 )  plabel.DrawLatexNDC(0.25, 0.5/yFactor, fpid[j]);
    else
      plabel.DrawLatexNDC(0.25, 0.7, fpid[j]);

    //      }
  }
}

void Draw_Indiv_v1SystemD() 
{
  LOG(INFO) << "Draw_Indiv_v1SystemD ... " << FairLogger::endl;

  LoadData();

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 850); iccv++;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const Int_t Ny = 5;

  Float_t lMargin = 0.02;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.05;

  CanvasPartitionY(ccv,Ny,lMargin,rMargin,bMargin,tMargin);

  TPad *pad[Ny];
  
  for(Int_t j = 0; j < Ny; j++ ) {
      
    ccv->cd(0);
    TString pname = Form("pad_%i",j);
    pad[j] = (TPad*)gROOT->FindObject(pname);
    pad[j] -> Draw();
    pad[j] -> SetFillStyle(4000);
    pad[j] -> SetFrameFillStyle(4000);
    pad[j] -> cd();

    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[j]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[j]->GetAbsHNDC();


    // Format for y axis
    g_v1slp[j] -> GetYaxis()->SetLabelFont(43);
    g_v1slp[j] -> GetYaxis()->SetLabelSize(20);
    g_v1slp[j] -> GetYaxis()->SetLabelOffset(0.02);
    g_v1slp[j] -> GetYaxis()->SetTitleFont(43);
    g_v1slp[j] -> GetYaxis()->SetTitleSize(20);
    g_v1slp[j] -> GetYaxis()->SetTitleOffset(3.);
    g_v1slp[j] -> GetYaxis()->CenterTitle();
    g_v1slp[j] -> GetYaxis()->SetNdivisions(504);
    g_v1slp[j] -> GetYaxis()->SetTickLength(xFactor*0.04/yFactor);

    Float_t Ymin = g_v1slp[j]->GetYaxis()->GetXmin();
    Ymin = Int_t(Ymin/0.01);
    Float_t Ymax = g_v1slp[j]->GetYaxis()->GetXmax();
    Ymax = Int_t(Ymax/0.01);
    //    g_v1slp[j] -> GetYaxis() -> SetRangeUser((Ymin*0.01)-0.003, (Ymax*0.01+0.003));

    // Format for x axis
    g_v1slp[j] -> GetXaxis()->SetLabelFont(43);
    g_v1slp[j] -> GetXaxis()->SetLabelSize(20);
    g_v1slp[j] -> GetXaxis()->SetLabelOffset(0.02);
    g_v1slp[j] -> GetXaxis()->SetTitleFont(43);
    g_v1slp[j] -> GetXaxis()->SetTitleSize(16);
    g_v1slp[j] -> GetXaxis()->SetTitleOffset(5);
    g_v1slp[j] -> GetXaxis()->SetNdivisions(505);
 
    // TICKS X Axis
    g_v1slp[j] -> GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
      
    g_v1slp[j] -> Draw("AP");

    cout << yFactor << endl;
    plabel.SetTextAlign(13);
    plabel.SetTextSize(0.09*yFactor);

    if( j == 0 )  plabel.DrawLatexNDC(0.3, 0.5/yFactor, fpid[j]);
    else
      plabel.DrawLatexNDC(0.3, 0.7, fpid[j]);

    //      }
  }
}


void Draw_Indiv_v2SystemD() 
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 850); iccv++;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const Int_t Ny = 5;

  Float_t lMargin = 0.02;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.05;

  CanvasPartitionY(ccv,Ny,lMargin,rMargin,bMargin,tMargin);

  TPad *pad[Ny];
  
  for(Int_t j = 0; j < Ny; j++ ) {
      
    ccv->cd(0);
    char pname[16];
    sprintf(pname,"pad_%i",j);
    pad[j] = (TPad*)gROOT->FindObject(pname);
    pad[j] -> Draw();
    pad[j] -> SetFillStyle(4000);
    pad[j] -> SetFrameFillStyle(4000);
    pad[j] -> cd();

    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[j]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[j]->GetAbsHNDC();


    // Format for y axis
    g_v2max[j] -> GetYaxis()->SetLabelFont(43);
    g_v2max[j] -> GetYaxis()->SetLabelSize(20);
    g_v2max[j] -> GetYaxis()->SetLabelOffset(0.02);
    g_v2max[j] -> GetYaxis()->SetTitleFont(43);
    g_v2max[j] -> GetYaxis()->SetTitleSize(20);
    g_v2max[j] -> GetYaxis()->SetTitleOffset(3);
    g_v2max[j] -> GetYaxis()->CenterTitle();
    g_v2max[j] -> GetYaxis()->SetNdivisions(504);
    g_v2max[j] -> GetYaxis()->SetTickLength(xFactor*0.04/yFactor);

    cout << j << " min " << g_v2max[j]->GetEXlow() << " max " << g_v2max[j]->GetEX() << endl; 
    Float_t Ymin = g_v2max[j]->GetYaxis()->GetXmin();
    Ymin = Int_t(Ymin/0.001);
    Float_t Ymax = g_v2max[j]->GetYaxis()->GetXmax();
    Ymax = Int_t(Ymax/0.001);
    g_v2max[j] -> GetYaxis() -> SetRangeUser((Ymin*0.001)-0.0001, (Ymax*0.001+0.0006));

    // Format for x axis
    g_v2max[j] -> GetXaxis()->SetLabelFont(43);
    g_v2max[j] -> GetXaxis()->SetLabelSize(20);
    g_v2max[j] -> GetXaxis()->SetLabelOffset(0.02);
    g_v2max[j] -> GetXaxis()->SetTitleFont(43);
    g_v2max[j] -> GetXaxis()->SetTitleSize(16);
    g_v2max[j] -> GetXaxis()->SetTitleOffset(5);
    g_v2max[j] -> GetXaxis()->SetNdivisions(505);
 
    // TICKS X Axis
    g_v2max[j] -> GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
      
    g_v2max[j] -> Draw("AP");

    cout << yFactor << endl;
    plabel.SetTextAlign(13);
    plabel.SetTextSize(0.09*yFactor);

    if( j == 0 )  plabel.DrawLatexNDC(0.3, 0.5/yFactor, fpid[j]);
    else
      plabel.DrawLatexNDC(0.3, 0.5, fpid[j]);

    //      }
  }
}

void Draw_v_y(UInt_t vn = 1) 
{
  if( vn > 2 ) vn = 1;

  LOG(INFO) << "Draw_v" << vn << "_y"  << FairLogger::endl;

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 1000, 800); iccv++;
  gStyle->SetOptTitle(0);

  std::vector<UInt_t> sqSys   = {0,1};
  std::vector<UInt_t> sqPart  = {0,1,2,3};
  const Int_t Nx = (Int_t)sqSys.size();
  const Int_t Ny = (Int_t)sqPart.size();

  Double_t vRange[2][5][2] ={ { {-0.24, 0.24}, {-0.4 , 0.4},  {-0.6 , 0.6 }, {-0.8 , 0.8},  {-1.0 , 1.0}},
			       { {-0.05, 0.01}, {-0.08, 0.01}, {-0.10, 0.01}, {-0.10, 0.01}, {-0.10, 0.01} }};

  Float_t lMargin = 0.08;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.10;
  Float_t mMargin = 0.08;

  TMultiGraph *mv;
  TGraphErrors* v_y;
  
  CanvasPartitionTwoColumn(ccv,Nx,Ny,lMargin,rMargin,bMargin,tMargin,mMargin);

  TPad *pad[Nx][Ny];


  for(auto i: ROOT::TSeqI(Nx) )for(auto j: ROOT::TSeqI(Ny) ) {
      ccv->cd(0);


      TString pname = Form("pad_%i_%i",i,j);
      pad[i][j] = (TPad*)gROOT->FindObject(pname);
      cout << " pad " << pname << endl;
      if( pad[i][j] == NULL ) {
        cout << " pad is not found " << pname << endl;
        continue;
      }
      pad[i][j] -> Draw();
      pad[i][j] -> SetFillStyle(4000);
      pad[i][j] -> SetFrameFillStyle(4000);
      pad[i][j] -> cd();
      Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
      Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();


      TMultiGraph *mv = new TMultiGraph();
      mv -> SetName(Form("mv%d_%d",i,j));
      mv -> SetTitle(";y; v1");

      TLegend *lg = new TLegend(0.6,0.3/yFactor,0.9,0.55/yFactor,lpid[sqPart[j]]);

      v_y = LoadData(sqSys[i], sqPart[j], Form("gu_v%d",vn));

      if( v_y != NULL ) {
	v_y -> SetMarkerStyle(20);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor(icol[sqPart[j]]);
	v_y -> SetMarkerColor(icol[sqPart[j]]);

	mv -> Add( v_y, "pl" );
	lg -> AddEntry(v_y, "Data : "+lsys[sqSys[i]]);
      }
      else
	LOG(ERROR) << " gu_v1  is not found. " << sqSys[i] << " : " << sqPart[j] << FairLogger::endl;
      

      v_y = LoadAMD(sqSys[i], sqPart[j], Form("v%d_y",vn));
      
      if( v_y != NULL ) {
	v_y -> SetMarkerStyle(25);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor(icol[sqPart[j]]);
	v_y -> SetMarkerColor(icol[sqPart[j]]);

	mv -> Add( v_y, "pl" );
	lg -> AddEntry(v_y, "AMD");
      }

      //      gROOT->Add(mv);
      

      mv->GetYaxis()->SetRangeUser(vRange[vn][j][0], vRange[vn][j][0]);

      mv->Draw("AP");
      //      lg->Draw();
    }  
}

void CanvasPartitionY(TCanvas *C,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
  if (!C) return;
 
  // Setup Pad layout:
  Float_t vSpacing = 0.0;
  Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
 
  Float_t hStep  = (1.- lMargin - rMargin );;

  Float_t vposd,vposu,vmard,vmaru,vfactor;
  Float_t hposl,hposr,hmarl,hmarr,hfactor;
 
  hposl = lMargin;
  hposr = 1.-rMargin;
  hfactor = hposr-hposl;
  hmarl = 0.18;
  hmarr = 0.01;

  for (Int_t j=0;j<Ny;j++) {
 
    if (j==0) {
      vposd = 0.0;
      vposu = bMargin + vStep;
      vfactor = vposu-vposd;
      vmard = bMargin / vfactor;
      vmaru = 0.0;
    } else if (j == Ny-1) {
      vposd = vposu + vSpacing;
      vposu = vposd + vStep + tMargin;
      vfactor = vposu-vposd;
      vmard = 0.0;
	vmaru = tMargin / (vposu-vposd);
    } else {
      vposd = vposu + vSpacing;
      vposu = vposd + vStep;
      vfactor = vposu-vposd;
      vmard = 0.0;
      vmaru = 0.0;
    }
 
    C->cd(0);
 
    char name[16];
    sprintf(name,"pad_%i",j);
    TPad *pad = (TPad*) gROOT->FindObject(name);
    if (pad) delete pad;
    pad = new TPad(name,"",hposl,vposd,hposr,vposu);
    pad->SetLeftMargin(hmarl);
    pad->SetRightMargin(hmarr);
    pad->SetBottomMargin(vmard);
    pad->SetTopMargin(vmaru);
 
    pad->SetFrameBorderMode(0);
    pad->SetBorderMode(0);
    pad->SetBorderSize(0);
    
    pad->Draw();
    
  }
}




 
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
  if (!C) return;
 
  // Setup Pad layout:
  Float_t vSpacing = 0.0;
  Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
 
  Float_t hSpacing = 0.0;
  Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
 
  Float_t vposd,vposu,vmard,vmaru,vfactor;
  Float_t hposl,hposr,hmarl,hmarr,hfactor;
 
  for (Int_t i=0;i<Nx;i++) {
 
    if (i==0) {
      hposl = 0.0;
      hposr = lMargin + hStep;
      hfactor = hposr-hposl;
      hmarl = lMargin / hfactor;
      hmarr = 0.0;
    } else if (i == Nx-1) {
      hposl = hposr + hSpacing;
      hposr = hposl + hStep + rMargin;
      hfactor = hposr-hposl;
      hmarl = 0.0;
      hmarr = rMargin / (hposr-hposl);
    } else {
      hposl = hposr + hSpacing;
      hposr = hposl + hStep;
      hfactor = hposr-hposl;
      hmarl = 0.0;
      hmarr = 0.0;
    }
 
    for (Int_t j=0;j<Ny;j++) {
 
      if (j==0) {
	vposd = 0.0;
	vposu = bMargin + vStep;
	vfactor = vposu-vposd;
	vmard = bMargin / vfactor;
	vmaru = 0.0;
      } else if (j == Ny-1) {
	vposd = vposu + vSpacing;
	vposu = vposd + vStep + tMargin;
	vfactor = vposu-vposd;
	vmard = 0.0;
	vmaru = tMargin / (vposu-vposd);
      } else {
	vposd = vposu + vSpacing;
	vposu = vposd + vStep;
	vfactor = vposu-vposd;
	vmard = 0.0;
	vmaru = 0.0;
      }
 
      C->cd(0);
 
      char name[16];
      sprintf(name,"pad_%i_%i",i,j);
      TPad *pad = (TPad*) gROOT->FindObject(name);
      if (pad) delete pad;
      pad = new TPad(name,"",hposl,vposd,hposr,vposu);
      pad->SetLeftMargin(hmarl);
      pad->SetRightMargin(hmarr);
      pad->SetBottomMargin(vmard);
      pad->SetTopMargin(vmaru);
 
      pad->SetFrameBorderMode(0);
      pad->SetBorderMode(0);
      pad->SetBorderSize(0);
 
      pad->Draw();
    }
  }
}


void CanvasPartitionTwoColumn(TCanvas *C,const Int_t Nx,const Int_t Ny,
			      Float_t lMargin, Float_t rMargin,
			      Float_t bMargin, Float_t tMargin,
			      Float_t midMargin)
{
  if (!C) return;
 
  // Setup Pad layout:
  Float_t vSpacing = 0.0;
  Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
 
  Float_t hSpacing = midMargin;
  Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
 
  Float_t vposd,vposu,vmard,vmaru,vfactor;
  Float_t hposl,hposr,hmarl,hmarr,hfactor;
 
  Float_t midl = 0.01;
  Float_t midr = midMargin - midl;


  cout << " hStep " << hStep << " vStep " << vStep << endl;

  for (Int_t i=0;i<Nx;i++) {
 
    if (i==0) {
      hposl = 0.0;
      hposr = lMargin + hStep;
      hfactor = hposr-hposl;
      hmarl = lMargin / hfactor;
      hmarr = midl;
    } else if (i == Nx-1) {
      hposl = hposr;
      hposr = hposl + hStep + rMargin + hSpacing;
      hfactor = hposr-hposl;
      hmarl = hSpacing / hfactor;
      hmarr = rMargin / (hposr-hposl);
    } else {
      hposl = hposr + hSpacing;;
      hposr = hposl + hStep;
      hfactor = hposr-hposl;
      hmarl = 0.0;
      hmarr = 0.0;
    }
 
    for (Int_t j=0;j<Ny;j++) {
 
      if (j==0) {
	vposd = 0.0;
	vposu = bMargin + vStep;
	vfactor = vposu-vposd;
	vmard = bMargin / vfactor;
	vmaru = 0.0;
      } else if (j == Ny-1) {
	vposd = vposu + vSpacing;
	vposu = vposd + vStep + tMargin;
	vfactor = vposu-vposd;
	vmard = 0.0;
	vmaru = tMargin / (vposu-vposd);
      } else {
	vposd = vposu + vSpacing;
	vposu = vposd + vStep;
	vfactor = vposu-vposd;
	vmard = 0.0;
	vmaru = 0.0;
      }
 
      C->cd(0);
 
      char name[16];
      sprintf(name,"pad_%i_%i",i,j);
      TPad *pad = (TPad*) gROOT->FindObject(name);
      if (pad) delete pad;
      pad = new TPad(name,"",hposl,vposd,hposr,vposu);
      pad->SetLeftMargin(hmarl);
      pad->SetRightMargin(hmarr);
      pad->SetBottomMargin(vmard);
      pad->SetTopMargin(vmaru);
 
      pad->SetFrameBorderMode(1);
      pad->SetBorderMode(0);
      pad->SetBorderSize(0);
 
      pad->Draw();
    }
  }
}
