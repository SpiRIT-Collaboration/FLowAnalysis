#include "DoFlow.h"
#include "SetStyle.C"
#include "FlowFunction.C"
#include "SetColor.C"
#include "CanvasPartition.C"

struct gplot{
  TString Version;
  TString fileHeader;
  TString config1;
  TString config2;
  TString config3;
};

gplot gnames[] = {
  {".v52.15.84"        ,"finYPt_","2to5fm","|#phi|<45" ,"acr"},
  //  {".2to5fm"       ,"ut_finYPt_" ,"2to5fm" ,"ave","y_nn"},
  //  {".v52.15.85"        ,"finYPt_","2to5fm","|#phi|>135" ,"acr"},
  // {".v52.15.70.v52.15.71" ,"ut_finYPt_" ,"m55to80","ave","y_nn"},
  // {".v52.15.72.v52.15.73" ,"ut_finYPt_" ,"m50to65","ave","y_nn"},
  // {".v52.15.74.v52.15.75" ,"ut_finYPt_" ,"m40to55","ave","y_nn"},
  // {".v52.15.76.v52.15.77" ,"ut_finYPt_" ,"m35to45","ave","y_nn"},
  // {".v52.15.78.v52.15.79" ,"ut_finYPt_" ,"m20to40","ave","y_nn"},
  // {".v52.15.80.v52.15.81" ,"ut_finYPt_" ,"m0to35" ,"ave","y_nn"},
  // {".v52.15.45.v52.15.46" ,"finYPt_","m40to55","ave","ndf20"},
  // {".v52.15.47.v52.15.48" ,"finYPt_","m30to40","ave","ndf20"},
  // {".v52.15.49.v52.15.50" ,"finYPt_","m0to30","ave","ndf20"},
  //{".v52.15.51" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  //{".v52.15.55" ,"finYPt_" ,"m40to55","|#phi|<45","y_nn"},
};
const UInt_t ndata = (UInt_t)sizeof(gnames)/sizeof(gplot);
TString vfitfname = "PlotFigure.v52.15.51.v52.15.52.root";

const UInt_t npart = 5;
std::vector<UInt_t> sqPart = {4,3,2,1,0};
std::vector<UInt_t> sqSys =  {1,3,0};

std::vector<UInt_t> sqv1sel = {10, 8, 2, 1};
std::vector<UInt_t> sqv2sel = { 4, 3, 2, 1};

UInt_t ix = 0;

TString FOPI_data_sys[] = {"Au","Ru","Ca"};

Double_t FOPI_AuAu_v11x[4]={0,1,2,4};
Double_t FOPI_AuAu_v11y[4]={0.384, 0.641, 0.800, 1.032};

Double_t FOPI_AuAu_v20[4]={0.048, 0.105, 0.170, 0.247720};

TFile *outFile;
TCanvas *ccv; UInt_t iccv = 0;
TString xlabel;
TLatex  plabel;
TLatex  clabel;

Double_t *syslabel;


TGraphErrors *v1Ut;;
TGraphErrors *v2Ut;

FairLogger *logger = FairLogger::GetLogger();
void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer);
void Draw_Indiv_v1SystemD(UInt_t igname);
void Draw_Indiv_v2SystemD(UInt_t igname);
void Draw_SystemD(UInt_t igname, TString gname);

void Draw_Ut_ParticleSystemD(UInt_t igname, TString gname);
void Draw_Ut_SystemD(UInt_t igname);
void Draw_Ut_Comparison(UInt_t igname, UInt_t isys);
void Draw_Ut_Comparison_FOPI(UInt_t igname, UInt_t isys);
void Draw_Ut_Particle(Int_t igname, UInt_t isys);
void Draw_v_y(UInt_t igname, UInt_t vn = 1);
void Draw_MultiplicityD(TString gname);
void Draw_ParticleD(TString gname, UInt_t igname);
void Draw_v20_Edependence();
void Draw_Ut_Ratio(UInt_t igname, TString gname);
void Draw_Ut_RatioToOne(UInt_t igname, TString gname);

TGraphErrors* LoadAMD(UInt_t isys, UInt_t ipart, TString grname, TString eos);
TGraphErrors* LoadFOPI(UInt_t ipart, TString fdir, TString sysname );
TGraphErrors* LoadTextGraph(TString fname);


void PlotFigure()
{


  gStyle->SetOptStat(0);

  SetStyle();
  SetColor();


  // Draw function
  for( auto igname : ROOT::TSeqI( ndata ) ) {

    //--- v10 and v20 system dependece

    if( 0 ) 
      Draw_Indiv_v1SystemD(igname);
    if( 0 ) 
      Draw_Indiv_v2SystemD(igname);

     if( 0 )
      Draw_SystemD(igname, "gu_v2");
    if( 0 )
      Draw_SystemD(igname, "gu_v1");
    //--------------------

    // //edit
    if( 0 ) 
      Draw_Ut_ParticleSystemD(igname,"g_utv1_0");
    if( 0 ) 
      Draw_Ut_ParticleSystemD(igname,"g_utv2");

    if( 0 )
      Draw_Ut_Ratio(igname, "g_utv1_0");
    if( 0 )
      Draw_Ut_Ratio(igname, "g_utv2");

    if( 1 )
      Draw_Ut_RatioToOne(igname, "g_utv1_0");
    if( 1 )
      Draw_Ut_RatioToOne(igname, "g_utv2");


    if( 0 ) //OK
      Draw_Ut_SystemD(igname); // Final

    if( 0 ) 
      Draw_Ut_Comparison(igname,0); // Comparison with AMD (0:system) 

    if( 0 ) 
      Draw_Ut_Comparison_FOPI(igname, 0); // Comparison with FOPI

    // if( 0 )
    //   Draw_Ut_Particle();

    if( 0 ) {
      Draw_v_y(igname, 1);
      Draw_v_y(igname, 2);
    }
  }


  
  //---  Multiplicioty dependence
  if( 0 && ndata > 1) {
    if( outFile ) outFile->Close();
    outFile = new TFile("data/PlotFigureMD_"+gnames[0].Version+".root","recreate");
    if( 1 )
      Draw_MultiplicityD("gu_v1");
    if( 1  )
      Draw_MultiplicityD("gu_v2");
  }

  if( 0 ) {
    Draw_ParticleD("gu_v1", 0);
    Draw_ParticleD("gu_v2", 0);
    // Draw_ParticleD("gu_v1", 2);
    // Draw_ParticleD("gu_v2", 2);
  }


  if( 0 ) 
    Draw_v20_Edependence();
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

TGraphErrors* LoadAMD(UInt_t isys, UInt_t ipart, TString grname="v1_ut", TString eos="SLy4")
{
  TGraphErrors *grv = NULL;

  // TString filename[] = {"Sn132/SLy4_L108/flow_proton.root",
  // 			"Sn132/SLy4_L108/flow_deuteron.root",
  // 			"Sn132/SLy4_L108/flxoxw_triton.root"};

  TString dirname  = "../../TransportModel/AMD/Sn132/"+eos;

  TString filename[] = {"flow_proton.root",
			"flow_deuteron.root",
			"flow_triton.root",
			"",
			"flow_alpha.root"};

  TString tfile = filename[ipart];
  cout << " tfile " << tfile << endl;
  TFile *ifile = NULL;
  if( gSystem -> FindFile(dirname, tfile ) ){
    ifile = new TFile(tfile,"READ");
    if( !ifile ) return NULL;
    LOG(INFO) << tfile << " is opened. " << FairLogger::endl;
  }
  else
    return NULL;

  grv = (TGraphErrors*)ifile->Get(grname);
  if( grv == NULL ) return NULL;

  TString gname = grname + "_" + lpid[ipart];
  grv -> SetName(gname);

  delete ifile;

  return grv;
}


TGraphErrors* LoadData(UInt_t igname, UInt_t isys, UInt_t ipart, TString gname)
{
  TGraphErrors* grv = NULL;

  TFile *fOpen;
  TString fname = gnames[igname].fileHeader + bName[isys] + fpid[ipart] + gnames[igname].Version + ".root";

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

TH1I* LoadHistogram(UInt_t isys, UInt_t igname, UInt_t ipart, TString gname)
{
  TH1I* grv = NULL;

  TFile *fOpen;
  TString fname = gnames[igname].fileHeader + bName[isys] + fpid[ipart] + gnames[igname].Version + ".root";

  if( !gSystem->FindFile("data", fname) ) {
    LOG(ERROR) << fname << " is not found " << FairLogger::endl;
    return NULL;
  }
  else {
    fOpen = TFile::Open( fname );
    LOG(INFO) << fname << " is opened. " << FairLogger::endl;
  }

  grv =  (TH1I*)fOpen->Get(gname);
  if( grv == NULL ) return NULL;


  gname += "_" + rsys[isys] + "_" + lpid[ipart];
  grv -> SetName(gname);


  grv -> SetDirectory(gROOT);
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

    cout << "x = " << x << " y= " << y << endl;
    if( !std::isnan(x) ) {
      vut -> SetPoint(in, x, y);
      in++;
    }
  }
  
  if( in == 0 ) return NULL;

  vut -> RemovePoint(in-1);
  return vut;
}

Double_t *ReadVfit(TString gname, UInt_t isys, UInt_t iPart, UInt_t imult)
{
  Double_t *vfit = new Double_t[2]; 
  vfit[0] = 0.;
  vfit[1] = 0.;

  TFile* fin = TFile::Open("data/"+vfitfname);
  TIter next(fin->GetListOfKeys());
  while( TGraphErrors* obj = (TGraphErrors*)next() ) {

    if( obj->GetName() == gname ) {
      TGraphErrors *grp = (TGraphErrors*)fin->Get(obj->GetName());
      if( grp ) {
	Double_t x, xe;
	grp -> GetPoint(imult, x, vfit[0]);
	xe      = grp -> GetErrorX(imult);
	vfit[1] = grp -> GetErrorY(imult);
	return vfit;
      }
    }
  }
  
  LOG(ERROR) << "Failer to read " << gname << FairLogger::endl;
  return vfit;
}

Double_t* GetV11(UInt_t igname, UInt_t isys, UInt_t ipart)
{
  Double_t *vfit = new Double_t[2];
  auto yv1 = (TGraphErrors*)LoadData(igname, isys, ipart, "gu_v1");
  
  auto yv1_rev = new TGraphErrors();
  if( isys == 3 ) {
    for( UInt_t i = 0; i < yv1->GetN(); i++ ) {
      Double_t x = 0, y = 0;
      yv1->GetPoint( yv1->GetN()-1-i, x, y );
      yv1_rev->SetPoint( i, -x, -y );
    }

    yv1 = yv1_rev;
  }

  vfit[0] = 0.;
  vfit[1] = 0.;
  if( yv1 ) {
    yv1->Fit("fv1fit","","",v1fit[0], v1fit[1]); //"Q0","");     

    vfit[0] = fv1fit->GetParameter(1);
    vfit[1] = fv1fit->GetParError(1); 
  }

  return vfit;
}

Double_t* GetV20(UInt_t igname, UInt_t isys, UInt_t ipart)
{
  Double_t *vfit = new Double_t[2];
  vfit[0] = 0.;
  vfit[1] = 0.;

  auto yv2 = (TGraphErrors*)LoadData(igname, isys, ipart, "gu_v2");

  if( yv2 ) { 
    fv2fit->SetParameter(0,v2para0[ipart][0]);
    fv2fit->SetParameter(1,v2para0[ipart][1]);
    fv2fit->SetParameter(2,v2para0[ipart][2]);

    auto ptr = yv2->Fit("fv2fit","","",v2para0[ipart][3],v2para0[ipart][4]); //"Q0","");         


    if( fv2fit -> GetParameter(1) > 0 ) {
      vfit[0] = -fv2fit->GetParameter(0);
      vfit[1] =  fv2fit->GetParError(0); 
    }
    else {
      std::vector< Double_t > sqPara = {0.05, 0.1, 0.15, 0.2};
      for( auto ck : sqPara ) {
	fv2fit->SetParameter(0, -0.12);

	auto lpar = v2para0[ipart][3] + ck;
	auto hpar = v2para0[ipart][4] + ck;
	ptr = yv2->Fit("fv2fit","","",lpar, hpar); //"Q0","");    
	//@	yv2->Draw("ALP");

	cout << " refit v2  " << ck << ", " << lpar << " , "  << endl;
	
	if( fv2fit -> GetParameter(1) > 0 ) {
	  vfit[0] = -fv2fit->GetParameter(0);
	  vfit[1] =  fv2fit->GetParError(0); 
	  break;
	}
      }
    }
  }
  return vfit;
}

TGraphErrors* GetSDGraph(TString gname, UInt_t igname, UInt_t ipart)
{
  Double_t sysBN[]   = {82./50., 58./50.,   74./50.,62./50.,   0,     50./50.};
  Double_t sysNA[]   = {156./100.,110./100.,136./100.,136./100., 0,  156./100.};
  switch(ix) {
  case 0:
    xlabel = "(N-P)/A";
    syslabel = sysdlt;
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

  auto grp = new TGraphErrors();
  grp -> SetName(gname+Form("_%d",ipart));
  

  UInt_t isyss = 0;
  for(auto isys: sqSys ) {
    Double_t *vfit;
    if( gname == "gu_v1" ) {
      vfit = GetV11(igname,isys,ipart);
      grp -> SetTitle(";"+gnames[igname].config1+" "+xlabel+"; v10");
    }
    else {
      vfit = GetV20(igname,isys,ipart);
      grp -> SetTitle(";"+gnames[igname].config1+" "+xlabel+"; v20");
    }

    if( vfit[1] != 0. ) {
      grp -> SetPoint     (isyss, *(syslabel+isys), vfit[0] );
      grp -> SetPointError(isyss, 0,                vfit[1] );
      isyss++;
    }
  }

  if( grp->GetN() > 0 )
    return grp;
  else
    return NULL;
} 


//------------------------------------------------
//------------------------------------------------
//------------------------------------------------

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
//void Draw_
void Draw_ParticleD(TString gname, UInt_t igname)
{
  LOG(INFO) << "Draw_PartileD " << FairLogger::endl;  

  auto mv = new TMultiGraph("mv","");
  auto lg = new TLegend(0.18,0.6,0.3,0.9,gnames[igname].config1);

  Double_t y ;
  for( auto isys : sqSys ) {      
    UInt_t ipp = 0;
    UInt_t iparts = 0;
    auto grp = new TGraphErrors();
    grp -> SetName(gname+Form("_%d",isys) );
    for( auto ipart : sqPart ){

      Double_t *vfit;

      if( gname == "gu_v1") {
	vfit = GetV11(igname, isys, ipart);
	mv -> SetTitle(";Particle; v11");
      }
      else {
	vfit = GetV20(igname,isys,ipart);
	mv -> SetTitle(";Particle; v20");
      }

      // TString ggname = gname + "_" + bName[isys]+lpid[ipart];
      // vfit = ReadVfit(ggname, isys, ipart, igname);

      grp -> SetPoint     (ipp, (Float_t)ipart, vfit[0] );
      grp -> SetPointError(ipp, 0, vfit[1] );
      ipp++;

      y = vfit[0];
    }
      
    grp -> SetMarkerStyle( imark[isys] );
    grp -> SetMarkerColor( icol[isys] );
    grp -> SetLineColor( icol[isys] );
    grp -> GetXaxis() -> SetRangeUser(-1.,(Float_t)sqSys.size()+2.);      

    mv -> Add( grp, "pl");
    lg -> AddEntry(grp, lsys[isys]);
  }

  auto grp = new TGraphErrors();
  grp -> SetPoint( 0,-0.5,y+0.01);
  grp -> SetPoint( 1, 4.5,y+0.01);
  grp -> SetLineColor(10);
  mv -> Add(grp,"");
  
  for( auto k : sqPart ) {
    auto fbin = mv -> GetXaxis() -> FindBin(k);
    mv -> GetXaxis() -> SetBinLabel( fbin+1, lpid[k]);
  }	

  if( 0 ) {  // FOPI
    if( gname == "gu_v1" ) {
      auto FOPI = new TGraph(4, FOPI_AuAu_v11x, FOPI_AuAu_v11y);    
      FOPI -> SetMarkerStyle(29);
      FOPI -> SetMarkerSize(1);
      mv -> Add(FOPI,"p");
      lg -> AddEntry(FOPI,"(FOPI)Au0.25 Ut>0.8");
    }
    else {
      auto FOPI = new TGraph(4, FOPI_AuAu_v11x, FOPI_AuAu_v20);    
      FOPI -> SetMarkerStyle(29);
      FOPI -> SetMarkerSize(1);
      mv -> GetYaxis() -> SetRangeUser(0.02,0.14);
      mv -> Add(FOPI,"p");
      lg -> AddEntry(FOPI,"(FOPI)Au0.25 Ut>0.8");
    }
  }


  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  mv -> Draw("ALP");
  lg -> Draw();
}

void Draw_MultiplicityD(TString gname) 
{
  LOG(INFO) << "Draw_MultiplicityD " << FairLogger::endl;  

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;
  const Int_t Nx = 1;
  const Int_t Ny = 3;

  Float_t lMargin = 0.12;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.03;
  Float_t mMargin = 0.08;

  CanvasPartitionTwoColumn(ccv,Nx,Ny,lMargin,rMargin,bMargin,tMargin,mMargin);

  TPad *pad[Nx][Ny];

  for( auto l : ROOT::TSeqI(Nx) )for( auto m : ROOT::TSeqI(Ny) ) {
      ccv->cd(0);

      TString pname = Form("pad_%i_%i",l,m);
      pad[l][m] = (TPad*)gROOT->FindObject(pname);

      if( pad[l][m] == NULL ) {
	cout << " pad is not found " << pname << endl;
	continue;
      }
      pad[l][m] -> Draw();
      pad[l][m] -> SetFillStyle(4000);
      pad[l][m] -> SetFrameFillStyle(4000);

      Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[l][m]->GetAbsWNDC();
      Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[l][m]->GetAbsHNDC();

      auto mv = new TMultiGraph(Form("mv_%d",sqSys[l]),"");
      auto lg = new TLegend(0.15,0.05,0.3,0.25,""); 
      auto ltx= new TLatex(55,1.0, fsys[sqSys[m]]);
      ltx -> SetTextSize( ltx->GetTextSize() * yFactor );
      for( auto ipart : sqPart ) {
	UInt_t isyss = 0;
	TGraphErrors* grp = new TGraphErrors();
	grp -> SetName(gname+"_"+bName[sqSys[m]]+lpid[ipart]);

	for( auto igname : ROOT::TSeqI(ndata) ){
      
	  Double_t *vfit;
	  if( gname == "gu_v1" ) {
	    vfit = GetV11(igname,sqSys[m],ipart);
	    mv -> SetTitle(";Multiplicity; v11"); 
	  }
	  else {
	    vfit = GetV20(igname,sqSys[m],ipart);
	    mv -> SetTitle(";multiplicity; v20");
	  }

	  TH1I *hmult  = (TH1I*)LoadHistogram(sqSys[m],igname,ipart,"hmult");
	  if( hmult == NULL ) continue;
	  Double_t mmean = hmult->GetMean();
	  Double_t mstd  = hmult->GetStdDev();///sqrt(hmult->GetEntries());
	  
	  if( vfit[1] != 0. ) {
	    grp -> SetPoint     (isyss, mmean, vfit[0] );
	    grp -> SetPointError(isyss, mstd,  vfit[1] );
	    isyss++;
	  }
	}

	grp -> SetMarkerStyle(imark[ipart]);
	grp -> SetMarkerColor(pcolor[ipart]);
	grp -> SetLineColor(pcolor[ipart]);
	
	mv -> Add(grp,"p");
	lg -> AddEntry(grp, lpid[ipart] );

	if(outFile) {
	  outFile -> cd();
	  grp -> Write();
	  gROOT -> cd();
	}
	  
      }
      
      mv -> GetYaxis()->SetTitleOffset(1.0 / yFactor);
      mv -> GetYaxis()->SetTitleSize(0.06 * yFactor);
      mv -> GetXaxis()->SetNdivisions(505);

      if( gname == "gu_v1" ) 
	mv -> GetYaxis() -> SetRangeUser(0.,1.1);
      else {
	mv -> GetYaxis() -> SetRangeUser(-0.041,0.181);
	ltx -> SetY(0.15);

      }

      pad[l][m] -> cd();
      mv -> Draw("AP");
      ltx-> Draw();

      if( m == Ny-1 )
	lg -> Draw();
    }
}    

void Draw_Ut_Comparison(UInt_t igname, UInt_t isys)
{
  LOG(INFO) << "Draw_Ut_Comparison " << FairLogger::endl;

  auto mv1 = new TMultiGraph(Form("mv1_%d",isys),";U_{t0}; v1");
  auto mv2 = new TMultiGraph(Form("mv2_%d",isys),";U_{t0}; v2");
  auto lg1 = new TLegend(0.64,0.15,0.8,0.5,"");
  auto lg2 = new TLegend(0.18 ,0.15,0.4,0.55,"");

  TString utv1name[] = {"g_utv1_0","g_utv1_1"};
  TString gtitle;

  // --------- DATA _______
  if( 1 ) {

    for(auto ipart: sqPart ) {
      TString gname = utv1name[1] ;
      v1Ut = (TGraphErrors*)LoadData(igname, isys, ipart, gname);
      v1Ut -> SetName(gname);

      if( v1Ut != NULL ) {
	LOG(INFO) << " v1ut " << gname << " is registered." << FairLogger::endl;
	v1Ut -> SetMarkerStyle(20);
	v1Ut -> SetMarkerSize(1.);
	v1Ut -> SetLineColor(icol[ipart]);
	v1Ut -> SetMarkerColor(icol[ipart]);

	gtitle = v1Ut->GetTitle();
	mv1 -> Add( v1Ut, "pl");
	lg1 -> AddEntry(v1Ut,lpid[ipart]);
      }
      else
	LOG(ERROR) << " v1ut " << gname << " is not found." << FairLogger::endl;

      gname = "g_utv2";
      v2Ut = (TGraphErrors*)LoadData(igname,isys, ipart, gname);
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
  if( 0 ) {
    TString FOPI_v1Ut[][3] = { {"PRC89/Fig8_v1Ut_" , "0.25","Au"},  //[2]
			       {"PRC89/Fig8_v1Ut_" , "0.4" ,"Au"}, //[1]
			       {"NPA876/Fig9_v1Ut_", "0.4","Au"},   //[3]
			       {"NPA876/Fig11_v1Ut_","0.4","Au"},   //[4]
			       {"NPA876/Fig11_v1Ut_","0.4","Ru"},   //[5]
			       {"NPA876/Fig11_v1Ut_","0.4","Ca"}};  //[6]


    TString FOPI_v2Ut[][3] = { {"PRC89/Fig8_v2Ut_","0.25","Au"},
			       {"PRC89/Fig8_v2Ut_","0.4" ,"Au"}};


    //			       {"NPA876/Fig7_v1y_", "0.4" ,"Au"}, //[2]

    // system comparison
    UInt_t ipart = 1;
    std::vector< UInt_t > sqFOPIfig1  = {1,0};
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
  
    std::vector< UInt_t > sqFOPIfig2  = {1,0};
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
  if( 1 ) {
    for(auto ipart: sqPart ) {
      v1Ut = (TGraphErrors*)LoadAMD(isys, ipart, "v1_ut","SLy4");
      if( v1Ut != NULL ) {
	v1Ut -> SetMarkerStyle(25);
	v1Ut -> SetMarkerSize(1.);
	v1Ut -> SetLineColor(icol[ipart]);
	v1Ut -> SetMarkerColor(icol[ipart]);


	mv1 -> Add( v1Ut, "p");
	lg1 -> AddEntry(v1Ut,lpid[ipart]+"(AMD)SLy4");
      }

      v1Ut = (TGraphErrors*)LoadAMD(isys, ipart, "v1_ut","SLy4_L108");
      if( v1Ut != NULL ) {
	v1Ut -> SetMarkerStyle(26);
	v1Ut -> SetMarkerSize(1.);
	v1Ut -> SetLineColor(icol[ipart]);
	v1Ut -> SetMarkerColor(icol[ipart]);

	mv1 -> Add( v1Ut, "p");
	lg1 -> AddEntry(v1Ut,lpid[ipart]+"(AMD)SLy4_L108");
      }

      v2Ut = (TGraphErrors*)LoadAMD(isys, ipart, "v2_ut","SLy4");    
      if( v2Ut != NULL ) {
	v2Ut -> SetMarkerStyle(25);
	v2Ut -> SetMarkerSize(1.);
	v2Ut -> SetLineColor(icol[ipart]);
	v2Ut -> SetMarkerColor(icol[ipart]);


	mv2 -> Add( v2Ut, "p" );
	lg2 -> AddEntry(v2Ut,lpid[ipart]+"(AMD)SLy4");
      }
      v2Ut = (TGraphErrors*)LoadAMD(isys, ipart, "v2_ut","SLy4_L108");    
      if( v2Ut != NULL ) {
	v2Ut -> SetMarkerStyle(26);
	v2Ut -> SetMarkerSize(1.);
	v2Ut -> SetLineColor(icol[ipart]);
	v2Ut -> SetMarkerColor(icol[ipart]);


	mv2 -> Add( v2Ut, "p" );
	lg2 -> AddEntry(v2Ut,lpid[ipart]+"(AMD)L108");
      }
    }
  }
  //----

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 800, 1000); iccv++;
  ccv -> Divide(1, 2);

  ccv -> cd(1);

  mv1->GetXaxis()->SetRangeUser( 0.,  2.2);
  mv1->GetYaxis()->SetRangeUser(-0.05,0.62);
  mv1->GetXaxis()->SetNdivisions(505);
  mv1 -> Draw("AP");
  //  lg1 -> Draw();
  TLatex *v1label = new TLatex(0.3,0.62,gnames[igname].config1+" "+gtitle(0,gtitle.First(";")));
  v1label -> Draw();
  auto labelv1 = new TLatex(1.6,  0.55, fsys[isys]);
  labelv1 -> Draw();


  ccv -> cd(2);
  mv2->GetXaxis()->SetRangeUser( 0.,  2.2);
  mv2->GetYaxis()->SetRangeUser(-0.32,0.02);
  mv2->GetXaxis()->SetNdivisions(505);
  mv2 -> Draw("AP");
  lg2 -> Draw();
  TLatex *v2label = new TLatex(0.3,0.02,gnames[igname].config1+" "+" -0.4 < y_{0} < 0.4");
  v2label -> Draw();
  auto labelv2 = new TLatex(1.6, -0.02, fsys[isys]);
  labelv2 -> Draw();
}
void Draw_Ut_Comparison_FOPI(UInt_t igname, UInt_t isys)
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
    v1Ut = (TGraphErrors*)LoadData(igname, isys, ipart, gname);
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
      v1Ut = (TGraphErrors*)LoadFOPI( ipart,"NPA876/Fig7_v1Ut_0.4", FOPI_data_sys[isys]);
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

void Draw_Ut_Particle(UInt_t igname, UInt_t isys)
{

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
	  v1Ut = (TGraphErrors*)LoadData(igname, isys, sqPart[j], "g_utv1");
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
	  v2Ut = (TGraphErrors*)LoadData(igname, isys, sqPart[j], "g_utv2");

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


void Draw_Ut_SystemD(UInt_t igname)
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
      TString gtitle;

      if( i == 0 ) {

	auto mv1 = new TMultiGraph(Form("mv1_%d",sqSys[j]),";U_{t0}; v1");
	auto lg1 = new TLegend(0.75,0.3,0.9,0.55,"");
	auto labelv1 = new TLatex(1.6,  0.50, fsys[sqSys[j]]);

	for(auto ipart: sqPart ) {

	  auto v1Ut = (TGraphErrors*)LoadData(igname, sqSys[j], ipart, "g_utv1_0");
	  if( v1Ut != NULL ) {
	    
	    gtitle = v1Ut->GetTitle();
	    v1Ut -> SetMarkerStyle(20);
	    v1Ut -> SetMarkerSize(1.);
	    v1Ut -> SetLineColor(pcolor[ipart]);
	    v1Ut -> SetMarkerColor(pcolor[ipart]);
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
	  TLatex *v1label = new TLatex(0.4,0.65,gnames[igname].config1+" "+gtitle);
	  v1label -> SetTextSize(0.08);
	  v1label -> Draw();; 
	}
      }
      else {  // Right side panel

	auto mv2 = new TMultiGraph(Form("mv2_%d",sqSys[j]),";U_{t0}; v2");
	auto lg2 = new TLegend(0.25,0.3,0.4,0.55,"");
	auto labelv2 = new TLatex(1.6,  0.50, fsys[sqSys[j]]);

	for(auto ipart: sqPart ) {
	  auto v2Ut = (TGraphErrors*)LoadData(igname, sqSys[j], ipart, "g_utv2");

	  if( v2Ut != NULL ) {
	    gtitle = v2Ut->GetTitle();
	    v2Ut -> SetMarkerStyle(20);
	    v2Ut -> SetMarkerSize(1.);
	    v2Ut -> SetLineColor(pcolor[ipart]);
	    v2Ut -> SetMarkerColor(pcolor[ipart]);

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
	  TLatex *v2label = new TLatex(0.4,0.02,gnames[igname].config1+" "+gtitle);
	  v2label -> SetTextSize(0.08);
	  v2label -> Draw();; 
	}
	
      }
    }
}

void Draw_Ut_RatioToOne(UInt_t igname, TString gname)
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 800); iccv++;
  gStyle->SetOptTitle(0);

  const Int_t Ny = 2;

  Float_t lMargin = 0.02;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.05;
  CanvasPartitionY(ccv,Ny,lMargin,rMargin,bMargin,tMargin);

  TPad *pad[Ny];

  TLegend   *lgv = new TLegend(0.6, 0.05, 0.9, 0.3);
  TString gtitle;

  UInt_t ii = 0;
  std::vector<UInt_t> sqRSys = {3,0};
  for(auto isys : sqRSys) {

    auto mvUt = new TMultiGraph(Form("mvUt_%d",isys),"");
    auto tlabel = new TLatex(1.2,  1.4, fsys[isys]); 

    std::vector<UInt_t> sqPart = {0, 1, 2, 3, 4};
    for(auto ipart: sqPart ) {

    TGraphErrors *vUt = (TGraphErrors*)LoadData(igname,1,ipart,gname);
    if( vUt != NULL ) {
      cout << " vut " << gname << " is registered." << endl;

      std::vector<Double_t> dx, dy, dye;
      for(UInt_t i = 0; i < vUt->GetN(); i++){
	Double_t x, y, ye;
	vUt -> GetPoint(i,x,y);
	ye  =  vUt->GetErrorY(i);
	dx.push_back(x);
	dy.push_back(y);
	dye.push_back(ye/y);
      }

	TGraphErrors *vUtR = (TGraphErrors*)LoadData(igname,isys,ipart,gname);    
	if( vUtR == NULL ) {
	  LOG(ERROR) << gname << " is not found. " << FairLogger::endl;
	  continue;
	}


	for(UInt_t i = 0; i < vUtR->GetN(); i++){
	  Double_t nx, ny, ye;
	  vUtR -> GetPoint(i, nx, ny);
	  ye = vUtR->GetErrorY(i);

	  if( dy.at(i) != 0 ) {
	    Double_t ry = ny/dy.at(i);
	    Double_t nye = sqrt( pow(ry,2) * (pow(ye/ny,2) + pow(dye.at(i),2)) );
	    vUtR -> SetPoint(i, nx, ry);
	    vUtR -> SetPointError(i, 0., nye);
	  }
	  else
	    vUtR -> SetPoint(i, ny, 0.);
	}

	gtitle = vUt->GetTitle();
	vUtR -> SetMarkerStyle(20);
	vUtR -> SetMarkerSize(1.);
	vUtR -> SetLineColor(icol[ipart]);
	vUtR -> SetMarkerColor(icol[ipart]);

	mvUt -> Add( vUtR, "pl" );
	if( ii == Ny - 1 )
	  lgv  -> AddEntry( vUtR, lpid[ipart]);
      }
    }
    
    auto vUt = new TGraphErrors();
    vUt -> SetPoint( 0, 0. , 0.);
    vUt -> SetPoint( 1, 2.2, 0.);
    vUt -> SetLineColor(10);
    mvUt-> Add(vUt,"");

    ccv->cd(0);
    TString pname = Form("pad_%i",ii);
    pad[ii] = (TPad*)gROOT->FindObject(pname);
    pad[ii] -> Draw();
    pad[ii] -> SetFillStyle(4000);
    pad[ii] -> SetFrameFillStyle(4000);
    pad[ii] -> cd();
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ii]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ii]->GetAbsHNDC();

    mvUt -> SetTitle( gtitle );
    mvUt -> GetYaxis()->SetTitle("R(v1)");
    mvUt ->GetYaxis()->SetLabelFont(43);
    mvUt ->GetYaxis()->SetLabelSize(16);
    mvUt ->GetYaxis()->SetLabelOffset(0.02);
    mvUt ->GetYaxis()->SetTitleFont(43);
    mvUt ->GetYaxis()->SetTitleSize(18);
    mvUt ->GetYaxis()->SetTitleOffset(3.5);
    mvUt ->GetYaxis()->CenterTitle();
    mvUt ->GetYaxis()->SetNdivisions(505);
    mvUt ->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
    
    mvUt ->GetXaxis()->SetLabelFont(43);
    mvUt ->GetXaxis()->SetLabelSize(16);
    mvUt ->GetXaxis()->SetLabelOffset(0.02);
    mvUt ->GetXaxis()->SetTitleFont(43);
    mvUt ->GetXaxis()->SetTitleSize(18);
    mvUt ->GetXaxis()->SetTitleOffset(4);
    mvUt ->GetXaxis()->CenterTitle();
    mvUt ->GetXaxis()->SetTickLength(xFactor*0.04/yFactor);
    mvUt ->GetXaxis()->SetNdivisions(505);

    Double_t yrange[2][2] = {{0.5,1.5}, {0.5,1.5}};
    UInt_t iyrange = 0;
    if( gname == "g_utv2"){
      iyrange = 1;
    }

    mvUt ->GetYaxis()->SetRangeUser(yrange[iyrange][0], yrange[iyrange][1]);
      // lgv -> SetX1(1.8);
      // lgv -> SetX2(2.0);
      // lgv -> SetY1(-0.05);
      // lgv -> SetY2(-0.01);

    mvUt -> Draw("AP");
    if(ii == Ny-1) lgv  -> Draw();

    tlabel -> SetTextSize(0.07*yFactor);
    tlabel -> Draw();

    if( ii == Ny - 1 ){
      lgv  -> Draw();

      TLatex *vlabel = new TLatex(0.3,1.55,gnames[igname].config1+" "+gtitle);
      vlabel -> SetTextSize(0.08);
      vlabel -> Draw();; 
    }

    ii++;
  }
}

void Draw_Ut_Ratio(UInt_t igname, TString gname)
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 800); iccv++;
  gStyle->SetOptTitle(0);

  const Int_t Ny = 5;

  Float_t lMargin = 0.02;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.05;
  CanvasPartitionY(ccv,Ny,lMargin,rMargin,bMargin,tMargin);

  TPad *pad[Ny];

  TLegend   *lgv = new TLegend(0.8, 0.1, 0.9, 0.3);
  TString gtitle;

  UInt_t ii = 0;
  std::vector<UInt_t> sqPart = {0, 1, 2, 3, 4};
  for(auto ipart: sqPart ) {

    auto mvUt = new TMultiGraph(Form("mvUt_%d",ipart),"");
    auto tlabel = new TLatex(1.8,  0.5, lpid[ipart]); 

    TGraphErrors *vUt = (TGraphErrors*)LoadData(igname,1,ipart,gname);
    if( vUt != NULL ) {
      cout << " vut " << gname << " is registered." << endl;

      std::vector<Double_t> dx, dy, dye;
      for(UInt_t i = 0; i < vUt->GetN(); i++){
	Double_t x, y, ye;
	vUt -> GetPoint(i,x,y);
	ye  =  vUt->GetErrorY(i);
	dx.push_back(x);
	dy.push_back(y);
	dye.push_back(ye/y);
      }

      std::vector<UInt_t> sqRSys = {0,3};
      for(auto isys : sqRSys) {

	TGraphErrors *vUtR = (TGraphErrors*)LoadData(igname,isys,ipart,gname);    
	if( vUtR == NULL ) {
	  LOG(ERROR) << gname << " is not found. " << FairLogger::endl;
	  continue;
	}


	for(UInt_t i = 0; i < vUtR->GetN(); i++){
	  Double_t nx, ny, ye;
	  vUtR -> GetPoint(i, nx, ny);
	  ye = vUtR->GetErrorY(i);

	  if( dy.at(i) != 0 ) {
	    Double_t ry = ny/dy.at(i);
	    Double_t nye = sqrt( pow(ry,2) * (pow(ye/ny,2) + pow(dye.at(i),2)) );
	    vUtR -> SetPoint(i, nx, ry);
	    vUtR -> SetPointError(i, 0., nye);
	  }
	  else
	    vUtR -> SetPoint(i, ny, 0.);
	}

	gtitle = vUt->GetTitle();
	vUtR -> SetMarkerStyle(20);
	vUtR -> SetMarkerSize(1.);
	vUtR -> SetLineColor(icol[isys]);
	vUtR -> SetMarkerColor(icol[isys]);

	mvUt -> Add( vUtR, "pl" );

      }
    }
    
    vUt = new TGraphErrors();
    vUt -> SetPoint( 0, 0. , 0.);
    vUt -> SetPoint( 1, 2.2, 0.);
    vUt -> SetLineColor(10);
    mvUt-> Add(vUt,"");

    ccv->cd(0);
    TString pname = Form("pad_%i",ii);
    pad[ii] = (TPad*)gROOT->FindObject(pname);
    pad[ii] -> Draw();
    pad[ii] -> SetFillStyle(4000);
    pad[ii] -> SetFrameFillStyle(4000);
    pad[ii] -> cd();
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ii]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ii]->GetAbsHNDC();

    mvUt -> SetTitle( gtitle );
    mvUt ->GetYaxis()->SetLabelFont(43);
    mvUt ->GetYaxis()->SetLabelSize(16);
    mvUt ->GetYaxis()->SetLabelOffset(0.02);
    mvUt ->GetYaxis()->SetTitleFont(43);
    mvUt ->GetYaxis()->SetTitleSize(18);
    mvUt ->GetYaxis()->SetTitleOffset(3.5);
    mvUt ->GetYaxis()->CenterTitle();
    mvUt ->GetYaxis()->SetNdivisions(505);
    mvUt ->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
    
    mvUt ->GetXaxis()->SetLabelFont(43);
    mvUt ->GetXaxis()->SetLabelSize(16);
    mvUt ->GetXaxis()->SetLabelOffset(0.02);
    mvUt ->GetXaxis()->SetTitleFont(43);
    mvUt ->GetXaxis()->SetTitleSize(18);
    mvUt ->GetXaxis()->SetTitleOffset(4);
    mvUt ->GetXaxis()->CenterTitle();
    mvUt ->GetXaxis()->SetTickLength(xFactor*0.04/yFactor);
    mvUt ->GetXaxis()->SetNdivisions(505);

    Double_t yrange[2][2] = {{0.5,1.5}, {0.5,1.5}};
    UInt_t iyrange = 0;
    if( gname == "g_utv2"){
      iyrange = 1;
    }

    mvUt ->GetYaxis()->SetRangeUser(yrange[iyrange][0], yrange[iyrange][1]);
      lgv -> SetX1(1.8);
      lgv -> SetX2(2.0);
      lgv -> SetY1(-0.05);
      lgv -> SetY2(-0.01);

    mvUt -> Draw("AP");

    tlabel -> SetTextSize(0.07*yFactor);
    tlabel -> Draw();

    if( ii == Ny - 1 ){
      lgv  -> Draw();

      TLatex *vlabel = new TLatex(0.4,0.65,gnames[igname].config1+" "+gtitle);
      vlabel -> SetTextSize(0.08);
      vlabel -> Draw();; 
    }

    ii++;
  }
}

void Draw_Ut_ParticleSystemD(UInt_t igname, TString gname)
{

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 800); iccv++;
  gStyle->SetOptTitle(0);

  const Int_t Ny = 5;

  Float_t lMargin = 0.02;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.05;
  CanvasPartitionY(ccv,Ny,lMargin,rMargin,bMargin,tMargin);

  TPad *pad[Ny];

  TLegend   *lgv = new TLegend(0.8, 0.1, 0.9, 0.3);
  TString gtitle;

  UInt_t ii = 0;
  std::vector<UInt_t> sqPart = {0, 1, 2, 3, 4};
  for(auto ipart: sqPart ) {
    auto mvUt = new TMultiGraph(Form("mvUt_%d",ipart),"");
    auto tlabel = new TLatex(1.8,  0.5, lpid[ipart]); 
    TGraphErrors *vUt = NULL;
    
    for(auto isys : sqSys) {
      vUt = (TGraphErrors*)LoadData(igname,isys,ipart,gname);

      if( vUt != NULL ) {
	cout << " vut " << gname << " is registered." << endl;

	gtitle = vUt->GetTitle();
	vUt -> SetMarkerStyle(20);
	vUt -> SetMarkerSize(1.);
	vUt -> SetLineColor(icol[isys]);
	vUt -> SetMarkerColor(icol[isys]);

	mvUt -> Add( vUt, "pl" );
	if( ii == 0)
	  lgv  -> AddEntry( vUt,  fsys[isys], "lp");
      }
      else 
	cout << " vUt " << gname << " not be found. " << endl;
    }
    
    vUt = new TGraphErrors();
    vUt -> SetPoint( 0, 0. , 0.);
    vUt -> SetPoint( 1, 2.2, 0.);
    vUt -> SetLineColor(10);
    mvUt -> Add(vUt,"");

    ccv->cd(0);
    TString pname = Form("pad_%i",ii);
    pad[ii] = (TPad*)gROOT->FindObject(pname);
    pad[ii] -> Draw();
    pad[ii] -> SetFillStyle(4000);
    pad[ii] -> SetFrameFillStyle(4000);
    pad[ii] -> cd();
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ii]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ii]->GetAbsHNDC();

    mvUt -> SetTitle( gtitle );
    mvUt ->GetYaxis()->SetLabelFont(43);
    mvUt ->GetYaxis()->SetLabelSize(16);
    mvUt ->GetYaxis()->SetLabelOffset(0.02);
    mvUt ->GetYaxis()->SetTitleFont(43);
    mvUt ->GetYaxis()->SetTitleSize(18);
    mvUt ->GetYaxis()->SetTitleOffset(3.5);
    mvUt ->GetYaxis()->CenterTitle();
    mvUt ->GetYaxis()->SetNdivisions(505);
    mvUt ->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
    
    mvUt ->GetXaxis()->SetLabelFont(43);
    mvUt ->GetXaxis()->SetLabelSize(16);
    mvUt ->GetXaxis()->SetLabelOffset(0.02);
    mvUt ->GetXaxis()->SetTitleFont(43);
    mvUt ->GetXaxis()->SetTitleSize(18);
    mvUt ->GetXaxis()->SetTitleOffset(4);
    mvUt ->GetXaxis()->CenterTitle();
    mvUt ->GetXaxis()->SetTickLength(xFactor*0.04/yFactor);
    mvUt ->GetXaxis()->SetNdivisions(505);

    Double_t yrange[2][2] = {{-0.05,0.62}, {-0.26,0.05}};
    UInt_t iyrange = 0;
    if( gname == "g_utv2"){
      iyrange = 1;
    }

    mvUt ->GetYaxis()->SetRangeUser(yrange[iyrange][0], yrange[iyrange][1]);
      lgv -> SetX1(1.8);
      lgv -> SetX2(2.0);
      lgv -> SetY1(-0.05);
      lgv -> SetY2(-0.01);

    mvUt -> Draw("AP");

    tlabel -> SetTextSize(0.07*yFactor);
    tlabel -> Draw();

    if( ii == Ny - 1 ){
      lgv  -> Draw();

      TLatex *vlabel = new TLatex(0.4,0.65,gnames[igname].config1+" "+gtitle);
      vlabel -> SetTextSize(0.08);
      vlabel -> Draw();; 
    }

    ii++;
  }
}


void Draw_SystemD(UInt_t igname, TString gname)
{
  
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  auto mv = new TMultiGraph();
  auto lg = new TLegend(0.3,0.6,0.6,0.9,"");  

  for(auto j: sqPart ) {
    TGraphErrors* gv = GetSDGraph(gname,igname,j);
    if( gv == NULL ) continue;

    gv -> SetMarkerStyle(imark[j]);
    gv -> SetMarkerColor(pcolor[j]);
    gv -> SetLineColor(pcolor[j]);

    mv -> Add(gv);
    mv -> SetTitle(gv -> GetTitle() );
    lg -> AddEntry(gv, lpid[j]);
  }


  mv -> Draw("ALP");
  //  lg -> Draw();
}




void Draw_Indiv_v1SystemD(UInt_t igname)
{

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 840); iccv++;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const Int_t Ny = sqPart.size();

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

    //@@@
    // Format for y axis
    TGraphErrors* gv1slp = GetSDGraph("gu_v1",igname,j);
    gv1slp->Print();
    if( gv1slp == NULL ) continue;
    
    gv1slp -> SetMarkerStyle(20);
    gv1slp -> SetMarkerColor(icol[j]);
    gv1slp -> SetLineColor(icol[j]);

    gv1slp -> GetYaxis()->SetLabelFont(43);
    gv1slp -> GetYaxis()->SetLabelSize(20);
    gv1slp -> GetYaxis()->SetLabelOffset(0.04);
    gv1slp -> GetYaxis()->SetTitleFont(43);
    gv1slp -> GetYaxis()->SetTitleSize(20);
    gv1slp -> GetYaxis()->SetTitleOffset(2.5);
    gv1slp -> GetYaxis()->CenterTitle();
    gv1slp -> GetYaxis()->SetNdivisions(504);
    gv1slp -> GetYaxis()->SetTickLength(xFactor*0.04/yFactor);

    cout << j << " min " << gv1slp->GetEXlow() << " max " << gv1slp->GetEX() << endl; 
    Float_t Ymin = gv1slp->GetYaxis()->GetXmin();
    //    Ymin = Int_t(Ymin/0.01);
    Float_t Ymax = gv1slp->GetYaxis()->GetXmax();
    //    Ymax = Int_t(Ymax/0.01);
    gv1slp -> GetYaxis() -> SetRangeUser(Ymin*(1-0.003), Ymax*(1+0.003));

    // Format for x axis
    gv1slp -> GetXaxis()->SetLabelFont(43);
    gv1slp -> GetXaxis()->SetLabelSize(20);
    gv1slp -> GetXaxis()->SetLabelOffset(0.02);
    gv1slp -> GetXaxis()->SetTitleFont(43);
    gv1slp -> GetXaxis()->SetTitleSize(16);
    gv1slp -> GetXaxis()->SetTitleOffset(5);
    gv1slp -> GetXaxis()->SetNdivisions(505);
 
    // TICKS X Axis
    gv1slp -> GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
      
    gv1slp -> Draw("AP");

    cout << yFactor << endl;
    plabel.SetTextAlign(13);
    plabel.SetTextSize(0.09*yFactor);

    if( j == 0 )  plabel.DrawLatexNDC(0.25, 0.7/yFactor, fpid[j]);
    else
      plabel.DrawLatexNDC(0.25, 0.7, fpid[j]);

    //      }
  }
}

void Draw_Indiv_v2SystemD(UInt_t igname)
{
  LOG(INFO) << "Draw_Indiv_v2SystemD " << igname << FairLogger::endl;

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 850); iccv++;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const Int_t Ny = sqPart.size();;

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


    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[j]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[j]->GetAbsHNDC();


    // Format for y axis
    TGraphErrors* gv2max = GetSDGraph("gu_v2",igname,j);
    if( gv2max == NULL ) continue;

    gv2max -> SetMarkerStyle(20);
    gv2max -> SetMarkerColor(icol[j]);
    gv2max -> SetLineColor(icol[j]);
    
    gv2max -> GetYaxis()->SetLabelFont(43);
    gv2max -> GetYaxis()->SetLabelSize(20);
    gv2max -> GetYaxis()->SetLabelOffset(0.02);
    gv2max -> GetYaxis()->SetTitleFont(43);
    gv2max -> GetYaxis()->SetTitleSize(20);
    gv2max -> GetYaxis()->SetTitleOffset(3);
    gv2max -> GetYaxis()->CenterTitle();
    gv2max -> GetYaxis()->SetNdivisions(504);
    gv2max -> GetYaxis()->SetTickLength(xFactor*0.04/yFactor);

    cout << j << " min " << gv2max->GetEXlow() << " max " << gv2max->GetEX() << endl; 
    Float_t Ymin = gv2max->GetYaxis()->GetXmin();
    Ymin = Int_t(Ymin/0.001);
    Float_t Ymax = gv2max->GetYaxis()->GetXmax();
    Ymax = Int_t(Ymax/0.001);
    gv2max -> GetYaxis() -> SetRangeUser((Ymin*0.001)-0.0001, (Ymax*0.001+0.0006));

    // Format for x axis
    gv2max -> GetXaxis()->SetLabelFont(43);
    gv2max -> GetXaxis()->SetLabelSize(20);
    gv2max -> GetXaxis()->SetLabelOffset(0.02);
    gv2max -> GetXaxis()->SetTitleFont(43);
    gv2max -> GetXaxis()->SetTitleSize(16);
    gv2max -> GetXaxis()->SetTitleOffset(5);
    gv2max -> GetXaxis()->SetNdivisions(505);
 
    // TICKS X Axis
    gv2max -> GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
      
    pad[j] -> cd();
    gv2max -> Draw("AP");

    cout << yFactor << endl;
    plabel.SetTextAlign(13);
    plabel.SetTextSize(0.09*yFactor);

    if( j == 0 )  plabel.DrawLatexNDC(0.3, 0.5/yFactor, fpid[j]);
    else
      plabel.DrawLatexNDC(0.3, 0.5, fpid[j]);

    //      }
  }
}

void Draw_v_y(UInt_t igname, UInt_t vn = 1) 
{
  if( vn > 2 ) vn = 1;

  LOG(INFO) << "Draw_v" << vn << "_y"  << FairLogger::endl;

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;
  gStyle->SetOptTitle(0);

  const Int_t Nx = 1;//(Int_t)sqSys.size();
  const Int_t Ny = (Int_t)sqPart.size();

  Double_t vRange[2][5][2] ={ { {-1.0 , 1.0},
				{-0.8 , 0.8},  
				{-0.6 , 0.6 }, 
				{-0.4 , 0.4},  
				{-0.24, 0.24} }, 
			      {	{-0.10, 0.02}, 
				{-0.10, 0.02}, 
				{-0.10, 0.02}, 
				{-0.08, 0.02}, 
				{-0.05, 0.02} }};
  
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

      TLegend *lg = new TLegend(0.6,0.4/yFactor,0.9,0.65/yFactor,lpid[sqPart[j]]);

      UInt_t isys = 0;
      v_y = LoadData(igname, isys, sqPart[j], Form("gu_v%d",vn));
      LOG(INFO) << v_y -> GetName() << " is registred. " << FairLogger::endl;

      if( v_y != NULL ) {
	v_y -> SetMarkerStyle(20);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor(icol[sqPart[j]]);
	v_y -> SetMarkerColor(icol[sqPart[j]]);

	mv -> Add( v_y, "p" );
	lg -> AddEntry(v_y, "Data : "+lsys[isys]);
      }
      else
	LOG(ERROR) << " gu_v1  is not found. " << isys << " : " << sqPart[j] << FairLogger::endl;
      

      v_y = LoadAMD(sqSys[i], sqPart[j], Form("v%d_y",vn), "SLy4");
      
      if( v_y != NULL ) {
	v_y -> SetMarkerStyle(25);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor(icol[sqPart[j]]);
	v_y -> SetMarkerColor(icol[sqPart[j]]);

	mv -> Add( v_y, "pl" );
	lg -> AddEntry(v_y, "AMD-SLy4");
      }

      v_y = LoadAMD(sqSys[i], sqPart[j], Form("v%d_y",vn), "SLy4_L108");
      
      if( v_y != NULL ) {
	v_y -> SetMarkerStyle(26);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor(icol[sqPart[j]]);
	v_y -> SetMarkerColor(icol[sqPart[j]]);

	mv -> Add( v_y, "pl" );
	lg -> AddEntry(v_y, "AMD-L108");
      }

      //      gROOT->Add(mv);
      

      mv->GetYaxis()->SetRangeUser(vRange[vn-1][j][0], vRange[vn-1][j][1]);
      //      mv->Print();

      mv->Draw("AP");
      if( j == 0 )
	lg->Draw();
    }  
}


void Draw_v20_Edependence()
{
  auto v20E = (TGraphErrors*)LoadFOPI(6, "NPA876/Fig29_v2E", "Au");
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;
  ccv -> SetLogx(1);
  v20E -> Draw("ALP");

}
