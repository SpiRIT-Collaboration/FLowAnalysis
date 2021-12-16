#include "DoRPRes.C"
// #include "DoFlow.h"
#include "SetStyle.C"
//#include "FlowFunction.C"
#include "SetColor.C"
#include "CanvasPartition.C"

struct gplot{
  TString Version;
  TString fileHeader;
  TString config1;
  TString config2;
  TString config3;
};

std::vector<gplot> gnames;

std::vector<gplot> gnamesMidCentral = {
  //  {".v53.0.0" ,"finYPt_" ,"-20<#phi<30","(M3)46to55","y_nn"},
  //  {".v53.0.1" ,"finYPt_" ,"(M3)49to57","|#phi|<45","y_nn"}, //<<- best
  //  {".v53.0.2" ,"finYPt_" ,"(M3)>55","|#phi|<45","y_nn"}, 
  //  {".v53.0.3" ,"finYPt_" ,"(M3)46to55","|#phi|<45","y_nn"},     //0<<- best
  //  {".v53.0.4" ,"finYPt_" ,"1to3fm ","|#phi|<45","y_nn"},
  // {".v54.0.0" ,"finYPt_" ,"1to3fm ","|#phi|<45","y_nn"}, 
  // {".v54.0.1" ,"finYPt_" ,"1to3fm ","|#phi|<45","y_nn"}, 
  // {".v54.0.2" ,"finYPt_" ,"1to3fm ","|#phi|<45","y_nn"}, 
  // {".v54.0.3" ,"finYPt_" ,"1to3fm ","|#phi|<45","y_nn"}, 
  // {".v54.0.4" ,"finYPt_" ,"b/b0<0.15 ","|#phi|<45","y_nn"}, 
  // {".v54.0.5" ,"finYPt_" ,"b<4fm ","|#phi|<45","y_nn"},   
  // {".v54.0.6" ,"finYPt_" ,"b<1fm ","|#phi|<45","y_nn"},   
 //-----
  //  {".v54.0.7" ,"finYPt_" ,"b<1fm ","-20<#phi<30","y_nn"},
  //  {".v54.0.8" ,"finYPt_" ,"b<1fm ","|#phi|<45","y_nn"},  
 //-----
  {".v55.1.1" ,"finYPt_" ,"b<4fm "  ,"|#phi|<45"  ,"y_nn"},  
};

std::vector<gplot> gnamesCentral = {
  {".v55.1.0" ,"finYPt_" ,"b0<0.15 ","-20<#phi<30","y_nn"},   
};
//@gdata

struct mplot{
  TString path;
  TString fname;
  TString config;
  Color_t fColor;
  UInt_t  fStyle;
};

std::vector<mplot> pBUUname = {
  {"ModelData/pBUU/","Soft_s030.2L38.7"  ,"pBUU s030.2L38.7"      ,kMagenta-9, 3001},
  {"ModelData/pBUU/","Soft_s032.9L64"    ,"pBUU s032.9L64"        ,kMagenta-3, 3001},
  {"ModelData/pBUU/","Soft_s037.4L105.5" ,"pBUU s037.4L105.5"     ,kMagenta-1, 3001},
  {"ModelData/pBUU/","Stiff_s037.4L105.5","pBUU s037.4L105.5(STF)",kMagenta,   3001}
};

std::vector<mplot> AMDnameText = {
  //  {"E270-124-112-sly4-l055-mdcorr20-cxAg","rapflow-","sly4-l055-md20-cxAg" , kCyan-7,   3001},
  //{"E270-124-112-sly4-l055-mdcorr50-cxAg","rapflow-","sly4-l055-md50-cxAg" , kOrange+8, 3001},
  //{"E270-124-112-sly4-mdcorr50-cxAq"     ,"rapflow-","sly4-l046-md50-cxAq" , kOrange-2, 3001},
  //{"Sn124Sn112-sly4-l108-mdcorr50-cxAq"  ,"rapflow-","sly4-l108-md50-cxAq" , kGreen+2,  3001},
  //{"Sn124Sn112-sly4-l152-mdcorr50-cxAq"  ,"rapflow-","sly4-l152-md50-cxAq" , kViolet-1, 3001},
  //
  //{"Sn132Sn124-sly4-l152-mdcorr50-cxAq"  ,"sly4-l152-md50-cxAq" , kViolet-1, 3001},
  //  {"Sn132Sn124-sly4-mdcorr50-cxAq"       ,"sly4-l046-md50-cxAq" , kOrange-2, 3001}
  {"Sn124Sn112-sly4-l46-mdcorr50-cxAg"   ,"rapflow-","sly4-l46-md50-cxAg" , kViolet-1, 3001}
};

std::vector<mplot> AMDname202111 = {
  //  {"20211111.cxAg/","amd_Sn124Sn112_SLy4_L46_mdcorr20_cxAq","SLy4_L46_md20_cxAq", kCyan-8, 3001},
  {"20211111.cxAg/","amd_Sn124Sn112_SLy4_L46_mdcorr50_cxAq","SLy4_L46_md50_cxAq", kOrange+7, 3001},
  {"20211111.cxAg/","amd_Sn124Sn112_SLy4_L46_mdcorr50_cxAq.06sgm","SLy4_L46_md50_cxAq.06s", kViolet-1, 3001},
  {"20211111.cxAg/","amd_Sn124Sn112_SLy4_cxAq.06sgm","SLy4_L46_cxAq.06s", kOrange-1, 3001},
  {"20211111.cxAg/","amd_Sn124Sn112_SLy4_cxAq.09sgm","SLy4_L46_cxAq.09s", kGreen-1, 3001},
};

std::vector<mplot> AMDnameCentral = {
  {"bimp0_1.5_highStat/","amd_Sn108Sn112_SLy4_b0.15",      "SLy4 b0.15"      ,kMagenta, 3001},
  {"bimp0_1.5_highStat/","amd_Sn132Sn124_SLy4_b0.15",      "SLy4_b0.15"      ,kMagenta, 3001},
  {"bimp0_1.5_highStat/","amd_Sn132Sn124_SLy4_sgm12_b0.15","SLy4_sgm12_b0.15",kGreen+2, 3001},
  {"bimp0_1.5_highStat/","amd_Sn108Sn112_SLy4_sgm12_b0.15","SLy4_sgm12_b0.15",kGreen+2, 3001},
  // {"ModelData/AMD/bimp0_1.5_highStat/","h_amd_SLy4_sigmul2_Sn108",         "SLy4_sigmul2"         ,kViolet-1,  3001},
  // {"ModelData/AMD/bimp0_1.5_highStat/","h_amd_SLy4_sigmul2_Sn132",         "SLy4_sigmul2"         ,kViolet-1,  3001},
  // {"h_amd_SLy4_gfg_sigmul2_Sn108",     "SLy4_gfg_sigmul2"     ,kMagenta-3, 3001},
  // {"h_amd_SLy4_gfg_sigmul2_Sn108_ss",  "SLy4_gfg_sigmul2_ss"  ,kMagenta-3, 3001},
  // {"h_amd_SLy4_gfg_sigmul2_Sn132",     "SLy4_gfg_sigmul2"     ,kMagenta-3, 3001},
  // {"h_amd_SLy4_L108_Sn108",            "SLy4_L108"            ,kOrange+8,  3001},
  // {"h_amd_SLy4_L108_Sn132_ss",         "SLy4_L108_ss"         ,kOrange+8,  3001},
  // {"h_amd_SLy4_L108_sigmul2_Sn108",    "SLy4_L108_sigmul2"    ,kViolet-1,  3001},  
  // {"h_amd_SLy4_L108_sigmul2_Sn132",    "SLy4_L108_sigmul2"    ,kViolet-1,  3001},  
  // {"h_amd_SLy4_L108_gfg_sigmul2_Sn108","SLy4_L108_gfg_sigmul2",kGreen+2,   3001},  
  // {"h_amd_SLy4_L108_gfg_sigmul2_Sn132","SLy4_L108_gfg_sigmul2",kGreen+2,   3001}    
};

// AMD data selection
Bool_t bCentral = 0;
std::vector<mplot>  AMDnames;

Bool_t bAMDRoot = 1;
Bool_t bReverse = 0;

const UInt_t ndata = gnames.size();
TString vfitfname = "PlotFigure.v52.15.51.v52.15.52.root";


Double_t ycut[5]= {2.0, 0.5, 0.4, 0.3, 0.2};

const UInt_t npart = 5;
std::vector<UInt_t> sqPart = {4,3,2,1,0};
std::vector<UInt_t> sqSys =  {0,1,3};

std::vector<UInt_t> sqv1sel = {10, 8, 2, 1};
std::vector<UInt_t> sqv2sel = { 4, 3, 2, 1};

UInt_t ix = 0;
Double_t sysBN[]   = {82./50., 58./50.,   74./50.,62./50.,   0,     50./50.};
Double_t sysNA[]   = {156./100.,110./100.,136./100.,136./100., 0,  156./100.};

TString FOPI_data_sys[] = {"Au","Ru","Ca"};

Double_t FOPI_AuAu_v11x[4]={0,1,2,4};
Double_t FOPI_AuAu_v11y[4]={0.384, 0.641, 0.800, 1.032};

Double_t FOPI_AuAu_v20[4]={0.048, 0.105, 0.170, 0.247720};

TFile *outFile;
TCanvas *ccv; UInt_t iccv = 0;
TLatex  plabel;
TLatex  clabel;

TString xlabel;


TGraphErrors *v1Ut;;
TGraphErrors *v2Ut;

FairLogger *logger = FairLogger::GetLogger();
void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer);
void Draw_Indiv_v1SystemD(UInt_t igname);
void Draw_Indiv_v2SystemD(UInt_t igname);
void Draw_DeltaDOne(UInt_t igname, TString gname);
void Draw_DeltaDIndiv(UInt_t igname, TString gname);

void Draw_Ut_ParticleSystemD(UInt_t igname, TString gname);
void Draw_Ut_SystemD(UInt_t igname);
void Draw_Ut_Comparison(UInt_t igname, UInt_t isys);
void Draw_Ut_Comparison_FOPI(UInt_t igname, UInt_t isys);
void Draw_Ut_Particle(Int_t igname, UInt_t isys);
void Draw_v_y(UInt_t igname, UInt_t vn = 1);
void Draw_v_y_System(UInt_t igname, UInt_t vn=1, TString para="y"); 
void Draw_MultiplicityD(TString gname);
void Draw_MultiplicityRatio(TString gname);
void Draw_ParticleD(TString gname, UInt_t igname);
void Draw_v20_Edependence();
void Draw_Ut_Ratio(UInt_t igname, TString gname);
void Draw_Ut_RatioToOne(UInt_t igname, TString gname);
void Draw_Ut_ProtonRatio(UInt_t igname, TString gname);
void Draw_Ut_ProtonRatioSystemD(UInt_t igname, TString gname);
void Draw_RPResolutionDM();
Double_t* GetRPResolution(UInt_t igname, UInt_t isys);
void Draw_v_y_Model(UInt_t igname, UInt_t isys, UInt_t vn, TString para);
void Draw_dndy(UInt_t isys, UInt_t ipart, TString opt="");
void Draw_dndyRatio(UInt_t isys, UInt_t ipart, TString gname);
void Draw_Ratio(UInt_t igname, TString gname);

TGraphErrors* LoadTextGraph(TString fname);


//------------------------------------------------
//------------------------------------------------
//------------------------------------------------
//###LoadData()
TObject* ReversePlot(TObject* obj, TString opt="")
{
  LOG(INFO) << " ReversePlot is called. " << obj->GetName() << FairLogger::endl;

  if( obj->InheritsFrom("TH1D") ) {
    TH1D* h_nrm = (TH1D*)obj->Clone();
    TH1D* h_rev = (TH1D*)obj->Clone();
    for( auto ibin: ROOT::TSeqI( h_nrm->GetNbinsX() ) ){
      Double_t x,y,ye;
      x  = h_nrm->GetBinCenter(ibin);
      y  = h_nrm->GetBinContent(ibin);
      ye = h_nrm->GetBinError(ibin);

      auto rbin = h_rev->GetXaxis()->FindBin(-x);
      h_rev->SetBinContent( rbin, y );
      h_rev->SetBinError(   rbin, ye );
    }
    
    delete h_nrm;
    return h_rev;
  }
  else if( obj->InheritsFrom("TGraphErrors") ) {
    TGraphErrors* g_nrm = (TGraphErrors*)obj->Clone();
    TGraphErrors* g_rev = new TGraphErrors();
    for( auto nx : ROOT::TSeqI(g_nrm->GetN()) ) {
      Double_t x, y, xe, ye;
      g_nrm -> GetPoint( g_nrm->GetN()-nx-1, x, y);
      xe = g_nrm -> GetErrorX( g_nrm->GetN()-nx-1);
      ye = g_nrm -> GetErrorY( g_nrm->GetN()-nx-1);
      
      if( opt == "v1" )
	g_rev->SetPoint(nx, -x, -y);
      else
	g_rev->SetPoint(nx, -x, y);
      g_rev->SetPointError(nx, xe, ye);
    }
    delete g_nrm;
    return g_rev;
  }

  return NULL;
}


TObject* LoadData(UInt_t isys, UInt_t ipart, TString grname, gplot gname, Int_t ylim = -1.)
{
  TGraphErrors* grv = NULL;
  TH1D*         hst = NULL;

  TFile *fOpen;
  TString fname = gname.fileHeader + bName[isys] + fpid[ipart] + gname.Version + ".root";

  cout << fname << endl;
  if( !gSystem->FindFile("data", fname) ) {
    LOG(ERROR) << fname << " is not found " << FairLogger::endl;
    return NULL;
  }
  else {
    fOpen = TFile::Open( fname );
    LOG(INFO) << fname << " is opened. " << FairLogger::endl;
  }

  TH1I* hmult = (TH1I*)fOpen->Get("hmult");
  if( !hmult )
    return NULL;
  Double_t total = (Double_t)hmult->GetEntries();

  cout << " LoadData: "<< grname << endl;

  if( grname.Contains("dndy") ) {
    hst = new TH1D();

    if( grname == "g_dndy" || grname == "h_dndy" ) {
      TH2D* hacpcorr = (TH2D*)fOpen->Get("hyutacpcrr");
      if( hacpcorr == NULL ) return NULL;
      cout << hacpcorr->GetName() << endl;

      hst -> SetName("h_dndy_"+rsys[isys]+"_"+lpid[ipart]+"_DATA");
    
      hst = hacpcorr->ProjectionX(hst->GetName());
      hst -> Scale(1./( total * hst->GetXaxis()->GetBinWidth(0)));
      hst -> SetDirectory(0);

      fOpen -> Close();

      if( isys == 2 )
	hst = (TH1D*)ReversePlot((TH1D*)hst);

      if( grname.Contains("h_") ) {
	return hst;
      }
      else  if( grname.Contains("g_") ) {
	grv = new TGraphErrors();
	grname += "_" + rsys[isys] + "_" + lpid[ipart];
	grv -> SetName(grname);
	
	UInt_t iibin = 0;
	for( auto ibin: ROOT::TSeqI(hst->GetXaxis()->GetNbins())) {
	  if( hst->GetBinError(ibin) != 0. ) {
	    grv -> SetPoint(iibin, hst->GetXaxis()->GetBinCenter(ibin), 
			    hst->GetBinContent(ibin));
	    grv -> SetPointError(iibin, 0., hst->GetBinError(ibin));
	    iibin++;
	  }
	}
	return grv;
      }
      return NULL;
    }
    else if( grname == "h_dndyx" || grname == "g_dndyx")  {
      TH2D* hacpcorr = (TH2D*)fOpen->Get("hyyxacpcrr");
      if( hacpcorr == NULL ) return NULL;
      cout << hacpcorr->GetName() << endl;

      hst -> Reset();
      hst -> SetName("h_dndyx_"+rsys[isys]+"_"+lpid[ipart]+"_DATA");

      Double_t lbin = hacpcorr->GetXaxis()->FindBin(-1.*ycut[ylim]);
      Double_t ubin = hacpcorr->GetXaxis()->FindBin(ycut[ylim]);
      cout << ycut[ylim] << " lbin " << lbin << " ubin " << ubin << endl;
      hst = hacpcorr->ProjectionY(hst->GetName(),lbin,ubin);
      cout << hst->GetName() << endl;
      Double_t rnorm = 2.*ycut[ylim] * hst->GetXaxis()->GetBinWidth(0);
      hst -> Scale(1./( rnorm * total));
      hst -> SetDirectory(0);

      fOpen->Close();
      
      if( grname.Contains("h_") )
	return hst;
      else if( grname.Contains("g_") ) {
	grv = new TGraphErrors();
	grname += "_" + rsys[isys] + "_" + lpid[ipart];
	grv -> SetName(grname);
	UInt_t iibin = 0;
	for( auto ibin: ROOT::TSeqI(hst->GetXaxis()->GetNbins())) {
	  if( hst->GetBinError(ibin) != 0. ) {
	    grv -> SetPoint(iibin, hst->GetXaxis()->GetBinCenter(ibin), 
			    hst->GetBinContent(ibin));
	    grv -> SetPointError(iibin, 0., hst->GetBinError(ibin));
	    iibin++;
	  }
	}
	return grv;
      }
      return NULL;
    }
  } 
  else if( grname.Contains("g") ){
    grv =  (TGraphErrors*)fOpen->Get(grname);
    grname += "_" + rsys[isys] + "_" + lpid[ipart]+"_DATA";
    grv -> SetName(grname);
    if( grv == NULL ) {
      fOpen->Close();
      return NULL;
    }
    fOpen->Close();
    return grv;
  }
  else if( grname.Contains("h_") ) {
    hst = (TH1D*)fOpen->Get(grname);
    if( hst == NULL ) return NULL;
    
    grname += "_" + rsys[isys] + "_" + lpid[ipart]+"_DATA";
    hst -> SetName(grname);
    hst -> SetDirectory(0);
    fOpen->Close();
    
    if( isys == 2 )
      hst = (TH1D*)ReversePlot((TH1D*)hst);
      
    return hst;
  }
  return NULL;
}

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


TGraphErrors* LoadAMDText(UInt_t isys, UInt_t ipart, mplot mdata, TString grname="v1")
{
  TString column[] = {"y","h_dndy","h_dndye","tmp1","tmp1e","tmp2","tmp2e","g_v1y","g_v1ye","g_v2y","g_v2ye"};
  TString pname[] = {"p", "d", "t", "h", "a", "n"};

  auto  gvt = new TGraphErrors();


  TString dirname = "ModelData/AMD/" + mdata.path;
  TString ifile = mdata.fname + pname[ipart] + ".dat"; 

  LOG(INFO) << dirname << "/" << ifile  << FairLogger::endl;

  if( !gSystem -> FindFile(dirname, ifile ) ) {
    LOG(ERROR) << ifile << " is not found. " << FairLogger::endl;
    return NULL;
  }
  else {
    LOG(INFO) << ifile << " is accessed. " << FairLogger::endl;
    
    std::fstream fread;
    fread.open(ifile, std::fstream::in);

    TString sget;
    UInt_t in = 0;
    
    while( !fread.eof() ) {
      
      for( auto ict : ROOT::TSeqI(3*5) )
	fread >> sget;

      while( !fread.eof() ) {

	Double_t x, y; 
	Double_t ye = 0.;
	for( auto iss : column ) {

	  fread >> sget;
	  if( iss == "y" ) 
	    x = (Double_t)atof(sget);

	  if( iss == grname )
	    y = (Double_t)atof(sget);

	  if( iss == grname+"e" )
	    ye = (Double_t)atof(sget);
	}
	//	cout  << grname << " " << in << " "  << x << " vs " << y << " +- " << ye << endl;
      
	///@@@@
      
	if( ye != 0. || y != 0 ) {
	  x = x/y_cm[10];

	  if( grname == "h_dndy" )
	    y *= y_cm[10]; 


	  gvt -> SetPoint(in, x, y);
	  if( ye == 0 ) ye = y*0.05;
	  gvt -> SetPointError(in, 0, ye);
	  in++;
	}
      }
    }

    return gvt;
  }

  return NULL;

}


//#LoadAMDRoot
TObject* LoadAMDRoot(UInt_t isys, UInt_t ipart, mplot mdata, TString grname="h_dndy")
{
  TString dirname  = "ModelData/AMD/" + mdata.path;  
  TString tfile = mdata.fname + ".root";
  LOG(INFO) << "AMD " << tfile << " is selected." << FairLogger::endl;

  if( !tfile.Contains(Form("%d",Asys[isys])) ) return NULL;

  TFile *ifile = NULL;
  if( gSystem -> FindFile(dirname, tfile ) ) {

    ifile = new TFile(tfile,"READ");
    if( !ifile ) return NULL;
    LOG(INFO) << tfile << " is opened. " << FairLogger::endl;
  }
  else{
    LOG(ERROR) << tfile << " is not found." << FairLogger::endl; 
    return NULL;
  }
  
  if( !grname.Contains("g_vartl") )
    grname += "_"+lpid[ipart];
  cout << grname << endl;


  TH1D* hist = NULL;
  TGraphErrors* grv = NULL;

  if( grname.Contains("dndy") ) {
    hist = (TH1D*)ifile->Get("h_"+grname(2, grname.Length()));
    if( hist == NULL ) return NULL;
  }
  else { 
    grv = (TGraphErrors*)ifile->Get(grname);
    if( grv == NULL ) return NULL;
  }

  if( grname.Contains("h_") && hist ) {
    if( bReverse && !grname.Contains("dndyx") ) 
      hist = (TH1D*)ReversePlot((TH1D*)hist);

    return hist;
  }
  else if( hist && grname.Contains("g_" )){
    grv = new TGraphErrors();
    UInt_t iibin = 0;
    for( auto ibin: ROOT::TSeqI(hist->GetXaxis()->GetNbins())) {
      if( hist->GetBinError(ibin) != 0. ) {
	grv -> SetPoint(iibin, hist->GetXaxis()->GetBinCenter(ibin), 
			hist->GetBinContent(ibin));
	grv -> SetPointError(iibin, 0., hist->GetBinError(ibin));
	iibin++;
      }
    }

    if( bReverse && !grname.Contains("g_vartl") ) {
      auto grv_rev = new TGraphErrors();
      for( auto nx : ROOT::TSeqI(grv->GetN()) ) {
	Double_t x, y, xe, ye;
	grv -> GetPoint( grv->GetN()-nx-1, x, y);
	xe = grv -> GetErrorX( grv->GetN()-nx-1);
	ye = grv -> GetErrorY( grv->GetN()-nx-1);
	
	if( grname.Contains("v1") )
	  grv_rev->SetPoint(nx, -x, -y);
	else
	  grv_rev->SetPoint(nx, -x, y);
	grv_rev->SetPointError(nx, xe, ye);
      }
      grv = grv_rev;
    }
    grv -> SetName(grname+"_AMD"+mdata.config);
    return grv;
  }

  return grv;
}

//###LoadAMD
TObject* LoadAMD(UInt_t isys, UInt_t ipart, mplot mdata, TString grname="v1")
{

  if( isys == 3 && (mdata.fname.Contains("Sn124Sn112") || mdata.path.Contains("Sn124Sn112")) ) bReverse = 1;
  else if( !mdata.fname.Contains(asys[isys]) ) return NULL;

  if( bAMDRoot )
    return (TObject*)LoadAMDRoot(isys, ipart, mdata, grname);
  else
    return (TObject*)LoadAMDText(isys, ipart, mdata, grname);
}


TGraphErrors* LoadpBUU(UInt_t isys, UInt_t ipart, mplot mdata, TString grname="v1_Ut0")
{
  TGraphErrors  *gvt = NULL;

  if( ipart == 4 ) return gvt;
  
  TString dirname = "ModelData/pBUU";
  TString ifile = grname + "_pd3Het_" + rsys[isys] + mdata.fname + ".txt";
  TString infile = ifile;

  if( !gSystem -> FindFile(dirname, infile ) ) {
    LOG(ERROR) << ifile << " is not found. " << FairLogger::endl;
    return NULL;
  }
  else {
    LOG(INFO) << infile << " is accessed. " << FairLogger::endl;
    

    std::fstream fread;
    fread.open(infile, std::fstream::in);
    gvt = new TGraphErrors();

    Double_t x, y[4], ye[4];
    TString sget;
    UInt_t  in = 0;
    
    while( !fread.eof() ) {
      fread >> sget;
      if( sget != "t" ) continue;

      while( !fread.eof() ) {
	fread >> sget;

	x = (Double_t)atof(sget);
      
	for( auto i : {0,1,3,2} ) {
	  fread >> sget;
	  y[i] = (Double_t)atof(sget);
	  fread >> sget;
	  ye[i] = (Double_t)atof(sget);
	}
	///@@@@
      
	if( ye[ipart] != 0. ) {
	  gvt -> SetPoint(in, x, y[ipart]);
	  gvt -> SetPointError(in, 0, ye[ipart]);

	  in++;
	}
      }
    }
  }

  return gvt;
}
TGraphErrors* LoadImQMD(UInt_t isys, UInt_t ipart, TString grname="v1", TString eos="Soft") 
{
  TString partname[] = {"p","d","t","He3","He4"};

  TGraphErrors  *gvt = new TGraphErrors();
  
  TString ifile = "imqmd.root";
  TString dirname = "ModelData/ImQMD";

      
  if( !gSystem -> FindFile(dirname, ifile ) ) {
    LOG(ERROR) << ifile << " is not found. " << FairLogger::endl;
    return NULL;
  }
  else {
    LOG(INFO) << ifile << " is accessed. " << FairLogger::endl;
  }

  TFile *fread = TFile::Open( ifile );

  TString gname = partname[ipart]+"_"+grname+"_"+"sn"+rsys[isys];
  LOG(INFO) << gname << " will be opened. " << FairLogger::endl;
  TH1F *hget = (TH1F*)fread->Get(gname);
  if( hget ) 
    LOG(INFO) << hget->GetName() << " is opened. " << FairLogger::endl;
  else
    return NULL;

  UInt_t in = 0;
  for( auto ibin : ROOT::TSeqI( hget->GetNbinsX() ) ) {
    
    auto x = hget->GetBinCenter(ibin);
    auto y = hget->GetBinContent(ibin);
    auto ye= hget->GetBinError(ibin);

    if( ye != 0. && ye < 0.5) {

      if( grname == "v1_Pt" )
	x = x/mass[ipart]/u_p;

      gvt -> SetPoint(in, x, y);
      gvt -> SetPointError(in, 0, ye);
      
      in++;
    }
  }


  return gvt;
}


TGraphErrors* LoadKanekoData(UInt_t isys, UInt_t ipart, TString gname="h1dndy")
{  
  TString partname[] = {"Proton","Deuteron","Triton","He3","He4"};

  auto *ofile = TFile::Open("data/dndy_Kaneko.root");
  
  if( ofile == NULL ) {
    LOG(ERROR) << "data/dndy_Kaneko.root is not found " << FairLogger::endl;
    return NULL;
  }


  TString hname = gname + "_" + bName[isys] + partname[ipart];
  TH1D* hist = (TH1D*)ofile->Get(hname);

  if( hist == NULL )
    return NULL;

  cout << " LoadKanekoData " << hname << " is found."<< endl;

  auto gvr = new TGraphErrors();
  Double_t x, y, ye;

  UInt_t iin = 0;
  for( auto ibin : ROOT::TSeqI(hist->GetNbinsX())) {
    x = hist->GetBinCenter(ibin);
    y = hist->GetBinContent(ibin);
    ye= hist->GetBinError(ibin);

    if( ye != 0 ) {
      gvr -> SetPoint(iin, x, y);
      gvr -> SetPointError(iin, 0, ye);
      iin++;
    }
  }
  
  return gvr;
}

TGraphErrors* LoadTommyData(UInt_t isys, UInt_t ipart, TString gname)
{
  TString partname[] = {"pim","pip","p","d","t","He3","He4"};
  TString sysname[]  = {"_Sn132","_Sn108","_Sn124","_Sn112"};

  TGraphErrors* grv = NULL;
  
  TString fname = "Content_Mult_" + gname  + sysname[isys] + ".txt";
  TString ffind = fname;

  if( !gSystem->FindFile("data/TommyPlots", ffind ) ) {
    LOG(ERROR) << fname << " is not found." << FairLogger::endl;
    return NULL;
  }

  
  LOG(INFO) << ffind << " is found." << FairLogger::endl;
  
  std::fstream fread;
  fread.open(ffind, std::fstream::in);
  
  grv = new TGraphErrors();
  
  Double_t x, y, ye;
  TString sget;
  UInt_t iin = 0;
  
  TString search = partname[ipart+2]+"_"+gname + sysname[isys]; 
  cout << " seraching for " << search << endl;

  Bool_t beof = kTRUE;
  while( beof ) {
    fread >> sget;
    if( sget == search ) {
      fread >> sget; fread >> sget; fread >> sget;

      while( !fread.eof() ) {
	fread >> sget;
	if( sget.Contains( gname ) || fread.eof()) {
	  beof = kFALSE;
	  break;
	}

	x = (Double_t)atof(sget);
	fread >> sget;
	y = (Double_t)atof(sget);
	fread >> sget;
	ye= (Double_t)atof(sget);

	if( ye != 0. ) {
	  grv -> SetPoint(iin, x, y);
	  grv -> SetPointError(iin, 0, ye);
	  iin++;
	}
      }
    }
  }


  return grv;
}


Double_t GetTommyDataIntegral(UInt_t isys, UInt_t ipart, TString gname="rapHist_M40_52")
{
  TString partname[] = {"pim","pip","p","d","t","He3","He4"};
  TString sysname[]  = {"_Sn132","_Sn108","_Sn124","_Sn112"};

  
  TString fname = "Content_Mult_" + gname  + sysname[isys] + ".txt";
  TString ffind = fname;

  if( !gSystem->FindFile("data/TommyPlots", ffind ) ) {
    LOG(ERROR) << fname << " is not found." << FairLogger::endl;
    return 0.;
  }

  
  LOG(INFO) << ffind << " is found." << FairLogger::endl;
  
  std::fstream fread;
  fread.open(ffind, std::fstream::in);
  
  auto hist = new TH1D(Form("hist_%d",ipart),"",100,-2.,2.);
  auto grv  = new TGraphErrors();
  
  Double_t x, y, ye;
  TString sget;
  UInt_t iin = 0;
  
  TString search = partname[ipart+2]+"_"+gname + sysname[isys]; 
  cout << " seraching for " << search << endl;

  Bool_t beof = kTRUE;
  while( beof ) {
    fread >> sget;
    if( sget == search ) {
      fread >> sget; fread >> sget; fread >> sget;

      while( !fread.eof() ) {
	fread >> sget;
	if( sget.Contains( gname ) || fread.eof()) {
	  beof = kFALSE;
	  break;
	}

	x = (Double_t)atof(sget);
	fread >> sget;
	y = (Double_t)atof(sget);
	fread >> sget;
	ye= (Double_t)atof(sget);

	if( ye != 0. ) {
	  hist -> Fill(x,y);
	  grv -> SetPoint(iin, x, y);
	  grv -> SetPointError(iin, 0, ye);
	  iin++;
	}
      }
    }
  }


  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  ccv -> Divide(1,2);
  ccv -> cd(1);
  grv -> Draw("ALP");
  ccv -> cd(2);
  hist -> Draw();

  return hist -> Integral();

}
 
Double_t* GetRPResolution(UInt_t igname, UInt_t isys)
{

  TFile *fOpen;
  TString fname = gnames[igname].fileHeader + bName[isys] + fpid[0] + gnames[igname].Version + ".root";

  if( !gSystem->FindFile("data", fname) ) {
    LOG(ERROR) << fname << " is not found " << FairLogger::endl;
    return NULL;
  }
  else {
    fOpen = TFile::Open( fname );
    LOG(INFO) << fname << " is opened. " << FairLogger::endl;
  }

  auto res = new Double_t[4];

  TH1D* hcos0  = (TH1D*)fOpen->Get("hdphi0_180");
  TH1D* hcos90 = (TH1D*)fOpen->Get("hdphi90_180");
  GetRPResolutionwChi(res, hcos0, hcos90, 1);

  TH1I* mult = (TH1I*)fOpen->Get("hmult");
  *(res+2) = mult->GetMean();
  *(res+3) = mult->GetStdDev();

  fOpen->Close();

  return res;
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

Double_t* GetV11(UInt_t igname, UInt_t isys, UInt_t ipart, TString sData, mplot mdata)
{
  gStyle->SetOptFit(1111);
  gStyle->SetStatColor(10);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.5);
  fv1fit->GetChisquare();

  Double_t *vfit = new Double_t[8];
  for( auto i : ROOT::TSeqI(8) )
    vfit[i] = 0.;

  TGraphErrors* yv1;
  Double_t ft_low  = -0.3;
  Double_t ft_high =  0.9;
  if( sData == "pBUU" ) {
    yv1 = LoadpBUU(isys, ipart, mdata, "v1_yn");
    ft_low  = -0.7;
    ft_high =  0.7;
  }
  else if( sData == "ImQMD" ) {
    yv1 = LoadImQMD(isys, ipart, "v1", "");
    ft_low  = -0.8;
    ft_high =  0.8;
  }
  else if( sData == "AMD"){
    yv1 = (TGraphErrors*)LoadAMD(isys, ipart, mdata, "v1");

    ft_low = -0.8;
    ft_high = 0.8;
  }
  else
    yv1 = (TGraphErrors*)LoadData(isys, ipart, "gu_v1",gnames[igname]);

  if( yv1 == NULL )
    return vfit;
  
  if( sData=="DATA" && kFALSE ) {  // reverse
    //  if( isys == 1 || isys == 3 ) {
    auto yv1_rev = new TGraphErrors();
    for( UInt_t i = 0; i < yv1->GetN(); i++ ) {
      Double_t x = 0, y = 0;
      yv1->GetPoint( yv1->GetN()-1-i, x, y );
      yv1_rev->SetPoint( i, -x, -y );

    }
    yv1 = yv1_rev;
  }
  

  if( yv1 ) {
    fv1fit -> SetLineColor(2);
    fv1fit -> SetParameter(1, v1para[ipart]);

    ft_low  = -0.5;
    ft_high =  0.5;
	
    yv1->Fit("fv1fit","QM","", ft_low, ft_high); //"Q0","");     
	
    vfit[0] = fv1fit->GetParameter(1);
    vfit[1] = fv1fit->GetParError(1); 
    vfit[2] = fv1fit->GetParameter(2);
    vfit[3] = fv1fit->GetParError(2);
    vfit[4] = fv1fit->GetChisquare();
    vfit[5] = fv1fit->GetNDF();
  

    if( sData != "DATA" ) {
      TCanvas *ccm =(TCanvas*)gROOT->FindObject("ccm_"+sData);
      if( !ccm ) {
	ccm = new TCanvas("ccm_"+sData,"ccm_"+sData,1000,600);//Form("ccm_%d-%d",isys,ipart),"v1fit_"+lsys[isys]+";"+lpid[ipart]);
	ccm -> Divide(2,5);
      }

      UInt_t idp = isys == 1 ? 2*(isys+ipart) : 2*(isys+ipart)+1 ;
      ccm->cd( idp );

      yv1 -> SetTitle(lsys[isys]+"_"+lpid[ipart]+";y/y_{nn}-1;v1");
      yv1 -> SetMarkerStyle(20);
      yv1 -> SetMarkerColor(3);
      yv1 -> SetLineColor(0);
      yv1 -> Draw("ALP");
      TLatex tlabel;
      tlabel.DrawLatex(0., yv1->GetYaxis()->GetXmax(), lsys[isys]+"_"+lpid[ipart]);    
      tlabel.DrawLatex(-0.4, 0.2, Form("%4.2f to %4.2f",ft_low, ft_high));
      
    }
  }

  return vfit;
}

Double_t* GetV20(UInt_t igname, UInt_t isys, UInt_t ipart, TString sData, mplot mdata)
{
  Double_t *vfit = new Double_t[2];
  vfit[0] = 0.;
  vfit[1] = 0.;

  auto yv2 = (TGraphErrors*)LoadData(isys, ipart, "gu_v2", gnames[igname]);

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

TGraphErrors* GetSDGraph(TString gname, UInt_t igname, UInt_t ipart, TString sData="Data", 
			 mplot mdata=AMDnames[0], UInt_t vn=1)
{
  Double_t *syslabel;

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


  UInt_t inn = 0;
  for( auto isys : {0, 3, 1} ){
    Double_t *fitval;

    if( vn == 2 )
      fitval = GetV20(igname, isys, ipart, sData, mdata);
    else
      fitval = GetV11(igname, isys, ipart, sData, mdata);

    if( *fitval == 0 ) continue;

    UInt_t ipara = 0;
    if( vn==1 && gname.Contains("_v13") )
	ipara += 2;
	
    grp -> SetPoint( inn, *(sysdlt+isys), *(fitval+ipara) );
    grp -> SetPointError( inn, 0, *(fitval+ipara+1) );
    inn++;
  }

  return grp;
}

void Draw_DeltaDIndiv(UInt_t igname, TString gname) 
{
  Bool_t bImQMD = 0;
  Bool_t bpBUU = 0;
  Bool_t bAMD  = 1;

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 840); iccv++;
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
  
    auto mv = new TMultiGraph();
    auto lg = new TLegend(0.3,0.3,0.6,0.6,"");  
    lg -> SetFillColor(0);

    if( gname.Contains("_v2") ) {
      mv -> SetTitle(gnames[igname].config1+"  "+gnames[igname].config2+
		     ";"+xlabel+"; v20");
    }
    else if( gname.Contains("_v13") && !gname.Contains("_v13/v11") ) {
      mv -> SetTitle(gnames[igname].config1+"  "+gnames[igname].config2+
		     ";"+xlabel+"; v13");
    }
    else if( gname.Contains("_v13/v11") ) {
      mv -> GetXaxis()->CenterTitle();
    }
    else {
      mv -> SetTitle(gnames[igname].config1+"  "+gnames[igname].config2+   
		     ";"+xlabel+"; v11"); 
    }
  
    TGraphErrors *grp;
    
    // PBUU    
    if( bpBUU ) {
      for( auto i : {1} ) { //ROOT::TSeqI(sizeof(pBUUname)/sizeof(mplot)) ) {
	grp = GetSDGraph(gname,igname,j, "pBUU",pBUUname[i]);
	if( grp == NULL ) continue;
      
	grp -> SetLineWidth(20);
	grp -> SetLineColor(CStyle[4+i].fColor);
	grp -> SetMarkerStyle(0);
	mv -> Add(grp, "3");
	lg -> AddEntry(grp, pBUUname[i].config);
      }
    }

    // ImQMD
    if( bImQMD ) {
      grp = GetSDGraph(gname,igname,j, "ImQMD");
      if( grp ) {
	grp -> SetLineWidth(1);
	grp -> SetLineColor(CStyle[6].fColor);
	grp -> SetFillColor(CStyle[6].fColor);
	grp -> SetMarkerStyle(0);
	mv -> Add(grp, "3");
	lg -> AddEntry(grp, "ImQMD");
      }
    }

    if( bAMD ) {
      for( auto samd : AMDnames ) {
	grp = GetSDGraph(gname, igname, j, "AMD", samd);
	if( grp ) {
	  grp -> SetLineWidth( 1 );
	  grp -> SetLineColor( samd.fColor);
	  grp -> SetFillColor( samd.fColor);
	  mv -> Add(grp, "3");
	  lg -> AddEntry(grp, "AMD"+ samd.config);
	}
      }
    }


    grp = GetSDGraph(gname,igname,j, "DATA");
    if( grp == NULL ) continue;
    Double_t ymean = grp -> GetMean(2);

    grp -> SetMarkerStyle(20);
    grp -> SetMarkerSize(1);
    grp -> SetMarkerColor(2);
    grp -> SetLineColor(2);

    //    grp -> Print();

    mv -> Add(grp, "P");
    lg -> AddEntry(grp, "Data");

    

    Double_t alpha = 0.4;

    if( gname.Contains("v2") )    alpha = 0.2;

    Double_t ylow  = Double_t( Int_t((1.-alpha)*ymean*1000.) / 1000. );
    Double_t yhigh = Double_t( Int_t((1.+alpha)*ymean*1000.) / 1000. );
    // mv -> SetMaximum( yhigh );
    // mv -> SetMinimum( ylow );

    // mv -> GetYaxis()->SetLabelFont(43);
    // mv -> GetYaxis()->SetLabelSize(20);
    // mv -> GetYaxis()->SetLabelOffset(0.04);
    // mv -> GetYaxis()->SetTitleFont(43);
    // mv -> GetYaxis()->SetTitleSize(20);
    // mv -> GetYaxis()->SetTitleOffset(2.5);
    // mv -> GetYaxis()->CenterTitle();
    // mv -> GetYaxis()->SetNdivisions(504);
    // mv -> GetYaxis()->SetTickLength(xFactor*0.04/yFactor);


    pad[j] -> cd(0);
    mv -> Draw("ALP");

    Double_t gx,gy;
    grp -> GetPoint(0, gx, gy);
    auto Line108 = new TLine(mv->GetXaxis()->GetXmin(), gy, 
			     mv->GetXaxis()->GetXmax(), gy);
    Line108 -> SetLineColor(2);
    Line108 -> SetLineStyle(3);
    //    Line108 -> Draw();
    //    Line108 -> Print();

    
    if( j == 0 ) {
      if( gname.Contains("v2") ) {
	lg -> SetX1(0.2);
	lg -> SetX2(0.45);
      }
      
      lg -> SetFillColor(0);
      lg -> SetTextSize(0.08*yFactor);
      lg -> Draw();
    }

    TLatex tlabel;
    tlabel.SetTextSize(0.08*yFactor);
      
    auto labely = (1.+alpha)*ymean-0.05*(1.+alpha)*ymean;
    
    tlabel.DrawLatex(0.2, labely, lpid[j]);

  }
}


void Draw_DeltaDOne(UInt_t igname, TString gname)
{  
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  auto mv = new TMultiGraph();
  auto lg = new TLegend(0.3,0.7,0.6,0.9,"");  
  lg -> SetFillColor(0);

  for(auto j: {0,1,2,3,4} ) {
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
  if( gname.Contains("v2") ) {
    lg -> SetX1(0.7);
    lg -> SetX2(0.9);
  }

  lg -> Draw();
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

    plabel.SetTextAlign(13);
    plabel.SetTextSize(0.09*yFactor);

    if( j == 0 )  plabel.DrawLatexNDC(0.25, 0.7/yFactor, lpid[j]);
    else
      plabel.DrawLatexNDC(0.25, 0.7, lpid[j]);

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

void Draw_v_y_System(UInt_t igname, UInt_t vn=1, TString para="y") 
{
  if( vn > 2 ) vn = 1;
  LOG(INFO) << "Draw_v" << vn << "_y_System"  << FairLogger::endl;

  Bool_t bData = 1;
  Bool_t bAMD = 0;
  Bool_t bImQMD = 0;
  Bool_t bpBUU = 0;

  TGraphErrors* v_y;
  TString ypara[4];
  TString header;
  Double_t xaxis[2];
  Double_t yaxis[2];

  for(auto ipart : {0,1,2,3,4} ) {

    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 600, 500); iccv++;
    gStyle -> SetOptTitle(0);

    auto mv = new TMultiGraph();
    mv -> SetName(Form("mv_v%d",vn));
    mv -> SetTitle(";"+para+Form("; v%d", vn) );
    auto lg = new TLegend(0.2, 0.2, 0.6, 0.5, lpid[ipart]);
    lg -> SetFillColor(0);

    if( para == "y" ) {
      ypara[0] = Form("gu_v%d",vn);
      xaxis[0] = -0.5;
      xaxis[1] =  0.5;
      yaxis[0] = -0.5;
      yaxis[1] =  0.5;

      if( vn == 2 ){
	yaxis[0] = -0.1;
	yaxis[1] =  0.02;
      }

    }
    else if( vn == 1 ) { 
      xaxis[0] =-0.01;
      xaxis[1] = 2.;
      yaxis[0] =-0.7;
      yaxis[1] = 0.7;

      if( ipart == 0 || ipart == 1 ) {
	yaxis[0] = - 0.1;
	yaxis[1] =   0.3;
      }      
      ypara[0] =  "g_utv1_0";
    }
    else if( vn == 2 ) {
      ypara[0] = "g_utv2";
    }
    
    for( auto isys : {1,3,0} ) {

      if( vn == 2 ){
	lg -> SetY1(0.15);
	lg -> SetY2(0.4);
      }

    
      v_y = (TGraphErrors*)LoadData(isys, ipart, ypara[0], gnames[igname]);
      if( v_y != NULL ) 
	LOG(INFO) << v_y -> GetName() << " is registred. " << FairLogger::endl;
      else {
	LOG(ERROR) << "Failed to get " << ypara[0] << FairLogger::endl;
	continue;
      }


      v_y -> SetMarkerStyle(imark[isys]);
      v_y -> SetMarkerSize(1.2);
      v_y -> SetLineColor(icol[isys]);
      v_y -> SetMarkerColor(icol[isys]);

      if( para == "y" ) {
	auto f1 = v_y -> GetFunction(Form("fv%dfit",vn));
	if( f1 != NULL )
	  f1 -> SetLineColor(icol[ipart]);
	// else {
	//   v_y -> Fit(Form("fv%dfit",vn),"","",v1fit[0], v1fit[1]);
	//   f1 = v_y -> GetFunction(Form("fv%dfit",vn));
	  
	//   f1 -> SetLineColor(icol[ipart]); 
	// }
      }


      header = v_y -> GetTitle();
      header = header(0,header.First(";"));
      lg -> SetHeader(lpid[ipart]+" "+header);

      if( bData ) {
	mv -> Add( v_y, "p" );
	lg -> AddEntry(v_y, lsys[isys]); 
      }

      if( bpBUU ) {
	
	for( auto i : {0} ){ //ROOT::TSeqI(sizeof(pBUUname)/sizeof(mplot)) ) {
	  if( para == "y" )
	    ypara[1] = Form("v%d_yn",vn);
	  else 
	    ypara[1] = Form("v%d_Ut0",vn);

	  v_y = LoadpBUU(isys, ipart, pBUUname[i], ypara[1]);

	
	  if( v_y == NULL ) continue;
      
	  v_y -> SetMarkerStyle( CStyle[4+i].mStyle );
	  v_y -> SetMarkerSize(1.);
	  v_y -> SetLineColor( icol[isys] );
	  v_y -> SetMarkerColor( icol[isys] );
	
	  mv -> Add( v_y, "p" );
	  lg -> AddEntry(v_y, lsys[isys]+pBUUname[i].config);
	}
      }

      if( bImQMD ) {
      
	if( para == "y" )
	  ypara[2] = Form("v%d",vn);
	else
	  ypara[2] = Form("v%d_Pt",vn);

	v_y = LoadImQMD(isys, ipart, ypara[2] );
      
	if( v_y != NULL ){
	  v_y -> SetMarkerStyle(CStyle[5].mStyle);
	  v_y -> SetMarkerSize(1.);
	  v_y -> SetLineColor( icol[isys] );
	  v_y -> SetMarkerColor( icol[isys] );
	
	  mv -> Add( v_y, "p" );
	  lg -> AddEntry(v_y, lsys[isys]+"ImQMD");
	}
      }
    }

    mv -> GetXaxis()->SetRangeUser(xaxis[0], xaxis[1]);
    mv -> GetYaxis()->SetRangeUser(yaxis[0], yaxis[1]);
    mv -> GetXaxis()->SetNdivisions(505);
    mv -> GetYaxis()->SetNdivisions(505);
    mv -> Draw("AP");
    lg -> Draw();


    if( para == "y" ) {
      auto Ymin = mv->GetYaxis()->GetXmin();
      auto Ymax = mv->GetYaxis()->GetXmax();
      auto Xmin = mv->GetXaxis()->GetXmin();
      auto Xmax = mv->GetXaxis()->GetXmax();
    
      auto aLineX1 = new TLine(Xmin, 0., Xmax, 0.);
      aLineX1->SetLineColor(1);
      aLineX1->SetLineStyle(3);
      aLineX1->Draw();
    
      auto aLineY1 = new TLine(0., Ymin, 0., Ymax);
      aLineY1->SetLineColor(1);
      aLineY1->SetLineStyle(3);
      aLineY1->Draw();
    }
  }
}

void Draw_v_y(UInt_t igname, UInt_t vn = 1) 
{
  if( vn > 2 ) vn = 1;

  Bool_t bAMD = 1;
  Bool_t bpBUU = 0;
  Bool_t bImQMD = 0;

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
      lg -> SetFillColor(0);

      UInt_t isys = 0;
      v_y = (TGraphErrors*)LoadData(isys, sqPart[j], Form("gu_v%d",vn), gnames[igname]);
      LOG(INFO) << v_y -> GetName() << " is registred. " << FairLogger::endl;

      if( v_y != NULL ) {
	v_y -> SetMarkerStyle(CStyle[0].mStyle);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor( CStyle[0].fColor );
	v_y -> SetMarkerColor( CStyle[0].fColor );

	mv -> Add( v_y, "p" );
	lg -> AddEntry(v_y, "Data : "+lsys[isys]);
      }
      else
	LOG(ERROR) << " gu_v1  is not found. " << isys << " : " << sqPart[j] << FairLogger::endl;
      
      if( bAMD ) {
	v_y = (TGraphErrors*)LoadAMD(sqSys[i], sqPart[j], AMDnames[0], Form("v%d_y",vn));
	
	if( v_y != NULL ) {
	  v_y -> SetMarkerStyle(CStyle[2].mStyle);
	  v_y -> SetMarkerSize(1.);
	  v_y -> SetLineColor(CStyle[2].fColor);
	  v_y -> SetMarkerColor(CStyle[2].fColor);
	
	  mv -> Add( v_y, "pl" );
	  lg -> AddEntry(v_y, "AMD-SLy4");
	}
      }

      if( bpBUU ) {
	
	for( auto ipbuu : pBUUname ) {
	  if( vn == 1 )
	    v_y = LoadpBUU(sqSys[i], sqPart[j], ipbuu , "v1_yn");
	  if( vn == 2 )
	    v_y = LoadpBUU(sqSys[i], sqPart[j], ipbuu , "v2_yn");
	  if( v_y == NULL ) continue;
      
	  v_y -> SetMarkerStyle(CStyle[4+i].mStyle);
	  v_y -> SetMarkerSize(1.);
	  v_y -> SetLineColor(ipbuu.fColor);
	  v_y -> SetMarkerColor(ipbuu.fColor);
	  
	  mv -> Add( v_y, "pl" );
	  lg -> AddEntry(v_y, ipbuu.config);
	}
      }

      if( bImQMD ) {
	if( vn == 1)
	  v_y = LoadImQMD(sqSys[i], sqPart[j], "v1" );
	else if( vn == 2)
	  v_y = LoadImQMD(sqSys[i], sqPart[j], "v2" );
      if( v_y != NULL ){
	v_y -> SetMarkerStyle(CStyle[6].mStyle);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor(CStyle[6].fColor);
	v_y -> SetMarkerColor(CStyle[6].fColor);

	mv -> Add( v_y, "pl" );
	lg -> AddEntry(v_y, "ImQMD");
	
      }


      //      mv->GetYaxis()->SetRangeUser(vRange[vn-1][j][0], vRange[vn-1][j][1]);
      //      mv->Print();

      mv->Draw("AP");
      if( j == 0 )
	lg->Draw();
      }
    }
}

void Draw_v_y_Model(UInt_t igname, UInt_t isys, UInt_t vn=1, TString para="y_{cm}") 
{}

void Draw_v_y_Model(UInt_t isys, UInt_t ipart, std::vector<gplot> gdata, UInt_t vn=1, TString para="y_{cm}") 
{
  if( vn > 2) vn = 1;

  Bool_t bAMD = 1;
  Bool_t bImQMD = 0;
  Bool_t bpBUU = 0;


  TGraphErrors* v_y;
  TString ypara;
  TString header;


  TMultiGraph *mv = new TMultiGraph();
  mv -> SetName(Form("mv_%d",ipart));
  mv -> SetTitle(";"+para+Form("; v%d",vn));

  TLegend *lg = new TLegend(0.2,0.5,0.6,0.9,lpid[ipart]);
  lg -> SetFillColor(0);

  if( para.Contains("y") )
    ypara = Form("gu_v%d",vn);
  else if( vn == 1 )
    ypara =  "g_utv1_0";
  else if( vn == 2 )
    ypara = "g_utv2";
	     

  // DATA
  std::vector< UInt_t > ivsys;
  ivsys.push_back(isys);
  if( isys == 3 ) 
    ivsys.push_back(2);

  UInt_t i = 0;
  for( auto gname : gdata )for(auto iisys : ivsys) {
    
      v_y = (TGraphErrors*)LoadData(iisys, ipart, ypara, gname); 
      if( v_y != NULL ) 
	LOG(INFO) << v_y -> GetName() << " is registred. " << FairLogger::endl;
      else {
	LOG(ERROR) << "Failed to get " << ypara << FairLogger::endl;
	continue;
      }

      header = v_y -> GetTitle();
      header = lpid[ipart]+":"+header(0,header.First(";")); 


      v_y -> SetMarkerStyle(CStyle[i].mStyle);
      v_y -> SetMarkerSize(1.);
      v_y -> SetLineColor( CStyle[i].fColor );
      v_y -> SetMarkerColor( CStyle[i].fColor );
    
      mv -> Add( v_y, "pl" );
      lg -> AddEntry(v_y, "Data : "+fsys[iisys]+" "+gname.config1);

      i++;
    }
  //---------

  if( bAMD && para.Contains("y") ) {
    
    for( auto samd : AMDnames ) {
      v_y = (TGraphErrors*)LoadAMD(isys, ipart, samd, Form("g_v%dy",vn));
      if( v_y != NULL ) {

	v_y -> SetFillStyle( samd.fStyle );
	v_y -> SetLineColor( samd.fColor );
	v_y -> SetFillColorAlpha( samd.fColor, 0.4 );
	
	mv -> Add( v_y, "3" );

	lg -> AddEntry(v_y, "AMD_"+samd.config);
      }
    }
  }

  if( bpBUU ) {
	
    UInt_t i = 0;
    for( auto spbuu : pBUUname ) {
      if( para == "y" )
	ypara = Form("v%d_yn",vn);
      else 
	ypara = Form("v%d_Ut0",vn);

      //      v_y = LoadpBUU(isys, ipart, spbuu, ypara);
      if( v_y == NULL ) continue;
      
      v_y -> SetMarkerStyle( CStyle[4+i].mStyle );
      v_y -> SetMarkerSize(1.);
      v_y -> SetLineColor( CStyle[4+i].fColor );
      v_y -> SetMarkerColor( CStyle[4+i].fColor );
      
      mv -> Add( v_y, "p" );
      lg -> AddEntry(v_y, spbuu.config);

      i++;
    }
  }

  if( bImQMD ) {
      
    if( para == "y" )
      ypara = Form("v%d",vn);
    else
      ypara = Form("v%d_Pt",vn);
    
    v_y = LoadImQMD(isys, ipart, ypara );
    
    if( v_y != NULL ){
      v_y -> SetMarkerStyle(CStyle[5].mStyle);
      v_y -> SetMarkerSize(1.);
      v_y -> SetLineColor( CStyle[5].fColor );
      v_y -> SetMarkerColor( CStyle[5].fColor );
      
      mv -> Add( v_y, "p" );
      lg -> AddEntry(v_y, "ImQMD");
    }
  }

  if( para.Contains("y") ) {
    if( vn == 1 ) {
      mv -> GetYaxis() -> SetRangeUser(-0.6, 0.6);
      mv -> GetXaxis() -> SetRangeUser(-1.2 , 1.2);
    }
    else if( vn == 2 ) {
      mv -> GetYaxis() -> SetRangeUser(-0.15, 0.15);
      mv -> GetXaxis() -> SetRangeUser(-1.2 , 1.2);
    }
  }
  else {
    if( vn == 1 ) {
      mv -> GetYaxis() -> SetRangeUser(-0.2, 0.8);
      mv -> GetXaxis() -> SetRangeUser( 0. , 2.);
    }
    else if( vn == 2 ) {
      mv -> GetYaxis() -> SetRangeUser(-0.4, 0.05);
      mv -> GetXaxis() -> SetRangeUser( 0.   ,2.);
      lg -> SetHeader( header );
    }
  }


  if(ipart == 0 ) {
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
    lg  -> Draw();
  }


  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 600); iccv++;
  gStyle->SetOptTitle(0);
  ccv->SetGrid(1);
  
  mv->Draw("AP");
    
  plabel.SetTextAlign(13);
  plabel.DrawLatexNDC(0.15,0.9, fsys[isys]+"_"+lpid[ipart]);


}

void Draw_RPResolutionDM()
{
  LOG(INFO) << "Draw_RPResolutionDM " << FairLogger::endl;  

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;
  ccv->Divide(1,3);

  auto mv = new TMultiGraph("mdep",";Multiplicity; <cos(#Psi_{M}-#Psi_{R})>");
  auto lg = new TLegend(0.55,0.15,0.9,0.45,""); 
  lg -> SetFillColor(0);


  Double_t mmean[4];
  Double_t msigm[4];

  UInt_t id = 1;
  TGraphErrors *mres;
  for( auto isys : {0,3,1}) {
    ccv->cd(id); id++;
    mres = new TGraphErrors();
    mres -> SetName(Form("mres_%d",isys));

    for( auto igname : ROOT::TSeqI(ndata) ) {
      Double_t *res = GetRPResolution(igname, isys );
      
      // cout << " m " << *(res+2) << " +- " << *(res+3)
      // 	   << " cos " << *res << " +- " << *(res+1) << endl;

      mres -> SetPoint(igname, *(res+2), *res);
      mres -> SetPointError(igname, *(res+3), *(res+1));
    }

    mres -> SetMarkerStyle(20);
    mres -> SetMarkerColor(icol[isys]);
    mres -> SetLineColor(icol[isys]);
    mv -> Add(mres, "p");
    lg -> AddEntry(mres, fsys[isys]);

    mres -> Draw("ALP");
    mres -> Fit("gaus");
    
    auto ffit = mres -> GetFunction("gaus");
    ffit -> SetLineColor(icol[isys]);

    mmean[isys] = ffit -> GetParameter(1);
    msigm[isys] = ffit -> GetParameter(2);

  }

  // Float_t lMargin = 0.12;
  // Float_t rMargin = 0.05;
  // Float_t bMargin = 0.10;
  // Float_t tMargin = 0.03;
  // Float_t mMargin = 0.08;
  
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  mv -> Draw("AP");
  lg -> Draw();

  Double_t fct = 0.2;
  for( auto isys : {0,1,3})
    LOG(INFO) << fsys[isys] << " m : " << mmean[isys] << " ->> " 
	      << mmean[isys] - fct*msigm[isys]
	      << " to " 
	      << mmean[isys] + fct*msigm[isys] << FairLogger::endl;
  
}

void Draw_v20_Edependence()
{
  auto v20E = (TGraphErrors*)LoadFOPI(6, "NPA876/Fig29_v2E", "Au");
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;
  ccv -> SetLogx(1);
  v20E -> Draw("ALP");

}


Double_t* Getvartxz(UInt_t isys, UInt_t ipart, Int_t ylim=-1) 
{
  gStyle -> SetLegendFillColor(0);
  
  Double_t *vart = new Double_t[2];
  THStack*  mhst = new THStack("hs","");
  
  auto lg  = new TLegend(0.15, 0.2, 0.5, 0.7,"");
  lg -> SetFillStyle(0);

  TH1D *h_dndy = NULL;
  TH1D *h_dndyx = NULL;

  for( auto gdata : {gnames[0]} ){
    h_dndy  = (TH1D*)LoadData(isys,ipart,"h_dndy" ,gdata,-1); 
    h_dndyx = (TH1D*)LoadData(isys,ipart,"h_dndyx",gdata, ylim); 

    if( h_dndy != NULL && h_dndyx != NULL) {

      Double_t vv[2];
      Double_t vve[2];
      TH1D* hists[] = {h_dndy,h_dndyx};     
      mhst -> SetTitle(";y_{cm}/y_{prj};dN/dy");

      UInt_t i = 0;
      for( auto hst : hists ) {

	TH1D* h_dndy_lim = (TH1D*)hst->Clone();
	auto x1 = h_dndy_lim->GetXaxis()->FindBin(-1.);
	auto x2 = h_dndy_lim->GetXaxis()->FindBin( 1.);
	auto mid= h_dndy_lim->GetXaxis()->FindBin( 0.);
	
	UInt_t nbin = h_dndy_lim->GetXaxis()->GetNbins();
	for( auto ibin : ROOT::TSeqI( nbin ) ) {
	  if( hst == hists[0] && ibin > mid ) {
	    auto rap  = h_dndy_lim -> GetBinCenter(ibin);
	    auto cont = h_dndy_lim -> GetBinContent(ibin);
	    auto conte= h_dndy_lim -> GetBinError(ibin);
	    auto rbin = h_dndy_lim -> GetXaxis()->FindBin( -rap );
	    h_dndy_lim->SetBinContent( rbin, cont );
	    h_dndy_lim->SetBinError( rbin, conte );
	  }

	  // if ( ibin < x1 || ibin > x2 )
	  //   h_dndy_lim->SetBinContent(ibin,0.);
	}

	h_dndy_lim -> SetMarkerStyle(CStyle[i].mStyle);
	h_dndy_lim -> SetMarkerSize(1.);
	h_dndy_lim -> SetLineColor( CStyle[i].fColor );
	h_dndy_lim -> SetMarkerColor( CStyle[i].fColor );
	h_dndy_lim -> SetMarkerSize(1.);

	mhst -> Add( h_dndy_lim);
	//	mhst -> Add( hst);
	lg  -> AddEntry( h_dndy_lim, fsys[i]+" "+gdata.config1);

	vv[i] = h_dndy_lim->GetStdDev();
	i++;
      }

      vart[0] = vv[1]*vv[1] / (vv[0]*vv[0]);
      vart[1] = 2.* abs(vve[1]/vv[1]) +  2.* abs(vve[0]/vv[0]);
      LOG(INFO) << " vartxz = : " << vart[0] << " +- " << vart[1]  << FairLogger::endl;
    }
      
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
    ccv -> SetGrid(1);
    //    mhst -> GetXaxis() -> SetRangeUser(-2.5, 2.5);
    mhst -> Draw("nostack");
  }
  
  return vart;
}

TGraphErrors* Get_gvartxz(UInt_t isys, Int_t ylim=-1)
{
  auto g_vartxz = new TGraphErrors();
  Double_t* vart;

  UInt_t in = 0;
  for( UInt_t ipart : {0,1,2,3,4} ){
    vart = Getvartxz(isys, ipart, ylim);
    
    if( vart[0] != 0 ) {
      g_vartxz->SetPoint( in, ipart+1, vart[0]);
      g_vartxz->SetPointError( in, 0,  vart[1]);
      in++;
    }
  }
  
  if( g_vartxz->GetN() == 0 ) return NULL;
  
  for(auto ipart : {0,1,2,3,4} )
    g_vartxz->GetXaxis()->SetBinLabel( g_vartxz->GetXaxis()->FindBin( ipart+1 ), lpid[ipart] );


  return g_vartxz;
}

void Draw_vartxz(UInt_t isys, Int_t ylim=2)
{
  Bool_t bAMD = 1;

  auto m_vartl = new TMultiGraph();
  auto legend = new TLegend(0.18,0.15,0.55,0.4,"");

  TGraphErrors* g_vartl;

  g_vartl = Get_gvartxz(isys,3);
  g_vartl-> SetMarkerStyle(20);
  g_vartl-> SetMarkerColor(2);
  g_vartl-> SetLineColor(2);

  if( g_vartl ) {
    m_vartl->Add(g_vartl,"alp");
    g_vartl->GetYaxis()->SetRangeUser(0.3,1.15);
    legend ->AddEntry(g_vartl,"DATA"+fsys[isys]);

  if( bAMD ){
    for( auto samd : AMDnames ) {
      g_vartl = (TGraphErrors*)LoadAMD(isys,0,samd, Form("g_vartl1_y%d",(UInt_t)(ycut[ylim]*10)));
      if( g_vartl ) {
	g_vartl -> SetFillStyle( samd.fStyle );
	g_vartl -> SetLineColor( samd.fColor );
	g_vartl -> SetFillColorAlpha( samd.fColor, 0.4 );
	m_vartl->Add(g_vartl,"3");
	legend ->AddEntry(g_vartl,samd.config);
      }
    }
  }

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  auto gtmp = new TGraph();
  gtmp -> SetPoint(0,0,0.2);
  gtmp -> SetPoint(1,6,0.2);
  m_vartl->Add(gtmp,"");
  for(auto ipart : {0,1,2,3,4} ) {
    auto fbin = m_vartl -> GetXaxis() -> FindBin(ipart+1);
    if( fbin != 101 )
      m_vartl->GetXaxis()->SetBinLabel( fbin, lpid[ipart] );
  }
  m_vartl -> Draw("LP");
  legend -> Draw();
  
  }
}
  
void Draw_dndy(UInt_t isys, UInt_t ipart, std::vector<gplot> gname)
{
  Bool_t bAMD  = 1;
  Bool_t bImQMD = 0;
  Bool_t bpBUU = 0;
  Bool_t bTommy = 1;
  Bool_t bKaneko = 1;
  if( !bCentral ) bKaneko = 0;

  auto mgr = new TMultiGraph();
  mgr -> SetTitle(";y_{cm}/y_{prj};dN/dy");

  auto lg  = new TLegend(0.15, 0.2, 0.5, 0.7,"");
  lg -> SetHeader(lpid[ipart]);
  lg -> SetFillStyle(0);

  TGraphErrors *g_dndy = NULL;

  std::vector< UInt_t > ivsys;
  ivsys.push_back(isys);
  if( isys == 3 ) 
    ivsys.push_back(2);

  UInt_t i = 0;
  for( auto igname : gname ) {
    for( auto iisys : ivsys ) {
      g_dndy =(TGraphErrors*)LoadData(iisys,ipart,"g_dndy",igname); 
      if( g_dndy != NULL ) {
	g_dndy -> SetMarkerStyle(CStyle[i].mStyle);
	g_dndy -> SetMarkerSize(1.);
	g_dndy -> SetLineColor( CStyle[i].fColor );
	g_dndy -> SetMarkerColor( CStyle[i].fColor );
	g_dndy -> SetMarkerSize(1.);
	
	if( iisys == 2 ) {
	  g_dndy -> SetMarkerStyle(25);
	  g_dndy -> SetMarkerSize(1.);
	  g_dndy -> SetLineColor( CStyle[i].fColor );
	  g_dndy -> SetMarkerColor( CStyle[i].fColor );
	  g_dndy -> SetMarkerSize(1.);
	}

	mgr -> Add( g_dndy, "P");
	lg  -> AddEntry( g_dndy, fsys[iisys]+" "+igname.config1+" "+igname.Version);
      }
      i++;
    }
  }  

  if( bTommy ) {
    TString opt = "";
    if( bCentral ) opt = "M55_55";
    else if( !bCentral ) {
      if( isys == 0 ) opt = "M46_46";
      else if( isys == 1 ) opt = "M43_52";
      else if( isys == 3 ) opt = "M46_52";
    }


    UInt_t iclsys[] = {2,2,2,2};
    TString glbl[]  = {"Data","Data","^{124}Sn(Reverse)","Data"};
    

    for( auto iisys : ivsys ) {

      TString hname = "rapHist_" + opt;
      if( iisys == 2 ) 
	hname = "rapHistflipped_" + opt;

      g_dndy = LoadTommyData(iisys, ipart, hname);
      if( g_dndy != NULL ) {
	  
	g_dndy -> SetMarkerStyle(CStyle[iclsys[iisys]].mStyle);
	g_dndy -> SetMarkerSize(1.);
	g_dndy -> SetLineColor( CStyle[iclsys[iisys]].fColor);
	g_dndy -> SetMarkerColor( CStyle[iclsys[iisys]].fColor );

	if( iisys == 2 ) {
	  g_dndy -> SetMarkerStyle(25);
	  g_dndy -> SetMarkerSize(1.);
	  g_dndy -> SetLineColor( CStyle[iclsys[iisys]].fColor );
	  g_dndy -> SetMarkerColor( CStyle[iclsys[iisys]].fColor );
	}	  

	mgr -> Add( g_dndy, "P");
	lg  -> AddEntry( g_dndy, fsys[iisys]+"_Tmy" );
      }
    }
  }


  if( bKaneko ) {
    g_dndy = LoadKanekoData(isys, ipart);
    if( g_dndy != NULL ) {
      g_dndy -> SetMarkerStyle(CStyle[8].mStyle);
      g_dndy -> SetMarkerSize(1.);
      g_dndy -> SetLineColor( CStyle[8].fColor );
      g_dndy -> SetMarkerColor( CStyle[8].fColor );
      mgr -> Add( g_dndy, "P");
      lg  -> AddEntry( g_dndy, CStyle[1].comment );
    }
  }

  // AMD
  if( bAMD ) {

    UInt_t iin = 0;
    for( auto samd : AMDnames ) {

      g_dndy = (TGraphErrors*)LoadAMD(isys, ipart, samd, "g_dndy");

      if( g_dndy != NULL ) {
	g_dndy -> SetFillStyle( samd.fStyle );
	g_dndy -> SetLineColorAlpha( samd.fColor, 0.4 );
	g_dndy -> SetFillColorAlpha( samd.fColor, 0.4 );
	//	g_dndy -> SetFillColor( samd.fColor );

	mgr -> Add( g_dndy, "3");
       	lg  -> AddEntry( g_dndy, "AMD "+ samd.config);
      } 
      // else
      // 	LOG(ERROR) << "AMD g_dndy is not found." << samd.config << FairLogger::endl;
    }
  }


  if(ipart == 0 ) {
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
    lg  -> Draw();
  }

  
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  ccv -> SetGrid(1);
  mgr -> Draw("AP");
  plabel.SetTextAlign(13);
  //  plabel.DrawLatexNDC(0.15,0.9, fsys[isys]+"_"+lpid[ipart]);
  plabel.DrawLatexNDC(0.15,0.9, lpid[ipart]); 
}

void Draw_dndydyx(UInt_t isys, UInt_t ipart, std::vector<gplot> gdata, Int_t ylim=-1.) 
{
  LOG(INFO) << "Draw_dndydyx .......... "<< FairLogger::endl;
  Bool_t bAMD  = 1;

  TString hists[] = {"h_dndy","h_dndyx"};
  TString label[] = {" Long"," Trans"};

  THStack*  mhsty[2];
  mhsty[0]  = new THStack("mdndy","");
  mhsty[0] -> SetName(Form("dndydy_%d_%d",isys,ipart));
  mhsty[0] -> SetTitle(";y_{cm}/y_{prj};dN/dy");

  mhsty[1] = new THStack("mdndyx","");
  mhsty[1] -> SetName(Form("dndydyx_%d_%d",isys,ipart));
  mhsty[1] -> SetTitle(";yx_{cm}/y_{prj};dN/dyx");

  auto lg  = new TLegend(0.15, 0.2, 0.5, 0.7,"");
  lg -> SetHeader(Form("|y|<%3.1f",ycut[ylim]));

  TH1D* h_dndy = NULL;

  std::vector< UInt_t > ivsys;
  ivsys.push_back(isys);
  if( isys == 3 ) 
    ivsys.push_back(2);

  UInt_t i = 0;
  for( auto hist : hists ){

    UInt_t ii = 0;
    for( auto igname : gdata )for( auto iisys : ivsys) {

	//    h_dndy->Reset();
	h_dndy = (TH1D*)LoadData(iisys,ipart,hist,igname,ylim); 
	if( h_dndy != NULL ) {
	  cout<< hist << " --> " << h_dndy -> GetName()<< endl;
	  h_dndy -> SetMarkerStyle(CStyle[ii].mStyle);
	  h_dndy -> SetMarkerSize(1.);
	  h_dndy -> SetLineColor( CStyle[ii].fColor );
	  h_dndy -> SetMarkerColor( CStyle[ii].fColor );
	  h_dndy -> SetMarkerSize(1.);
	  mhsty[i] -> Add( h_dndy);
	  if( i == 0 )
	    lg  -> AddEntry( h_dndy, "DATA"+fsys[iisys]+" "+igname.config1);
	  //	  lg  -> AddEntry( h_dndy, "DATA"+igname.config1);
	}
	ii++;
      }
    i++;
  }
    // AMD
  if( bAMD ) {    

    i = 0;      
    for( auto hist : hists ) {
      for( auto samd : AMDnames ) {
	
	TString hname = hist;
	if( i == 1 )
	  hname = hist + Form("_y%d",(Int_t)(ycut[ylim]*10));

	cout <<" AMD " << samd.config << " " << hist << endl;
	h_dndy = (TH1D*)LoadAMD(isys, ipart, samd, hname);
	if( h_dndy != NULL ) {
	  h_dndy -> SetFillStyle( samd.fStyle );
	  h_dndy -> SetLineColor( samd.fColor );
	  h_dndy -> SetFillColor( samd.fColor );
	  mhsty[i] -> Add( h_dndy,"E4");
	  if( i == 0 )
	  lg  -> AddEntry( h_dndy, "AMD"+ samd.config);
	}   
      }
      i++;
    }
  }

  if( ipart == 0 ) {
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
    lg  -> Draw();
  }

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),600,800); iccv++;
  const Int_t Ny = 2;
  Float_t lMargin = 0.02;
  Float_t rMargin = 0.02;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.02;

  CanvasPartitionY(ccv,Ny,lMargin,rMargin,bMargin,tMargin);
  TPad* pad[Ny];

  for( auto ipad : ROOT::TSeqI(Ny) ) {
    ccv->cd(0);
    TString pname = Form("pad_%d",ipad);
    pad[ipad] = (TPad*)gROOT->FindObject(pname);
    pad[ipad] -> Draw();
    pad[ipad] -> SetFillStyle(4000);
    pad[ipad] -> SetFrameFillStyle(4000);
    pad[ipad] -> SetGrid(1);
    pad[ipad] -> cd();

    mhsty[Ny-ipad-1] -> Draw("nostack");
    if( ipad == Ny-1 ) {
      plabel.SetTextAlign(13);
      plabel.DrawLatexNDC(0.2,0.9, fsys[isys]+"_"+lpid[ipart]);
      //      plabel.DrawLatexNDC(0.2,0.8, lpid[ipart]); 
    }
  }
}


void Draw_dndyRatio(UInt_t isys, UInt_t ipart, TString gname) 
{
  Bool_t bAMD = 0;
  Bool_t bImQMD = 0;
  Bool_t bpBUU = 0;
  Bool_t bTommy = 1;
  Bool_t bKaneko = 1;
  

  auto mgr = new TMultiGraph();
  mgr -> SetTitle(gname+";y_{cm};dN/dy");
  auto lg  = new TLegend(0.7, 0.7, 0.95, 0.95,"");
  lg -> SetHeader(bName[isys]+lpid[ipart]);

  TGraphErrors *g_dndy[2] = {NULL, NULL};

  
  if( bTommy ) {
    g_dndy[0] = LoadTommyData(isys, ipart, gname);
    if( g_dndy[0] != NULL ) {
      g_dndy[0] -> SetMarkerStyle(CStyle[0].mStyle);
      g_dndy[0] -> SetMarkerSize(1.);
      g_dndy[0] -> SetLineColor( CStyle[0].fColor );
      g_dndy[0] -> SetMarkerColor( CStyle[0].fColor );
      mgr -> Add( g_dndy[0], "P");
      lg  -> AddEntry( g_dndy[0], CStyle[0].comment );
    }
  }

  if( bKaneko ) {
    g_dndy[1] = LoadKanekoData(isys, ipart);
    if( g_dndy[1] != NULL ) {
      g_dndy[1] -> SetMarkerStyle(CStyle[1].mStyle);
      g_dndy[1] -> SetMarkerSize(1.);
      g_dndy[1] -> SetLineColor( CStyle[1].fColor );
      g_dndy[1] -> SetMarkerColor( CStyle[1].fColor );
      mgr -> Add( g_dndy[1], "P");
      lg  -> AddEntry( g_dndy[1], CStyle[1].comment );
    }
  }

  // AMD
  if( bAMD ) {
    for( auto samd : AMDnames ) {
      
      auto g_dndy = (TGraphErrors*)LoadAMD(isys, ipart,  samd, "dndy");
      if( g_dndy != NULL ) {
	g_dndy -> SetMarkerStyle(samd.fStyle);
	g_dndy -> SetLineColor( samd.fColor);
	g_dndy -> SetMarkerColor( samd.fColor );

	//	mgr -> Add( g_dndy, "4");
	//	lg  -> AddEntry( g_dndy, "AMD+"+samd.config);
      }   
    }
  }

  auto g_dndyr = new TGraphErrors();
  UInt_t iin = 0;
  auto kx = g_dndy[0]->GetX();
  for( auto ix : ROOT::TSeqI(g_dndy[0]->GetN()) ) {
    //    cout << *(kx+ix) << " -> " << g_dndy[0]->Eval(*(kx+ix)) << endl;

    Double_t n1 = g_dndy[1]->Eval(*(kx+ix));
    if( n1 > 0 ) {
      Double_t n_ratio = n1/g_dndy[0]->Eval(*(kx+ix));
      
      g_dndyr -> SetPoint(iin, *(kx+ix), n_ratio);
      iin++;
    }
  }
    
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;

  g_dndyr -> GetYaxis()->SetRangeUser(0.8, 1.2);
  g_dndyr -> Draw("ALP");

}


void Draw_Ratio(UInt_t igname, TString gname)
{
  TGraphErrors *grp[3];
  TMultiGraph  *mv;
  TLegend      *lg;

  for( auto ipart : {0,1,2,3,4} ) {

    grp[0] = (TGraphErrors*)LoadData( 1, ipart, gname, gnames[igname]);

    mv = new TMultiGraph();
    mv -> SetName(Form("mv_%d",ipart));
    mv -> SetTitle(";U_{t0};R(v2)");
    lg = new TLegend(0.2, 0.2, 0.5, 0.5,lpid[ipart]);

    for( auto isys : {0,3} ) {

      grp[2] = new TGraphErrors();
      grp[2] -> SetName(lsys[isys]+gname);



      grp[1] = (TGraphErrors*)LoadData( isys, ipart, gname, gnames[igname]);

      UInt_t inn = 0;
      for( auto i : ROOT::TSeqI(grp[0]->GetN()) ) {
      
	Double_t x0,y0,y0e;
	grp[0] -> GetPoint( i, x0, y0);
	y0e = grp[0] -> GetErrorY(i);
	
	Double_t x, y, ye;
	grp[1] -> GetPoint( i, x, y);
	ye = grp[1] -> GetErrorY(i);
	
	Double_t ratio, ratioe;
	if( y0 != 0  ) {
	  ratio  = abs(y / y0);
	  ratioe = GetError( y, y0, ye, y0e);

	  grp[2] -> SetPoint( inn, x, ratio );
	  grp[2] -> SetPointError( inn, 0, ratioe);
	  inn++;
	}  
      }

      grp[2] -> SetMarkerStyle(imark[isys]);
      grp[2] -> SetMarkerColor(icol[isys]);
      grp[2] -> SetLineColor(icol[isys]);
      
      grp[0] ->  Print();
      grp[1] -> Print();
      grp[2] -> Print();

      mv -> Add( grp[2], "pl" );
      lg -> AddEntry( grp[2], lsys[isys]+"/"+lsys[1]);
    }
     
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++; 
    mv -> GetYaxis()->SetRangeUser(0.5,1.5);
    mv -> Draw("AP");
    lg -> Draw();
    
    
  }
}

//----------------------------------------------------------------------
//--------------------

void PlotFigure(Bool_t nplot=kFALSE)
{
  if( ndata == 0 ) {
    LOG(ERROR) << " No data is selected. " << FairLogger::endl;
    exit(0);
  }

  gStyle -> SetLegendFillColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatStyle(0);
  gStyle->SetOptFit(1111);
  gStyle->SetStatColor(10);

  SetStyle();
  SetColor();

  if( nplot )
    return 0;
 
//##main##
  bCentral = 1;

  if( bCentral ) {
    AMDnames = AMDnameCentral;
    gnames   = gnamesCentral;
  }
  else{ 
    AMDnames = AMDname202111;
    gnames   = gnamesMidCentral;
  }

  cout << " gnames[0] -> " <<gnames[0].Version << endl;

  if(1) {
    std::vector<gplot> gname = {gnames[0]};
    for( auto isys : {0})for( auto ipart : {0,1,2,3,4} ){
	Draw_dndy(isys, ipart, gname);
      }
  }
  
  //Double_t ycut[5]= {2.0, 0.5, 0.4, 0.3, 0.2};
  UInt_t ycutid = 3;
  if(0) {
    for( auto isys : {1} )
      Draw_vartxz(isys,ycutid); //<-- The last index is ycut[id]
  }

  if(0) {
    std::vector<gplot> gname = {gnames[0]};
    for( auto isys : {1})for( auto ipart : {0,1,2,3,4} ) {
	Draw_dndydyx(isys,ipart,gname,ycutid); //<-- The last index is ycut[id]
      }
  }

  if(0) {
    std::vector<gplot> gname = {gnames[1]};//, gnames[3]};
    for( auto isys : {3})for( auto ipart : {0,1,2,3,4} ) {
	  Draw_v_y_Model(isys,ipart,gname,1,"y/y_{nn}-1"); //112
	}
  }

  // if(0)
  //   Draw_DeltaDIndiv(igname, "gu_v1");


  // Draw function
  UInt_t igname = 0;
  for( auto gname : gnames ) {

    //--- v11 and v20 system dependece

    if( ndata == 1 ){

      if( kFALSE )//@@@@
	Draw_Ratio(igname, "g_utv1_0" );
      if( kFALSE )//@@@@
	Draw_Ratio(igname, "g_utv2" );

      // if( kFALSE )//@@@@
      // 	Draw_Ut_ProtonRatioSystemD(igname, "g_utv1_0");
      // if( kFALSE )
      // 	Draw_Ut_ProtonRatioSystemD(igname, "g_utv2");

      // if( kFALSE )//@@@@ 
      // 	Draw_Ut_ParticleSystemD(igname,"g_utv1_0");
      // if( kFALSE ) 
      // 	Draw_Ut_ParticleSystemD(igname,"g_utv2");

      if( kFALSE )
	Draw_DeltaDIndiv(igname, "gu_v1_v13");
      if( kFALSE )
	Draw_DeltaDIndiv(igname, "gu_v1_v13/v11");

      if( kFALSE )
	Draw_DeltaDIndiv(igname, "gu_v2");

      //v1
      if(0)  //JPS2021
	Draw_v_y_Model(igname,0,1,"y/y_{nn}-1"); //132
      if(0)  //JPS2021
	Draw_v_y_Model(igname,1,1,"y/y_{nn}-1"); //108


      //v2
      if(0)  //JPS2021
	Draw_v_y_Model(igname,0,2,"y/y_{nn}-1"); //132
      if(0)  //JPS2021
	Draw_v_y_Model(igname,1,2,"y/y_{nn}-1"); //108
      if(0)  //JPS2021
	Draw_v_y_Model(igname,3,2,"y/y_{nn}-1"); //112

      
      //      Draw_v_y_Model(igname,0,1,"y"); //132

      if( kFALSE ) {
	Draw_v_y_Model(igname,0,1,"ut0");  //132 v1
	Draw_v_y_Model(igname,1,1,"ut0");  //108 v1
	Draw_v_y_Model(igname,3,1,"ut0");  //132 v1
      }

      if( kFALSE ) {
	Draw_v_y_Model(igname,0,2,"ut0");  //132 v2
	Draw_v_y_Model(igname,1,2,"ut0");  //108 v2
      }





      if( kFALSE ) 
	Draw_Ut_Comparison(igname,0); // Comparison with AMD and pBUU (0:system) 

      if( kFALSE ) {
	for( auto ii : {0} )
	  Draw_dndyRatio(0, ii, "rapHist_M55_55");
      }

      if( kFALSE ) {
	for( auto ii : {0,1,2,3,4} )
	  Draw_dndy(0, ii, "rapHist_M55_55");
      }

      if( kFALSE ) {
	for( auto ii : {0,1,2} )
	  Draw_dndy(1, ii, "rapHist_M55_55");
      }


      if( kFALSE ) {
	for( auto ii : {0,1,2,3,4} )
	  Draw_dndy(3, ii, "rapHist_M43_52");
      }
      //@@@@@@@


      if( kFALSE ) 
	Draw_Indiv_v1SystemD(igname);
      if( kFALSE ) 
	Draw_Indiv_v2SystemD(igname);


      if( kFALSE )
	Draw_DeltaDOne(igname, "gu_v1");
      if( kFALSE )
	Draw_DeltaDOne(igname, "gu_v2");
      //--------------------


      // if( kFALSE )
      // 	Draw_Ut_RatioToOne(igname, "g_utv1_0");
      // if( kFALSE )
      // 	Draw_Ut_RatioToOne(igname, "g_utv2");


      if( kFALSE ) //OK
	Draw_Ut_SystemD(igname); // Final

      if( kFALSE ) 
	Draw_Ut_Comparison_FOPI(igname, 0); // Comparison with FOPI


      // if( kFALSE )
      //   Draw_Ut_Particle();

      if( kFALSE ) 
	Draw_v_y(igname, 1);
      if( kFALSE ) 
	Draw_v_y(igname, 2);

      if( kFALSE ) //@@@@@ fig0
	Draw_ParticleD("gu_v1", 0);
      if( kFALSE ) 
	Draw_ParticleD("gu_v2", 0);

      if(0) 
	Draw_v_y_System(igname, 1, "y");

      if(0) 
	Draw_v_y_System(igname, 2, "y");
            
    }

    if( kFALSE ) 
      Draw_DeltaDIndiv(igname, "gu_v1");

    igname++;
  }

  if( kFALSE ) {
    for( auto ii : {0,1,2,3,4} )
      Draw_dndy(3, ii);
  }
  
  //---  Multiplicity dependence
  if( ndata > 3) {
    if( outFile ) outFile->Close();
    outFile = new TFile("data/PlotFigureMD_"+gnames[0].Version+".root","recreate");
    if( kFALSE )
      Draw_MultiplicityD("gu_v1");
    if( kFALSE )
      Draw_MultiplicityD("gu_v2");
    if( kFALSE )
      Draw_MultiplicityRatio("gu_v1");
    if( kFALSE )
      Draw_MultiplicityRatio("gu_v2");

    if( kFALSE )
      Draw_RPResolutionDM();

  }


  if( kFALSE ) 
    Draw_v20_Edependence();
}
