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

gplot gnames[] = {
  // {".v52.15.100" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.102" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.104" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.106" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.108" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.110" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.112" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.114" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.116" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.118" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  // {".v52.15.120" ,"finYPt_" ,"m55to80","|#phi|<45","y_nn"},
  //-----
  {".v52.15.122" ,"finYPt_" ,"5fm","|#phi|<45","y_nn"},
  //  {".v52.15.123" ,"finYPt_" ,"5fm","|#phi|>135","y_nn"},
  //  {".v52.15.122.v52.15.123" ,"ut_finYPt_" ,"5fm","ave","y_nn"},
  //-----
};

struct mplot{
  TString Version;
  TString config;
};

mplot pBUUname[] = {
  {"Soft_s030.2L38.7"  ,"pBUU s030.2L38.7"      },
  {"Soft_s032.9L64"    ,"pBUU s032.9L64"     },
  {"Soft_s037.4L105.5" ,"pBUU s037.4L105.5"  },
  {"Stiff_s037.4L105.5","pBUU s037.4L105.5(STF)"}
};

Bool_t bAMD = 1;
Bool_t bImQMD = 1;
Bool_t bpBUU = 1;


const UInt_t ndata = (UInt_t)sizeof(gnames)/sizeof(gplot);
TString vfitfname = "PlotFigure.v52.15.51.v52.15.52.root";

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
void Draw_v_y_oneSystem(UInt_t igname, UInt_t vn = 1, UInt_t sSys = 0);
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

TGraphErrors* LoadTextGraph(TString fname);

void PlotFigure()
{
  if( ndata == 0 ) {
    LOG(ERROR) << " No data is selected. " << FairLogger::endl;
    exit(0);
  }

  gStyle->SetOptStat(0);
  gStyle->SetStatStyle(0);

  SetStyle();
  SetColor();


  // Draw function
  for( auto igname : ROOT::TSeqI( ndata ) ) {

    //--- v11 and v20 system dependece


    if( ndata == 1 ){
      if( kFALSE )//@@@@
	Draw_Ut_ProtonRatioSystemD(igname, "g_utv1_0");
      if( kFALSE )
	Draw_Ut_ProtonRatioSystemD(igname, "g_utv2");

      if( kFALSE )//@@@@ 
	Draw_Ut_ParticleSystemD(igname,"g_utv1_0");
      if( kFALSE ) 
	Draw_Ut_ParticleSystemD(igname,"g_utv2");

      if( kFALSE )
	Draw_DeltaDIndiv(igname, "gu_v1");
      if( kFALSE )
	Draw_DeltaDIndiv(igname, "gu_v1_v13");
      if( kFALSE )
	Draw_DeltaDIndiv(igname, "gu_v1_v13/v11");

      if( kFALSE )
	Draw_DeltaDIndiv(igname, "gu_v2");

      if( kFALSE ) {
	// Draw_v_y_oneSystem(igname, 1, 0);
	// Draw_v_y_oneSystem(igname, 1, 1);
	Draw_v_y_oneSystem(igname, 1, 3);
      }

      //v1
      if( kFALSE ) {
	Draw_v_y_Model(igname,0,1,"y"); //132
	Draw_v_y_Model(igname,1,1,"y"); //108
	Draw_v_y_Model(igname,3,1,"y"); //112
      }

      //v2
      if( kFALSE ) {
	Draw_v_y_Model(igname,0,2,"y"); //132
	Draw_v_y_Model(igname,1,2,"y"); //108
	Draw_v_y_Model(igname,3,2,"y"); //112
      }
      
      //      Draw_v_y_Model(igname,0,1,"y"); //132

      if( 1 ) {
	Draw_v_y_Model(igname,0,1,"ut0");
	Draw_v_y_Model(igname,1,1,"ut0");
      }

      if( kFALSE )
	Draw_Ut_Ratio(igname, "g_utv1_0");
      if( kFALSE )
	Draw_Ut_Ratio(igname, "g_utv2");


      if( kFALSE ) 
	Draw_Ut_Comparison(igname,0); // Comparison with AMD and pBUU (0:system) 
      //@@@@@@@@@@@@@@@@


      if( kFALSE ) 
	Draw_Indiv_v1SystemD(igname);
      if( kFALSE ) 
	Draw_Indiv_v2SystemD(igname);


      if( kFALSE )
	Draw_DeltaDOne(igname, "gu_v1");
      if( kFALSE )
	Draw_DeltaDOne(igname, "gu_v2");
      //--------------------


      if( kFALSE )
	Draw_Ut_RatioToOne(igname, "g_utv1_0");
      if( kFALSE )
	Draw_Ut_RatioToOne(igname, "g_utv2");


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

      if( kFALSE ) 
	Draw_v_y_oneSystem(igname, 2, 0);
      
      
    }

    if( kFALSE ) 
      Draw_DeltaDIndiv(igname, "gu_v1");

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

    if( 1 )
      Draw_RPResolutionDM();

  }


  if( kFALSE ) 
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


TGraphErrors* LoadAMDText(UInt_t isys, UInt_t ipart, UInt_t vn=1, TString grname="rapflow", TString eos="")
{
  if( isys != 3 ) return NULL;


  if( vn > 2 || vn == 0 ) vn = 1;

  TString pname[] = {"p", "d", "t", "h", "a", "n"};
  TGraphErrors *gvt = NULL;

  TString dirname = "ModelData/AMD/" + eos;
  TString ifile = grname + "-" + pname[ipart] + ".dat"; 

  //  LOG(INFO) << dirname << "/" << ifile  << FairLogger::endl;

  if( !gSystem -> FindFile(dirname, ifile ) ) {
    LOG(ERROR) << ifile << " is not found. " << FairLogger::endl;
    return gvt;
  }
  else {
    LOG(INFO) << ifile << " is accessed. " << FairLogger::endl;
    
    std::fstream fread;
    fread.open(ifile, std::fstream::in);
    gvt = new TGraphErrors();

    Double_t x, y, ye;
    TString sget;
    TString seget;
    UInt_t  in = 0;
    
    while( !fread.eof() ) {
      
      for( auto ict : ROOT::TSeqI(3*5) )
	fread >> sget;


      while( !fread.eof() ) {
	fread >> sget;
	x = (Double_t)atof(sget);
      
	for( auto i : ROOT::TSeqI(5) ) {
	  fread >> sget;
	  fread >> seget;
	  if( (vn==1 && i == 3) || (vn==2 && i==4) ){
	    y = (Double_t)atof(sget);
	    ye = (Double_t)atof(seget);
	  }
	}
	///@@@@
      
	//	cout << " x " << x << " y " << y << " +- " << ye << endl;
	if( ye != 0. ) {
	  x = x/y_cm[10];
	  // if( isys == 3 ) {
	  //   x *= -1.;
	  //   if( vn == 1 )
	  //     y *= -1;
	  // }
	  
	  gvt -> SetPoint(in, x, y);
	  gvt -> SetPointError(in, 0, ye);

	  in++;
	}
      }
    }
  }

  return gvt;

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


TGraphErrors* LoadpBUU(UInt_t isys, UInt_t ipart, TString grname="v1_Ut0", TString eos="Soft") 
{
  TGraphErrors  *gvt = NULL;

  if( ipart == 4 ) return gvt;
  
  TString dirname = "ModelData/pBUU";
  TString ifile = grname + "_pd3Het_" + rsys[isys] + eos + ".txt";
  //  LOG(INFO) << dirname << ifile << " will be opened. " << FairLogger::endl;

  if( !gSystem -> FindFile(dirname, ifile ) ) 
    LOG(ERROR) << ifile << " is not found. " << FairLogger::endl;
  else {
    LOG(INFO) << ifile << " is accessed. " << FairLogger::endl;
    

    std::fstream fread;
    fread.open(ifile, std::fstream::in);
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



TGraphErrors* FindData(UInt_t igname, UInt_t isys, UInt_t ipart, TString gname)
{
  TString fname = gnames[igname].fileHeader + bName[isys] + fpid[ipart] + gnames[igname].Version + ".root";
  return (TGraphErrors*)gROOT->FindObject(fname);
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
  
  //  grv -> SetDirectory(gROOT);

  fOpen->Close();

  return grv;
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

Double_t* GetV11(UInt_t igname, UInt_t isys, UInt_t ipart, TString sData="DATA", TString eos = "Soft")
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
    yv1 = LoadpBUU(isys, ipart, "v1_yn", eos);
    ft_low  = -0.7;
    ft_high =  0.7;
  }
  else if( sData == "ImQMD" ) {
    yv1 = LoadImQMD(isys, ipart, "v1", "");
    ft_low  = -0.8;
    ft_high =  0.8;
  }
  else
    yv1 = LoadData(igname, isys, ipart, "gu_v1");
  //    yv1 = FindData(igname, isys, ipart, "gu_v1");

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

TGraphErrors* GetSDGraph(TString gname, UInt_t igname, UInt_t ipart, TString sData="DATA", TString eos="Soft",
			 Double_t fit_low=v1fit[0], Double_t fit_high=v1fit[1])
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
    fitval = GetV11(igname, isys, ipart, sData, eos);
    if( *fitval == 0 ) continue;

    UInt_t ipara = 0;
    if( gname.Contains("_v13") )
	ipara += 2;
	
    grp -> SetPoint( inn, *(sysdlt+isys), *(fitval+ipara) );
    grp -> SetPointError( inn, 0, *(fitval+ipara+1) );
    inn++;
  }

  return grp;
}

void Draw_DeltaDIndiv(UInt_t igname, TString gname) 
{


  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 512, 840); iccv++;
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
	grp = GetSDGraph(gname,igname,j, "pBUU",pBUUname[i].Version);
	if( grp == NULL ) continue;
      
	grp -> SetLineWidth(20);
	grp -> SetLineColor(mdlcol[3+i]);
	grp -> SetMarkerStyle(0);
	mv -> Add(grp, "AL");
	lg -> AddEntry(grp, pBUUname[i].config);
      }
    }

    // ImQMD
    if( bImQMD ) {
      grp = GetSDGraph(gname,igname,j, "ImQMD");
      if( grp ) {
	grp -> SetLineWidth(20);
	grp -> SetLineColor(mdlcol[6]);
	grp -> SetMarkerStyle(0);
	mv -> Add(grp, "AL");
	lg -> AddEntry(grp, "ImQMD");
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

    mv -> Add(grp, "AP");
    lg -> AddEntry(grp, "Data");

    

    Double_t alpha = 0.4;

    if( gname.Contains("v2") )    alpha = 0.2;

    Double_t ylow  = Double_t( Int_t((1.-alpha)*ymean*1000.) / 1000. );
    Double_t yhigh = Double_t( Int_t((1.+alpha)*ymean*1000.) / 1000. );
    mv -> SetMaximum( yhigh );
    mv -> SetMinimum( ylow );

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

void Draw_v_y_oneSystem(UInt_t igname, UInt_t vn = 1, UInt_t sSys = 0 )
{
  if( vn > 2 ) vn = 1;
  LOG(INFO) << "Draw_v" << vn << "_y_one"  << FairLogger::endl;


  auto mv = new TMultiGraph(Form("mv%d",vn), Form(";y/y_{nn}-1;v%d",vn));

  TGraphErrors* vy;
  TLegend *lg;
  lg -> SetFillColor(0);

  for( auto isys : {sSys} ) {
    lg = new TLegend(0.6, 0.15, 0.9, 0.5,fsys[isys]); 
    if( vn == 2 ){
      lg -> SetY1(0.6);
      lg -> SetY2(0.9);
    }
    
    for( auto ipart : sqPart ){
    
      vy = LoadData(igname, isys, ipart, Form("gu_v%d",vn));
      LOG(INFO) << vy -> GetName() << " is registred. " << FairLogger::endl;


      if( vy != NULL ) {
	vy -> SetMarkerStyle(20);
	vy -> SetMarkerSize(1.);
	vy -> SetLineColor(icol[ipart]);
	vy -> SetMarkerColor(icol[ipart]);

	auto f1 = vy -> GetFunction(Form("fv%dfit",vn));
	if( f1 != NULL )
	  f1 -> SetLineColor(icol[ipart]);
	else {
	  vy -> Fit(Form("fv%dfit",vn),"","",v1fit[0], v1fit[1]);
	  f1 = vy -> GetFunction(Form("fv%dfit",vn));
	  f1 -> SetLineColor(icol[ipart]); 
	}

	mv -> Add( vy, "p" );
	lg -> AddEntry(vy, fpid[ipart]);
      }
      else
	LOG(ERROR) << " gu_v1  is not found. " << isys << " : " << ipart << FairLogger::endl;
    }
  }

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 600, 500); iccv++;

  mv -> GetXaxis()->SetRangeUser(-0.65,0.8);
  mv -> GetXaxis()->SetNdivisions(505);
  mv -> GetYaxis()->SetNdivisions(505);
  mv -> Draw("AP");
  lg -> Draw();

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

void Draw_v_y(UInt_t igname, UInt_t vn = 1) 
{
  if( vn > 2 ) vn = 1;

  Bool_t bAMD = 0;
  Bool_t bImQMD = 1;
  Bool_t bpBUU = 1;


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
      v_y = LoadData(igname, isys, sqPart[j], Form("gu_v%d",vn));
      LOG(INFO) << v_y -> GetName() << " is registred. " << FairLogger::endl;

      if( v_y != NULL ) {
	v_y -> SetMarkerStyle(mstyle[0]);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor( mdlcol[0] );
	v_y -> SetMarkerColor( mdlcol[0] );

	mv -> Add( v_y, "p" );
	lg -> AddEntry(v_y, "Data : "+lsys[isys]);
      }
      else
	LOG(ERROR) << " gu_v1  is not found. " << isys << " : " << sqPart[j] << FairLogger::endl;
      
      if( bAMD ) {
	v_y = LoadAMD(sqSys[i], sqPart[j], Form("v%d_y",vn), "SLy4");
	
	if( v_y != NULL ) {
	  v_y -> SetMarkerStyle(mstyle[1]);
	  v_y -> SetMarkerSize(1.);
	  v_y -> SetLineColor(mdlcol[1]);
	  v_y -> SetMarkerColor(mdlcol[1]);
	
	  mv -> Add( v_y, "pl" );
	  lg -> AddEntry(v_y, "AMD-SLy4");
	}

	v_y = LoadAMD(sqSys[i], sqPart[j], Form("v%d_y",vn), "SLy4_L108");
      
	if( v_y != NULL ) {
	  v_y -> SetMarkerStyle(mstyle[1]);
	  v_y -> SetMarkerSize(1.);
	  v_y -> SetLineColor(mdlcol[2]);
	  v_y -> SetMarkerColor(mdlcol[2]);

	  mv -> Add( v_y, "pl" );
	  lg -> AddEntry(v_y, "AMD-L108");
	}
      }

      if( bpBUU ) {
	
	for( auto i : ROOT::TSeqI(sizeof(pBUUname)/sizeof(mplot)) ) {
	  if( vn == 1 )
	    v_y = LoadpBUU(sqSys[i], sqPart[j], "v1_yn", pBUUname[i].Version);
	  if( vn == 2 )
	    v_y = LoadpBUU(sqSys[i], sqPart[j], "v2_yn", pBUUname[i].Version);
	  if( v_y == NULL ) continue;
      
	  v_y -> SetMarkerStyle(mstyle[2]);
	  v_y -> SetMarkerSize(1.);
	  v_y -> SetLineColor(mdlcol[3+i]);
	  v_y -> SetMarkerColor(mdlcol[3+i]);
	  
	  mv -> Add( v_y, "pl" );
	  lg -> AddEntry(v_y, pBUUname[i].config);
	}
      }

      if( bImQMD ) {
	if( vn == 1)
	  v_y = LoadImQMD(sqSys[i], sqPart[j], "v1" );
	else if( vn == 2)
	  v_y = LoadImQMD(sqSys[i], sqPart[j], "v2" );
      if( v_y != NULL ){
	v_y -> SetMarkerStyle(mstyle[3]);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor(mdlcol[6]);
	v_y -> SetMarkerColor(mdlcol[6]);

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

void Draw_v_y_Model(UInt_t igname, UInt_t isys, UInt_t vn=1, TString para="y") 
{
  if( vn > 2) vn = 1;
  TGraphErrors* v_y;
  TString ypara;

  for(auto ipart: {0,1,2,3,4} ) {

    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 600); iccv++;
    gStyle->SetOptTitle(0);

    TMultiGraph *mv = new TMultiGraph();
    mv -> SetName(Form("mv_%d",ipart));
    mv -> SetTitle(";"+para+Form("; v%d",vn));

    TLegend *lg = new TLegend(0.2,0.7,0.6,0.9,lpid[ipart]);
    lg -> SetFillColor(0);

    if( para == "y")
      ypara = Form("gu_v%d",vn);
    else 
      ypara = Form("g_utv%d_0",vn);
    
    v_y = LoadData(igname, isys, ipart, ypara); 

    if( v_y != NULL ) {
      v_y -> SetMarkerStyle(mstyle[0]);
      v_y -> SetMarkerSize(1.);
      v_y -> SetLineColor( mdlcol[0] );
      v_y -> SetMarkerColor( mdlcol[0] );
      
      mv -> Add( v_y, "p" );
      lg -> AddEntry(v_y, "Data : "+lsys[isys]);
    }
    else
      LOG(ERROR) << "Failed to get " << ypara << FairLogger::endl;

    if( bAMD && para == "y" ) {
      
      v_y = LoadAMDText(isys, ipart, vn,  "rapflow", "E270-124-112-sly4-l055-mdcorr50-cxAg");
      
      if( v_y != NULL ) {
	v_y -> SetMarkerStyle(mstyle[1]);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor( mdlcol[1] );
	v_y -> SetMarkerColor( mdlcol[1] );

	mv -> Add( v_y, "p" );
	lg -> AddEntry(v_y, "AMD mdcr50");
      }

      v_y = LoadAMDText(isys, ipart, vn, "rapflow", "E270-124-112-sly4-l055-mdcorr20-cxAg");
      //v_y = LoadAMD(isys, ipart, Form("v%d_y",vn), "SLy4_L108");
      
      if( v_y != NULL ) {
	v_y -> SetMarkerStyle(mstyle[1]);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor( mdlcol[2] );
	v_y -> SetMarkerColor( mdlcol[2] );
	  
	mv -> Add( v_y, "p" );
	lg -> AddEntry(v_y, "AMD mdcr20");
      }
    }

    if( bpBUU ) {
	
      for( auto i : ROOT::TSeqI(sizeof(pBUUname)/sizeof(mplot)) ) {
	if( para == "y" )
	  ypara = Form("v%d_yn",vn);
	else 
	  ypara = Form("v%d_Ut0",vn);

	v_y = LoadpBUU(isys, ipart, ypara, pBUUname[i].Version);

	
	if( v_y == NULL ) continue;
      
	v_y -> SetMarkerStyle( mstyle[2] );
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor( mdlcol[3+i] );
	v_y -> SetMarkerColor( mdlcol[3+i] );
	
	mv -> Add( v_y, "p" );
	lg -> AddEntry(v_y, pBUUname[i].config);
      }
    }

    if( bImQMD ) {
      
      if( para == "y" )
	ypara = Form("v%d",vn);
      else
	ypara = Form("v%d_Pt",vn);

      v_y = LoadImQMD(isys, ipart, ypara );
      
      if( v_y != NULL ){
	v_y -> SetMarkerStyle(mstyle[3]);
	v_y -> SetMarkerSize(1.);
	v_y -> SetLineColor( mdlcol[5] );
	v_y -> SetMarkerColor( mdlcol[5] );
	
	mv -> Add( v_y, "p" );
	lg -> AddEntry(v_y, "ImQMD");
      }
    }

    if( para == "y" ) {
      if( vn == 1 ) {
	mv -> GetYaxis() -> SetRangeUser(-1.2, 1.2);
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
	mv -> GetYaxis() -> SetRangeUser(-0.15, 0.15);
	mv -> GetXaxis() -> SetRangeUser(-1.2 , 1.2);
      }
    }

    mv->Draw("AP");
    lg->Draw();
  }
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
      
      cout << " m " << *(res+2) << " +- " << *(res+3)
	   << " cos " << *res << " +- " << *(res+1) << endl;

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
