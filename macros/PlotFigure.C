#include "DoRPRes.C"
#include "SetStyle.C"
//#include "FlowFunction.C"
#include "SetColor.C"
#include "CanvasPartition.C"
#include "ImageProcessing.h"

//@@@@@Parameters==================================================

struct gplot{
  TString Version;
  TString fileHeader;
  TString config1;
  TString config2;
  TString config3;
  UInt_t  mtconfig;  //MNTRK
  UInt_t  phiconfig; //PHICUT
  UInt_t  centrality;
};

//std::vector<gplot> gnames;

std::vector<gplot> gnamesMidCentral = {
  //  {".v1.0.1"      ,"phys_" ,"b0 2~4fm" ,"-20<#phi<30" ," 46<M<55", 1, 5, 0}
  {".v1.0.3"      ,"phys_" ,"b0 2~4fm" ,"-20<#phi<30" ," 46<M<55", 1, 5, 0}
};

std::vector<gplot> gnamesCentral = {
  //  {".v1.0.0"      ,"phys_" ,"b0<0.15" ,"-20<#phi<30" ," M>55", 0, 5, 1},
  {".v1.0.2"      ,"phys_" ,"b0 0-2fm" ,"-20<#phi<30" ," M>55", 0, 5, 1}
};
std::vector<gplot> gnames = {
  //  {".v1.0.3"      ,"phys_" ,"b0 2~4fm" ,"-20<#phi<30" ," 46<M<55", 1, 5, 0},
  //  {".v1.0.2"      ,"phys_" ,"b0 0-2fm" ,"-20<#phi<30" ," M>55",    0, 5, 1},
  // {".v1.0.7"      ,"phys_" ,"b0 2~4fm" ,"-20<#phi<30" ," 46<M<55", 1, 5, 0},
  // {".v1.0.6"      ,"phys_" ,"b0 0-2fm" ,"-20<#phi<30" ," M>55",    0, 5, 1}
  //  {".v1.0.10"        ,"phys_" ,"b0 5fm" ,"-20<#phi<30" ," 28<M<55", 1, 5, 0},
  //  {".v1.0.8"      ,"phys_" ,"b0 0-2fm" ,"-20<#phi<30" ," M>55",    0, 5, 1},
  //  {".v1.0.9"      ,"phys_" ,"b0 2~4fm" ,"-20<#phi<30" ," 46<M<55", 1, 5, 0},
  //
  //=--
  //  {".v1.0.24"        ,"phys_" ,"b0 2-4fm" ,"-Tommy" ," Tommy",    1, 6, 0},
  //  {".v1.0.23"        ,"phys_" ,"b0 2-4fm" ,"-Tommy" ," Tommy",    1, 6, 0},
  //  {".v1.0.22"        ,"phys_" ,"b0 2-4fm" ,"-Tommy" ," Tommy",    1, 6, 0},

  // {".v1.0.16"      ,"phys_" ,"b0 2~4fm" ,"|phi|>140"   ," 46<M<55", 1, 2, 0},
  // {".v1.0.15"      ,"phYs_" ,"B0 0-2fm" ,"|phi|>140"   ," M>55",    0, 2, 1},
  // {".v1.0.14"      ,"phys_" ,"b0 2~4fm" ,"-20<#phi<30&&|phi|>140" ," 46<M<55", 1, 6, 0},
  // {".v1.0.13"      ,"phys_" ,"b0 0-2fm" ,"-20<#phi<30&&|phi|>140" ,"M>55"    , 0, 6, 1}
  //
  //  {".v1.0.17"     ,"phys_" ,"b0 2~4fm" ,"-20<#phi<30" ," 46<M<55", 1, 5, 0}, //_left
  //  {".v1.0.18"     ,"phys_" ,"b0 0~2fm" ,"-20<#phi<30" ,"    M>55", 0, 5, 1}, //_left
  //  {".v2.0.0"      ,"phys_" ,"b0 2~4fm" ,"-20<#phi<30" ," 46<M<55", 1, 5, 0}, //_left
  {".v2.0.1"      ,"phys_" ,"b0 2-4fm" ,"-20<#phi<30&&|phi|>140" ,"46<M<55" , 1, 6, 0},
    {".v2.0.2"      ,"phys_" ,"b0 2-4fm" ,"-20<#phi<30&&|phi|>140" ," Tommy",    1, 6, 0},
};
//@gdata

std::vector< std::pair< Double_t, TString> > fphysdatafile = {
  {  0.,"_fit0to05"},
  {-0.1,"_fit01to05"},
  {-0.2,"_fit02to05"},
  {-0.3,"_fit03to05"},
  {-0.4,"_fit04to05"},
  {-0.5,"_fit05to05"}
};
UInt_t fphysdataid = 3;


Double_t phiacp[] = {180./90., 180./45., 180./45., 180./90., 1., 360./50. };

struct mplot{
  UInt_t  id;
  TString path;
  TString fname;
  TString config;
  TString fullconfig;
  Color_t fColor;
  Color_t fColor2;
  UInt_t  fStyle;
  UInt_t  category;
};
 
std::vector<mplot> pBUUname = {
  {0,"ModelData/pBUU/","Soft_s030.2L38.7"  ,"pBUU s030.2L38.7"      ,"pBUU s030.2L38.7"      ,kMagenta-9, kMagenta-9 ,3001, 1},
  {1,"ModelData/pBUU/","Soft_s032.9L64"    ,"pBUU s032.9L64"        ,"pBUU s032.9L64"        ,kMagenta-3, kMagenta-3 ,3001, 1},
  {2,"ModelData/pBUU/","Soft_s037.4L105.5" ,"pBUU s037.4L105.5"     ,"pBUU s037.4L105.5"     ,kMagenta-1, kMagenta-1 ,3001, 1},
  {3,"ModelData/pBUU/","Stiff_s037.4L105.5","pBUU s037.4L105.5(STF)","pBUU s037.4L105.5(STF)",kMagenta,   kMagenta   ,3001, 1}
};


// category: 0 no plot
//         :>0 Draw_compCorrelation, Draw_compIntegral, Draw_compParameter
//         : 1 Draw_meanpx, Draw_dndbtgm, Draw_dndptRatio, Draw_ParticleDependence, 
//         : 2 Draw_compCorrelation(1)
//         : 3 Draw_compCorrelation(2) 
//         : 4 Draw_meanpx, Draw_dndydX
//----  
//AMDDATA 
std::vector<mplot> AMD2022={
  // { 0,"2022.cxBs/","amd_Sn108Sn112_SLy4_gfg_sgml2",                    "bxXp_gfg"          ,"SLy4_gfg_sgm12"            ,kOrange   ,kOrange   ,3101, 14},
  // { 1,"2022.cxBs/","amd_Sn132Sn124_SLy4_gfg_sgml2",                    "bxXp_gfg"          ,"SLy4_gfg_sgm12"            ,kOrange   ,kOrange   ,3101, 14},
  { 0,"2022.cxBs/","amd_Sn132Sn124-cxBs-rhoxdep10-gd",                 "md50_rhx10_gd"     ,"L46_md50_rhx10_gd"    ,kGray+2   ,kGray+2   ,3001, 142},
  { 1,"2022.cxBs/","amd_Sn124Sn112-cxBs-rhoxdep10-gd",                 "md50_rhx10_gd"     ,"L46_md50_rhx10_gd"    ,kGray+2   ,kGray+2   ,3001, 142},
  { 2,"2022.cxBs/","amd_Sn108Sn112-cxBs-rhoxdep10-gd",                 "md50_rhx10_gd"     ,"L46_md50_rhx10_gd"    ,kGray+2   ,kGray+2   ,3001, 142},
  //
  { 3,"2022.cxBs/","amd_Sn132Sn124-sly4-mdcorr50-cxBs-skms",           "md50_rhx10_skms"   ,"L108_sksm"  ,kYellow   ,kYellow   ,3001, 14},
  { 4,"2022.cxBs/","amd_Sn124Sn112-sly4-mdcorr50-cxBs-skms",           "md50_rhx10_skms"   ,"L108_sksm"  ,kYellow   ,kYellow   ,3001, 14},
  { 5,"2022.cxBs/","amd_Sn108Sn112-sly4-mdcorr50-cxBs-skms",           "md50_rhx10_skms"   ,"L108_sksm"  ,kYellow   ,kYellow   ,3001, 14},
  //
  { 6,"2022.cxBs/","amd_Sn132Sn124-sly4-mdcorr50-b04-cxBs",            "md50"         ,"md50"             ,kOrange   ,kOrange   ,3001, 41},
  { 7,"2022.cxBs/","amd_Sn108Sn112-sly4-mdcorr50-b04-cxBs",            "md50"         ,"md50"             ,kOrange   ,kOrange   ,3001, 41},
  //														           
  { 8,"2022.cxBs/","amd_Sn132Sn124-cxBs-b04-sly4-mdcorr20",            "md20"         ,"md20"             ,kGreen    ,kGreen    ,3014, 42},
  { 9,"2022.cxBs/","amd_Sn124Sn112-cxBs-b04-sly4-mdcorr20",            "md20"         ,"md20"             ,kGreen    ,kGreen    ,3014, 42},
  {10,"2022.cxBs/","amd_Sn108Sn112-cxBs-b04-sly4-mdcorr20",            "md20"         ,"md20"             ,kGreen    ,kGreen    ,3014, 42},
  //
  {11,"2022.cxBs/","amd_Sn132Sn124-cxBs-rhoxdep10",                    "md50"    ,"md50_L46"       ,kPink-2   ,kPink-2   ,3002, 42}, 
  {12,"2022.cxBs/","amd_Sn124Sn112-cxBs-rhoxdep10",                    "md50"    ,"md50_L46"       ,kPink-2   ,kPink-2   ,3014, 42},
  {13,"2022.cxBs/","amd_Sn108Sn112-cxBs-rhoxdep10",                    "md50"    ,"md50_L46"       ,kPink-2   ,kPink-2   ,3014, 42},
  //										   
  {14,"2022.cxBs/","amd_Sn132Sn124-cxBs-rhoxdep10-l108",               "md50_rhx10_L108"   ,"md50_L108"      ,kBlue+2   ,kBlue     ,3003, 40},
  {15,"2022.cxBs/","amd_Sn124Sn112-cxBs-rhoxdep10-l108",               "md50_rhx10_L108"   ,"md50_L108"      ,kBlue+2   ,kBlue     ,3003, 40},
  {16,"2022.cxBs/","amd_Sn108Sn112-cxBs-rhoxdep10-l108",               "md50_rhx10_L108"   ,"md50_L108"      ,kBlue+2   ,kBlue     ,3003, 40},
  //														        
};

UInt_t   bCentral = 1;
TString lbCentral[] = {"2to4fm","0to2fm","mix","152to200M","0to152M"};
// AMD data selection
std::vector<mplot>  AMDnames = AMD2022;
//std::vector<mplot>  AMDnames = AMDKaneko;

const UInt_t nDATA = gnames.size();
const UInt_t nAMD  = AMDnames.size();
TFile* fileEffCor[2];
std::vector< std::vector< std::vector< TFile*>>> fileDATA;
std::vector< std::vector<TFile*>> fileAMD;
  
Bool_t bAMDRoot = 1;
Bool_t bReverse = 0;

TString vfitfname = "PlotFigure.v52.15.51.v52.15.52.root";

Double_t ycut[5]= {2.0, 0.5, 0.4, 0.3, 0.2};

std::vector<UInt_t> sqPart  = {0,1,2,3,4,5,6,7};
std::vector<UInt_t> sqSys   = {0,1,2,3};
std::vector<UInt_t> sqv1sel = {10, 8, 2, 1};
std::vector<UInt_t> sqv2sel = { 4, 3, 2, 1};

UInt_t ix = 0;
Double_t sysBN[]   = {82./50., 58./50.,   74./50.,62./50.,   0,     50./50.};
Double_t sysNA[]   = {156./100.,110./100.,136./100.-0.011,136./100., 0,  156./100.};
Double_t sysN2[]   = {2.43,    1.21,      1.85,   1.85   ,  2.43};
Double_t sysdlt[]  = {(32+24.)/(132+124), (8+12.)/(108+112.), (24+112.)/(124+112.),  (24+112.)/(124+112.), 0.21875}; 

TString FOPI_data_sys[] = {"Au","Ru","Ca"};

Double_t FOPI_AuAu_v11x[4]={0,1,2,4};
Double_t FOPI_AuAu_v11y[4]={0.384, 0.641, 0.800, 1.032};
Double_t FOPI_AuAu_v20[4]={0.048, 0.105, 0.170, 0.247720};

TFile *outFile;
TCanvas *ccv; UInt_t iccv = 0;
TLatex  plabel;
TString xlabel;
Bool_t bsaveData    = 0;

#include "PhysData.C"

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
void Draw_Ut_Particle(Int_t igname, UInt_t isys);
void Draw_v_y(UInt_t igname, UInt_t vn = 1);
void Draw_v_y_System(UInt_t igname, UInt_t vn=1, TString para="y"); 
void Draw_MultiplicityD(TString gname);
void Draw_MultiplicityRatio(TString gname);
void Draw_v20_Edependence();
void Draw_Ut_Ratio(UInt_t igname, TString gname);
void Draw_Ut_RatioToOne(UInt_t igname, TString gname);
void Draw_Ut_ProtonRatio(UInt_t igname, TString gname);
void Draw_Ut_ProtonRatioSystemD(UInt_t igname, TString gname);
void Draw_RPResolutionDM();
Double_t* GetRPResolution(UInt_t igname, UInt_t isys);
void Draw_v_y_Model(UInt_t igname, UInt_t isys, UInt_t vn, TString para);
void Draw_dndy(UInt_t isys, UInt_t ipart, TString opt="");
void Draw_Ratio(UInt_t igname, TString gname);

TGraphErrors* LoadTextGraph(TString fname);


//------------------------------------------------
//------------------------------------------------
//------------------------------------------------
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
      
      TString objname = (TString)obj->GetName();
      if( objname.Contains("v1") || objname.Contains("px") )
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

//@loaddata
TObject* LoadData(UInt_t isys, UInt_t ipart, TString grname, gplot gname, Int_t ylim = -1.)
{
  LOG(INFO) << "[LoadData] Opening.. " << grname << " in " << gname.Version << FairLogger::endl;

  TFile* fOpen;

  // Acceptance correction file
  if( fileEffCor[gname.mtconfig] == NULL )
    fileEffCor[gname.mtconfig] = GetAcceptanceCorrectionFile(gname.mtconfig, gname.phiconfig, mcphicut);

  if( (bsaveData && fileDATA[isys][ipart][gname.centrality] == NULL) || !bsaveData ) {

    TString fname = gname.fileHeader + bName[isys] + ncls[ipart].name + gname.Version + ".root";
    LOG(INFO) << "[LoadData] Opening file is " << fname << FairLogger::endl;

    if( !gSystem->FindFile("data", fname) ) {
      LOG(ERROR) << fname << " is not found " << FairLogger::endl;
      return NULL;
    }
    else {
      fileDATA[isys][ipart][gname.centrality] = TFile::Open( fname, "READ" );
      fOpen = fileDATA[isys][ipart][gname.centrality];
      LOG(INFO) << fname << " is opened. " << FairLogger::endl;
    }
  }
  else 
    fOpen = fileDATA[isys][ipart][gname.centrality]; 
  
  LOG(INFO) << "[LoadData] " << fOpen->GetName() << " is connected." << FairLogger::endl;

  TGraphErrors* grv = NULL;
  TH1D*         hst = new TH1D();
  hst -> Sumw2();
  TString h2dname;
  TString acpcorrName;  

  //---->>> Projection on X
  if( grname == "@@h_dndy1" ) { 
    h2dname = "hypt";
    acpcorrName ="h2PtYEff_108Sn_"+ncls[ipart].Name+"_mbin0_iter1_GaussBlur";    
    if( ipart == 5 || ipart == 6)   
      acpcorrName ="h2PtYEff_108Sn_"+ncls[0].Name+"_mbin0_iter1_GaussBlur"; 
  }
  else if( grname == "h_dndy" ){ 
    h2dname = "hyyx";
    acpcorrName ="h2yyxEff_108Sn_"+ncls[ipart].Name+"_mbin0_iter1_GaussBlur";    
    if( ipart == 5 || ipart == 6)   
      acpcorrName ="h2yyxEff_108Sn_"+ncls[0].Name+"_mbin0_iter1_GaussBlur"; 
  }


  if( h2dname != "" ) {
    TH2D* hacp = (TH2D*)fOpen->Get(h2dname);

    if( hacp == NULL ) {
      LOG(ERROR) << h2dname<< " is not found " << FairLogger::endl;
      return NULL;
    }    
    LOG(INFO) << hacp->GetName() << " is found " << FairLogger::endl;

    TH2D* hcorrected;
    if(  fileEffCor[gname.mtconfig] ) {
      TH2D* hacpcorr = (TH2D*)fileEffCor[gname.mtconfig]->Get(acpcorrName);
      
      if( hacpcorr == NULL ) {
	LOG(ERROR) << " Efficiency correction : " << acpcorrName << " is not found " << FairLogger::endl;
	hcorrected = hacp;
      }
      else {
	LOG(INFO) << hacpcorr->GetName() << " is loaded." << FairLogger::endl;
	
	TH2D* h_smear = (TH2D*)ImageProcessing::GaussianBlur(hacp, 3,1.);
	hcorrected = GetAcpCorrection(*h_smear, hacpcorr);
	
	if( hcorrected ) 
	  LOG(INFO) << h2dname << ": Efficiency correciton have been done." << FairLogger::endl;      
	
	else {
	  LOG(ERROR) << h2dname << ": Efficiency correciton is failed(A)." << FairLogger::endl;
	  hcorrected = hacp;
	}
      }
    }
    else {
      LOG(ERROR) << h2dname << ": Efficiency correciton file is not found.." << FairLogger::endl;
      hcorrected = hacp;
    }

    hst -> Reset();
    hst -> SetName(grname+"_"+rsys[isys]+"_"+ncls[ipart].sName+"_DATA");
    hst = hcorrected->ProjectionX(hst->GetName());
    hst -> Scale(hcorrected->GetYaxis()->GetBinWidth(0));
    
    // cut out y<-0.5 because of poor acceptance
    for( auto iy : ROOT::TSeqI(hst->GetNbinsX()) ) {
      if( hst -> GetBinCenter(iy) < -0.5 ) {
	hst -> SetBinContent(iy,0.);
	hst -> SetBinError(iy,0.);
      }
    }

    hst -> SetTitle(";;");
    
    return hst;
  }


  //-------->>>> Projection on Y
  if( grname == "h_dndpx" ) { 
    h2dname = "hypxrnd";
    acpcorrName ="h2ypxrndEff_108Sn_"+ncls[ipart].Name+"_mbin0_iter1_GaussBlur";    
    if( ipart == 5 || ipart == 6)   
      acpcorrName ="h2ypxrndEff_108Sn_"+ncls[0].Name+"_mbin0_iter1_GaussBlur"; 
  }
  else if( grname == "h_dndyx") {
    h2dname = "hyyxrnd";
    acpcorrName ="h2yyxEff_108Sn_"+ncls[ipart].Name+"_mbin0_iter1_GaussBlur";
    if( ipart == 5 || ipart == 6)   
      acpcorrName ="h2yyxEff_108Sn_"+ncls[0].Name+"_mbin0_iter1_GaussBlur"; //<<<-- need to  be considered
  }
  else if( grname == "h_dndbtgm" ){
    h2dname = "hybtgm";
    acpcorrName = "h2ybgEff_108Sn_"+ncls[ipart].Name+"_mbin0_iter1_GaussBlur";
  }
  else if( grname == "h_dndEt" ){
    h2dname = "hyEt";
    acpcorrName = "h2yetEff_108Sn_"+ncls[ipart].Name+"_mbin0_iter1_GaussBlur";
  }
  else if( grname == "h_dndpt" ){
    h2dname = "hypt";
    acpcorrName ="h2PtYEff_108Sn_"+ncls[ipart].Name+"_mbin0_iter1_GaussBlur";    
  }
  

  if( h2dname != "" ) {
    TH2D* hacp = (TH2D*)fOpen->Get(h2dname);
    
    if( hacp == NULL ) {
      LOG(ERROR) << h2dname<< " is not found " << FairLogger::endl;
      return NULL;
    }    
    LOG(INFO) << hacp->GetName() << " is found " << FairLogger::endl;
    
    TH2D* hcorrected;
    if(  fileEffCor[gname.mtconfig] ) {
      TH2D* hacpcorr = (TH2D*)fileEffCor[gname.mtconfig]->Get(acpcorrName);
      
      if( hacpcorr == NULL ) {
	LOG(ERROR) << "Efficiency correction : " << acpcorrName << " is not found " << FairLogger::endl;
	hcorrected = hacp;
      }
      else {
	LOG(INFO) << hacpcorr->GetName() << " is loaded." << FairLogger::endl;
	
	TH2D* h_smear = (TH2D*)ImageProcessing::GaussianBlur(hacp, 3,1.);
	hcorrected = GetAcpCorrection(*h_smear, hacpcorr);
	
	if( hcorrected ) 
	  LOG(INFO) << h2dname << ": Efficiency correciton have been done." << FairLogger::endl;      
        else {
          LOG(ERROR) << h2dname << ": Efficiency correciton is failed(A)." << FairLogger::endl;
          hcorrected = hacp;
        }
     }
    }
    else {
      LOG(ERROR) << h2dname << ": Efficiency correciton is failed(B)." << FairLogger::endl;
      hcorrected = hacp;
    }


    Double_t lbin = hcorrected->GetXaxis()->FindBin(-1.*ycut[ylim]);
    Double_t ubin = hcorrected->GetXaxis()->FindBin(ycut[ylim]);

    hst -> Reset();
    hst -> SetName(grname+"_"+rsys[isys]+"_"+ncls[ipart].sName+"_DATA");
    hst = hcorrected->ProjectionY(hst->GetName(),lbin,ubin);
    Double_t rnorm = hcorrected->GetXaxis()->GetBinWidth(0);
    hst -> Scale(rnorm/(2.*ycut[ylim]));
    hst -> SetDirectory(0);

    return hst;
  }


  ////---> corrected <px>
  else if( grname == "hypx" ) {
    Double_t rpresall[2] ={ 1.,0.};
    auto hdphi0_180 = (TH1D*)fOpen->Get("hdphi0_180");
    auto hdphi90_180 = (TH1D*)fOpen->Get("hdphi90_180");
    
    if( hdphi0_180 && hdphi90_180 ) 
      GetRPResolutionwChi(rpresall, hdphi0_180, hdphi90_180, 1.);

    TH2D* hypx = (TH2D*)fOpen->Get("hypx");
    if( hypx == NULL ) {
      LOG(ERROR) << "hypx  is not found " << FairLogger::endl;
      return NULL;
    }
    LOG(INFO) << hypx->GetName() << " is found " << FairLogger::endl;

    TString acpcorrName ="h2ypxrndEff_108Sn_"+ncls[ipart].Name+"_mbin0_iter1_GaussBlur";
    if( ipart == 5 || ipart == 6 )
      acpcorrName ="h2ypxrndEff_108Sn_"+ncls[0].Name+"_mbin0_iter1_GaussBlur";
    TH2D* hacpcorr = (TH2D*)fileEffCor[gname.mtconfig]->Get(acpcorrName);
    
    if( hacpcorr == NULL ) {
      LOG(ERROR) << acpcorrName << " is not found " << FairLogger::endl;
    }
    else {
      LOG(INFO) << hacpcorr->GetName() << " is loaded." << FairLogger::endl;
      TH2D* hcorrected = GetAcpCorrection(*hypx, hacpcorr);
      hypx = hcorrected;
      LOG(INFO) << hypx->GetName() << " is corrected." << FairLogger::endl;
    }

    hypx->RebinX(2);
    grv = new TGraphErrors();

    for(auto iy : ROOT::TSeqI(hypx->GetXaxis()->GetNbins()) ) {
      if( iy == 0 ) continue;
      
      auto rap = hypx->GetXaxis()->GetBinCenter(iy);
      auto hpx = hypx->ProjectionY("hpx",iy,iy);
      auto mean = hpx->GetMean();;
      mean /= rpresall[0];
      auto meanerr = hpx->GetMeanError();
      meanerr = GetError(mean, rpresall[0], meanerr, rpresall[1]);
      
      if( meanerr > 0 ) {
	grv -> SetPoint     (grv->GetN(), rap, mean);
	grv -> SetPointError(grv->GetN()-1, 0, meanerr );
      }
    }
    
    delete hacpcorr;
    delete hypx;

    return grv;
  }
  else if( grname == "gyv1" || grname == "gyv1A" ){
    Double_t rpresall[2] ={ 1.,0.};
    auto hdphi0_180 = (TH1D*)fOpen->Get("hdphi0_180");
    auto hdphi90_180 = (TH1D*)fOpen->Get("hdphi90_180");
    
    if( hdphi0_180 && hdphi90_180 ) 
      GetRPResolutionwChi(rpresall, hdphi0_180, hdphi90_180, 1.);

    TH2D* hyv1 = (TH2D*)fOpen->Get("hycos1");
    if( hyv1 ) {
      TGraphErrors* gyv1 = new TGraphErrors();
      UInt_t nbin = 5;
      hyv1->RebinX(nbin);
      TH1D* hyprj;
      for( auto ix : ROOT::TSeqI(hyv1->GetXaxis()->GetNbins() ) ){
	if( ix == 0 ) continue;

	Double_t rap = hyv1 -> GetXaxis()->GetBinCenter(ix);
	Double_t rape= hyv1 -> GetXaxis()->GetBinWidth(ix)/sqrt(12);
	hyprj = hyv1->ProjectionY(Form("hyv1_%d",ix), ix, ix);
	auto v1 = hyprj -> GetMean();

	v1 /= rpresall[0];
	auto v1error = hyprj -> GetMeanError();
	v1error = GetError(v1, rpresall[0], v1error, rpresall[1]);

	if( grname == "gyv1A" ) {
	  v1 /= ncls[ipart].A;
	  v1error /= ncls[ipart].A;
	}

	gyv1->SetPoint     ( gyv1->GetN()  ,  rap, v1 );
	gyv1->SetPointError( gyv1->GetN()-1, rape, v1error );
      }
      return gyv1;
    }
    delete hyv1;
  }
  else if( grname == "gyv2" || grname == "gyv2A" ){
    Double_t rpresall[2] ={ 1.,0.};
    auto hdphi0_180 = (TH1D*)fOpen->Get("hdphi0_180");
    auto hdphi90_180 = (TH1D*)fOpen->Get("hdphi90_180");
    
    if( hdphi0_180 && hdphi90_180 ) 
      GetRPResolutionwChi(rpresall, hdphi0_180, hdphi90_180, 2.);

    TH2D* hyv2 = (TH2D*)fOpen->Get("hycos2");
    if( hyv2 ) {
      TGraphErrors* gyv2 = new TGraphErrors();
      UInt_t nbin = 10;
      if( ipart == 4 )nbin = 5;
      hyv2->RebinX(nbin);
      TH1D* hyprj;
      for( auto ix : ROOT::TSeqI( hyv2->GetXaxis()->GetNbins() ) ){
	if( ix == 0 ) continue;

	Double_t rap = hyv2 -> GetXaxis()->GetBinCenter(ix);
       	Double_t rape= hyv2 -> GetXaxis()->GetBinWidth(ix)/sqrt(12);
	hyprj = hyv2->ProjectionY(Form("hyv2_%d",ix), ix, ix);
	auto v2 = hyprj -> GetMean();
	v2 /= rpresall[0];
	auto v2error = hyprj -> GetMeanError();//hyprj -> GetStdDev()/sqrt(hyprj -> GetEntries();
	v2error = GetError(v2, rpresall[0], v2error, rpresall[1]);

	if( grname == "gyv2A" ) {
	  v2 /= ncls[ipart].A;
	  v2error /= ncls[ipart].A;
	}
	
	gyv2->SetPoint     ( gyv2->GetN()  ,  rap, v2 );
	gyv2->SetPointError( gyv2->GetN()-1, rape, v2error );
      }
      return gyv2;
    }
    delete hyv2;
  }
  else if( grname=="gyv2ave"){
    Double_t rpresall[2] ={ 1.,0.};
    auto hdphi0_180 = (TH1D*)fOpen->Get("hdphi0_180");
    auto hdphi90_180 = (TH1D*)fOpen->Get("hdphi90_180");
    
    if( hdphi0_180 && hdphi90_180 ) 
      GetRPResolutionwChi(rpresall, hdphi0_180, hdphi90_180, 2.);

    TH2D* hyv2 = (TH2D*)fOpen->Get("hycos2am");
    TH2D* hyanum = (TH2D*)fOpen->Get("hyanum");
    if( hyv2 && hyanum ) {
      TGraphErrors* gyv2 = new TGraphErrors();
      UInt_t nbin = 10;
      
      hyv2->RebinX(nbin);
      hyanum->RebinX(nbin);

      TH1D* hyprj;
      TH1D* hynprj;
      for( auto ix : ROOT::TSeqI( hyv2->GetXaxis()->GetNbins() ) ){
	if( ix == 0 ) continue;

	Double_t rap = hyv2 -> GetXaxis()->GetBinCenter(ix);
       	Double_t rape= hyv2 -> GetXaxis()->GetBinWidth(ix)/sqrt(12);
	hyprj = hyv2->ProjectionY(Form("hyv2_%d",ix), ix, ix);
	auto v2 = hyprj -> GetMean();
	//rpcorre
	v2 /= rpresall[0];
	auto v2error = hyprj -> GetMeanError();//hyprj -> GetStdDev()/sqrt(hyprj -> GetEntries();
	v2error = GetError(v2, rpresall[0], v2error, rpresall[1]);
	
	hynprj = hyanum->ProjectionY(Form("hyanum_%d",ix), ix, ix);
	auto anum  = hynprj -> GetMean();
	auto anume = hynprj -> GetMeanError();
	
	v2 /= anum;
	v2error = GetError(v2, anum, v2error, anume);

	gyv2->SetPoint     ( gyv2->GetN()  ,  rap, v2 );
	gyv2->SetPointError( gyv2->GetN()-1, rape, v2error );
      }
      return gyv2;
    }
  }
  else if( grname == "gu_v1" ){
    TGraphErrors* gv1 = (TGraphErrors*)fOpen->Get("gu_v1");
    if( gv1 )
      return gv1;
  } 
  else if( grname == "gu_v2" ){
    TGraphErrors* gv2 = (TGraphErrors*)fOpen->Get("gu_v2");
    if( gv2 )
      return gv2;
  } 
  else if( grname == "gy_v1" ){
    TGraphErrors* gv1 = (TGraphErrors*)fOpen->Get("gy_v1");
    if( gv1 )
      return gv1;
  } 
  else if( grname == "gy_v2" ){
    TGraphErrors* gv2 = (TGraphErrors*)fOpen->Get("gy_v2");
    if( gv2 )
      return gv2;
  } 


  return NULL;
}

TGraphErrors* LoadFOPI(UInt_t ipart=0, TString fdir="PRC89/Fig8_v1Ut_0.25", TString sysname="" )
{
  TString FOPI_dir =  "data/FOPI/";

  FOPI_dir += fdir + "_" + ncls[ipart].name + "_" + sysname + ".txt";

  TGraphErrors* grph = LoadTextGraph(FOPI_dir); 
  grph -> SetName("g_utv1_FOPI_" +  ncls[ipart].sName + sysname );

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
TObject* LoadAMDRoot(UInt_t isys, UInt_t ipart, mplot mdata, TString grname="h_dndy", UInt_t icent=1, Int_t ylim=-1.)
{
  TFile *ifile;
  if( !mdata.fname.Contains(Form("%d",Asys[isys])) ) return NULL;

  if( fileAMD[mdata.id][icent] == NULL ) {
    TString dirname  = "ModelData/AMD/" + mdata.path;  
    TString tfile = mdata.fname + "_2to0fm.root";
    if( icent == 0 )
      tfile = mdata.fname + "_4to2fm.root";
    else if( icent == 3 )
      tfile = mdata.fname + "_152to200M.root";
    else if( icent == 4 )
      tfile = mdata.fname + "_0to152M.root";

    LOG(INFO) << "[LoadAMD] --> " << grname << " " << dirname << tfile  << FairLogger::endl;

    if( gSystem -> FindFile(dirname, tfile ) ) {
      fileAMD[mdata.id][icent] = TFile::Open(tfile,"READ");
      ifile = fileAMD[mdata.id][icent];
      LOG(INFO) << ifile->GetName() << FairLogger::endl;
      if( !ifile ) return NULL;
    }
    else{
      LOG(ERROR) << tfile << " is not found." << FairLogger::endl; 
      return NULL;
    }
  }
  else
    ifile = fileAMD[mdata.id][icent];


  if( grname == "gy_px" || grname == "g_v1y"  || grname == "g_v2y" || grname == "gyv1" || grname == "gyv2" ||
      grname == "gyv1A" || grname == "gyv2A" )  {

    TString pgrname = grname;
    if( grname ==  "gyv1A" || grname == "gyv2A" ) 
      pgrname = grname(0,4) +  "_"+ncls[ipart].sName;
    else
      pgrname += "_"+ncls[ipart].sName;
 
    LOG(INFO) << "[AMD searching] " << pgrname << " in " << ifile  << FairLogger::endl;

    auto obj = ifile->Get(pgrname);
    if( obj ) {
      LOG(INFO) << "[AMD] " << pgrname << " is found. " << FairLogger::endl;

      if( grname ==  "gyv1A" || grname == "gyv2A" ) {
	TGraphErrors* gr = (TGraphErrors*)obj;
	Double_t x, y, xe, ye;
	for( auto nx : ROOT::TSeqI(gr->GetN()) ) {
	  gr -> GetPoint( nx, x, y );
	  ye  = gr->GetErrorY( nx );
	  xe  = gr->GetErrorX( nx );
	  auto ya = y / ncls[ipart].A;
	  auto yea= ye/ ncls[ipart].A;
	  gr -> SetPoint( nx, x, ya );
	  gr -> SetPointError(nx, 0, yea );
	}
      }


      if( bReverse )
	return ReversePlot(obj);
      else
	return obj;
    }
    return  NULL;
  }

  TH1D* hist = NULL;
  TGraphErrors* grv = NULL;

  //---> Projection on X
  TString h2dname;
  if( grname == "h_dndy" ) 
    h2dname = "hypx_" + ncls[ipart].sName;

  if( h2dname != "" ) {
    TH2D* hypx = (TH2D*)ifile->Get(h2dname);
    if( hypx == NULL ) {
      LOG(ERROR) << h2dname << " is not found " << FairLogger::endl;
      return NULL;
    }
    
    LOG(INFO) << hypx->GetName() << " will be projected on X axis." << FairLogger::endl;
    hypx -> SetName(grname + "_2d_" + rsys[isys]+"_AMD");
    auto hypx_px = hypx->ProjectionX(Form("%s_px",hypx->GetName()));
    Double_t rnorm = hypx->GetYaxis()->GetBinWidth(0);
    hypx_px -> Scale(rnorm);
    hypx_px -> SetDirectory(0);

    if( bReverse ) {
      LOG(INFO) << "[AMD] " << hypx_px->GetName() << " is reversed. "  << FairLogger::endl;
      hypx_px = (TH1D*)ReversePlot(hypx_px);
    }

    return hypx_px;
  }

  //---> Projection on Y
  h2dname = "";
  if( grname == "h_dndyx" )
    //    h2dname = "hyyx_" + ncls[ipart].sName;
    h2dname = "hyyxrnd_" + ncls[ipart].sName;
  else if( grname == "h_dndpx" )
    h2dname = "hypxrnd_" + ncls[ipart].sName;
  else if( grname == "h_dndEt" )
    h2dname = "hyEt_" + ncls[ipart].sName;
  else if( grname == "h_dndbtgm" )
    h2dname = "hybtgm_" + ncls[ipart].sName;
  else if( grname == "h_dndpt" )
    h2dname = "hypt_" + ncls[ipart].sName;

  if( h2dname != "" ) {
    TH2D* hypy = (TH2D*)ifile->Get(h2dname);
    if( hypy == NULL ) {
      LOG(ERROR) << h2dname << " is not found " << FairLogger::endl;
      return NULL;
    }
    
    LOG(INFO) << hypy->GetName() << " is found and cut in |y|< " << ycut[ylim] << FairLogger::endl;
      
    Double_t lbin = hypy->GetXaxis()->FindBin(-1.*ycut[ylim]);
    Double_t ubin = hypy->GetXaxis()->FindBin(ycut[ylim]);
    hypy -> SetName(grname + "_2dy_" + rsys[isys]+"_AMD");
    auto hypy_py = hypy->ProjectionY(Form("%s_py",hypy->GetName()), lbin, ubin);
    Double_t rnorm = hypy->GetXaxis()->GetBinWidth(0);
    hypy_py -> Scale(rnorm/(2.*ycut[ylim]));
    hypy_py -> SetDirectory(0);

    return hypy_py;
  }

  if( grname == "hnuclei" ) {
    TH2D *hist = (TH2D*)ifile->Get(grname);
    return hist;
  }
  if( grname == "hb" ) {
    TH1D *hist = (TH1D*)ifile->Get(grname);
    return hist;
  }
  


  if( grname == "g_dndy" ) {
    grname = "h_dndy_"+ncls[ipart].sName;

    hist = (TH1D*)ifile->Get(grname);
    if( hist == NULL ) return NULL;

    if( bReverse ) {
      LOG(INFO) << "[AMD] reversed." << FairLogger::endl;
      hist = (TH1D*)ReversePlot(hist);
    }

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

    return grv;
  }

    

  return NULL;
}

//###LoadAMD
TObject* LoadAMD(UInt_t isys, UInt_t ipart, mplot mdata, TString grname="v1", UInt_t icent=1,  Int_t ylim = -1.)
{
  bReverse = 0;
  if( isys == 3 && (mdata.fname.Contains("Sn124Sn112") || mdata.path.Contains("Sn124Sn112")) ) bReverse = 1;
  else if( !mdata.fname.Contains(asys[isys]) ) return NULL;

  if( bAMDRoot )
    return (TObject*)LoadAMDRoot(isys, ipart, mdata, grname, icent, ylim);
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


TObject* LoadKanekoData(UInt_t isys, UInt_t ipart, TString gname="g_dndy")
{  
  TString partname[] = {"Proton","Deuteron","Triton","He3","He4"};

  auto *ofile = TFile::Open("data/dndy_Kaneko.root");
  
  if( ofile == NULL ) {
    LOG(ERROR) << "data/dndy_Kaneko.root is not found " << FairLogger::endl;
    return NULL;
  }

  TString hname = "h1dndy_" + bName[isys] + partname[ipart];
  TH1D* hist = (TH1D*)ofile->Get(hname);

  if( hist == NULL )
    return NULL;

  hist->SetDirectory(gROOT);

  if( gname.Contains("h_") && hist ) {
    LOG(INFO) << hist->GetName() << FairLogger::endl;
    ofile->Close();
    return hist;
  }


  LOG(INFO) << " LoadKanekoData " << hname << " is found."<< FairLogger::endl;

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
  
  ofile->Close();
  return gvr;
}

TObject* LoadTommyData(UInt_t isys, UInt_t ipart, TString gtype)
{
  TString partname[] = {"pim","pip","p","d","t","He3","He4"};
  TString sysname[]  = {"_Sn132","_Sn108","_Sn124","_Sn112"};

  TGraphErrors* grv = NULL;
  
  TString fname = "Content_Mult_" + gtype  + sysname[isys] + ".txt";
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
  
  TString search = partname[ipart+2]+"_"+gtype + sysname[isys]; 
  LOG(INFO) << " seraching for " << search << FairLogger::endl;

  Bool_t beof = kTRUE;
  while( beof ) {
    fread >> sget;
    if( sget == search ) {
      fread >> sget; fread >> sget; fread >> sget;

      while( !fread.eof() ) {
	fread >> sget;
	if( sget.Contains( gtype ) || fread.eof()) {
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
  LOG(INFO) << " seraching for " << search << FairLogger::endl;

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
  TString fname = gnames[igname].fileHeader + bName[isys] + ncls[0].name + gnames[igname].Version + ".root";

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
  TString fname = gnames[igname].fileHeader + bName[isys] + ncls[ipart].name + gnames[igname].Version + ".root";

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


  gname += "_" + rsys[isys] + "_" + ncls[ipart].sName;
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

    LOG(INFO) << "x = " << x << " y= " << y << FairLogger::endl;
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

Double_t* GetSlopeParameter(TGraphErrors *grph, Double_t ft_low=-0.5, Double_t ft_high=0.8, Double_t par0=0, Double_t par1=0.)
{
  Double_t *vfit = new Double_t[4];
  for( auto i : ROOT::TSeqI(4) )
    vfit[i] = 0.;

  fv1fit->SetParameter(0,par0);
  fv1fit->SetParameter(1,par1);

  grph -> Fit("fv1fit","Q0","", ft_low, ft_high); //"Q0","");     
  
  vfit[0] = fv1fit->GetParameter(1);
  vfit[1] = fv1fit->GetParError(1); 
  if( vfit[0] < 0.1 ) {vfit[0] = 0.; vfit[1] = 0.;}
  vfit[2] = fv1fit->GetChisquare();
  vfit[3] = fv1fit->GetNDF();

  return vfit;
}

Double_t* Getv2FitParameters(TGraphErrors *grph, Double_t ft_low=-0.5, Double_t ft_high=0.8)
{
  Double_t *vfit = new Double_t[8];
  for( auto i : ROOT::TSeqI(8) )
    vfit[i] = 0.;

  fv2fit->SetParameter(0,-0.03);
  fv2fit->SetParameter(1,0.1);
  fv2fit->SetParameter(2,0.);

  grph -> Fit("fv2fit","Q0","", ft_low, ft_high); //"Q0","");     
  
  vfit[0] = fv2fit->GetParameter(0);
  vfit[1] = fv2fit->GetParError(0); 
  vfit[2] = fv2fit->GetParameter(1);
  vfit[3] = fv2fit->GetParError(1); 
  vfit[4] = fv2fit->GetParameter(2);
  vfit[5] = fv2fit->GetParError(2); 

  vfit[6] = fv1fit->GetChisquare();
  vfit[7] = fv1fit->GetNDF();

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

	LOG(INFO) << " refit v2  " << ck << ", " << lpar << " , "  << FairLogger::endl;
	
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
    //   else
      //      fitval = GetV11(igname, isys, ipart, sData, mdata);

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
    
    tlabel.DrawLatex(0.2, labely, ncls[j].sName);

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
    lg -> AddEntry(gv, ncls[j].Name);
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

    if( j == 0 )  plabel.DrawLatexNDC(0.25, 0.7/yFactor, ncls[j].sName);
    else
      plabel.DrawLatexNDC(0.25, 0.7, ncls[j].sName);

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

    plabel.SetTextAlign(13);
    plabel.SetTextSize(0.09*yFactor);

    if( j == 0 )  plabel.DrawLatexNDC(0.3, 0.5/yFactor, ncls[j].name);
    else
      plabel.DrawLatexNDC(0.3, 0.5, ncls[j].name);

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
    auto lg = new TLegend(0.2, 0.2, 0.6, 0.5, ncls[ipart].sName);
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
      lg -> SetHeader(ncls[ipart].sName+" "+header);

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

void Draw_v_y(UInt_t igname, std::vector<UInt_t> sqpart, UInt_t vn = 1) 
{
  if( vn > 2 ) vn = 1;

  Bool_t bAMD = 1;
  Bool_t bpBUU = 0;
  Bool_t bImQMD = 0;

  LOG(INFO) << "Draw_v" << vn << "_y"  << FairLogger::endl;

  const Int_t Nx = 1;//(Int_t)sqSys.size();
  const Int_t Ny = (Int_t)sqpart.size();

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 300*Ny); iccv++;
  gStyle->SetOptTitle(0);


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
      if( pad[i][j] == NULL ) {
        LOG(ERROR) << " pad is not found " << pname << FairLogger::endl;
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

      TLegend *lg = new TLegend(0.6,0.4/yFactor,0.9,0.65/yFactor,ncls[sqpart[j]].sName);
      lg -> SetFillColor(0);

      UInt_t isys = 0;
      v_y = (TGraphErrors*)LoadData(isys, sqpart[j], Form("gu_v%d",vn), gnames[igname]);
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
	LOG(ERROR) << " gu_v1  is not found. " << isys << " : " << sqpart[j] << FairLogger::endl;
      
      if( bAMD ) {
	v_y = (TGraphErrors*)LoadAMD(sqSys[i], sqpart[j], AMDnames[0], Form("v%d_y",vn));
	
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
	    v_y = LoadpBUU(sqSys[i], sqpart[j], ipbuu , "v1_yn");
	  if( vn == 2 )
	    v_y = LoadpBUU(sqSys[i], sqpart[j], ipbuu , "v2_yn");
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
	  v_y = LoadImQMD(sqSys[i], sqpart[j], "v1" );
	else if( vn == 2)
	  v_y = LoadImQMD(sqSys[i], sqpart[j], "v2" );
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

  TLegend *lg = new TLegend(0.2,0.5,0.6,0.9,ncls[ipart].sName);
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
      header = ncls[ipart].sName+":"+header(0,header.First(";")); 


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
  plabel.DrawLatexNDC(0.15,0.9, fsys[isys]+"_"+ncls[ipart].sName);


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

    for( auto igname : ROOT::TSeqI(nDATA) ) {
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

void GetModelMultiplicity(TH1D* hist, Double_t var[])
{
  var[0] = hist->GetMean();
  var[1] = hist->GetMeanError();
}

void GetIntegralandSigma(TH1D* hist, Double_t var[])
{
  TString hname = hist->GetName();
  Bool_t blim = kFALSE;

  if( hname.Contains("h_dndy_") ){
    auto x1 = hist->FindBin( 0.); 
    auto x2 = hist->FindBin( 2.); 
    var[0] = 2. * hist -> IntegralAndError(x1,x2,var[1],"width");
    blim = kTRUE;
    cout << " hist --> " << hist->GetName() << " " << var[0] << endl;

  }    
  else if( hname.Contains("h_dndyx_") ){
    var[4] = hist->IntegralAndError(-1,-1,var[5],"width");
    var[6] = hist->GetStdDev();
    var[7] = hist->GetStdDevError();
    blim = kTRUE;
  }
  else if( hname.Contains("h_dndpx_") ){
    var[8] = hist->GetStdDev();
    var[9] = hist->GetStdDevError();    
  }
  else if( hname.Contains("h_dndEt")){
    var[10] = hist->GetMean();
    var[11] = hist->GetMeanError();
  }
  else if( hname.Contains("h_dndbtgm")){
    var[12] = hist->GetMean();
    var[13] = hist->GetMeanError();
    var[14] = hist->IntegralAndError(-1,-1,var[15],"width");
  }


  if( blim ) {
    TH1D* h_dndy_lim = (TH1D*)hist->Clone();
    auto x1 = h_dndy_lim->GetXaxis()->FindBin(-1.);
    auto x2 = h_dndy_lim->GetXaxis()->FindBin( 1.);
    auto mid= h_dndy_lim->GetXaxis()->FindBin( 0.);
    
    UInt_t nbin = h_dndy_lim->GetXaxis()->GetNbins();
    
    for( auto ibin : ROOT::TSeqI( nbin ) ) { // reflect y>0 to y<0
      if ( ibin < x1 || ibin > x2 ) {
	h_dndy_lim->SetBinContent(ibin,0.);
	h_dndy_lim->SetBinError  (ibin,0.);
      }
      
      if( ibin > mid ) {
	auto rap  = h_dndy_lim -> GetBinCenter(ibin);
	auto cont = h_dndy_lim -> GetBinContent(ibin);
	auto conte= h_dndy_lim -> GetBinError(ibin);
	auto rbin = h_dndy_lim -> GetXaxis()->FindBin( -rap );
	h_dndy_lim->SetBinContent( rbin, cont );
	h_dndy_lim->SetBinError( rbin, conte );
      }
    }
  
    var[2] = h_dndy_lim->GetStdDev();
    var[3] = h_dndy_lim->GetStdDevError();
    delete h_dndy_lim;
  }
}


void Getvarxz(UInt_t isys, UInt_t ipart, std::vector<gplot> gname, Double_t var[], Int_t ylim=-1) 
{
  gStyle -> SetLegendFillColor(0);
  
  THStack*  mhst;
  if( mhst ) mhst = new THStack("hs","");
  
  auto lg  = new TLegend(0.15, 0.2, 0.5, 0.7,ncls[ipart].sName);
  lg -> SetFillStyle(0);

  TH1D *h_dndy = NULL;
  TH1D *h_dndyx = NULL;

  for( auto gdata : gname ){
    h_dndy  = (TH1D*)LoadData(isys,ipart,"h_dndy" ,gdata,-1); 
    h_dndyx = (TH1D*)LoadData(isys,ipart,"h_dndyx",gdata, ylim); 

    if( h_dndy != NULL && h_dndyx != NULL) {

      Double_t vv[2];
      Double_t vve[2];
      Double_t mtot[3];
      Double_t err[3];

      TH1D* hists[] = {h_dndy,h_dndyx};     
      mhst -> SetTitle(";y_{cm}/y_{prj};dN/dy(dN/dyx)");

      UInt_t i = 0;
      for( auto hst : hists ) {
	if( i == 0 ) {
	  auto x1 = hst->FindBin(-ycut[ylim]);
	  auto x2 = hst->FindBin( ycut[ylim]);
	  mtot[2] = hst->IntegralAndError(x1,x2,err[2],"width");//*hst->GetXaxis()->GetBinWidth(0);
	}
	
	mtot[i] = hst->IntegralAndError(-1,-1,err[i],"width");//*hst->GetXaxis()->GetBinWidth(0);

	TH1D* h_dndy_lim = (TH1D*)hst->Clone();
	auto x1 = h_dndy_lim->GetXaxis()->FindBin(-1.);
	auto x2 = h_dndy_lim->GetXaxis()->FindBin( 1.);
	auto mid= h_dndy_lim->GetXaxis()->FindBin( 0.);
	
	UInt_t nbin = h_dndy_lim->GetXaxis()->GetNbins();

	for( auto ibin : ROOT::TSeqI( nbin ) ) { // reflect y>0 to y<0
	  if ( ibin < x1 || ibin > x2 ) {
	    h_dndy_lim->SetBinContent(ibin,0.);
	    h_dndy_lim->SetBinError(ibin,0.);
	  }

	  if( hst == hists[0] && ibin > mid ) {
	    auto rap  = h_dndy_lim -> GetBinCenter(ibin);
	    auto cont = h_dndy_lim -> GetBinContent(ibin);
	    auto conte= h_dndy_lim -> GetBinError(ibin);
	    auto rbin = h_dndy_lim -> GetXaxis()->FindBin( -rap );
	    h_dndy_lim->SetBinContent( rbin, cont );
	    h_dndy_lim->SetBinError( rbin, conte );
	  }
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
	vve[i] = h_dndy_lim->GetStdDevError();
	i++;
      }

      var[0] = vv[1]*vv[1] / (vv[0]*vv[0]);
      var[1] = 2.* abs(vve[1]/vv[1]) +  2.* abs(vve[0]/vv[0]);
      var[2] = vv[0];
      var[3] = vve[0];
      var[4] = vv[1];
      var[5] = vve[1];
      var[6] = mtot[0];
      var[7] = err[0];
      var[8] = mtot[1];
      var[9] = err[1];

      LOG(INFO) << " varxz = : " << var[0] << " +- " << var[1]  << FairLogger::endl;
      LOG(INFO) << " Integral in Longitudinal -> " << mtot[0]
		<< " Long in " << -ycut[ylim] << " to " << ycut[ylim] << " -> " << mtot[2]  
		<< " Transverse -> " << mtot[1] << FairLogger::endl;
    }
      
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
    ccv -> SetGrid(1);
    // //    mhst -> GetXaxis() -> SetRangeUser(-2.5, 2.5);
    mhst -> Draw("nostack");
    lg->Draw();
  }
  
  return;
}


void Draw_varxz(UInt_t isys, std::vector<gplot> gname, Int_t ylim=2)
{
  Bool_t bAMD = 1;


  Double_t* var = new Double_t[10];

  TString ypara[] = {"varxz","varz","varx","multz","multx"};
  TGraphErrors* g_var[5];
  for( auto i : ROOT::TSeqI(5) ) {
    g_var[i] = new TGraphErrors();
    g_var[i] -> SetTitle(";PID;"+ypara[i]);
  }

  UInt_t in = 0;
  for( UInt_t ipart : {0,1,2,3,4} ){
    Getvarxz(isys, ipart, gname, var, ylim);

    if( var[0] != 0 ) {
      for( auto i : {0,1,2,3,4} ) {
	g_var[i]->SetPoint( in, ipart+1, var[2*i]);
	g_var[i]->SetPointError( in, 0,  var[2*i+1]);
      }
      in++;
    }
  }
 
  TMultiGraph* m_varxz;
  TLegend* legend;
 
  for( auto i : {0,1,2,3,4} ) {    
    m_varxz = new TMultiGraph();
    legend = new TLegend(0.2,0.5,0.9,0.9,"");

    UInt_t loop = 1;
    //    if( i == 1 || i == 3 ) loop = 2;
    for( auto j : ROOT::TSeqI(loop) ) {
      if( g_var[i+j] -> GetN() > 0 ) {
	g_var[i+j] -> SetMarkerStyle(20+j);
	g_var[i+j] -> SetMarkerColor(2+j);
	g_var[i+j] -> SetLineColor(2+j);

	if( j == 20 ) {
	  m_varxz->Add(g_var[i+j],"alp");
	  for(auto ipart : {0,1,2,3,4} )
	    g_var[i]->GetXaxis()->SetBinLabel( g_var[i]->GetXaxis()->FindBin( ipart+1 ), ncls[ipart].sName);
	}
	else
	  m_varxz->Add(g_var[i+j],"lp");

	m_varxz -> SetTitle( g_var[i+j] -> GetTitle() );

	legend ->AddEntry(g_var[i+j],"DATA"+fsys[isys]+" "+gname[0].Version);
      }
      //      g_var[i]->GetYaxis()->SetRangeUser(0.3,1.15);
    }
  
    if( bAMD ){
      for( auto samd : AMDnames ) {

	loop = 1;
	//	if( i == 1 || i == 3 ) loop = 2;
	for( auto j : ROOT::TSeqI(loop) ) {
	  g_var[i+j] = (TGraphErrors*)LoadAMD(isys,0,samd, "g_"+ypara[i+j]+"1");
	  if( g_var[i+j] ) {
	    g_var[i+j] -> SetFillStyle( samd.fStyle );
	    g_var[i+j] -> SetLineColor( samd.fColor+j );
	    g_var[i+j] -> SetFillColorAlpha( samd.fColor+j, 0.4 );
	    m_varxz  ->Add(g_var[i+j],"3");
	    legend ->AddEntry(g_var[i+j],samd.config);
	  }
	}
      }
    }
  
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
    auto gtmp = new TGraph();
    gtmp -> SetPoint(0,0,0.2);
    gtmp -> SetPoint(1,6,0.2);
    m_varxz->Add(gtmp,"");
    for(auto ipart : {0,1,2,3,4} ) {
      auto fbin = m_varxz -> GetXaxis() -> FindBin(ipart+1);
      if( fbin != 101 )
	m_varxz->GetXaxis()->SetBinLabel( fbin, ncls[ipart].sName );
    }
    m_varxz -> Draw("ALP");
    // legend -> Draw();
  }

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  legend -> Draw();

}

void Draw_meanpx(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, std::vector<gplot> gname, TString grname="hypx", UInt_t icateg=4, Bool_t bsaveData=0)
{
  LOG(INFO) << "[Draw_meanpx] .......... "<< fsys[ivsys.at(0)] 
	    << " bsaveData "   << bsaveData 
	    << " bCentrality " << bCentral << FairLogger::endl;

  gStyle->SetOptFit(0);
  Bool_t bAMD  = 1;

  if( ivsys.size() > 1 && bsaveData == 0 ) bAMD = 0;
  //  if( bsaveData ) grname = "hypx";

  const UInt_t npart = (UInt_t)sqpart.size(); 
  TMultiGraph* mgr[npart];
  for( auto ipart : sqpart ) {
    mgr[ipart] = new TMultiGraph();
    if( grname == "hypx" )
      mgr[ipart] -> SetTitle(";y_{cm}/y_{prj};<px>/A");
    else if( grname == "gyv1" )
      mgr[ipart] -> SetTitle(";y_{cm}/y_{prj}; v1");
    else if( grname == "gyv2" )
      mgr[ipart] -> SetTitle(";y_{cm}/y_{prj}; v2");
    else if( grname == "gyv1A" )
      mgr[ipart] -> SetTitle(";y_{cm}/y_{prj}; v1/A");
    else if( grname == "gyv2A" )
      mgr[ipart] -> SetTitle(";y_{cm}/y_{prj}; v2/A");

  }

  auto lg  = new TLegend(0.15, 0.2, 0.5, 0.7,"");
  lg -> SetFillStyle(0);

  TGraphErrors *grph = NULL;


  for( auto iisys : ivsys ) {
    for( auto ipart : sqpart) {


      // AMD
      if( bAMD ) {

	for( auto ic : {0,1} ) { // midcentral andl central

	  if( !bsaveData and bCentral != ic) continue;

	  for( auto samd : AMDnames ) {
	  
	    if( !bsaveData && samd.category != icateg ) continue;

	    TString agrname = grname;
	    if( grname == "hypx" ) agrname = "gy_px";
	    else if( grname == "gu_v1" || grname == "gy_v1" ) agrname = "g_v1y";
	    else if( grname == "gu_v2" || grname == "gy_v2" ) agrname = "g_v2y";

	    grph = (TGraphErrors*)LoadAMD(iisys, ipart, samd, agrname, ic);

 	    if( grph != NULL ) {
	      UInt_t iin = samd.id;

	      grph -> SetLineColor( samd.fColor);
	      grph -> SetFillColorAlpha( samd.fColor, 0.35);
	      grph -> SetFillStyle( 3001 );
	      grph -> SetLineWidth( 5002);

	      if(!bsaveData && grname == "gyv1") {
		GetSlopeParameter(grph, -0.5, 0.5,0,0.4);
	      }
	      else if(!bsaveData &&  grname == "gyv2") {
		Getv2FitParameters(grph, -1.0, 1.0);
	      }

	      
	      mgr[ipart] -> Add( grph, "3E");
	      if( ipart == 0 )
		//		lg  -> AddEntry( grph, "AMD "+ samd.config);
		lg  -> AddEntry( grph, samd.fullconfig);
     

	      if( bsaveData ) {
		PhysParameter *sAMD = &physAMD[iisys][ipart][iin][ic];
		sAMD->version = samd.config;
		sAMD->system  = iisys;
		sAMD->centrality = ic;
		
		if( grname == "hypx" ) {
		  Double_t* para = GetSlopeParameter(grph, fphysdatafile[fphysdataid].first, 0.5 );
		  sAMD->meanpxSlope = *(para);
		  sAMD->meanpxSlopeError = *(para+1);
		}
		else if( grname == "gu_v1" || grname == "gy_v1" || grname == "gyv1") {
		  Double_t* para = GetSlopeParameter(grph, fphysdatafile[fphysdataid].first, 0.5, 0., 0.4 );
		  sAMD->v1Slope = *(para);
		  sAMD->v1SlopeError = *(para+1);
		}
		else if( grname == "gyv2" ) {
		  Double_t* para = Getv2FitParameters(grph, -1.0, 1.0);		  
		  sAMD->v2minimum     =*(para);
		  sAMD->v2minimumError=*(para+1);
		  sAMD->v2width       =*(para+2);
		  sAMD->v2widthError  =*(para+3);
		  sAMD->v2offset      =*(para+4);
		  sAMD->v2offsetError =*(para+5);
		}
	      }
	    }
	  }
	}	
      }

      // Data
      UInt_t ig = 0;
      for( auto igname : gname ) {
	if( !bsaveData && igname.centrality != bCentral ) continue; 

	grph =(TGraphErrors*)LoadData(iisys,ipart,grname,igname); // acceptance corrected
	
	if( grph != NULL ) {
	  grph -> SetMarkerStyle(DStyle[iisys].mStyle+ig);
	  grph -> SetMarkerSize( DStyle[iisys].mSize);
	  grph -> SetLineColor(  DStyle[iisys].fColor+2*ig );
	  grph -> SetMarkerColor(DStyle[iisys].fColor+2*ig );

	  if( 0 ) {
	    if(!bsaveData && grname == "hypx" ) {
	      auto ffit = (TF1*)gROOT->FindObject("fv1fit");
	      if( ffit ) {
		gStyle->SetOptFit(0);
		ffit->SetLineColor(DStyle[iisys].fColor );
	      }
	      GetSlopeParameter(grph, fphysdatafile[fphysdataid].first, 0.5,0.,0.6);
	    }
	    else if(!bsaveData && grname == "gyv1" ) {
	      auto ffit = (TF1*)gROOT->FindObject("fv1fit");
	      if( ffit ) {
		gStyle->SetOptFit(0);
		ffit->SetLineColor(DStyle[iisys].fColor );
	      }
	      GetSlopeParameter(grph, fphysdatafile[fphysdataid].first, 0.5,0.,0.4);
	    }
	    else if(!bsaveData && grname == "gyv2" ) {
	      auto ffit = (TF1*)gROOT->FindObject("fv2fit");
	      if( ffit ) {
		gStyle->SetOptFit(0);
		ffit->SetLineColor(DStyle[iisys].fColor );
	      }
	      Getv2FitParameters(grph, -1., 1.);
	    }
	  }

	  mgr[ipart] -> Add( grph, "P");
	  if( ipart == 0 )
	    // lg  -> AddEntry( grph, fsys[iisys]+" "+igname.config1);
	    lg  -> AddEntry( grph, fsys[iisys]+" "+igname.Version);
	    
	  //@@@ 
	  if( bsaveData ) {	    

	    PhysParameter *sData = &physData[iisys][ipart][ig];
	    if( sData ) {
	      sData->version=igname.Version;
	      sData->system=iisys;
	      sData->centrality=ig;
	      sData->partID=ipart;
	      if( grname == "hypx" ) {
		Double_t* para = GetSlopeParameter(grph, fphysdatafile[fphysdataid].first, 0.5);
		sData->meanpxSlope=*(para);
		sData->meanpxSlopeError = *(para+1);
	      }
	      else if( grname == "gu_v1" || grname == "gy_v1" || grname == "gyv1") {
		Double_t* para = GetSlopeParameter(grph, fphysdatafile[fphysdataid].first, 0.5,0,0.6);
		sData->v1Slope=*(para);
		sData->v1SlopeError=*(para+1);
		cout << " v1 slope " << sData->v1Slope << " +- " << sData->v1SlopeError << endl;
	      }
	      else if( grname == "gyv2" ) {
		Double_t* para = Getv2FitParameters(grph, -1.0, 1.0);
		sData->v2minimum     =*(para);
		sData->v2minimumError=*(para+1);
		sData->v2width       =*(para+2);
		sData->v2widthError  =*(para+3);
		sData->v2offset      =*(para+4);
		sData->v2offsetError =*(para+5);
		cout << " v20 " << sData->v2minimum << " +- " << sData->v2minimumError << endl;
	      }
	    }
	  }
	}
	ig++;
      }
      //-- endof data

    }
  }

  // auto glist = mgr[0]->GetListOfGraphs();
  // glist->Print();


  if( grname == "gyv2ave" ){
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
    auto pad = ccv->cd(1);
    pad -> SetGrid(1);

    mgr[6] -> GetYaxis()->SetRangeUser(-0.05,0.05);
    mgr[6] -> GetXaxis()->SetLimits(-0.9, 1.3);

    mgr[6] -> GetXaxis()->SetNdivisions(5,5,0,kTRUE);
    mgr[6] -> GetXaxis()->SetNdivisions(5,5,0,kTRUE);


    mgr[6] -> GetYaxis()->SetTitleFont(43);
    mgr[6] -> GetYaxis()->SetTitleSize(18);
    mgr[6] -> GetYaxis()->SetTitleOffset(2.8);
    mgr[6] -> GetYaxis()->SetLabelFont(43);
    mgr[6] -> GetYaxis()->SetLabelSize(20);

    mgr[6] -> GetXaxis()->SetTitleFont(43);
    mgr[6] -> GetXaxis()->SetTitleSize(18);
    mgr[6] -> GetXaxis()->SetTitleOffset(1.8);
    mgr[6] -> GetXaxis()->SetLabelFont(43);
    mgr[6] -> GetXaxis()->SetLabelSize(20);


    mgr[6] -> SetTitle(";y/y_{b}; v2");
    mgr[6] -> Draw("AP");    
    plabel.SetTextSize(0.08);
    plabel.DrawLatexNDC(0.55,0.2, "averaged H");

    mgr[6]->Print();

  }

  else if( !bsaveData ) {  
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),250*npart,450); iccv++;
    //    ccv->Divide(npart+1,1);

    const Int_t Nx = npart+1;

    Float_t lMargin = 0.041;
    Float_t rMargin = 0.05;
    Float_t bMargin = 0.05;
    Float_t tMargin = 0.05;

    CanvasPartitionX(ccv,Nx,lMargin,rMargin,bMargin,tMargin);

    TPad* pad[Nx];

    std::vector< UInt_t >::iterator ipx = sqpart.begin();      
    for( auto i : ROOT::TSeqI(Nx) ) {

      ccv->cd(0);
      char pname[16];
      sprintf(pname,"pad_%i",i);
      pad[i] = (TPad*)gROOT->FindObject(pname);
      pad[i] -> Draw();
      pad[i] -> SetFillStyle(4000);
      pad[i] -> SetFrameFillStyle(400);
      pad[i] -> SetGrid(1);
      pad[i] -> cd();

      Float_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
      Float_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

      if( i < npart ) {
	mgr[*ipx] -> GetXaxis()->SetLimits(-1.,1.5);

	if( grname == "gu_v1"      || grname == "gy_v1" || grname == "gyv1" ) {
	  mgr[*ipx] -> GetYaxis()->SetRangeUser(-0.78,0.82);
	  mgr[*ipx] -> GetXaxis()->SetLimits(-0.9, 1.3);
	}
	else if( grname == "gu_v2" || grname == "gy_v2" || grname  == "gyv2" ) {
	  mgr[*ipx] -> GetYaxis()->SetRangeUser(-0.085,0.075);
	  mgr[*ipx] -> GetXaxis()->SetLimits(-0.9, 1.3);
	}
	else if( grname == "gyv1A" ){
	  mgr[*ipx] -> GetYaxis()->SetRangeUser(-0.4,0.4);
	  mgr[*ipx] -> GetXaxis()->SetLimits(-0.9, 1.3);
	}
	else if( grname == "gyv2A" ){
	  mgr[*ipx] -> GetYaxis()->SetRangeUser(-0.04,0.04);
	  mgr[*ipx] -> GetXaxis()->SetLimits(-0.9, 1.3);
	}
	else
	  mgr[*ipx] -> GetYaxis()->SetRangeUser(-100.,100.);	  

	mgr[*ipx] -> GetXaxis()->SetNdivisions(5,5,0,kTRUE);
	mgr[*ipx] -> GetXaxis()->SetNdivisions(5,5,0,kTRUE);


	mgr[*ipx] -> GetYaxis()->SetTitleFont(43);
	mgr[*ipx] -> GetYaxis()->SetTitleSize(18);
	mgr[*ipx] -> GetYaxis()->SetTitleOffset(2.8);
	mgr[*ipx] -> GetYaxis()->SetLabelFont(43);
	mgr[*ipx] -> GetYaxis()->SetLabelSize(20);

	mgr[*ipx] -> GetXaxis()->SetTitleFont(43);
	mgr[*ipx] -> GetXaxis()->SetTitleSize(18);
	mgr[*ipx] -> GetXaxis()->SetTitleOffset(1.8);
	mgr[*ipx] -> GetXaxis()->SetLabelFont(43);
	mgr[*ipx] -> GetXaxis()->SetLabelSize(20);

	mgr[*ipx] -> Draw("AP");
	//	mgr[*ipx] -> Paint("");

	plabel.SetTextSize(0.08);
	plabel.DrawLatexNDC(0.85,0.2, ncls[sqpart.at(i)].sName);
      }
      else
	lg->Draw();

      ipx++;
    }
  }

  if( ccv ) {
    TString ssys = "";
    for( auto iisys : ivsys ) {
      ssys += rsys[iisys]+"Sn";
    }
    ssys += lbCentral[bCentral];
    ccv->SaveAs(Form("%s_%s.png",grname.Data(),ssys.Data()));  
  }
}


void Draw_dndy(UInt_t isys, UInt_t ipart, std::vector<gplot> gname)
{
  Bool_t bAMD    = 1;
  Bool_t bImQMD  = 0;
  Bool_t bpBUU   = 0;
  Bool_t bTommy  = 1;
  Bool_t bKaneko = 1;
  if( !bCentral ) bKaneko = 0;

  auto mgr = new TMultiGraph();
  mgr -> SetTitle(";y_{cm}/y_{prj};dN/dy");

  auto lg  = new TLegend(0.15, 0.2, 0.5, 0.7,"");
  lg -> SetHeader(ncls[ipart].sName);
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
	g_dndy -> SetLineColor( CStyle[i].fColor );
	g_dndy -> SetMarkerColor( CStyle[i].fColor );
	g_dndy -> SetMarkerSize(1.);
	
	if( iisys == 2 ) {
	  g_dndy -> SetMarkerStyle(24);
	  g_dndy -> SetLineColor( CStyle[0].fColor );
	  g_dndy -> SetMarkerColor( CStyle[0].fColor );
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

      g_dndy = (TGraphErrors*)LoadTommyData(iisys, ipart, hname);
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
    g_dndy = (TGraphErrors*)LoadKanekoData(isys, ipart);
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
  //  plabel.DrawLatexNDC(0.15,0.9, fsys[isys]+"_"+ncls[ipart].sName);
  plabel.DrawLatexNDC(0.15,0.9, fsys[isys]+" "+ncls[ipart].sName); 
}

void Draw_dndpt(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, std::vector<gplot> gdata) 
{
  LOG(INFO) << "[Draw_dndpt] .......... "<< ivsys.at(0)  << FairLogger::endl;  
  Bool_t bAMD  = 1;
  if( ivsys.size() > 1 ) bAMD = 0;

  THStack*  mhst[5];
  for( auto ipart : sqpart ) {
    mhst[ipart] = new THStack(Form("dndpt_%s",ncls[ipart].sName.Data()),";Pt;1/NdN/dpt");
    //    mhst[ipart] -> SetDirectory(0);
  }

  auto lg  = new TLegend(0.15, 0.2, 0.5, 0.7,"");

  TH1D* hst = NULL;

  for( auto iisys : ivsys){
    
    for( auto ipart : sqpart) {

	UInt_t ig = 0;
	for( auto igname : gdata ) {

	  hst = (TH1D*)LoadData(iisys,ipart,"h_dndpt",igname,3); 

	  if( hst != NULL ) {
	    hst -> SetMarkerStyle(DStyle[iisys].mStyle);
	    hst -> SetMarkerSize( DStyle[iisys].mSize-0.6);
	    //	    hst -> SetFillStyle( 4001);
	    hst -> SetFillColor(  DStyle[iisys].fColor);
	    hst -> SetLineColor(  DStyle[iisys].fColor);
	    hst -> SetMarkerColor(DStyle[iisys].fColor);

	    mhst[ipart] -> Add( hst );
	    if( ipart == 0 )
	      lg  -> AddEntry( hst, "DATA"+fsys[iisys]+" "+igname.config1);
	  }
	}
    }
  }
  const Int_t Ny = 1;
  const Int_t Nx = 6;
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1800,500); iccv++;
  ccv->Divide(Nx,Ny);

  ccv->cd(1);
  lg->Draw();

  for( auto ipart : sqpart ){
    auto pad = ccv->cd(ipart+2);
    pad->SetLogy();
    pad->SetGrid();
    pad->SetLeftMargin(0.12);
    pad->SetRightMargin(0.02);
    //    mhst[ipart]->GetXaxis()->SetLimits(0.,1200.);
    //    mhst[ipart]->GetXaxis()->SetMaxDigits(2); 
    mhst[ipart]->Draw("nostack");
  }
  
  TString ssys = "";
  for( auto iisys : ivsys ) {
    ssys += rsys[iisys]+"Sn";
  }
  ssys += lbCentral[bCentral];
  ccv->SaveAs(Form("dndpt_%s.png",ssys.Data()));  

}

void Draw_dnIntegral(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, std::vector<gplot> gname, Int_t ylim=-1.)
{
  LOG(INFO) << "[Draw_dnIntegral] ......" << endl;

  const UInt_t npart = 1;
  auto iisys  = ivsys.at(0);
  auto ipart  = 6;
  auto igname = gname.at(0);

  THStack*  mhsty[npart][3];
  mhsty[0][0] = new THStack();
  mhsty[0][1] = new THStack();
  mhsty[0][2] = new THStack();
  TH1D*       h_dndy;

  Double_t var[8];

  h_dndy = (TH1D*)LoadData(iisys,ipart,"h_dndy",igname); 
  if( h_dndy != NULL ) {

    h_dndy -> SetMarkerStyle(DStyle[2].mStyle);
    h_dndy -> SetMarkerSize( DStyle[2].mSize-0.6);
    h_dndy -> SetLineColor(  DStyle[2].fColor);
    h_dndy -> SetMarkerColor(DStyle[2].fColor);
    
    mhsty[0][0] -> Add( h_dndy );

    auto x1 = h_dndy -> FindBin(0.);
    auto x2 = h_dndy -> FindBin(1.);
    var[0]  = 2.* h_dndy -> IntegralAndError(x1, x2, var[1],"width");

    x1 = h_dndy -> FindBin(-0.3);
    x2 = h_dndy -> FindBin( 0.3);
    var[6]  = h_dndy -> IntegralAndError(x1, x2, var[7],"width");
    //    var[6] /= 0.6;
  }

  h_dndy = (TH1D*)LoadData(iisys,ipart,"h_dndyx",igname,3); 
  if( h_dndy != NULL ) {

    h_dndy -> SetMarkerStyle(DStyle[0].mStyle);
    h_dndy -> SetMarkerSize( DStyle[0].mSize-0.6);
    h_dndy -> SetLineColor(  DStyle[0].fColor);
    h_dndy -> SetMarkerColor(DStyle[0].fColor);
    
    mhsty[0][1] -> Add( h_dndy );

    auto x1 = h_dndy -> FindBin(-1.);
    auto x2 = h_dndy -> FindBin( 1.);
    var[2]  = h_dndy -> IntegralAndError(x1, x2, var[3],"width");
  }

  h_dndy = (TH1D*)LoadData(iisys,ipart,"h_dndpx",igname,3); 
  if( h_dndy != NULL ) {

    h_dndy -> SetMarkerStyle(DStyle[1].mStyle);
    h_dndy -> SetMarkerSize( DStyle[1].mSize-0.6);
    h_dndy -> SetLineColor(  DStyle[1].fColor);
    h_dndy -> SetMarkerColor(DStyle[1].fColor);
    
    mhsty[0][2] -> Add( h_dndy );

    auto x1 = h_dndy -> FindBin(-1.);
    auto x2 = h_dndy -> FindBin( 1.);
    var[4]  = h_dndy -> IntegralAndError(x1, x2, var[5],"width");
  }

  for( auto i : {0, 2, 4, 6} ) {
    cout << " Int[" << i << "] " << var[i] << " +- " << var[i+1] << endl;
  }
  

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  ccv -> Divide(2,2);
  ccv -> cd(1);
  mhsty[0][0] -> Draw();
  ccv -> cd(3);
  mhsty[0][1] -> Draw();
  ccv -> cd(4);
  mhsty[0][2] -> Draw();

}

void Get_ClusterMultiplicity(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, mplot samd, Double_t numClust[9][2])
{
  LOG(INFO) << "[Draw_ClusterMultiplicity] .......... " << FairLogger::endl;
  std::vector<TString> partname  = {"1H", "2H", "3H", "3He", "4He", "H", "HHe", "n"};

  for( auto iisys : ivsys){

    for( auto ic : {0,1} ) { // midcentral andl central
    
      if( bCentral != ic) continue;
      
      for( auto ipart : sqpart ) {

	auto hist = (TH1D*)LoadAMD(iisys, 0, samd, "hmult_"+partname[ipart], ic);
	//	auto hb   = (TH1D*)LoadAMD(iisys, 0, samd, "hb",   ic);

	for( auto i : ROOT::TSeqI(9) ) {
	  numClust[i][0] = 0.;
	  numClust[i][1] = 0.;
	}
	if( hist ) {
	  numClust[ipart][0] = hist->GetMean();
	  numClust[ipart][1] = hist->GetMeanError();  ///nevt;
	}
      }

      for( auto i : ROOT::TSeqI(5) ) {
	if( i < 3 ) {
	  numClust[5][0] += numClust[i][0];
	  numClust[5][1] += pow(numClust[i][1],2);
	}
	numClust[6][0] += numClust[i][0] *  ncls[i].Z;
	numClust[6][1] += pow(numClust[i][1], 2) * ncls[i].Z;
	numClust[8][0] += numClust[i][0] * (ncls[i].A - ncls[i].Z);
	numClust[8][1] += pow(numClust[i][1], 2) * (ncls[i].A - ncls[i].Z);
      }
      numClust[5][1] = sqrt(numClust[5][1]);
      numClust[6][1] = sqrt(numClust[6][1]);
      numClust[8][1] = sqrt(numClust[8][1]);

    }
    
  }
  return &numClust[0][0];
}

void Draw_cluster(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, mplot samd, Double_t numClust[9][2])
{
  LOG(INFO) << "[Draw_cluster] .......... " << FairLogger::endl;

  for( auto iisys : ivsys){

    for( auto ic : {0,1} ) { // midcentral andl central
    
      if( bCentral != ic) continue;
      
      auto hist = (TH1D*)LoadAMD(iisys, 0, samd, "hnuclei", ic);
      auto hb   = (TH1D*)LoadAMD(iisys, 0, samd, "hb",   ic);

      for( auto i : ROOT::TSeqI(9) ) {
	numClust[i][0] = 0.;
	numClust[i][1] = 0.;
      }
      if( hist && hb ) {
	Double_t nevt = hb->GetEntries();
	//proton	
	numClust[0][0] = hist->GetBinContent(1,2); //nevt;
	numClust[0][1] = hist->GetBinError(1,2);  ///nevt;
	numClust[0][1] = numClust[0][1]/nevt * sqrt(1/numClust[0][0] + 1/nevt); 
	numClust[0][0] /= nevt;

	//deuteron
	numClust[1][0] = hist->GetBinContent(2,2);
	numClust[1][1] = hist->GetBinError(2,2);
	numClust[1][1] = numClust[1][1]/nevt* sqrt(1/numClust[1][0] + 1/nevt);
	numClust[1][0] /= nevt;

	//triton
	numClust[2][0] = hist->GetBinContent(3,2);
	numClust[2][1] = hist->GetBinError(3,2);
	numClust[2][1] = numClust[2][1]/nevt* sqrt(1/numClust[2][0] + 1/nevt);
	numClust[2][0] /= nevt;

	//3He
	numClust[3][0] = hist->GetBinContent(2,3);
	numClust[3][1] = hist->GetBinError(2,3);
	numClust[3][1] = numClust[3][1]/nevt* sqrt(1/numClust[3][0] + 1/nevt);
	numClust[3][0] /= nevt;

	//4He
	numClust[4][0] = hist->GetBinContent(3,3);
	numClust[4][1] = hist->GetBinError(3,3);
	numClust[4][1] = numClust[4][1]/nevt* sqrt(1/numClust[4][0] + 1/nevt);
	numClust[4][0] /= nevt;

	//n
	numClust[7][0] = hist->GetBinContent(2,1);
	numClust[7][1] = hist->GetBinError(2,1);
	numClust[7][1] = numClust[7][1]/nevt* sqrt(1/numClust[7][0] + 1/nevt);
	numClust[7][0] /= nevt;

	for( auto i : ROOT::TSeqI(5) ) {
	  if( i < 3 ) {
	    numClust[5][0] += numClust[i][0];
	    numClust[5][1] += pow(numClust[i][1],2);
	  }
	  numClust[6][0] += numClust[i][0] *  ncls[i].Z;
	  numClust[6][1] += pow(numClust[i][1], 2) * ncls[i].Z;
	  numClust[8][0] += numClust[i][0] * (ncls[i].A - ncls[i].Z);
	  numClust[8][1] += pow(numClust[i][1], 2) * (ncls[i].A - ncls[i].Z);
	}
	numClust[5][1] = sqrt(numClust[5][1]);
	numClust[6][1] = sqrt(numClust[6][1]);
	numClust[8][1] = sqrt(numClust[8][1]);

      }
    }
  }
  return &numClust[0][0];
}

void Draw_dndbtgm(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, std::vector<gplot> gname, Int_t ylim=-1.,
		 Bool_t bsaveData=0) 
{
  LOG(INFO) << "[Draw_dndbtgm] .......... "<< ivsys.at(0) << endl;
  Bool_t bAMD  = 1;
  if( ivsys.size() > 1 && bsaveData == 0) bAMD = 0;
  LOG(INFO) << " bAMD : " << bAMD << FairLogger::endl;

  Bool_t bKaneko = 0;

  UInt_t hplt[]    = {0,1,2};
  std::vector< std::pair< UInt_t, TString >> phist = {{1,"h_dndEt"},{2,"h_dndbtgm"},{4,"h_dndpx"}};
  TString Ylabel, Xlabel;

  const UInt_t npart = (UInt_t)sqpart.size();
  const UInt_t ndata = (UInt_t)gname.size();
  const UInt_t nhist = (UInt_t)phist.size();
  THStack*  mhsty[8][nhist];

  for( auto ipart : sqpart ) for( auto hist : phist ) {
      mhsty[ipart][hist.first] = new THStack(Form(hist.second+"_%d",ipart),"");
    }

  auto lg  = new TLegend(0.15, 0.2, 0.5, 0.7,"");

  TH1D* h_dndy = NULL;
  TH1D* h_dndyH[ndata];
  TH1D* h_dndyHHe[ndata];
  // TH1D* hratiot3He = new TH1D("hr_t3He");
  // Th1D* hrationp   = new TH1D("hr_np"};

  Double_t ymax[3][2] = {{35.,100},{0.3,0.8},{100.,210.}};

  Double_t var[16];

  for( auto iisys : ivsys){

    Bool_t bfirst = kTRUE;
    for( auto ihist : phist )for( auto ipart : sqpart ) {
	if( !bsaveData && ihist.first > 2 ) break;

	UInt_t ig = 0;
	for( auto igname : gname ) {
	  if( !bsaveData && igname.centrality != bCentral && bCentral != 2 ) continue;

	  h_dndy = (TH1D*)LoadData(iisys,ipart,ihist.second,igname,ylim); 
	  //	  h_dndy = NULL;

	  //	  cout << igname.config1 << " " << ipart << " " << ihist.second << " ->  ";

	  if( h_dndy != NULL ) {
	    TString nname = (TString)h_dndy->GetName()+Form("%d",ig);
	    h_dndy -> SetName(nname);

	    if( ig == 0 ) {
	      Ylabel = h_dndy->GetYaxis()->GetTitle();
	      Xlabel = h_dndy->GetXaxis()->GetTitle();
	    }

	    h_dndy -> Scale( ncls[ipart].A );

	    h_dndy -> SetMarkerStyle(DStyle[iisys].mStyle);
	    h_dndy -> SetMarkerSize( DStyle[iisys].mSize-0.6);
	    h_dndy -> SetLineColor(  DStyle[iisys].fColor+ig*5);
	    h_dndy -> SetMarkerColor(DStyle[iisys].fColor+ig*5);

	    if( bfirst ){ 
	      lg  -> AddEntry( h_dndy, "DATA"+fsys[iisys]+" "+igname.config1);
	      bfirst = kFALSE;
	    }

	    if( ipart == 0 ) {
	      //	      ymax[ihist.first][0] = h_dndy->GetYaxis()->GetXmax();
	      mhsty[ipart][ihist.first] -> Add( h_dndy );

	      h_dndyH[ig]   = (TH1D*)h_dndy->Clone();
	      h_dndyHHe[ig] = (TH1D*)h_dndy->Clone();
 	    }
	    else if( ipart < 3 ) {
	      mhsty[ipart][ihist.first] -> Add( h_dndy );

	      h_dndyH[ig]  ->Add(h_dndy,1);
	      h_dndyHHe[ig]->Add(h_dndy,1);
	    }
	    else if( ipart < 5 ) {
	      mhsty[ipart][ihist.first] -> Add( h_dndy );

	      h_dndyHHe[ig]->Add(h_dndy,2);
	    }
	    else if( ipart == 5 ) {
	      h_dndyH[ig] -> SetMarkerColor(DStyle[iisys].fColor+ig*5);
	      mhsty[ipart][ihist.first] -> Add( h_dndyH[ig] );
	      //	      if( hist.first == 0 )   lg  -> AddEntry( h_dndyH[ig], "DATA_H"); 
	    }
	    else if( ipart == 6 ) {
	      cout << ipart << " " << h_dndy->GetName() << " -> ";
	      //	      ymax[ihist.first][1] = h_dndyHHe[ig]->GetYaxis()->GetXmax();

	      h_dndyHHe[ig] -> SetMarkerColor(DStyle[iisys].fColor+ig*5);
	      mhsty[ipart][ihist.first] -> Add( h_dndyHHe[ig] );
	      //	      if( hist.first == 0 )   lg  -> AddEntry( h_dndyH[ig], "DATA_HHe"); 
	    }
	    cout << ig << endl;


	    if( bsaveData ){
	      if( ipart < 5 )
		GetIntegralandSigma(h_dndy, var);
	      else if( ipart == 5 ) 
		GetIntegralandSigma(h_dndyH[ig], var);
	      else if( ipart == 6 ) 
		GetIntegralandSigma(h_dndyHHe[ig], var);

	      PhysParameter *sData = &physData[iisys][ipart][ig];

	      if( ihist.second == "h_dndy" ) {
 		sData->version = igname.Version;
		sData->system = iisys;
		sData->centrality=igname.centrality;
		sData->partID=ipart;

		sData->integralM = var[0];
		sData->integralMError = var[1];
		sData->RapiditystdDev = var[2];
		sData->RapiditystdDevError = var[3];
	      }
	      else if( ihist.second == "h_dndyx" ){
		sData->integralyx = var[4];
		sData->integralyxError = var[5];
		sData->xRapiditystdDev = var[6];
		sData->xRapiditystdDevError = var[7];
		cout << ncls[ipart].sName << " yx_std = " << var[6] << " +- " << var[7] << endl;
	      }
	      else if( ihist.second == "h_dndpx" ){
		sData->pxstdDev = var[8];
		sData->pxstdDevError = var[9];
		cout << ncls[ipart].sName << " Px_std = " << var[8] << " +- " << var[9] << endl;
	      }
	      else if( ihist.second == "h_dndEt"){
		sData->meanEt = var[10];
		sData->meanEtError = var[11];
		cout << ncls[ipart].sName << " <Et> = " << var[10] << " +- " << var[11] << endl;
	      }
	      else if( ihist.second == "h_dndbtgm"){
		sData->meanbtgm = var[12];
		sData->meanbtgmError = var[13];
		sData->integralbtgm = var[14];
		sData->integralbtgmError = var[15];
		cout << ncls[ipart].sName << " integral<btgm> = " << var[14] << " +- " << var[15] << endl;
	      }
	    }
	  }	 
	  ig++;
	}

	
    
	if( bKaneko && ihist.second == "h_dndy" ) {
	  h_dndy = (TH1D*)LoadKanekoData(iisys, ipart, "h_dndy");
	  if( h_dndy != NULL ) {
	    h_dndy -> SetFillStyle( 3001 );
	    h_dndy -> SetLineColor( CStyle[8].fColor );
	    h_dndy -> SetMarkerColor( CStyle[8].fColor );
	    mhsty[ipart][0] -> Add( h_dndy, "E4" );
	    if( ipart == 0 )
	      lg  -> AddEntry( h_dndy, CStyle[1].comment );
	  }
	}
  
	//@@@@ AMD
	if( bAMD ) {    

	  for( auto ic : {0,1} ) { // midcentral andl central

	    if( !bsaveData and bCentral != ic) continue;

	    UInt_t iin = 0;
	    for( auto samd : AMDnames ) {

	      if( !bsaveData and samd.category != 1 ) continue;
	  
	      h_dndy = (TH1D*)LoadAMD(iisys, ipart, samd, ihist.second, ic, ylim);
	      LOG(INFO) <<" [AMD] " << samd.config << " " << ihist.second << " iisys " << iisys <<  FairLogger::endl;
	      
	      if( h_dndy != NULL ) {
		h_dndy -> SetName((TString)h_dndy->GetName()+Form("_%d",iin));
		LOG(INFO) <<" [AMD] open " <<  h_dndy->GetName() <<  FairLogger::endl;

		h_dndy -> Scale( ncls[ipart].A );
		
		h_dndy -> SetFillStyle( samd.fStyle );
		h_dndy -> SetLineColor( samd.fColor );
		h_dndy -> SetFillColor( samd.fColor );
		

		mhsty[ipart][ihist.first] -> Add( h_dndy,"aE4");
		if( ipart == 0 && ihist.first == 0 ) 
		  lg  -> AddEntry( h_dndy, "AMD"+ samd.config);

		if( bsaveData ){

		  Double_t var[16];
		  GetIntegralandSigma(h_dndy, var);
	      
		  //@@@physamd
		  PhysParameter *sAMD = &physAMD[iisys][ipart][iin][ic];

		  if( ihist.second == "h_dndy" ) {
		    sAMD->version = samd.config;
		    sAMD->system  = iisys;
		    sAMD->centrality=ic;
		    sAMD->partID    =ipart;
		    
		    sAMD->integralM = var[0];
		    sAMD->integralMError = var[1];
		    sAMD->RapiditystdDev = var[2];
		    sAMD->RapiditystdDevError = var[3];
		  }
		  else if( ihist.second == "h_dndyx" ){
		    sAMD->integralyx      = var[4];
		    sAMD->integralyxError = var[5];
		    sAMD->xRapiditystdDev      = var[6];
		    sAMD->xRapiditystdDevError = var[7];
		  }
		  else if( ihist.second == "h_dndpx" ){
		    sAMD->pxstdDev      = var[8];
		    sAMD->pxstdDevError = var[9];
		  }
		  else if( ihist.second == "h_dndEt"){
		    sAMD->meanEt      = var[10];
		    sAMD->meanEtError = var[11];
		  }
		  else if( ihist.second == "h_dndbtgm"){
		    sAMD->meanbtgm      = var[12];
		    sAMD->meanbtgmError = var[13];
		    sAMD->integralbtgm  = var[14];
		    sAMD->integralbtgmError  = var[15];
		  }
		}
	      }
	      iin++;
	    }
	  }
	}
      }
  }

  //--->
  // if( 0 ) {
  //   const Int_t Ny = 2;
  //   const Int_t Nx = 3;
  //   ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),300*Nx,400*Ny); iccv++;
 
  //   ccv-> Divide(Nx, Ny);
  //   id = 1;

  //   for( auto j :  {1,2} ) {
  //     auto mh = new THStack();
  //     auto lg1 = new TLegend();

  //     auto pad = ccv->cd(id); id++;
  //     pad->SetLogy();

  //     for( auto i : {0,1,2,3,4} ) {
  // 	h_ratio[i][j]->SetMarkerStyle( PStyle[i].mStyle );
  // 	h_ratio[i][j]->SetMarkerColor( PStyle[i].fColor );
  // 	mh->Add(h_ratio[i][j],"p");
  // 	lg1->AddEntry(h_ratio[i][j], PStyle[i].comment );
  //     }
  //     mh->Draw("nostack");
  //     lg1->Draw();
    
  //     ccv->cd(id); id++;
  //     auto h1 =  (TH1D*)h_ratio[2][j]->Clone();
  //     h1 -> Divide( h_ratio[3][j]);
  //     h1 -> Draw();
    
  //     ccv->cd(id); id++;
  //     auto h2 =  (TH1D*)h_ratio[1][j]->Clone();
  //     h2 -> Divide( h_ratio[0][j]);
  //     h2 ->Draw();
  //   }
  // }
  

  if( !bsaveData) {
    //--> plot all histogram
    const Int_t Ny = 3;
    const Int_t Nx = npart+1;
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),300*Nx,800); iccv++;
    ccv->Divide(Nx,Ny);

    Float_t lMargin = 0.15;
    Float_t rMargin = 0.02;
    Float_t bMargin = 0.10;
    Float_t tMargin = 0.02;

    //  mhsty[ipart][2] -> GetXaxis()->SetLimits(-800.,800.);
    
    id = 1;

    for( auto ipady : ROOT::TSeqI(Ny) ) {
      std::vector< UInt_t >::iterator ipx = sqpart.begin(); 
      for( auto ipadx : ROOT::TSeqI(Nx) ) {
	auto pad = ccv -> cd(id); id++;
	pad->SetGrid();
	//      pad->SetLogy();
	pad->SetLeftMargin  (lMargin);
	pad->SetRightMargin (rMargin);
	pad->SetTopMargin   (tMargin);
	pad->SetBottomMargin(bMargin);


	mhsty[*ipx][hplt[ipady]] -> SetTitle();

	if( ipadx < npart ) {
	  mhsty[*ipx][hplt[ipady]] -> Draw("nostack");      

	  if( ipadx < 5 ) {
	    mhsty[*ipx][hplt[ipady]] -> SetMaximum( ymax[ipady][0] );	    
	  }
	  else {
	    mhsty[*ipx][hplt[ipady]] -> SetMaximum( ymax[ipady][1] );	    
	  }

	  mhsty[*ipx][hplt[ipady]] -> GetYaxis()->SetTitle( Ylabel );

	  if( ipady == 0 ) {
	    plabel.SetTextAlign(15);
	    plabel.DrawLatexNDC(0.2,0.9, ncls[ipadx].sName);
	  }
	  // else
	  //   plabel.DrawLatexNDC(0.2,0.9,Form("|y|<%3.1f",ycut[ylim]));
	}
	else if( ipadx == npart && ipady == 1)
	  lg->Draw();
	  
	ipx++;
      }
    }

    TString ssys = "";
    for( auto iisys : ivsys ) {
      ssys += rsys[iisys]+"Sn";
    }
    ssys += lbCentral[bCentral];
    ccv->SaveAs(Form("dndydyx_%s.png",ssys.Data()));  

    if( 0 ) {
      ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1200,800); iccv++;
      ccv -> Divide(3,2);
      id = 1;
      for( auto i : {5,6} )for( auto j: {0,1,2} ) {
	  ccv -> cd(id); id++;
	  mhsty[i][j] -> Draw("nostack");
	}
    }
  }
}

void Draw_dndptRatio(std::vector<UInt_t> ivsys, std::vector<gplot> gname )
{
  LOG(INFO) << "[Draw_dndptRatio] .......... Centrality " << bCentral << " sysID " << ivsys.at(0) ;

  
  TH1D* hptRatioData_t3He;
  TH1D* hptRatioAMD_t3He;


  for( auto iisys : ivsys ) {

    for( auto samd : AMDnames ){
      if(samd.category != 1 ) continue;


      auto h_dndpt_t   = (TH1D*)LoadAMD(iisys, 2, samd, "h_dndpt", 1, 3);
      auto h_dndpt_3He = (TH1D*)LoadAMD(iisys, 3, samd, "h_dndpt", 1, 3);
      if( h_dndpt_t && h_dndpt_3He ){

	hptRatioAMD_t3He = (TH1D*)h_dndpt_t->Clone();
	
	//hptRatioAMD_t3He -> Divide( h_dndpt_3He );

	ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
	ccv -> Divide(1, 2);
	ccv -> cd(1);
	h_dndpt_t -> SetLineColor(2);
	h_dndpt_t -> Draw();
	h_dndpt_3He -> SetLineColor(4);
	//	h_dndpt_3He -> Draw("same");
      
	ccv -> cd(2);
	hptRatioAMD_t3He -> Draw();

      }

    }
  }
}

void Draw_dndydX(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, std::vector<gplot> gname,  Int_t ylim=-1.,
		 Bool_t bsaveData=0, std::vector<UInt_t> hpltorder={0,3}, Int_t icateg=1) 
{
  LOG(INFO) << "[Draw_dndydX] .......... Centrality " << bCentral << " sysID " << ivsys.at(0) <<  " bsaveData " << bsaveData ;
  Bool_t bAMD  = 1;
  if( ivsys.size() > 2 && bsaveData == 0) bAMD = 0;
  LOG(INFO) << " bAMD : " << bAMD << FairLogger::endl;


  Bool_t bKaneko = 0;

  std::vector< std::pair< UInt_t, TString >> phist = {{0, "h_dndy"},{1,"h_dndEt"},{2,"h_dndbtgm"},{3,"h_dndyx"},{4,"h_dndpx"},{5,"h_dndpt"}};
  TString ylabel[] = {"A*dN/dy","A*dN/Et","A*dN/d(#beta#gamma_{T})","A*dN/dy_{x}","dN/dP_{x}","dN/dP_{t}"};
  TString xlabel[] = {"y_{cm}/y_{prj}","E_{t}","#beta#gamma_{T}","y_{x}","P_{x}","P_{t}"};

  Double_t umax[][2] = {{48.,-1.95},{0.4,0},{140.,1.},{35.0,-1.95},{0.,0.},{0.,0.}};
  Double_t xlim[][2] = {{-1.8,1.8},{0.,500.},{0.,1.2},{-1.8,1.8},{0.,1000.},{0.,1000.}};

  std::vector< UInt_t > hbgpt = {2,5};
  std::vector< UInt_t > hdsp  = {0,1,2,5};
  std::vector< UInt_t > hsave = {0,1,2,3,4};
  std::vector< UInt_t > hyyx  = {0,3};
  std::vector< UInt_t > hyyxbg= {0,3,2};


  //  std::vector< UInt_t > hplt =  hyyx;//hbgpt; //hyyx;
  //  std::vector< UInt_t > hplt =  {2};
  std::vector< UInt_t > hplt =  hpltorder;
  if( bsaveData )
    hplt = hsave;
  std::vector< UInt_t >::iterator itplt;

  TString Ylabel, Xlabel;
  
  if( !bAMD && sqpart.at(sqpart.size()-1) == 7 )
    sqpart.pop_back();

  const UInt_t npart = (UInt_t)sqpart.size();
  const UInt_t ndata = (UInt_t)gname.size();
  const UInt_t nhist = (UInt_t)hplt.size();
  THStack*  mhsty[8][nhist];

  for( auto ipart : sqpart ) for( itplt = hplt.begin(); itplt != hplt.end(); itplt++ ) {
      Int_t j = itplt - hplt.begin();
      TString histname = phist[j].second;
      cout << j << " p " << ipart << " h " << phist[*itplt].second  << " -> " << histname  <<endl; 
      mhsty[ipart][j] = new THStack(Form(histname+"_%d",ipart),"");
    }

  auto lg  = new TLegend(0.15, 0.2, 0.5, 0.7,fsys[ivsys[0]]);

  TH1D* h_dndy = NULL;
  TH1D* h_dndyH[ndata];
  TH1D* h_dndyHHe[ndata];

  Double_t var[16];

  for( auto ic : {0,1} ) { // midcentral andl central
  //  for( auto ic : {4} ) { // midcentral andl central
    
    if( !bsaveData && bCentral != ic  ) continue;

    for( auto iisys : ivsys){

      for( itplt = hplt.begin(); itplt != hplt.end(); itplt++ ) {
	Int_t j = itplt - hplt.begin();
	TString histname = phist[*itplt].second;
	
	for( auto ipart : sqpart ) {
	  
	  for( auto igname : gname ) {
	    
	    if( igname.centrality != ic ) continue;

	    h_dndy = (TH1D*)LoadData(iisys,ipart,histname,igname,ylim); 
	    //	  h_dndy = NULL;

	    cout << igname.config1 << " " << ipart << " " << histname << " ->  ";

	    if( h_dndy != NULL ) {
	      TString nname = (TString)h_dndy->GetName()+Form("%d",ic);
	      h_dndy -> SetName(nname);

	      Ylabel = h_dndy->GetYaxis()->GetTitle();
	      Xlabel = h_dndy->GetXaxis()->GetTitle();
	      
	      if( !bsaveData  )
		h_dndy -> Scale( ncls[ipart].A );
	      
	      h_dndy -> SetMarkerStyle(DStyle[iisys].mStyle);
	      h_dndy -> SetMarkerSize( DStyle[iisys].mSize-0.6);
	      h_dndy -> SetLineColor(  DStyle[iisys].fColor);//+ic*5);
	      h_dndy -> SetMarkerColor(DStyle[iisys].fColor);//+ic*5);
	      
	      if( j == 0 && ipart == 0 )
		//		lg  -> AddEntry( h_dndy, "DATA");//+fsys[iisys]+" "+igname.config1);
		lg  -> AddEntry( h_dndy, "DATA");//+fsys[iisys]);//+" "+igname.config1);
	      

	      mhsty[ipart][j] -> Add( h_dndy );
	      
	      if( histname == "h_dndyx" && j == 1 && kFALSE){
		h_dndy->SetMarkerStyle(25);
		h_dndy->SetMarkerColor(4);
		mhsty[ipart][0] -> Add( h_dndy );
	      }

	      if( ipart == 0 ) {
		h_dndyH[ic]   = (TH1D*)h_dndy->Clone();
		h_dndyHHe[ic] = (TH1D*)h_dndy->Clone();
	      }
	      else if( ipart < 3 ) {
		h_dndyH[ic]  ->Add(h_dndy,1);
		h_dndyHHe[ic]->Add(h_dndy,1);
	      }
	      else if( ipart < 5 ) {
		h_dndyHHe[ic]->Add(h_dndy,2);
	      }
	      else if( ipart == 5 ) {
		h_dndyH[ic] -> SetMarkerColor(DStyle[iisys].fColor+ic*5);
		mhsty[ipart][j] -> Add( h_dndyH[ic] );
		//	      if( hist.first == 0 )   lg  -> AddEntry( h_dndyH[ic], "DATA_H"); 
	      }
	      else if( ipart == 6 ) {
		cout << ipart << " " << h_dndy->GetName() << " -> ";
		
		h_dndyHHe[ic] -> SetMarkerColor(DStyle[iisys].fColor+ic*5);
		mhsty[ipart][j] -> Add( h_dndyHHe[ic] );
		//	      if( hist.first == 0 )   lg  -> AddEntry( h_dndyH[ic], "DATA_HHe"); 
	      }

	      if( bsaveData ){
		if( ipart < 5 )
		  GetIntegralandSigma(h_dndy, var);
		else if( ipart == 5 ) 
		  GetIntegralandSigma(h_dndyH[ic], var);
		else if( ipart == 6 ) 
		  GetIntegralandSigma(h_dndyHHe[ic], var);
		
		PhysParameter *sData = &physData[iisys][ipart][ic];

		if( histname == "h_dndy" ) {
		  sData->version = igname.Version;
		  sData->system = iisys;
		  sData->centrality=igname.centrality;
		  sData->partID=ipart;

		  sData->integralM = var[0];
		  sData->integralMError = var[1];
		  sData->RapiditystdDev = var[2];
		  sData->RapiditystdDevError = var[3];
		}
		else if( histname == "h_dndyx" ){
		  sData->integralyx = var[4];
		  sData->integralyxError = var[5];
		  sData->xRapiditystdDev = var[6];
		  sData->xRapiditystdDevError = var[7];
		  cout << ncls[ipart].sName << " yx_std = " << var[6] << " +- " << var[7] << endl;
		}
		else if( histname == "h_dndpx" ){
		  sData->pxstdDev = var[8];
		  sData->pxstdDevError = var[9];
		  cout << ncls[ipart].sName << " Px_std = " << var[8] << " +- " << var[9] << endl;
		}
		else if( histname == "h_dndEt"){
		  sData->meanEt = var[10];
		  sData->meanEtError = var[11];
		  cout << ncls[ipart].sName << " <Et> = " << var[10] << " +- " << var[11] << endl;
		}
		else if( histname == "h_dndbtgm"){
		  sData->meanbtgm = var[12];
		  sData->meanbtgmError = var[13];
		  sData->integralbtgm = var[14];
		  sData->integralbtgmError = var[15];
		  cout << ncls[ipart].sName << " integral<btgm> = " << var[14] << " +- " << var[15] << endl;
		}
	      }
	    }	 
	  }
	
    
	  if( bKaneko && histname == "h_dndy" ) {
	    h_dndy = (TH1D*)LoadKanekoData(iisys, ipart, "h_dndy");
	    if( h_dndy != NULL ) {
	      h_dndy -> SetFillStyle( 3001 );
	      h_dndy -> SetLineColor( CStyle[8].fColor );
	      h_dndy -> SetMarkerColor( CStyle[8].fColor );
	      mhsty[ipart][j] -> Add( h_dndy, "E4" );
	      if( ipart == 0 )
		lg  -> AddEntry( h_dndy, CStyle[1].comment );
	    }
	  }
  
	  //@@@@ AMD
	  if( bAMD ) {    

	    for( auto samd : AMDnames ) {

	      if( !bsaveData and samd.category != icateg ) continue;
	  
	      h_dndy = (TH1D*)LoadAMD(iisys, ipart, samd, histname, ic, ylim);
	      LOG(INFO) <<" [AMD] " << samd.config << " " << histname << " iisys " << iisys <<  FairLogger::endl;
	      
	      if( h_dndy != NULL ) {
		
		UInt_t iin = samd.id;

		h_dndy -> SetName((TString)h_dndy->GetName()+Form("_%d",iin));
		LOG(INFO) <<" [AMD] open " <<  h_dndy->GetName() <<  FairLogger::endl;

		if( !bsaveData )
		  h_dndy -> Scale( ncls[ipart].A );
		
		h_dndy -> SetFillStyle( samd.fStyle );
		h_dndy -> SetLineColor( samd.fColor );
		h_dndy -> SetFillColor( samd.fColor );
		// h_dndy -> SetFillStyle( 3001 );
		// h_dndy -> SetLineWidth(-5002 );
		
		mhsty[ipart][j] -> Add( h_dndy,"aE4");
		if( ipart == 0 && j == 0 ) 
		  lg  -> AddEntry( h_dndy, "AMD"+ samd.config);

		if( bsaveData ){

		  Double_t var[16];
		  GetIntegralandSigma(h_dndy, var);
	      
		  //@@@physamd
		  PhysParameter *sAMD = &physAMD[iisys][ipart][iin][ic];

		  if( histname == "h_dndy" ) {
		    sAMD->version = samd.config;
		    sAMD->system  = iisys;
		    sAMD->centrality=ic;
		    sAMD->partID    =ipart;
		    
		    sAMD->integralM = var[0];
		    sAMD->integralMError = var[1];
		    sAMD->RapiditystdDev = var[2];
		    sAMD->RapiditystdDevError = var[3];
		  }
		  else if( histname == "h_dndyx" ){
		    sAMD->integralyx      = var[4];
		    sAMD->integralyxError = var[5];
		    sAMD->xRapiditystdDev      = var[6];
		    sAMD->xRapiditystdDevError = var[7];
		  }
		  else if( histname == "h_dndpx" ){
		    sAMD->pxstdDev      = var[8];
		    sAMD->pxstdDevError = var[9];
		  }
		  else if( histname == "h_dndEt"){
		    sAMD->meanEt      = var[10];
		    sAMD->meanEtError = var[11];
		  }
		  else if( histname == "h_dndbtgm"){
		    sAMD->meanbtgm      = var[12];
		    sAMD->meanbtgmError = var[13];
		    sAMD->integralbtgm  = var[14];
		    sAMD->integralbtgmError  = var[15];
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }


  if( !bsaveData) {
    //--> plot all histogram
    const Int_t Ny = nhist;
    const Int_t Nx = npart+1;
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),250*Nx,300*Ny); iccv++;
    //    ccv->Divide(Nx,Ny);
    Float_t lMargin = 0.04;
    Float_t rMargin = 0.01;
    Float_t bMargin = 0.15;
    Float_t tMargin = 0.02;
    Float_t midMargin = 0.00;

    CanvasPartitionTwoColumn( ccv, Nx, Ny, lMargin, rMargin, bMargin, tMargin, midMargin );

    //  mhsty[ipart][2] -> GetXaxis()->SetLimits(-800.,800.);
    
    id = 1;

    for( auto ipady : ROOT::TSeqI(Ny) ) {
      std::vector< UInt_t >::iterator ipx = sqpart.begin(); 
      for( auto ipadx : ROOT::TSeqI(Nx) ) {
	//	auto pad = ccv -> cd(id); id++;
	auto pad = (TPad*)gROOT->FindObject(Form("pad_%d_%d",ipadx,Ny-ipady-1));
	if( !pad ) return;

	pad->cd();
	pad->SetGrid();
	//      pad->SetLogy();
	// pad->SetLeftMargin  (lMargin);
	// pad->SetRightMargin (rMargin);
	// pad->SetTopMargin   (tMargin);
	// pad->SetBottomMargin(bMargin);
	
	if( phist.at( hplt[ipady] ).second == "h_dndpt") // ||
	  //	    phist.at( hplt[ipady] ).second == "h_dndbtgm" )
	  pad->SetLogy();


	if( ipadx < npart && (ivsys.size()==1 ||  (ivsys.size()>2&&ipadx>0) ) ){
	  mhsty[*ipx][ipady] -> SetTitle();
	  mhsty[*ipx][ipady] -> Draw("nostack");      
	  if( mhsty[*ipx][ipady] -> GetXaxis() ){
	    mhsty[*ipx][ipady] -> GetXaxis()->SetLimits(xlim[hplt[ipady]][0],xlim[hplt[ipady]][1]);
	    mhsty[*ipx][ipady] -> GetXaxis()->SetNdivisions(4,5,0,kTRUE);
	  }
	  if( umax[hplt[ipady]][0] != 0 ) {
	    mhsty[*ipx][ipady] -> SetMaximum( umax[hplt[ipady]][0] );	    

	    if( umax[hplt[ipady]][1] < 0 ) {
	      mhsty[*ipx][ipady] -> GetXaxis() -> SetRangeUser(umax[hplt[ipady]][1], abs(umax[hplt[ipady]][1]));
	    }
	    else
	      mhsty[*ipx][ipady] -> GetXaxis() -> SetRangeUser(0., abs(umax[hplt[ipady]][1]));
	  }

	  if( ipady == 0 ) {
	    plabel.SetTextAlign(15);
	    plabel.SetTextSize(0.08);
	    plabel.DrawLatexNDC(0.8,0.85, ncls[sqpart.at(ipadx)].sName);
	  }
	  // else
	  //   plabel.DrawLatexNDC(0.2,0.9,Form("|y|<%3.1f",ycut[ylim]));

	  if( ipadx == 0 )
	    mhsty[*ipx][ipady] -> GetYaxis()->SetTitle(ylabel[hplt[ipady]]);

	  //	  if( ipady == 1 ) 
	  mhsty[*ipx][ipady] -> GetXaxis()->SetTitle(xlabel[hplt[ipady]]);

	}
	else if( ipadx == npart && ipady == Ny-1)
	  lg->Draw();
	  
	ipx++;
      }
    }

    TString ssys = "";
    for( auto iisys : ivsys ) {
      ssys += rsys[iisys]+"Sn";
    }
    ssys += lbCentral[bCentral];
    ccv->SaveAs(Form("cmp%s_%s.png",phist[hplt[0]].second.Data(),ssys.Data()));  

    if( 0 ) {
      ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1200,800); iccv++;
      ccv -> Divide(3,2);
      id = 1;
      for( auto i : {5,6} )for( auto j: {0,1,2} ) {
	  ccv -> cd(id); id++;
	  mhsty[i][j] -> Draw("nostack");
	}
    }
  }
}


void Draw_ypt(UInt_t isys, std::vector< gplot >gname)
{
  std::vector< UInt_t > sqpart = {0,1,2,3,4};
  TH2D* h2ypt;
  for( auto igname : gname ) {

    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1000,700); iccv++;
    ccv->Divide(3,2);
    
    UInt_t id = 1;
    for( auto ipart : sqpart ) {
      ccv->cd(id); id++;
    
      h2ypt = (TH2D*)LoadData(isys, ipart, "hyptacp", igname);
      if( !h2ypt ) continue;
      h2ypt -> Draw("colz");

    }

    ccv->cd(id);
    plabel.SetTextAlign(13);
    plabel.SetTextSize(0.1);
    plabel.DrawLatexNDC(0.2,0.7, fsys[isys]);
    plabel.DrawLatexNDC(0.2,0.5, igname.fileHeader+igname.Version);
    plabel.DrawLatexNDC(0.2,0.3, lbCentral[igname.centrality]);
    
  }
}

void Draw_dndyRatio(UInt_t isys, UInt_t ipart, std::vector<gplot> gname,
		    std::vector<std::pair<UInt_t, UInt_t>> ipair )
{
  Bool_t bKaneko = 0; //1
  Bool_t bTommy  = 0; //2
  Bool_t bAMD    = 0; //3
  Bool_t bImQMD  = 0; //4
  Bool_t bpBUU   = 0; //5
 
  for( auto spair : ipair ) {
    if( spair.first == 1 || spair.second == 1 )
      bKaneko = 1;
    if( spair.first == 2 || spair.second == 2 )
      bTommy = 1;
    if( spair.first == 3 || spair.second == 3 )
      bAMD = 1;
  }

  TClonesArray* dataArray = new TClonesArray("TGraphErrors",5);
  TGraphErrors* g_dndy[3];
  TGraphErrors* g_data;

  UInt_t i = 0;
  for( auto igname : gname ) {
    
    g_data =(TGraphErrors*)LoadData(isys,ipart,"g_dndy",igname); 
    if( g_data == NULL ) {
      LOG(ERROR) << " No data is found. " << FairLogger::endl;
      return;
    }
    else {
      new( (*dataArray)[i] ) TGraphErrors(*g_data);
      i++;
    }
  }
  
  if( bKaneko ) {
    g_dndy[0] = (TGraphErrors*)LoadKanekoData(isys, ipart);
    if( g_dndy[0] == NULL ) 
      LOG(ERROR) << " Kaneko's data are not found. " << FairLogger::endl;
    else
      g_dndy[0]->SetName("KanekoAna");
  }

  if( bTommy ) {
    g_dndy[1] = (TGraphErrors*)LoadTommyData(isys, ipart, "g_dndy");
    if( g_dndy[1] == NULL ) 
      LOG(ERROR) << " Tommy's data are not found. " << FairLogger::endl;
    else
      g_dndy[1]->SetName("TommyAna");
  }


  TGraphErrors* g_dndyr[5];
  UInt_t ir = 0;
  TIter next( dataArray );
  while( (g_data = (TGraphErrors*)next() ) ) {
    for( auto fpair : ipair ) {
      if( g_dndy[fpair.second-1] ) {
	g_dndyr[ir] = new TGraphErrors();
	
	UInt_t iin = 0;
	for( auto ix : ROOT::TSeqI(g_data->GetN()) ) {
	  Double_t x, y;
	  g_data->GetPoint(ix, x, y);
	  Double_t numerator = g_dndy[fpair.second-1]->Eval(x);
	  if( y > 0 ) {
	    Double_t n_ratio  = numerator / y;
	    g_dndyr[ir] -> SetPoint(iin, x, n_ratio);	  
	    iin++;
	  }
	}

	if( g_dndyr[ir]->GetN() > 0 ) {
	  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),700,1000); iccv++;
	  ccv->Divide(1,2);
	  ccv->cd(2);
	  g_dndyr[ir] -> SetMarkerStyle(CStyle[0].mStyle);
	  g_dndyr[ir] -> SetMarkerSize(1.);
	  g_dndyr[ir] -> SetMarkerColor(CStyle[0].fColor);
	  g_dndyr[ir] -> GetYaxis()->SetRangeUser(0.8, 1.2);
	  g_dndyr[ir] -> Draw("AP");

	  TLatex tlabel;
	  tlabel.DrawLatex(0.7,1.15,fsys[isys]+" "+ncls[ipart].sName);
	  tlabel.DrawLatex(0.7,1.1 ,(TString)g_dndy[fpair.second-1]->GetName()+"/MizukiAna");
	  
	  
	  ccv->cd(1);
	  g_data->SetMarkerStyle(20);
	  g_data->SetMarkerSize(1);
	  g_data->SetMarkerColor(2);
	  g_data->SetLineColor(2);
	  g_data->Draw("AP");

	  g_dndy[fpair.second-1]->SetMarkerStyle(CStyle[8].mStyle); 
	  g_dndy[fpair.second-1]->SetMarkerSize(1.);
	  g_dndy[fpair.second-1]->SetLineColor(CStyle[8].fColor);
	  g_dndy[fpair.second-1]->SetMarkerColor(CStyle[8].fColor);
	  g_dndy[fpair.second-1]->Draw("P");
		  
	  
	}
      }
    }
  }
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
    lg = new TLegend(0.2, 0.2, 0.5, 0.5,ncls[ipart].sName);

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


void Draw_systematicError(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, gplot gname)
{
  LOG(INFO) << "[Draw_systematicError]... " ;
  for( auto iisys : ivsys )
    LOG(INFO) << rsys[iisys] << " : " ;
  LOG(INFO)<< FairLogger::endl; 



  TMultiGraph*  mgrph; 
  TGraphErrors*  grph;

  for( auto iisys : ivsys ) {
    mgrph = new TMultiGraph();
    mgrph -> SetTitle(fsys[iisys]);
			 
    for( auto ipart : sqpart) {
      grph = new TGraphErrors();
      grph -> SetName(Form("sysError_%d_%d",iisys,ipart));
      UInt_t iin = 0;
      for( auto idata : fphysdatafile ) {
	TString fname = Form("physData_cnt%d%s%s.dat",bCentral,gname.Version.Data(),idata.second.Data());
	Load_physDATA(fname, ipart);
	
	grph -> SetPoint     (iin,  idata.first, physData[iisys][ipart][0].meanpxSlope);
	grph -> SetPointError(iin,            0, physData[iisys][ipart][0].meanpxSlopeError);
	iin++;
      }
      grph->SetMarkerStyle(20);
      grph->SetMarkerSize(1.);
      grph->SetMarkerColor( CStyle[ipart].fColor );
      

      mgrph -> Add(grph,"P");
    }

    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++; 
    mgrph -> Draw("AP");
  }

  // auto  mgrpha = new TMultiGraph();
  // for( auto iisys : ivsys )for( auto ipart : sqpart) {
  //     grph = new TGraphErrors();
  //     grph -> SetName(Form("AMDsysError_%d_%d",iisys,ipart));
  //     UInt_t iin = 0;
  //     for( auto idata : fphysdatafile ) {
  // 	if( idata.second == "_fit0to05" ) continue;

  // 	TString fname = Form("physAMD_cnt%d%s.dat",bCentral,idata.second.Data());
  // 	if( Load_physAMD(fname, ipart) ) {
  // 	  grph -> SetPoint     (iin,  iin, physAMD[iisys][ipart][5].meanpxSlope);
  // 	  grph -> SetPointError(iin,    0, physAMD[iisys][ipart][5].meanpxSlopeError);
  // 	  iin++;
  // 	}
  //     }
  //     grph->SetMarkerStyle(20);
  //     grph->SetMarkerSize(1.);
  //     grph->SetMarkerColor( CStyle[ipart].fColor );
      

  //     mgrpha -> Add(grph,"P");
  //   }

  // ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1000,700); iccv++; 
  // mgrpha -> Draw("AP");
}

  //@@@@comp
void Draw_compCorrelation(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, std::vector<gplot> gname, UInt_t isg=2, TString phys="pxstdDev",
			  UInt_t bnrm=1, UInt_t AMDid=2)
{
  Bool_t bAMD = 1;
  if( ivsys.size() > 1 ) bAMD=0;

  LOG(INFO) << "[Draw_compCorrelation]... " ;

  auto f1func = new TF1("f1func","[0]/(x-[1])",0.,500.);

  for( auto iisys : ivsys ) 
    LOG(INFO) << rsys[iisys] << " : " ;
  LOG(INFO)<< FairLogger::endl; 

  UInt_t npart = 8;
 
  TString fname;
  Bool_t fload = kTRUE;

  fname = Form("physData%s.dat",fphysdatafile[fphysdataid].second.Data());
  fload *= Load_physDATA(fname);

  if( !fload ) return;

  auto lg  = new TLegend(0.,0.1,0.8,0.9,Form("Slope fit in %4.2f<y<0.5",fphysdatafile[fphysdataid].first));
  auto lgt = new TLegend(0.6,0.7,0.9,0.95,Form("Slope fit in %4.2f<y<0.5",fphysdatafile[fphysdataid].first));

  TMultiGraph*  mtotal = new TMultiGraph();
  TMultiGraph*  mgrph[npart];
  TGraphErrors* grph[npart][4];
  for(auto i: ROOT::TSeqI(npart) ) {
    mgrph[i] = new TMultiGraph();
    for(auto j: ROOT::TSeqI(4) ) 
      grph[i][j]  = new TGraphErrors();
  }

  auto msys  = new TMultiGraph();

  UInt_t isg0 = isg;
  UInt_t isg1 = isg;
  if( isg == 2 ) {
    isg0 = 1;
    isg1 = 0;
  }


  for( auto ipart : sqpart ) {
    if( ipart == 7 ) continue;
    
    for( auto iisys : ivsys ) {

      if( phys == "xRapiditystdDev" ) {
	grph[ipart][iisys] -> SetPoint     ( 0, physData[iisys][ipart][isg0].xRapiditystdDev,      physData[iisys][ipart][isg1].meanpxSlope );
	grph[ipart][iisys] -> SetPointError( 0, physData[iisys][ipart][isg0].xRapiditystdDevError, physData[iisys][ipart][isg1].meanpxSlopeError );
      }
      else if( phys == "pxstdDev" ) {
	Double_t nrmx  = physData[iisys][ipart][isg0].pxstdDev;
	Double_t nrmy  = physData[iisys][ipart][isg1].meanpxSlope;
	Double_t nrmxe = physData[iisys][ipart][isg0].pxstdDevError;
	Double_t nrmye = physData[iisys][ipart][isg1].meanpxSlopeError;

	if( bnrm == 1 && physData[iisys][6][isg0].pxstdDev != 0 && physData[iisys][6][isg1].meanpxSlope != 0) {
	  nrmx  = physData[iisys][ipart][isg0].pxstdDev / physData[iisys][6][isg0].pxstdDev;
	  nrmxe = GetError( physData[iisys][ipart][isg0].pxstdDev,      physData[iisys][6][isg0].pxstdDev,
			    physData[iisys][ipart][isg0].pxstdDevError, physData[iisys][6][isg0].pxstdDevError);

	  nrmy  = physData[iisys][ipart][isg1].meanpxSlope / physData[iisys][6][isg1].meanpxSlope;
	  nrmye = GetError(physData[iisys][ipart][isg1].meanpxSlope,      physData[iisys][6][isg1].meanpxSlope, 
			   physData[iisys][ipart][isg1].meanpxSlopeError, physData[iisys][6][isg1].meanpxSlopeError); 
	}

	//	cout << "DATA nrmx " << nrmx << " nrmy " << nrmy << endl;
	grph[ipart][iisys] -> SetPoint     ( 0, nrmx,  nrmy );
	grph[ipart][iisys] -> SetPointError( 0, nrmxe, nrmye );

      }
      else if( phys == "xRaid_px" ){
	grph[ipart][iisys] -> SetPoint     ( 0, physData[iisys][ipart][isg1].pxstdDev,      physData[iisys][ipart][isg1].xRapiditystdDev );
	grph[ipart][iisys] -> SetPointError( 0, physData[iisys][ipart][isg1].pxstdDevError, physData[iisys][ipart][isg1].xRapiditystdDevError );
      }
      else if( phys == "meanbtgm" ) {
	grph[ipart][iisys] -> SetPoint     ( 0, physData[iisys][ipart][isg0].meanbtgm,      physData[iisys][ipart][isg1].meanpxSlope );
	grph[ipart][iisys] -> SetPointError( 0, physData[iisys][ipart][isg0].meanbtgmError, physData[iisys][ipart][isg1].meanpxSlopeError );
      }

      if( bAMD ) 
	grph[ipart][iisys] -> SetMarkerStyle( PStyle[ipart+8].mStyle );//45 );
      else
	grph[ipart][iisys] -> SetMarkerStyle( PStyle[ipart].mStyle );//45 );

      grph[ipart][iisys] -> SetMarkerSize(  PStyle[ipart].mSize+0.5 );
      grph[ipart][iisys] -> SetMarkerColor( DStyle[iisys].fColor );
      grph[ipart][iisys] -> SetLineColor(   DStyle[iisys].fColor );

      mtotal -> Add(grph[ipart][iisys],"p");
      mgrph[ipart] -> Add(grph[ipart][iisys],"P");
      if( ipart == 0 ) {
	lg ->AddEntry(grph[ipart][iisys], "DATA "+fsys[iisys]+" "+lbCentral[bCentral]);
	lgt->AddEntry(grph[ipart][iisys], fsys[iisys]);
      }
      if( ipart < 5 && bAMD )
	msys -> Add(grph[ipart][iisys], "");

    }
  }

  //------ AMD ------

  TMultiGraph*  mtotalAMD[nAMD]; 
  TGraphErrors* grphAMD[npart][nAMD];
  TMultiGraph*  msysAMD[nAMD];
  auto lgpart = new TLegend(0.1,0.4,0.7,0.9,"");

  UInt_t sg = 99;

  for(auto j: ROOT::TSeqI(nAMD) ) {
    mtotalAMD[j] = new TMultiGraph();
    msysAMD[j]   = new TMultiGraph();

    for(auto i: ROOT::TSeqI(npart) ) {
      grphAMD[i][j]  = new TGraphErrors();
    }
  }
 

  if( bAMD ) {
    fname = Form("physAMD%s.dat",fphysdatafile[fphysdataid].second.Data());
    fload *= Load_physAMD(fname, sqpart);
    if( !fload ) return;      

    UInt_t isg0 = isg;
    UInt_t isg1 = isg;
    if( isg == 2 ) {
      isg0 = 1;
      isg1 = 0;
    }

    for( auto ipart : sqpart ) {
      UInt_t is[npart][nAMD];
      for( auto i : ROOT::TSeqI(nAMD) )
	is[ipart][i] = 0;
      
      for( auto ig : ROOT::TSeqI( nAMD ) ) {
	for( auto iisys : ivsys ){
	  if( physAMD[iisys][0][ig][isg0].system != 0 && physAMD[iisys][0][ig][isg1].system != 0 ) {
	    //	      AMDnames.at(ig).category >= 1) {
	    if( sg == 99 ) sg = ig;
	    
	    if( ipart == 0 ) {
	      lg->AddEntry(grphAMD[ipart][ig],  Form("%d %s",physAMD[iisys][ipart][ig][isg0].system,AMDnames.at(ig).config.Data()) ); 
	    }

	    if( phys == "xRapiditystdDev" ) {
	      grphAMD[ipart][ig] -> SetPoint     ( is[ipart][ig], physAMD[iisys][ipart][ig][isg0].xRapiditystdDev,      physAMD[iisys][ipart][ig][isg1].meanpxSlope );
	      grphAMD[ipart][ig] -> SetPointError( is[ipart][ig], physAMD[iisys][ipart][ig][isg0].xRapiditystdDevError, physAMD[iisys][ipart][ig][isg1].meanpxSlopeError );
	    }
	    else if( phys == "pxstdDev" ) {
	      Double_t nrmx  = physAMD[iisys][ipart][ig][isg0].pxstdDev;
	      Double_t nrmy  = physAMD[iisys][ipart][ig][isg1].meanpxSlope;
	      Double_t nrmxe = physAMD[iisys][ipart][ig][isg0].pxstdDevError;
	      Double_t nrmye = physAMD[iisys][ipart][ig][isg1].meanpxSlopeError;

	      if( bnrm == 1 && physAMD[iisys][6][ig][isg0].pxstdDev != 0 && physAMD[iisys][6][ig][isg1].meanpxSlope != 0 ) {
		nrmx  = physAMD[iisys][ipart][ig][isg0].pxstdDev / physAMD[iisys][6][ig][isg0].pxstdDev;
		nrmxe = GetError( physAMD[iisys][ipart][ig][isg0].pxstdDev,      physAMD[iisys][6][ig][isg0].pxstdDev,
				  physAMD[iisys][ipart][ig][isg0].pxstdDevError, physAMD[iisys][6][ig][isg0].pxstdDevError);
		
		nrmy  = physAMD[iisys][ipart][ig][isg1].meanpxSlope / physAMD[iisys][6][ig][isg1].meanpxSlope;
		nrmye = GetError(physAMD[iisys][ipart][ig][isg1].meanpxSlope,   physAMD[iisys][6][ig][isg1].meanpxSlope, 
				 physAMD[iisys][ipart][ig][isg1].meanpxSlopeError,   physAMD[iisys][6][ig][isg1].meanpxSlopeError); 
	      }

	      //	      cout << iisys << " : " << ipart << " : " << ig  <<" AMD nrmx " << nrmx << " nrmy " << nrmy << endl;	
	      grphAMD[ipart][ig] -> SetPoint     ( is[ipart][ig], nrmx,    nrmy );
	      grphAMD[ipart][ig] -> SetPointError( is[ipart][ig], nrmxe, nrmye );

	    }
	    else if( phys == "xRaid_px" ){
	      grphAMD[ipart][ig] -> SetPoint     ( is[ipart][ig], physAMD[iisys][ipart][ig][isg0].pxstdDev,      physAMD[iisys][ipart][ig][isg0].xRapiditystdDev );
	      grphAMD[ipart][ig] -> SetPointError( is[ipart][ig], physAMD[iisys][ipart][ig][isg0].pxstdDevError, physAMD[iisys][ipart][ig][isg0].xRapiditystdDevError );
	    }
	    else if( phys == "meanbtgm" ) {
	      grphAMD[ipart][ig] -> SetPoint     ( is[ipart][ig], physAMD[iisys][ipart][ig][isg0].meanbtgm,     physAMD[iisys][ipart][ig][isg1].meanpxSlope );
	      grphAMD[ipart][ig] -> SetPointError( is[ipart][ig], physAMD[iisys][ipart][ig][isg0].meanbtgmError,physAMD[iisys][ipart][ig][isg1].meanpxSlopeError );     
	    }
	  

	    is[ipart][ig]++;


	    grphAMD[ipart][ig] -> SetMarkerStyle( PStyle[ipart].mStyle);//45 );
	    grphAMD[ipart][ig] -> SetMarkerSize(  PStyle[ipart].mSize+0.5 );

	    if( ivsys.size() == 1 ) {
	      grphAMD[ipart][ig] -> SetMarkerColor( AMDnames.at(ig).fColor );
	      grphAMD[ipart][ig] -> SetLineColor(   AMDnames.at(ig).fColor );
	    }
	    else {
	      grphAMD[ipart][ig] -> SetMarkerColor( AMDnames.at(ig).fColor2 );
	      grphAMD[ipart][ig] -> SetLineColor(   AMDnames.at(ig).fColor2 );
	    }

	    mgrph[ipart]  -> Add(grphAMD[ipart][ig],"P");
	    mtotalAMD[ig] -> Add(grphAMD[ipart][ig],"p");

	    if( ipart < 5 )
	      msysAMD[ig] -> Add(grphAMD[ipart][ig],"");
	  }
	}
      }
      lgpart -> AddEntry(grphAMD[ipart][sg],ncls[ipart].Name);
    }
  }

  TString ssys = "";
  for( auto iisys : ivsys ) {
    ssys += rsys[iisys]+"Sn";
  }
  ssys += lbCentral[bCentral];
  
  if( 0 ) {
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1000,700); iccv++; 
    ccv -> Divide(3,2);

    for( auto ipart : sqpart ) {
      ccv->cd(1+ipart);
      if( phys == "xRapiditystdDev" ) {
	mgrph[ipart] -> SetTitle(ncls[ipart].sName+";yx_std; Slope of <px>/A");
	mgrph[ipart] -> GetXaxis()->SetLimits(0.36, 0.58); 
	mgrph[ipart] -> GetXaxis()->SetNdivisions(405); 
	mgrph[ipart] -> GetYaxis()->SetRangeUser(0.,   165.); 
      }
      else if( phys == "pxstdDev" ) {
	mgrph[ipart] -> SetTitle(ncls[ipart].sName+";px_std; Slope of <px>/A");
	mgrph[ipart] -> GetXaxis()->SetLimits(130, 265); 
	//      mgrph[ipart] -> GetXaxis()->SetNdivisions(405); 
	mgrph[ipart] -> GetYaxis()->SetRangeUser(0.,   165.); 
      }
      else if( phys == "xRaid_px" ){
	mgrph[ipart] -> SetTitle(ncls[ipart].sName+";px_std; yx_std");
	mgrph[ipart] -> GetXaxis()->SetLimits(130, 265); 
	mgrph[ipart] -> GetXaxis()->SetNdivisions(405); 
	mgrph[ipart] -> GetYaxis()->SetRangeUser(0.36, 0.58); 
      }
      else if( phys == "meanbtgm" ){
	mgrph[ipart] -> SetTitle(ncls[ipart].sName+";<#beta #gamma>; Slope of <px>/A");
	mgrph[ipart] -> GetXaxis()->SetLimits(130, 265); 
	mgrph[ipart] -> GetXaxis()->SetNdivisions(405); 
	mgrph[ipart] -> GetYaxis()->SetRangeUser(0.36, 0.58); 
      }

      mgrph[ipart] -> Draw("A");
    }

    ccv->cd(6);
    lg->Draw();
  }

  if( !bAMD ) {
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
    mtotal->SetTitle(";StdDev(px/A); Slope of <px>/A");

    f1func->SetLineColor( DStyle[3].fColor );
    f1func->SetLineStyle( 9);

    if( bnrm ) {
      mtotal->Draw("A");
      //      msys[3]->Fit("f1func","","");
      ccv->SetGrid();
    }
    else {
      mtotal->GetXaxis()->SetLimits(125.,270.);    mtotal->GetYaxis()->SetRangeUser(60.,160.);
      mtotal->Draw("A");
      gStyle->SetOptFit(0);
      mtotal->Fit("f1func","","", 130., 260.);  
      //      msys[3]->Fit("f1func","","", 130., 260.);  
      f1func->Draw("same");
      plabel.DrawText(245., 75.,"Proton");
      plabel.DrawText(180., 76.,"Deuteron");
      plabel.DrawText(140.,100.,"Triton");
      plabel.DrawText(170.,130.,"3He");
      plabel.DrawText(145.,150.,"4He");
      plabel.DrawText(220.,100.,"H");
      plabel.DrawText(195.,115.,"H+2He");
      lgt->Draw();
    }
    ccv->SaveAs(Form("cmpCorrelation%d_%s_%s.png",AMDid,phys.Data(),ssys.Data()));
  }


  if( 0 ) {
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv_%d",iccv),1500,500); iccv++;
    ccv->Divide(4,1);

    id = 1;  
    for(auto j: ROOT::TSeqI(nAMD) ) {
      //      if(  AMDnames.at(j).category == 1 && mtotalAMD[j]->GetListOfGraphs()) {
      if( mtotalAMD[j]->GetListOfGraphs()) {
	mtotalAMD[j]->SetTitle(";StdDev(px/A); Slope of <px>/A");
	mtotalAMD[j]->GetXaxis()->SetLimits(130.,270.);
	mtotalAMD[j]->GetYaxis()->SetRangeUser(20.,160.);
	ccv->cd(id); id++;
	mtotalAMD[j]->Draw("A");
	mtotal->Draw();
	plabel.DrawLatexNDC(0.28,0.85,AMDnames.at(j).fname(6,3)+" "+AMDnames.at(j).config);
      }
    }
    ccv->SaveAs(Form("cmpCorreAMD_%s_%s.png",phys.Data(),ssys.Data()));
  }

  if( bAMD ) {
    const UInt_t ngr = 1;
    TMultiGraph* mcmp[ngr];
    TLegend*     lgm[ngr];

    for( auto i : ROOT::TSeqI(ngr) ) {
      mcmp[i] = new TMultiGraph();
      lgm[i]  = new TLegend(0.25,0.75,0.9,0.9,"");//lg->GetHeader());

      lgm[i]->AddEntry(mtotal->GetListOfGraphs()->First(),"DATA "+fsys[ivsys[0]]);
      
      mcmp[i]->Add(mtotal,"p");
    }

    UInt_t k = 0;
    for(auto j: ROOT::TSeqI(nAMD) ) {
      if( mtotalAMD[j]->GetListOfGraphs() ) {

	TString sysname = AMDnames.at(j).fname(6,3);
	//	sysname = "";

	if( j == AMDid ) {
	  lgm[k]->AddEntry(mtotalAMD[j]->GetListOfGraphs()->At(1),AMDnames.at(j).config);
	  mcmp[k]->Add(mtotalAMD[j],"P");
	  msysAMD[j]->Fit("f1func","","", 130., 260.); 
	  f1func->SetLineColor( AMDnames.at(j).fColor );
	  f1func->SetLineStyle( 9);
	  mcmp[k]->Add(msysAMD[j],"l");
	}      
      }
    }

    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv_%d",iccv),400*(ngr+1),500); iccv++;
    ccv->Divide(ngr+1,1);

    Double_t rng[2];
    for(auto k: ROOT::TSeqI(ngr) ) {
      ccv->cd(k+1); 

      if( phys == "xRaid_px" ) {
	mcmp[k]->SetTitle(";StdDev(px/A); Slope of <px>/A");
	if( !bnrm ) {
	  mcmp[k]->GetXaxis()->SetLimits(130.,270.);
	  mcmp[k]->GetYaxis()->SetRangeUser(40.,200.);
	  rng[0] = 130.;
	  rng[1] = 260.;
	}
      }
      else if ( phys == "meanbtgm" ) {
	mcmp[k]->SetTitle(";<#beta #gamma>; Slope of <px>/A");
	if( !bnrm ) {
	  mcmp[k] -> GetXaxis()->SetLimits(0.17, 0.36);
	  mcmp[k] -> GetXaxis()->SetNdivisions(405); 
	  mcmp[k] -> GetYaxis()->SetRangeUser(70., 165.);
	  rng[0] = 0.1;
	  rng[1] = 0.35;
	}
      }

      mcmp[k]->Draw("A");
      msys->Fit("f1func","","", rng[0], rng[1]); 
      f1func->SetLineColor( 1);
      f1func->SetLineStyle( 9);
      f1func->Draw("same");

      lgm[k] ->Draw();
    }
    ccv->cd(ngr+1);
    lgpart->Draw();
    
    ccv->SaveAs(Form("cmpCorreAMD_%s_%s.png",phys.Data(),ssys.Data()));
  }
}


void Draw_ParticleDependence(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, gplot gname, UInt_t isg=0, TString phys="xRapiditystdDev", UInt_t idiv = 0 , UInt_t ndiv = 2, Double_t min=-1, Double_t max=-1, UInt_t icateg=41, UInt_t icateg1=41)
{
  Bool_t bRATIO = 0;

  if( ivsys.size() == 2 )
    bRATIO = 1;

  LOG(INFO) << "[Draw_ParticleDependence]... " << phys << " isg -> "  << isg << " RATIO " << bRATIO << " " ;

  for( auto iisys : ivsys )
    LOG(INFO) << rsys[iisys] << " : " ;
  LOG(INFO)<< FairLogger::endl; 

  auto lgd  = new TLegend(0.2,0.75, 0.9,0.90,"");//lbCentral[isg]);
  auto lg   = new TLegend(0.2,0.65, 0.9,0.78,"");//lbCentral[isg]);
  auto lgr  = new TLegend(0.2,0.2 ,0.4,0.8,"");//lbCentral[isg]);
  lg  -> SetFillStyle(0);
  lgd -> SetFillStyle(0);
  lgr -> SetFillStyle(0);

  //TString fname = Form("physData%s_both.dat",fphysdatafile[fphysdataid].second.Data());
  TString fname = Form("physData%s_left.dat",fphysdatafile[fphysdataid].second.Data());
  //  TString fname = Form("physData%s_3.dat",fphysdatafile[fphysdataid].second.Data());
  Load_physDATA(fname);


  fname = Form("physAMD%s.dat",fphysdatafile[fphysdataid].second.Data());
  Load_physAMD(fname, sqpart);

  const UInt_t ngrph = 3;
  TMultiGraph* mgr[ngrph];
  for( auto i : ROOT::TSeqI(ngrph) ) {
    mgr[i] = new TMultiGraph();
    mgr[i] -> SetName(Form("mgr_%d",i));
  }
  
  TGraphErrors* grph[ngrph];

  auto sDire = new TDirectory("sDire","sDire");
  TGraphErrors* gr132A;
  TGraphErrors* gr108A;
  TGraphErrors* gr132D;
  TGraphErrors* gr108D;
  Int_t  vig132 = -1;
  Int_t  vig108 = -1;
  auto mgrRatio = new TMultiGraph();

  // AMD
  if( ivsys.size() >= 1 || ivsys.size() != 3 ) {
    for( auto ig : ROOT::TSeqI( nAMD ) ) {
      Bool_t bfill = kFALSE;

      for( auto iisys : ivsys ){
	if( physAMD[iisys][0][ig][isg].system != 0 && (AMDnames.at(ig).category == icateg || AMDnames.at(ig).category == icateg1 ) ) {
	  LOG(INFO) << "[AMD compare] " << physAMD[iisys][0][ig][isg].version << " " << physAMD[iisys][0][ig][isg].system << FairLogger::endl;


	  for( auto i : ROOT::TSeqI(3) ) {
	    grph[i] = new TGraphErrors();
	    grph[i] -> SetName(Form("phAMD_%s_%d_%d_%d",phys.Data(),iisys,i,ig) );
	  }
	  
	  if( bRATIO && (iisys == 0 || iisys == 1 ) ) {
	    sDire -> GetList() -> Add( grph[0] );
	  }

	  UInt_t iin = 0;
	  for( auto ipart : sqpart ){
	    if( phys == "xRapiditystdDev" ) {
	      grph[0] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].xRapiditystdDev );
	      grph[0] -> SetPointError(iin,   0, physAMD[iisys][ipart][ig][isg].xRapiditystdDevError );
	    }
	    else if( phys == "RapiditystdDev" ){
	      grph[0] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].RapiditystdDev );
	      grph[0] -> SetPointError(iin,   0, physAMD[iisys][ipart][ig][isg].RapiditystdDevError );
	    }
	    else if( phys == "pxstdDev" ) {
	      grph[0] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].pxstdDev );
	      grph[0] -> SetPointError(iin,   0, physAMD[iisys][ipart][ig][isg].pxstdDevError );

	      grph[1] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].pxstdDev / physAMD[iisys][0][ig][isg].pxstdDev );
	      Double_t error = GetError( physAMD[iisys][ipart][ig][isg].pxstdDev,      physAMD[iisys][0][ig][isg].pxstdDev, 
					 physAMD[iisys][ipart][ig][isg].pxstdDevError, physAMD[iisys][0][ig][isg].pxstdDevError );
	      grph[1] -> SetPointError(iin,   0, error);

	      if( physAMD[1][ipart][ig][isg].pxstdDev != 0 ) {
		bfill = kTRUE;
		grph[2] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].pxstdDev / physAMD[1][ipart][ig][isg].pxstdDev );
		error = GetError( physAMD[iisys][ipart][ig][isg].pxstdDev,      physAMD[1][ipart][ig][isg].pxstdDev, 
				  physAMD[iisys][ipart][ig][isg].pxstdDevError, physAMD[1][ipart][ig][isg].pxstdDevError );
		grph[2] -> SetPointError(iin,  0, error);
	      }
	    }
	    else if( phys == "meanpxSlope" ) {
	      grph[0] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].meanpxSlope );
	      grph[0] -> SetPointError(iin,   0, physAMD[iisys][ipart][ig][isg].meanpxSlopeError );
	      mgr[0]->SetTitle(";;slope<px>");

	      cout << " AMD slope<px> " << iin << " " << physAMD[iisys][ipart][ig][isg].meanpxSlope << endl;
	
	      grph[1] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].meanpxSlope / physAMD[iisys][0][ig][isg].meanpxSlope );
	      Double_t error = GetError( physAMD[iisys][ipart][ig][isg].meanpxSlope,      physAMD[iisys][0][ig][isg].meanpxSlope, 
					 physAMD[iisys][ipart][ig][isg].meanpxSlopeError, physAMD[iisys][0][ig][isg].meanpxSlopeError );
	      grph[1] -> SetPointError(iin,   0, error);
	      mgr[1]->SetTitle(";;R(slope<px>/slope<px>_{proton})");


	      if ( physAMD[1][ipart][ig][isg].meanpxSlope != 0 ) {
		bfill = kTRUE;
		grph[2] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].meanpxSlope / physAMD[1][ipart][ig][isg].meanpxSlope );
		error = GetError( physAMD[iisys][ipart][ig][isg].meanpxSlope,      physAMD[1][ipart][ig][isg].meanpxSlope, 
				  physAMD[iisys][ipart][ig][isg].meanpxSlopeError, physAMD[1][ipart][ig][isg].meanpxSlopeError );
		grph[2] -> SetPointError(iin,   0, error);	
		mgr[2]->SetTitle(";;R(slope<px>/slope<px>_{108Sn+112Sn})");
	      }
	    }
	    else if( phys == "meanbtgm" ) {
	      grph[0] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].meanbtgm );
	      grph[0] -> SetPointError(iin,   0, physAMD[iisys][ipart][ig][isg].meanbtgmError );

	      grph[1] -> SetPoint     (iin, iin+1, 0.5*pow(physAMD[iisys][ipart][ig][isg].meanbtgm,2)*mass[ipart] );
	      grph[1] -> SetPointError(iin,   0, physAMD[iisys][ipart][ig][isg].meanbtgmError );
	    }
	    else if( phys == "meanKt") {
	      grph[0] -> SetPoint     (iin, iin+1, 0.5*pow(physAMD[iisys][ipart][ig][isg].meanbtgm,2)*mass[ipart] );
	      grph[0] -> SetPointError(iin,   0, physAMD[iisys][ipart][ig][isg].meanbtgmError*mass[ipart] );
	      cout << "AMD  meanKt " << iin << " "  << 0.5*pow(physAMD[iisys][ipart][ig][isg].meanbtgm,2)*mass[ipart] << endl;
	    }
	    else if( phys == "v1Slope") {
	      grph[0] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].v1Slope );
	      grph[0] -> SetPointError(iin,     0, physAMD[iisys][ipart][ig][isg].v1SlopeError );
	      mgr[0]->SetTitle(";;v11");
	    }
	    else if( phys == "v2minimum") {
	      grph[0] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].v2minimum );
	      grph[0] -> SetPointError(iin,     0, physAMD[iisys][ipart][ig][isg].v2minimumError );
	      mgr[0]->SetTitle(";;v20");
	    }
	    else if( phys == "v2width") {
	      grph[0] -> SetPoint     (iin, iin+1, physAMD[iisys][ipart][ig][isg].v2width );
	      grph[0] -> SetPointError(iin,     0, physAMD[iisys][ipart][ig][isg].v2widthError );
	      mgr[0]->SetTitle(";;v21");
	    }
	    
	    iin++;
	  }

	  Color_t amdcol = AMDnames.at(ig).fColor;
	  for( auto i : ROOT::TSeqI(ngrph) ) {
	    grph[i]->SetFillStyle( AMDnames.at(ig).fStyle );
	    grph[i]->SetLineColor( amdcol );
	    grph[i]->SetFillColorAlpha( amdcol, 0.01 );

	    mgr[i] ->Add(grph[i], "3");
	  }


	  if( !bRATIO )
	    lg ->AddEntry(grph[0],Form("%s",AMDnames.at(ig).fullconfig.Data())); 

	  if( bfill )
	    lgr ->AddEntry(grph[0],Form("%d %s",physAMD[iisys][0][ig][isg].system,AMDnames.at(ig).fullconfig.Data())); 
	}
      
      
	//	sDire->ls();

	if( bRATIO && iisys == 0 ) {
	  gr132A = (TGraphErrors*)sDire->FindObject(Form("phAMD_%s_0_0_%d",phys.Data(),ig) );
	  if( gr132A )
	    vig132 = ig;
	}
	else if( bRATIO && iisys == 1 ) {
	  gr108A = (TGraphErrors*)sDire->FindObject(Form("phAMD_%s_1_0_%d",phys.Data(),ig) );	
	  if( gr108A )
	    vig108 = ig;
 	}
      }


      if( bRATIO && (vig132 > -1 && vig108 > -1 ) ) {
	gr132A = (TGraphErrors*)sDire->FindObject(Form("phAMD_%s_0_0_%d",phys.Data(),vig132) );
	gr108A = (TGraphErrors*)sDire->FindObject(Form("phAMD_%s_1_0_%d",phys.Data(),vig108) );

	if( !gr132A || !gr108A ) continue;

	auto grRatio = new TGraphErrors();

	Double_t x0,x1,y0,y1,x0e,x1e,y0e,y1e;
	for( auto np : ROOT::TSeqI(gr132A->GetN()) ) {
	  gr132A->GetPoint(np, x0, y0);
	  gr108A->GetPoint(np, x1, y1);
	  y0e = gr132A->GetErrorY(np);
	  y1e = gr108A->GetErrorY(np);

	  if( y0 != 0 ) {
	    grRatio -> SetPoint(np, x0, y1/y0);
	    auto err = GetError(y1, y0, y1e, y0e);
	    grRatio -> SetPointError(np, 0, err);
	  }
	}

	grRatio -> SetFillStyle( gr132A->GetFillStyle() );
	grRatio -> SetFillColor( gr132A->GetFillColor() );
	grRatio -> SetFillColor( gr132A->GetLineColor() );
	
	mgrRatio -> Add(grRatio, "3Y+");
	cout << AMDnames.at(vig132).fullconfig << endl;

	lg ->AddEntry(grRatio,Form("%s",AMDnames.at(vig132).fullconfig.Data())); 

	vig132 = -1;
	vig108 = -1;
      }	
    }
  }


  // data
  for( auto iisys : ivsys ) {

    for( auto i : ROOT::TSeqI(ngrph) ) {
      grph[i] = new TGraphErrors();
      grph[i] -> SetName(Form("phpara_%s_%d_%d",phys.Data(),iisys,i) );
    }

    if( bRATIO && (iisys == 0 || iisys == 1 ) )
      sDire -> GetList()->Add( grph[0] );

    UInt_t iin = 2;
    for( auto ipart : sqpart ) {
      if( ipart == 7 ) continue;

      if( phys == "xRapiditystdDev" ) {
	grph[0] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].xRapiditystdDev );
	grph[0] -> SetPointError(iin,   0, physData[iisys][ipart][isg].xRapiditystdDevError );
	mgr[0]->SetTitle(";;#sigma(Transvese)");

	grph[1] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].xRapiditystdDev / physData[iisys][0][isg].xRapiditystdDev );
	Double_t error = GetError( physData[iisys][ipart][isg].xRapiditystdDev,      physData[iisys][0][isg].xRapiditystdDev, 
				   physData[iisys][ipart][isg].xRapiditystdDevError, physData[iisys][0][isg].xRapiditystdDevError );
	grph[1] -> SetPointError(iin,   0, error);
	mgr[1]->SetTitle(";;R(yx_std/yx_std_{proton})");

	grph[2] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].xRapiditystdDev / physData[1][ipart][isg].xRapiditystdDev );
	error = GetError( physData[iisys][ipart][isg].xRapiditystdDev,      physData[1][ipart][isg].xRapiditystdDev, 
			  physData[iisys][ipart][isg].xRapiditystdDevError, physData[1][ipart][isg].xRapiditystdDevError );
	grph[2] -> SetPointError(iin,   0, error);	
	mgr[2]->SetTitle(";;R(yx_std/yx_std_{108Sn+112Sn})");
      }

      else if( phys == "RapiditystdDev" ){
	grph[0] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].RapiditystdDev );
	grph[0] -> SetPointError(iin,   0, physData[iisys][ipart][isg].RapiditystdDevError );
	mgr[0]->SetTitle(";;#sigma(Longigudinal)");

	grph[1] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].RapiditystdDev / physData[iisys][0][isg].RapiditystdDev );
	Double_t error = GetError( physData[iisys][ipart][isg].RapiditystdDev,      physData[iisys][0][isg].RapiditystdDev, 
				   physData[iisys][ipart][isg].RapiditystdDevError, physData[iisys][0][isg].RapiditystdDevError );
	grph[1] -> SetPointError(iin,   0, error);
	mgr[1]->SetTitle(";;R(y_std/y_std_{proton})");

	grph[2] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].RapiditystdDev / physData[1][ipart][isg].RapiditystdDev );
	error = GetError( physData[iisys][ipart][isg].RapiditystdDev,      physData[1][ipart][isg].RapiditystdDev, 
			  physData[iisys][ipart][isg].RapiditystdDevError, physData[1][ipart][isg].RapiditystdDevError );
	grph[2] -> SetPointError(iin,   0, error);	
	mgr[2]->SetTitle(";;R(y_std/y_std_{108Sn+112Sn})");

      }
      else if( phys == "pxstdDev" ) {
	grph[0] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].pxstdDev );
	grph[0] -> SetPointError(iin,   0, physData[iisys][ipart][isg].pxstdDevError );
	mgr[0]->SetTitle(";;px_std");
	
	grph[1] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].pxstdDev / physData[iisys][0][isg].pxstdDev );
	Double_t error = GetError( physData[iisys][ipart][isg].pxstdDev,      physData[iisys][0][isg].pxstdDev, 
				   physData[iisys][ipart][isg].pxstdDevError, physData[iisys][0][isg].pxstdDevError );
	grph[1] -> SetPointError(iin,   0, error);
	mgr[1]->SetTitle(";;R(px_std/px_std_{proton})");

	if( iisys != 1 ){
	  grph[2] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].pxstdDev / physData[1][ipart][isg].pxstdDev );
	  error = GetError( physData[iisys][ipart][isg].pxstdDev,      physData[1][ipart][isg].pxstdDev, 
			    physData[iisys][ipart][isg].pxstdDevError, physData[1][ipart][isg].pxstdDevError );
	  grph[2] -> SetPointError(iin,   0, error);	
	  mgr[2]->SetTitle(";;R(px_std/px_std_{108Sn+112Sn})");
	}
      }
      else if( phys == "meanpxSlope" ) {
	cout << " meanpxslope " << ipart << " "<< iisys << " " << isg << " " << physData[iisys][ipart][isg].meanpxSlope << endl;

	grph[0] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].meanpxSlope );	grph[0] -> SetPointError(iin,   0, physData[iisys][ipart][isg].meanpxSlopeError );
	mgr[0]->SetTitle(";;slope<px>");
	
	grph[1] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].meanpxSlope / physData[iisys][0][isg].meanpxSlope );
	Double_t error = GetError( physData[iisys][ipart][isg].meanpxSlope,      physData[iisys][0][isg].meanpxSlope, 
				   physData[iisys][ipart][isg].meanpxSlopeError, physData[iisys][0][isg].meanpxSlopeError );
	grph[1] -> SetPointError(iin,   0, error);
	mgr[1]->SetTitle(";;R(slope<px>/slope<px>_{proton})");


	if( iisys != 1 ) {
	  grph[2] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].meanpxSlope / physData[1][ipart][isg].meanpxSlope );
	  error = GetError( physData[iisys][ipart][isg].meanpxSlope,      physData[1][ipart][isg].meanpxSlope, 
			    physData[iisys][ipart][isg].meanpxSlopeError, physData[1][ipart][isg].meanpxSlopeError );
	  grph[2] -> SetPointError(iin,   0, error);	
	  mgr[2]->SetTitle(";;R(slope<px>/slope<px>_{108Sn+112Sn})");
	}
      }
      else if( phys == "meanbtgm" ) {
	grph[0] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].meanbtgm );
	grph[0] -> SetPointError(iin,   0, physData[iisys][ipart][isg].meanbtgmError );
	mgr[0] -> SetTitle(";; <#beta#gamma>");

	grph[1] -> SetPoint     (iin, iin, 0.5*pow(physData[iisys][ipart][isg].meanbtgm,2)*mass[ipart] );
	grph[1] -> SetPointError(iin,   0, physData[iisys][ipart][isg].meanbtgmError );
	mgr[1]  -> SetTitle(";; <K_{T}>");
	mgr[1]  -> GetYaxis()->SetRangeUser(40.,100.);
      }
      else if( phys == "meanKt" ) {
	grph[0] -> SetPoint     (iin, iin, 0.5*pow(physData[iisys][ipart][isg].meanbtgm,2)*mass[ipart] );
	grph[0] -> SetPointError(iin,   0, physData[iisys][ipart][isg].meanbtgmError );
	mgr[0]  -> SetTitle(";; <K_{T}>");
      }
      else if( phys == "v1Slope") {
	grph[0] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].v1Slope );
	grph[0] -> SetPointError(iin,   0, physData[iisys][ipart][isg].v1SlopeError );
	mgr[0]->SetTitle(";;v11");
      }
      else if( phys == "v2minimum") {
	grph[0] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].v2minimum );
	grph[0] -> SetPointError(iin,   0, physData[iisys][ipart][isg].v2minimumError );
	mgr[0]->SetTitle(";;v20");
      }
      else if( phys == "v2width") {
	grph[0] -> SetPoint     (iin, iin, physData[iisys][ipart][isg].v2width );
	grph[0] -> SetPointError(iin,   0, physData[iisys][ipart][isg].v2widthError );
	mgr[0]->SetTitle(";;v21");
      }
      
      iin++;
    }

    for( auto i : ROOT::TSeqI(ngrph) ) {
      grph[i]->SetMarkerStyle( DStyle[iisys].mStyle );
      grph[i]->SetMarkerSize(  DStyle[iisys].mSize);
      grph[i]->SetLineColor(   DStyle[iisys].fColor );
      grph[i]->SetMarkerColor( DStyle[iisys].fColor );
      
      mgr[i] ->Add(grph[i], "p");
    }

    if( !bRATIO )
      lgd  ->AddEntry(grph[0],fsys[iisys]);//+"_"+physData[iisys][0][isg].version); 

    if( iisys != 1 )
      lgr ->AddEntry(grph[0],fsys[iisys]+"/"+fsys[1]+"_"+physData[iisys][0][isg].version); 
  }


  if( bRATIO ) {
    gr132D = (TGraphErrors*)sDire->FindObject(Form("phpara_%s_0_0",phys.Data()) );
    gr108D = (TGraphErrors*)sDire->FindObject(Form("phpara_%s_1_0",phys.Data()) );

    if( gr132D && gr108D ) {
      auto grRatio = new TGraphErrors();


      // gr132D->Print();
      // gr108D->Print();

      Double_t x0,x1,y0,y1,x0e,x1e,y0e,y1e;
      UInt_t iin = 0;
      for( auto np : ROOT::TSeqI(gr132D->GetN()) ) {
	gr132D->GetPoint(np, x0, y0);
	gr108D->GetPoint(np, x1, y1);
	y0e = gr132D->GetErrorY(np);
	y1e = gr108D->GetErrorY(np);
	auto err = GetError(y1, y0, y1e, y0e);
	
	if( y0 != 0 ) {
	  grRatio -> SetPoint(iin, x0, y1/y0);
	  auto err = GetError(y1, y0, y1e, y0e);
	  grRatio -> SetPointError(iin, 0, err);
	  iin++;

	  //	  cout << " data " << np << " " << x0 << " " << y1 << "/" << y0  << "=" << y1/y0 << endl;
	}
      }
      
      // data
      grRatio->Print();
      
      grRatio -> SetMarkerStyle( 25);
      grRatio -> SetMarkerColor( 1 );
      
      mgrRatio -> Add(grRatio, "pY+");
      TString slabel =";;R_"+(TString)mgr[0]->GetYaxis()->GetTitle()+"("+fsys[1]+"/"+fsys[0]+")";
      mgrRatio -> SetTitle( slabel );

      lgd  ->AddEntry(grRatio,"DATA");//fsys[1]+"/"+fsys[0]);//+"_"+physData[iisys][0][isg].version); 
  
    }
  }
  

  SetXTextLabel(*mgr[0], sqpart);
  if( bRATIO )
    SetXTextLabel(*mgrRatio, sqpart);
  if( max != -1 && min != -1) {
    mgr[0]->GetYaxis()->SetRangeUser(min,max);
    mgr[0]->GetYaxis()->SetNdivisions(5,5,0,kTRUE);

    if( bRATIO )
      mgrRatio->GetYaxis()->SetRangeUser(min, max);
  }


  //----------------- plot --------------------
  //-- draft 

  std::vector< std::pair<Int_t,Int_t>> pmap = {{0,1},{1,1},{0,0},{1,0}};
  if( idiv == 0) {
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),800,1000); iccv++;
    
    const Int_t Ny = 2;
    const Int_t Nx = 2;
    Float_t lMargin = 0.13;
    Float_t rMargin = 0.02;
    Float_t bMargin = 0.1;
    Float_t tMargin = 0.02;
    Float_t midMargin = 0.08;
    CanvasPartitionTwoColumn( ccv, Nx, Ny, lMargin, rMargin, bMargin, tMargin, midMargin );
  }
  
  TPad* pad0 = (TPad*)gROOT->FindObject("pad_0_0");
  
  TString pname = Form("pad_%i_%i",pmap[idiv].first,pmap[idiv].second); 
  TPad* pad = (TPad*)gROOT->FindObject(pname);
  
  if( !pad || !pad0 ) return;
  
  pad -> SetFillStyle(4000);
  pad -> SetFrameFillStyle(4000);
  pad -> cd();
  Float_t xFactor = pad0->GetAbsWNDC()/pad->GetAbsWNDC();
  Float_t yFactor = pad0->GetAbsHNDC()/pad->GetAbsHNDC();
  
  
  

  // division
  // 0 | 1
  //-------
  // 2 | 3
  
  if( bRATIO ) {
    if( idiv == 1 || idiv == 3 ) {
      //format for x axis
      mgrRatio ->GetXaxis()->SetLabelFont(43);
      mgrRatio ->GetXaxis()->SetLabelSize(20);
      mgrRatio ->GetXaxis()->SetLabelOffset(0.01);
      mgrRatio ->GetXaxis()->SetTitleFont(43);
      mgrRatio ->GetXaxis()->SetTitleSize(20);
      mgrRatio ->GetXaxis()->SetTitleOffset(4);
      mgrRatio ->GetXaxis()->CenterTitle();
      mgrRatio ->GetXaxis()->SetNdivisions(505);
      mgrRatio ->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
  
      //format for y axis
      mgrRatio ->GetYaxis()->SetLabelFont(43);
      mgrRatio ->GetYaxis()->SetLabelSize(20);
      mgrRatio ->GetYaxis()->SetLabelOffset(0.02);
      mgrRatio ->GetYaxis()->SetTitleFont(43);
      mgrRatio ->GetYaxis()->SetTitleSize(20);
      mgrRatio ->GetYaxis()->SetTitleOffset(3.5);
      mgrRatio ->GetYaxis()->CenterTitle();
      mgrRatio ->GetYaxis()->SetNdivisions(505);
      mgrRatio ->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
  
      mgrRatio -> GetXaxis() -> SetNdivisions(5,5,0,kTRUE);
      mgrRatio -> GetYaxis() -> SetNdivisions(5,5,0,kTRUE);

      mgrRatio -> Draw("A");      
      cout << " bratio " << bRATIO << " idiv " << idiv << endl;
    }
  }
  else {
    //format for x axis
    mgr[0] ->GetXaxis()->SetLabelFont(43);
    mgr[0] ->GetXaxis()->SetLabelSize(20);
    mgr[0] ->GetXaxis()->SetLabelOffset(0.01);
    mgr[0] ->GetXaxis()->SetTitleFont(43);
    mgr[0] ->GetXaxis()->SetTitleSize(20);
    mgr[0] ->GetXaxis()->SetTitleOffset(4);
    mgr[0] ->GetXaxis()->CenterTitle();
    mgr[0] ->GetXaxis()->SetNdivisions(505);
    mgr[0] ->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
  
    //format for y axis
    mgr[0] ->GetYaxis()->SetLabelFont(43);
    mgr[0] ->GetYaxis()->SetLabelSize(18);
    mgr[0] ->GetYaxis()->SetLabelOffset(0.02);
    mgr[0] ->GetYaxis()->SetTitleFont(43);
    mgr[0] ->GetYaxis()->SetTitleSize(20);
    mgr[0] ->GetYaxis()->SetTitleOffset(3.5);
    mgr[0] ->GetYaxis()->CenterTitle();
    mgr[0] ->GetYaxis()->SetNdivisions(505);
    mgr[0] ->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
  
    mgr[0] -> GetXaxis() -> SetNdivisions(5,5,0,kTRUE);
    mgr[0] -> GetYaxis() -> SetNdivisions(5,5,0,kTRUE);

    mgr[0] -> Draw("alp");      
  }

  if( idiv == 0 ) {
    lgd -> SetX1(0.3);
    lgd -> SetX2(0.8);
    lgd -> SetY1(0.8);
    lgd -> SetY2(0.9);

    lg  -> SetX1(0.3);
    lg  -> SetX2(0.8);
    lg  -> SetY1(0.6);
    lg  -> SetY2(0.8);

    lgd  -> Draw();
    lg   -> Draw();

    TString ssys = "";
    for( auto iisys : ivsys ) {
      ssys += rsys[iisys]+"Sn";
    }
    ssys += lbCentral[isg];
    ccv->SaveAs(Form("partdep_%s_%s.png",phys.Data(),ssys.Data()));  

  }
  else if( idiv == 1 ) {
    lgd -> SetX1(0.25);
    lgd -> SetX2(0.8);
    lgd -> SetY1(0.8);
    lgd -> SetY2(0.9);
    lg  -> SetX1(0.25);
    lg  -> SetX2(0.8);
    lg  -> SetY1(0.6);
    lg  -> SetY2(0.8);
      
    lgd -> Draw();
    lg  -> Draw();
  }


  // else if( ivsys.size() <= 2) {

  //   SetXTextLabel(*mgr[0], sqpart);

  //   if( ndiv == 1 ){
  //     ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),500,400); iccv++;
  //     mgr[0]->Draw("ALP");


  //     lg -> Draw();
  //     lgd-> Draw();

  //     //      plabel.DrawLatexNDC(0.25,0.85, lbCentral[isg]);
  //     TString ssys = "";
  //     for( auto iisys : ivsys ) {
  // 	ssys += rsys[iisys]+"Sn";
  //     }
  //     ssys += lbCentral[isg];
  //     ccv->SaveAs(Form("pdep1_%s_%s.png",phys.Data(),ssys.Data()));
  //   }

  //   else if( ndiv == 2) {

  //     if( idiv == 0) {
  // 	ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),700,500); iccv++;
  // 	ccv -> Divide(ndiv,1);
  //     }

  //     ccv -> cd(idiv+1);
  //     mgr[0]->Draw("ALP");


  //     if( idiv == 0) { 
  // 	lg -> Draw();
  // 	lgd-> Draw();
  //     }
  //     else if( idiv == 1) {
  // 	//	plabel.DrawLatexNDC(0.25,0.85, lbCentral[isg]);
  // 	TString ssys = "";
  // 	for( auto iisys : ivsys ) {
  // 	  ssys += rsys[iisys]+"Sn";
  // 	}
  // 	ssys += lbCentral[isg];
  // 	ccv->SaveAs(Form("pdep2_%s_%s.png",phys.Data(),ssys.Data()));
  //     }
  //   }
  

  //   //-- presen
  //   else if( ndiv == 4 ) {
      
  //     std::vector< std::pair<Int_t,Int_t>> pmap = {{0,1},{1,1},{0,0},{1,0}};
  //     if( idiv == 0) {
  // 	ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),500,600); iccv++;
      
  // 	const Int_t Ny = 2;
  // 	const Int_t Nx = 2;
  // 	Float_t lMargin = 0.13;
  // 	Float_t rMargin = 0.01;
  // 	Float_t bMargin = 0.10;
  // 	Float_t tMargin = 0.02;
  // 	Float_t midMargin = 0.0;
  // 	CanvasPartitionTwoColumn( ccv, Nx, Ny, lMargin, rMargin, bMargin, tMargin, midMargin );
  //     }
          
  //     TPad* pad0 = (TPad*)gROOT->FindObject("pad_0_0");
      
  //     TString pname = Form("pad_%i_%i",pmap[idiv].first,pmap[idiv].second); 
  //     TPad* pad = (TPad*)gROOT->FindObject(pname);
      
  //     if( !pad || !pad0 ) return;
      
  //     pad -> SetFillStyle(4000);
  //     pad -> SetFrameFillStyle(4000);
  //     pad -> cd();
  //     Float_t xFactor = pad0->GetAbsWNDC()/pad->GetAbsWNDC();
  //     Float_t yFactor = pad0->GetAbsHNDC()/pad->GetAbsHNDC();
      

  //     //format for x axis
  //     mgr[0] ->GetXaxis()->SetLabelFont(43);
  //     mgr[0] ->GetXaxis()->SetLabelSize(20);
  //     mgr[0] ->GetXaxis()->SetLabelOffset(0.01);
  //     mgr[0] ->GetXaxis()->SetTitleFont(43);
  //     mgr[0] ->GetXaxis()->SetTitleSize(20);
  //     mgr[0] ->GetXaxis()->SetTitleOffset(4);
  //     mgr[0] ->GetXaxis()->CenterTitle();
  //     mgr[0] ->GetXaxis()->SetNdivisions(505);
  //     mgr[0] ->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
      
  //     //format for y axis
  //     mgr[0] ->GetYaxis()->SetLabelFont(43);
  //     mgr[0] ->GetYaxis()->SetLabelSize(18);
  //     mgr[0] ->GetYaxis()->SetLabelOffset(0.02);
  //     mgr[0] ->GetYaxis()->SetTitleFont(43);
  //     mgr[0] ->GetYaxis()->SetTitleSize(20);
  //     mgr[0] ->GetYaxis()->SetTitleOffset(3.5);
  //     mgr[0] ->GetYaxis()->CenterTitle();
  //     mgr[0] ->GetYaxis()->SetNdivisions(505);
  //     mgr[0] ->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
    
  //     mgr[0] -> GetXaxis() -> SetNdivisions(5,5,0,kTRUE);
  //     mgr[0] -> GetYaxis() -> SetNdivisions(5,5,0,kTRUE);


  //     // division
  //     // 0 | 1
  //     //-------
  //     // 2 | 3
      
  //     if( idiv == 0 )
  // 	mgr[0] -> Draw("alp");      
  //     else if( idiv == 1 ) {
  // 	mgr[0] -> Draw("alp");
  // 	//	mgrRatio -> Draw("alp");      

  //     }      
  //     else if( idiv == 3 ) {
  // 	mgr[0] -> Draw("alp");
  // 	//	mgrRatio -> Draw("alp");
  // 	lgd -> SetX1(0.1);
  // 	lgd -> SetX2(0.8);
  // 	lgd -> SetY1(0.35);
  // 	lgd -> SetY2(0.45);
  // 	lg  -> SetX1(0.1);
  // 	lg  -> SetX2(0.8);
  // 	lg  -> SetY1(0.2);
  // 	lg  -> SetY2(0.35);

  // 	lgd  -> Draw();
  // 	lg   -> Draw();
  //     }
  //     else if( idiv == 2 ) {
  // 	mgr[0] -> Draw("alp");      

  // 	lgd -> SetX1(0.25);
  // 	lgd -> SetX2(0.8);
  // 	lgd -> SetY1(0.35);
  // 	lgd -> SetY2(0.5);
  // 	lg  -> SetX1(0.25);
  // 	lg  -> SetX2(0.8);
  // 	lg  -> SetY1(0.2);
  // 	lg  -> SetY2(0.35);

  // 	lgd -> Draw();
  // 	lg  -> Draw();
  //     }

	
  //     if( idiv == 3 ) {
  // 	TString ssys = "";
  // 	for( auto iisys : ivsys ) {
  // 	  ssys += rsys[iisys]+"Sn";
  // 	}
  // 	ssys += lbCentral[isg];
  // 	ccv->SaveAs(Form("partdep_%s_%s.png",phys.Data(),ssys.Data()));  
  //     }
  //   }
  // }

  // if( ivsys.size() > 2 ) {
  //   ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1000,500); iccv++;
  //   ccv->Divide(3,1);
  //   ccv->cd(1);
  //   lg ->Draw();

  //   for( auto i : {0,1} ){
  //     ccv->cd(i+2);
  //     SetXTextLabel(*mgr[i], sqpart);
  //     mgr[i]->Draw("ALP");
  //   }

  //   TString ssys = "";
  //   for( auto iisys : ivsys ) {
  //     ssys += rsys[iisys]+"Sn";
  //   }
  //   ssys += lbCentral[isg];
  //   ccv->SaveAs(Form("partratio_%s_%s_1.png",phys.Data(),ssys.Data()));  
    
  //   ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),700,500); iccv++;
  //   ccv->Divide(2,1);
  //   ccv->cd(1);
  //   lgr ->Draw();
    
  //   ccv->cd(2);
  //   SetXTextLabel(*mgr[2], sqpart);

  //   mgr[2]->Draw("ALP");
  //   ccv->SaveAs(Form("partratio_%s_%s_2.png",phys.Data(),ssys.Data()));  
  // }
}

void Draw_compAMD(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart)
{
  LOG(INFO) << "[Draw_compAMD]... " ;
  for( auto iisys : ivsys )
    LOG(INFO) << rsys[iisys] << " : " ;
  for( auto ipart : sqpart )
    LOG(INFO) << " >> " << ncls[ipart].sName ;
  LOG(INFO)<< FairLogger::endl; 

  TString fname = Form("physData%s.dat",fphysdatafile[fphysdataid].second.Data());
  Load_physDATA(fname);
  fname = Form("physAMD%s.dat",fphysdatafile[fphysdataid].second.Data());
  Load_physAMD(fname, sqpart);


  const UInt_t nsys  = (UInt_t)ivsys.size();
  const UInt_t npart = (UInt_t)sqpart.size();

  const UInt_t nplot = 2;
  TMultiGraph*  mgrph[nsys][npart][nplot];
  TGraphErrors* grphAMD[nsys][npart][nplot];
  TGraphErrors* grphData[nsys][npart][nplot];
  std::vector< TString > AMDLabel;
  //TMultiGraph*  mgrphHHe[nsys];
  //{"meanpxSlope","meanbtgm"},{"meanpxSlope","pxstdDev"},{"meanbtgm","meanEt"}, {"integralbtgm" ,"integralyx"}};
  TString physpara[nplot] = {"numClust","integralM"};
    //{"meanpxSlope","integralM"};
  TGraphErrors* grphAMDSub[nsys][npart];


  UInt_t jsys = 0;
  for( auto isys : ivsys ){

    UInt_t jpart = 0;
    for( auto ipart: sqpart ){
      for( auto jpara : ROOT::TSeqI(nplot) ) {
	grphAMD[jsys][jpart][jpara]  = new TGraphErrors();
	grphData[jsys][jpart][jpara] = new TGraphErrors();
	mgrph[jsys][jpart][jpara]    = new TMultiGraph();
      }
      grphAMDSub[jsys][jpart] = new TGraphErrors();

      mgrph[jsys][jpart][0] -> SetTitle(lsys[isys]+"("+ncls[ipart].sName+");;"+ physpara[0]);
      mgrph[jsys][jpart][1] -> SetTitle(";; "+ physpara[1]);

      AMDLabel.clear();
      Double_t value, valuee;
      Double_t max[2] = {-999.,-999.};
      Double_t min[2] = { 999., 999.};

      for( auto ig : ROOT::TSeqI( AMDnames.size() ) ) {
	if( physAMD[isys][ipart][ig][0].system != 0 ){

	  for( auto j : ROOT::TSeqI(nplot) ) {
	    value = 0.; valuee = 0.;


	    // upper panel
	    LOG(INFO) << "[AMD compare [" << physpara[j] << "] " << physAMD[isys][ipart][ig][0].version << " " << physAMD[isys][ipart][ig][0].system << FairLogger::endl;

	    if( physpara[j] == "numClust" ) {
	      Double_t numClust[9][2];
	      Draw_cluster({isys},{ipart}, AMDnames.at(ig) ,numClust);

	      value  = numClust[ipart][0];
	      valuee = numClust[ipart][1];
	    }
	    else if( physpara[j] == "meanpxSlope") {
	      value  = physAMD[isys][ipart][ig][0].meanpxSlope;
	      valuee = physAMD[isys][ipart][ig][0].meanpxSlopeError ;
	    }
	    else if( physpara[j] == "meanbtgm") {
	      value  = physAMD[isys][ipart][ig][1].meanbtgm;
	      valuee = physAMD[isys][ipart][ig][1].meanbtgmError; 
	    }
	    else if( physpara[j] == "integralbtgm") {
	      value  = physAMD[isys][ipart][ig][1].integralbtgm;
	      valuee = physAMD[isys][ipart][ig][1].integralbtgmError; 
	    }
	    else if( physpara[j] == "integralM") {
	      value  = physAMD[isys][ipart][ig][1].integralM;
	      valuee = physAMD[isys][ipart][ig][1].integralMError; 
	    }
	    else if( physpara[j] == "meanEt") {
	      value = physAMD[isys][ipart][ig][1].meanEt;
	      valuee = physAMD[isys][ipart][ig][1].meanEtError;
	    }	  
	    else if( physpara[j] == "pxstdDev") {
	      value  = physAMD[isys][ipart][ig][1].pxstdDev;
	      valuee = physAMD[isys][ipart][ig][1].pxstdDevError;
	    }
	    else if( physpara[j] == "integralyx") {
	      value  = physAMD[isys][ipart][ig][1].integralyx;
	      valuee = physAMD[isys][ipart][ig][1].integralyxError;
	    }
	    else if( physpara[j] == "meanbtgm") {
	      value  = physAMD[isys][ipart][ig][1].meanbtgm;
	      valuee = physAMD[isys][ipart][ig][1].meanbtgmError; 
	    }


	    grphAMD[jsys][jpart][j] -> SetPoint(grphAMD[jsys][jpart][j]->GetN(), grphAMD[jsys][jpart][j]->GetN()+1,  value );
	    grphAMD[jsys][jpart][j] -> SetPointError(grphAMD[jsys][jpart][j]->GetN()-1,                          0,  valuee );

	    if( max[j] < value )
	      max[j] =   value + valuee *2.;
	    if( min[j] > value )
	      min[j] =   value - valuee *2.;

	    //	  AMDLabel.push_back( physAMD[isys][ipart][ig][0].version(5,20) );
	    if( j == 0 )
	      AMDLabel.push_back( physAMD[isys][ipart][ig][0].version );
	  
	    grphAMD[jsys][jpart][j] -> SetMarkerStyle(20);
	    grphAMD[jsys][jpart][j] -> SetMarkerColor(2);
	    grphAMD[jsys][jpart][j] -> SetLineColor(2);
	    
	    mgrph[jsys][jpart][j] -> Add(grphAMD[jsys][jpart][j],"P" );
	  	    
	  }
	}
      }
 
      if( kFALSE && physpara[0] == "numClust" && (grphAMD[jsys][jpart][0]->GetN() ==  grphAMD[jsys][jpart][1]->GetN()) ) {
	grphAMD[jsys][jpart][0] -> SetMarkerColor(4);
	mgrph[jsys][jpart][1] -> Add( grphAMD[jsys][jpart][0], "P");
      }

      if(physpara[0] == "numClust" && (grphAMD[jsys][jpart][0]->GetN() ==  grphAMD[jsys][jpart][1]->GetN()) ) {

	for( auto k : ROOT::TSeqI(grphAMD[jsys][jpart][0]->GetN()) ) {
	  Double_t x0, y0, x0e, y0e;
	  Double_t x1, y1, x1e, y1e;
	  
	  grphAMD[jsys][jpart][0]->GetPoint(k, x0, y0);
	  grphAMD[jsys][jpart][1]->GetPoint(k, x1, y1);
	  y0e = grphAMD[jsys][jpart][0]->GetErrorY(k);
	  y1e = grphAMD[jsys][jpart][1]->GetErrorY(k);
	  auto err = sqrt(y0e*y0e+y1e+y1e);
	  ///	   cout << " k " << k << " = " << y0 << " " << y1 << " y0-y1 " << y0-y1 << endl; 
	  grphAMDSub[jsys][jpart]->SetPoint(k, x0, y0-y1);
	  grphAMDSub[jsys][jpart]->SetPointError(k, 0, err);

	}	
      }
      

      //++++ DATA ++++++++++++++
      if( jpart < 5 ) {

	for( auto j : ROOT::TSeqI(nplot) ) {
	  value = 0.; valuee = 0.;
	  if( physpara[j] == "meanpxSlope") {
	    value  = physData[isys][ipart][0].meanpxSlope;
	    valuee = physData[isys][ipart][0].meanpxSlopeError;
	  }
	  else if( physpara[j] == "integralM") {
	    value  = physData[isys][ipart][1].integralM;
	    valuee = physData[isys][ipart][1].integralMError;
	  }
	  else if( physpara[j] == "meanbtgm") {
	    value  = physData[isys][ipart][1].meanbtgm;
	    valuee = physData[isys][ipart][1].meanbtgmError;
	  }
	  else if( physpara[j] == "integralyx") {
	    value  = physData[isys][ipart][1].integralyx;
	    valuee = physData[isys][ipart][1].integralyxError;
	  }
	  else if( physpara[j] == "meanEt") {
	    value  = physData[isys][ipart][1].meanEt;
	    valuee = physData[isys][ipart][1].meanEtError;
	  }
	  else if( physpara[j] == "pxstdDev") {
	    value  = physData[isys][ipart][1].pxstdDev;
	    valuee = physData[isys][ipart][1].pxstdDevError;
	  }
	  else if( physpara[j] == "integralbtgm") {
	    value  = physData[isys][ipart][1].integralbtgm;
	    valuee = physData[isys][ipart][1].integralbtgmError;
	  }
	

	  grphData[jsys][jpart][j] -> SetPoint(0, 0.5,    value );
	  grphData[jsys][jpart][j] -> SetPointError(0, 0, valuee );
	  grphData[jsys][jpart][j] -> SetPoint(1, AMDLabel.size()+1.5, value );
	  grphData[jsys][jpart][j] -> SetPointError(1, 0, valuee );
	
	  if( max[j] < value )
	    max[j] =   value + valuee *2.;
	  if( min[j] > value )
	    min[j] =   value - valuee *2.;

	  grphData[jsys][jpart][j] -> SetFillStyle(3001);
	  grphData[jsys][jpart][j] -> SetFillColor(kGreen-3);
	  mgrph[jsys][jpart][j]    -> Add(grphData[jsys][jpart][j], "A3" );

	  mgrph[jsys][jpart][j]    -> GetYaxis() -> SetRangeUser(min[j], max[j]);
	}
      }
      jpart++;
    }
    jsys++;
  }

  jsys = 0;
  for( auto isys : ivsys ) {

    UInt_t jpart = 0;
    for( auto ipart: sqpart ){

      mgrph[jsys][jpart][0] -> GetXaxis()->SetLimits(0.5,AMDLabel.size()+0.5);
      mgrph[jsys][jpart][1] -> GetXaxis()->SetLimits(0.5,AMDLabel.size()+0.5);
      for( auto i : ROOT::TSeqI(AMDLabel.size()) ) { 
	cout << " AMDLabel " << AMDLabel.at(i) << endl;
	mgrph[jsys][jpart][0]->GetXaxis()->SetBinLabel( mgrph[jsys][jpart][0]->GetXaxis()->FindBin( i+1 ), AMDLabel.at(i) );
	mgrph[jsys][jpart][1]->GetXaxis()->SetBinLabel( mgrph[jsys][jpart][1]->GetXaxis()->FindBin( i+1 ), AMDLabel.at(i) );

	if( grphAMDSub[jsys][jpart]->GetN()>0 )
	  grphAMDSub[jsys][jpart]->GetXaxis()->SetBinLabel( grphAMDSub[jsys][jpart]->GetXaxis()->FindBin( i+ 1 ), AMDLabel.at(i) );
      }

      ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),700, 600); iccv++;
      CanvasPartitionY(ccv, 2, 0.02, 0.02, 0.2, 0.05);

      ccv->cd(0);
      TString pname = Form("pad_%i",0);
      TPad* pad0 = (TPad*)gROOT->FindObject(pname);
      if( pad0 ) {
	pad0->Draw();
	pad0->SetFillStyle(4000);
	pad0->SetFrameFillStyle(4000);
	pad0->SetGrid();
	pad0->cd();
	mgrph[jsys][jpart][1] -> GetYaxis() -> SetTitleSize(0.06);

	mgrph[jsys][jpart][1] -> Draw("AP");
	plabel.SetTextSize(0.08);
	plabel.DrawLatexNDC(0.8, 0.5, "b:"+lbCentral[1] );
      }

      ccv->cd(0);
      TPad* pad1 = (TPad*)gROOT->FindObject("pad_1");
      if( pad1 ) {
	pad1->Draw();
	pad1->SetFillStyle(4000);
	pad1->SetFrameFillStyle(4000);
	pad1->SetGrid();
	pad1->cd();

	Float_t yFactor = pad0->GetAbsHNDC()/pad1->GetAbsHNDC();
	//	cout << yFactor << endl;
	mgrph[jsys][jpart][0] -> GetYaxis() -> SetTitleSize(0.06*yFactor);
	mgrph[jsys][jpart][0] -> Draw("AP");
	plabel.SetTextSize(0.08*yFactor);
	plabel.DrawLatexNDC(0.8, 0.1, "b:"+lbCentral[0] );
      }

      ccv->SaveAs(Form("cmpAMD_%dSn_%s_%s_%s.png",Asys[isys],ncls[ipart].sName.Data(),physpara[0].Data(),physpara[1].Data()));


      ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
      grphAMDSub[jsys][jpart] -> SetTitle(Form("%u %s",Asys[isys],ncls[ipart].sName.Data()));
      grphAMDSub[jsys][jpart] -> SetMarkerStyle(20);
      grphAMDSub[jsys][jpart] -> SetMarkerColor(2);
      grphAMDSub[jsys][jpart] -> Draw("AP");
      ccv->SaveAs(Form("compdeltM_%dSn_%s_%s_%s.png",Asys[isys],ncls[ipart].sName.Data(),physpara[0].Data(),physpara[1].Data()));



      jpart++;
    }
    jsys++;
  }
}

void Draw_compParameter(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, UInt_t isg=0, TString para = "slope")
{
  LOG(INFO) << "[Draw_compParameter]... " << para << " central " << isg << " " ;
  for( auto iisys : ivsys )
    LOG(INFO) << rsys[iisys] << " : " ;
  LOG(INFO)<< FairLogger::endl; 
  
  Bool_t bAMD = 1;
  if( ivsys.size() > 1 ) bAMD=0;

  TString fname = Form("physData%s.dat",fphysdatafile[fphysdataid].second.Data());
  Load_physDATA(fname);
  if( bAMD ) {
    fname = Form("physAMD%s.dat",fphysdatafile[fphysdataid].second.Data());
    Load_physAMD(fname,sqpart);
  }

  //+++ <px>/M
  auto mgr = new TMultiGraph();
  auto lg  = new TLegend(0.1,0.2,0.8,0.8,lbCentral[isg]);

  if( para == "meanpxSlope" ) {
    mgr -> SetTitle(";Particle;Slope of <px>/A");
  }
  else if( para == "meanEt" ) {
    mgr -> SetTitle(";Particle; <Et>");
  }
  else if( para == "meanbtgm" ) {
    mgr -> SetTitle(";Particle; <#beta #gamma>");
  }
  else if( para == "pxstdDev" ) {
    mgr -> SetTitle(";Particle; pxstd_Dev");
  }
  else if( para == "integralM" ) {
    mgr -> SetTitle(";Particle; Multiplicity");
  }
  else if( para == "RapiditystdDev" ) {
    mgr -> SetTitle(";Particle; Rapidity_Dev");
  }
  else if( para == "xRapiditystdDev" ) {
    mgr -> SetTitle(";Particle; xRapidity_Dev");
  }


  TGraphErrors* grph;

  // AMD
  if( bAMD ) {
    for( auto ig : ROOT::TSeqI( AMDnames.size() ) ) {
      for( auto iisys : ivsys ){
	if( physAMD[iisys][0][ig][isg].system != 0 && AMDnames.at(ig).category == 1) {
	  LOG(INFO) << "[AMD compare] " << physAMD[iisys][0][ig][isg].version << " " << physAMD[iisys][0][ig][isg].system << FairLogger::endl;
	  grph = new TGraphErrors();
	  grph -> SetName(Form("phAMD_%d_%d",iisys,ig) );
	  
	  std::vector< UInt_t >::iterator itpart = sqpart.begin();
	  for( auto iin : ROOT::TSeqI(sqpart.size()) ){
	    Double_t x, xe;
	    if( para == "meanpxSlope" ) {
	      x  = physAMD[iisys][*itpart][ig][isg].meanpxSlope;
	      xe = physAMD[iisys][*itpart][ig][isg].meanpxSlopeError;
	    }
	    else if( para == "meanEt" ) {
	      x  = physAMD[iisys][*itpart][ig][isg].meanEt;
	      xe = physAMD[iisys][*itpart][ig][isg].meanEtError;
	    }
	    else if( para == "meanbtgm" ) {
	      x  = physAMD[iisys][*itpart][ig][isg].meanbtgm;
	      xe = physAMD[iisys][*itpart][ig][isg].meanbtgmError;
	    }
	    else if( para == "pxstdDev" ) {
	      x  = physAMD[iisys][*itpart][ig][isg].pxstdDev;
	      xe = physAMD[iisys][*itpart][ig][isg].pxstdDevError;
	    }
	    else if( para == "integralM" ) {
	      x  = physAMD[iisys][*itpart][ig][isg].integralM;
	      xe = physAMD[iisys][*itpart][ig][isg].integralMError;
	    }
	    else if( para == "RapiditystdDev" ) {
	      x  = physAMD[iisys][*itpart][ig][isg].RapiditystdDev;
	      xe = physAMD[iisys][*itpart][ig][isg].RapiditystdDevError;
	    }
	    else if( para == "xRapiditystdDev" ) {
	      x  = physAMD[iisys][*itpart][ig][isg].xRapiditystdDev;
	      xe = physAMD[iisys][*itpart][ig][isg].xRapiditystdDevError;
	    }




	    grph -> SetPoint     (iin, iin+1, x  );
	    grph -> SetPointError(iin,     0, xe );
	    itpart++;
	  }

	  grph->SetFillStyle( AMDnames.at(ig).fStyle );
	  grph->SetLineColor( AMDnames.at(ig).fColor);
	  grph->SetFillColor( AMDnames.at(ig).fColor);

	  mgr->Add(grph, "3");
	  lg ->AddEntry(grph,"AMD "+AMDnames.at(ig).config); 
	} 
      }
    }
  }

  // data
  for( auto iisys : ivsys ) {
    grph = new TGraphErrors();
    grph -> SetName(Form("phpara_%d",iisys) );

    std::vector< UInt_t >::iterator itpart = sqpart.begin();
    for( auto iin : ROOT::TSeqI(sqpart.size()) ){
      if( *itpart == 7 ) continue;
      Double_t x, xe;
      if( para == "meanpxSlope" ) {
	x  = physData[iisys][*itpart][isg].meanpxSlope ;
	xe = physData[iisys][*itpart][isg].meanpxSlopeError ;
      }
      else if( para == "meanEt" ) {
	x  = physData[iisys][*itpart][isg].meanEt ;
	xe = physData[iisys][*itpart][isg].meanEtError ;
      }
      else if( para == "meanbtgm" ) {
	x  = physData[iisys][*itpart][isg].meanbtgm ;
	xe = physData[iisys][*itpart][isg].meanbtgmError ;
      }
      else if( para == "meanKt" ) {
	x  = 0.5*pow(physData[iisys][*itpart][isg].meanbtgm,2)*mass[*itpart] ;
	xe = physData[iisys][*itpart][isg].meanbtgmError ;
      }
      else if( para == "pxstdDev" ) {
	x  = physData[iisys][*itpart][isg].pxstdDev ;
	xe = physData[iisys][*itpart][isg].pxstdDevError ;
      }
      else if( para == "integralM" ) {
	x  = physData[iisys][*itpart][isg].integralM;
	xe = physData[iisys][*itpart][isg].integralMError;
      }
      else if( para == "RapiditystdDev" ) {
	x  = physData[iisys][*itpart][isg].RapiditystdDev;
	xe = physData[iisys][*itpart][isg].RapiditystdDevError;
      }
      else if( para == "xRapiditystdDev" ) {
	x  = physData[iisys][*itpart][isg].xRapiditystdDev;
	xe = physData[iisys][*itpart][isg].xRapiditystdDevError;
      }
      
      grph -> SetPoint     (iin, iin+1, x );
      grph -> SetPointError(iin,     0, xe);
      itpart++;
    }

    grph->SetMarkerStyle( DStyle[iisys].mStyle );
    grph->SetMarkerSize(  DStyle[iisys].mSize );
    grph->SetLineColor(   DStyle[iisys].fColor );
    grph->SetMarkerColor( DStyle[iisys].fColor );

    mgr->Add(grph, "p");
    lg ->AddEntry(grph,fsys[iisys]+"_"+physData[iisys][0][isg].version); 
  }  
  
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),800,500); iccv++; 
  ccv->Divide(2,1);
  ccv->cd(1);
  lg ->Draw();

  ccv->cd(2);
  SetXTextLabel(*mgr, sqpart);
  mgr->Draw("A");

  TString ssys = "";
  for( auto iisys : ivsys ) {
    ssys += rsys[iisys]+"Sn";
  }
  ssys += lbCentral[isg];
  ccv->SaveAs(Form("cmp%s_%s.png",para.Data(),ssys.Data()));
}

void Draw_compIntegral(std::vector<UInt_t> ivsys, std::vector<UInt_t> sqpart, UInt_t isg, UInt_t icateg )
{
  LOG(INFO) << "[Draw_compIntegral]... " ;
  Double_t neutAlpha[4] = {1.2, 1.11, 1.15, 1.15};
  Double_t neutAlphae   = 0.05;

  const UInt_t npart = (UInt_t)sqpart.size();
  TString sfsys;
  for( auto iisys : ivsys ) {
    LOG(INFO) << rsys[iisys] << " : " ;
    sfsys += fsys[iisys] + " ";
  }
  LOG(INFO)<< FairLogger::endl; 

  TString fname;
  Bool_t fload = kTRUE;

  fname = Form("physData%s_left.dat",fphysdatafile[fphysdataid].second.Data());
  fload *= Load_physDATA(fname);
  if( !fload ) return;

  fname = Form("physAMD%s.dat",fphysdatafile[fphysdataid].second.Data());
  fload *= Load_physAMD(fname,sqpart);
  if( !fload ) return;      

  const UInt_t nAMDsys = 3;
  std::vector< std::vector< UInt_t >> amdsys(nAMDsys);
  UInt_t MarkerStyle[nAMDsys];
  UInt_t jamd = 0; 
  amdsys[jamd] = {14,15,16};    MarkerStyle[jamd]=21;  jamd++;
  amdsys[jamd] = {3,4,5};       MarkerStyle[jamd]=29;  jamd++;
  amdsys[jamd] = {0,1,2};       MarkerStyle[jamd]=25;  jamd++;
 
  TGraphErrors* grpht3heAMDsys[nAMDsys]; 
  TGraphErrors* grphnpAMDsys[nAMDsys];   
  TGraphErrors* grphDRAMDsys[nAMDsys];   
  TGraphErrors* grpht3heAMDsys0[nAMDsys]; 
  TGraphErrors* grphnpAMDsys0[nAMDsys];   
  TGraphErrors* grphDRAMDsys0[nAMDsys];   
  for( auto i : ROOT::TSeqI(nAMDsys) ) {
    grpht3heAMDsys[i]  = new  TGraphErrors();
    grphnpAMDsys[i]    = new  TGraphErrors();
    grphDRAMDsys[i]    = new  TGraphErrors();
    grpht3heAMDsys0[i] = new  TGraphErrors();
    grphnpAMDsys0[i]   = new  TGraphErrors();
    grphDRAMDsys0[i]   = new  TGraphErrors();
  }

  auto grphnpDatasys   = new  TGraphErrors();
  auto grpht3heDatasys = new  TGraphErrors();
  auto grphDRDatasys   = new  TGraphErrors();

  auto mgr   = new TMultiGraph();
  auto mgrA  = new TMultiGraph();
  auto mgrZ  = new TMultiGraph();
  auto mgrx  = new TMultiGraph();
  auto mgrxA = new TMultiGraph();
  auto mgnp   = new TMultiGraph();
  auto mgt3he = new TMultiGraph();
  auto mgDR   = new TMultiGraph();
  auto mgnpsys   = new TMultiGraph();
  auto mgt3hesys = new TMultiGraph();
  auto mgDRsys   = new TMultiGraph();
  std::vector< TString > AMDLabel;  

  auto grphnpAMD    = new  TGraphErrors();
  auto grpht3heAMD  = new  TGraphErrors();
  auto grphnpAMD0   = new  TGraphErrors();
  auto grpht3heAMD0 = new  TGraphErrors();
  auto grphnpData   = new  TGraphErrors();
  auto grpht3heData = new  TGraphErrors();

  auto grphDRAMD     = new  TGraphErrors();
  auto grphDRAMD0    = new  TGraphErrors();
  auto grphDRData    = new  TGraphErrors();


  auto lg  = new TLegend(0.15,0.3,0.9,0.9, sfsys+" "+lbCentral[isg]);
  mgr  -> SetTitle(";Particle; Multiplicity");
  mgrA -> SetTitle(";Particle; A*Multiplicity");
  mgrZ -> SetTitle(";Particle; Z*Multiplicity");
  mgrx -> SetTitle("|y|<0.3;Particle; Multiplicity");
  mgrxA-> SetTitle("|y|<0.3;Particle; A*Multiplicity");


  TGraphErrors* g_mH[4];
  g_mH[0] = new TGraphErrors();
  g_mH[0] ->SetName("g_mZ");
  g_mH[0] ->SetTitle(";AMD;M(p+d+t+2*(^{3}He+^{4}He))");
  g_mH[1] = new TGraphErrors();
  g_mH[1] ->SetName("g_mZx");
  g_mH[1] ->SetTitle("|y|<0.3;AMD;M(p+d+t+2*(^{3}He+^{4}He))");
  g_mH[2] = new TGraphErrors();
  g_mH[2] ->SetName("g_mH");
  g_mH[2] ->SetTitle(";AMD;M(p+d+t)");
  g_mH[3] = new TGraphErrors();
  g_mH[3] ->SetName("g_mHx");
  g_mH[3] ->SetTitle("|y|<0.3;AMD;M(p+d+t)");

  std::vector< std::vector<Double_t> > AMDInteglralM(2);
  TGraphErrors* grph;
  TGraphErrors* grphN;
  TGraphErrors* grphA;
  TGraphErrors* grphZ;
  TGraphErrors* grphx;
  TGraphErrors* grphxA;
  TGraphErrors* grphO;
  TGraphErrors* grphOA;
  TBox* aBox[4][5];

  // AMD

  UInt_t inx = 0;
  for( auto iisys : ivsys ){
    AMDLabel.clear();
    Double_t value, valuee;
    Double_t max[2] = {-999.,-999.};
    Double_t min[2] = { 999., 999.};

    UInt_t iamd = 0;
    for( auto igg : ROOT::TSeqI( AMDnames.size() ) ) {
      UInt_t ig = AMDnames.at(igg).id;

      if( physAMD[iisys][0][ig][isg].system != 0 ) {
	LOG(INFO) << "[AMD comparison] " << ig << " " << physAMD[iisys][0][ig][isg].version << FairLogger::endl;

	grph   = new TGraphErrors();
	grph   -> SetName(Form("g_phAMD_%d_%d",iisys,ig) );
	grphN  = new TGraphErrors();
	grphN  -> SetName(Form("g_phAMDN_%d_%d",iisys,ig) );
	grphA  = new TGraphErrors();
	grphA  -> SetName(Form("g_phAMDZ_%d_%d",iisys,ig) );
	grphZ  = new TGraphErrors();
	grphZ  -> SetName(Form("g_phAMDZ_%d_%d",iisys,ig) );
	grphx  = new TGraphErrors();
	grphx  -> SetName(Form("g_phAMDx_%d_%d",iisys,ig) );
	grphxA = new TGraphErrors();
	grphxA -> SetName(Form("g_phAMDxA_%d_%d",iisys,ig) );


	grphO = new TGraphErrors();
	grphO -> SetName(Form("g_AMDO_%d_%d",iisys,ig));
	grphOA= new TGraphErrors();
	grphOA-> SetName(Form("g_AMDOA_%d_%d",iisys,ig));

	Double_t numClust[9][2];
	Draw_cluster({iisys},sqpart, AMDnames.at(ig) ,numClust);
	

	// M(p+d+t+2(3He+4He)
	//

	Double_t multH[npart], multHe[npart];
	for( auto i : ROOT::TSeqI(npart) ) {multH[i] = 0.; multHe[i]=0.;}

	Double_t multn = 0, multne = 0.;
	Double_t multp = 0, multpe = 0.;
	for( auto ipart : ROOT::TSeqI(5) ){
	  if( ipart < 3 ) {
	    multH[2]  += physAMD[iisys][ipart][ig][isg].integralM*ncls[ipart].Z;
	    multHe[2] += pow( physAMD[iisys][ipart][ig][isg].integralMError*ncls[ipart].Z ,2);
	    multH[3]  += physAMD[iisys][ipart][ig][isg].integralyx*ncls[ipart].Z;
	    multHe[3] += pow( physAMD[iisys][ipart][ig][isg].integralyxError*ncls[ipart].Z ,2);
	  }

	  if( ipart < 5 ) {
	    multn  += physAMD[iisys][ipart][ig][isg].integralM* (ncls[ipart].A - ncls[ipart].Z);
	    multne += pow( physAMD[iisys][ipart][ig][isg].integralMError * (ncls[ipart].A - ncls[ipart].Z), 2);
	    
	    multp  += physAMD[iisys][ipart][ig][isg].integralM* ncls[ipart].Z;
	    multpe += pow( physAMD[iisys][ipart][ig][isg].integralMError* ncls[ipart].Z, 2);
	  }
	}

	UInt_t iin = 0;
	for( auto ipart : sqpart ){
	  if( ipart == 5 || ipart == 6 ) continue;

	  grph   -> SetPoint     (iin,   iin+1, physAMD[iisys][ipart][ig][isg].integralM );
	  grph   -> SetPointError(iin,     0.5, physAMD[iisys][ipart][ig][isg].integralMError );

	  grphA  -> SetPoint     (iin,   iin+1, physAMD[iisys][ipart][ig][isg].integralM*ncls[ipart].A );
	  grphA  -> SetPointError(iin,     0.5, physAMD[iisys][ipart][ig][isg].integralMError*ncls[ipart].A );

	  grphZ  -> SetPoint     (iin,   iin+1, physAMD[iisys][ipart][ig][isg].integralM*ncls[ipart].Z );
	  grphZ  -> SetPointError(iin,     0.5, physAMD[iisys][ipart][ig][isg].integralMError*ncls[ipart].Z );

	  grphx  -> SetPoint     (iin,   iin+1, physAMD[iisys][ipart][ig][isg].integralyx );
	  grphx  -> SetPointError(iin,     0.5, physAMD[iisys][ipart][ig][isg].integralyxError);

	  grphxA -> SetPoint     (iin,   iin+1, physAMD[iisys][ipart][ig][isg].integralyx*ncls[ipart].A );
	  grphxA -> SetPointError(iin,     0.5, physAMD[iisys][ipart][ig][isg].integralyxError*ncls[ipart].A);

	  grphO  -> SetPoint     (iin,   iin+1, numClust[ipart][0]);
	  grphO  -> SetPointError(iin,     0.5, numClust[ipart][1]);
	  grphOA -> SetPoint     (iin,   iin+1, numClust[ipart][0] * ncls[ipart].A);
	  grphOA -> SetPointError(iin,     0.5, numClust[ipart][1] * ncls[ipart].A);
	  iin++;

	  multH[0]  += physAMD[iisys][ipart][ig][isg].integralM*ncls[ipart].Z;
	  multHe[0] += pow( physAMD[iisys][ipart][ig][isg].integralMError*ncls[ipart].Z ,2);
	  multH[1]  += physAMD[iisys][ipart][ig][isg].integralyx*ncls[ipart].Z;
	  multHe[1] += pow( physAMD[iisys][ipart][ig][isg].integralyxError*ncls[ipart].Z ,2);

	  if( ipart == 7 ) {
	    grphN  -> SetPoint (0, ipart+1, sysN[iisys]-(100.- multp) - multn );
	    grphO  -> SetPoint     (iin,   iin+1, numClust[ipart][0]);
	    grphO  -> SetPointError(iin,     0.5, numClust[ipart][1]);
	    grphOA -> SetPoint     (iin,   iin+1, numClust[ipart][0] * ncls[ipart].A);
	    grphOA -> SetPointError(iin,     0.5, numClust[ipart][1] * ncls[ipart].A);
	  }
	}
	
	multne = sqrt( multne );
	multpe - sqrt( multpe );

	Double_t z2 = 0., z2e = 0.;
	Double_t n2 = 0., n2e = 0.;
	for( auto i : {0,1,2,3,4,7} ){
	  cout << setw(4) << ncls[i].sName << setw(8)<<setprecision(4) << numClust[i][0] << " +- " << setw(6)<< setprecision(4) << numClust[i][1];

	  z2  += numClust[i][0] * ncls[i].Z;
	  z2e += pow(numClust[i][1] * ncls[i].Z, 2);

	  if( i != 7 ) {
	    n2  += numClust[i][0] * (ncls[i].A - ncls[i].Z);
	    n2e += pow(numClust[i][1] * (ncls[i].A - ncls[i].Z), 2);
	  }
	}
	cout << endl;

	z2e = sqrt(z2e);
	n2e = sqrt(n2e);

	Double_t YnAmt4 = sysN[iisys] - (n2 + numClust[7][0]);
	Double_t YpAmt4 = 100 - z2;

	cout << " ---- alpha ---- " << endl;	
	cout << physAMD[iisys][0][ig][isg].version  << " ===> " ;
	cout << " Yn(free) "   << setw(6) << setprecision(4) << numClust[7][0] ;
	cout << " Yn(1<A<=4) " << setw(6) << setprecision(4) << n2 ;
	cout << " Yn(A<4)    " << setw(6) << setprecision(4) << YnAmt4 ;
	cout << " Yp(1<A<=4) " << setw(6) << setprecision(4) << z2 ;
	cout << " Yp(A<4)    " << setw(6) << setprecision(4) << YpAmt4 ;
	cout << " Yn/Yp  ->  " << setw(6) << setprecision(4) << YnAmt4/YpAmt4 << endl;


	Double_t np    = numClust[7][0]/numClust[0][0];
	Double_t npe   = np * sqrt( pow(numClust[7][1]/numClust[7][0],2) + pow(numClust[0][1]/numClust[0][0],2) );
	Double_t nest  = sysN[iisys]- (n2 + 100.-z2 );
	Double_t neste = sqrt(n2e*n2e + z2e+z2e);
	Double_t nestp = nest/numClust[0][0];
	Double_t nestpe= nestp * sqrt( pow(neste/nest,2) + pow(numClust[0][1]/numClust[0][0],2) );
	Double_t t3He  = numClust[2][0]/numClust[3][0];
	Double_t t3Hee = numClust[2][0]/numClust[3][0] * sqrt( pow(numClust[2][1]/numClust[2][0],2) + pow(numClust[3][1]/numClust[3][0],2) ); 

	/// Y(t/3He)/Y(n/p)
	grphnpAMD0   -> SetPoint( AMDLabel.size(),   AMDLabel.size()+1, np  );
	grphnpAMD0   -> SetPointError( AMDLabel.size(),              0, npe );
	grpht3heAMD0 -> SetPoint( AMDLabel.size(), AMDLabel.size()+1, t3He   );
	grpht3heAMD0 -> SetPointError( AMDLabel.size(),            0, t3Hee  );


	// Double_t dr  = t3He / nestp;
	// Double_t dre = dr * sqrt( pow(t3Hee/t3He,2) + pow(nestpe/nestp,2) ); 
	
	cout << " AMD DR is direct output " << endl; 
	Double_t dr  = t3He / np;
	Double_t dre = dr * sqrt( pow(t3Hee/t3He,2) + pow(npe/np,2) );

	grphDRAMD0   -> SetPoint( AMDLabel.size(),   AMDLabel.size()+1, dr );
	grphDRAMD0   -> SetPointError( AMDLabel.size(),              0, dre);


	//	if( ivsys.size() >= 2 && (ig==amdsys[0]||ig==amdsys[1]||ig==amdsys[2] || ig==amdsys[3]||ig==amdsys[4]) ) {
	if( ivsys.size() > 2 ){

	  Int_t isys = -1;
	  for( auto jsys : ROOT::TSeqI(nAMDsys) )for( auto iamd : ROOT::TSeqI(amdsys[jsys].size()) ){
	    
	      if( ig==amdsys[jsys][iamd] ) {
		isys = jsys;
		break;
	      }
	    }

	  cout << ig << " --- " << isys << endl;
	  if( isys > -1 ) { 
	    grphnpAMDsys0[isys]   -> SetPoint( grphnpAMDsys0[isys]->GetN(),    sysNA[iisys], np  );
	    grphnpAMDsys0[isys]   -> SetPointError(grphnpAMDsys0[isys]->GetN(),           0., npe ); 
	    grpht3heAMDsys0[isys] -> SetPoint( grpht3heAMDsys0[isys]->GetN(),  sysNA[iisys], t3He  );
	    grpht3heAMDsys0[isys] -> SetPointError(grpht3heAMDsys0[isys]->GetN(),           0., t3Hee ); 
	    grphDRAMDsys0[isys]   -> SetPoint( grphDRAMDsys0[isys]->GetN(),    sysNA[iisys], dr  );
	    grphDRAMDsys0[isys]   -> SetPointError(grphDRAMDsys0[isys]->GetN(),           0., dre ); 
	  }    
	}


	// OUTPUT
	cout << setw(25)  << physAMD[iisys][0][ig][isg].version << " : ";
	cout << "n/p "      << setw(6) << setprecision(4) << numClust[7][0]/numClust[0][0] << "," << npe << ", ";
	cout << "n(est)/p " << setw(6) << setprecision(4) << nestp << "," << nestpe << ", ";
	cout << "t/3He "   << setw(6)   << setprecision(4) << t3He << "+-" << t3Hee << endl;
	cout << "Integ --- " << endl;

	nest  = sysN[iisys]-(multn + 100.- multp); 
	neste = sqrt( multne*multne + multpe*multpe );
	nestp = nest/physAMD[iisys][0][ig][isg].integralM;
	nestpe= nestp * sqrt( pow(neste/nest,2) + pow(physAMD[iisys][0][ig][isg].integralMError/physAMD[iisys][0][ig][isg].integralM,2) );
	t3He  = physAMD[iisys][2][ig][isg].integralM / physAMD[iisys][3][ig][isg].integralM;
	t3Hee = t3He * sqrt( pow(physAMD[iisys][2][ig][isg].integralMError/physAMD[iisys][2][ig][isg].integralM,2) +
			     pow(physAMD[iisys][3][ig][isg].integralMError/physAMD[iisys][3][ig][isg].integralM,2) );
	
	cout << "n(est) " << setw(6)  << setprecision(4) << nest << "+-" << neste << ", ";
	cout << "p      " << setw(6)  << setprecision(4) << physAMD[iisys][0][ig][isg].integralM << "+-" 
	     << physAMD[iisys][0][ig][isg].integralMError << ", ";
	cout << "n(est)/p " << setw(6) << setprecision(4) << nestp << "+-" << nestpe << ", ";
	cout << "t/3He " << setw(6)   << setprecision(4) << t3He << "+-" << t3Hee << endl; 

	cout << " ---" << endl;

	cout << " n "   << setw(6)   << setprecision(4) << physAMD[iisys][7][ig][isg].integralM 
	     << " p "   << setw(6)   << setprecision(4) << physAMD[iisys][0][ig][isg].integralM
	     << " t "   << setw(6)   << setprecision(4) << physAMD[iisys][2][ig][isg].integralM
	     << " 3He " << setw(6)   << setprecision(4) << physAMD[iisys][3][ig][isg].integralM;
	  
	cout << " ---" << endl;
	
	cout << " n "   << setw(6)   << setprecision(4) << numClust[7][0]
	     << " p "   << setw(6)   << setprecision(4) << numClust[0][0]
	     << " t "   << setw(6)   << setprecision(4) << numClust[2][0]
	     << " 3He " << setw(6)   << setprecision(4) << numClust[3][0]
	     << endl;

	grphnpAMD   -> SetPoint( AMDLabel.size(),   AMDLabel.size()+1, nestp  );
	grphnpAMD   -> SetPointError( AMDLabel.size(),              0, nestpe );
	grpht3heAMD -> SetPoint( AMDLabel.size(), AMDLabel.size()+1, t3He   );
	grpht3heAMD -> SetPointError( AMDLabel.size(),            0, t3Hee  );

	dr  = t3He / nestp;
	dre = dr * sqrt( pow(t3Hee/t3He,2) + pow(nestpe/nestp,2) ); 
	grphDRAMD   -> SetPoint( AMDLabel.size(),   AMDLabel.size()+1, dr );
	grphDRAMD   -> SetPointError( AMDLabel.size(),              0, dre);
	///@@@@@@@@@@@2
	AMDLabel.push_back( physAMD[iisys][0][ig][isg].version(5,20) );


	if( ivsys.size() > 2 ){

	  Int_t isys = -1;
	  for( auto jsys : ROOT::TSeqI(nAMDsys) )for( auto iamd : ROOT::TSeqI(amdsys[jsys].size()) ){
	    
	      if( ig==amdsys[jsys][iamd] ) {
		isys = jsys;
		break;
	      }
	    }

	  cout << ig << " --- " << isys << endl;
	  if( isys > -1 ) { 
	    grphnpAMDsys[isys]   -> SetPoint( grphnpAMDsys[isys]->GetN(),    sysNA[iisys], np  );
	    grphnpAMDsys[isys]   -> SetPointError(grphnpAMDsys[isys]->GetN(),           0., npe ); 
	    grpht3heAMDsys[isys] -> SetPoint( grpht3heAMDsys[isys]->GetN(), sysNA[iisys], t3He  );
	    grpht3heAMDsys[isys] -> SetPointError(grpht3heAMDsys[isys]->GetN(),           0., t3Hee ); 
	    grphDRAMDsys[isys]   -> SetPoint( grphDRAMDsys[isys]->GetN(),   sysNA[iisys], dr  );
	    grphDRAMDsys[isys]   -> SetPointError(grphDRAMDsys[isys]->GetN(),           0., dre ); 
	  }    
	}


	//	grph  ->SetFillStyle( AMDnames.at(ig).fStyle);
	grph  ->SetLineColor( AMDnames.at(ig).fColor+5);
	grph  ->SetFillColor( AMDnames.at(ig).fColor+5);

	grphN  ->SetLineColor( AMDnames.at(ig).fColor+10);
	grphN  ->SetFillColor( AMDnames.at(ig).fColor+10);

	grphA ->SetLineColor( AMDnames.at(ig).fColor);
	grphA ->SetFillColor( AMDnames.at(ig).fColor);

	//	grphZ ->SetFillStyle( AMDnames.at(ig).fStyle);
	grphZ ->SetLineColor( AMDnames.at(ig).fColor);
	grphZ ->SetFillColor( AMDnames.at(ig).fColor);

	//	grphx ->SetFillStyle( AMDnames.at(ig).fStyle);
	grphx ->SetLineColor( AMDnames.at(ig).fColor);
	grphx ->SetFillColor( AMDnames.at(ig).fColor);

	//	grphxA->SetFillStyle( AMDnames.at(ig).fStyle);
	grphxA->SetLineColor( AMDnames.at(ig).fColor);
	grphxA->SetFillColor( AMDnames.at(ig).fColor);

	grphO ->SetLineColor( AMDnames.at(ig).fColor );
	grphO ->SetFillColor( AMDnames.at(ig).fColor );

	grphOA->SetLineColor( AMDnames.at(ig).fColor );
	grphOA ->SetLineWidth( 2);
	//	grphOA->SetFillColor( AMDnames.at(ig).fColor );

	if( ivsys.size()<=2 && AMDnames.at(ig).category == icateg){
	  LOG(INFO) << "[AMD] Adding to plot (" << ig << ") "  <<AMDnames.at(ig).config  << endl;

	  // mgr->Add(grph,  "pA");
	  // mgr->Add(grphN, "pA");
	  mgr->Add(grphO, "p");
	  
	  //	  lg ->AddEntry(grphA,AMDnames.at(ig).config+"(y>0)");	  mgrA ->Add(grphA, "p");
	  lg ->AddEntry(grphO,AMDnames.at(ig).config);	  mgrA ->Add(grphOA,"pl");
	  mgrZ ->Add(grphZ, "3");
	  mgrx ->Add(grphx, "p");
	  mgrxA->Add(grphxA,"3");
	}

	for( auto i : ROOT::TSeqI(4) ) {
	  g_mH[i] -> SetPoint  (iamd,    iamd+1, multH[i]);
	  g_mH[i] -> SetPointError(iamd,    0.5, sqrt(multHe[i]));
	}
	iamd++;
      }
    }
    

    grphnpAMD0 -> SetMarkerStyle(20);
    grphnpAMD0 -> SetMarkerColor(3);
    grphnpAMD0 -> SetLineColor(3);
    mgnp -> Add(grphnpAMD0, "p");

    grphnpAMD -> SetMarkerStyle(20);
    grphnpAMD -> SetMarkerColor(2);
    grphnpAMD -> SetLineColor(2);
    mgnp -> Add(grphnpAMD, "p");

    grpht3heAMD0 -> SetMarkerStyle(20);
    grpht3heAMD0 -> SetMarkerColor(3);
    grpht3heAMD0 -> SetLineColor(3);
    mgt3he -> Add(grpht3heAMD0, "p");

    grpht3heAMD -> SetMarkerStyle(20);
    grpht3heAMD -> SetMarkerColor(2);
    grpht3heAMD -> SetLineColor(2);
    mgt3he    -> Add(grpht3heAMD, "p");

    grphDRAMD0 -> SetMarkerStyle(20);
    grphDRAMD0 -> SetMarkerColor(3);
    grphDRAMD0 -> SetLineColor(3);
    mgDR     -> Add(grphDRAMD0, "p");

    grphDRAMD -> SetMarkerStyle(20);
    grphDRAMD -> SetMarkerColor(2);
    grphDRAMD -> SetLineColor(2);
    mgDR     -> Add(grphDRAMD, "p");

    inx++;
  }

		 
  //---------> data
  inx = 0;
  for( auto iisys : ivsys ) {

    grph = new TGraphErrors();
    grph -> SetName(Form("phpara_%d",iisys) );

    grphA = new TGraphErrors();
    grphA -> SetName(Form("phpara_%d",iisys) );

    grphZ = new TGraphErrors();
    grphZ -> SetName(Form("phparz_%d",iisys) );

    grphx = new TGraphErrors();
    grphx -> SetName(Form("phparx_%d",iisys) );

    grphxA = new TGraphErrors();
    grphxA -> SetName(Form("phparxa_%d",iisys) );

    Double_t multH[npart], multHe[npart];
    for( auto i : ROOT::TSeqI(npart) ) {multH[i] = 0.; multHe[i]=0.;}
    Double_t multn = 0, multne = 0.;
    Double_t multp = 0, multpe = 0.;
    Double_t nest  = 0.,neste  = 0.;

    for( auto ipart : ROOT::TSeqI(5) ) {
      if( ipart < 3 ) {
	multH[2]  += physData[iisys][ipart][isg].integralM*ncls[ipart].Z;
	multHe[2] += pow( physData[iisys][ipart][isg].integralMError*ncls[ipart].Z ,2);
	multH[3]  += physData[iisys][ipart][isg].integralyx*ncls[ipart].Z;
	multHe[3] += pow( physData[iisys][ipart][isg].integralyxError*ncls[ipart].Z ,2);
      }

      if( ipart < 5 ){ 
	multn += physData[iisys][ipart][isg].integralM * (ncls[ipart].A - ncls[ipart].Z);
	multne+= pow( physData[iisys][ipart][isg].integralMError * (ncls[ipart].A - ncls[ipart].Z), 2 );
	multp += physData[iisys][ipart][isg].integralM * ncls[ipart].Z;
	multpe+= pow( physData[iisys][ipart][isg].integralMError * ncls[ipart].Z, 2 );
      }
    }

    slabel.clear();
    UInt_t iin = 0;
    for( auto ipart : sqpart ) {
      if( ipart == 5 || ipart == 6 ) continue;

      slabel.push_back(ncls[ipart].pLabel);

      grph  -> SetPoint     (iin,   iin+1, physData[iisys][ipart][isg].integralM );
      grph  -> SetPointError(iin,       0, physData[iisys][ipart][isg].integralMError );
      grphA -> SetPoint     (iin,   iin+1, physData[iisys][ipart][isg].integralM*ncls[ipart].A );
      grphA -> SetPointError(iin,       0, physData[iisys][ipart][isg].integralMError*ncls[ipart].A );
      grphZ -> SetPoint     (iin,   iin+1, physData[iisys][ipart][isg].integralM*ncls[ipart].Z );
      grphZ -> SetPointError(iin,       0, physData[iisys][ipart][isg].integralMError*ncls[ipart].Z );

      grphx  -> SetPoint     (iin,   iin+1, physData[iisys][ipart][isg].integralyx );
      grphx  -> SetPointError(iin,       0, physData[iisys][ipart][isg].integralyxError );
      grphxA -> SetPoint     (iin,   iin+1, physData[iisys][ipart][isg].integralyx*ncls[ipart].A );
      grphxA -> SetPointError(iin,       0, physData[iisys][ipart][isg].integralyxError*ncls[ipart].A );

      if( ipart == 7 ) {
	//	cout  << "neutron " << setw(12)  << (sysN[iisys]-(100.-multp) - multn)/ physData[iisys][0][isg].integralM << endl;
	nest  = sysN[iisys] - neutAlpha[iisys]*(100.-multp) - multn;
	neste = sqrt( pow(neutAlpha[iisys]*multp * ( pow(neutAlphae/neutAlpha[iisys],2) + pow(multpe/multp,2) ),2) + multne*multne );

	grph  -> SetPoint     (iin, iin+1, nest  );
	grph  -> SetPointError(iin,     0, neste );
	grphA -> SetPoint     (iin, iin+1, nest  );
	grphA -> SetPointError(iin,     0, neste );
      }
      iin++;

      multH[0]  += physData[iisys][ipart][isg].integralM*ncls[ipart].Z;
      multHe[0] += pow( physData[iisys][ipart][isg].integralMError*ncls[ipart].Z ,2);
      multH[1]  += physData[iisys][ipart][isg].integralyx*ncls[ipart].Z;
      multHe[1] += pow( physData[iisys][ipart][isg].integralyxError*ncls[ipart].Z ,2);

    }

    multne = sqrt( multne );
    multpe = sqrt( multpe );

    SetXTextLabel(*mgr,  sqpart);
    SetXTextLabel(*mgrA, sqpart);
    SetXTextLabel(*mgr,  sqpart);
    SetXTextLabel(*mgrZ, sqpart);


    // Double_t nest  = sysN[iisys] - neutAlpha[iisys]*(100. - multp)- multn; 
    // Double_t neste = sqrt( multne*multne + pow(1.2*multpe,2) );
    Double_t nestp = nest/physData[iisys][0][isg].integralM;
    Double_t nestpe= nestp * sqrt( pow(neste/nest,2) + pow( physData[iisys][0][isg].integralMError / physData[iisys][0][isg].integralM , 2) );
    Double_t t3He  = physData[iisys][2][isg].integralM / physData[iisys][3][isg].integralM;
    Double_t t3Hee = t3He * sqrt( pow(physData[iisys][2][isg].integralMError/physData[iisys][2][isg].integralM, 2) +
				  pow(physData[iisys][3][isg].integralMError/physData[iisys][3][isg].integralM, 2) );


    grphnpData -> SetPoint     ( 0, 0., nestp );
    grphnpData -> SetPointError( 0, 0., nestpe);
    grphnpData -> SetPoint     ( 1, AMDLabel.size()+1, nestp );
    grphnpData -> SetPointError( 1, 0,  nestpe);

    cout << " nestp " << endl;
    grphnpData -> Print();

    grpht3heData -> SetPoint     ( 0, 0., t3He );
    grpht3heData -> SetPointError( 0, 0., t3Hee);
    grpht3heData -> SetPoint     ( 1, AMDLabel.size()+1, t3He );
    grpht3heData -> SetPointError( 1, 0,  t3Hee);
    
    cout << " t3he " << endl;
    grpht3heData -> Print();

    Double_t dr  = t3He / nestp;
    Double_t dre = dr * sqrt( pow(t3Hee/t3He,2) + pow(nestpe/nestp,2) ); 
    grphDRData   -> SetPoint( 0,      0, dr );
    grphDRData   -> SetPointError( 0, 0, dre);
    grphDRData   -> SetPoint( 1,      AMDLabel.size()+1, dr );
    grphDRData   -> SetPointError( 1, 0, dre);

    cout << " dr " << endl;
    grphDRData->Print(); 


    if( ivsys.size() >= 2 ) {
      grphnpDatasys   -> SetPoint( inx, sysNA[iisys], nestp  );
      grphnpDatasys   -> SetPointError(inx,       0., nestpe ); 
      grpht3heDatasys -> SetPoint( inx, sysNA[iisys], t3He  );
      grpht3heDatasys -> SetPointError(inx,       0., t3Hee ); 
      grphDRDatasys   -> SetPoint( inx, sysNA[iisys], dr  );
      grphDRDatasys   -> SetPointError(inx,       0., dre ); 
    }

    

    grphnpData -> SetFillStyle(3001);
    grphnpData -> SetFillColor( kGreen );
    grpht3heData -> SetFillStyle(3001);
    grpht3heData -> SetFillColor( kGreen );
    grphDRData -> SetFillStyle(3001);
    grphDRData -> SetFillColor( kGreen );
    mgnp   -> Add( grphnpData,   "3" );
    mgt3he -> Add( grpht3heData, "3" );
    mgDR   -> Add( grphDRData,   "3" );
    

    cout << "DATA : " << Asys[iisys] << " " << sysN[iisys] << " " ;
    cout << setw(6)  << setprecision(4) << physData[iisys][0][isg].integralM << "+- " << physData[iisys][0][isg].integralMError << " , ";
    cout << setw(6)  << setprecision(4) << physData[iisys][1][isg].integralM << "+- " << physData[iisys][1][isg].integralMError << " , ";
    cout << setw(6)  << setprecision(4) << physData[iisys][2][isg].integralM << "+- " << physData[iisys][2][isg].integralMError << " , ";
    cout << setw(6)  << setprecision(4) << physData[iisys][3][isg].integralM << "+- " << physData[iisys][3][isg].integralMError << " , ";
    cout << setw(6)  << setprecision(4) << physData[iisys][4][isg].integralM << "+- " << physData[iisys][4][isg].integralMError << " , ";
    cout << " n(est)  " << setw(6) << setprecision(4) << nest  << " +- " << neste << " , ";
    cout << " n(est)/p "<< setw(6) << setprecision(4) << nestp << " +- " << nestpe << " , ";
    cout << " t/3He   " << setw(6) << setprecision(4) << t3He  << " +- " << t3Hee << endl;
      

    // cout << setw(6) << setprecision(4) << physData[iisys][0][isg].integralM
    // 	 << setw(6) << " " 
    // 	 << setw(6) << setprecision(4) << physData[iisys][2][isg].integralM
    // 	 << setw(6) << setprecision(4) << physData[iisys][3][isg].integralM
    // 	 << endl;

    for( auto i : ROOT::TSeqI(4)) {
      aBox[iisys][i] = new TBox( g_mH[i]->GetXaxis()->GetXmin(), multH[i]-multHe[i]*0.5,
				 g_mH[i]->GetXaxis()->GetXmax(), multH[i]+multHe[i]*0.5);
      aBox[iisys][i] -> SetFillColorAlpha( CStyle[iisys].fColor, 0.6 );
    }

    grph->SetMarkerStyle( DStyle[iisys].mStyle );
    grph->SetMarkerSize(  DStyle[iisys].mSize);
    grph->SetLineColor(   DStyle[iisys].fColor );
    grph->SetMarkerColor( DStyle[iisys].fColor );
    mgr->Add(grph,"p");


    grphA->SetMarkerStyle( DStyle[iisys].mStyle );
    grphA->SetMarkerSize(  DStyle[iisys].mSize+0.5 );
    grphA->SetLineColor(   DStyle[iisys].fColor );
    grphA->SetMarkerColor( DStyle[iisys].fColor );
    mgrA->Add(grphA,"p");

    grphZ->SetMarkerStyle( DStyle[iisys].mStyle );
    grphZ->SetMarkerSize(  DStyle[iisys].mSize );
    grphZ->SetLineColor(   DStyle[iisys].fColor );
    grphZ->SetMarkerColor( DStyle[iisys].fColor );
    mgrZ->Add(grphZ,"p");

    grphx->SetMarkerStyle( DStyle[iisys].mStyle );
    grphx->SetMarkerSize(DStyle[iisys].mSize);
    grphx->SetLineColor( DStyle[iisys].fColor );
    grphx->SetMarkerColor( DStyle[iisys].fColor );
    mgrx->Add(grphx,"p");

    grphxA->SetMarkerStyle( DStyle[iisys].mStyle );
    grphxA->SetMarkerSize(  DStyle[iisys].mSize );
    grphxA->SetLineColor(   DStyle[iisys].fColor );
    grphxA->SetMarkerColor( DStyle[iisys].fColor );
    mgrxA ->Add(grphxA, "p");

    ///    lg ->AddEntry(grph,fsys[iisys]+"_"+physData[iisys][0][isg].version); 
    lg ->AddEntry(grph,"DATA"); 
  

    Double_t yrange[4][2] = {{2.3, 3.2},{1.3, 2.2},{1.6, 2.4},{1.6, 2.4}};
    mgnp   -> GetYaxis()->SetRangeUser(yrange[iisys][0], yrange[iisys][1]);
    mgt3he -> GetYaxis()->SetRangeUser(yrange[iisys][0], yrange[iisys][1]);
    mgnp   -> GetXaxis()->SetLimits(0.5,AMDLabel.size()+0.5);
    mgt3he -> GetXaxis()->SetLimits(0.5,AMDLabel.size()+0.5);
    
    inx++;
  }

  auto lgsys = new TLegend(0.18,0.7, 0.6,0.9,"");
  if( ivsys.size() >= 2 ) {

    grphnpDatasys  -> SetFillColor(kBlack);
    grphnpDatasys  -> SetFillStyle(3001);
    grphnpDatasys  -> SetLineColor(kBlack);
    grphnpDatasys  -> SetMarkerColor(kBlack);
    grphnpDatasys  -> SetMarkerStyle(20);
    grphnpDatasys  -> SetMarkerSize(1.5);

    grpht3heDatasys  -> SetFillColor(kBlack);
    grpht3heDatasys  -> SetFillStyle(3001);
    grpht3heDatasys  -> SetLineColor(kBlack);
    grpht3heDatasys  -> SetMarkerColor(kBlack);
    grpht3heDatasys  -> SetMarkerStyle(20);
    grpht3heDatasys  -> SetMarkerSize(1.5);

    grphDRDatasys  -> SetFillColor(kBlack);
    grphDRDatasys  -> SetFillStyle(3001);
    grphDRDatasys  -> SetLineColor(kBlack);
    grphDRDatasys  -> SetMarkerColor(kBlack);
    grphDRDatasys  -> SetMarkerStyle(20);
    grphDRDatasys  -> SetMarkerSize(1.5);

    cout << " np " << endl;
    grphnpDatasys->Print();
    cout << " t3he " << endl;
    grpht3heDatasys -> Print();
    cout << " DR " << endl;
    grphDRDatasys -> Print();
    cout << endl;
    mgnpsys   -> Add(grphnpDatasys,  "p3");
    mgt3hesys -> Add(grpht3heDatasys,"p3");
    mgDRsys   -> Add(grphDRDatasys,  "p3");

    lgsys -> AddEntry(grphnpDatasys,"DATA");

    for( auto isys : ROOT::TSeqI(nAMDsys) ) {
      Color_t fcol = AMDnames.at(amdsys[isys][0]).fColor;
      grphnpAMDsys[isys]  -> SetFillColor(fcol);
      grphnpAMDsys[isys]  -> SetLineColor(fcol);
      grphnpAMDsys[isys]  -> SetMarkerColor(fcol);
      grphnpAMDsys[isys]  -> SetMarkerStyle(MarkerStyle[isys]);
      grphnpAMDsys[isys]  -> SetMarkerSize(1.8);

      grpht3heAMDsys[isys]  -> SetFillColor(fcol);
      grpht3heAMDsys[isys]  -> SetLineColor(fcol);
      grpht3heAMDsys[isys]  -> SetMarkerColor(fcol);
      grpht3heAMDsys[isys]  -> SetMarkerStyle(MarkerStyle[isys]);
      grpht3heAMDsys[isys]  -> SetMarkerSize(1.8);

      grphDRAMDsys[isys]   -> SetFillColor(fcol);
      grphDRAMDsys[isys]   -> SetLineColor(fcol);
      grphDRAMDsys[isys]   -> SetMarkerColor(fcol);
      grphDRAMDsys[isys]   -> SetMarkerStyle(MarkerStyle[isys]);
      grphDRAMDsys[isys]   -> SetMarkerSize(1.8);

      mgnpsys   -> Add(grphnpAMDsys[isys], "p");
      mgt3hesys -> Add(grpht3heAMDsys[isys],"p");
      mgDRsys   -> Add(grphDRAMDsys[isys], "p");
      lgsys -> AddEntry(grphnpAMDsys[isys], AMDnames.at( amdsys[isys][0] ).fullconfig, "lep");

    }
  }

  //++++ plot
  for( auto i : ROOT::TSeqI(AMDLabel.size()) ) { 
    mgnp  ->GetXaxis()->SetBinLabel( mgnp  ->GetXaxis()->FindBin( i+1 ), AMDLabel.at(i) );
    mgt3he->GetXaxis()->SetBinLabel( mgt3he->GetXaxis()->FindBin( i+1 ), AMDLabel.at(i) );
    mgDR  ->GetXaxis()->SetBinLabel( mgDR  ->GetXaxis()->FindBin( i+1 ), AMDLabel.at(i) );
  }     

  TString ssys = "";
  for( auto iisys : ivsys ) {
    ssys += rsys[iisys]+"Sn";
  }
  ssys += lbCentral[bCentral];


  if( 1 ) {
    id = 1;
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),700,500); iccv++; 
    ccv ->  Divide(2,1);

    ccv -> cd(id); id++;
    lg -> Draw();

    ccv->cd(id); id++;
    SetXTextLabel(*mgrA, sqpart);
    mgrA->GetYaxis()->SetRangeUser(0.,80.);
    mgrA->Draw("ALP");

    ccv->SaveAs(Form("cmpItg0_%s.png",ssys.Data()));
  }

  if( 0 ) {
    id = 1;
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1000,500); iccv++; 
    ccv ->  Divide(3,1);

    ccv -> cd(id); id++;
    lg -> Draw();

    ccv->cd(id); id++;
    mgr->GetYaxis()->SetRangeUser(0.,80.);
    mgr->Draw("ALP");

    ccv->cd(id); id++;
    SetXTextLabel(*mgrA, sqpart);
    mgrA->GetYaxis()->SetRangeUser(0.,80.);
    mgrA->Draw("ALP");

    ccv->SaveAs(Form("cmpItg2_%s.png",ssys.Data()));
  }


  if( ivsys.size() >= 2 ) {
    id = 1;
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1000,500); iccv++; 
    ccv ->  Divide(3,1);
    auto pad = ccv -> cd(id); id++;
    mgnpsys -> SetTitle(";System(n/p); R(n/p)");
    mgnpsys -> GetXaxis() -> SetLimits(1.0, 1.68);
    // mgnpsys -> GetXaxis() -> SetLimits(1.0, 2.5);
    mgnpsys -> GetYaxis() -> SetRangeUser(1., 3.2);
    mgnpsys -> Draw("AP");
    mgnpsys -> Print();
    lgsys   -> Draw();

    pad = ccv -> cd(id); id++;
    mgt3hesys -> SetTitle(";System(n/p); R(t/^{3}He)");
    mgt3hesys -> GetXaxis() -> SetLimits(1.0, 1.68);
    //mgt3hesys -> GetXaxis() -> SetLimits(1.0, 2.5);
    mgt3hesys -> GetYaxis() -> SetRangeUser(1., 3.2);
    mgt3hesys -> Draw("AP");

    pad = ccv -> cd(id); id++;
    mgDRsys -> SetTitle(";System(n/p); DR(R(t/^{3}He) / R(n/p))");
    //    mgDRsys -> GetXaxis() -> SetLimits(1.0, 2.5);
    mgDRsys -> GetXaxis() -> SetLimits(1.0, 1.68);
    mgDRsys -> GetYaxis() -> SetRangeUser(0.7, 1.4);
    mgDRsys -> Draw("AP");
    TLine line1(1., 1., 1.68, 1.);
    line1.SetLineStyle(2);
    line1.Draw();

    ccv->SaveAs(Form("cmpnpDR_%s.png",ssys.Data()));
  }


  if( ivsys.size() == 2 ) {
    id = 1;
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1000,500); iccv++; 
    ccv ->  Divide(3,1);

    auto pad = ccv -> cd(id); id++;
    pad -> SetBottomMargin(5.);
    pad -> Draw();
    mgnp -> Draw("AP");

    pad = ccv -> cd(id); id++;
    pad -> SetBottomMargin(5.);
    pad -> Draw();
    mgt3he -> Draw("AP");

    pad = ccv -> cd(id); id++;
    pad -> SetBottomMargin(5.);
    pad -> Draw();
    mgDR -> Draw("AP");
    ccv->SaveAs(Form("cmpnp_%s.png",ssys.Data()));
  }

  if( 0 ) {
    id = 1;
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1400,800); iccv++; 
    ccv->Divide(4,2);
    ccv->cd(id); id++;
    lg ->Draw();

    //  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++; 
    ccv->cd(id); id++;
    mgr -> Draw("ALP");

    //  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++; 
    ccv->cd(id); id++;
    mgrA->Draw("A");

    ccv->cd(id); id++;
    mgrZ->Draw("A");

    UInt_t imt = 0;
    ccv->cd(id); id++;
    g_mH[imt]->GetYaxis()->SetRangeUser(50., 100.);
    g_mH[imt]->Draw("A*");
    for( auto iisys : ivsys )    
      aBox[iisys][imt]->Draw();

    //  ccv->cd(id); id++;
    //  for( auto iisys : ivsys ) {
    //   LOG(INFO) << " Total M " << rsys[iisys] << " ->  " << multH[iisys] << " +- " << multHe[iisys] 
    // 	      << " ->  " << multxH[iisys] << " +- " << multxHe[iisys] << FairLogger::endl;
    
    //   plabel.DrawLatexNDC(0.2,0.9-0.1*iisys,rsys[iisys]+": "+Form("M %8.3f +- %6.3f",multH[iisys],multHe[iisys] ));    
    //   plabel.DrawLatexNDC(0.2,0.4-0.1*iisys,rsys[iisys]+Form(" |y|<0.3 : %8.3f +- %6.3f",multxH[iisys],multxHe[iisys] ));    
    //  }

    ccv->cd(id); id++;
    mgrx->SetTitle("|y|<0.3;;M");
    SetXTextLabel(*mgrx, sqpart);
    mgrx->Draw("A PLC");

    ccv->cd(id); id++;
    mgrx->SetTitle("|y|<0.3;;A*M");
    SetXTextLabel(*mgrxA, sqpart);
    mgrxA->Draw("A");

    TString ssys = "";
    for( auto iisys : ivsys ) {
      ssys += rsys[iisys]+"Sn";
    }
    ssys += lbCentral[bCentral];
    ccv->SaveAs(Form("cmpItg_%s.png",ssys.Data()));
  }
}
//----------------------
// paper plots
//-----

void Draw_PaperHistgram(std::vector<gplot> gname,  Int_t ylim=-1.,  Bool_t bAMD=1, UInt_t hplt=0, Int_t icateg=1) 
{
  LOG(INFO) << "[Draw_PaperHistgram] .......... Centrality " << bCentral ;

  std::vector<UInt_t> ivsys  = {0, 3, 1};
  std::vector<UInt_t> sqpart = {0,1,2,3,4};;

  std::vector< std::pair< UInt_t, TString >> phist = {{0, "h_dndy"},{1,"h_dndEt"},{2,"h_dndbtgm"},{3,"h_dndyx"},{4,"h_dndpx"},{5,"h_dndpt"},
						      {6,"hypx"}   ,{7,"gy_v1"}  ,{8,"gy_v2"}};
  TString ylabel[] = {"A*dN/dy","A*dN/Et","A*dN/d(#beta#gamma_{T})","A*dN/dy_{x}","dN/dP_{x}","dN/dP_{t}",
		      "<px>/A","v1","v2"};
  TString xlabel[] = {"y/y^{cm}_{NN}-1","E_{t}","#beta#gamma_{T}","y_{x}","P_{x}","P_{t}",
		      "y/y^{cm}_{NN}-1","y/y^{cm}_{NN}-1","y/y^{cm}_{NN}-1"};

  Double_t umax[][2] = {{-2.,35.},{0, 0.4},{0.,0.4},{0,35.0},{0.,0.},{0.,0.},{-100.,100.},{-0.8,0.8},{-0.07,0.07}};

  TString Ylabel, Xlabel;
  
  const UInt_t nsys  = (UInt_t)ivsys.size();
  const UInt_t npart = (UInt_t)sqpart.size();
  const UInt_t ndata = (UInt_t)gname.size();

  THStack*     mplt[nsys][npart];
  TLegend*     lg[nsys];

  TString histname = phist[hplt].second;
  LOG(INFO) << phist[hplt].second  << " -> " << histname  <<endl; 
  for( auto ksys : ROOT::TSeqI(nsys) ) {
    lg[ksys]  = new TLegend(0.2, 0.55, 0.9, 0.9,fsys[ ivsys[ksys] ] );

    for( auto ipart : sqpart )  {
      mplt[ksys][ipart] = new THStack(Form(histname+"_%d",ipart),"");
    }
  }



  TH1D*         hPlot = NULL;

  for( auto ic : {0,1} ) { // midcentral andl central
    
    if( bCentral != ic  ) continue;

    UInt_t ksys = 0;
    for( auto iisys : ivsys){

      for( auto ipart : sqpart ) {
	  
	for( auto igname : gname ) {
	    
	  if( igname.centrality != ic ) continue;

	  hPlot = (TH1D*)LoadData(iisys,ipart,histname,igname,ylim); 

	  LOG(INFO) << "[LoadDATA] opening... "<< igname.config1 << " " << ipart << " " << histname << " ->  ";

	  if( hPlot != NULL ) {
	    TString nname = (TString)hPlot->GetName()+Form("%d",ic);
	    hPlot -> TH1D::SetName(nname);

	    hPlot -> SetMarkerStyle(DStyle[iisys].mStyle);
	    hPlot -> SetMarkerSize( DStyle[iisys].mSize-0.6);
	    hPlot -> SetLineColor(  DStyle[iisys].fColor);//+ic*5);
	    hPlot -> SetMarkerColor(DStyle[iisys].fColor);//+ic*5);
	      
	    if( ipart == 0 )
	      lg[ksys]  -> AddEntry( hPlot, "DATA");
	      
	    if( hplt < 6 ) 
	      hPlot -> Scale( ncls[ipart].A );

	    mplt[ksys][ipart] -> Add( hPlot );
	      
	  }
	}
	
    
	//@@@@ AMD
	if( bAMD ) {    

	  for( auto samd : AMDnames ) {

	    if( samd.category != icateg ) continue;

	    hPlot = (TH1D*)LoadAMD(iisys, ipart, samd, histname, ic, ylim);

	    LOG(INFO) <<" [AMD] " << samd.config << " " << histname << " iisys " << iisys <<  FairLogger::endl;
	      
	    if( hPlot != NULL ) {
		
	      UInt_t iin = samd.id;

	      hPlot -> SetName((TString)hPlot->GetName()+Form("_%d",iin));
	      LOG(INFO) <<" [AMD] open " <<  hPlot->GetName() <<  FairLogger::endl;

	      hPlot -> Scale( ncls[ipart].A );
		
	      hPlot -> SetFillStyle( samd.fStyle );
	      hPlot -> SetLineColor( samd.fColor );
	      hPlot -> SetFillColor( samd.fColor );
	      // hPlot -> SetFillStyle( 3001 );
	      // hPlot -> SetLineWidth(-5002 );
		
	      mplt[ksys][ipart] -> Add( hPlot,"aE4");
		  
	      if( ipart == 0 && ksys == 0 ) 
		lg[ksys]  -> AddEntry( hPlot, "AMD"+ samd.config);

	    }
	  }
	}
      }

      ksys++;
    }
  }


  //--> plot all histogram
  const Int_t Ny = nsys;
  const Int_t Nx = npart;
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),300*Nx,300*Ny); iccv++;
  //    ccv->Divide(Nx,Ny);
  Float_t lMargin = 0.04;
  Float_t rMargin = 0.01;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.02;
  Float_t midMargin = 0.005;

  CanvasPartitionTwoColumn( ccv, Nx, Ny, lMargin, rMargin, bMargin, tMargin, midMargin );

  //  mhsty[ipart][2] -> GetXaxis()->SetLimits(-800.,800.);
    
  id = 1;

  TPad* pad0 = (TPad*)gROOT->FindObject("pad_0_0");

  for( auto ipady : ROOT::TSeqI(Ny) ) {

    std::vector< UInt_t >::iterator ipx = sqpart.begin(); 
    for( auto ipadx : ROOT::TSeqI(Nx) ) {
      //	auto pad = ccv -> cd(id); id++;
      auto pad = (TPad*)gROOT->FindObject(Form("pad_%d_%d",ipadx,Ny-ipady-1));
      if( !pad ) return;

      pad->cd();
      pad->SetGrid();

      Float_t xFactor = pad0->GetAbsWNDC()/pad->GetAbsWNDC();
      Float_t yFactor = pad0->GetAbsHNDC()/pad->GetAbsHNDC();

      cout << " xFactor " << xFactor << " yFactor " << yFactor << endl;

      // pad->SetLogy();
      // pad->SetLeftMargin  (lMargin);
      // pad->SetRightMargin (rMargin);
      // pad->SetTopMargin   (tMargin);
      // pad->SetBottomMargin(bMargin);
      
      if( phist.at( hplt ).second == "h_dndpt") // ||
	pad->SetLogy();



      if( ipadx < npart ){
	mplt[ipady][*ipx] -> SetTitle();

	mplt[ipady][*ipx] -> Draw("nostack");      

	if( mplt[ipady][*ipx] -> GetXaxis() )
	  mplt[ipady][*ipx] -> GetXaxis()->SetNdivisions(510);

	if( umax[hplt][0] != 0 ) {
	  mplt[ipady][*ipx] -> SetMaximum( umax[hplt][1] );	    
	  
	  if( umax[hplt][0] < 0 ) {
	    mplt[ipady][*ipx] -> GetXaxis() -> SetRangeUser(umax[hplt][0], abs(umax[hplt][0]));
	  }
	  else
	    mplt[ipady][*ipx] -> GetXaxis() -> SetRangeUser(0., abs(umax[hplt][0]));
	}

	if( ipady == 0 ) {
	  plabel.SetTextAlign(15);
	  plabel.SetTextSize(0.08);
	  plabel.DrawLatexNDC(0.8,0.85, ncls[sqpart.at(ipadx)].sName);
	}

	if( ipadx == 0 )
	  mplt[ipady][*ipx] -> GetYaxis()->SetTitle(ylabel[hplt]);

	//	  if( ipady == 1 ) 
	mplt[ipady][*ipx] -> GetXaxis()->SetTitle(xlabel[hplt]);

      }
      
      if( ipadx == 0 ){
	lg[ipady]->SetFillStyle(0);
	lg[ipady]->SetTextSize(0.06*yFactor);	
	lg[ipady]->Draw();
      }

      ipx++;
    }
  }

  TString ssys = "";
  for( auto iisys : ivsys ) {
    ssys += rsys[iisys]+"Sn";
  }
  ssys += lbCentral[bCentral];
  ccv->SaveAs(Form("cmp%s_%s.png",phist[hplt].second.Data(),ssys.Data()));  

  if( 0 ) {
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1200,800); iccv++;
    ccv -> Divide(3,2);
    id = 1;
    for( auto i : {5,6} )for( auto j: {0,1,2} ) {
	ccv -> cd(id); id++;
	mplt[i][j] -> Draw("nostack");
      }
  }

}
void Draw_PaperGraph(std::vector<gplot> gname,  Int_t ylim=-1.,  Bool_t bAMD=1, UInt_t hplt=0, Int_t icateg=1) 
{
  LOG(INFO) << "[Draw_PaperGraph] .......... Centrality " << bCentral ;

  std::vector<UInt_t> ivsys  = {1, 0};
  std::vector<UInt_t> sqpart = {0,1,2,3,4};;

  std::vector< std::pair< UInt_t, TString >> phist = {{0, "h_dndy"},{1,"h_dndEt"},{2,"h_dndbtgm"},{3,"h_dndyx"},{4,"h_dndpx"},{5,"h_dndpt"},
						      {6,"hypx"}   ,{7,"gyv1"}  ,{8,"gyv2"}};
  TString ylabel[] = {"A*dN/dy","A*dN/Et","A*dN/d(#beta#gamma_{T})","A*dN/dy_{x}","dN/dP_{x}","dN/dP_{t}",
		      "<px>/A","v1","v2"};
  TString xlabel[] = {"y_{cm}/y_{proj}","E_{t}","#beta#gamma_{T}","y_{x}","P_{x}","P_{t}",
		      "y_{cm}/y_{proj}","y_{cm}/y_{proj}","y_{cm}/y_{proj}"};

  Double_t umax[][2] = {{-2.,35.},{0, 0.4},{0.,0.4},{0,35.0},{0.,0.},{0.,0.},{-100.,100.},{-0.8,0.8},{-0.08,0.08}};

  TString Ylabel, Xlabel;
  
  const UInt_t nsys  = (UInt_t)ivsys.size();
  const UInt_t npart = (UInt_t)sqpart.size();
  const UInt_t ndata = (UInt_t)gname.size();

  TMultiGraph*     mplt[nsys][npart];
  TLegend*         lg[nsys];
  TString          sysName[nsys];
  auto lgAMD = new TLegend();

  TString histname = phist[hplt].second;
  LOG(INFO) << phist[hplt].second  << " -> " << histname  <<endl; 
  for( auto ksys : ROOT::TSeqI(nsys) ) {
    lg[ksys]  = new TLegend();

    for( auto ipart : sqpart )  {
      mplt[ksys][ipart] = new TMultiGraph();
    }
  }



  TGraphErrors*     hPlot = NULL;

  for( auto ic : {0,1} ) { // midcentral andl central
    
    if( bCentral != ic  ) continue;

    UInt_t ksys = 0;
    for( auto iisys : ivsys){

      for( auto ipart : sqpart ) {
	  
	
    
	//@@@@ AMD
	if( bAMD ) {    

	  for( auto samd : AMDnames ) {

	    if( samd.category != icateg ) continue;

	    TString ahistname = histname;
	    if( histname == "hypx" ) ahistname = "gy_px";
	    else if( histname == "gu_v1" || histname == "gy_v1" ) ahistname = "g_v1y";
	    else if( histname == "gu_v2" || histname == "gy_v2" ) ahistname = "g_v2y";
	    hPlot = (TGraphErrors*)LoadAMD(iisys, ipart, samd, ahistname, ic, ylim);

	    LOG(INFO) <<" [AMD] " << samd.config << " " << histname << " iisys " << iisys <<  FairLogger::endl;
	      
	    if( hPlot != NULL ) {
		
	      UInt_t iin = samd.id;

	      hPlot -> SetName((TString)hPlot->GetName()+Form("_%d",iin));
	      LOG(INFO) <<" [AMD] open " <<  hPlot->GetName() <<  FairLogger::endl;

	      hPlot -> SetFillStyle( samd.fStyle );
	      hPlot -> SetLineColor( samd.fColor );
	      hPlot -> SetFillColor( samd.fColor );
	      // hPlot -> SetFillStyle( 3001 );
	      // hPlot -> SetLineWidth(-5002 );

	      
	      mplt[ksys][ipart] -> Add( hPlot,"3");
		  
	      if( ipart == 0 && ksys == 0 ) 
		lgAMD  -> AddEntry( hPlot, "AMD"+ samd.config);

	    }
	  }
	}

	// DATA
	for( auto igname : gname ) {
	    
	  if( igname.centrality != ic ) continue;

	  hPlot = (TGraphErrors*)LoadData(iisys,ipart,histname,igname,ylim); 

	  LOG(INFO) << "[LoadDATA] opening... "<< igname.config1 << " " << ipart << " " << histname << " ->  ";

	  if( hPlot != NULL ) {
	    TString nname = (TString)hPlot->GetName()+Form("%d",ic);
	    hPlot -> SetName(nname);

	    hPlot -> SetMarkerStyle(DStyle[iisys].mStyle);
	    hPlot -> SetMarkerSize( 0.8);
	    hPlot -> SetLineColor(  DStyle[iisys].fColor);//+ic*5);
	    hPlot -> SetMarkerColor(DStyle[iisys].fColor);//+ic*5);

	    if( hplt == 8 ) {
	      auto fv2fit = (TF1*)gROOT->FindObject("fv2fit");
	      if( fv2fit ) {
		gStyle->SetOptFit(0);
		fv2fit->SetLineColor(DStyle[iisys].fColor );
	      }
	      Getv2FitParameters(hPlot, -1., 1.);
	    }
	    
	      
	    if( ipart == 0 ) {
	      lg[ksys]  -> AddEntry( hPlot, "DATA" );//+fsys[iisys]);
	      sysName[ksys] = fsys[iisys];

	    }
	    mplt[ksys][ipart] -> Add( hPlot, "P" );
	      
	  }
	}

	//--------
      }

      ksys++;
    }
  }


  mplt[0][0] -> SetTitle(";;"+ylabel[hplt]);
  mplt[1][0] -> SetTitle(";;"+ylabel[hplt]);
  mplt[0][2] -> SetTitle(";"+xlabel[hplt]+";");


  //--> plot all histogram
  const Int_t Ny = nsys;
  const Int_t Nx = npart;
  //  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),100*Nx,200*Ny); iccv++;
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),150*Nx,300*Ny); iccv++;
  //    ccv->Divide(Nx,Ny);
  Float_t lMargin = 0.08;
  Float_t rMargin = 0.01;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.02;
  Float_t midMargin = 0.0;

  CanvasPartitionTwoColumn( ccv, Nx, Ny, lMargin, rMargin, bMargin, tMargin, midMargin );

  //  mhsty[ipart][2] -> GetXaxis()->SetLimits(-800.,800.);
    
  id = 1;

  TPad* pad[Nx][Ny];

  for( auto ipady : ROOT::TSeqI(Ny) ) {

    std::vector< UInt_t >::iterator ipx = sqpart.begin(); 
    for( auto ipadx : ROOT::TSeqI(Nx) ) {
      //	auto pad = ccv -> cd(id); id++;
      TString pname = Form("pad_%i_%i",ipadx,ipady); 
      pad[ipadx][ipady] = (TPad*)gROOT->FindObject(pname);
      if( !pad[ipadx][ipady] ) return;

      pad[ipadx][ipady]->SetGrid();

      pad[ipadx][ipady] -> SetFillStyle(4000);
      pad[ipadx][ipady] -> SetFrameFillStyle(4000);
      pad[ipadx][ipady] -> cd();
      Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[ipadx][ipady]->GetAbsWNDC();
      Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[ipadx][ipady]->GetAbsHNDC();
      
      if( phist.at( hplt ).second == "h_dndpt") // ||
	pad[ipadx][ipady]->SetLogy();

      if( ipadx < npart ){

	//format for x axis
	mplt[ipady][*ipx] ->GetXaxis()->SetLabelFont(43);
	mplt[ipady][*ipx] ->GetXaxis()->SetLabelSize(18);
	mplt[ipady][*ipx] ->GetXaxis()->SetLabelOffset(0.02);
	mplt[ipady][*ipx] ->GetXaxis()->SetTitleFont(43);
	mplt[ipady][*ipx] ->GetXaxis()->SetTitleSize(20);
	mplt[ipady][*ipx] ->GetXaxis()->SetTitleOffset(3);
	mplt[ipady][*ipx] ->GetXaxis()->CenterTitle();
	mplt[ipady][*ipx] ->GetXaxis()->SetNdivisions(505);
	mplt[ipady][*ipx] ->GetXaxis()->SetTickLength(yFactor*0.05/xFactor);

	//format for y axis
	mplt[ipady][*ipx] ->GetYaxis()->SetLabelFont(43);
	mplt[ipady][*ipx] ->GetYaxis()->SetLabelSize(18);
	mplt[ipady][*ipx] ->GetYaxis()->SetLabelOffset(0.02);
	mplt[ipady][*ipx] ->GetYaxis()->SetTitleFont(43);
	mplt[ipady][*ipx] ->GetYaxis()->SetTitleSize(20);
	mplt[ipady][*ipx] ->GetYaxis()->SetTitleOffset(3.5);
	mplt[ipady][*ipx] ->GetYaxis()->CenterTitle();
	mplt[ipady][*ipx] ->GetYaxis()->SetNdivisions(505);

	mplt[ipady][*ipx] ->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);

	mplt[ipady][*ipx] -> GetXaxis() -> SetRangeUser(-0.9,1.3);
	mplt[ipady][*ipx] -> GetYaxis() -> SetRangeUser(umax[hplt][0], umax[hplt][1]);
	mplt[ipady][*ipx] -> GetXaxis() -> SetNdivisions(5,5,0,kTRUE);
	mplt[ipady][*ipx] -> GetYaxis() -> SetNdivisions(5,5,0,kTRUE);

	mplt[ipady][*ipx] -> Draw("alp");      

	if( ipady == Ny-1 ) {
	  plabel.SetTextAlign(15);
	  plabel.SetTextSize(0.1*xFactor);
	  plabel.DrawLatexNDC(XtoPad(0.1),YtoPad(0.85), ncls[sqpart.at(ipadx)].pLabel);
	}

	if( ipadx == 0 && ipady == 0) {
	  lgAMD -> SetFillStyle(0);
	  lgAMD -> SetTextFont(43);
	  lgAMD -> SetTextSize(18);	
	  lgAMD -> SetX1(0.28);
	  lgAMD -> SetX2(0.9);
	  lgAMD -> SetY1(YtoPad(0.05));
	  lgAMD -> SetY2(YtoPad(0.18));
	  lgAMD -> Draw();
	  lgAMD -> Draw();
	}

	if( ipadx == 0 ){
          plabel.SetTextAlign(15);
          plabel.SetTextSize(0.1*xFactor);
	  plabel.DrawLatexNDC(XtoPad(0.1),YtoPad(0.7),sysName[ipady]);
	  
	  lg[ipady]->SetFillStyle(0);
	  lg[ipady]->SetTextFont(43);
	  lg[ipady]->SetTextSize(18);	
	  lg[ipady]->SetX1(0.28);
	  lg[ipady]->SetX2(0.9);
	  lg[ipady]->SetY1(YtoPad(0.2));
	  lg[ipady]->SetY2(YtoPad(0.25));
	  lg[ipady]->Draw();
	}
      }

      ipx++;
    }
  }

  TString ssys = "";
  for( auto iisys : ivsys ) {
    ssys += rsys[iisys]+"Sn";
  }
  ssys += lbCentral[bCentral];
  ccv->SaveAs(Form("cmp%s_%s.png",phist[hplt].second.Data(),ssys.Data()));  
  ccv->SaveAs(Form("cmp%s_%s.eps",phist[hplt].second.Data(),ssys.Data()));  





  if( 0 ) {
    ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv),1200,800); iccv++;
    ccv -> Divide(3,2);
    id = 1;
    for( auto i : {5,6} )for( auto j: {0,1,2} ) {
	ccv -> cd(id); id++;
	mplt[i][j] -> Draw("nostack");
      }
  }

}

void Draw_PaperPlots(std::vector<gplot> gname,  Int_t ylim=-1.,  Bool_t bAMD=1, UInt_t hplt=0, Int_t icateg=1) 
{
  if( hplt < 6 )
    Draw_PaperHistgram( gname, ylim,  bAMD, hplt, icateg);  
  else {
    bCentral = 0;
    Draw_PaperGraph(gname,  ylim, bAMD, hplt, icateg);  
    bCentral = 1;
  }
}

 //----------------------------------------------------------------------

void PlotFigure(Bool_t nplot=kTRUE)
{
  if( !nplot )
    return 0;

  beamAToSys[132] = 0;
  beamAToSys[108] = 1;
  beamAToSys[124] = 2;
  beamAToSys[112] = 3;


  const UInt_t npart = (UInt_t)sqPart.size();
  fileDATA.resize(5);
  for( auto i : ROOT::TSeqI(5) ) {
    fileDATA[i].resize(npart);
    for( auto j : ROOT::TSeqI(npart) ) {
      fileDATA[i][j].resize(nDATA);
    }
  }

  fileAMD.resize(nAMD);
  for( auto i : ROOT::TSeqI(nAMD) )
    fileAMD[i].resize(5);

  if( nDATA == 0 ) {
    LOG(ERROR) << " No data is selected. " << FairLogger::endl;
    exit(0);
  }


  gStyle->SetLegendFillColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatStyle(0);
  gStyle->SetOptFit(1111);
  gStyle->SetStatColor(10);

  SetStyle();
  SetColor();

  //##main#
  //  std::vector< UInt_t > sqsys = {0,1};
  //std::vector< UInt_t > sqsys = {0,3,1};
  //std::vector< UInt_t > sqsys = {0};
  std::vector< UInt_t > sqsys = {1};
  //std::vector< UInt_t > sqsys = {3};
  bCentral = 1;
  //  std::vector<gplot>    gname = {gnames[bCentral]};
  std::vector<gplot>    gname =  gnames;
  TString sVersion    = gname[0].Version;
  //  std::vector< UInt_t > sqpart = {0};
  std::vector< UInt_t > sqpart = {0,1,2,3,4};
  //std::vector< UInt_t > sqpart = {7,0,1,2,3,4};
  UInt_t ycutid = 3;

  //SAVEDATA
  //main
  bsaveData    = 0;
  //++------------------
  Bool_t bcorrelation = 0;
  //++-------------------
  //  bCentral = 1;
  Bool_t bdndydx      = 0; //dN/dy, dN/dpt, dN/d(beta*gamma) 
  Bool_t bdndy        = 0;
  Bool_t bdndydyx_one = 0;
  Bool_t bvertxz      = 0;
  Bool_t bgetvertxz   = 0;
  Bool_t bypt         = 0;
  Bool_t bdndyratio   = 0;
  Bool_t bdndpt       = 0;
  Bool_t bdnintegral  = 0;
  Bool_t bcluster     = 0; // Text ouput for multipliicity of clusters in AMD
  Bool_t bdndptratio  = 0;
  Bool_t bv1y         = 1;
  Bool_t bv2y         = 1;
  Bool_t bympx        = 0;
  Bool_t bpaper       = 0;
  Bool_t bTommy       = 0;
  //++------------------

  if( bsaveData )
    bcorrelation = 0;

  if( bcorrelation ) {
    bsaveData  = 0;
    bdndy      = 0;
    bdndydx    = 0;
    bvertxz    = 0;
    bgetvertxz = 0;
    bv1y       = 0;
    bv2y       = 0;
    bympx      = 0;
    bypt       = 0;
    bdndyratio = 0;
    bdndpt     = 0;
    bdnintegral= 0;
    bcluster   = 0;
    bdndptratio= 0;
    bpaper     = 0;
    bTommy     = 0;
    bCentral   = 1;
  }

  for( auto i : ROOT::TSeqI(4) ) {
    physData[i].resize(npart);
    for( auto j : ROOT::TSeqI(npart) ) {
      physData[i][j].resize((Int_t)(gnames.size()));
    }
  }
  for( auto i : ROOT::TSeqI(4) ) {
    physAMD[i].resize(npart);
    for( auto j : ROOT::TSeqI(npart) ) {
      physAMD[i][j].resize((Int_t)(nAMD));
      for( auto k : ROOT::TSeqI(nAMD) )  
	physAMD[i][j][k].resize(2);
    }
  }
  


  //--------------------------------------------------
  if( !bsaveData ) {
    if( bTommy ) {
      bCentral = 0;
      std::vector< UInt_t > sqsys1={0};
      std::vector< UInt_t > sqpart1={6};
      Draw_meanpx(sqsys1, sqpart1, gname, "gyv2ave",1,bsaveData); //corrected.
      bCentral = 1;
    }


    if( bpaper ) {
      //   {{0, "h_dndy"},{1,"h_dndEt"},{2,"h_dndbtgm"},{3,"h_dndyx"},{4,"h_dndpx"},{5,"h_dndpt"},
      //    {6,"hypx"}   ,{7,"gyv1"}  ,{8,"gyv2"}};

      Draw_PaperPlots(gname, -1, 1, 7, 42);  
      Draw_PaperPlots(gname, -1, 1, 8, 42); //v2

    }

    if( bdndy ) {
      for( auto isys : sqsys )for( auto ipart : sqpart ){
	  Draw_dndy(isys, ipart, gname);
	}
    }
  
    if( bdndyratio ) { 
      // 0: Mizuki, 1: Kaneko, 3: TTommy, 4: AMD
      std::vector<std::pair<UInt_t, UInt_t>> fpair={ {0,1} };
      for( auto isys : sqsys )for( auto ipart : sqpart ){
	  Draw_dndyRatio(isys, ipart, gname, fpair);
	}
    }
  
    //Double_t ycut[5]= {2.0, 0.5, 0.4, 0.3, 0.2};
    if( bvertxz ) {
      for( auto isys : sqsys )
	Draw_varxz(isys, gname, ycutid); //<-- The last index is ycut[id]
    }

    if( bgetvertxz ) {
      auto para = new Double_t[10];
      for( auto isys : sqsys )for( auto ipart : sqpart ) {
	  Getvarxz(isys,ipart,gname,para,ycutid); //<-- The last index is ycut[id]
	  for( auto ivar : ROOT::MakeSeq(8) ) {
	    LOG(INFO) << para[ivar] <<  FairLogger::endl;
	  }
	}
    }

    if( bympx ) {
      bCentral = 0;
      Draw_meanpx(sqsys, sqpart, gname, "hypx",41,bsaveData); //corrected.
      bCentral = 1;
    }

    if( bv1y ) {
      bCentral = 0;
      Draw_meanpx(sqsys, sqpart, gname, "gyv1",42,bsaveData); //corrected.
      //      Draw_meanpx(sqsys, spart, gname, "gu_v1",42,bsaveData); //corrected.
      //    Draw_meanpx(ssys, sqpart, gname, "gy_v1",41,bsaveData); //corrected.
      Draw_meanpx(sqsys, sqpart, gname, "hypx",42,bsaveData); //corrected.
      ////      Draw_meanpx(sqsys, sqpart, gname, "gyv1",42,bsaveData); //corrected.
      //Draw_meanpx(sqsys, sqpart, gname, "gyv1A",42,bsaveData); //corrected.
      bCentral = 1;

    }
    if( bv2y ) {
      bCentral = 0;
      Draw_meanpx(sqsys, sqpart, gnames, "gyv2",42,bsaveData); //corrected.
      //      Draw_meanpx(sqsys, sqpart, gname, "gyv2A",42,bsaveData); //corrected.
      bCentral = 1;
    }

    if( bdndydx ) {
      std::vector< UInt_t > spart = {0,1,2,3,4};
      Draw_dndydX(sqsys,sqpart,gname,ycutid, bsaveData,{0},42); 
      Draw_dndydX(sqsys,sqpart,gname,ycutid, bsaveData,{3},42); 
      Draw_dndydX(sqsys,sqpart,gname,ycutid, bsaveData,{2},42); 
    }

    if( bdnintegral ) 
      Draw_dnIntegral(sqsys,sqpart,gname,ycutid); //<-- The last index is ycut[id]

    if( bcluster ) {
      for( auto samd : AMDnames ) {
	if( samd.category != 1 ) continue;

	Double_t numClust[9][2]; 
	Draw_cluster(sqsys,sqpart,samd, numClust);
	if( numClust[0][0] != 0 ) {
	  cout << setw(25) << samd.config << "," ;
	  for( auto i : {0,1,2,3,4,7,5,6,8} ){
	    //	  cout << setw(4) << ncls[i].sName 
	    cout << setw(6)<<setprecision(4) << numClust[i][0] << "," << setw(6)<< setprecision(4) << numClust[i][1] << "," ;
	  }
	  cout << endl;
	}
      }
    }

    if( bdndpt ) 
      Draw_dndpt(sqsys,sqpart,gname);

    if( bdndptratio )
      Draw_dndptRatio(sqsys, gname );

    if( bypt ) {
      for( auto isys : {0,1,2,3} )
	Draw_ypt(isys, gname);
    }
  }
  //----- SAVE DATA ------------------
  else if( bsaveData ) {
    sqsys  = {0,1,2,3};
    gname  = gnames;
    sqpart = sqPart;

    bCentral = 0;
    Draw_meanpx(sqsys, sqpart, gname, "gyv1",42,bsaveData); //corrected.
    Draw_meanpx(sqsys, sqpart, gname, "gyv2",42,bsaveData); //corrected.
    bCentral = 1;
    Draw_dndydX(sqsys,sqpart,gname,ycutid, bsaveData,{0},42); 


    TString fname = Form("physData%s_left.dat",fphysdatafile[fphysdataid].second.Data());
    SavePhysData( fname, sqpart );

    fname = Form("physAMD%s.dat",fphysdatafile[fphysdataid].second.Data());
    SavePhysData( fname, sqpart );
  }
 

  //----- CORRELATION -----------------
  if( bcorrelation ) {

    if( 1 ) { //paper2022
      std::vector< UInt_t > ssys = {0};
      Draw_ParticleDependence(ssys, sqpart, gnames[0], 0, "v1Slope"        ,0,4,   0., 1.2, 42, 40);
      Draw_ParticleDependence(ssys, sqpart, gnames[0], 0, "v2minimum"      ,2,4, -0.1, 0.,  42, 40);
      ssys = {0,1};
      Draw_ParticleDependence(ssys, sqpart, gnames[0], 0, "v1Slope"        ,1,4,   0.45, 1.6, 40, 1);
      Draw_ParticleDependence(ssys, sqpart, gnames[0], 0, "v2minimum"      ,3,4,   0.,  1.8,  40, 1);
    }

    UInt_t ig = 0;
    for( auto igname : gnames ){

      if( 0 ) {
	bCentral = 1;
	if( igname.centrality == bCentral ) {
	  Draw_ParticleDependence(sqsys, sqpart, igname, bCentral, "RapiditystdDev" ,0,4,0.3 ,0.8);
	  Draw_ParticleDependence(sqsys, sqpart, igname, bCentral, "xRapiditystdDev",1,4,0.3, 0.8);
	  Draw_ParticleDependence(sqsys, sqpart, igname, 0,        "meanpxSlope"    ,2,4, 50.,165.);

	}
      }
      ig++;
    }

    //CORRELATION
    if( 0 ) { //paper2022
      std::vector< UInt_t > ssys = {0};
      Draw_ParticleDependence(sqsys, sqpart, gnames[1], 1, "RapiditystdDev" ,0,1,0.4 ,0.75, 41);
    }

    if( 0 ) { //paper2022
      std::vector< UInt_t > ssys = {0};
      Draw_ParticleDependence(ssys, sqpart, gnames[1], 1, "RapiditystdDev" ,0,2,0.36 ,0.74, 42);
      Draw_ParticleDependence(ssys, sqpart, gnames[1], 1, "xRapiditystdDev",1,2,0.36, 0.74, 42);
    }

    if( 0 ) { //paper2022
      std::vector< UInt_t > ssys = {0};
      Draw_ParticleDependence(ssys, sqpart, gnames[0], 0, "v1Slope"        ,0,2,   0., 1.2, 42);
      Draw_ParticleDependence(ssys, sqpart, gnames[0], 0, "v2minimum"      ,1,2, -0.1, 0.,  42);
    }


    if( 0 ) { //paper2022
      bCentral = 2;
      //      Draw_compCorrelation(sqsys, sqpart, gnames, 2,"pxstdDev", 0);
      sqsys={0};
      Draw_compCorrelation(sqsys, sqpart, gnames, 2,"meanbtgm", 0, 3);
      sqsys={1};
      Draw_compCorrelation(sqsys, sqpart, gnames, 2,"meanbtgm", 0, 4);
    } 

    if( 0 ) { //paper2022
      std::vector< UInt_t > ssys = {0,3,1};
      Draw_compIntegral(ssys, sqpart, 1, 42);
    }

    if( 0 ) { //paper2022
      std::vector< UInt_t > ssys = {3};
      Draw_compIntegral(ssys, sqpart, 1, 42);
    }

    if( 0 ) {
      // for( auto ipart : sqpart ) {
      std::vector< UInt_t > ssqpart = sqpart;
      Draw_compAMD(sqsys, ssqpart);
      //      }
    }

    if( 0 ) {
      Draw_compParameter(sqsys, sqpart, 0, "meanpxSlope");
    }
    if( 0 ) {
      Draw_compParameter(sqsys, sqpart, 1, "meanEt");
    }
    if( 0 ) {
      Draw_compParameter(sqsys, sqpart, 1, "meanbtgm");
    }
    if( 0 ) {
      Draw_compParameter(sqsys, sqpart, 1, "pxstdDev");
    }
    if( 0 ) {
      Draw_compParameter(sqsys, sqpart, 1, "integralM");
    }

    if( 0 ) {
      bCentral = 0;
      Draw_systematicError(sqsys,sqpart, gnames[0]);
    }
  }


  // Draw function
  UInt_t igname = 0;
  for( auto gname : gnames ) {

    //--- v11 and v20 system dependece

    if( nDATA == 1 ){

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
	Draw_v_y(igname, sqpart, 1);
      if( kFALSE ) 
	Draw_v_y(igname, sqpart, 2);


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
    for( auto ii : sqpart )
      Draw_dndy(3, ii);
  }
  
  //---  Multiplicity dependence
  if( nDATA > 3) {
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


