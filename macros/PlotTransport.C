
UInt_t    ndata;
Double_t  yz[50];
Double_t  rms_yz[50];
Double_t  v1[50];
Double_t  rms_v1[50];
Double_t  v2[50];
Double_t  rms_v2[50];
TGraphErrors *gu_v1; 
TGraphErrors *gu_v2; 

Bool_t bdebug = 0;

TString FindTransportData(UInt_t i=0, UInt_t j=0, UInt_t k=0) //model, y OR pt, particle
{
  TString fdir[] = {"AMDPt","ChiBUUPt","IBUUPt","IQMDPt","SMFPt","UrQMDPt","dcQMDPt"};
  TString dname[] = {"rapflow-","ptdist-"};
  TString pname[] = {"p","d","t"};

  TString fname = "data/All_models/" + fdir[i] ;

  if( i == 3 )
    fname +=  "/IQMD-BNU-" + dname[j] + pname[k] + ".dat";
  else 
    fname += "/" + dname[j] + pname[k] + ".dat"; 

 
    std::cout << fname << " is not found." << std::endl;


  if( gSystem -> FindFile(".",fname ) ) return fname;
  
  return "";
  
}

TGraphErrors *ReadTransportData(TString fname, UInt_t iv)
{

  std::fstream TData;
  TData.open(fname, std::fstream::in);

  TString sget;
  UInt_t iy = 0;
  while( !TData.eof() ) {
    TData >> sget;
    if( sget == "yz" || sget == "v1" || sget == "rms_v1" || sget == "v2" || sget == "rms_v2" || sget == "0") continue;


    yz[iy] = atof(sget) ;
    if( bdebug )    cout <<" yz " << yz[iy] ;

    TData >> sget;
    v1[iy] = atof(sget) ;
    if( bdebug ) cout << " v1 " << v1[iy] ;

    TData >> sget;
    rms_v1[iy] = atof(sget) ;
    if( bdebug ) cout << " rms_v1 " << rms_v1[iy] ;

    TData >> sget;
    v2[iy] = atof(sget) ;
    if( bdebug ) cout << " v2 " << v2[iy] ;

    TData >> sget;
    rms_v2[iy] = atof(sget) ;
    if( bdebug ) cout << " rms_v2 " << rms_v2[iy] ;
    if( bdebug ) cout << endl;

    rms_yz[iy] = 0. ;



    iy++;
  }

  ndata=iy;

  if( fname.Contains("dcQMD") || fname.Contains("IBUU") ){
    for(UInt_t i = 0; i < ndata; i++ ) 
      yz[i] = -yz[i];

  }


  auto gu_v1 = new TGraphErrors( ndata, yz, v1, rms_yz, rms_v1 );
  auto gu_v2 = new TGraphErrors( ndata, yz, v2, rms_yz, rms_v2 );

 

  TData.close();

  if( iv == 2 )
    return gu_v2;
  else
    return gu_v1;

}

void ReadTransportPtData(TString fname, UInt_t iv)
{
  std::fstream TData;
  TData.open( fname, std::fstream::in );
  
  std::cout << fname << " is opened. " << std::endl;

  TString sget;
  UInt_t ipt = 0;
  Double_t pt[40];
  Double_t pt_v1[40][80];

  TData >> sget;
  cout << sget;

  for(UInt_t iy = 0; iy < 40; iy++ ){

    TData >> sget;
    cout << sget << " ";
    yz[iy] = atof(sget);

  }
  cout << "-----" << endl;

  for(UInt_t ipt = 0; ipt < 40; ipt++) {

    TData >> sget;
    pt[ipt] = atof(sget);
    cout << " (" << ipt << ") " << sget << " pt " << pt[ipt] <<endl;
      
    for(UInt_t iyy = 0; iyy < 80; iyy++ ){
      TData >> sget;
      cout << sget;
      
      pt_v1[iyy][ipt] = atof(sget);
      cout << " <" << iyy << "> " << sget << " v1 " << pt_v1[iyy][ipt] <<endl;
    }
  }    


  std::cout << endl;

  for( UInt_t k = 0; k < 2; k++ ) {
    std::cout << k << " : " << " rap " << yz[k] << std::endl;

    for( UInt_t j = 0; j < 40; j++ ) {
      std::cout << j << " pt " << pt[j] << " : v1 " << pt_v1[k][j] << std::endl;
    }
  }
  
  

  // auto gu_v1 = new TGraphErrors( ndata, yz, v1, rms_yz, rms_v1 );
  // auto gu_v2 = new TGraphErrors( ndata, yz, v2, rms_yz, rms_v2 );

  TData.close();
  
}


void PlotTransport(UInt_t i=0, UInt_t j=1, UInt_t k=0)
{
  //  ReadTransportPtData("/cache/scr/spirit/mizuki/SpiRITAnalysis/macros/./data/All_models/AMDPt/ptdist-p.dat",1);
}

 
// void temp() 
// {
//   TGraphErrors *gu_v1 = new TGraphErrors();
//   TGraphErrors *gu_v2 = new TGraphErrors();

//   gu_v1 = ReadTransportData( FindTransportData(i, j, k) , 1);
//   gu_v1 -> SetName("gu_v1");

//   gu_v2 = ReadTransportData( FindTransportData(i, j, k) , 2);
//   gu_v2 -> SetName("gu_v2");


//   TCanvas *cc;
//   UInt_t icc = 0;

//   cc = new TCanvas(Form("cc%d",icc), Form("cc%d",icc) ); icc++;
//   gu_v1->Draw("ALP");

//   cc = new TCanvas(Form("cc%d",icc), Form("cc%d",icc) ); icc++;
//   gu_v2->Draw("ALP");

// }
