#include "STFlowCorrection.hh"

void STFlowCorrection::Initialize(TChain *chele, UInt_t ival1, UInt_t ival2)
{
  ChEle = chele;
  Initialize(ival1, ival2);
}

void STFlowCorrection::Initialize(UInt_t ival1, UInt_t ival2)
{
  harm  = ival1;
  charm = harm;

  //irm   = ival2;
  iver  = ival2;

  Init();
}

void STFlowCorrection::Init()
{
  indx   = new UInt_t[harm];
  An     = new Double_t[harm];
  Bn     = new Double_t[harm];
  An_rms = new Double_t[harm];
  Bn_rms = new Double_t[harm];

  for(UInt_t i = 0; i < harm; i++){
    indx[i]   = 0;
    An[i]     = 0.;
    Bn[i]     = 0.;
    An_rms[i] = 0.;
    Bn_rms[i] = 0.; 
  }

  for(UInt_t i = 0; i < 2; i++){
    binmax[i] = 0.;;
    binmin[i] = 0.;;
    binpara[i]= "";
  }

  constX   = 0.;
  meanX    = 0.;
  sigX     = 0.;
  constY   = 0.;
  meanY    = 0.;
  sigY     = 0.;


  if(ChEle != NULL)  SetFileName();  

  clear();
}

void STFlowCorrection::clear()
{
  vphi.clear();
  bphi.clear();
  vtheta.clear();
  vmtrack.clear();

  vvec.clear();
  bvec.clear();
  rcphi.clear();

}

void STFlowCorrection::SetFileName()
{
  TIter nnext((TCollection*)ChEle->GetListOfFiles());

  TString sval = ((TFile)nnext()->GetTitle()).GetName();


  Ssiz_t ifnd = sval.First("_")-7;
  fname = sval(ifnd,sval.Length()-ifnd-4);
  
}

void STFlowCorrection::SetFileName(TString sval)
{
  fname = sval;
  LOG(INFO) << fname << " in defined. " << FairLogger::endl;
}

void STFlowCorrection::SetReCenteringParameter(TString cprm, Double_t *val)
{
  if( cprm == "X"){
    constX = val[0];
    meanX  = val[1];
    sigX   = val[2];
  }
  else if( cprm == "Y"){
    constY = val[0];
    meanY  = val[1];
    sigY   = val[2];
  }
  else
    LOG(ERROR) << "Failed:  SetReCenteringParameter(X or Y, " << FairLogger::endl;
}


UInt_t STFlowCorrection::LoadCorrectionFactor(UInt_t val)
{

  LOG(INFO) << "STFlowCorrection::GetCorrectionFactor : "<< FairLogger::endl;

  TString header = "->,";

  std::fstream fin;
  fin.open(fname, std::fstream::in);
  if(fin == NULL) {
    LOG(ERROR) << "A file " << fname << " was not found " << FairLogger::endl;
    exit(0);
  }


  TString sget;
  
  fin >> sget;
  harm = atoi(sget);
  charm = harm;  // number of harmonics


  Init();

  Int_t j = 0;
  while(!fin.eof()){
    fin >> sget;

    if(sget == "X:"){
      fin >> sget; constX = atof(sget);
      fin >> sget; meanX  = atof(sget);
      fin >> sget; sigX   = atof(sget);

      fin >> sget;
      fin >> sget; constY = atof(sget);
      fin >> sget; meanY  = atof(sget);
      fin >> sget; sigY   = atof(sget);

      LOG(DEBUG) << " ReCentering correction factors > " ;
      LOG(DEBUG) << " X: " 
		 << constX << ", "
		 << meanX  << ", "
		 << sigX   << ", " ;
      LOG(DEBUG) << " Y: " 
		 << constY << ", "
		 << meanY  << ", "
		 << sigY   << ", "
		 << FairLogger::endl;

    }
    
    if(sget == "mtrack>"){
      fin >> sget;
      binmin[0] = atof(sget);
      binpara[0]= "mtrack";
    }
    else if(sget == "mtrack<"){
      fin >> sget;
      binmax[0] = atof(sget);
      binpara[0]= "mtrack";
    }

    if(sget == "ntrack>="){
      fin >> sget;
      binmin[0] = atof(sget);
      binpara[0]= "ntrack";
    }
    else if(sget == "ntrack<"){
      fin >> sget;
      binmax[0] = atof(sget);
      binpara[0]= "ntrack";
    }

    if(sget == "pz<"){
      fin >> sget;
      binmax[0] = atof(sget);
      binpara[0] = "pz";
    }
    else if(sget =="theta<"){
      fin >> sget;
      binmax[1] = atof(sget);
      binpara[1] = "theta";
    }
    else if(sget =="theta>="){
      fin >> sget;
      binmin[1] = atof(sget);
      binpara[1] = "theta";
    }

    
    else if(sget == header){
      fin >> sget;
      indx[j] = atoi(sget);

      fin >> sget;
      Bn[j] = atof(sget);

      fin >> sget;
      Bn_rms[j] = atof(sget);

      fin >> sget;
      An[j] = atof(sget);

      fin >> sget;
      An_rms[j] = atof(sget);

      j++;
    }
  }

  if( val == 1 ) ShowParameters();

  fin.close();
  return 0;
}

void STFlowCorrection::ShowBinInformation()
{
  LOG(DEBUG) << binpara[0] << " > " << binmin[0] << " && < " << binmax[0] << FairLogger::endl; 
  if(binpara[1] != "")
    LOG(DEBUG) << binpara[1] << " > " << binmin[1] << " && < " << binmax[1] << FairLogger::endl; 

}

void STFlowCorrection::ShowParameters()
{  
  ShowBinInformation();

  for(UInt_t k = 0; k < charm; k++){ 
    LOG(DEBUG) << "->, " << std::setw(5) << indx[k] << ", "
	       << std::scientific << std::setprecision(5) << std::right
	       << std::setw(20) << Bn[k] << ", " << Bn_rms[k] << ", "
	       << std::setw(20) << An[k] << ", " << An_rms[k]
	       << FairLogger::endl;
  }
   
  LOG(INFO) << fname << " was loaded." << FairLogger::endl;
}


TVector3 STFlowCorrection::ReCentering(TVector3 val)
{
  return TVector3((val.X() - meanX)/sigX, (val.Y() - meanY)/sigY, val.Z());
}

TVector3 STFlowCorrection::ReCenteringFourierCorrection(TVector3 val) // ReCentering and Fourier correction for stored bvec
{
  TVector3 rc = ReCentering(val);
  rc.SetPhi( GetCorrection(rc.Phi()) );

  return rc;
}

UInt_t STFlowCorrection::ReCenteringFourierCorrection() // ReCentering and Fourier correction for stored bvec
{
  std::vector<TVector3>::iterator itr;
  for(itr = bvec.begin(); itr != bvec.end(); itr++){    
    if(sigX > 0 && sigY > 0){
      TVector3 rc = ReCentering( *itr );
      vvec.push_back(rc);
      vphi.push_back(rc.Phi());
      rcphi.push_back(rc.Phi());
    }
  }

  return FourierCorrection();
}


void STFlowCorrection::GetCorrection(std::vector<Double_t> &val)
{
  for(Int_t i = 0; i < (Int_t)val.size(); i++) {
    val.at(i) = GetCorrection(val.at(i)) ;
  }
}

Double_t STFlowCorrection::GetCorrection(Double_t val)
{
  Double_t cpphi = val;

  for(UInt_t k = 0; k < charm; k++){
    cpphi += An[k]*cos((Double_t)(k+1) * val);
    cpphi += Bn[k]*sin((Double_t)(k+1) * val);
  } 
  return cpphi;  
}

UInt_t STFlowCorrection::FourierCorrection()
{
  std::vector<TVector3>::iterator itr;

  if(vphi.size() == 0){
    for(itr = bvec.begin(); itr != bvec.end(); itr++)
      vphi.push_back(itr->Phi());
  }

  FourierCorrection(vphi);
  
  for(itr = bvec.begin(); itr != bvec.end(); itr++) 
    itr->SetPhi(vphi.at( itr - bvec.begin() ) );
  
  return (UInt_t)vphi.size();
}

void STFlowCorrection::FourierCorrection(std::vector<Double_t> &val)
{
  LOG(DEBUG) << " harm = " << charm <<  "val.size " << val.size() << FairLogger::endl;
  std::vector<Double_t> fvCos;
  std::vector<Double_t> fvSin;

  for(UInt_t ihm = 0; ihm < charm; ihm++){ 

    fvCos.clear();
    fvSin.clear();

    Double_t findx = (Double_t)(ihm+1);

    for(Int_t ival = 0; ival < (Int_t)val.size(); ival++){

      fvCos.push_back(cos(findx * val.at(ival)));
      fvSin.push_back(sin(findx * val.at(ival)));

    }

    if( fvCos.size() > 0 ) {
      std::vector<Double_t>::iterator ibgn = fvCos.begin();
      std::vector<Double_t>::iterator iend = fvCos.end();

      Bn[ihm]     =  2./findx * TMath::Mean(ibgn, iend); 
      Bn_rms[ihm] =  2./findx * TMath::RMS(ibgn, iend); 

      ibgn = fvSin.begin();
      iend = fvSin.end();

      An[ihm]     = -2./findx * TMath::Mean(ibgn, iend);  
      An_rms[ihm] =  2./findx * TMath::RMS(ibgn, iend);  
    }
    else{
      An[ihm] = 0.;
      Bn[ihm] = 0.;
      An_rms[ihm] = 0.;
      Bn_rms[ihm] = 0.;
    }
  }

  for(UInt_t k = 0; k < charm; k++){
    Double_t findx = (Double_t)(k+1);
    
    LOG(DEBUG) << std::setw(5) << std::noshowpos << k+1 
	       << std::scientific << std::setprecision(5)  << std::right //<< showpos
	       << std::setw(6) << " Bn<cos> : " <<  std::setw(15) << Bn[k]  << " rms " << Bn_rms[k]
	       << "    " 
	       << std::setw(6) << " An<sin> : " <<  std::setw(15) << An[k]  << " rms " << An_rms[k]
	       << FairLogger::endl;
    
  }
  
  if(val.size()>0) 
    GetCorrection(val);

}

std::vector<Double_t> STFlowCorrection::GetTheta()
{
  std::vector<TVector3>::iterator itr;
  for(itr = vvec.begin(); itr != vvec.end(); itr++)
    vtheta.push_back(itr->Theta());
						   
  return vtheta;
}

Double_t  *STFlowCorrection::GetAverageCosin(Int_t ival, std::vector<Double_t> &val)
{
  std::vector<Double_t> fvCos;

  fvCos.clear();

  Double_t findx = (Double_t)ival;
  for(Int_t i = 0; i < (Int_t)val.size(); i++)
    fvCos.push_back(cos(findx * val.at(i)));
  
  std::vector<Double_t>::iterator ibgn = fvCos.begin();
  std::vector<Double_t>::iterator iend = fvCos.end();

  Double_t *vcos;
  vcos = new Double_t[2];
  vcos[0]  =  TMath::Mean(ibgn, iend);
  vcos[1]  =  TMath::RMS(ibgn, iend);

  return vcos;
}

void STFlowCorrection::SetDirectory()
{
  gSystem->cd("db");
}

UInt_t STFlowCorrection::SaveCorrectionFactor(TString comm1, TString comm2)
{ 
  SetDirectory();
  
  std::fstream fout;

  Ssiz_t ifnd = comm1.First(":");
  fname = comm1(0,ifnd);
  fname += ".txt";

  fout.open(fname,std::fstream::out);
  fout << charm << std::endl;

  fout << "comment : " << comm1 << std::endl;

  if( comm2 != "")
    fout << "ReCentering parameter : " << comm2 << std::endl;

  TString comm3 = Form("X: %f, %f, %f Y: %f, %f, %f :",constX,meanX,sigX,constY,meanY,sigY);
  fout << comm3 << std::endl;


  TChainElement *ele = 0;
  UInt_t ith = 0;
  TIter nnext((TCollection*)ChEle->GetListOfFiles());

  while(( ele = (TChainElement*)nnext() )){
    TFile f(ele->GetTitle());
    fout << ith << " : " << f.GetName() << std::endl;
    ith++;
  }

  fout << "nth  "  << "Bn<cos>        " << "rms        " <<"An<sin>        " << "rms" << std::endl;

  for(UInt_t k = 0; k < charm; k++){
    fout << "->, " << std::setw(5) << k+1 << ", "
   	 << std::scientific << std::setprecision(5) << std::right
   	 << std::setw(20) << Bn[k] << ", " << Bn_rms[k] << ", "
   	 << std::setw(20) << An[k] << ", " << An_rms[k]
   	 << std::endl;
  }  

  fout.close();
  LOG(INFO) << fname << " is created."<< FairLogger::endl;


  gSystem->cd("..");
  return 1;
}


void STFlowCorrection::PrintContents()
{
  PrintRange();

  for(UInt_t i = 0; i < (UInt_t)vphi.size(); i++){
    LOG(INFO) << " mtrack " << std::setw(5) << vmtrack.at(i) << " theta " << std::setw(8) << vtheta.at(i) 
   	      << " phi " << vphi.at(i) << FairLogger::endl;
  } 
}

void STFlowCorrection::PrintRange()
{
  std::cout << fname << FairLogger::endl;

  if(binpara[0] != "")
    LOG(INFO) << binmin[0] << " < " << binpara[0] << " < " << binmax[0] << FairLogger::endl; 
  if(binpara[1] != "")
    LOG(INFO) << binmin[1] << " < " << binpara[1] << " < " << binmax[1] << FairLogger::endl; 
}
#if !defined(__CINT__)
ClassImp(STFlowCorrection);
#endif





