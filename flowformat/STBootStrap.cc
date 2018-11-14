#include "STBootStrap.hh"
STBootStrap::STBootStrap(UInt_t ival=1)
{
  clear();
  
  nboot  = ival;

  //  hbsphi = new TH1D("hbsphi","BootStrap distribution",100, -1. * TMath::Pi(), TMath::Pi() );
		    
}
STBootStrap::STBootStrap(UInt_t ival1, UInt_t ival2, Double_t *sample)
{
  clear();

  nboot       = ival1;
  numElements = ival2;


  for(UInt_t i = 0; i < numElements; i++)
    elements.push_back( sample[i] );

}


STBootStrap::STBootStrap(UInt_t ival1, UInt_t ival2, TVector2 *sample)
{
  clear();
  
  nboot       = ival1;
  numElementsTV2 = ival2;


  for(UInt_t i = 0; i < numElementsTV2; i++)
    elementsTV2.push_back( sample[i] );

  BootStrapping();
}

STBootStrap::STBootStrap(UInt_t ival1, std::vector<TVector2> *sample)
{
  clear();
  
  nboot       = ival1;
  numElementsTV2 = sample->size();

  elementsTV2 = *sample;

  BootStrapping();
}

void STBootStrap::Add(TVector2 sample)
{
  elementsTV2.push_back(sample);
  org_sum += sample.Unit();
  numElementsTV2 = (Int_t)elementsTV2.size();
}

void STBootStrap::Add(Double_t sample)
{
  elements.push_back(sample);
  r_sum += sample;
  numElements = (Int_t)elements.size();
}

UInt_t STBootStrap::BootStrapping(UInt_t nbt)
{
  if ( numElements > 0 ) return BootStrappingDouble(nbt);

  if ( numElementsTV2 > 0 ) return BootStrappingTVector2(nbt);

  if ( numElements < 0 && numElementsTV2 < 0 ) return -1;


  
  return 0;
}


UInt_t STBootStrap::BootStrappingDouble(UInt_t nbt)
{
  if(nbt > 0 ) nboot = nbt;
  else nboot = 100;
  
  ordn_Mean   = TMath::Mean(  elements.begin(), elements.end() );
  ordn_StdDev = TMath::StdDev(elements.begin(), elements.end()); 

  replace.clear();
  
  for(UInt_t i = 0; i < nboot; i++) {
    std::vector< UInt_t > rep = Resampling(numElements);

    Double_t sum_double = 0.;
    std::vector< Double_t > resampling;
    
    for( std::vector< UInt_t >::iterator it = rep.begin(); it != rep.end(); it++ ) 
      resampling.push_back( elements.at( *it ) );
    
    Double_t mean   = TMath::Mean( resampling.begin(), resampling.end() );    

    if( kFALSE ){
      for( std::vector< Double_t >::iterator it = resampling.begin(); it != resampling.end(); it++ )
	std::cout  << std::setw(10) << *it ;
      std::cout  << std::setw(10) << mean << std::endl;
    }

    replace.push_back( mean );

    // for median 
    // std::sort( resampling.begin(), resampling.end() );
    // Double_t median = *(resampling.begin() + resampling.size()/2);
    // if( resampling.size()%2 == 1 ) median = (median + *(resampling.begin() + resampling.size()/2 + 1))/2.;
    // replace.push_back( median );

    //    StoreResults;
    cnvMean = TMath::Mean(replace.begin(), replace.end()) ;
    cnvStdv = TMath::StdDev(replace.begin(), replace.end());

    resMean.push_back( cnvMean );
    //    cnvStdv = TMath::StdDev(resMean.begin(), resMean.end());
    resStdv.push_back( cnvStdv );

  }

  StoreConfideneLevel();

  return 1;
}

UInt_t STBootStrap::BootStrappingTVector2(UInt_t nbt)
{
  if(nbt > 0 ) nboot = nbt;
  else nboot = 1000;

  // Ordinary method
  ordn_Mean = TVector2::Phi_mpi_pi(org_sum.Phi());

  std::vector< Double_t > org_phiv;
  for(std::vector< TVector2 >::iterator it = elementsTV2.begin(); it != elementsTV2.end(); it++) 
    org_phiv.push_back( TVector2::Phi_mpi_pi( (*it).DeltaPhi( org_sum ) ) );

  
  ordn_StdDev = TMath::StdDev( org_phiv.begin(), org_phiv.end() );

  replace.clear();

  for(UInt_t i = 0; i < nboot; i++) {
    std::vector< UInt_t> rep = Resampling(numElementsTV2);

    TVector2 sum_vec = TVector2(0,0);
    TVector2 bs_vec  = TVector2(0,0);

    for(std::vector< UInt_t >::iterator it = rep.begin(); it != rep.end(); it++)  {
      TVector2 rotElements = elementsTV2.at( *it ).Rotate( -1.* org_sum.Phi() );
      sum_vec += rotElements.Unit();
      bs_vec  += elementsTV2.at( *it ).Unit();
    }      

    replace.push_back( TVector2::Phi_mpi_pi( sum_vec.Phi() ) );

    //    StoreResults();
    cnvMean = TVector2::Phi_mpi_pi( TMath::Mean(replace.begin(), replace.end()) );
    cnvStdv = TMath::StdDev(replace.begin(), replace.end()); ;

    resMean.push_back( cnvMean );
    resStdv.push_back( cnvStdv );

    cnvMean += org_sum.Phi() ;

    cnvMod = org_sum.Mod();
  }  

  StoreConfideneLevel();

  return 1;
}


void STBootStrap::StoreConfideneLevel()
{
  std::sort(replace.begin(), replace.end());

  UInt_t CL_25  = UInt_t( (Double_t)replace.size() * 0.025);
  UInt_t CL_975 = UInt_t( (Double_t)replace.size() * 0.975);

  std::vector< Double_t >::iterator itr;
  itr = replace.begin();

  CL_25lw  = *(itr + CL_25 );

  CL_975up = *(itr + CL_975 - 1);


  Error = (CL_975up - CL_25lw) * cnvStdv;
}

void STBootStrap::clear()
{
  numElements = 0;
  numElementsTV2 = 0;
  elements.clear();
  elementsTV2.clear();

  replace.clear();
  resMean.clear();
  resStdv.clear();
  //  replaceMod.clear();

  cnvMean    = 0.;
  cnvStdv    = 0.;
  cnvStdv2   = 0.;
  cnvCosMean = 0.;
  cnvMod     = 0.;

  org_sum = TVector2(0.,0.);
  r_sum    = 0.;
  
  ordn_StdDev = 0.;
  
  CL_25lw = 0.;
  CL_975up= 0.;
  Error   = 0.;
}

std::vector< UInt_t > STBootStrap::Resampling(UInt_t ival)
{
  std::vector< UInt_t > rep;
  
  for(UInt_t j = 0; j < ival; j++)
    rep.push_back( rnd.Integer(ival) );
  
  return rep;
}   

Double_t STBootStrap::GetReplace(UInt_t ival)
{
  if( replace.size() > 0 && ival < replace.size() ) return replace.at(ival);

  else return -999.;
}


Double_t STBootStrap::GetResidualMean(UInt_t ival)
{

  if( resMean.size() > 0 && ival < resMean.size()) return resMean.at(ival) ; 

  else return -999.;
}

Double_t STBootStrap::GetResidualStdDev(UInt_t ival)
{

  if( resStdv.size() > 0 && ival < resStdv.size()) return resStdv.at(ival); 

  else return -999.;
}


Double_t STBootStrap::GetStdDevError()
{
  UInt_t vsize = 1.;
  return cnvStdv / sqrt( Double_t(vsize) ); 
}

