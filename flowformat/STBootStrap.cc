#include "STBootStrap.hh"
STBootStrap::STBootStrap(UInt_t ival)
{
  clear();
  
  nboot  = ival;
}

STBootStrap::STBootStrap(UInt_t ival1, UInt_t ival2, Double_t *sample)
{
  clear();

  nboot       = ival1;
  numElements = ival2;


  for(UInt_t i = 0; i < numElements; i++)
    elements.push_back( sample[i] );

  for(UInt_t i = 0; i < nboot; i++) {

    std::vector< UInt_t> rep = Resampling(numElements);

    replace.clear();
    for(std::vector< UInt_t >::iterator it = rep.begin(); it != rep.end(); it++)      
      replace.push_back( elements.at( *it ) );
			 
    StoreResults();
  }
}


STBootStrap::STBootStrap(UInt_t ival1, UInt_t ival2, TVector2 *sample)
{
  clear();
  
  nboot       = ival1;
  numElements = ival2;


  for(UInt_t i = 0; i < numElements; i++)
    elementsTV2.push_back( sample[i] );

  BootStrapping();
}

STBootStrap::STBootStrap(UInt_t ival1, std::vector<TVector2> *sample)
{
  clear();
  
  nboot       = ival1;
  numElements = sample->size();

  elementsTV2 = *sample;

  BootStrapping();
}

void STBootStrap::Add(TVector2 sample)
{
  elementsTV2.push_back(sample);

  numElements = (Int_t)elementsTV2.size();
}


UInt_t STBootStrap::BootStrapping(UInt_t nbt)
{
  if ( numElements <= 0 ) return -1;


  std::cout << " num elements " << numElements << std::endl;

  if(nbt > 0 ) nboot = nbt;
  else nboot = 100;

  std::vector<TVector2> rvec;
  TVector2 rvec_sum = TVector2(0,0);

  for(UInt_t ielm = 0; ielm < numElements; ielm++) { 
    rvec.push_back( elementsTV2.at(ielm) );
    rvec_sum += elementsTV2.at(ielm);
  }
  
  phi_off = rvec_sum.Phi();


  replace.clear();

  for(UInt_t i = 0; i < nboot; i++) {
    std::vector< UInt_t> rep = Resampling(numElements);
    

    TVector2 sum_vec = TVector2(0,0);

    for(std::vector< UInt_t >::iterator it = rep.begin(); it != rep.end(); it++)  {

      sum_vec += rvec.at( *it ).Unit(); 
    }      

    Double_t vec_delt = rvec_sum.DeltaPhi( sum_vec );
    replace.push_back( TVector2::Phi_mpi_pi( vec_delt ) );


    StoreResults();

  }  

  
  return numElements;
}


void STBootStrap::StoreResults()
{

  std::vector< Double_t >::iterator ibgn;
  std::vector< Double_t >::iterator iend;
    
  ibgn = replace.begin();
  iend = replace.end();

  resMean.push_back( TMath::Mean(ibgn, iend) );
  resStdv.push_back( TMath::StdDev(ibgn, iend) );
  
  std::cout << " res " << iend - ibgn <<  " std " << TMath::StdDev(ibgn, iend) <<  std::endl;

  //  ibgn = resMean.begin();
  //  iend = resMean.end();

  cnvMean = TVector2::Phi_mpi_pi( TMath::Mean(ibgn, iend) + phi_off ) ;
  cnvStdv = TMath::StdDev(ibgn, iend) ;
  std::cout << " mean std " << cnvMean << " " << cnvStdv << std::endl;

  cnvCosMean = cos(TMath::Mean(ibgn, iend) ) ;
}

void STBootStrap::clear()
{
  elements.clear();
  elementsTV2.clear();

  replace.clear();
  resMean.clear();
  resStdv.clear();
  cnvMean = 0.;
  cnvStdv = 0.;
  cnvStdv2 = 0.;
  cnvCosMean = 0.;

}

std::vector< UInt_t > STBootStrap::Resampling(UInt_t ival)
{
  std::vector< UInt_t > rep;
  
  for(UInt_t j = 0; j < ival; j++)
    rep.push_back( rnd.Integer(ival) );
  
  return rep;
}   


Double_t STBootStrap::GetResidualMean(UInt_t ival)
{

  if( resMean.size() > 0 && ival < resMean.size()) return TVector2::Phi_mpi_pi( resMean.at(ival) + phi_off); 

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


