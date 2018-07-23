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

    replace[0].clear();
    for(std::vector< UInt_t >::iterator it = rep.begin(); it != rep.end(); it++)      
      replace[0].push_back( elements.at( *it ) );
			 
    StoreResults(0);
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


void STBootStrap::BootStrapping(UInt_t nbt)
{
  Double_t hrm = 1.;

  if(nbt != 0 ) nboot = nbt;

  std::vector<TVector2> rvec;
  TVector2 rvec_sum = TVector2(0,0);
  for(UInt_t ielm = 0; ielm < numElements; ielm++) { 
    rvec.push_back( elementsTV2.at(ielm) );
    rvec_sum += elementsTV2.at(ielm);

  }


  Double_t phi_rot = 0;

  for(UInt_t i = 0; i < nboot; i++) {
    std::vector< UInt_t> rep = Resampling(numElements);
    
    for(UInt_t j = 0; j < nidx; j++)
      replace[j].clear();

    TVector2 vsum = TVector2(0,0);

    for(std::vector< UInt_t >::iterator it = rep.begin(); it != rep.end(); it++)      
      vsum += rvec.at( *it );


      

    // Double_t tsumX = 0;
    // Double_t tsumY = 0;
    
    // for(std::vector< UInt_t >::iterator it = rep.begin(); it != rep.end(); it++)      
    //   tsumX += cos(hrm * rvec.at( *it ).Phi() );
    // //    tsumX /= (Double_t)numElements;

    // replace[0].push_back( tsumX );
    // StoreResults(0);

    // for(std::vector< UInt_t >::iterator it = rep.begin(); it != rep.end(); it++)      
    //   tsumY += sin(hrm * rvec.at( *it ).Phi() );

    // //  tsumY /= (Double_t)numElements;
   
    // replace[1].push_back( tsumY );

    // StoreResults(1);

    // TVector2 vt2 = TVector2(tsumX, tsumY);
    //    replace[2].push_back( vt2.DeltaPhi( rvec_sum )); 

 
    phi_mean =  vsum.Phi() ;
    replaceVec.push_back( vsum );
    StoreVectors();
    
  }  
}


void STBootStrak::StoreVectors(TVector2 vec)
{
  std::vector< Double_t >::iterator ibgn;
  std::vector< Double_t >::iterator iend;

  ibgn = replaceVec.begin();
  iend = replaceVec.end();

  


  resVec.push_back( vec );

  
}


void STBootStrap::StoreResults(UInt_t idx, Double_t off)
{

  std::vector< Double_t >::iterator ibgn;
  std::vector< Double_t >::iterator iend;
    
  ibgn = replace[idx].begin();
  iend = replace[idx].end();
  resMean[idx].push_back( TMath::Mean(ibgn, iend) );
  resStdv[idx].push_back( TMath::StdDev(ibgn, iend) );
  
  //  final results
  ibgn = resMean[idx].begin();
  iend = resMean[idx].end();

  cnvMean[idx].push_back(        TMath::Mean(ibgn, iend) + off );
  cnvStdv[idx].push_back(        TMath::StdDev(ibgn, iend) );
  cnvCosMean[idx].push_back( cos(TMath::Mean(ibgn, iend) ) );

  ibgn = resStdv[idx].begin();
  iend = resStdv[idx].end();
  cnvStdv2[idx].push_back( TMath::StdDev(ibgn, iend) );



}

void STBootStrap::clear()
{
  elements.clear();
  elementsTV2.clear();

  Initialize();
  for(UInt_t i = 0; i < nidx; i++){
    replace[i].clear();
    resMean[i].clear();
    resStdv[i].clear();
    cnvMean[i].clear();
    cnvStdv[i].clear();
    cnvStdv2[i].clear();
    cnvCosMean[i].clear();
  }
}

std::vector< UInt_t > STBootStrap::Resampling(UInt_t ival)
{
  std::vector< UInt_t > rep;
  
  for(UInt_t j = 0; j < ival; j++)
    rep.push_back( rnd.Integer(ival) );
  
  return rep;
}   


Double_t STBootStrap::GetMean(UInt_t idx)
{
  if( cnvMean[idx].size() > 0 ) return TVector2::Phi_mpi_pi( cnvMean[idx].at( cnvMean[idx].size() - 1 ) ); 

  else    return -999.;
}

Double_t STBootStrap::GetCosMean(UInt_t idx)
{
  if( cnvCosMean[idx].size() > 0 ) return cnvCosMean[idx].at( cnvCosMean[idx].size() - 1  ); 

  else    return -999.;
}

Double_t STBootStrap::GetStdDev(UInt_t idx)
{
  if( cnvStdv[idx].size() > 0 ) return cnvStdv[idx].at( cnvStdv[idx].size() - 1 ); 

  else    return -999.;
}

Double_t STBootStrap::GetStdDevError(UInt_t idx)
{
  UInt_t vsize = cnvStdv[idx].size();
  if( cnvStdv[idx].size() > 0 ) return cnvStdv[idx].at( vsize - 1 ) / sqrt( Double_t(vsize) ); 

  else    return -999.;
}

Double_t STBootStrap::GetStdDev2(UInt_t idx)
{
  if( cnvStdv2[idx].size() > 0 ) return cnvStdv2[idx].at( cnvStdv2[idx].size() - 1 ); 

  else    return -999.;

}


void STBootStrap::Initialize()
{
  replace.resize(nidx);
  resMean.resize(nidx);
  resStdv.resize(nidx);
  cnvMean.resize(nidx);
  cnvCosMean.resize(nidx);
  cnvStdv.resize(nidx);
  cnvStdv2.resize(nidx);
}
