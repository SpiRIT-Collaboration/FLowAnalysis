#include "STSource.hh"

STSource* STSource::fgRinstance = 0;

STSource* STSource::Instance()
{
  return fgRinstance;
}

STSource::STSource()
{
  fgRinstance = this;
}


UInt_t STSource::SetInputChain(TChain *fchain)
{
  auto chname = (std::string)fchain->GetName();

  fInputChain.insert(std::make_pair(chname,  fchain));
  
  return (UInt_t)fInputChain.size();
}

TChain* STSource::GetInputChain(TString fcname)
{
  std::map<std::string, TChain*>::iterator iiter = fInputChain.begin();

  // while( iiter != fInputChain.end() ){
  //   if( *(iiter->first) == fcname.c_str() )
  //     return iiter->second;
  //   else
  //     iiter++;
  // }

  return NULL;
}

