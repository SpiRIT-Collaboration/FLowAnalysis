#include "STNeuLANDCluster.hh"

#include <TMath.h>
#include <iostream>

ClassImp(STNeuLANDCluster)

STNeuLANDCluster::STNeuLANDCluster()
{
  Clear();
  distance_to_center =  8561.05 + 280; // Mizuki's photogrametry to the blue frame + to the 1st plane surface
  angle = 29.579 * TMath::Pi()/180.; // deg.
  
  global_offset = TVector3(distance_to_center*TMath::Sin(angle), 0., 
			   distance_to_center*TMath::Cos(angle));


  c = 29.979245; //cm/ns
  nmass = 939.565; // MeV/cc
}

void STNeuLANDCluster::Clear(Option_t *)
{
  nhit = -1;
  edep = -1;
  tof = -1;
  tof_last = -1;

  localx = -500; localy = -500; localz = -1;
  localx_last = -500; localy_last = -500; localz_last = -1;
  distance = -1;
  localPos.SetX(localx);
  localPos.SetY(localy);
  localPos.SetZ(localz);
  
  localPos_last.SetX(localx_last);
  localPos_last.SetY(localy_last);
  localPos_last.SetZ(localz_last);

  globalPos.SetX(-10000);
  globalPos.SetY(-10000);
  globalPos.SetZ(-10000);

  veto_all = -1;
  veto_bar = -1;
  veto_mid = -1;
  veto_loose = -1;
}

void STNeuLANDCluster::SetLocalPos()
{
  localPos.SetX(localx);
  localPos.SetY(localy);
  localPos.SetZ(localz);
  
  if(localx > -500){
    globalPos = localPos + global_offset;

    std::cout << " gx " << globalPos.X() 
	      << " lx " << localPos.X()
	      << " off " << global_offset.X()
	      << std::endl;
    //    globalPos.RotateY(-angle);

    

  }
  else{
    globalPos.SetX(-500);
    globalPos.SetY(-500);
    globalPos.SetZ(-1);
  }
}


void STNeuLANDCluster::SetLocalPosLast()
{
  localPos_last.SetX(localx);
  localPos_last.SetY(localy);
  localPos_last.SetZ(localz);
  
  if(localx > -500){
    globalPos_last = localPos + global_offset;
    globalPos_last.RotateY(-angle);
  }
  else{
    globalPos_last.SetX(-500);
    globalPos_last.SetY(-500);
    globalPos_last.SetZ(-1);
  }
}

Double_t STNeuLANDCluster::GetGlobalX()
{
  return localx*TMath::Cos(angle) + (localz+distance_to_center)*TMath::Sin(angle);
}

Double_t STNeuLANDCluster::GetGlobalY()
{
  return localy;
}

Double_t STNeuLANDCluster::GetGlobalZ()
{
  return -localx*TMath::Sin(angle) + (localz+distance_to_center)*TMath::Cos(angle);
}

Double_t STNeuLANDCluster::GetPx()
{
  if(0>distance) distance = TMath::Sqrt((distance_to_center + localz)*(distance_to_center + localz) + localx*localx + localy*localy);
  return GetMom() * GetGlobalX()/distance;
}

Double_t STNeuLANDCluster::GetPy()
{
  if(0>distance) distance = TMath::Sqrt((distance_to_center + localz)*(distance_to_center + localz) + localx*localx + localy*localy);
  return GetMom() * GetGlobalY()/distance;
}

Double_t STNeuLANDCluster::GetPz()
{
  if(0>distance) distance = TMath::Sqrt((distance_to_center + localz)*(distance_to_center + localz) + localx*localx + localy*localy);
  return GetMom() * GetGlobalZ()/distance;
}


Double_t STNeuLANDCluster::GetBeta()
{
  if(0 > distance) distance = TMath::Sqrt((distance_to_center + localz)*(distance_to_center + localz) + localx*localx + localy*localy);
  Double_t v = distance/tof;
  return v/c;
}

Double_t STNeuLANDCluster::GetGamma()
{
  Double_t mybeta = GetBeta();

  if(0!=mybeta){
    return 1/TMath::Sqrt(1-mybeta*mybeta);
  }
  else{
    return 0;
  }
 
}
