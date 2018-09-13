#include "STNeuLANDCluster.hh"
#include <TMath.h>
#include <iostream>

ClassImp(STNeuLANDCluster)

STNeuLANDCluster::STNeuLANDCluster()
{
  Clear();
  angle = 29.579*TMath::Pi()/180. ; // deg.

  target_offset = -8.9; //mm
  targetPos = TVector3(0,0,target_offset);

  distance_to_center =  8561.05 + 280; //mm Mizuki's photogrametry to the blue frame + to the 1st plane surface
  
  ncglobal_offset = TVector3(distance_to_center*TMath::Sin(angle), 0., 
			   distance_to_center*TMath::Cos(angle));

  c = 299.79245; //mm/ns

  //  tof_offset = 37.5; //ns temporaly.
  tof_offset = 0.; // after 1 Feb. 2018
}


void STNeuLANDCluster::SetLocalPos(TVector3 posv)
{
  nclocalPos = 10.*posv;
  nclocalx = nclocalPos.X();
  nclocaly = nclocalPos.Y();
  nclocalz = nclocalPos.Z();

  ncglobalPos = nclocalPos;
  ncglobalPos.RotateY(angle);
  ncglobalPos += ncglobal_offset;
  //  ncglobalPos -= targetPos;

  //  globalPos.SetX(distance_to_center*TMath::Sin(angle) +
  // 	      localx*TMath::Cos(angle) + localz*TMath::Sin(angle) );
  // globalPos.SetZ(distance_to_center*TMath::Cos(angle) -
  // 	      localx*TMath::Sin(angle) + localz*TMath::Cos(angle) + target_offset);
  // globalPos.SetY(localy);


  ncdistance = ncglobalPos.Mag();

}

void STNeuLANDCluster::SetLocalPosLast(TVector3 posv)
{
  nclocalPos_last = 10.*posv;

  nclocalx_last = nclocalPos_last.X();
  nclocaly_last = nclocalPos_last.Y();
  nclocalz_last = nclocalPos_last.Z();

  ncglobalPos_last = nclocalPos_last;
  ncglobalPos_last.RotateY(angle);
  ncglobalPos_last += ncglobal_offset;
  //  ncglobalPos_last += targetPos;
}

void STNeuLANDCluster::SetMass(TString v)
{
  if(v == "neutron")
    SetMass(2112);

  else if(v == "proton")
    SetMass(2212);
  
  else if(v == "deuteron")
    SetMass(1000010020);
  
  else if(v == "triton")
    SetMass(1000010030);

  else
    SetMass(0);
}

void STNeuLANDCluster::SetMass(UInt_t v)
{
  // neutron = "fPID==2112";
  // proton  = "fPID==2212";
  // deuteron= "fPID==1000010020";
  // triton  = "fPID==1000010030";
  
  ncPID = v;
 
  auto mp   =    938.2720813;
  auto mn   =    939.565346;
  auto md   =   1875.612762;
  auto mt   =   2808.921112;
 

  Double_t mass = 0;
  switch (ncPID) {
  case 2112:
    mass = mn;
    break;

  case 2212:
    mass = mp;
    break;

  case 1000010020:
    mass = md;
    break;

  case 1000010030:
    mass = mt;
    break;

  default:
    mass = 0;
    break;
  }

  ncmass = mass;
  SetMomentum();
}


void STNeuLANDCluster::SetMomentum() 
{

  nctof -= tof_offset;

  ncbeta  = (ncdistance/nctof)/c;

  if(ncbeta < 1.)
    ncgamma = 1/TMath::Sqrt(1. - ncbeta*ncbeta);
  else
    ncgamma = 0.;
  
  if(nctof < 5){ 
    ncPID = 22;
    ncmass = 0.;
  }

  if(ncmass*ncbeta*ncgamma > 0){
    ncP = ncglobalPos;

    ncP.SetMag(ncmass*ncbeta*ncgamma);

    ncE = ncmass * ncgamma;

    Double_t p_para = ncP.Mag()*TMath::Cos(ncP.Theta());
    ncRapidity = 0.5 * log( (ncE + p_para)/(ncE - p_para) );
  }
  else
    ncPID = 0;

}

void STNeuLANDCluster::SetBeamAngle(Double_t va, Double_t vb)
{
  ncbeamAngleA = va; 
  ncbeamAngleB = vb;

  if(ncmass*ncbeta*ncgamma > 0){

    ncP.RotateY(-ncbeamAngleA);
    ncP.RotateX(-ncbeamAngleB);

    Double_t p_para = ncP.Mag()*TMath::Cos(ncP.Theta());
    ncRapidity = 0.5 * log( (ncE + p_para)/(ncE - p_para) );
  }
}

void STNeuLANDCluster::Clear(Option_t *)
{
  ncnhit = -1;
  ncedep = -1;
  nctof = -1;
  nctof_last = -1;

  nclocalx = -500; nclocaly = -500; nclocalz = -1;
  nclocalx_last = -500; nclocaly_last = -500; nclocalz_last = -1;
  
  ncdistance = -1;
  nclocalPos.SetX(nclocalx);
  nclocalPos.SetY(nclocaly);
  nclocalPos.SetZ(nclocalz);

  nclocalPos_last.SetX(nclocalx_last);
  nclocalPos_last.SetY(nclocaly_last);
  nclocalPos_last.SetZ(nclocalz_last);

  ncglobalPos.SetX(nclocalx);
  ncglobalPos.SetY(nclocaly);
  ncglobalPos.SetZ(nclocalz);

  ncvetot_all = -1;
  ncvetoq_all = -1;
  ncvetot_bar = -1;
  ncvetoq_bar = -1;
  ncvetot_mid = -1;
  ncvetoq_mid = -1;
  ncvetot_loose = -1;
  ncvetoq_loose = -1;
  ncvetot_loose1 = -1;
  ncvetoq_loose1 = -1;
  ncvetot_loose2 = -1;
  ncvetoq_loose2 = -1;

  ncP = TVector3(0.,0.,1.);
  ncbeta  = 0.;
  ncgamma = 0.;
  ncmass = 0.;

  ncE = 0.;
  ncRapidity = 0.;

  ncbeamAngleA = 0.;
  ncbeamAngleB = 0.;
}
