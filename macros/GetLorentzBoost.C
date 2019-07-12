TVector3 LorentzBoost(UInt_t m)
{
  Double_t amu  = 931.4940954; //MeV/c2  
  
  TString system[6];
  Double_t AB[6];
  Double_t mB[6];
  Double_t eB_lb[6];
  Double_t mT[6];     

  system[0] =  "(132Sn + 124Sn)";
  AB[0]     =  132.;
  mB[0]     =  131.8906 ; //amu                                                                 
  eB_lb[0]  =  268.9 ;    //MeV/amu; incident energy    
  mT[0]     =  123.8773895;   //amu Target mass

  system[1] = "(108Sn + 112Sn)";
  AB[1]     =  108.;
  mB[1]     =  107.8844964; //amu
  eB_lb[1]  =  268.9;
  mT[1]     =  111.8773895;;

  system[2] = "(124Sn + 112Sn)";
  AB[2]     =  124.;
  mB[2]     =  123.8778449; //amu
  eB_lb[2]  =  270.2;
  mT[2]     =  111.8773895;;

  system[3] = "(112Sn + 124Sn)";
  AB[3]     =  112.;
  mB[3]     =  111.8773895; //amu
  eB_lb[3]  =  270.2;
  mT[3]     =  123.8773895;

  system[4] = "(p + p)";
  AB[4]     = 1.;
  mB[4]     = 1.00727646688;
  eB_lb[4]  = 268.9;
  mT[4]     = 1.00727646688;

  system[5] = "(100Sn + 100Sn)";
  AB[5]     = 108.;
  mB[5]     = 105.;
  eB_lb[5]  = 268.9;
  mT[5]     = 108.8;

  


  UInt_t sysid = 4;
  if( m < 5)
    sysid = m;


  Double_t EkB_lb    = eB_lb[sysid]  * mB[sysid];
  mB[sysid] *= amu;
  mT[sysid] *= amu;

  // Beam      
  Double_t EB_lb  = EkB_lb + mB[sysid];
  Double_t PB_lb     = sqrt(EB_lb*EB_lb - mB[sysid]*mB[sysid]);


  //PB_lb = sqrt(EkB_lb*EkB_lb + 2.*mB[sysid]*EkB_lb);
  Double_t YB_lb  = 0.5 * log( (EB_lb + PB_lb) / (EB_lb - PB_lb) );

  // CM system     
  auto E_cm    = sqrt( mB[sysid]*mB[sysid] + mT[sysid]*mT[sysid] + 2.* mT[sysid] * EB_lb);
  auto Beta_cm = PB_lb/(EB_lb + mT[sysid]);
  auto Gamm_cm = (EB_lb + mT[sysid]) / E_cm;

  auto y_cm = 0.5 * log( (1+Beta_cm)/ (1-Beta_cm) );


  auto EB_cm = Gamm_cm * EB_lb - Gamm_cm*Beta_cm * PB_lb;
  auto PB_cm = -Gamm_cm* Beta_cm * EB_lb + Gamm_cm * PB_lb;
  auto YB_cm =  0.5 * log( (EB_cm + PB_cm) / (EB_cm - PB_cm) );

  auto ET_cm = Gamm_cm * mT[sysid];
  auto PT_cm = -Gamm_cm * Beta_cm * mT[sysid];
  auto YT_cm =  0.5 * log( (ET_cm + PT_cm) / (ET_cm - PT_cm) );


  cout << " --- Reaction  " << system[sysid] << endl
       << " Beam Mass   " << mB[sysid]/amu << " [u] "
       << " Beam Enery  " << eB_lb[sysid] << " [MeV/u] "
       << " Target Mass " << mT[sysid]/amu
       << " --- "
       << endl;

  cout << " Beam --- " 
       << " Ekine  " << EkB_lb << " [MeV] "
       << " Energy " << EB_lb  << " [MeV] "
       << " P      " << PB_lb  << " [MeV/c] "
       << " Y      " << YB_lb  
       << endl;

  cout << " Ecm  --- " << E_cm << " [MeV] "
       << " Beta_cm " << Beta_cm 
       << " Gamm_cm " << Gamm_cm 
       << " Y_cm " << y_cm
       << endl;


  cout << " E_beam   " << EB_cm
       << " P_beam   " << PB_cm
       << " Y_beam   " << YB_cm
       << endl;
  cout << " E_target " << ET_cm
       << " P_target " << PT_cm
       << " Y_target " << YT_cm 
       << endl;

  cout << " dY " << YB_cm - YT_cm << endl;

  cout << " ------------------------ " << endl;
  auto bmVec = new TLorentzVector( TVector3(0., 0., PB_lb), EB_lb );
  auto tgVec = new TLorentzVector( TVector3(0., 0., 0.), mT[sysid] );
  auto totalVec = *bmVec + *tgVec;
  TVector3 boostVec = totalVec.BoostVector();

  cout << " before " << endl;
  cout << " beam g : " << bmVec->Gamma()
       << "      Y : " << bmVec->Y()
       << "      Z : " << bmVec->Z()
       << "      E : " << bmVec->E()
       << endl;

  cout << " Targ g : " << tgVec->Gamma()
       << "      Y : " << tgVec->Y()
       << "      Z : " << tgVec->Z()
       << "      E : " << tgVec->E()
       << endl;

  bmVec->Boost(-boostVec);
  tgVec->Boost(-boostVec);

  cout << " after " << endl;
  cout << " beam g : " << bmVec->Gamma()
       << "      Y : " << bmVec->Y()
       << "      Z : " << bmVec->Z()
       << "      E : " << bmVec->E()
       << endl;

  cout << " Targ g : " << tgVec->Gamma()
       << "      Y : " << tgVec->Y()
       << "      Z : " << tgVec->Z()
       << "      E : " << tgVec->E()
       << endl;

  cout << " Ecm " << (*bmVec + *tgVec).E() << endl;

  cout << " Beam   Rapidity " << 0.5 * log( ( bmVec->E() + bmVec->Z())/(bmVec->E() - bmVec->Z()) )  << endl;
  cout << " Target Rapidity " << 0.5 * log( ( tgVec->E() + tgVec->Z())/(tgVec->E() - tgVec->Z()) )  << endl;

  cout << " beam   Beta " << bmVec->Beta() << " gammma " << bmVec->Gamma() << endl;
  cout << " target Beta " << tgVec->Beta() << " gammma " << tgVec->Gamma() << endl;

  cout << " =====-------=====---------- " << endl;


  return boostVec;  
}


Double_t GetRapidity_cm(TVector3 P, Double_t mass, TVector3 bvec)
{

  Double_t Etot = sqrt( P.Mag2() + pow(mass,2) );
  TLorentzVector LrnzVec( P, Etot);

  LrnzVec.Boost(bvec);

  Double_t PZcm = LrnzVec.Z();
  Double_t Ecm  = LrnzVec.E();
  //  auto rapidity = 0.5*log( (Ecm + PZcm)/(Ecm - PZcm) );
  auto rapidity = LrnzVec.Rapidity();

  return rapidity;
}

void GetLorentzBoost(UInt_t m = 0)
{
  LorentzBoost(m);
}
