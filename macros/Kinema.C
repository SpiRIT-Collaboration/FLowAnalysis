TVector3 GetLorentzBoost(UInt_t isel);

void RapidityShift()
{
  TGraph *grph = new TGraph();
  grph->SetName("grph");

  TGraph *grphypt = new TGraph();
  grphypt->SetName("grphypt");

  Double_t mass_p = 1.00727646688*931.4940954;


  TVector3 p(0.,0.,356.);


  Double_t Etot = sqrt( p.Mag2() + mass_p*mass_p);
  TLorentzVector pp(p, Etot);
  
  TVector3 boost = GetLorentzBoost(4);


  UInt_t iy = 0;
  for(UInt_t i = 0; i < 180; i++) {
    Double_t theta = i*TMath::Pi()/180. ;


    pp.SetTheta(theta);
    Double_t rapidity = 0.5 * log( (Etot + pp.Z())/ (Etot - pp.Z() ) );

    pp.Boost(boost);    
    theta = pp.Theta();
    
    grph->SetPoint(i, theta*180./TMath::Pi(), rapidity);

    if( abs( rapidity ) < 0.12 ) {
      grphypt->SetPoint(iy, theta*180./TMath::Pi(), rapidity);
      iy++;

      if( abs( rapidity ) < 0.001 )
	cout << " Mid Rapidity " << theta << " " << theta*180./TMath::Pi() << endl;

      if( rapidity < -0.1 )
	cout << " Shift Rapidity " << theta << " " << theta*180./TMath::Pi() << endl;
    }

    pp.Boost(-boost);
  }

  auto c1 = new TCanvas("c1","c1");
  grph->Draw("ALP");

  auto c2 = new TCanvas("c2","c2");
  grphypt->Draw("ALP");
  
  
}

void Kinema(UInt_t isel = 0)
{
  

  //  RapidityShift();

  GetLorentzBoost(0);
  GetLorentzBoost(1);
  GetLorentzBoost(2);
  GetLorentzBoost(3);
  GetLorentzBoost(4);

  return;
}



TVector3 GetLorentzBoost(UInt_t isel = 4)
{
  if( isel > 4 ) exit(0);

  Double_t amu  = 931.4940954; //MeV/c2  
  Double_t c    = 299792458.; //m/s                                                                 

  TString system[5];
  Double_t AB[5];
  Double_t mB[5];
  Double_t eB_lb[5];
  Double_t mT[5];     

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

  // system[3] = "(208Pb + 208Pb";
  // AB[3]     =  208.;
  // mB[3]     =  208; //amu
  // eB_lb[3]  =  158000.;
  // mT[3]     =  208;

  system[4] = "(p + p)";
  AB[4]     = 1.;
  mB[4]     = 1.00727646688;
  eB_lb[4]  = 268.9;
  mT[4]     = 1.00727646688;
  


  UInt_t isys = isel;


  Double_t EkB_lb    = eB_lb[isys]  * mB[isys];
  mB[isys] *= amu;
  mT[isys] *= amu;

  // Beam      
  Double_t EB_lb  = EkB_lb + mB[isys];
  Double_t PB_lb     = sqrt(EB_lb*EB_lb - mB[isys]*mB[isys]);


  //PB_lb = sqrt(EkB_lb*EkB_lb + 2.*mB[isys]*EkB_lb);
  Double_t YB_lb  = 0.5 * log( (EB_lb + PB_lb) / (EB_lb - PB_lb) );

  // CM system     
  auto E_cm    = sqrt( mB[isys]*mB[isys] + mT[isys]*mT[isys] + 2.* mT[isys] * EB_lb);
  // Jon's 
  auto Gamm_cm = E_cm/(mB[isys] + mT[isys]);
  auto Beta_cm = TMath::Sqrt(1. -1./Gamm_cm/Gamm_cm);
  
  // be out since 27 June 2019
  // auto Beta_cm = PB_lb/(EB_lb + mT[isys]);
  // auto Gamm_cm = (EB_lb + mT[isys]) / E_cm;  

  auto y_cm = 0.5 * log( (1+Beta_cm)/ (1-Beta_cm) );

  auto EB_cm = Gamm_cm * EB_lb - Gamm_cm*Beta_cm * PB_lb;
  auto PB_cm = -Gamm_cm* Beta_cm * EB_lb + Gamm_cm * PB_lb;
  auto YB_cm =  0.5 * log( (EB_cm + PB_cm) / (EB_cm - PB_cm) );

  auto ET_cm = Gamm_cm * mT[isys];
  auto PT_cm = -Gamm_cm * Beta_cm * mT[isys];
  auto YT_cm =  0.5 * log( (ET_cm + PT_cm) / (ET_cm - PT_cm) );

  cout << " ------------------------ " << endl;

  cout << " --- Reaction  " << system[isys] << endl;
   //      << " Beam Mass   " << mB[isys]/amu << " [u] "
   //      << " Beam Enery  " << eB_lb[isys] << " [MeV/u] "
   //      << " Target Mass " << mT[isys]/amu
   //      << " --- "
   //      << endl;

   // cout << " Beam --- " 
   //      << " Ekine  " << EkB_lb << " [MeV] "
   //      << " Energy " << EB_lb  << " [MeV] "
   //      << " P      " << PB_lb  << " [MeV/c] "
   //      << " Y      " << YB_lb  
   //      << endl;

  cout << " Ecm  --- " << E_cm << " [MeV] "
       << " Beta_cm " << Beta_cm 
       << " Gamm_cm " << Gamm_cm 
       << " Y_cm " << y_cm
       << endl;


  cout << " E_beam   " << EB_cm
    //      << " P_beam   " << PB_cm
       << " Y_beam   " << YB_cm
       << endl;
  cout << " E_target " << ET_cm
   //      << " P_target " << PT_cm
       << " Y_target " << YT_cm 
       << endl;

  // cout << " dY " << YB_cm - YT_cm << endl;

  auto bmVec = new TLorentzVector( TVector3(0., 0., PB_lb), EB_lb );
  auto tgVec = new TLorentzVector( TVector3(0., 0., 0.), mT[isys] );

  auto totalVec = *bmVec + *tgVec;
  TVector3 boostVec = totalVec.BoostVector();

  cout << " Beta_cm "     << totalVec.Beta(); 
  cout << " Gamma_cm "    << totalVec.Gamma(); 
  cout << " Rapidity_cm " << totalVec.Rapidity(); 
  cout << " boostVec Z "     << boostVec.Z(); 
  cout << endl;

  // cout << " ---> before " << endl;
  // cout << " beam g : " << bmVec->Gamma()
  //      << "      Y : " << bmVec->Y()
  //      << "      Z : " << bmVec->Z()
  //      << "      E : " << bmVec->E()
  //      << endl;

  // cout << " Targ g : " << tgVec->Gamma()
  //      << "      Y : " << tgVec->Y()
  //      << "      Z : " << tgVec->Z()
  //      << "      E : " << tgVec->E()
  //      << endl;

  bmVec->Boost(-boostVec);
  tgVec->Boost(-boostVec);

  // cout << " ---> after " << endl;
  // cout << " beam g : " << bmVec->Gamma()
  //      << "      Y : " << bmVec->Y()
  //      << "      Z : " << bmVec->Z()
  //      << "      E : " << bmVec->E()
  //      << endl;

  // cout << " Targ g : " << tgVec->Gamma()
  //      << "      Y : " << tgVec->Y()
  //      << "      Z : " << tgVec->Z()
  //      << "      E : " << tgVec->E()
  //      << endl;

  //  cout << " Ecm " << (*bmVec + *tgVec).E() << endl;

  // cout << " Beam   Rapidity " << 0.5 * log( ( bmVec->E() + bmVec->Z())/(bmVec->E() - bmVec->Z()) ) ;
  // cout << " Target Rapidity " << 0.5 * log( ( tgVec->E() + tgVec->Z())/(tgVec->E() - tgVec->Z()) );

  cout << " y_beam   = " << bmVec->Rapidity()<< endl;
  cout << " y_target = " << tgVec->Rapidity()<<endl;




  cout << " =====-------=====---------- " << endl;
  
  return boostVec;


  Double_t mass = mT[4] * amu;

  Double_t initR = 0;

  for(UInt_t i = 0; i < 3; i++) {
  
    Double_t p    = 423. - i*50.;

    auto Etot = sqrt( mass*mass + p*p );
  
    auto Rapp = 0.5 * log(  (Etot + p) / (Etot - p) );

    auto Beta = p / mass;


    if( i == 0 )
      initR = Rapp;


    cout << " mass " << mass 
	 << " p " << p << endl;
    cout << " Etot " << Etot << endl;
    cout << " beta " << Beta << endl;
    cout << " TOF  " << Beta*88.60 << endl;
    cout << "rapditiy " << Rapp << endl;
    cout << " rap diff " << Rapp - initR
	 << " P diff " << 373 - p
	 << endl;
  }

}

