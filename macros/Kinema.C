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
  
  TGraph *grph = new TGraph();
  //  RapidityShift();
  TVector3 vec;

  vec = GetLorentzBoost(4);
  TLine *ln = new TLine(100., vec.Z(), 135., vec.Z());

  vec = GetLorentzBoost(1);
  grph->SetPoint(0, 108., vec.Z() );

  vec = GetLorentzBoost(3);
  grph->SetPoint(1, 112., vec.Z() );

  vec = GetLorentzBoost(2);
  grph->SetPoint(2, 124., vec.Z() );


  vec = GetLorentzBoost(0);
  grph->SetPoint(3, 132., vec.Z() );


  TCanvas *cc;
  grph->SetMarkerStyle(20);
  grph->SetMarkerSize(0.4);
  grph->Draw("AP");
  ln->Draw();




  return;
}



TVector3 GetLorentzBoost(UInt_t isel = 4)
{
  if( isel >= 6 ) exit(0);

  Double_t amu  = 931.4940954; //MeV/c2  
  Double_t c    = 299792458.; //m/s                                                                 

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

  TVector3 nnBoost(0.,0., 0.355151);

  system[5] = "(112Sn + XSn)";
  AB[5]     = 112.;
  mB[5]     = 111.8773895;
  eB_lb[5]  = 270.2;
  mT[5]     = 123.8773895 * 112./124.;;

  //Sn-108    107.911892833   +- 0.000005900    0+                        10.30M 8 
  //Sn-112    111.904821807   +- 0.000001813    0+               0.97     stable  
  //Sn-124    123.905273581   +- 0.000001492    0+               5.79     stable 
  //Sn-132    131.917821719   +- 0.000006866    0+                        39.7S 8    

  // mB[0]     = 131.917821719;
  // mB[1]     = 107.911892833;
  // mB[2]     = 123.905273581;
  // mB[3]     = 111.904821807;

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

  cout << " Beam   Boost " << bmVec->BoostVector().Z() << endl;

  bmVec->Boost(-boostVec);
  tgVec->Boost(-boostVec);


  cout << " y_beam   = " << bmVec->Rapidity()<< endl;
  cout << " y_target = " << tgVec->Rapidity()<<endl;


  bmVec->Boost(boostVec-nnBoost);
  cout << " y_bam - y_nn = " << bmVec->Rapidity()<<endl;


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

