{

  const UInt_t th_bin = 4;

  Float_t mass  =   938.2720813; // 1875.612762; // 
  Float_t maxmom   = 900.;
  TLorentzVector lv;

  // lv.SetXYZM(0.,0.,maxmom,mass);
  // lv.SetTheta(TMath::Pi()/4.);



  TMultiGraph *mgr = new TMultiGraph("ytheta",";Rapidity; Pt");
  const UInt_t mom_bin = 10;
  TGraph grph[th_bin];

  for(UInt_t j = 0; j < mom_bin; j++ )
    {
      Float_t mom = (j+1)*maxmom/(Double_t)mom_bin; 
      lv.SetXYZM(0.,0.,mom, mass);

      for(UInt_t i = 0; i < th_bin; i++ ) 
	{
	  Float_t theta = (i+1)*TMath::Pi()/(2.*(Double_t)th_bin);
	  lv.SetTheta(theta);
	  std::cout << " theta " << theta*TMath::RadToDeg() 
		    << " mass " << mass << " mom " << mom << std::endl;
	  grph[i].SetPoint(j, lv.Rapidity(), lv.Pt());
	  
	}
    }

  for(UInt_t j = 0; j < th_bin; j++)
    {
      grph[j].SetName(Form("gr_%d",j));
      //      std::cout << "---" << j << std::endl;
      //grph[j].Print();
      
      mgr->Add(&grph[j],"lp");
    }

  mgr->Draw("ALP");

}
