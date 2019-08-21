{
  

  gSystem->Load("libMathMore");


  auto mbs1 = new TF1("mbs1","sqrt(TMath::Pi())/(2.*sqrt(2))*x*exp(-x*x/4.)*(ROOT::Math::cyl_bessel_i(0,x*x/4.)+ROOT::Math::cyl_bessel_i(1,x*x/4.))",0,3);
  auto mbs2 = new TF1("mbs2","sqrt(TMath::Pi())/(2.*sqrt(2))*x*exp(-x*x/4.)*(ROOT::Math::cyl_bessel_i(0.5,x*x/4.)+ROOT::Math::cyl_bessel_i(1.5,x*x/4.))",0,3);
  auto mbs3 = new TF1("mbs3","sqrt(TMath::Pi())/(2.*sqrt(2))*x*exp(-x*x/4.)*(ROOT::Math::cyl_bessel_i(1,x*x/4.)+ROOT::Math::cyl_bessel_i(2.,x*x/4.))",0,3);

  
  auto mbs4 = new TF1("mbs4","0.626657*x-0.09694*pow(x,3) + 0.02754 * pow(x,4)-0.002283*pow(x,5)",0,3);
  auto mbs5 = new TF1("mbs5","0.25*pow(x,2)-0.011414*pow(x,3)-0.034726*pow(x,4)+0.006815*pow(x,5)",0,3);

  mbs1->SetMaximum(1.05);
  mbs1->Draw("LP");
  mbs2->Draw("same");
  mbs3->Draw("same");

  mbs4->SetLineColor(4);
  mbs4->Draw("same");
  mbs5->SetLineColor(4);
  mbs5->Draw("same");


  Double_t x = 1.5;
  cout << "mbs 1" << sqrt(TMath::Pi())/(2.*sqrt(2))*x*exp(-x*x/4.)*(ROOT::Math::cyl_bessel_i(0,x*x/4.)+ROOT::Math::cyl_bessel_i(1,x*x/4.)) << endl;
  cout << "mbs 4" << 0.626657*x-0.09694*pow(x,3) + 0.02754 * pow(x,4)-0.002283*pow(x,5) << endl;

}
