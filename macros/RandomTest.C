void RandomTest(UInt_t is = 123)
{
  TF1 *rnd = new TF1("rnd","1+2.*[0]*cos(x-[2]) + 2.*[1]*cos(2.* (x-[2]))",-3.14, 3.14);
  rnd->SetParameter(0, 0.2);
  rnd->SetParameter(1, 0.);
  rnd->SetParameter(2, 0.);
  TString sget;

  while(1){
    TRandom pp(0);
    gRandom->SetSeed(is);
    
    cout << rnd->GetRandom() << endl;

    cin >> sget;
    if( sget == "y" ) break;


  }


}
