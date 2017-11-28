Double_t *vMean(vector<Double_t> &vec)
{
  //  cout << "vMean" << vec.size() << endl;

  vector<Double_t>::iterator ibgn = vec.begin();
  vector<Double_t>::iterator iend = vec.end();

  auto *vc = new Double_t[2];
  vc[0]  =  TMath::Mean(ibgn, iend);
  vc[1]  =  TMath::RMS(ibgn, iend);

  if(vec.size() < 2) {
    vc[0] = -999.;
    vc[1] = 0.;
  }

  cout << " vc[0] " << vc[0] << endl; 
  return vc;
}
Double_t *vn(UInt_t hm, vector<Double_t> &vphi)
{

  vector<Double_t> fvCos;
  
  //  cout << " vphi size " << vphi.size() << endl;
  Double_t findx = (Double_t)hm;
  for(UInt_t i = 0; i < (UInt_t)vphi.size(); i++)
    fvCos.push_back(cos(findx * vphi.at(i)));

  vector<Double_t>::iterator ibgn = fvCos.begin();
  vector<Double_t>::iterator iend = fvCos.end();

  Double_t *vcos = new Double_t[2];
  vcos[0]       =  TMath::Mean(ibgn, iend);
  Double_t rms  =  TMath::RMS(ibgn, iend);

  rms /= sqrt(vphi.size());
  vcos[1] = rms;

  if(vphi.size() < 2) {
    vcos[0] = -999.;
    vcos[1] = 0.;
  }

  //  cout << hm << " Mean " << vcos[0] << " RMS " << vcos[1] << endl; 
  return vcos;
}


