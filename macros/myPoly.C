#include "DoRPRes.C"

void myPoly()
{

  double pp1[] = {0., 0.626657,    0.,   -0.09694,  0.02754,  -0.002283};
  double pp2[] = {0.,       0.,  0.25,  -0.011414, -0.034726,  0.006815};

  pp1[0] = -0.415801;


  Complex xx[5];
  Complex aa[6];
  for(UInt_t i = 0; i < 6; i++ )
    aa[i] = tocomplex( pp1[5-i], 0.);

  dka(aa, xx, 5, 0.001, 100);

  std::vector<Double_t> xsolve;

  for(UInt_t i = 0; i < 5; i++ ) {
    cout << " xm[" << i+1 << "] " << xx[i].r << " + " << xx[i].i << "i" << endl;

    if( abs(xx[i].i) < 0.01 && xx[i].r > 0 && xx[i].r < 3.)
      {
	//	cout << "--->  xm[" << i+1 << "] " << xx[i].r << " + " << xx[i].i << "i" << endl;
	xsolve.push_back( xx[i].r );
      }
  }

  auto res = new Double_t[4];
  for(UInt_t k = 0; k < (UInt_t)xsolve.size(); k++ ) {
    Double_t tt = sqrt(2. * xsolve.at(k) ) ;
    
    res = GetRPResolutionwChi(res, tt, 0., 2);
  }

  cout << " res " << res[0] << endl;
}


