#include "STBetheBlochFittingFunction.hh"


Double_t BetheBlochFitter::fddedx(Double_t *x, Double_t *p)
{
  auto mass = x[0];
  // p[2] = mom
  // p[3] = Z
  // p[4] = dedx
  auto ZA = 17.2/37.6;
  auto I2 = 0.9*I_Ar + 0.1*I_CH4; I2 = I2*I2;
  auto pZ = p[2]*p[3];
  auto b2 = pZ*pZ/(pZ*pZ+mass*mass);
  auto g2 = 1./(1.-b2);
  auto Wx = 2*kme*b2*g2*mass*mass/((kme+mass)*(kme+mass)+2*kme*mass*(TMath::Sqrt(g2)-1));
  auto ddedx = p[4]-p[0]*ZA*p[3]*p[3]/b2*(0.5*TMath::Log(2*kme*b2*g2*Wx/I2)-b2-p[1]);
  return ddedx;
}

Double_t BetheBlochFitter::fdedx(Double_t Z, Double_t m, Double_t *x, Double_t *p)
{
  auto ZA = 17.2/37.6;
  auto I2 = 0.9*I_Ar + 0.1*I_CH4; I2 = I2*I2;
  auto pZ = x[0]*Z;
  auto b2 = pZ*pZ/(pZ*pZ+m*m);
  auto g2 = 1./(1.-b2);
  auto Wx = 2*kme*b2*g2*m*m/((kme+m)*(kme+m)+2*kme*m*(TMath::Sqrt(g2)-1));
  auto dedx = p[0]*ZA*Z*Z/b2*(0.5*TMath::Log(2*kme*b2*g2*Wx/I2)-b2-p[1]);
  return dedx;
}

Double_t BetheBlochFitter::CalcMass(Double_t *p)
{
  TF1 f("f",fddedx,-500.,5000.,5,1);
  f.SetParameters(p);
  ROOT::Math::RootFinder finderp(ROOT::Math::RootFinder::kBRENT);
  ROOT::Math::WrappedTF1 funcp(f);
  finderp.SetFunction(funcp, 0.1, 10000.);
  finderp.Solve();
  return finderp.Root();
}

