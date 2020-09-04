
auto fv1fitorg = new TF1("fv1fitorg","[0]+[1]*x+[2]*x^3",-1.3,1.3); //-0.5,1.3);
auto fv1fit    = new TF1("fv1fit","[1]*(x-[3])+[2]*(x-[3])^3",-1.3,1.3); //-0.5,1.3);
auto fv2fit    = new TF1("fv2fit","[0]+[1]*(x-[2])^2",-1.3,1.3);

auto lslope    = new TF1("lslope","[0]+[1]*x",-1., 1.);

TF1 *fv1y;
TF1 *fv2y;
void FlowFunction()
{
  //---- Flow parameter -----
  //-----------------------------------------
  fv1y = new TF1("fv1y","[0]+[1]*x+[2]*x^3"  ,-1.,1.);
  fv1y->SetParameter(0,0);
  fv1y->SetParameter(1,5.18056e-01);
  fv1y->SetParameter(1,5.18056e-01);
  fv1y->SetParameter(2,-1.84025e-01);

  fv2y = new TF1("fv2y","[0]+[1]*x^2+[2]*x^4",-1.,1.);
  fv2y->SetParameter(0, -0.08);
  fv2y->SetParameter(1,  0.1);
  fv2y->SetParameter(2, -0.02);

}

