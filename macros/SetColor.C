
void SetColor()
{
  TColor *fcol;
  
  fcol = gROOT->GetColor(3);
  fcol -> SetRGB(0., 0., 1.);

  fcol = gROOT->GetColor(4);
  fcol -> SetRGB(0., 256/256., 0.);

  fcol = gROOT->GetColor(5);
  fcol -> SetRGB(256./256., 215./256., 0.);

  


  UInt_t i = 11;
  fcol =  gROOT->GetColor(i); i++;
  fcol -> SetRGB(255./256, 165./355., 0./256.);

  fcol =  gROOT->GetColor(i); i++;
  fcol -> SetRGB(34./256, 139./355., 134./256.);

  fcol =  gROOT->GetColor(i); i++;
  fcol -> SetRGB(128./256, 128./355., 15./256.);

  fcol =  gROOT->GetColor(i); i++;
  fcol -> SetRGB(245./256, 153./355., 153./256.);

  fcol =  gROOT->GetColor(i); i++;
  fcol -> SetRGB(204./256, 102./355., 255./256.);

  fcol =  gROOT->GetColor(i); i++;
  fcol -> SetRGB(153./256, 51./355., 51./256.);

  fcol =  gROOT->GetColor(i); i++;
  fcol -> SetRGB(102./256, 204./355., 102./256.);

  fcol =  gROOT->GetColor(i); i++;
  fcol -> SetRGB(102./256, 51./355., 256./256.);

  fcol =  gROOT->GetColor(i); i++;
  fcol -> SetRGB(128./256, 128./355., 128./256.);



  // TCanvas *c = new TCanvas();
  // c->DrawColorTable();

}
