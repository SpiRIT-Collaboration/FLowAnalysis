
void FOPI_v2_energy()
{
  
  Double_t FOPI_v2new[][2] = {
    {0.09009744708060832, 0.07133565040405002},
    {0.1200079863689576, 0.022777667577971186},
    {0.14960577368289288, -0.0048473349374591},
    {0.25067600880901303, -0.0541932367983027},
    {0.39985520708461514, -0.0706093666002593},
    {0.5968201457318721, -0.07402138259552166},
    {0.8026823561286831, -0.07154907389040771},
    {1.005707917589353, -0.06583171627095352 },
    {1.2086168169606657, -0.05910243694894107},
    {1.5083468193747376, -0.04484675998194637}};




  Double_t FOPI_v2and05[][2] = {};

  auto grh_1 = new TGraph();
  auto grh_2 = new TGraph();

  UInt_t siz = (UInt_t)sizeof(FOPI_v2new)/sizeof(Double_t)/2.;
  for( UInt_t i = 0; i < siz; i++) 
    grh_2->SetPoint( i, FOPI_v2new[i][0], FOPI_v2new[i][1] );
  

  grh_2->SetLineColor(4);
  grh_2->SetMarkerStyle(25);
  grh_2->SetMarkerColor(4);


  auto mgrh = new TMultiGraph();

  mgrh->Add(grh_2, "p");

  auto cc1 = new TGraph("cc1","cc1");
  mgrh->Draw("ALP");
  cc1->SetLogx();
}

