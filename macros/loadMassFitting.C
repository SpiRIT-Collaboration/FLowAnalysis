void loadMassFitting() {


  auto massFitter = new STMassFunction();

  auto fstatus = massFitter->SetMassFitFunction("/cache/scr/spirit/mizuki/Bethe-Bloch_Fitter/mk_NewFitter_20190111/",
					   "MassFitter.root",
					   "f1IterMassFitRotate_132Sn");   

  massFitter->SetMassRegion();


}
