{
  gStyle->SetOptStat(0);

  TString partName ;
  partName = "proton";
  partName = "3He";
  partName = "deuteron";
  partName = "triton";
  partName = "4He";

  auto c1 = new TCanvas("c1","c1",500,1000);
  c1->Divide(1,3);

  // TString fversion[][2] = { {".v50.2.19.","|phi| >150&&M 0-30"}, 
  // 			    {".v50.2.20.","|phi| >150&&M 30-40"}, 
  // 			    {".v50.2.21.","|phi| >150&&M 40-50"}, 
  // 			    {".v50.2.22.","|phi| >150&&M 50-80"}};

  TString fversion[][3] = { {"132",".v52.10.16.","132Sn"}, 
			    {"112",".v52.10.16.","112Sn"},
			    {"108",".v52.10.16.","108Sn"}};

  TDirectory *gDir = new TDirectory();
  TFile *_ffin;

  for(UInt_t i = 0; i < 3; i++ ) {
    _ffin = TFile::Open("data/finYPt_"+fversion[i][0]+"Sn_"+ partName + fversion[i][1] +"root");
    c1->cd(i+1);
    TString hname = Form("hyptacp_%d",i );
    hyptacp->SetName( hname );
    hyptacp->SetTitle(fversion[i][2]+":"+partName);
    hyptacp->Draw("colz");
    hyptacp->SetDirectory( gDir );
    _ffin->Close();
  }

}
