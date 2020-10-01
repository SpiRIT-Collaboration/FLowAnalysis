TCanvas* c1;
TFile* _ffin;
TH2D*  hyptacp;

void PlotAcc(UInt_t ipart=0){
  gStyle->SetOptStat(0);

  TString partName[] = {"proton","deuteron","triton","3He","4He"};

  c1 = new TCanvas("c1","c1",500,1000);
  c1->Divide(1,3);
  hyptacp = new TH2D();

  TString fversionC[][3] = { {"132",".v52.10.28.","132Sn"}, 
			    {"112",".v52.10.28.","112Sn"},
			    {"108",".v52.10.28.","108Sn"}};

  TString fversion[][3] = { {"132",".v52.10.16.","132Sn"}, 
			    {"112",".v52.10.16.","112Sn"},
			    {"108",".v52.10.16.","108Sn"}};


  for(UInt_t i = 0; i < 3; i++ ) {
    _ffin = TFile::Open("data/finYPt_"+fversion[i][0]+"Sn_"+ partName[ipart] + fversion[i][1] +"root");
    TString hname = "hyptacp" ;
    hyptacp = (TH2D*)_ffin->Get(hname);

    hname = Form("hyptacp_%d",i );
    hyptacp->SetName( hname );
    hyptacp->SetTitle(fversion[i][2]+":"+partName[ipart]);

    c1->cd(i+1);
    hyptacp->Draw("colz");
    hyptacp->SetDirectory( 0 );
    _ffin->Close();
  }

  TString sName = "acc_"+partName[ipart]+fversion[0][1]+"png";
  c1->SaveAs(sName);

}
