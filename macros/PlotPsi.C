{

    //  TCutCUTntrack[4]>=0 && ntrack[4]< 20
    // X 21424.1 mean : -0.30243 sig : 1.81192
    // Y 36257.8 mean : 0.0885154 sig : 1.07187
    // OBJ: TCutCUTntrack[4]>=20 && ntrack[4]< 40
    // X 75937 mean : -0.311857 sig : 2.71388
    // Y 129269 mean : 0.0515195 sig : 1.59689
    // OBJ: TCutCUTntrack[4]>=40 && ntrack[4]< 60
    // X 8702.31 mean : -0.277663 sig : 3.30879
    // Y 14671.6 mean : -0.0376367 sig : 1.96842

  UInt_t icolor[] = {2, 3, 4, 7, 6};
  Double_t meanX[]={-0.30243, -0.311857, -0.277663, -3.19874e-01};
  Double_t sigX[] ={1.81192, 2.71388, 3.30879,    2.67219e+00};
  Double_t meanY[]={0.0885154, 0.0515195, -0.0376367, 4.06381e-02 };
  Double_t sigY[] ={1.07187, 1.59689, 1.96842, 1.55774e+00};


  auto rChain0 = (TChain*)gROOT->FindObject("rChain0");

  TVector3 *unitP_rot = NULL;
  TBranch  *bunitP_rot;

  rChain0->SetBranchAddress("ntrack"  , ntrack);
  rChain0->SetBranchAddress("unitP_rot"  ,&unitP_rot,&bunitP_rot);
  

  TH1D *htf[4];
  for(UInt_t i = 0; i < 4; i++ )
    htf[i] = new TH1D(Form("htf_%d",i),";#Psi; 1/N dN/d#Psi", 100, -3.2, 3.2);

  auto aleg2 = new TLegend(0.24,0.7,0.75,0.98,"");


  UInt_t nevt = rChain0->GetEntries();
  for(UInt_t i = 0; i < nevt; i++){

    rChain0->GetEntry(i);

    Double_t vec_cc[5];

    if( ntrack[4] >= 0 && ntrack[4] < 20 ) {
      auto tvec = TVector2( (unitP_rot->X()-meanX[0])/sigX[0], (unitP_rot->Y()-meanY[0])/sigY[0] );
      vec_cc[0] = TVector2::Phi_mpi_pi(unitP_rot->Phi()); //tvec.Phi());
      htf[0]->Fill(vec_cc[0]);
    }
    else if( ntrack[4] >= 20 && ntrack[4] < 40 ) {
      auto tvec = TVector2( (unitP_rot->X()-meanX[1])/sigX[1], (unitP_rot->Y()-meanY[1])/sigY[1] );
      vec_cc[1] = TVector2::Phi_mpi_pi(unitP_rot->Phi()); //tvec.Phi());  
      htf[1]->Fill(vec_cc[1]);
    }
    else if( ntrack[4] >= 40 ) {
      auto tvec = TVector2( (unitP_rot->X()-meanX[2])/sigX[2], (unitP_rot->Y()-meanY[2])/sigY[2] );
      vec_cc[2] = TVector2::Phi_mpi_pi(unitP_rot->Phi()); //tvec.Phi());  
      htf[2]->Fill(vec_cc[2]);
    }

    auto tvec = TVector2( (unitP_rot->X()-meanX[3])/sigX[3], (unitP_rot->Y()-meanY[3])/sigY[3] );
    vec_cc[3] = TVector2::Phi_mpi_pi(unitP_rot->Phi()); //tvec.Phi());  
    htf[3]->Fill(vec_cc[3]);
  }

  aleg2->AddEntry(htf[0], "ntrack[4] >= 0 && ntrack[4] < 20");  
  aleg2->AddEntry(htf[1], "ntrack[4] >= 20 && ntrack[4] < 40");
  aleg2->AddEntry(htf[2], "ntrack[4] >= 40 ");
  aleg2->AddEntry(htf[3], "All events");

  auto ccc = new TCanvas("ccc1","ccc1");

  for(UInt_t i = 0; i < 4; i++ ){
    htf[i]->SetLineColor(icolor[i]);
    htf[i]->SetNormFactor(100);

    if( i == 0 )
      htf[i]->Draw("e");
    else
      htf[i]->Draw("same e");
  }
  aleg2->Draw();

}


void vv(){
  TH1D *ht;
  auto aleg = new TLegend(0.24,0.7,0.75,0.98,"");  

  
  auto htt = new TH1D("htt",";Q_xy",100,-15., 15.);
  auto fgX  = new TF1("fgX","gaus",-30,30);

  auto ccc = new TCanvas("ccc0","ccc0");

  for(UInt_t i = 0; i < 3; i++ ){
    ht = new TH1D(Form("ht_%d",i),";#Psi", 100, -3.2, 3.2);

    TCut ntcut = Form("ntrack[4]>=%d && ntrack[4]< %d",20*i, 20*(i+1));
    
    rChain0->Project("htt", "unitP_rot.X()",ntcut);
    htt->Fit("fgX","Q0");
    
    Double_t constX= fgX->GetParameter(0);
    Double_t meanX = fgX->GetParameter(1);
    Double_t sigX  = fgX->GetParameter(2);

    ntcut.Print();
    cout << " X " << constX << " mean : " << meanX << " sig : " << sigX << endl;

    rChain0->Project("htt", "unitP_rot.Y()",ntcut);
    htt->Fit("fgX","Q0");
    
    Double_t constY= fgX->GetParameter(0);
    Double_t meanY = fgX->GetParameter(1);
    Double_t sigY  = fgX->GetParameter(2);

    cout << " Y " << constY << " mean : " << meanY << " sig : " << sigY << endl;

    rChain0->Project(Form("ht_%d",i),"unitP_rot.Phi()",ntcut);

    ht->SetLineColor(icolor[i]);
    //    ht->Set
    ht->SetNormFactor(100);
    
    aleg->AddEntry(ht, ntcut.GetTitle());

    if(i == 0)
      ht->Draw("e");
    else
      ht->Draw("same e");    
    
  }

  auto htall = new TH1D("htall",";#Psi", 100, -3.2, 3.2);
  rChain0->Project("htall","unitP_rot.Phi()");
  
  htall->SetLineColor(6);
  htall->SetNormFactor(100);
  htall->Draw("same e");
  aleg->AddEntry(htall,"ALL");
  
  aleg->Draw();

}

