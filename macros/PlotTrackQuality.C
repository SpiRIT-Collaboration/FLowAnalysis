{
  if( kFALSE ){

    TCanvas cvx1("cvx1","cvx1",1200,1500);
    cvx1.Divide(3,4);
    UInt_t idd = 1;

    cvx1.cd(idd); idd++;
    rChain->Draw("fdEdx:fP>>h(200,0.,2500,200,0,800)","fReactionPlanef%2==1","colz");

    cvx1.cd(idd); 
    rChain->Draw("fBBMass:fP>>h2(200,0.,3500,200,0,8000)","fReactionPlanef%2==1","colz");
    cvx1.GetPad(idd)->SetLogz(); idd++;

    cvx1.cd(idd); idd++;
    rChain->Draw("fRotatedP3.Phi():fRotatedP3.Theta()","fReactionPlanef%2==1","colz");

    cvx1.cd(idd); 
    rChain->Draw("fRotatedP3.Phi():fRotatedP3.Theta()","fReactionPlanef%2==0","colz");
    cvx1.GetPad(idd)->SetLogz(); idd++;

    cvx1.cd(idd); idd++;
    rChain->Draw("mtrack4>>h3(65,0,65)");

    cvx1.cd(idd); idd++;
    rChain->Draw("TVector2::Phi_mpi_pi(unitP_2fc.Phi()-unitP_1fc.Phi()):unitP_fc.Phi()","","colz");

    TH1D *h4 = new TH1D("h4","h4",100,-3.14,3.14);
    rChain->Project("h4","unitP_fc.Phi()");
    h4->SetNormFactor(100);
    cvx1.cd(idd); idd++;
    h4->Draw("e");
  }

  TCanvas cvx2("cvx2","cvx2",1200,1500);
  cvx2.Divide(3,4);
  idd = 1;

  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid11(200,-1.,1.3,200,0.,800)","fPID_norm==2212","colz");
  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid12(200,-1.,1.3,200,0.,800)","fPID_tight==2212","colz");
  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid13(200,-1.,1.3,200,0.,800)","fPID_loose==2212","colz");


  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid21(200,-0.5,1.,200,0.,800)","fPID_norm==1000010020","colz");
  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid22(200,-0.5,1.,200,0.,800)","fPID_tight==1000010020","colz");
  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid23(200,-0.5,1.,200,0.,800)","fPID_loose==1000010020","colz");

  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid31(200,-0.5,1.,200,0.,800)","fPID_norm==1000010030","colz");
  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid32(200,-0.5,1.,200,0.,800)","fPID_tight==1000010030","colz");
  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid33(200,-0.5,1.,200,0.,800)","fPID_loose==1000010030","colz");

  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid41(200,-0.5,1.,200,0.,800)","fPID_norm==1000020030","colz");
  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid42(200,-0.5,1.,200,0.,800)","fPID_tight==1000020030","colz");
  cvx2.cd(idd); idd++;
  rChain->Draw("fRotatedP3.Pt():fRapiditycm>>hpid43(200,-0.5,1.,200,0.,800)","fPID_loose==1000020030","colz");




}

