

void ProcessInfo(Long64_t nevt, Long64_t ievt)
{
  TDatime dtime;
  TDatime beginTime;

  UInt_t norm =  nevt > 50 ? 50: 1;

  beginTime.Copy(dtime);
  if(ievt%(UInt_t)(nevt/norm) == 0) {
    dtime.Set();
    Int_t ptime = dtime.Get() - beginTime.Get();

    std::cout << "Processing .... "
	      << setw(4) << Int_t(((Double_t)ievt/(Double_t)nevt)*100.) << " % = "
	      << setw(8) << ievt << "/"<< nevt
	      << "--->"
	      << dtime.AsString() << " ---- "
	      << std::endl;
  }
}
