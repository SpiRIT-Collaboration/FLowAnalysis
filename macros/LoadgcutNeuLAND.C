{
  TFile *_file0 = TFile::Open("data/gcutNLNeutron.root");
  TCutG *gcutNLNeutron =(TCutG*)_file0->Get("gcutNLNeutron");
  _file0->Close();
  gcutNLNeutron->Print();

  TFile *_file1 = TFile::Open("data/gcutNLProton.root");
  TCutG *gcutNLProton =(TCutG*)_file1->Get("gcutNLProton");
  _file1->Close();
  gcutNLProton->Print();

  TFile *_file2 = TFile::Open("data/gcutNLDeuteron.root");
  TCutG *gcutNLDeuteron =(TCutG*)_file2->Get("gcutNLDeuteron");
  _file2->Close();
  gcutNLDeuteron->Print();

  TFile *_file3 = TFile::Open("data/gcutNLTriton.root");
  TCutG *gcutNLTriton =(TCutG*)_file3->Get("gcutNLTriton");
  _file3->Close();
  gcutNLTriton->Print();

}
