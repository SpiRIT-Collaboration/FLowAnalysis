#include <TFile.h>
#include <TParameter.h>
#include <TH1D.h>
#include <TMath.h>

#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cRED "\033[1;31m"

std::map<Int_t, Double_t> xsAbs;
std::map<Int_t, Double_t> xsAbsSigma;
std::map<Int_t, Double_t> bMax;
std::map<Int_t, Double_t> bMaxSigma;
std::map<Int_t, TH1D *> hBhat;
std::map<Int_t, TH1D *> hMult;
std::map<Int_t, TH1D *> hMultNorm;
std::map<Int_t, TH1D *> hMultWeighted;

Double_t CalculateNReturn(Int_t returnIndex, Int_t system, Int_t minM, Int_t maxM = 100) {
  if (!(system == 132 || system == 124 || system == 112 || system == 108)) {
    cout << " = You specified wrong system!" << endl;
    cout << " = Available systems = {108, 112, 124, 132}" << endl;

    return -9999;
  }
  auto minBin = hMultWeighted[system] -> GetBin(minM) + 1;
  auto maxBin = hMultWeighted[system] -> GetBin(maxM) + 1;
  auto numeratorError = 0.;
  auto denominatorError = 0.;
  auto weightedSum = hMultWeighted[system] -> IntegralAndError(minBin, maxBin, numeratorError);
  auto probInRange = hMultNorm[system] -> IntegralAndError(minBin, maxBin, denominatorError);

  auto bHatAvg = weightedSum/probInRange;
  auto bHatAvgError = TMath::Sqrt(TMath::Power(numeratorError/probInRange, 2) + TMath::Power(denominatorError*weightedSum, 2)/TMath::Power(probInRange, 4));
  auto bAvg = bHatAvg*bMax[system];
  auto bAvgError = TMath::Sqrt(TMath::Power(bHatAvgError*bMax[system], 2) + TMath::Power(bHatAvg*bMaxSigma[system], 2));

  auto bHatMin = hBhat[system] -> GetBinContent(maxM + 1);
  auto bHatMinError = hBhat[system] -> GetBinError(maxM + 1);
  auto bHatMax = hBhat[system] -> GetBinContent(minM + 1);
  auto bHatMaxError = hBhat[system] -> GetBinError(minM + 1);

  switch (returnIndex) {
    case 0:
      cout << "============================================================================================" << endl;
      cout << " Cuts for multiplicity distribution: Beam 2-sigma + On-target + dist<20 cuts" << endl;
      cout << endl;
      cout << cRED << " Ignore error values." << cNORMAL << endl;
      cout << cRED << " Do not use average b yet." << cNORMAL << endl;
      cout << endl;
      cout << " bMax of " << system << "Sn system = " << bMax[system] << " fm" << endl;
      cout << " Multiplicity range = [" << minM << ", " << maxM << "] (Max. value meaningless if you didn't put it in)" << endl;
      cout << " Average bHat = " << bHatAvg << " +- " << bHatAvgError << " (bHat=[" << bHatMin << ", " << bHatMax << "])"<< endl;
      cout << " Average b [fm] = " << bAvg << " +- " << bAvgError << " (% Error=" << bAvgError/bAvg*100. << "%, b=[" << bHatMin*bMax[system] << ", " << bHatMax*bMax[system] << "])" << endl;
      cout << " #events in multiplicity cut = " << (Long64_t) hMult[system] -> Integral(minBin, maxBin) << " / " << (Long64_t) hMult[system] -> Integral() << " (" << probInRange*100. << "%)" << endl;
      cout << "============================================================================================" << endl;

      return 0;

    case 1: return bHatMax*bMax[system];
    case 2: return bHatMin*bMax[system];
    case 3: return bHatMaxError*bMax[system];
    case 4: return bHatMinError*bMax[system];
    case 5: return hMult[system] -> Integral();
    case 6: return hMult[system] -> Integral(minBin, maxBin);
    case 7: return bAvg;
    case 8: return bAvgError;
  }

  return 0;
}

void Print(Int_t system, Int_t minM, Int_t maxM = 100) { CalculateNReturn(0, system, minM, maxM); }
Double_t GetMaxB(Int_t system, Int_t minM, Int_t maxM = 100) { return CalculateNReturn(1, system, minM, maxM); }
Double_t GetMinB(Int_t system, Int_t minM, Int_t maxM = 100) { return CalculateNReturn(2, system, minM, maxM); }
Double_t GetMaxBError(Int_t system, Int_t minM, Int_t maxM = 100) { return CalculateNReturn(3, system, minM, maxM); }
Double_t GetMinBError(Int_t system, Int_t minM, Int_t maxM = 100) { return CalculateNReturn(4, system, minM, maxM); }
Long64_t GetNumTotalEvents(Int_t system) { return CalculateNReturn(5, system, 0, 100); }
Long64_t GetNumEvents(Int_t system, Int_t minM, Int_t maxM = 100) { return CalculateNReturn(6, system, minM, maxM); }
Double_t GetAvgB(Int_t system, Int_t minM, Int_t maxM = 100) { return CalculateNReturn(7, system, minM, maxM); }
Double_t GetAvgBError(Int_t system, Int_t minM, Int_t maxM = 100) { return CalculateNReturn(8, system, minM, maxM); }

void help() {
  cout << "" << endl;
  cout << "  = You can omit the value in [ ]." << endl;
  cout << "  = Systems = {108, 112, 124, 132}" << endl;
  cout << "" << endl;
  cout << "  = Print(system, minMult[, maxMult])" << endl;
  cout << "  = GetMinB(system, minMult[, maxMult])" << endl;
  cout << "  = GetMaxB(system, minMult[, maxMult])" << endl;
  cout << "  = GetMinBError(system, minMult[, maxMult])" << endl;
  cout << "  = GetMaxBError(system, minMult[, maxMult})" << endl;
  cout << "  = GetNumTotalEvents(system)" << endl;
  cout << "  = GetNumEvents(system, minMult[, maxMult])" << endl;
  cout << "  = GetAvgB(system, minMult[, maxMult])" << endl;
  cout << "  = GetAvgBError(system, minMult[, maxMult])" << endl;
  cout << "" << endl;
  cout << "  = help() to see this menu." << endl;
  cout << "" << endl;
  cout << cYELLOW << "  = All error values are not reliable in this version!!!!" << cNORMAL << endl;
  cout << "" << endl;
}

void calculateImpactParameter() {
  auto db = new TFile("data/bDB.root");
  if (!db -> IsOpen()) {
    cout << endl;
    cout << " = No bDB.root file is found!" << endl;
    cout << endl;

    exit(0);
  }

//  for (auto system : {108, 112, 124, 132}) {
  for (auto system : {108, 132, 112, 124}) {
    // [barn]
    xsAbs[system] = ((TParameter<Double_t> *) db -> Get(Form("xsAbs%d", system))) -> GetVal();
    xsAbsSigma[system] = ((TParameter<Double_t> *) db -> Get(Form("xsAbsSigma%d", system))) -> GetVal();

    // [fm] barn = 100fm^2
    bMax[system] = ((TParameter<Double_t> *) db -> Get(Form("bMax%d", system))) -> GetVal();

    // First order approx. The values are the same up to 1E-2
    bMaxSigma[system] = ((TParameter<Double_t> *) db -> Get(Form("bMaxSigma%d", system))) -> GetVal();

    hMult[system] = (TH1D *) db -> Get(Form("hMult%d", system));
    hMultNorm[system] = (TH1D *) db -> Get(Form("hMultNorm%d", system));
    hBhat[system] = (TH1D *) db -> Get(Form("hBhat%d", system)); // Each bin of the histogram is the sum up to that bin from infinity
    hMultWeighted[system] = (TH1D *) db -> Get(Form("hMultWeighted%d", system));
  }

  if (!gROOT -> IsBatch())
    help();
}
