#include "STAnaTaskManager.hh"

#include "FairEventHeader.h"            // for FairEventHeader
#include "FairFileHeader.h"             // for FairFileHeader
#include "FairLogger.h"                 // for FairLogger, MESSAGE_ORIGIN
#include "FairParIo.h"                  // for FairParIo
#include "FairParSet.h"                 // for FairParSet
#include "FairRootManager.h"            // for FairRootManager
#include "FairRunIdGenerator.h"         // for FairRunIdGenerator
#include "FairRuntimeDb.h"              // for FairRuntimeDb
#include "FairTask.h"                   // for FairTask
#include "FairTrajFilter.h"             // for FairTrajFilter

#include "FairFileSource.h"             // ONLY TEMPORARILY, FOR COMPABILITY
#include "FairMixedSource.h"            // ONLY TEMPORARILY, FOR COMPABILITY

#include <iosfwd>                       // for ostream
#include "TChain.h"                     // for TChain
#include "TCollection.h"                // for TIter
#include "TDirectory.h"                 // for TDirectory, gDirectory
#include "TFile.h"                      // for TFile, gFile
#include "TGeoManager.h"                // for gGeoManager, TGeoManager
#include "TKey.h"                       // for TKey
#include "TList.h"                      // for TList
#include "TNamed.h"                     // for TNamed
#include "TObjArray.h"                  // for TObjArray
#include "TObject.h"                    // for TObject
#include "TROOT.h"                      // for TROOT, gROOT
#include "TSeqCollection.h"             // for TSeqCollection
#include "TSystem.h"                    // for TSystem, gSystem
#include "TTree.h"                      // for TTree

#include <stdlib.h>                     // for NULL, exit
#include "signal.h"
#include <string.h>                     // for strcmp
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <list>                         // for list

using std::cout;
using std::endl;
using std::list;

Bool_t gFRAIsInterrupted;

//_____________________________________________________________________________
STAnaTaskManager* STAnaTaskManager::Instance()
{
  static thread_local STAnaTaskManager instance;
  return &instance;
}
//_____________________________________________________________________________
void FRA_handler_ctrlc(int)
{
  LOG(INFO) << "*********** CTRL C PRESSED *************" << FairLogger::endl;
  gFRAIsInterrupted = kTRUE;
}
//_____________________________________________________________________________

STAnaTaskManager* STAnaTaskManager::fgRinstance= 0;

//_____________________________________________________________________________
STAnaTaskManager::STAnaTaskManager()
  :FairRun(),
   fRunInfo(),
   fIsInitialized(kFALSE),
   fStatic(kFALSE),
   fTimeStamps(kFALSE),
   fInFileIsOpen(kFALSE),
   fEventTimeMin(0),
   fEventTimeMax(0),
   fEventTime(0),
   fEventMeanTime(0),
   fTimeProb(0),
   fFinishProcessingLMDFile(kFALSE),
   fFileSource(0),
   fMixedSource(0),
   fStoreEventHeader(kTRUE)

{

  fgRinstance=this;
  fAna=kTRUE;
}
//_____________________________________________________________________________

//_____________________________________________________________________________
STAnaTaskManager::~STAnaTaskManager()
{
}


//_____________________________________________________________________________

void STAnaTaskManager::Init()
{

  if (fIsInitialized) {
    LOG(FATAL) << "Error Init is already called before!" << FairLogger::endl;
    exit(-1);
  } else {
    fIsInitialized=kTRUE;
  }
  fRtdb= GetRuntimeDb();

  // Check if we have an input file to be used
  //  fInFileIsOpen = fRootManager->InitSource();

  gROOT->GetListOfBrowsables()->Add(fTask);


  /**Set the IO Manager to run with time stamps*/
  if (fTimeStamps) {
   fRootManager->RunWithTimeStamps();
  }

  std::cout << " run id ___> " << fRunId << std::endl;

  fRtdb->initContainers(fRunId);
  fFileHeader->SetRunId(fRunId);


  // Now call the User initialize for Tasks
  fTask->InitTask();
  // if the vis manager is available then initialize it!
  FairTrajFilter* fTrajFilter = FairTrajFilter::Instance();
  if (fTrajFilter) {
    fTrajFilter->Init();
  }
  // Create a list of time based branches (if any).

  fRootManager->UpdateListOfTimebasedBranches();



  // create the output tree after tasks initialisation
  fOutFile->cd();
  //TTree* outTree =new TTree("flw", "/cbmout", 99);
  TTree* outTree =new TTree(FairRootManager::GetTreeName(), "/cbmout", 99);
  LOG(INFO) << "TREENAME : " << FairRootManger::GetTreeName() << FairLogger::endl;
  fRootManager->TruncateBranchNames(outTree, "/flow");
  fRootManager->SetOutTree(outTree);
  fRootManager->CreatePersistentBranchesAny();
  fRootManager->WriteFolder();
  fRootManager->WriteFileHeader(fFileHeader);
}
//_____________________________________________________________________________

//_____________________________________________________________________________
void STAnaTaskManager::Run(Int_t Ev_start, Int_t Ev_end)
{
  gFRAIsInterrupted = kFALSE;

  if (fTimeStamps) {
    RunTSBuffers();
  } else {
    UInt_t tmpId =0;
    //  if (fInputFile==0) {
    if (!fInFileIsOpen) {
      DummyRun(Ev_start,Ev_end);
      return;
    }

   Int_t MaxAllowed=fRootManager->CheckMaxEventNo(Ev_end);
    if ( MaxAllowed != -1 ) {
      if (Ev_end==0) {
        if (Ev_start==0) {
          Ev_end=MaxAllowed;
        } else {
          Ev_end =  Ev_start;
          if ( Ev_end > MaxAllowed ) {
            Ev_end = MaxAllowed;
          }
          Ev_start=0;
        }
      } else {
        if (Ev_end > MaxAllowed) {
          cout << "-------------------Warning---------------------------" << endl;
          cout << " -W STAnaTaskManager : File has less events than requested!!" << endl;
          cout << " File contains : " << MaxAllowed  << " Events" << endl;
          cout << " Requested number of events = " <<  Ev_end <<  " Events"<< endl;
          cout << " The number of events is set to " << MaxAllowed << " Events"<< endl;
          cout << "-----------------------------------------------------" << endl;
          Ev_end = MaxAllowed;
        }
      }
      LOG(INFO) << "STAnaTaskManager::Run() After checking, the run will run from event " << Ev_start << " to " << Ev_end << "." << FairLogger::endl;
    }
    else {
      LOG(INFO) << "STAnaTaskManager::Run() continue running without stop" << FairLogger::endl;
    }

    if (fGenerateRunInfo) {
      fRunInfo.Reset();
    }

    Int_t readEventReturn = 0;

    for (int i=Ev_start; i< Ev_end || MaxAllowed==-1 ; i++) {

      gSystem->IgnoreInterrupt();
      gFRAIsInterrupted = kFALSE;
      signal(SIGINT, FRA_handler_ctrlc);

      if ( gFRAIsInterrupted ) {
        LOG(WARNING) << "STAnaTaskManager::Run() Event loop was interrupted by the user!" << FairLogger::endl;
        break;
      }

      readEventReturn = fRootManager->ReadEvent(i);

      if ( readEventReturn != 0 ) {
        LOG(WARNING) << "STAnaTaskManager::Run() fRootManager->ReadEvent(" << i << ") returned " << readEventReturn << ". Breaking the event loop" << FairLogger::endl;
        break;
      }

      fRootManager->FillEventHeader(fEvtHeader);

      tmpId = fEvtHeader->GetRunId();
      if ( tmpId != fRunId ) {
        fRunId = tmpId;
        if ( !fStatic ) {
          Reinit( fRunId );
          fTask->ReInitTask();
        }
      }
      //std::cout << "WriteoutBufferData with time: " << fRootManager->GetEventTime();
      fRootManager->StoreWriteoutBufferData(fRootManager->GetEventTime());
      fTask->ExecuteTask("");
      Fill();
      fRootManager->DeleteOldWriteoutBufferData();
      fTask->FinishEvent();

      if (fGenerateRunInfo) {
        fRunInfo.StoreInfo();
      }
      if (NULL !=  FairTrajFilter::Instance()) {
        FairTrajFilter::Instance()->Reset();
      }

    }

    fRootManager->StoreAllWriteoutBufferData();
    fTask->FinishTask();
    if (fGenerateRunInfo) {
      fRunInfo.WriteInfo();
    }
    fRootManager->LastFill();
    fRootManager->Write();
  }
}
//_____________________________________________________________________________

//_____________________________________________________________________________
void STAnaTaskManager::RunEventReco(Int_t Ev_start, Int_t Ev_end)
{
  UInt_t tmpId =0;

  Int_t MaxAllowed=fRootManager->CheckMaxEventNo(Ev_end);
  if ( MaxAllowed != -1 ) {
    if (Ev_end==0) {
      if (Ev_start==0) {
	Ev_end=MaxAllowed;
      } else {
	Ev_end =  Ev_start;
	if ( Ev_end > MaxAllowed ) {
	  Ev_end = MaxAllowed;
	}
	Ev_start=0;
      }
    } else {
      if (Ev_end > MaxAllowed) {
	cout << "-------------------Warning---------------------------" << endl;
	cout << " -W STAnaTaskManager : File has less events than requested!!" << endl;
	cout << " File contains : " << MaxAllowed  << " Events" << endl;
	cout << " Requested number of events = " <<  Ev_end <<  " Events"<< endl;
	cout << " The number of events is set to " << MaxAllowed << " Events"<< endl;
	cout << "-----------------------------------------------------" << endl;
	Ev_end = MaxAllowed;
      }
    }
    LOG(INFO) << "STAnaTaskManager::Run() After checking, the run will run from event " << Ev_start << " to " << Ev_end << "." << FairLogger::endl;
  }
  else {
    LOG(INFO) << "STAnaTaskManager::Run() continue running without stop" << FairLogger::endl;
  }

  if (fGenerateRunInfo) {
    fRunInfo.Reset();
  }

  for (int i=Ev_start; i< Ev_end; i++) {
    fRootManager->ReadEvent(i);
    /**
     * if we have simulation files then they have MC Event Header and the Run Id is in it, any way it
     * would be better to make FairMCEventHeader a subclass of FairEvtHeader.
     */
    if ( tmpId != fRunId ) {
      fRunId = tmpId;
      if ( !fStatic ) {
        Reinit( fRunId );
        fTask->ReInitTask();
      }
    }
    //FairMCEventHeader* header = dynamic_cast<FairMCEventHeader*>(fRootManager->GetObject("MCEventHeader.");
    //    std::cout << "WriteoutBufferData with time: " << fRootManager->GetEventTime();
    fRootManager->StoreWriteoutBufferData(fRootManager->GetEventTime());
    fTask->ExecuteTask("");

    fRootManager->FillEventHeader(fEvtHeader);
    // Fill();
    fTask->FinishEvent();

    if (fGenerateRunInfo) {
      fRunInfo.StoreInfo();
    }
    if (NULL !=  FairTrajFilter::Instance()) {
      FairTrajFilter::Instance()->Reset();
    }

  }

  fTask->FinishTask();
  if (fGenerateRunInfo) {
    fRunInfo.WriteInfo();
  }
  fRootManager->LastFill();
  fRootManager->Write();
}
//_____________________________________________________________________________

//_____________________________________________________________________________
void STAnaTaskManager::Run(Double_t delta_t)
{
  while (fRootManager->ReadNextEvent(delta_t)==kTRUE) {
    fTask->ExecuteTask("");
    fRootManager->FillEventHeader(fEvtHeader);
    Fill();
    fRootManager->DeleteOldWriteoutBufferData();
    fTask->FinishEvent();
    if (NULL !=  FairTrajFilter::Instance()) {
      FairTrajFilter::Instance()->Reset();
    }
  }

  fRootManager->StoreAllWriteoutBufferData();
  fTask->FinishTask();
  fRootManager->LastFill();
  fRootManager->Write();

}
//_____________________________________________________________________________

//_____________________________________________________________________________
void STAnaTaskManager::Run(Long64_t entry)
{
  UInt_t tmpId =0;
  fRootManager->ReadEvent(entry);
  tmpId = fEvtHeader->GetRunId();
  if ( tmpId != fRunId ) {
    fRunId = tmpId;
    if ( !fStatic ) {
      Reinit( fRunId );
      fTask->ReInitTask();
    }
  }
  fTask->ExecuteTask("");
  fRootManager->FillEventHeader(fEvtHeader);
  fTask->FinishTask();
  Fill();
  fRootManager->DeleteOldWriteoutBufferData();
  fRootManager->LastFill();
  fRootManager->Write();
}
//_____________________________________________________________________________

//_____________________________________________________________________________
void STAnaTaskManager::RunTSBuffers()
{
  Int_t globalEvent = 0;

  bool firstRun = true;
  while (firstRun || fRootManager->AllDataProcessed() == kFALSE) {
    firstRun = false;
    if (globalEvent < fRootManager->CheckMaxEventNo(0) ) { //this step is necessary to load in all data which is not read in via TSBuffers
      fRootManager->ReadNonTimeBasedEventFromBranches(globalEvent++);
    }
    fTask->ExecuteTask("");
    fRootManager->FillEventHeader(fEvtHeader);
    Fill();
    fRootManager->DeleteOldWriteoutBufferData();
    fTask->FinishEvent();
    if (NULL !=  FairTrajFilter::Instance()) {
      FairTrajFilter::Instance()->Reset();
    }
  }
  fRootManager->StoreAllWriteoutBufferData();
  fTask->FinishTask();
  fRootManager->LastFill();
  fRootManager->Write();
}
//_____________________________________________________________________________
//_____________________________________________________________________________

void STAnaTaskManager::RunOnLmdFiles(UInt_t NStart, UInt_t NStop)
{
  if(NStart==0 && NStop==0) {
    NStart=0;
    NStop=1000000000;
    LOG(INFO) << " Maximum number of event is set to 1E9" << FairLogger::endl;
  }
  for (UInt_t i=NStart; i< NStop; i++) {
    if ( fFinishProcessingLMDFile ) {
      i = NStop; ///Same result like break

    }

    fTask->ExecuteTask("");
    fRootManager->FillEventHeader(fEvtHeader);
    Fill();
  }

  fTask->FinishTask();
  fRootManager->Write();

}
//_____________________________________________________________________________
void STAnaTaskManager::RunOnTBData() {
      std::cout << "STAnaTaskManager::RunOnTBData " << std::endl;
        while (fRootManager->FinishRun() != kTRUE) {
		fTask->ExecuteTask("");
            Fill();
            fTask->FinishEvent();
        }

        fTask->FinishTask();
        fRootManager->LastFill();
        fRootManager->Write();
}
//_____________________________________________________________________________
void STAnaTaskManager::DummyRun(Int_t Ev_start, Int_t Ev_end)
{

  /** This methode is just for testing, if you are not sure about what you do, don't use it */
  for (int i=Ev_start; i< Ev_end; i++) {
    fTask->ExecuteTask("");
    fRootManager->FillEventHeader(fEvtHeader);
    Fill();
  }
  fTask->FinishTask();
  fRootManager->Write();

}
//_____________________________________________________________________________

//_____________________________________________________________________________
void STAnaTaskManager::TerminateRun()
{
  fRootManager->StoreAllWriteoutBufferData();
  fTask->FinishTask();
  gDirectory->SetName(fRootManager->GetOutFile()->GetName());
  //  fRunInfo.WriteInfo(); // CRASHES due to file ownership i guess...
  //   cout << ">>> SlaveTerminate fRootManager->GetInChain()->Print()" << endl;
  //   fRootManager->GetInChain()->Print();
  //   cout << ">>>------------------------------------------------<<<" << endl;
  fRootManager->LastFill();
  fRootManager->Write();
  fRootManager->CloseOutFile();
}
//_____________________________________________________________________________

void STAnaTaskManager::Reinit(UInt_t runId)
{
  // reinit procedure
  fRtdb->initContainers( runId );
}
//_____________________________________________________________________________

void  STAnaTaskManager::RunWithTimeStamps()
{
  if (fIsInitialized) {
    LOG(WARNING) << "RunWithTimeStamps has to be set before Run::Init !" << FairLogger::endl;
    exit(-1);
  } else {
    fTimeStamps=kTRUE;
    fRootManager->RunWithTimeStamps();
  }
}
//_____________________________________________________________________________

//_____________________________________________________________________________
void  STAnaTaskManager::SetContainerStatic(Bool_t tempBool)
{
  fStatic=tempBool;
  if ( fStatic ) {
    LOG(INFO) << "Parameter Cont. initialisation is static" << FairLogger::endl;
  } else {
    LOG(INFO) << "Parameter Cont. initialisation is NOT static" << FairLogger::endl;
  }
}

// BELOW FUNCTIONS SHOULD BE DELETED AND MOVED TO FairFileSource ONLY
//_____________________________________________________________________________
void STAnaTaskManager::SetInputFile(TString name)
{
  LOG(WARNING) << "STAnaTaskManager::SetInputFile is obsolete. Set it by FairFileSource" << FairLogger::endl;
  if ( fMixedSource )
    {
      LOG(ERROR) << "Mixed input already set!" << FairLogger::endl;
      return;
    }
  if ( !fFileSource )
    {
      fFileSource = new FairFileSource(name);
      SetSource(fFileSource);
      return;
    }
  fFileSource->SetInputFile(name);
}
//_____________________________________________________________________________
void STAnaTaskManager::AddFriend (TString name)
{
  LOG(WARNING) << "STAnaTaskManager::AddFriend is obsolete. Set it by FairFileSource" << FairLogger::endl;
  if ( fMixedSource )
    {
      LOG(ERROR) << "Mixed input already set!" << FairLogger::endl;
      return;
    }
  if ( !fFileSource )
    {
      LOG(ERROR) << "Input file not yet set!" << FairLogger::endl;
      return;
    }
  fFileSource->AddFriend(name);
}
//_____________________________________________________________________________
void STAnaTaskManager::AddFile(TString name)
{
  LOG(WARNING) << "STAnaTaskManager::AddFile is obsolete. Set it by FairFileSource" << FairLogger::endl;
  if ( fMixedSource )
    {
      LOG(ERROR) << "Mixed input already set!" << FairLogger::endl;
      return;
    }
  if ( !fFileSource )
    {
      LOG(ERROR) << "Input file not yet set!" << FairLogger::endl;
      return;
    }
  fFileSource->AddFile(name);
}
//_____________________________________________________________________________
// ABOVE FUNCTIONS SHOULD BE DELETED AND MOVED TO FairFileSource ONLY

// BELOW FUNCTIONS SHOULD BE DELETED AND MOVED TO FairMixedSource ONLY
//_____________________________________________________________________________
void STAnaTaskManager::SetSignalFile(TString name, UInt_t identifier )
{
  LOG(WARNING) << "STAnaTaskManager::SetSignalFile is obsolete. Set it by FairMixedSource" << FairLogger::endl;
  if (identifier==0) {
    LOG(FATAL) << " ----- Identifier 0 is reserved for background files! please use other value ------ " << FairLogger::endl;
  }
  if ( fFileSource )
    {
      LOG(ERROR) << "Standard input already set!" << FairLogger::endl;
      return;
    }
  if ( !fMixedSource )
    {
      fMixedSource = new FairMixedSource(name,identifier);
      SetSource(fMixedSource);
      return;
    }
  fMixedSource->AddSignalFile(name, identifier);
}
//_____________________________________________________________________________
void STAnaTaskManager::AddSignalFile(TString name, UInt_t identifier )
{
  LOG(WARNING) << "STAnaTaskManager::AddSignalFile is obsolete. Set it by FairMixedSource" << FairLogger::endl;
  if (identifier==0) {
    LOG(FATAL) << " ----- Identifier 0 is reserved for background files! please use other value ------ " << FairLogger::endl;
  }
  if ( fFileSource )
    {
      LOG(ERROR) << "Standard input already set!" << FairLogger::endl;
      return;
    }
  if ( !fMixedSource )
    {
      fMixedSource = new FairMixedSource(name,identifier);
      SetSource(fMixedSource);
      return;
    }
  fMixedSource->AddSignalFile(name, identifier);
}
//_____________________________________________________________________________
void STAnaTaskManager::SetBackgroundFile(TString name)
{
  LOG(WARNING) << "STAnaTaskManager::SetBackgroundFile is obsolete. Set it by FairMixedSource" << FairLogger::endl;
  if ( fFileSource )
    {
      LOG(ERROR) << "Standard input already set!" << FairLogger::endl;
      return;
    }
  if ( !fMixedSource )
    {
      fMixedSource = new FairMixedSource(name,0);
      SetSource(fMixedSource);
      return;
    }
  fMixedSource->SetBackgroundFile(name);
}
//_____________________________________________________________________________
void STAnaTaskManager::AddBackgroundFile(TString name)
{
  LOG(WARNING) << "STAnaTaskManager::AddBackgroundFile is obsolete. Set it by FairMixedSource" << FairLogger::endl;
  if ( fFileSource )
    {
      LOG(ERROR) << "Standard input already set!" << FairLogger::endl;
      return;
    }
  if ( !fMixedSource )
    {
      LOG(ERROR) << "Background file not yet set!" << FairLogger::endl;
      return;
    }
  fMixedSource->AddBackgroundFile(name);
}
//_____________________________________________________________________________
void  STAnaTaskManager::BGWindowWidthNo(UInt_t background, UInt_t Signalid)
{
  LOG(WARNING) << "STAnaTaskManager::BGWindowWidthNo is obsolete. Set it by FairMixedSource" << FairLogger::endl;
  if ( fFileSource )
    {
      LOG(ERROR) << "Standard input already set!" << FairLogger::endl;
      return;
    }
  if ( !fMixedSource )
    {
      LOG(ERROR) << "Background file not yet set!" << FairLogger::endl;
      return;
    }
  fMixedSource->BGWindowWidthNo(background, Signalid);
}
//_____________________________________________________________________________
void  STAnaTaskManager::BGWindowWidthTime(Double_t background, UInt_t Signalid)
{
  LOG(WARNING) << "STAnaTaskManager::BGWindowWidthTime is obsolete. Set it by FairMixedSource" << FairLogger::endl;
  if ( fFileSource )
    {
      LOG(ERROR) << "Standard input already set!" << FairLogger::endl;
      return;
    }
  if ( !fMixedSource )
    {
      LOG(ERROR) << "Background file not yet set!" << FairLogger::endl;
      return;
    }
  fMixedSource->BGWindowWidthTime(background, Signalid);
}
//_____________________________________________________________________________
// ABOVE FUNCTIONS SHOULD BE DELETED AND MOVED TO FairMixedSource ONLY

// BELOW FUNCTIONS SHOULD BE DELETED AND MOVED TO FairFileSource AND FairMixedSource ONLY
//_____________________________________________________________________________
void STAnaTaskManager::SetEventTimeInterval(Double_t min, Double_t max)
{
  LOG(WARNING) << "STAnaTaskManager::SetEventTimeInterval is obsolete. Set it by FairSource" << FairLogger::endl;
  if ( fFileSource )
    {
      fFileSource->SetEventTimeInterval(min,max);
      return;
    }
  if ( fMixedSource )
    {
      fMixedSource->SetEventTimeInterval(min,max);
      return;
    }
  LOG(ERROR) << "SetEventTimeInterval only by input source!" << FairLogger::endl;
}
//_____________________________________________________________________________
void  STAnaTaskManager::SetEventMeanTime(Double_t mean)
{
  LOG(WARNING) << "STAnaTaskManager::SetEventMeanTime is obsolete. Set it by FairSource" << FairLogger::endl;
  if ( fFileSource )
    {
      fFileSource->SetEventMeanTime(mean);
      return;
    }
  if ( fMixedSource )
    {
      fMixedSource->SetEventMeanTime(mean);
      return;
    }
  LOG(ERROR) << "SetEventMeanTime only by input source!" << FairLogger::endl;
}
//_____________________________________________________________________________
void STAnaTaskManager::SetBeamTime(Double_t beamTime, Double_t gapTime)
{
  LOG(WARNING) << "STAnaTaskManager::SetBeamTime is obsolete. Set it by FairSource" << FairLogger::endl;
  if ( fFileSource )
    {
      fFileSource->SetBeamTime(beamTime, gapTime);
      return;
    }
  if ( fMixedSource )
    {
      fMixedSource->SetBeamTime(beamTime, gapTime);
      return;
    }
  LOG(ERROR) << "SetBeamTime only by input source!" << FairLogger::endl;
}
//_____________________________________________________________________________

//_____________________________________________________________________________
Int_t  FairRootManager::AddBranchToList(const char* name)
{
  if(fBranchNameList->FindObject(name)==0) {
    fBranchNameList->AddLast(new TObjString(name));
    fBranchSeqId++;
  }
  return fBranchSeqId;
}



//_____________________________________________________________________________
void STAnaTaskManager::Fill()
{
  if(fMarkFill)
  {
    fRootManager->Fill();
  }
  else
  {
    fMarkFill = kTRUE;
  }
}
//_____________________________________________________________________________


// void  STAnaTaskManager::SetMixAllInputs(Bool_t Status)
// {
//    fLogger->Info(MESSAGE_ORIGIN, "Mixing for all input is choosed, in this mode one event per input file is read per ste

p");
//    fRootManager->SetMixAllInputs(Status);
// }
//_____________________________________________________________________________
// ABOVE FUNCTIONS SHOULD BE DELETED AND MOVED TO FairFileSource AND FairMixedSource ONLY


ClassImp(STAnaTaskManager)
