#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <g4centrality/PHG4CentralityReco.h>

#include <HIJetReco.C>
#include <JetValidation.h>
#include<Calo_Calib.C>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libJetValidation.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libg4dst.so)


#endif
void Fun4All_JetVal(const char *filelistcalo = "dst_calo_cluster.list",
		    const int jobID = 0)
{

  //Process_Calo_Calib();

  char outname[80];
  sprintf(outname, "/gpfs/mnt/gpfs02/sphenix/user/mstojano/jetPt0MB/joboutput_%d.root", jobID);
  //sprintf(outname, "outputjob.root");
 
  Fun4AllServer *se = Fun4AllServer::instance();
  int verbosity = 0;


  ifstream file(filelistcalo);
  string first_file;
  getline(file, first_file);

  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(first_file);
  int runnumber = runseg.first;
  cout << "run number = " << runnumber << endl;

  se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();
  //rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  //rc->set_uint64Flag("TIMESTAMP", runnumber);

  PHG4CentralityReco *cent = new PHG4CentralityReco();
  cent->Verbosity(verbosity);
  cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
  se->registerSubsystem( cent );
  
  Enable::VERBOSITY = verbosity;
  HIJetReco();

    //Process_Calo_Calib();

  JetValidation *myJetVal = new JetValidation("AntiKt_Tower_r04_Sub1", "AntiKt_Truth_r04", outname);

  myJetVal->setPtRange(5, 100);
  myJetVal->setEtaRange(-1.1, 1.1);
  myJetVal->doUnsub(1);
  myJetVal->doTruth(0);
  myJetVal->doSeeds(1);
  se->registerSubsystem(myJetVal);
 
  //Fun4AllInputManager *intrue = new Fun4AllDstInputManager("DSTtruth");
  //intrue->AddListFile(filelisttruth,1);
  //se->registerInputManager(intrue);

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTcalo");
  in2->AddListFile(filelistcalo,1);
  se->registerInputManager(in2);

  //Fun4AllInputManager *in3 = new Fun4AllDstInputManager("DSTglobal");
  //in3->AddListFile(filelistglobal,1);
  //se->registerInputManager(in3);
  
  se->run(-1);
  se->End();

  gSystem->Exit(0);
  return 0;

}
