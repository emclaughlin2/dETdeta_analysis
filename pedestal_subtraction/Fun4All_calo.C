#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
//#include <jetbase/FastJetAlgo.h>
//#include <jetbase/JetReco.h>
//#include <jetbase/TowerJetInput.h>
//#include <g4jets/TruthJetInput.h>
#include <fstream>
#include <phool/recoConsts.h>
#include <TSystem.h>
#include "mdctreemaker/MDCTreeMaker.h"
#include <caloreco/CaloTowerCalib.h>
#include <g4mbd/MbdDigitization.h>
#include <mbd/MbdReco.h>
#include <frog/FROG.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <centrality/CentralityReco.h>
#include <calotrigger/MinimumBiasClassifier.h>
//#include <G4Setup_sPHENIX.C>
//#include <energycorrection/EnergyCorrection.h>

using namespace std;

R__LOAD_LIBRARY(libg4centrality.so)
//R__LOAD_LIBRARY(libEnergyCorrection.so)
R__LOAD_LIBRARY(libFROG.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libg4vertex.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libmdctreemaker.so)
R__LOAD_LIBRARY(libcalo_io.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4mbd.so)
R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcalotrigger.so)

void Fun4All_calo()
{
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();

  MDCTreeMaker *tt = new MDCTreeMaker("hcal_pedestal_data_output.root", 0, 0, 0);
  se->registerSubsystem( tt );

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTcalo");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0000.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0001.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0002.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0003.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0004.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0005.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0006.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0007.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0008.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0009.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0010.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0011.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0012.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0013.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0014.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0015.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0016.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0017.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0018.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0019.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0020.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0021.root");
  in->AddFile("/sphenix/user/egm2153/calib_study/zs_testing/pedestal_subtraction_hcal/Run30146_Processing/DST_CALOR-00030146-0022.root");
  se->registerInputManager(in);

  se->run(-1);
  se->End();
  se->PrintTimer();
  gSystem->Exit(0);
  return 0;
}
