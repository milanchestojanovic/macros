//module for producing a TTree with jet information for doing jet validation studies
// for questions/bugs please contact Virginia Bailey vbailey13@gsu.edu
#include <fun4all/Fun4AllBase.h>
#include <JetValidation.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <jetbase/JetMap.h>
#include <jetbase/JetContainer.h>
#include <jetbase/Jetv2.h>
#include <jetbase/Jetv1.h>
#include <centrality/CentralityInfo.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <ffarawobjects/Gl1Packet.h>
#include <jetbackground/TowerBackground.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerDefs.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeomContainer.h>

#include <TTree.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include <TTree.h>
#include <numeric>  // For std::accumulate

//____________________________________________________________________________..
JetValidation::JetValidation(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename):
  SubsysReco("JetValidation_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-1, 1)
  , m_ptRange(5, 100)
  , m_doTruthJets(0)
  , m_doSeeds(0)
  , m_doUnsubJet(0)
  , m_T(nullptr)
  , m_event(-1)
  , m_nTruthJet(-1)
  , m_nJet(-1)
  , m_id()
  , m_nComponent()
  , m_eta()
  , m_phi()
  , m_e()
  , m_pt()
  , m_cleta()
  , m_clphi()
  , m_cle()
  , m_clecore()
  , m_clpt()
  , m_clprob()
  , m_sub_et()
  , m_truthID()
  , m_truthNComponent()
  , m_truthEta()
  , m_truthPhi()
  , m_truthE()
  , m_truthPt()
  , m_eta_rawseed()
  , m_phi_rawseed()
  , m_pt_rawseed()
  , m_e_rawseed()
  , m_rawseed_cut()
  , m_eta_subseed()
  , m_phi_subseed()
  , m_pt_subseed()
  , m_e_subseed()
  , m_subseed_cut()
{
  std::cout << "JetValidation::JetValidation(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
JetValidation::~JetValidation()
{
  std::cout << "JetValidation::~JetValidation() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
//____________________________________________________________________________..
int JetValidation::Init(PHCompositeNode *topNode)
{
  std::cout << "JetValidation::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  std::cout << "JetValidation::Init - Output to " << m_outputFileName << std::endl;

  m_T = new TTree("T", "MyJetAnalysis Tree");
    m_T->Branch("jetPt",&jetPt);
    m_T->Branch("jetEta",&jetEta);
    m_T->Branch("jetPhi",&jetPhi);
    m_T->Branch("triggerVector", &m_triggerVector);
    m_T->Branch("cetmMult",&m_cetmMult);
    m_T->Branch("cemMult",&m_cemMult);
    m_T->Branch("z_vertex",&m_vertex);


/*
  // configure Tree
  m_T = new TTree("T", "MyJetAnalysis Tree");
  m_T->Branch("m_event", &m_event, "event/I");
  m_T->Branch("nJet", &m_nJet, "nJet/I");
  m_T->Branch("cent", &m_centrality);
  m_T->Branch("zvtx", &m_zvtx);
  m_T->Branch("b", &m_impactparam);
  m_T->Branch("id", &m_id);
  m_T->Branch("nComponent", &m_nComponent);
  m_T->Branch("triggerVector", &m_triggerVector);

  m_T->Branch("eta", &m_eta);
  m_T->Branch("phi", &m_phi);
  m_T->Branch("e", &m_e);
  m_T->Branch("pt", &m_pt);

  m_T->Branch("cleta", &m_cleta);
  m_T->Branch("clphi", &m_clphi);
  m_T->Branch("cle", &m_cle);
  m_T->Branch("clecore", &m_clecore);
  m_T->Branch("clpt", &m_clpt);
  m_T->Branch("clprob", &m_clprob);

  if(m_doUnsubJet)
    {
      m_T->Branch("pt_unsub", &m_unsub_pt);
      m_T->Branch("subtracted_et", &m_sub_et);
    }
  if(m_doTruthJets){
    m_T->Branch("nTruthJet", &m_nTruthJet);
    m_T->Branch("truthID", &m_truthID);
    m_T->Branch("truthNComponent", &m_truthNComponent);
    m_T->Branch("truthEta", &m_truthEta);
    m_T->Branch("truthPhi", &m_truthPhi);
    m_T->Branch("truthE", &m_truthE);
    m_T->Branch("truthPt", &m_truthPt);
  }

  if(m_doSeeds){
    m_T->Branch("rawseedEta", &m_eta_rawseed);
    m_T->Branch("rawseedPhi", &m_phi_rawseed);
    m_T->Branch("rawseedPt", &m_pt_rawseed);
    m_T->Branch("rawseedE", &m_e_rawseed);
    m_T->Branch("rawseedCut", &m_rawseed_cut);
    m_T->Branch("subseedEta", &m_eta_subseed);
    m_T->Branch("subseedPhi", &m_phi_subseed);
    m_T->Branch("subseedPt", &m_pt_subseed);
    m_T->Branch("subseedE", &m_e_subseed);
    m_T->Branch("subseedCut", &m_subseed_cut);
  }
 
  m_T->Branch("towerEta", &m_eta_tower);
  m_T->Branch("retowerEta", &m_eta_retower);

  m_T->Branch("toweriEta", &m_eta_itower);
  m_T->Branch("retoweriEta", &m_eta_iretower);

  m_T->Branch("towerPhi", &m_phi_tower);
  m_T->Branch("retowerPhi", &m_phi_retower);

  m_T->Branch("toweriPhi", &m_phi_itower);
  m_T->Branch("retoweriPhi", &m_phi_iretower);
*/
  //2pc
  doSymmetrisation_= true; doBackground_ = true; doEffCorr_ = false;

  ptTbins_= 5; ptNbins_= 5;
  triggerH.resize(ptTbins_);
  triggerH2D.resize(ptTbins_);
  fullSet.resize(12);
  event.resize(10);
  signal.resize(ptTbins_);
  background.resize(ptTbins_);
  for(int i=0; i<ptTbins_; i++) {
    signal[i].resize(ptNbins_);
    background[i].resize(ptNbins_);
  }      

  //histo_trigger_all = fs->make<TH2D>(Form("Trigger_all") , "trigger_all" , 10000 , 0 , 10, occBins_.size()-1, &occBins_[0] );

  for (Int_t j=0; j<ptTbins_; j++ ) {
    triggerH[j]=new TH1F(Form("trig_%d",j) , ";#eta;" , 96,-1.1,1.1);
    triggerH2D[j]=new TH2F(Form("trig2D_%d",j) , ";#eta;#phi" , 96,-1.1,1.1,62,-TMath::Pi()+TMath::Pi(), TMath::Pi());

    for (Int_t i=0; i<ptNbins_; i++) {
   signal[j][i]=new TH2F(Form("signal_trig_%d_%d",j,i) , ";#Delta#eta;#Delta#phi;" , 21,-2.2,2.2,  31, -0.5*TMath::Pi()+TMath::Pi()/32, 1.5*TMath::Pi()-TMath::Pi()/32);
   background[j][i]=new TH2F(Form("background_trig_%d_%d",j,i) , ";#Delta#eta;#Delta#phi;" , 21,-2.2,2.2,  31, -0.5*TMath::Pi()+TMath::Pi()/32, 1.5*TMath::Pi()-TMath::Pi()/32);
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetValidation::InitRun(PHCompositeNode *topNode)
{
  std::cout << "JetValidation::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetValidation::process_event(PHCompositeNode *topNode)
{
  //  std::cout << "JetValidation::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  ++m_event;

  if (m_event%1000==0) std::cout<<m_event<<std::endl;

  m_vertex = -99999;

  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

  if (!vertexmap)
  {
    std::cout << "GlobalVertexMap node is missing" << std::endl;
  }
  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx)
    {
      m_vertex = vtx->get_z();
    }
  }
  //std::cout<<"z vertex = "<<m_vertex<<std::endl;
  if (fabs(m_vertex)>50) return Fun4AllReturnCodes::ABORTEVENT;
/*
  TowerInfoContainer *ohcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  TowerInfoContainer *ihcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");

  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  RawTowerGeomContainer *tower_geomIN = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");

  unsigned int tower_range = 0;
*/
  for (int i=0; i<10; i++) event[i].resize(0);
  for (int i=0; i<10; i++) event[i].shrink_to_fit();

  Gl1Packet *gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (!gl1PacketInfo)
  {
    std::cout << PHWHERE << "caloTreeGen::process_event: " << "GL1Packet" << " node is missing. Output related to this node will be empty" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (gl1PacketInfo)
  {
    uint64_t triggervec = gl1PacketInfo->getScaledVector();
    for (int i = 0; i < 64; i++)
    {
      bool trig_decision = ((triggervec & 0x1U) == 0x1U);
      m_triggerVector.push_back(trig_decision);
      triggervec = (triggervec >> 1U) & 0xffffffffU;
    }
  }

  //if (m_triggerVector[21] || m_triggerVector[22] || m_triggerVector[23]) return Fun4AllReturnCodes::ABORTEVENT;

  // interface to reco jets
  //int Njets = 0;
  jetPt=0.;
  jetPhi=-999.;
  jetEta=-999.;

  //JetContainer* jets = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  JetContainer* jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Tower_r04_Sub1");
  //std::cout<<"m_recoJetName = "<<m_recoJetName<<std::endl;
  for (auto jet : *jets)
    {
      if(jet->get_pt() < 1) continue; // to remove noise jets
      //if(jet->get_pt()>20) Njets++;
        //std::cout<<"jet pt = "<<jet->get_pt()<<std::endl;
        //particle.SetPtEtaPhi(jet->get_pt(),jet->get_eta(),jet->get_phi());
        //if (jet->get_pt()>20) event[3].push_back(particle);
      if(jet->get_pt()>jetPt) {
        jetPt=jet->get_pt();
        jetEta=jet->get_eta();
        jetPhi=jet->get_phi();
      }
      	
    }
   //if (Njets==0)     return Fun4AllReturnCodes::ABORTEVENT;

  //fill the tree
  m_T->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
  //double mean_eta = 0.0;
  //return Fun4AllReturnCodes::EVENT_OK;

  //totalClusterEEMCal = 0;

  m_cetmMult=0;
  m_cemMult=0;

  RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
   RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
    RawClusterContainer::ConstIterator clusterIter;
    for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
    {
      RawCluster *recoCluster = clusterIter->second;

      CLHEP::Hep3Vector vertex(0, 0, 0);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);
      CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*recoCluster, vertex);

      //float clusE = E_vec_cluster_Full.mag();
      float clus_eta = E_vec_cluster.pseudoRapidity();
      float clus_phi = E_vec_cluster.phi();
      //float clus_pt = E_vec_cluster.perp();
      //if (clus_pt>0.5) totalClusterEEMCal += clusE;

      //particle.SetPtEtaPhi(clus_pt,clus_eta,clus_phi);
      particle.SetPtEtaPhi(E_vec_cluster.mag()/cosh(clus_eta),clus_eta,clus_phi);
      //if (clus_pt>0.5) event[0].push_back(particle);
      //if (clus_pt>0.5) event[5].push_back(particle);
      if (recoCluster->get_chi2() < 4 && E_vec_cluster.mag()/cosh(clus_eta)>.3) event[0].push_back(particle);
      if (recoCluster->get_chi2() < 4 && E_vec_cluster.mag()/cosh(clus_eta)>.3) event[5].push_back(particle);
      if (recoCluster->get_chi2() < 4 && E_vec_cluster.mag()/cosh(clus_eta)>.3) m_cetmMult++;
      if (recoCluster->get_chi2() < 4 && E_vec_cluster.mag()>.5) m_cemMult++;

      //if (recoCluster->get_chi2() < 4 && E_vec_cluster.mag()/clus_eta>.5) std::cout<<"E = "<<E_vec_cluster.mag()<< " eta = "<<clus_eta<<"  E_T = "<<E_vec_cluster.mag()/clus_eta<<std::endl;

   } 
  if (m_cetmMult<2)  return Fun4AllReturnCodes::ABORTEVENT;
     
  //m_T->Fill();
  //return Fun4AllReturnCodes::EVENT_OK;

  TowerInfoContainer *emcTowerContainer;
  emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");

  TowerInfoContainer *emcReTowerContainer;
  emcReTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");

  TowerInfoContainer *ohcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  TowerInfoContainer *ihcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");

  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  RawTowerGeomContainer *tower_geomIN = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");

  unsigned int tower_range = 0;
  if (emcTowerContainer)
  {
    tower_range = emcTowerContainer->size();
    for (unsigned int iter = 0; iter < tower_range; iter++)
    {
      unsigned int towerkey = emcTowerContainer->encode_key(iter);
      unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
      unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
      double energy = emcTowerContainer->get_tower_at_channel(iter)->get_energy();
      float chi2 = emcTowerContainer->get_tower_at_channel(iter)->get_chi2();

      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
                  float tower_phi = tower_geomEM->get_tower_geometry(key)->get_phi();
                  float tower_eta = tower_geomEM->get_tower_geometry(key)->get_eta();

      //const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
      //            float tower_phi = tower_geomIN->get_tower_geometry(key)->get_phi();
      //            float tower_eta = tower_geomIN->get_tower_geometry(key)->get_eta();

       particle.SetPtEtaPhi(energy,tower_eta,tower_phi);
       if (chi2<1e4 && energy/cosh(tower_eta)>0.35) 
       event[1].push_back(particle);
       if (chi2<1e4 && energy/cosh(tower_eta)>0.35) event[6].push_back(particle);

       //m_eta_tower.push_back(tower_eta);
       //m_phi_tower.push_back(tower_phi);

       //m_eta_itower.push_back(ieta);
       //m_phi_itower.push_back(iphi);

    }
  }

  if (emcReTowerContainer)
  {
    tower_range = emcReTowerContainer->size();
    for (unsigned int iter = 0; iter < tower_range; iter++)
    {
      unsigned int towerkey = emcReTowerContainer->encode_key(iter);
      unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
      unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
      double energy = emcReTowerContainer->get_tower_at_channel(iter)->get_energy();
      float chi2 = emcReTowerContainer->get_tower_at_channel(iter)->get_chi2();

      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
                  float tower_phi = tower_geomIN->get_tower_geometry(key)->get_phi();
                  float tower_eta = tower_geomIN->get_tower_geometry(key)->get_eta();

       particle.SetPtEtaPhi(energy,tower_eta,tower_phi);
       if (chi2<1e4 && energy/cosh(tower_eta)>0.35) event[2].push_back(particle);
       if (chi2<1e4 && energy/cosh(tower_eta)>0.35) event[7].push_back(particle);

       //m_eta_retower.push_back(tower_eta);
       //m_phi_retower.push_back(tower_phi);

       //m_eta_iretower.push_back(ieta);
       //m_phi_iretower.push_back(iphi);

    }
  }

//  TowerInfoContainer *ohcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
//  TowerInfoContainer *ihcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
//
//  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
//  RawTowerGeomContainer *tower_geomIN = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");

  if (ohcTowerContainer && ihcTowerContainer)
  {
    tower_range = ohcTowerContainer->size();
    for (unsigned int iter = 0; iter < tower_range; iter++)
    {
      unsigned int towerkey = ohcTowerContainer->encode_key(iter);
      unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
      unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
      double energy = ohcTowerContainer->get_tower_at_channel(iter)->get_energy();
      float chi2 = ohcTowerContainer->get_tower_at_channel(iter)->get_chi2();

      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
                  float tower_phi = tower_geomOH->get_tower_geometry(key)->get_phi();
                  float tower_eta = tower_geomOH->get_tower_geometry(key)->get_eta();
        particle.SetPtEtaPhi(energy,tower_eta,tower_phi);
        if (chi2>0 && energy/cosh(tower_eta)>0.1) 
		event[3].push_back(particle);
        if (chi2>0 && energy/cosh(tower_eta)>0.1) event[8].push_back(particle);
        if (chi2>0 && energy/cosh(tower_eta)>0.1) m_ohcTowEta.push_back(tower_eta);

    }
    //mean_eta = std::accumulate(m_ohcTowEta.begin(), m_ohcTowEta.end(), 0.0) / m_ohcTowEta.size();
    //if (fabs(mean_eta)>0.1) event[2].resize(0);
    //if (fabs(mean_eta)>0.1) event[6].resize(0);
    //std::cout<<std::accumulate(m_ohcTowEta.begin(), m_ohcTowEta.end(), 0.0) <<" "<<m_ohcTowEta.size()<<std::endl;

    tower_range = ihcTowerContainer->size();
    for (unsigned int iter = 0; iter < tower_range; iter++)
    {
      unsigned int towerkey = ihcTowerContainer->encode_key(iter);
      unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
      unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
      double energy = ihcTowerContainer->get_tower_at_channel(iter)->get_energy();
      float chi2 = ihcTowerContainer->get_tower_at_channel(iter)->get_chi2();

      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
                  float tower_phi = tower_geomIN->get_tower_geometry(key)->get_phi();
                  float tower_eta = tower_geomIN->get_tower_geometry(key)->get_eta();
        particle.SetPtEtaPhi(energy,tower_eta,tower_phi);
        if (chi2>0 && energy/cosh(tower_eta)>0.03) event[4].push_back(particle);
        if (chi2>0 && energy/cosh(tower_eta)>0.03) event[9].push_back(particle);

    }
  }
 
 //commenting out everything except clusters 
  //m_T->Fill();

  //if (event[2].size()==0)     return Fun4AllReturnCodes::ABORTEVENT;
  //std::cout<<mean_eta<<std::endl;
  //if (std::isnan(fabs(mean_eta)) || fabs(mean_eta)>1.1) return Fun4AllReturnCodes::ABORTEVENT;  

  //if (event[0].size()<20)     return Fun4AllReturnCodes::ABORTEVENT;
  
  //std::cout<<" Multiplicity = "<<event[0].size()<<std::endl; 
//Filling signal histograms for 2pc
  double trigEff = 1.0; 
  double assEff = 1.0;
  int Ntrig, Nass;
  double dEta=0, dPhi=0;
   for (int i_trigbin=0; i_trigbin<ptTbins_; i_trigbin++) {
     Ntrig = event[i_trigbin].size();     
     for (int j=0; j<Ntrig; j++) {
       trigger = event[i_trigbin][j];
       triggerH[i_trigbin]->Fill(trigger.Eta());
       triggerH2D[i_trigbin]->Fill(trigger.Eta(),trigger.Phi());
       //if (doEffCorr_) trigEff=getEff(trigger,occ);

       for (int i_assbin=0; i_assbin<ptNbins_; i_assbin++) {
         Nass = event[ptTbins_+i_assbin].size();
         for (int i=0; i<Nass; i++) {
           associated = event[ptTbins_+i_assbin][i];
           //if (doEffCorr_) assEff=getEff(associated,occ);
	   if (trigger.Pt()<=associated.Pt()) continue;
           dEta=deltaEta(trigger,associated); dPhi=deltaPhi(trigger,associated);
           if (dEta != 0  && dPhi != 0) signal[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
           if (dEta != 0  && dPhi != 0 && doSymmetrisation_){
              dEta=deltaEta(trigger,associated); dPhi=deltaPhi(associated,trigger);
              signal[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);       
              dEta=deltaEta(associated,trigger); dPhi=deltaPhi(associated,trigger);
              signal[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
              dEta=deltaEta(associated,trigger); dPhi=deltaPhi(trigger,associated);
              signal[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
           } 


         }
       }
     }
   }

   int vtxBin;
   vtxBin=(50+m_vertex)/10;
   //vtxBin=0;
   //std::cout<<"vtxBin = "<<vtxBin<<std::endl; 
   //vtxBin=(mean_eta+1.1)*5;
   //std::cout<<"mean eta = "<<mean_eta<<" (mean_eta+1.1)*5 = "<<(mean_eta+1.1)*5<<" Vtx bin = "<<vtxBin<<std::endl;
   //vtxBin = (vsorted[0].z() +15)/2;
   //if (!useCentrality_) centBin=int (multiplicity-500)/125;   
   //fullSet[15*centBin+vtxBin].push_back(event);
   fullSet[vtxBin].push_back(event);

/*
   if (doBackground_ && fullSet[0].size()==2000) {

      std::cout<<"background start = "<<std::endl;

   int Ntrig, Nass;
   int Nbgev, Nev;// j_bgev;
   double trigEff=1.0, assEff=1.0; double dEta, dPhi;
   int jev; //occ, centBin;

   for (int i_evtClass = 0; i_evtClass < 12; i_evtClass++) {
     Nev = fullSet[i_evtClass].size();
     Nbgev = 10;
     if ( Nev<11 ) Nbgev = Nev-1;

     for (int i_ev=0; i_ev<Nev; i_ev++) {
       event=fullSet[i_evtClass][i_ev];
       eventsBg=fullSet[i_evtClass];
       eventsBg.erase(eventsBg.begin() + i_ev);

       for (int i_bgev=0; i_bgev<Nbgev; i_bgev++)  {

         jev=i_bgev;//gRandom->Rndm()*eventsBg.size();         
         for (int i_trigbin=0; i_trigbin<ptTbins_; i_trigbin++) {
           Ntrig = event[i_trigbin].size();

           for (int i_assbin=0; i_assbin<ptNbins_; i_assbin++) {
             Nass = eventsBg[jev][ptTbins_+i_assbin].size();
             for (int j=0; j<Ntrig; j++) {
               trigger = event[i_trigbin][j];
               for (int i=0; i<Nass; i++) {
                 associated = eventsBg[jev][ptTbins_+i_assbin][i];
                 dEta=deltaEta(trigger,associated); dPhi=deltaPhi(trigger,associated);
                 if (dEta != 0  && dPhi != 0) background[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
                 if (dEta != 0  && dPhi != 0 && doSymmetrisation_){
                   dEta=deltaEta(trigger,associated); dPhi=deltaPhi(associated,trigger);
                   background[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
                   dEta=deltaEta(associated,trigger); dPhi=deltaPhi(associated,trigger);
                   background[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
                   dEta=deltaEta(associated,trigger); dPhi=deltaPhi(trigger,associated);
                   background[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
                 }
               }
             }
           }

         }

         eventsBg.erase (eventsBg.begin()+jev);
       }
     }
    }
    std::cout<<"cleaning background pool "<<std::endl;
    fullSet[0].resize(0);
   } //end background loop
*/

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetValidation::ResetEvent(PHCompositeNode *topNode)
{
  //std::cout << "JetValidation::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  m_ohcTowEta.clear();

  m_id.clear();
  m_nComponent.clear();
  m_eta.clear();
  m_phi.clear();
  m_e.clear();
  m_pt.clear();
  m_unsub_pt.clear();
  m_sub_et.clear();

  m_truthID.clear();
  m_truthNComponent.clear();
  m_truthEta.clear();
  m_truthPhi.clear();
  m_truthE.clear();
  m_truthPt.clear();
  m_truthdR.clear();

  m_eta_subseed.clear();
  m_phi_subseed.clear();
  m_e_subseed.clear();
  m_pt_subseed.clear();
  m_subseed_cut.clear();

  m_eta_rawseed.clear();
  m_phi_rawseed.clear();
  m_e_rawseed.clear();
  m_pt_rawseed.clear();
  m_rawseed_cut.clear();
  
  m_triggerVector.clear();

  m_cle.clear();
  m_clecore.clear();
  m_cleta.clear();
  m_clphi.clear();
  m_clpt.clear();
  m_clprob.clear();

  m_eta_tower.clear();
  m_eta_itower.clear();
  m_eta_retower.clear();
  m_eta_iretower.clear();
  m_phi_tower.clear();
  m_phi_itower.clear();
  m_phi_retower.clear();
  m_phi_iretower.clear();


  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetValidation::EndRun(const int runnumber)
{
  std::cout << "JetValidation::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetValidation::End(PHCompositeNode *topNode)
{
  std::cout << "Mixed event background"<< std::endl;

   doBackground_= false;

   if (doBackground_) {


   std::cout<<"full size = "<<fullSet[0].size()<<std::endl;
   int Ntrig, Nass;
   int Nbgev, Nev;// j_bgev;
   double trigEff=1.0, assEff=1.0; double dEta, dPhi;
   int jev; //occ, centBin;

   for (int i_evtClass = 0; i_evtClass < 10; i_evtClass++) {
     Nev = fullSet[i_evtClass].size();
     Nbgev = 10;
   //std::cout<<__LINE__<<std::endl;
     //occ = i_evtClass/3;
     //centBin = i_evtClass/15;
     if ( Nev<11 ) Nbgev = Nev-1;

     for (int i_ev=0; i_ev<Nev; i_ev++) {
       event=fullSet[i_evtClass][i_ev];
       eventsBg=fullSet[i_evtClass];
       eventsBg.erase(eventsBg.begin() + i_ev);
       //std::cout<<"bg size = "<<eventsBg.size()<<std::endl;      
       for (int i_bgev=0; i_bgev<Nbgev; i_bgev++)  {
   //std::cout<<__LINE__<<std::endl;

	 jev=i_bgev;//gRandom->Rndm()*eventsBg.size();         
         if (int(eventsBg.size()) < jev+1) continue; //condition for less than 10 mixed events
	//std::cout<<__LINE__<<" eventsBg size "<<eventsBg.size()<<" jev =  "<<jev<<std::endl;
         for (int i_trigbin=0; i_trigbin<ptTbins_; i_trigbin++) {
           Ntrig = event[i_trigbin].size();
   //std::cout<<__LINE__<<std::endl;

           for (int i_assbin=0; i_assbin<ptNbins_; i_assbin++) {
             //std::cout<<" i_ass bin = "<<ptTbins_+i_assbin<<std::endl;
             //std::cout<<" eventsBg size "<<eventsBg.size()<<" jev =  "<<jev<<std::endl;

             Nass = eventsBg[jev][ptTbins_+i_assbin].size();
   //std::cout<<__LINE__<<std::endl;
             for (int j=0; j<Ntrig; j++) {
               trigger = event[i_trigbin][j];
   //std::cout<<__LINE__<<std::endl;
   //std::cout<<"Nass = "<<Nass<<std::endl;
               for (int i=0; i<Nass; i++) {
                 associated = eventsBg[jev][ptTbins_+i_assbin][i];
                 //if (doEffCorr_) trigEff=getEff(trigger,occ);
                 //if (doEffCorr_) assEff=getEff(associated,occ);
                 if (trigger.Pt()<=associated.Pt()) continue;

                 dEta=deltaEta(trigger,associated); dPhi=deltaPhi(trigger,associated);
   //std::cout<<__LINE__<<std::endl;
                 if (dEta != 0  && dPhi != 0) background[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
                 if (dEta != 0  && dPhi != 0 && doSymmetrisation_){
                   dEta=deltaEta(trigger,associated); dPhi=deltaPhi(associated,trigger);
                   background[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
                   dEta=deltaEta(associated,trigger); dPhi=deltaPhi(associated,trigger);
                   background[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
                   dEta=deltaEta(associated,trigger); dPhi=deltaPhi(trigger,associated);
                   background[i_trigbin][i_assbin]->Fill(dEta,dPhi,1.0/trigEff/assEff);
                 }
               }
             }
           }

         }

	 eventsBg.erase (eventsBg.begin()+jev); 
       }
     }
    }

   } //end background loop

  std::cout << "JetValidation::End - Output to " << m_outputFileName << std::endl;
  PHTFileServer::get().cd(m_outputFileName);

  for (int i=0; i<ptTbins_; i++){
    triggerH[i]->Write();
    triggerH2D[i]->Write();
    for (int ii=0; ii<ptNbins_; ii++){
        signal[i][ii]->Write();
        background[i][ii]->Write();
    }
  } 
  m_T->Write();
  std::cout << "JetValidation::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetValidation::Reset(PHCompositeNode *topNode)
{
  std::cout << "JetValidation::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void JetValidation::Print(const std::string &what) const
{
  std::cout << "JetValidation::Print(const std::string &what) const Printing info for " << what << std::endl;
}

double JetValidation::deltaEta(const TVector3 a, const TVector3 b)
{

return a.Eta() - b.Eta();

}

double JetValidation::deltaPhi(const TVector3 a, const TVector3 b)
{

  double Dphi;

  Dphi = a.Phi() - b.Phi();
  while (Dphi < -0.5*TMath::Pi()) Dphi = Dphi + 2.0*TMath::Pi();
  while (Dphi > 1.5 *TMath::Pi()) Dphi = Dphi - 2.0*TMath::Pi();
  return Dphi;
}
    
