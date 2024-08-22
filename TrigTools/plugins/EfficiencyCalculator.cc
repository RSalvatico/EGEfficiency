#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/EgammaObject.h"
#include "L1Trigger/L1TGlobal/plugins/L1TGlobalProducer.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "HLTrigger/Egamma/plugins/HLTEgammaL1TMatchFilterRegional.cc"
#include "CondFormats/DataRecord/interface/L1TGlobalParametersRcd.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include <vector>
#include <string>
#include <iostream>
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"
#include <stdlib.h>

#define TWOPI 6.283185308
int globCounter = 0;
/****************************************************************************
 Simple example on how to access the L1 decision via  l1t::L1TGlobalUtil
 To make things easier it uses HLTPrescaleProvider to obtain the  
 l1t::L1TGlobalUtil object

 note, its very important that you run with the correct global tag as it
 retrieves the menu from this

 

 author Sam  Harper (RAL), 2019
*****************************************************************************/

class EfficiencyCalculator : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
 
private:
  
  std::string hltProcess_; //name of HLT process, usually "HLT"
  edm::EDGetTokenT<std::vector<pat::Electron> > eleToken_;
  edm::EDGetTokenT<edm::TriggerResults > hltToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsToken_;
  
  //edm::EDGetTokenT<trigger::TriggerEvent > hltsevtToken_;
  //edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  //edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs > L1Token_;
  //edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::Service<TFileService> fs;
  double endcap_end_ = 2.5;
  TH1D* num_ele_pt_EB;
  TH1D* num_ele_pt_EE;
  TH1D* num_ele_eta;
  TH1D* num_ele_phi;
  TH1D* num_ele_pu;

  TH1D* den_ele_pt_EB;
  TH1D* den_ele_pt_EE;
  TH1D* den_ele_eta;
  TH1D* den_ele_phi;
  TH1D* den_ele_pu;
  float barrel_end_ = 1.4442;
  TH2D* occupancy_phi_eta_all;

public:
  explicit EfficiencyCalculator(const edm::ParameterSet& iConfig);
  ~EfficiencyCalculator(){}
  
 private:
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& iSetup);
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override{}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob(){}
};


EfficiencyCalculator::EfficiencyCalculator(const edm::ParameterSet& iConfig):
  hltProcess_("HLTX")
{
  eleToken_     = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons","","RECO"));
  hltToken_     = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  //hltsevtToken_ = consumes<trigger::TriggerEvent>(edm::InputTag("hltTriggerSummaryAOD","","HLTX"));
  triggerObjectsToken_ = consumes<std::vector<pat::TriggerObjectStandAlone> > (edm::InputTag("slimmedPatTrigger","","RECO"));
  
  // pT
  num_ele_pt_EB = fs->make<TH1D>("num_ele_pt_EB",";pT (GeV);Events",10,0,100);
  num_ele_pt_EE = fs->make<TH1D>("num_ele_pt_EE",";pT (GeV);Events",10,0,100);

  den_ele_pt_EB = fs->make<TH1D>("den_ele_pt_EB",";pT (GeV);Events",10,0,100);
  den_ele_pt_EE = fs->make<TH1D>("den_ele_pt_EE",";pT (GeV);Events",10,0,100);

  // eta
  num_ele_eta = fs->make<TH1D>("num_ele_eta",";eta;Events",53,-2.65,2.65);
  den_ele_eta = fs->make<TH1D>("den_ele_eta",";eta;Events",53,-2.65,2.65);

  // phi
  num_ele_phi = fs->make<TH1D>("num_ele_phi",";phi;Events",63,-3.15,3.15);
  den_ele_phi = fs->make<TH1D>("den_ele_phi",";phi;Events",63,-3.15,3.15);
  occupancy_phi_eta_all = fs->make<TH2D>("occupancy_phi_eta_all",";#eta_{SC};#phi",50,-2.5,2.5,32,-3.2,3.2);
}


const int getFilterIndex(const trigger::TriggerEvent& trigEvt, const std::string filterName){
  for(auto i=0; i < trigEvt.sizeFilters(); i++){

    //if("hltEle32WPTightGsfTrackIsoFilter"==trigEvt.filterLabel(i))  {
    //std::cout << trigEvt.filterLabel(i) << std::endl;
    //globCounter++;}
    if(filterName==trigEvt.filterLabel(i)) return i;
  }
  return trigEvt.sizeFilters();  
}

std::vector<const trigger::TriggerObject*> getListFilterPassedObj(const std::string filterName,const trigger::TriggerEvent& hlts){
  std::vector<const trigger::TriggerObject*> eg_trig_objs;
  const int filterIndex = getFilterIndex(hlts,filterName);
  if(filterIndex < hlts.sizeFilters() ){
    for(auto filterKey : hlts.filterKeys(filterIndex)){
      auto& obj = hlts.getObjects()[filterKey];
      eg_trig_objs.push_back(&obj);
    }
  }
  return eg_trig_objs;
}


//std::vector<const trigger::TriggerObject*> matchTrigObjs(const float eta,const float phi,const float pT,std::vector<const trigger::TriggerObject*> trigObjs,const float maxDeltaR=0.1, const float maxDpT=0.1)
//{
//  std::vector<const trigger::TriggerObject*> matchedObjs;
//  const float maxDR2 = maxDeltaR*maxDeltaR;
//  for(auto& trigObj : trigObjs){
//    const float deltaR2 = reco::deltaR2(eta, phi, trigObj->eta(), trigObj->phi());
//    if(deltaR2 < maxDR2 && (fabs(pT - trigObj->pt())/pT) < maxDpT) matchedObjs.push_back(trigObj);
//  }
//  return matchedObjs;
//}

//std::vector<const pat::TriggerObjectStandAlone*> matchTrigObjs(const float eta,const float phi,std::vector<const pat::TriggerObjectStandAlone*> trigObjs,const float maxDeltaR=0.1)
//{
//  std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
//  const float maxDR2 = maxDeltaR*maxDeltaR;
//  for(auto& trigObj : trigObjs){
//    const float deltaR2 = reco::deltaR2(eta, phi, trigObj->eta(), trigObj->phi());
//    if(deltaR2 < maxDR2) matchedObjs.push_back(trigObj);
//  }
//  return matchedObjs;
//}

std::vector<pat::TriggerObjectStandAlone> matchTrigObjs(const float eta,const float phi,std::vector<pat::TriggerObjectStandAlone> trigObjs,const float maxDeltaR=0.1)
{
  std::vector<pat::TriggerObjectStandAlone> matchedObjs;
  const float maxDR2 = maxDeltaR*maxDeltaR;
  for(auto trigObj : trigObjs){
    const float deltaR2 = reco::deltaR2(eta, phi, trigObj.eta(), trigObj.phi());
    if(deltaR2 < maxDR2) matchedObjs.push_back(trigObj);
  }
  return matchedObjs;
}

const bool matchDRAndpT(const float eta1,const float phi1,const float pT1,const float eta2,const float phi2,const float pT2,const float maxDeltaR=0.1, const float maxDpT=0.1){
  
  bool isMatched = false;
  const float maxDR2 = maxDeltaR*maxDeltaR;
  const float deltaR2 = reco::deltaR2(eta1, phi1, eta2, phi2);
  const float deltaPt = fabs(pT1-pT2)/pT1;
  if(deltaR2 < maxDR2 && deltaPt < maxDpT) isMatched = true;
  return isMatched;
}

std::vector<const reco::GenParticle*> getGenparts(const std::vector<reco::GenParticle>& genparts,const int pid=11, bool antipart=true, const int status=1){

  std::vector<const reco::GenParticle*> selected;
  if(genparts.empty()) return selected;

  for(auto& part : genparts){
    const int pdg_id = part.pdgId();
    if(pdg_id == pid || (antipart && abs(pdg_id) == abs(pid))){
      if(part.isHardProcess() && status == 1){
	selected.push_back(&part);
      }
    }
  }
  return selected;
}

float matchToGen(const float eta,const float phi,const std::vector<reco::GenParticle>& genparts,const int pid=11,bool antipart=true,const float max_dr=0.1,const int status=1){

  float best_dr2 = max_dr*max_dr;
  float best_pt = -999999;
  auto selected_parts = getGenparts(genparts,pid,antipart,status);
  for(auto part : selected_parts){
    const float dr2 = reco::deltaR2(eta, phi, part->eta(), part->phi());
    if(dr2 < best_dr2){
      best_dr2 = dr2;
      best_pt = part->pt();
    }
  }
  return best_pt;
}


//float calculateInvMass(const trigger::EgammaObject tagElectron, const trigger::EgammaObject probeCandidate) {
//
//  TLorentzVector tag;
//  TLorentzVector probe;
//
//  tag.SetPxPyPzE(tagElectron.px(),tagElectron.py(),tagElectron.pz(),tagElectron.energy());
//  probe.SetPxPyPzE(probeCandidate.px(),probeCandidate.py(),probeCandidate.pz(),probeCandidate.energy());
//
//  float invMass = (tag + probe).M();
//
//  return invMass;
//}

float calculateInvMass(const pat::Electron tagElectron, const pat::Electron probeCandidate) {

  TLorentzVector tag;
  TLorentzVector probe;

  tag.SetPxPyPzE(tagElectron.px(),tagElectron.py(),tagElectron.pz(),tagElectron.energy());
  probe.SetPxPyPzE(probeCandidate.px(),probeCandidate.py(),probeCandidate.pz(),probeCandidate.energy());

  float invMass = (tag + probe).M();

  return invMass;
}



//we need to initalise the menu each run (menu can and will change on run boundaries)
void EfficiencyCalculator::beginRun(const edm::Run& run,const edm::EventSetup& setup)
{
  //bool changed=false;
  //hltPSProv_.init(run,setup,hltProcess_,changed);
}

void EfficiencyCalculator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<std::vector<pat::Electron> > ele;
  iEvent.getByToken(eleToken_,ele);

  edm::Handle<edm::TriggerResults > hlt;
  iEvent.getByToken(hltToken_,hlt);

  //edm::Handle<trigger::TriggerEvent > hltsevt;
  //iEvent.getByToken(hltsevtToken_,hltsevt);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
  iEvent.getByToken(triggerObjectsToken_, triggerObjects);


  //Create a list of trigger objects with unpacked filter names
  std::vector<pat::TriggerObjectStandAlone> unpackedTriggerObjects;
  for(auto& trigObj : *triggerObjects){
    unpackedTriggerObjects.push_back(trigObj);
    unpackedTriggerObjects.back().unpackFilterLabels(iEvent,*hlt);
  }


  //auto eg_trig_objs_filter = getListFilterPassedObj("hltEle32WPTightGsfTrackIsoFilter",*hltsevt.product());

  auto electrons = ele.product();
 
  if(electrons->size() < 2) return;

  std::vector<pat::Electron> listOfTags;



  for(auto& el : *electrons){

    if(fabs(el.eta()) > barrel_end_ || el.pt() < 30. || !(el.electronID("mvaEleID-RunIIIWinter22-iso-V1-wp80"))) continue;

    auto matchedTrigObjsTags = matchTrigObjs(el.eta(),el.phi(),unpackedTriggerObjects);
    
    //for(const auto& trigObj : *matchedTrigObjsTags){
  
    for(auto trigObj : matchedTrigObjsTags){
      if(trigObj.hasFilterLabel("hltEle30WPTightGsfTrackIsoFilter")){
       listOfTags.push_back(el);
      }
    }
  }

  bool isGoodPair = false;
  for(auto& el : *electrons){

    if(fabs(el.eta()) > endcap_end_ || !(el.electronID("mvaEleID-RunIIIWinter22-iso-V1-wp80"))) continue;

    for(auto tag : listOfTags){

      float invMass = calculateInvMass(tag, el);
      isGoodPair = (fabs(invMass - 91.1876) < 30. && tag.pdgId()*el.pdgId() < 0)? true : false;
      if(isGoodPair) break;
    }

    if(!isGoodPair) continue;

    //Fill denominators
    if (fabs(el.eta()) < 1.5 ) den_ele_pt_EB->Fill(el.pt());
    else den_ele_pt_EE->Fill(el.pt());

    if (el.pt() > 30.) {
    	den_ele_eta->Fill(el.eta());
    	den_ele_phi->Fill(el.phi());
    }

    auto matchedTrigObjsProbes = matchTrigObjs(el.eta(),el.phi(),unpackedTriggerObjects);
    //for(const auto& trigObj : *matchedTrigObjsProbes){
    for(auto trigObj : matchedTrigObjsProbes){
      if(trigObj.hasFilterLabel("hltEle30WPTightGsfTrackIsoFilter")){
        if (fabs(el.eta()) < 1.5 ) num_ele_pt_EB->Fill(el.pt());
        else num_ele_pt_EE->Fill(el.pt());

        if (el.pt() > 30.) {
    	    num_ele_eta->Fill(el.eta());
    	    num_ele_phi->Fill(el.phi());
        }
        break; //Avoid to fill the numerator more than once with the same object if more than one offline-online matching is found
      }
    }   
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(EfficiencyCalculator);
