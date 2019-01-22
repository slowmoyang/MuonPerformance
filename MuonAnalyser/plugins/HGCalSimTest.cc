// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include<map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GEMDigi/interface/ME0DigiCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0RecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class HGCalSimTest : public edm::EDAnalyzer {
public:
  explicit HGCalSimTest(const edm::ParameterSet&);
  ~HGCalSimTest();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<ME0DigiCollection> me0Digis_;
  edm::EDGetTokenT<ME0SegmentCollection> me0Segments_;
  edm::Service<TFileService> fs;


  TTree *t_event;
  int b_run, b_lumi, b_latency;
  int b_nME0Digis;

  TTree *t_digi;
  int b_firstStrip, b_nStrips, b_chamber, b_layer, b_etaPartition;
};

HGCalSimTest::HGCalSimTest(const edm::ParameterSet& iConfig)
{ 
  me0Digis_ = consumes<ME0DigiCollection>(iConfig.getParameter<edm::InputTag>("me0Digis"));
  me0Segments_ = consumes<ME0SegmentCollection>(iConfig.getParameter<edm::InputTag>("me0Segments"));

  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nME0Digis", &b_nME0Digis, "nME0Digis/I");

  t_digi = fs->make<TTree>("Hit", "Hit");
  t_digi->Branch("firstStrip", &b_firstStrip, "firstStrip/I");
  t_digi->Branch("nStrips", &b_nStrips, "nStrips/I");
  t_digi->Branch("chamber", &b_chamber, "chamber/I");
  t_digi->Branch("layer", &b_layer, "layer/I");
  t_digi->Branch("etaPartition", &b_etaPartition, "etaPartition/I");
}

HGCalSimTest::~HGCalSimTest()
{
}

void
HGCalSimTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<ME0Geometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  const ME0Geometry* ME0Geometry_ = &*hGeom;
  
  edm::Handle<ME0DigiCollection> me0Digis;  
  iEvent.getByToken(me0Digis_, me0Digis);
  
  edm::Handle<ME0SegmentCollection> me0Segments;  
  iEvent.getByToken(me0Segments_, me0Segments);

  b_nME0Digis = 0;

  for (auto ch : ME0Geometry_->chambers()) {
    for (auto ly : ch->layers()) {
      for (auto roll : ly->etaPartitions()) {
        ME0DetId rId = roll->id();

        auto digisRange = me0Digis->get(rId); 
        auto me0Digi = digisRange.first;
        int roll_ = rId.roll();
        int chamber = rId.chamber();
        int layer = rId.layer();
        for (auto hit = me0Digi; hit != digisRange.second; ++hit) {
          int strip = hit->strip();
          t_digi->Fill();
          b_nME0Digis++;
        }
      }
    }
  }
  t_event->Fill();
}

void HGCalSimTest::beginJob(){}
void HGCalSimTest::endJob(){}

void HGCalSimTest::beginRun(Run const& run, EventSetup const&){
}
void HGCalSimTest::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalSimTest);
