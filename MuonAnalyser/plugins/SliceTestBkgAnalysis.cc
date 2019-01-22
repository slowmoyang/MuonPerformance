// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
//#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
//#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
//#include "DataFormats/GEMDigi/interface/GEMAMCdataCollection.h"
//#include "DataFormats/GEMDigi/interface/GEMGEBStatusDigiCollection.h"
//#include "DataFormats/GEMDigi/interface/GEMVfatStatusDigiCollection.h"

//Libraries for Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/CommonTopologies/interface/RadialStripTopology.h"

//Libraries for tracking
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class SliceTestBkgAnalysis : public edm::EDAnalyzer {
public:
  explicit SliceTestBkgAnalysis(const edm::ParameterSet&);
  ~SliceTestBkgAnalysis();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
  virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
  
  const GEMEtaPartition* findEtaPartition(const GEMChamber*& chamber, GlobalPoint& tsosGP);

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  //edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;
  //edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  //edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  //edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalers_;
  //edm::EDGetTokenT<GEMAMCdataCollection> amcData_;
  //edm::EDGetTokenT<GEMGEBStatusDigiCollection> gebStatusCol_;
  //edm::EDGetTokenT<GEMVfatStatusDigiCollection> vfatStatusCol_;
  edm::Service<TFileService> fs;
  //edm::EDGetTokenT<CSCRecHit2DCollection> cscRecHits_;
  //edm::EDGetTokenT<GEMVfatStatusDigiCollection> gemDigisvfat_;
  
  MuonServiceProxy* theService_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_; 
  //edm::Handle<GEMAMCdataCollection> amcData;
  //edm::Handle<GEMGEBStatusDigiCollection> gebStatusCol;  
  //edm::Handle<GEMVfatStatusDigiCollection> vfatStatusCol;  

  TH2D* h_allStrips[36][2];
  TH2D* h_allStrips_[36][2];

  TH2D* h_activeLumi;
  TH2D* h_activeLumi_;
  TH1D* h_lumiStatus;
  
  TH1D* h_clusterSize, *h_totalStrips, *h_bxtotal;

  TH1D* h_res_x, *h_res_y, *h_pull_x, *h_pull_y;

  TTree *t_setup;
  float b_stripLength[2][8];
  float b_stripPitch[2][8];
  float b_etaPosR[2][8];
  //float b_cscArea;
  
  TTree *t_run;
  
  TTree *t_event;
  int b_run, b_lumi, b_latency;
  int b_nGEMHits, b_nCSCHits, b_nMuons;
  float b_instLumi, b_pileup;
  unsigned int b_timeLow, b_timeHigh;
  int b_vfatQuality[36][2][24];

  TTree *t_hit;
  int b_firstStrip, b_nStrips, b_chamber, b_layer, b_etaPartition;
  float b_x, b_y, b_z, b_eta, b_mag;

  TTree *t_muon;
  bool b_isValid;
  int b_id;
  float b_pt, b_phi;
  float b_resX, b_resY, b_resPhi;
  float b_pullX, b_pullY;
  
  TTree *t_csc;

};

SliceTestBkgAnalysis::SliceTestBkgAnalysis(const edm::ParameterSet& iConfig)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  //gemDigisvfat_ = consumes<GEMVfatStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigisvfat"));
  //vertexCollection_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  //amcData_ = consumes<GEMAMCdataCollection>(iConfig.getParameter<edm::InputTag>("amcData"));
  //gebStatusCol_ = consumes<GEMGEBStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("gebStatusCol"));
  //vfatStatusCol_ = consumes<GEMVfatStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("vfatStatusCol"));
  lumiScalers_ = consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);

  h_clusterSize=fs->make<TH1D>(Form("clusterSize"),"clusterSize",100,0,100);
  h_totalStrips=fs->make<TH1D>(Form("totalStrips"),"totalStrips",200,0,200);
  h_bxtotal=fs->make<TH1D>(Form("bx"),"bx",1000,0,1000);


  t_run = fs->make<TTree>("Run", "Run");
  t_run->Branch("run", &b_run, "run/I");

  t_setup = fs->make<TTree>("Setup", "Setup");
  t_setup->Branch("stripLength", &b_stripLength, "stripLength[2][8]/F");
  t_setup->Branch("stripPitch", &b_stripPitch, "stripPitch[2][8]/F");
  t_setup->Branch("etaPosR", &b_etaPosR, "etaPosR[2][8]/F");
  //t_setup->Branch("cscArea", &b_cscArea, "cscArea/F");

  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  t_event->Branch("nCSCHits", &b_nCSCHits, "nCSCHits/I");
  t_event->Branch("nMuons", &b_nMuons, "nMuons/I");
  t_event->Branch("run", &b_run, "run/I");
  t_event->Branch("lumi", &b_lumi, "lumi/I");
  t_event->Branch("latency", &b_latency, "latency/I");
  t_event->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_event->Branch("timeLow", &b_timeLow, "timeLow/i");
  t_event->Branch("timeHigh", &b_timeHigh, "timeHigh/i");
  t_event->Branch("vfatQuality", b_vfatQuality, "vfatQuality[36][2][24]/I");

  t_hit = fs->make<TTree>("Hit", "Hit");
  t_hit->Branch("firstStrip", &b_firstStrip, "firstStrip/I");
  t_hit->Branch("nStrips", &b_nStrips, "nStrips/I");
  t_hit->Branch("chamber", &b_chamber, "chamber/I");
  t_hit->Branch("layer", &b_layer, "layer/I");
  t_hit->Branch("etaPartition", &b_etaPartition, "etaPartition/I");
  t_hit->Branch("x", &b_x, "x/F");
  t_hit->Branch("y", &b_y, "y/F");
  t_hit->Branch("z", &b_z, "z/F");
  t_hit->Branch("eta", &b_eta, "eta/F");
  t_hit->Branch("mag", &b_mag, "mag/F");

  //t_csc = fs->make<TTree>("CSC", "CSC");
  //t_csc->Branch("chamber", &b_chamber, "chamber/I");
  //t_csc->Branch("layer", &b_layer, "layer/I");

  //t_muon = fs->make<TTree>("muon", "muon");
  //t_muon->Branch("firstStrip", &b_firstStrip, "firstStrip/I");
  //t_muon->Branch("nStrips", &b_nStrips, "nStrips/I");
  //t_muon->Branch("chamber", &b_chamber, "chamber/I");
  //t_muon->Branch("layer", &b_layer, "layer/I");
  //t_muon->Branch("etaPartition", &b_etaPartition, "etaPartition/I");
  //t_muon->Branch("id", &b_id, "id/I");
  //t_muon->Branch("pt", &b_pt, "pt/F");
  //t_muon->Branch("eta", &b_eta, "eta/F");
  //t_muon->Branch("phi", &b_phi, "phi/F");
  //t_muon->Branch("resX", &b_resX, "resX/F");
  //t_muon->Branch("resY", &b_resY, "resY/F");
  //t_muon->Branch("resPhi", &b_resPhi, "resPhi/F");
  //t_muon->Branch("isValid", &b_isValid, "isValid/O");


}

SliceTestBkgAnalysis::~SliceTestBkgAnalysis()
{
  t_setup->Fill();
}

void
SliceTestBkgAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_run = iEvent.run();
  b_lumi = iEvent.luminosityBlock();
  
  b_nGEMHits = 0;
  b_nCSCHits = 0;
  b_nMuons = 0;
  
  for (int ch = 0; ch < 36; ch++) {
    for (int ly = 0; ly < 2; ly++) {
      for (int vfat = 0; vfat < 24; vfat++) {
        b_vfatQuality[ch][ly][vfat] = -1;
      }
    }
  }
  
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;
  
  //edm::ESHandle<CSCGeometry> hCSCGeom;
  //iSetup.get<MuonGeometryRecord>().get(hCSCGeom);
  //const CSCGeometry* CSCGeometry_ = &*hCSCGeom;
  
  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);

  //edm::Handle<CSCSegmentCollection> cscSegments;  
  //iEvent.getByToken(cscSegments_, cscSegments);

  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  //theService_->update(iSetup);
  //auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  
  //edm::Handle<reco::VertexCollection> vertexCollection;
  //iEvent.getByToken( vertexCollection_, vertexCollection );
  //if(vertexCollection.isValid()) {
  //  vertexCollection->size();
  //}

  edm::Handle<LumiScalersCollection> lumiScalers;
  iEvent.getByToken( lumiScalers_, lumiScalers );

  //Handle<View<reco::Muon> > muons;
  //iEvent.getByToken(muons_, muons);

  //edm::Handle<GEMVfatStatusDigiCollection> gemDigisvfat;
  //iEvent.getByToken(gemDigisvfat_, gemDigisvfat);
  
  //edm::Handle<CSCRecHit2DCollection> cscRecHits;  
  //iEvent.getByToken(cscRecHits_, cscRecHits);
  //iEvent.getByToken(amcData_, amcData);
  //iEvent.getByToken(gebStatusCol_, gebStatusCol);
  //iEvent.getByToken(vfatStatusCol_, vfatStatusCol);

  int totalStrips = 0;
  auto instLumi = (lumiScalers->at(0)).instantLumi()/10000;
  h_lumiStatus->Fill(instLumi);
  b_instLumi = instLumi;
  b_pileup = (lumiScalers->at(0)).pileup();
  b_timeHigh = iEvent.time().unixTime();
  b_timeLow = iEvent.time().microsecondOffset();

  b_latency = -1;

  //for (auto g : *amcData) {
  //  for (auto a = g.second.first; a != g.second.second; ++a) {
  //    if (b_latency != -1 && b_latency != a->param1())
  //      std::cout << "CHANGING LATENCY - old: " << b_latency << " new: " << a->param1() << std::endl;
  //    b_latency = a->param1();
  //    if (b_latency == -1) std::cout << "-1 LATENCY VALUE - " << iEvent.id().event() << " " << iEvent.id().run() << std::endl;
  //  }
  //}

  //for (size_t i = 0; i < muons->size(); ++i) {
  //  edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
  //  const reco::Muon* mu = muRef.get();

  //  const reco::Track* muonTrack = 0;  
  //  if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
  //  else if ( mu->outerTrack().isNonnull() ) muonTrack = mu->outerTrack().get();
  //  if (muonTrack) {
  //    b_pt = muonTrack->pt();
  //    b_eta = muonTrack->eta();
  //    b_phi = muonTrack->phi();
  //    if (b_eta < -2.15 or b_eta > -1.6) continue;
  //    t_muon->Fill();
  //    b_nMuons++;
  //  }
  //}
  
  for (auto ch : GEMGeometry_->chambers()) {
    for(auto roll : ch->etaPartitions()) {
      GEMDetId rId = roll->id();

      auto recHitsRange = gemRecHits->get(rId); 
      auto gemRecHit = recHitsRange.first;
      b_chamber = rId.chamber();
      b_layer = rId.layer();
      b_etaPartition = rId.roll();
      
      //auto vfats = vfatStatusCol->get(rId); 
      //int index = 0;
      //for (auto vfat = vfats.first; vfat != vfats.second; ++vfat) {
      //  int vfatNu = (8-b_etaPartition)+8*(index);
      //  b_vfatQuality[b_chamber-1][b_layer-1][vfatNu] = vfat->quality();
      //  index++;
      //}
      
      for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {

        h_clusterSize->Fill(hit->clusterSize());
        h_bxtotal->Fill(hit->BunchX());

        for (int nstrip = hit->firstClusterStrip(); nstrip < hit->firstClusterStrip()+hit->clusterSize(); ++nstrip) {
          totalStrips++;
          h_allStrips[rId.chamber()-1][rId.layer()-1]->Fill(nstrip, rId.roll());
          h_allStrips_[rId.chamber()-1][rId.layer()-1]->Fill(nstrip, rId.roll());
        }

        b_firstStrip = hit->firstClusterStrip();
        b_nStrips = hit->clusterSize();

        int vfat = (8-b_etaPartition)+8*((b_firstStrip-1)/128);
        float vfatNu = b_chamber+(b_layer-1)/2.+vfat/48.;
        h_activeLumi->Fill(b_lumi, vfatNu);

        auto globalPosition = roll->toGlobal(hit->localPosition());
        b_x = globalPosition.x();
        b_y = globalPosition.y();
        b_z = globalPosition.z();
        b_eta = globalPosition.eta();
        b_mag = globalPosition.mag();

        t_hit->Fill();
        b_nGEMHits++;
      }

      //auto gemVfatRange = gemDigisvfat->get(rId);
      //
      //for (auto vfatStatus = gemVfatRange.first; vfatStatus != gemVfatRange.second; vfatStatus++) {
      //  auto ix = vfatStatus->position();
      //  int vfat = 8-b_etaPartition+8*ix;
      //  bool flag = int(vfatStatus->flag()) != 0;
      //  flag = int(vfatStatus->quality()) != 0;
      //  if (!flag) h_activeLumi_->Fill(b_lumi, b_chamber+(b_layer-1)/2.+vfat/48.);
      //}
    }
  }

  //for (size_t i = 0; i < muons->size(); ++i) {

  //  edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
  //  const reco::Muon* mu = muRef.get();

  //  // muon id
  //  int muonId = 0;
  //  if (mu->passed(reco::Muon::Selector::CutBasedIdTight)) muonId = 2;
  //  else if (mu->passed(reco::Muon::Selector::CutBasedIdLoose)) muonId = 1;

  //  // tight and pt > 20 muon only
  //  //if (muonId != 2) continue;
  //  //if (mu->pt() < 20) continue;
  //  
  //  const reco::Track* muonTrack = 0;  
  //  if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
  //  else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
  //  if (!muonTrack) continue;
  //  else {
  //    b_id = muonId;
  //    b_pt = muonTrack->pt();
  //    b_eta = muonTrack->eta();
  //    b_phi = muonTrack->phi();
  //  }

  //  reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
  //  for (auto chamber : GEMGeometry_->chambers()) {
  //    b_isValid = false;
  //    if (chamber->id().chamber() == 1) continue; // ignore chammber 1
  //    if (mu->eta() * chamber->id().region() < 0 ) continue;

  //    TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),chamber->surface());
  //    if (!tsos.isValid()) continue;

  //    GlobalPoint tsosGP = tsos.globalPosition();
  //    auto etaPart =  findEtaPartition(chamber, tsosGP);
  //    if (!etaPart) continue;

  //    auto gemid = etaPart->id();
  //    auto locPos = etaPart->toLocal(tsosGP);
  //    auto strip = (int) etaPart->strip(locPos);
  //    
  //    b_chamber = gemid.chamber();
  //    b_layer = gemid.layer();
  //    b_etaPartition = gemid.roll();

  //    //Find hit
  //    b_resX = 999;
  //    GEMRecHit closestHit;
  //    auto recHitsRange = gemRecHits->get(gemid);
  //    for (auto hit = recHitsRange.first; hit != recHitsRange.second; ++hit) {
  //      LocalPoint hitLocPos = hit->localPosition();
  //      if ( fabs(hitLocPos.x() - locPos.x()) < fabs(b_resX) ) {
  //        b_resX = hitLocPos.x() - locPos.x();
  //        closestHit = (*hit);
  //      }
  //    }
  //    if (b_resX != 999) { //We could not find hit associated with muon 
  //      auto hitLocPos = closestHit.localPosition();
  //      auto hitGlobPos = etaPart->toGlobal(hitLocPos);
  //      b_firstStrip = closestHit.firstClusterStrip();
  //      b_nStrips = closestHit.clusterSize();
  //      b_resY = hitLocPos.y() - locPos.y();
  //      b_resPhi = hitGlobPos.phi() - tsosGP.phi();

  //      LocalError && locErr = tsos.localError().positionError();
  //      LocalError && hitLocErr = closestHit.localPositionError();
  //      b_pullX = b_resX / std::sqrt(hitLocErr.xx() + locErr.xx());
  //      b_pullY = b_resY / std::sqrt(hitLocErr.yy() + locErr.yy());

  //      b_isValid = true;
  //    }
  //    
  //    t_muon->Fill();
  //    b_nMuons++;
  //  }
  //}
  //for (auto ch : CSCGeometry_->chambers()) {
  //  CSCDetId cId = ch->id();
  //  // Selection for CSC chambers are in same regime with GEMINI
  //  if (cId.station() != 1) continue;
  //  if (cId.endcap() != 2) continue;
  //  if (cId.chamber() < 27 or cId.chamber() > 30) continue;

  //  auto segmentsRange = cscSegments->get(cId);
  //  auto cscSegment = segmentsRange.first;
  //  
  //  for (auto hit = cscSegment; hit != segmentsRange.second; ++hit) {
  //    b_chamber = cId.chamber();

  //    t_csc->Fill();
  //    b_nCSCHits++;
  //  }
  //}

  h_totalStrips->Fill(totalStrips);

  t_event->Fill();
}

const GEMEtaPartition* SliceTestBkgAnalysis::findEtaPartition(const GEMChamber*& chamber, GlobalPoint& tsosGP){
  for (auto etaPart : chamber->etaPartitions()) {
    const LocalPoint locPos = etaPart->toLocal(tsosGP);
    const LocalPoint locPos2D(locPos.x(), locPos.y(), 0);
    const BoundPlane& bps(etaPart->surface());
    if (!bps.bounds().inside(locPos2D)) continue;
    return etaPart;
  }
  return nullptr;
}

void SliceTestBkgAnalysis::beginJob(){}
void SliceTestBkgAnalysis::endJob(){}

void SliceTestBkgAnalysis::beginRun(Run const& run, EventSetup const& iSetup){
  for (int ichamber=26; ichamber<30;++ichamber) {
  // for (int ichamber=27; ichamber<=30;++ichamber) {
    for (int ilayer=0; ilayer<2;++ilayer) {
      h_allStrips_[ichamber][ilayer] = fs->make<TH2D>(Form("%i ch %i lay %i", run.run(), ichamber+1, ilayer+1),"",384,0.5,384.5,8,0.5,8.5);
      h_allStrips_[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_allStrips_[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");
    }
  }


  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;
  
  //edm::ESHandle<CSCGeometry> hCSCGeom;
  //iSetup.get<MuonGeometryRecord>().get(hCSCGeom);
  //const CSCGeometry* CSCGeometry_ = &*hCSCGeom;

  h_activeLumi = fs->make<TH2D>(Form("%i active lumi", run.run()),Form("Run number %i", run.run()),3000, 0, 3000, 400, 27, 31);
  h_activeLumi_ = fs->make<TH2D>(Form("%i active lumi_", run.run()),Form("Run number %i", run.run()),3000, 0, 3000, 400, 27, 31);
  
  for (auto ch : GEMGeometry_->chambers()) {
    for(auto roll : ch->etaPartitions()) {
      GEMDetId rId = roll->id();

      const TrapezoidalStripTopology* top_(dynamic_cast<const TrapezoidalStripTopology*>(&(roll->topology())));
      const float striplength(top_->stripLength());
      const float pitch(roll->pitch());
      auto centrePosition = roll->toGlobal(roll->centreOfStrip(192));

      b_stripLength[rId.chamber()%2][rId.roll()-1] = striplength;
      b_stripPitch[rId.chamber()%2][rId.roll()-1] = pitch;
      b_etaPosR[rId.chamber()%2][rId.roll()-1] = sqrt(pow(centrePosition.x(), 2) + pow(centrePosition.y(), 2));
    }
  }

  //for (auto ch : CSCGeometry_->chambers()) {
  //  for (auto layer : ch->layers()) {
  //    CSCDetId lId = layer->id();
  //    if (lId.station() != 1) continue;
  //    if (lId.endcap() != 2) continue;
  //    if (lId.chamber() < 27 or lId.chamber() > 30) continue;

  //    const RadialStripTopology* top_(dynamic_cast<const RadialStripTopology*>(&(layer->topology())));
  //    const float striplength(top_->stripLength());
  //    const float pitch(layer->geometry()->stripPitch());
  //    b_cscArea = striplength * pitch * top_->nstrips();
  //  }
  //}
  
  b_run = run.run();
  t_run->Fill();
}
void SliceTestBkgAnalysis::endRun(Run const&, EventSetup const&){}

void SliceTestBkgAnalysis::beginLuminosityBlock(LuminosityBlock const& lumiBlock, EventSetup const& iSetup){
  for (int ichamber=26; ichamber<30;++ichamber) {
  // for (int ichamber=27; ichamber<=30;++ichamber) {
    for (int ilayer=0; ilayer<2;++ilayer) {
      h_allStrips[ichamber][ilayer] = fs->make<TH2D>(Form("%i %i ch %i lay %i",lumiBlock.run(), lumiBlock.luminosityBlock(), ichamber+1, ilayer+1),"",384,0.5,384.5,8,0.5,8.5);
      h_allStrips[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_allStrips[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");
    }
  }

  h_lumiStatus = fs->make<TH1D>(Form("%i %i status", lumiBlock.run(), lumiBlock.luminosityBlock()), "", 1000, 0, 2);
}
void SliceTestBkgAnalysis::endLuminosityBlock(LuminosityBlock const& lumiBlock, EventSetup const& iSetup){}

//define this as a plug-in
DEFINE_FWK_MODULE(SliceTestBkgAnalysis);
