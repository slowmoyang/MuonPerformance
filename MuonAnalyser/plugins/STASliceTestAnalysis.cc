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

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMCStatusDigi.h"
#include "DataFormats/GEMDigi/interface/GEMVfatStatusDigi.h"
#include "DataFormats/GEMDigi/interface/GEMGEBStatusDigi.h"
#include "EventFilter/GEMRawToDigi/interface/AMCdata.h"

#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMC13EventCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMCdataCollection.h"
#include "DataFormats/GEMDigi/interface/GEMGEBStatusDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMVfatStatusDigiCollection.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "TVector2.h"


using namespace std;
using namespace edm;

class STASliceTestAnalysis : public edm::EDAnalyzer {
public:
  explicit STASliceTestAnalysis(const edm::ParameterSet&);
  ~STASliceTestAnalysis();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  //  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<reco::TrackCollection> staTracks_;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalers_;

  edm::EDGetTokenT<GEMAMC13EventCollection> amc13Event_;
  edm::EDGetTokenT<GEMAMCdataCollection> amcData_;
  edm::EDGetTokenT<GEMGEBStatusDigiCollection> gebStatusCol_;
  edm::EDGetTokenT<GEMVfatStatusDigiCollection> vfatStatusCol_;

  bool checkEtaPartitionGood(const GEMEtaPartition* part);
  
  edm::Handle<GEMAMC13EventCollection> amc13Event;
  edm::Handle<GEMAMCdataCollection> amcData;
  edm::Handle<GEMGEBStatusDigiCollection> gebStatusCol;  
  edm::Handle<GEMVfatStatusDigiCollection> vfatStatusCol;  

  ULong64_t b_event;
  int b_run, b_lumi;
  float b_instLumi;
  int b_firstStrip, b_nStrips, b_chamber, b_layer, b_etaPartition, b_muonQuality, b_bx;
  int b_latency;
  vector<int> b_strips;
  float b_x, b_y, b_z;
  float b_mu_eta, b_mu_phi, b_mu_pt;
  float b_pull_x, b_pull_y, b_res_x, b_res_y;

  int nEvents, nMuonTotal, nGEMFiducialMuon, nGEMTrackWithMuon, nGEMTotal;
  int b_nMuons, b_nMuonsWithGEMHit, b_nMuonsInMuonTree, b_nHitsInHitTree;
  int b_valid;

  int b_nGEMHits;

  // muon branches
  int m_nGEMhits, m_nCSChits, m_nhits, m_nvalidhits;
  int m_nbounds;
  int m_quality, m_charge;
  float m_pt, m_eta, m_phi, m_chi2;
  // GEMHits included in Muon
  vector<int> m_roll, m_chamber, m_layer; // hit info
  vector<float> m_resx, m_resy, m_pullx, m_pully;
  // Propagation only information
  vector<float> m_in_vfat;
  vector<int> m_in_roll, m_in_chamber, m_in_layer, m_in_goodEta; // propagation bound info
  vector<float> m_in_globalPhi, m_in_globalEta, m_in_nearGemPhi, m_in_nearGemEta; // global info
  vector<float> m_in_x, m_in_y, m_in_local_x, m_in_local_y, m_in_gemx, m_in_gemy, m_in_local_gemx, m_in_local_gemy, m_in_pullx, m_in_pully, m_in_resx, m_in_resy, m_in_trkextdx, m_in_gem_dx;
  vector<float> m_in_local_x_closetsos, m_in_local_y_closetsos, m_in_trkextdx_closetsos, m_in_trkextdx_inner;
  vector<float> m_in_local_x_inner, m_in_local_y_inner;
  vector<int> m_in_gemNStrips, m_in_gemFirstStrip, m_in_strip;
  vector<bool> m_in_matchingGem;

  vector<int> m_rec_roll, m_rec_chamber, m_rec_layer;

  vector<vector<float>> m_in_resx_tests;
  
  TTree *t_hit;
  TTree *t_run;
  TTree *t_muon;
  TTree *t_event;
};

STASliceTestAnalysis::STASliceTestAnalysis(const edm::ParameterSet& iConfig) :
  nEvents(0),
  nMuonTotal(0),
  nGEMFiducialMuon(0),
  nGEMTrackWithMuon(0),
  nGEMTotal(0)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  staTracks_ = consumes<reco::TrackCollection>(iConfig.getParameter<InputTag>("muons"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);
  lumiScalers_ = consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers"));
  amc13Event_ = consumes<GEMAMC13EventCollection>(iConfig.getParameter<edm::InputTag>("amc13Event"));
  amcData_ = consumes<GEMAMCdataCollection>(iConfig.getParameter<edm::InputTag>("amcData"));
  gebStatusCol_ = consumes<GEMGEBStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("gebStatusCol"));
  vfatStatusCol_ = consumes<GEMVfatStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("vfatStatusCol"));
 
  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nMuons", &b_nMuons, "nMuons/I");
  t_event->Branch("nMuonsWithGEMHit", &b_nMuonsWithGEMHit, "nMuonsWithGEMHit/I");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  t_event->Branch("nMuonsInMuonTree", &b_nMuonsInMuonTree, "nMuonsInMuonTree/I");
  t_event->Branch("nHitsInHitTree", &b_nHitsInHitTree, "nHitsInHitTree/I");
  t_event->Branch("event", &b_event, "event/l");
  t_event->Branch("run", &b_run, "run/I");
  t_event->Branch("lumi", &b_lumi, "lumi/I");
  t_event->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_event->Branch("latency", &b_latency, "latency/I");

  t_muon = fs->make<TTree>("Muon", "Muon");
  t_muon->Branch("event", &b_event, "event/l");
  t_muon->Branch("run", &b_run, "run/I");
  t_muon->Branch("lumi", &b_lumi, "lumi/I");
  t_muon->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_muon->Branch("latency", &b_latency, "latency/I");
  t_muon->Branch("chi2", &m_chi2, "chi2/F");
  t_muon->Branch("nhits", &m_nhits, "nhits/I");
  t_muon->Branch("nCSChits", &m_nCSChits, "nCSChits/I");
  t_muon->Branch("nGEMhits", &m_nGEMhits, "nGEMhits/I")->SetTitle("n GEM hits associated to muon");
  t_muon->Branch("nvalidhits", &m_nvalidhits, "nvalidhits/I")->SetTitle("n GEM hits associated to muon, and muon can propagate to eta partition of hit");
  t_muon->Branch("nbounds", &m_nbounds, "nbounds/I")->SetTitle("times muon is in GEM eta partition bounds");
  t_muon->Branch("pt", &m_pt, "pt/F");
  t_muon->Branch("eta", &m_eta, "eta/F");
  t_muon->Branch("phi", &m_phi, "phi/F");
  t_muon->Branch("charge", &m_charge, "charge/I");
  t_muon->Branch("quality", &m_quality, "quality/I")->SetTitle("muon quality :: 0:noid 1:looseID 2:tightID");
  t_muon->Branch("in_strip", &m_in_strip);
  t_muon->Branch("in_vfat", &m_in_vfat);
  t_muon->Branch("in_roll", &m_in_roll);
  t_muon->Branch("in_chamber", &m_in_chamber);
  t_muon->Branch("in_layer", &m_in_layer);
  t_muon->Branch("in_goodEta", &m_in_goodEta);
  t_muon->Branch("in_resx", &m_in_resx);
  t_muon->Branch("in_resx_tests", &m_in_resx_tests);
  t_muon->Branch("in_resy", &m_in_resy);
  t_muon->Branch("in_trkextdx", &m_in_trkextdx);
  t_muon->Branch("in_pullx", &m_in_pullx);
  t_muon->Branch("in_pully", &m_in_pully);
  t_muon->Branch("in_phi", &m_in_globalPhi);
  t_muon->Branch("in_eta", &m_in_globalEta);
  t_muon->Branch("in_nearGemPhi", &m_in_nearGemPhi);
  t_muon->Branch("in_nearGemEta", &m_in_nearGemEta);
  t_muon->Branch("in_x", &m_in_x);
  t_muon->Branch("in_y", &m_in_y);

  t_muon->Branch("in_local_x_inner", &m_in_local_x_inner);
  t_muon->Branch("in_local_y_inner", &m_in_local_y_inner);
  t_muon->Branch("in_local_x_closetsos", &m_in_local_x_closetsos);
  t_muon->Branch("in_local_y_closetsos", &m_in_local_y_closetsos);
  t_muon->Branch("in_trkextdx_closetsos", &m_in_trkextdx_closetsos);
  t_muon->Branch("in_trkextdx_inner", &m_in_trkextdx_inner);

  t_muon->Branch("in_local_x", &m_in_local_x);
  t_muon->Branch("in_local_y", &m_in_local_y);
  t_muon->Branch("in_gemx", &m_in_gemx);
  t_muon->Branch("in_gem_dx", &m_in_gem_dx);
  t_muon->Branch("in_gemy", &m_in_gemy);
  t_muon->Branch("in_local_gemx", &m_in_local_gemx);
  t_muon->Branch("in_local_gemy", &m_in_local_gemy);
  t_muon->Branch("in_gemFirstStrip", &m_in_gemFirstStrip);
  t_muon->Branch("in_gemNStrips", &m_in_gemNStrips);
  t_muon->Branch("in_matchingGem", &m_in_matchingGem);

  t_muon->Branch("rec_roll", &m_rec_roll);
  t_muon->Branch("rec_chamber", &m_rec_chamber);
  t_muon->Branch("rec_layer", &m_rec_layer);
  
  t_hit = fs->make<TTree>("Hit", "Hit");
  t_hit->Branch("event", &b_event, "event/l");
  t_hit->Branch("run", &b_run, "run/I");
  t_hit->Branch("lumi", &b_lumi, "lumi/I");
  t_hit->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_hit->Branch("latency", &b_latency, "latency/I");
  t_hit->Branch("bx", &b_bx, "bx/I");
  t_hit->Branch("firstStrip", &b_firstStrip, "firstStrip/I");
  t_hit->Branch("nStrips", &b_nStrips, "nStrips/I");
  t_hit->Branch("strips", &b_strips);
  t_hit->Branch("chamber", &b_chamber, "chamber/I");
  t_hit->Branch("layer", &b_layer, "layer/I");
  t_hit->Branch("etaPartition", &b_etaPartition, "etaPartition/I");
  t_hit->Branch("muonQuality", &b_muonQuality, "muonQuality/I")->SetTitle("muonQuality -1:none 0:noid 1:looseID 2:tightID");
  t_hit->Branch("x", &b_x, "x/F");
  t_hit->Branch("y", &b_y, "y/F");
  t_hit->Branch("z", &b_z, "z/F");
}

STASliceTestAnalysis::~STASliceTestAnalysis()
{
  std::cout << "::: GEM Slice Test Results :::" << std::endl;
  std::cout << ": From " << nEvents << " events" << std::endl;
  std::cout << std::endl;
  std::cout << " # GEMHits " << nGEMTotal << std::endl;
  std::cout << " # Muons   " << nMuonTotal << std::endl;
  std::cout << " # FidMu   " << nGEMFiducialMuon << std::endl;
  std::cout << " # GEMMu   " << nGEMTrackWithMuon << std::endl;
  std::cout << std::endl;
}

void
STASliceTestAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<LumiScalersCollection> lumiScalers;
  iEvent.getByToken(lumiScalers_, lumiScalers);

  edm::ESHandle<GEMGeometry> GEMGeometry_;
  iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);

  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<reco::TrackCollection> muons;
  iEvent.getByToken(staTracks_, muons);

  iEvent.getByToken(amc13Event_, amc13Event);
  iEvent.getByToken(amcData_, amcData);
  iEvent.getByToken(gebStatusCol_, gebStatusCol);
  iEvent.getByToken(vfatStatusCol_, vfatStatusCol);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttrackBuilder_);
  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  b_latency = -1;

  for (auto g : *amcData) {
    for (auto a = g.second.first; a != g.second.second; ++a) {
      if (b_latency != -1 && b_latency != a->param1())
 	std::cout << "CHANGING LATENCY - old: " << b_latency << " new: " << a->param1() << std::endl;
      b_latency = a->param1();
      if (b_latency == -1) std::cout << "-1 LATENCY VALUE - " << iEvent.id().event() << " " << iEvent.id().run() << std::endl;
    }
  }

  nEvents++;

  b_nMuons = 0;
  b_nMuonsWithGEMHit = 0;
  b_nGEMHits = 0;
  b_nMuonsInMuonTree = 0;
  b_nHitsInHitTree = 0;
  
  b_run = iEvent.run();
  b_lumi = iEvent.luminosityBlock();
  if (lumiScalers->size() > 0)
    b_instLumi = (lumiScalers->at(0)).instantLumi();
  else
    b_instLumi = -1;
  b_event = iEvent.id().event(); // = event number

  
  //  cout << "---- MUONS --- " << endl;
  for (auto & mu : *muons) {
    b_nMuons++;
    // only consider muons going in the right direction (toward the gem slice test)
    if (mu.eta() > 0) continue;
    
    m_nGEMhits = 0;
    m_nbounds = 0;
    m_in_strip.clear(); m_in_vfat.clear(); m_in_roll.clear(); m_in_chamber.clear(); m_in_layer.clear(); m_in_goodEta.clear();
    m_in_x.clear(); m_in_y.clear(); m_in_gemx.clear(); m_in_gem_dx.clear(); m_in_gemy.clear();
    m_in_local_x_closetsos.clear(); m_in_local_y_closetsos.clear();
    m_in_local_x_inner.clear(); m_in_local_y_inner.clear(); 
    m_in_trkextdx_closetsos.clear(); m_in_trkextdx_inner.clear();
    m_in_local_gemx.clear(); m_in_local_gemy.clear();
    m_in_local_x.clear(); m_in_local_y.clear();
    m_in_resx.clear(); m_in_trkextdx.clear(); m_in_resy.clear(); m_in_pullx.clear(); m_in_pully.clear();
    m_in_resx_tests.clear();
    m_in_matchingGem.clear();
    m_in_gemFirstStrip.clear(); m_in_gemNStrips.clear();
    m_in_nearGemPhi.clear(); m_in_nearGemEta.clear();
    m_in_globalPhi.clear(); m_in_globalEta.clear();

    m_rec_layer.clear(); m_rec_roll.clear(); m_rec_chamber.clear();
    m_quality = mu.qualityMask();
    
  //   if (mu.passed(reco::Muon::Selector::CutBasedIdTight)) m_quality = 2;
  //   else if (mu.passed(reco::Muon::Selector::CutBasedIdLoose)) m_quality = 1;
  //   else m_quality = 0;

    m_charge = mu.charge();

    const reco::Track* muonTrack = &mu;
  //   if (mu.globalTrack().isNonnull()) muonTrack = mu.globalTrack().get();
  //   else if (mu.outerTrack().isNonnull()) muonTrack = mu.outerTrack().get();
    if (muonTrack) {
      bool fidMu = false;
      reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
      for (auto ch : GEMGeometry_->etaPartitions()) {

  	TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),
  							      ch->surface());
  	if (!tsos.isValid()) continue;

  	GlobalPoint tsosGP = tsos.globalPosition();
	
  	const LocalPoint pos = ch->toLocal(tsosGP);
  	const LocalPoint pos2D(pos.x(), pos.y(), 0);
  	const BoundPlane& bps(ch->surface());

  	if (bps.bounds().inside(pos2D)) {
  	  fidMu = true;
  	  m_nbounds++;
	  
  	  auto gemid = ch->id();

  	  m_in_strip.push_back(ch->strip(pos2D));
  	  m_in_vfat.push_back(ch->strip(pos2D)/128.);
  	  m_in_roll.push_back(gemid.roll());
  	  m_in_chamber.push_back(gemid.chamber());
  	  m_in_layer.push_back(gemid.layer());

	  m_in_goodEta.push_back(checkEtaPartitionGood(ch));
	  
  	  float in_x, in_y;
  	  in_x = tsosGP.x(); in_y = tsosGP.y();
  	  m_in_x.push_back(in_x); m_in_y.push_back(in_y);
  	  m_in_globalPhi.push_back(tsosGP.phi());
  	  m_in_globalEta.push_back(tsosGP.eta());

  	  // propagate from the "closest" tsos
  	  float in_local_x_closetsos = -9999, in_local_y_closetsos = -9999, in_trkextdx_closetsos = -99;
  	  GlobalPoint gemtagp1 = ch->position();
  	  TrajectoryStateOnSurface tsos1 = ttTrack.stateOnSurface(gemtagp1);
  	  if (!tsos1.isValid()) goto out1;	  
  	  tsos1 = propagator->propagate(tsos1, ch->surface());
  	  if (!tsos1.isValid()) goto out1;
  	  in_local_x_closetsos = tsos1.localPosition().x(); in_local_y_closetsos = tsos1.localPosition().y();	  
  	  in_trkextdx_closetsos = std::sqrt(tsos1.localError().positionError().xx());
  	out1:
  	  m_in_local_x_closetsos.push_back(in_local_x_closetsos); m_in_local_y_closetsos.push_back(in_local_y_closetsos);
  	  m_in_trkextdx_closetsos.push_back(in_trkextdx_closetsos);

  	  // Propagate from the innermost state
  	  float in_local_x_inner = -9999, in_local_y_inner = -9999, in_trkextdx_inner = -99;
  	  TrajectoryStateOnSurface tsos2 = propagator->propagate(ttTrack.innermostMeasurementState(),
  								ch->surface());
  	  if (!tsos2.isValid()) goto out2;	  
  	  in_local_x_inner = tsos2.localPosition().x(); in_local_y_inner = tsos2.localPosition().y();
  	  in_trkextdx_inner = std::sqrt(tsos2.localError().positionError().xx());
  	out2:
  	  m_in_local_x_inner.push_back(in_local_x_inner); m_in_local_y_inner.push_back(in_local_y_inner);
  	  m_in_trkextdx_inner.push_back(in_trkextdx_inner);


  	  bool matchingGem = false;
  	  int nstrips = -1, firststrip = -1;
  	  float in_resx = -99, in_trkextdx, in_resy = -99, in_pullx = -99, in_pully = -99, in_gem_dx = -99;
  	  float gemEta = +99.0, gemPhi = +99.0, in_gemx=-999, in_gemy=-999, in_local_gemx=-999, in_local_gemy=-999, in_local_x=-999, in_local_y=-999;
  	  auto recHitsRange = gemRecHits->get(gemid);
  	  auto gemRecHit = recHitsRange.first;
  	  vector<float> in_resx_tests;

  	  in_local_x = pos.x();
  	  in_local_y = pos.y();

  	  LocalPoint loc;

  	  for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {
  	    auto gemGlob = ch->toGlobal(hit->localPosition());
  	    if (fabs(gemGlob.phi() - tsosGP.phi()) < gemPhi) {
  	      in_gemx = gemGlob.x();
  	      in_gemy = gemGlob.y();

  	      in_local_gemx = hit->localPosition().x();
  	      in_local_gemy = hit->localPosition().y();
	      
  	      matchingGem = true;
  	      firststrip = hit->firstClusterStrip();
  	      nstrips = hit->clusterSize();
  	      gemPhi = gemGlob.phi(); gemEta = gemGlob.eta();

  	      LocalPoint && tsos_localpos = tsos.localPosition();
  	      loc = hit->localPosition();
  	      LocalError && tsos_localerr = tsos.localError().positionError();
  	      LocalPoint && dethit_localpos = hit->localPosition();     
  	      LocalError && dethit_localerr = hit->localPositionError();
  	      in_gem_dx = dethit_localerr.xx();
  	      in_resx = (dethit_localpos.x() - tsos_localpos.x());
  	      in_trkextdx = std::sqrt(tsos.localError().positionError().xx());
  	      in_resy = (dethit_localpos.y() - tsos_localpos.y()); 
  	      in_pullx = (dethit_localpos.x() - tsos_localpos.x()) / 
  		std::sqrt(dethit_localerr.xx() + tsos_localerr.xx());
  	      in_pully = (dethit_localpos.y() - tsos_localpos.y()) / 
  		std::sqrt(dethit_localerr.yy() + tsos_localerr.yy());
  	    }
  	  }

  	  m_in_resx_tests.push_back(in_resx_tests); 
  	  m_in_gemx.push_back(in_gemx); m_in_gemy.push_back(in_gemy);
  	  m_in_gem_dx.push_back(in_gem_dx);
  	  m_in_local_x.push_back(in_local_x); m_in_local_y.push_back(in_local_y);
  	  m_in_local_gemx.push_back(in_local_gemx); m_in_local_gemy.push_back(in_local_gemy);
  	  m_in_resx.push_back(in_resx); m_in_trkextdx.push_back(in_trkextdx); m_in_resy.push_back(in_resy);
  	  m_in_pullx.push_back(in_pullx); m_in_pully.push_back(in_pully);
  	  m_in_gemNStrips.push_back(nstrips); m_in_gemFirstStrip.push_back(firststrip);
  	  m_in_matchingGem.push_back(matchingGem);
  	  m_in_nearGemPhi.push_back(gemPhi); m_in_nearGemEta.push_back(gemEta);
  	  // cout << "   --- eta: " << tsosGP.eta() << " phi: " << tsosGP.phi()  << " xyz: " << in_x << " " << in_y << " GEMxyz: " << in_gemx << " " << in_gemy << " qual:" << m_quality << " eta:" << mu.eta() << " phi:" << mu.phi() << " pt:" << mu.pt() << " ch:" << gemid.chamber() << " l:" << gemid.layer() << " p:" << gemid.roll() << endl;
  	}
      }

      m_nhits = 0; m_nCSChits = 0;
      for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
	m_nhits++;
	if ((*hit)->geographicalId().det() == DetId::Muon &&
  	    (*hit)->geographicalId().subdetId() == MuonSubdetId::CSC) m_nCSChits++;
	
  	if ((*hit)->geographicalId().det() == DetId::Muon &&
  	    (*hit)->geographicalId().subdetId() == MuonSubdetId::GEM) {
  	  GEMDetId gemid((*hit)->geographicalId());

  	  m_rec_chamber.push_back(gemid.chamber());
  	  m_rec_layer.push_back(gemid.layer());
  	  m_rec_roll.push_back(gemid.roll());
  	}
      }

      if (fidMu) ++nGEMFiducialMuon;
    }

    if (m_nGEMhits > 0 or m_nbounds > 0) {
      m_pt = mu.pt();
      m_eta = mu.eta();
      m_phi = mu.phi();
      m_pt = mu.normalizedChi2();

      t_muon->Fill();
      b_nMuonsInMuonTree++;
    }
  } // Muon Loop

  //  cout << "--- GEM HITS ---" << endl;
  for (auto & gem : *gemRecHits) {
    ++b_nGEMHits;
    auto detId = gem.gemId();

    b_strips.clear();
    b_chamber = detId.chamber();
    b_layer = detId.layer();
    b_etaPartition = detId.roll();
    b_firstStrip = gem.firstClusterStrip();
    b_bx = gem.BunchX();
    b_nStrips = gem.clusterSize();
    for (int nstrip = gem.firstClusterStrip(); nstrip < (gem.firstClusterStrip() + gem.clusterSize()); ++nstrip) {
      b_strips.push_back(nstrip);
    }

    auto roll = GEMGeometry_->etaPartition(detId);
    auto globalPosition = roll->toGlobal(gem.localPosition());
    b_x = globalPosition.x();
    b_y = globalPosition.y();
    b_z = globalPosition.z();
    
    t_hit->Fill();
    b_nHitsInHitTree++;
    //    cout << "   --- eta: " << globalPosition.eta() << " phi: " << globalPosition.phi()  << " xyz: " << b_x << " " << b_y << " " << b_z << " ch:" << detId.chamber() << " l:" << detId.layer() << " p:" << detId.roll() << endl;
  } // GEM Loop
  //  cout << endl << endl;

  nMuonTotal += b_nMuons;
  nGEMTotal += b_nGEMHits;

  t_event->Fill();
}

bool STASliceTestAnalysis::checkEtaPartitionGood(const GEMEtaPartition* part)
{
  GEMDetId rId = part->id();
  int amcBx = -1;

  for (GEMAMCdataCollection::DigiRangeIterator amcsIt = amcData->begin(); amcsIt != amcData->end(); ++amcsIt){
    auto amcs = (*amcsIt).second;
    for (auto amc = amcs.first; amc != amcs.second; ++amc) {
      amcBx = amc->bx();
    }
  }
  
  auto gebs = gebStatusCol->get(rId.chamberId()); 
  for (auto geb = gebs.first; geb != gebs.second; ++geb) {
    if (int(geb->getInFu()) != 0 ) return false;
  }
  
  auto vfats = vfatStatusCol->get(rId); 
  for (auto vfat = vfats.first; vfat != vfats.second; ++vfat) {
    if (vfat->bc() != amcBx ) return false;
    if (int(vfat->quality()) != 0 ) return false;
    if (int(vfat->flag()) != 0 ) return false;
  }
  return true;
}

void STASliceTestAnalysis::beginJob(){}
void STASliceTestAnalysis::endJob(){}

void STASliceTestAnalysis::beginRun(Run const& run, EventSetup const&){
}
void STASliceTestAnalysis::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(STASliceTestAnalysis);
