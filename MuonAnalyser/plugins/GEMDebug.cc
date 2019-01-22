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

class GEMDebug : public edm::EDAnalyzer {
public:
  explicit GEMDebug(const edm::ParameterSet&);
  ~GEMDebug();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalers_;
  // edm::EDGetTokenT<MuonDigiCollection<unsigned short,GEMAMCStatusDigi>> gemDigis_;
  // // edm::EDGetTokenT<MuonDigiCollection<GEMDetId,GEMVfatStatusDigi>> vfatStatus_;
  // edm::EDGetTokenT<MuonDigiCollection<GEMDetId,GEMGEBStatusDigi>> gebStatus_;

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
  int m_nhits, m_nvalidhits;
  int m_nbounds;
  int m_quality, m_charge;
  float m_pt, m_eta, m_phi;
  // GEMHits included in Muon
  vector<int> m_roll, m_chamber, m_layer; // hit info
  vector<float> m_resx, m_resy, m_pullx, m_pully;
  // Propagation only information
  vector<float> m_in_vfat;
  vector<int> m_in_roll, m_in_chamber, m_in_layer; // propagation bound info
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

GEMDebug::GEMDebug(const edm::ParameterSet& iConfig) :
  nEvents(0),
  nMuonTotal(0),
  nGEMFiducialMuon(0),
  nGEMTrackWithMuon(0),
  nGEMTotal(0)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);
  lumiScalers_ = consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers"));
}

GEMDebug::~GEMDebug()
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

#include <sstream>

void
GEMDebug::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<LumiScalersCollection> lumiScalers;
  iEvent.getByToken(lumiScalers_, lumiScalers);

  edm::ESHandle<GEMGeometry> GEMGeometry_;
  iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);

  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);

  Handle<View<reco::Muon> > muons;
  iEvent.getByToken(muons_, muons);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttrackBuilder_);
  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  bool output = false;
  std::stringstream ss; 
  
  ss << "\n\n ---------- \n\n" << "\n";
  ss << "nEvent: " << nEvents++ << " : " << iEvent.id().event() << "\n\n";
  ss << "GEM REC HITS" << "\n";
  for (auto & gem : *gemRecHits) {
    auto detId = gem.gemId();

    auto roll = GEMGeometry_->etaPartition(detId);
    auto globalPosition = roll->toGlobal(gem.localPosition());
    b_x = globalPosition.x();
    b_y = globalPosition.y();
    b_z = globalPosition.z();

    ss << "   --- eta: " << globalPosition.eta() << " phi: " << globalPosition.phi()  << " xyz: " << b_x << " " << b_y << " " << b_z << " ch:" << detId.chamber() << " l:" << detId.layer() << " p:" << detId.roll() << " s:" << gem.firstClusterStrip() << " [" << gem.clusterSize() << "]" << "\n";
  }

  ss << "\n\nMUONS" << "\n";
  for (auto & mu : *muons) {
    ss << "  -- Mu e:" << mu.eta() << " p:" << mu.phi() << " pt:" << mu.pt() << "\n";
    // only consider muons going in the right direction (toward the gem slice test)
    if (mu.eta() > 0) continue;

    if (mu.passed(reco::Muon::Selector::CutBasedIdTight)) m_quality = 2;
    else if (mu.passed(reco::Muon::Selector::CutBasedIdLoose)) m_quality = 1;
    else m_quality = 0;

    
    const reco::Track* muonTrack = 0;
    if (mu.globalTrack().isNonnull()) muonTrack = mu.globalTrack().get();
    else if (mu.outerTrack().isNonnull()) muonTrack = mu.outerTrack().get();
    if (muonTrack) {
      reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
      for (auto ch : GEMGeometry_->etaPartitions()) {

	TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),
							      ch->surface());
	if (!tsos.isValid()) continue;
	auto gemid = ch->id();
	if (gemid.chamber() == 1) continue;

	GlobalPoint tsosGP = tsos.globalPosition();
	const LocalPoint pos = ch->toLocal(tsosGP);
	const LocalPoint pos2D(pos.x(), pos.y(), 0);
	const BoundPlane& bps(ch->surface());

	if (bps.bounds().inside(pos2D)) {
	  output = true;
	  ss << "    ONgem: " << gemid.chamber() << " " << gemid.layer() << " " << gemid.roll() << " " << ch->strip(pos2D) << " - xyz: " << tsosGP.x() << " " << tsosGP.y() << " " << tsosGP.z() << "\n";
	}
      }

      for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
	if ((*hit)->geographicalId().det() == DetId::Muon &&
	    (*hit)->geographicalId().subdetId() == MuonSubdetId::GEM) {
	  output = true;
	  GEMDetId gemid((*hit)->geographicalId());

	  ss << "      found hit " << gemid.chamber() << " " << gemid.layer() << " " << gemid.roll() << "\n";
	}
      }
    }
  }

  ss << " ---------- " << "\n";
  if (output) { cout << ss.str() << endl; }
}

void GEMDebug::beginJob(){}
void GEMDebug::endJob(){}

void GEMDebug::beginRun(Run const& run, EventSetup const&){
}
void GEMDebug::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(GEMDebug);
