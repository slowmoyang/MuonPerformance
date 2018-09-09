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

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMCStatusDigi.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

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

class SliceTestAnalysis : public edm::EDAnalyzer {
public:
  explicit SliceTestAnalysis(const edm::ParameterSet&);
  ~SliceTestAnalysis();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
  virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<CSCRecHit2DCollection> cscRecHits_;
  edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  edm::EDGetTokenT<MuonDigiCollection<unsigned short, GEMAMCStatusDigi>> gemDigisAMC_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalers_;
  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_; 
  
  TH2D* h_firstStrip[36][2];
  TH2D* h_allStrips[36][2];

  TH2D* h_activeLumi;
  TH1D* h_lumiStatus;
  
  bool vfatStatus[36][2][24];

  TH2D* h_globalPosOnGem;
  TH1D* h_instLumi[36][2][24];
  TH1D* h_pileup;
  TH1D* h_clusterSize, *h_totalStrips, *h_bxtotal;
  TH1D* h_inEta[36][2];
  TH1D* h_hitEta[36][2];
  TH1D* h_trkEta[36][2];

  TH1D* h_res_x, *h_res_y, *h_pull_x, *h_pull_y;

  TTree *t_setup;
  float b_stripLength[2][8];
  float b_stripPitch[2][8];
  float b_cscArea[2];
  
  TTree *t_run;
  
  TTree *t_event;
  int b_run, b_lumi, b_latency;
  int b_nGEMHits, b_nCSCHits;
  float b_instLumi;
  unsigned int b_timeLow, b_timeHigh;

  TTree *t_hit;
  int b_firstStrip, b_nStrips, b_chamber, b_layer, b_etaPartition;
  float b_x, b_y, b_z;
  
  TTree *t_csc;

};

SliceTestAnalysis::SliceTestAnalysis(const edm::ParameterSet& iConfig)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  cscRecHits_ = consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("cscRecHits"));
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  gemDigisAMC_ = consumes<MuonDigiCollection<unsigned short, GEMAMCStatusDigi>>(iConfig.getParameter<edm::InputTag>("gemDigisAMC"));
  vertexCollection_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  lumiScalers_ = consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);

  h_clusterSize=fs->make<TH1D>(Form("clusterSize"),"clusterSize",100,0,100);
  h_totalStrips=fs->make<TH1D>(Form("totalStrips"),"totalStrips",200,0,200);
  h_instLumi[0][0][0]=fs->make<TH1D>(Form("instLumi"),"instLumi",100,0,10);
  h_pileup=fs->make<TH1D>(Form("pileup"),"pileup",80,0,80);
  h_bxtotal=fs->make<TH1D>(Form("bx"),"bx",1000,0,1000);

  h_globalPosOnGem = fs->make<TH2D>(Form("onGEM"), "onGEM", 100, -100, 100, 100, -100, 100);

  t_run = fs->make<TTree>("Run", "Run");
  t_run->Branch("run", &b_run, "run/I");

  t_setup = fs->make<TTree>("Setup", "Setup");
  t_setup->Branch("stripLength", &b_stripLength, "stripLength[2][8]/F");
  t_setup->Branch("stripPitch", &b_stripPitch, "stripPitch[2][8]/F");
  t_setup->Branch("cscArea", &b_cscArea, "cscArea[2]/F");

  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  t_event->Branch("run", &b_run, "run/I");
  t_event->Branch("lumi", &b_lumi, "lumi/I");
  t_event->Branch("latency", &b_latency, "latency/I");
  t_event->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_event->Branch("timeLow", &b_timeLow, "timeLow/i");
  t_event->Branch("timeHigh", &b_timeHigh, "timeHigh/i");

  t_hit = fs->make<TTree>("Hit", "Hit");
  t_hit->Branch("firstStrip", &b_firstStrip, "firstStrip/I");
  t_hit->Branch("nStrips", &b_nStrips, "nStrips/I");
  t_hit->Branch("chamber", &b_chamber, "chamber/I");
  t_hit->Branch("layer", &b_layer, "layer/I");
  t_hit->Branch("etaPartition", &b_etaPartition, "etaPartition/I");
  t_hit->Branch("x", &b_x, "x/F");
  t_hit->Branch("y", &b_y, "y/F");
  t_hit->Branch("z", &b_z, "z/F");

  t_csc = fs->make<TTree>("CSC", "CSC");
  t_csc->Branch("chamber", &b_chamber, "chamber/I");
  t_csc->Branch("layer", &b_layer, "layer/I");
  t_csc->Branch("nStrips", &b_nStrips, "nStrips/I");

  for (int ichamber=26; ichamber<30;++ichamber) {
  // for (int ichamber=27; ichamber<=30;++ichamber) {
    for (int ilayer=0; ilayer<2;++ilayer) {
      h_firstStrip[ichamber][ilayer] = fs->make<TH2D>(Form("firstStrip ch %i lay %i",ichamber+1, ilayer+1),"firstStrip",384,1,385,8,0.5,8.5);
      h_firstStrip[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_firstStrip[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");
      
      h_allStrips[ichamber][ilayer] = fs->make<TH2D>(Form("allStrips ch %i lay %i",ichamber+1, ilayer+1),"allStrips",384,1,385,8,0.5,8.5);
      h_allStrips[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_allStrips[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");

      h_inEta[ichamber][ilayer] = fs->make<TH1D>(Form("inEta ch %i lay %i",ichamber+1, ilayer+1),"inEta",8,0.5,8.5);
      h_hitEta[ichamber][ilayer] = fs->make<TH1D>(Form("hitEta ch %i lay %i",ichamber+1, ilayer+1),"hitEta",8,0.5,8.5);
      h_trkEta[ichamber][ilayer] = fs->make<TH1D>(Form("trkEta ch %i lay %i",ichamber+1, ilayer+1),"trkEta",8,0.5,8.5);
      
      for (int vfat = 0; vfat < 24; vfat++) {
        h_instLumi[ichamber][ilayer][vfat] = fs->make<TH1D>(Form("instLumi ch %i lay %i", ichamber+1, ilayer+1),"instlumi", 50, 0, 2);
        h_instLumi[ichamber][ilayer][vfat]->GetXaxis()->SetTitle("instantaneous luminosity[cm^{-2}s^{-1}]");
        h_instLumi[ichamber][ilayer][vfat]->GetYaxis()->SetTitle("Event");
        vfatStatus[ichamber][ilayer][vfat] = false;
      }
    }
  }


}

SliceTestAnalysis::~SliceTestAnalysis()
{
  t_setup->Fill();
}

void
SliceTestAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_run = iEvent.run();
  b_lumi = iEvent.luminosityBlock();
  
  b_nGEMHits = 0;
  
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;
  
  //edm::ESHandle<CSCGeometry> hCSCGeom;
  //iSetup.get<MuonGeometryRecord>().get(hCSCGeom);
  //const CSCGeometry* CSCGeometry_ = &*hCSCGeom;
  
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  
  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<CSCRecHit2DCollection> cscRecHits;  
  iEvent.getByToken(cscRecHits_, cscRecHits);

  edm::Handle<MuonDigiCollection<unsigned short,GEMAMCStatusDigi>> gemDigisAMC;
  iEvent.getByToken(gemDigisAMC_, gemDigisAMC);

  edm::Handle<reco::VertexCollection> vertexCollection;
  iEvent.getByToken( vertexCollection_, vertexCollection );
  if(vertexCollection.isValid()) {
    vertexCollection->size();
  }

  edm::Handle<LumiScalersCollection> lumiScalers;
  iEvent.getByToken( lumiScalers_, lumiScalers );

  Handle<View<reco::Muon> > muons;
  iEvent.getByToken(muons_, muons);

  int totalStrips = 0;
  auto instLumi = (lumiScalers->at(0)).instantLumi()/10000;
  auto pileup = (lumiScalers->at(0)).pileup();
  h_instLumi[0][0][0]->Fill(instLumi);
  h_lumiStatus->Fill(instLumi);
  b_instLumi = instLumi;
  b_timeHigh = iEvent.time().unixTime();
  b_timeLow = iEvent.time().microsecondOffset();
  h_pileup->Fill(pileup);

  b_latency = -1;
  for (auto g : *gemDigisAMC) {
    for (auto a = g.second.first; a != g.second.second; a++) {
      b_latency = a->Param1();
    }
  }

  for (auto ch : GEMGeometry_->chambers()) {
    for(auto roll : ch->etaPartitions()) {
      GEMDetId rId = roll->id();

      auto recHitsRange = gemRecHits->get(rId); 
      auto gemRecHit = recHitsRange.first;
      
      for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {

	h_firstStrip[rId.chamber()-1][rId.layer()-1]->Fill(hit->firstClusterStrip(), rId.roll());
	h_clusterSize->Fill(hit->clusterSize());
	h_bxtotal->Fill(hit->BunchX());
	for (int nstrip = hit->firstClusterStrip(); nstrip < hit->firstClusterStrip()+hit->clusterSize(); ++nstrip) {
	  totalStrips++;
	  h_allStrips[rId.chamber()-1][rId.layer()-1]->Fill(nstrip, rId.roll());
	}

	b_firstStrip = hit->firstClusterStrip();
	b_nStrips = hit->clusterSize();
	b_chamber = rId.chamber();
	b_layer = rId.layer();
	b_etaPartition = rId.roll();
        int vfat = (8-b_etaPartition)+8*((b_firstStrip-1)/128);
        h_activeLumi->Fill(b_lumi, b_chamber+(b_layer-1)/2.+vfat/48.);
        vfatStatus[b_chamber-1][b_layer-1][vfat] = true;

	auto globalPosition = roll->toGlobal(hit->localPosition());
	b_x = globalPosition.x();
	b_y = globalPosition.y();
	b_z = globalPosition.z();

	t_hit->Fill();
	b_nGEMHits++;
      }
    }
  }

  for (int ichamber=26; ichamber<30;++ichamber) {
  // for (int ichamber=27; ichamber<=30;++ichamber) {
    for (int ilayer=0; ilayer<2;++ilayer) {
      for (int vfat = 0; vfat < 24; vfat++) {
        if (!vfatStatus[ichamber][ilayer][vfat]) continue;
        h_instLumi[ichamber][ilayer][vfat]->Fill(instLumi);
        vfatStatus[ichamber][ilayer][vfat] = false;
      }
    }
  }

  //for (auto ch : CSCGeometry_->chambers()) {
  //  for (auto layer : ch->layers()) {
  //    CSCDetId lId = layer->id();
  //    if (lId.station() != 1) continue;
  //    if (lId.endcap() != 2) continue;
  //    if (lId.chamber() < 27 or lId.chamber() > 30) continue;

  //    auto recHitsRange = cscRecHits->get(lId);
  //    auto cscRecHit = recHitsRange.first;
  //    
  //    for (auto hit = cscRecHit; hit != recHitsRange.second; ++hit) {
  //      b_nStrips = hit->nStrips();
  //      b_chamber = lId.chamber();
  //      b_layer = lId.layer();

  //      t_csc->Fill();
  //      b_nCSCHits++;
  //    }
  //  }
  //}

  h_totalStrips->Fill(totalStrips);

  t_event->Fill();
}

void SliceTestAnalysis::beginJob(){}
void SliceTestAnalysis::endJob(){}

void SliceTestAnalysis::beginRun(Run const& run, EventSetup const& iSetup){
  h_activeLumi = fs->make<TH2D>(Form("%i active lumi", run.run()),Form("Run number %i", run.run()),5000, 0, 5000, 400, 27, 31);
  
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;
  
  edm::ESHandle<CSCGeometry> hCSCGeom;
  iSetup.get<MuonGeometryRecord>().get(hCSCGeom);
  const CSCGeometry* CSCGeometry_ = &*hCSCGeom;
  
  for (auto ch : GEMGeometry_->chambers()) {
    for(auto roll : ch->etaPartitions()) {
      GEMDetId rId = roll->id();

      const TrapezoidalStripTopology* top_(dynamic_cast<const TrapezoidalStripTopology*>(&(roll->topology())));
      const float striplength(top_->stripLength());
      const float pitch(roll->pitch());

      b_stripLength[rId.chamber()%2][rId.roll()-1] = striplength;
      b_stripPitch[rId.chamber()%2][rId.roll()-1] = pitch;
    }
  }

  //for (auto ch : CSCGeometry_->chambers()) {
  //  for (auto layer : ch->layers()) {
  //    CSCDetId lId = layer->id();
  //    if (lId.station() != 1) continue;
  //    if (lId.endcap() != 2) continue;
  //    if (lId.chamber() < 27 or lId.chamber() > 30) continue;

  //    const OffsetRadialStripTopology* top_(dynamic_cast<const OffsetRadialStripTopology*>(&(layer->topology())));
  //    const float striplength(top_->stripLength());
  //    const float pitch(top_->pitch());
  //    b_cscArea[lId.chamber()%2] = striplength * pitch * top_->nstrips();
  //  }
  //}
  
  b_run = run.run();
  t_run->Fill();
}
void SliceTestAnalysis::endRun(Run const&, EventSetup const&){}

void SliceTestAnalysis::beginLuminosityBlock(LuminosityBlock const& lumiBlock, EventSetup const& iSetup){
  h_lumiStatus = fs->make<TH1D>(Form("%i %i status", lumiBlock.run(), lumiBlock.luminosityBlock()), "", 1000, 0, 2);
}
void SliceTestAnalysis::endLuminosityBlock(LuminosityBlock const& lumiBlock, EventSetup const& iSetup){}

//define this as a plug-in
DEFINE_FWK_MODULE(SliceTestAnalysis);
