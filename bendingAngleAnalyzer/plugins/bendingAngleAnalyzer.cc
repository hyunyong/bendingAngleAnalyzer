#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"


class bendingAngleAnalyzer : public edm::EDAnalyzer {
  public:
    explicit bendingAngleAnalyzer(const edm::ParameterSet&);
    ~bendingAngleAnalyzer(){};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

    void resetBr();

    edm::ESHandle<CSCGeometry> CSCGeometry_;
    edm::ESHandle<GEMGeometry> GEMGeometry_;

    edm::EDGetTokenT<edm::View<reco::Muon>> muonToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<GEMRecHitCollection> gemRecHitToken_;
    edm::EDGetTokenT<CSCRecHit2DCollection> cscRecHitToken_;

};

bendingAngleAnalyzer::bendingAngleAnalyzer(const edm::ParameterSet& iConfig)
{
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtx"));
  muonToken_ = consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
  cscRecHitToken_ = consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("cscRecHits"));
  gemRecHitToken_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
}

void bendingAngleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("vtx", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("gemRecHits", edm::InputTag("gemRecHits"));
  desc.add<edm::InputTag>("cscRecHits", edm::InputTag("csc2DRecHits"));
  descriptions.add("bendingAngleAnalyzer", desc);
}

void bendingAngleAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);
  iSetup.get<MuonGeometryRecord>().get(CSCGeometry_);
}

void bendingAngleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::VertexCollection> vertexCollection;
  if (!iEvent.getByToken( vtxToken_, vertexCollection )) return;
 
  edm::Handle<edm::View<reco::Muon>> muons;
  if (!iEvent.getByToken(muonToken_, muons)) return;

  edm::Handle<GEMRecHitCollection> gemRecHits;
  if (!iEvent.getByToken(gemRecHitToken_, gemRecHits)) return;

  edm::Handle<CSCRecHit2DCollection> cscRecHits;
  if (!iEvent.getByToken(cscRecHitToken_, cscRecHits)) return;

  for (auto muon = muons->begin(); muon != muons->end(); ++muon){
    if (!muon->isGEMMuon()) continue;
    std::cout << "muon eta: " << muon->eta() << ", muon phi: " << muon->phi() << ", muon pT:" << muon->pt() << ", muon charge: " << muon->charge() <<std::endl;
    for (const auto& gemHit:*(gemRecHits.product())){
      auto gemDetId = gemHit.gemId();
      if (!(gemDetId.station() == 1 and gemDetId.layer() ==2 and gemDetId.region()*muon->eta() > 0.0)) continue;
      const auto& gemDet = GEMGeometry_->idToDet(gemHit.geographicalId());
      double gemHitPhi = gemDet->toGlobal(gemHit.localPosition()).phi();
      std::cout << gemDetId << ", GEM recHit phi: " << gemHitPhi <<  std::endl;
        for (const auto& cscHit:*(cscRecHits.product())){
          auto cscDetId = cscHit.cscDetId();
          if (!(cscDetId.station() == 1 and cscDetId.ring() == 1 and cscDetId.layer() == 3 and cscDetId.chamber() == gemDetId.chamber() and (cscDetId.endcap() == 1 ? 1 : -1) == gemDetId.region())) continue;
          const auto& cscDet = CSCGeometry_->idToDet(cscHit.geographicalId());
          double cscHitPhi = cscDet->toGlobal(cscHit.localPosition()).phi();
          std::cout << cscDetId << ", CSC recHit phi:" << cscHitPhi << std::endl;
        }
    }
  }  
}

void bendingAngleAnalyzer::resetBr()
{
}

DEFINE_FWK_MODULE(bendingAngleAnalyzer);
