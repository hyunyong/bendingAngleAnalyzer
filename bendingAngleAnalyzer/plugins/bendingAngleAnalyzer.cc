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
#include "DataFormats/Math/interface/deltaPhi.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "TTree.h"
#include "TString.h"

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

    edm::Service<TFileService> fs;

    TTree* ttree_;
   
    double bMuEta, bMuPhi, bMuPt, bMuCharge;
    double bGEMPhi, bCSCPhi, bDelPhi;
};

bendingAngleAnalyzer::bendingAngleAnalyzer(const edm::ParameterSet& iConfig)
{
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtx"));
  muonToken_ = consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
  cscRecHitToken_ = consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("cscRecHits"));
  gemRecHitToken_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));


  ttree_ = fs->make<TTree>("bendingAngleAnalyzer", "bendingAngleAnalyzer");

  ttree_->Branch("muEta", &bMuEta);
  ttree_->Branch("muPhi", &bMuPhi);
  ttree_->Branch("muPt", &bMuPt);
  ttree_->Branch("muCharge", &bMuCharge);
  ttree_->Branch("GEMPhi", &bGEMPhi);
  ttree_->Branch("CSCPhi", &bCSCPhi);
  ttree_->Branch("DelPhi", &bDelPhi);
 
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
    resetBr();
    std::cout << "muon eta: " << muon->eta() << ", muon phi: " << muon->phi() << ", muon pT:" << muon->pt() << ", muon charge: " << muon->charge() <<std::endl;
    bMuEta = muon->eta(); bMuPhi = muon->phi(); bMuPt = muon->pt(); bMuCharge = muon->charge();
    for (const auto& gemHit:*(gemRecHits.product())){
      auto gemDetId = gemHit.gemId();
      if (!(gemDetId.station() == 1 and gemDetId.layer() ==2 and gemDetId.region()*muon->eta() > 0.0)) continue;
      const auto& gemDet = GEMGeometry_->idToDet(gemHit.geographicalId());
      double gemHitPhi = gemDet->toGlobal(gemHit.localPosition()).phi();
      std::cout << gemDetId << ", GEM recHit phi: " << gemHitPhi <<  std::endl;
      bGEMPhi = gemHitPhi;
        for (const auto& cscHit:*(cscRecHits.product())){
          auto cscDetId = cscHit.cscDetId();
          if (!(cscDetId.station() == 1 and cscDetId.ring() == 1 and cscDetId.layer() == 3 and cscDetId.chamber() == gemDetId.chamber() and (cscDetId.endcap() == 1 ? 1 : -1) == gemDetId.region())) continue;
          const auto& cscDet = CSCGeometry_->idToDet(cscHit.geographicalId());
          double cscHitPhi = cscDet->toGlobal(cscHit.localPosition()).phi();
          std::cout << cscDetId << ", CSC recHit phi:" << cscHitPhi << std::endl;
          bCSCPhi = cscHitPhi;
          bDelPhi = reco::deltaPhi(bGEMPhi, bCSCPhi);
          ttree_->Fill();
        }
    }
    //ttree_->Fill();
  }  
}

void bendingAngleAnalyzer::resetBr()
{
  bMuEta = -999.; bMuPhi = -999.; bMuPt = -999.; bMuCharge = -999.;
  bGEMPhi = -999.; bCSCPhi = -999.; bDelPhi = -999.;
}

DEFINE_FWK_MODULE(bendingAngleAnalyzer);
