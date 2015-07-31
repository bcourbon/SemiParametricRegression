//
// Class: GedPhotonRegTreeMaker.cc
//
// Info: Processes a track into histograms of delta-phis and such
//
// Author: L. Gray (FNAL)
//


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"

#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"

#include "GedPhotonRegTreeMaker/GedPhotonRegTreeMaker/plugins/GedPhotonRegTreeMaker.h"

#include "TTree.h"
#include "TVector2.h"

#include <map>


namespace MK = reco::MustacheKernel;

GedPhotonRegTreeMaker::GedPhotonRegTreeMaker(const PSet& p) 
{
    _calib.reset(new PFEnergyCalibration());
    _N_ECALClusters = 1;
    _N_PSClusters   = 1;
    _tree = _fs->make<TTree>("RegressionTree","Regression tree for GedPhotons");
    _tree->Branch("eventNumber"             , &_eventNumber            , "eventNumber/I");
    _tree->Branch("N_ECALClusters"          , &_N_ECALClusters         , "N_ECALClusters/I");
    _tree->Branch("nVtx"                    , &_nVtx                   , "nVtx/I");
    _tree->Branch("rho"                     , &_rho                    , "rho/F");
    _tree->Branch("hoe"                     , &_hoe                    , "hoe/F");
    _tree->Branch("N_PSClusters"            , &_N_PSClusters           , "N_PSClusters/I");
    _tree->Branch("scIndex"                 , &_scIndex                , "scIndex/I");
    _tree->Branch("scRawEnergy"             , &_scRawEnergy            , "scRawEnergy/F");
    _tree->Branch("scCalibratedEnergy"      , &_scCalibratedEnergy     , "scCalibratedEnergy/F");
    _tree->Branch("scPreshowerEnergy"       , &_scPreshowerEnergy      , "scPreshowerEnergy/F");
    _tree->Branch("scPreshowerEnergyPlane1" , &_scPreshowerEnergyPlane1, "scPreshowerEnergyPlane1/F");
    _tree->Branch("scPreshowerEnergyPlane2" , &_scPreshowerEnergyPlane2, "scPreshowerEnergyPlane2/F");
    _tree->Branch("scIsEB"                  , &_scIsEB                 , "scIsEB/I");
    _tree->Branch("scEta"                   , &_scEta                  , "scEta/F");
    _tree->Branch("scPhi"                   , &_scPhi                  , "scPhi/F");
    _tree->Branch("scR"                     , &_scR                    , "scR/F");
    _tree->Branch("scPhiWidth"              , &_scPhiWidth             , "scPhiWidth/F");
    _tree->Branch("scEtaWidth"              , &_scEtaWidth             , "scEtaWidth/F");
    _tree->Branch("scSeedRawEnergy"         , &_scSeedRawEnergy        , "scSeedRawEnergy/F");
    _tree->Branch("scSeedCalibratedEnergy"  , &_scSeedCalibratedEnergy , "scSeedCalibratedEnergy/F");
    _tree->Branch("scSeedEta"               , &_scSeedEta              , "scSeedEta/F");
    _tree->Branch("scSeedPhi"               , &_scSeedPhi              , "scSeedPhi/F");
    _tree->Branch("scSeedSize"              , &_scSeedSize             , "scSeedSize/I");
    _tree->Branch("scSeedR9"                , &_scSeedR9               , "scSeedR9/F");
    _tree->Branch("scSeedEmax"              , &_scSeedEmax             , "scSeedEmax/F");
    _tree->Branch("scSeedE2nd"              , &_scSeedE2nd             , "scSeedE2nd/F");
    _tree->Branch("scSeedLeftRightAsym"     , &_scSeedLeftRightAsym    , "scSeedLeftRightAsym/F");
    _tree->Branch("scSeedTopBottomAsym"     , &_scSeedTopBottomAsym    , "scSeedTopBottomAsym/F");
    _tree->Branch("scSeedE2x5max"           , &_scSeedE2x5max          , "scSeedE2x5max/F");
    _tree->Branch("scSeed2x5LeftRightAsym"  , &_scSeed2x5LeftRightAsym , "scSeed2x5LeftRightAsym/F");
    _tree->Branch("scSeed2x5TopBottomAsym"  , &_scSeed2x5TopBottomAsym , "scSeed2x5TopBottomAsym/F");
    _tree->Branch("scSeedE3x3"  , &_scSeedE3x3 , "scSeedE3x3/F");
    _tree->Branch("scSeedE5x5"  , &_scSeedE5x5 , "scSeedE5x5/F");
    _tree->Branch("scSeedETop"  , &_scSeedETop , "scSeedETop/F");
    _tree->Branch("scSeedEBottom"  , &_scSeedEBottom , "scSeedEBottom/F");
    _tree->Branch("scSeedELeft"  , &_scSeedELeft , "scSeedELeft/F");
    _tree->Branch("scSeedERight"  , &_scSeedERight , "scSeedERight/F");
    _tree->Branch("scSeedE2x5Top"  , &_scSeedE2x5Top , "scSeedE2x5Top/F");
    _tree->Branch("scSeedE2x5Bottom"  , &_scSeedE2x5Bottom , "scSeedE2x5Bottom/F");
    _tree->Branch("scSeedE2x5Left"  , &_scSeedE2x5Left , "scSeedE2x5Left/F");
    _tree->Branch("scSeedE2x5Right"  , &_scSeedE2x5Right , "scSeedE2x5Right/F");
    _tree->Branch("scSeedSigmaIetaIeta"     , &_scSeedSigmaIetaIeta    , "scSeedSigmaIetaIeta/F");
//    _tree->Branch("scSeedSigmaIetaIeta_full5x5"     , &_scSeedSigmaIetaIeta_full5x5    , "scSeedSigmaIetaIeta_full5x5/F");
    _tree->Branch("scSeedSigmaIetaIphi"     , &_scSeedSigmaIetaIphi    , "scSeedSigmaIetaIphi/F");
    _tree->Branch("scSeedSigmaIphiIphi"     , &_scSeedSigmaIphiIphi    , "scSeedSigmaIphiIphi/F");
    _tree->Branch("scSeedCryEta"            , &_scSeedCryEta           , "scSeedCryEta/F");
    _tree->Branch("scSeedCryPhi"            , &_scSeedCryPhi           , "scSeedCryPhi/F");
    _tree->Branch("scSeedCryIeta"           , &_scSeedCryIeta          , "scSeedCryIeta/F");
    _tree->Branch("scSeedCryIphi"           , &_scSeedCryIphi          , "scSeedCryIphi/F");
    _tree->Branch("scSeedCryX"              , &_scSeedCryX             , "scSeedCryX/F");
    _tree->Branch("scSeedCryY"              , &_scSeedCryY             , "scSeedCryY/F");
    _tree->Branch("scSeedCryIx"             , &_scSeedCryIx            , "scSeedCryIx/F");
    _tree->Branch("scSeedCryIy"             , &_scSeedCryIy            , "scSeedCryIy/F");
    _tree->Branch("scSeedCryIeta_glob"           , &_scSeedCryIeta_glob          , "scSeedCryIeta_glob/F");
    _tree->Branch("scSeedCryIphi_glob"           , &_scSeedCryIphi_glob          , "scSeedCryIphi_glob/F");
    _tree->Branch("scSeedCryIx_glob"             , &_scSeedCryIx_glob            , "scSeedCryIx_glob/F");
    _tree->Branch("scSeedCryIy_glob"             , &_scSeedCryIy_glob            , "scSeedCryIy_glob/F");

    // ecal cluster information
    _clusterRawEnergy     .reset(new float[1], array_deleter<float>());
    _clusterCalibEnergy   .reset(new float[1], array_deleter<float>());
    _clusterEta           .reset(new float[1], array_deleter<float>());
    _clusterPhi           .reset(new float[1], array_deleter<float>());
    _clusterDEtaToSeed    .reset(new float[1], array_deleter<float>());  
    _clusterDPhiToSeed    .reset(new float[1], array_deleter<float>());
    _clusterDPhiToCentroid.reset(new float[1], array_deleter<float>());
    _clusterDEtaToCentroid.reset(new float[1], array_deleter<float>());
    _clusterInMustache    .reset(new int[1],   array_deleter<int>());
    _clusterInDynDPhi     .reset(new int[1],   array_deleter<int>());
    _clusterLeakage       .reset(new int[1],   array_deleter<int>());
    _tree->Branch("clusterRawEnergy"      , _clusterRawEnergy.get()     , "clusterRawEnergy[N_ECALClusters]/F");
    _tree->Branch("clusterCalibEnergy"    , _clusterCalibEnergy.get()   , "clusterCalibEnergy[N_ECALClusters]/F");
    _tree->Branch("clusterEta"            , _clusterEta.get()           , "clusterEta[N_ECALClusters]/F");
    _tree->Branch("clusterPhi"            , _clusterPhi.get()           , "clusterPhi[N_ECALClusters]/F");
    _tree->Branch("clusterDPhiToSeed"     , _clusterDPhiToSeed.get()    , "clusterDPhiToSeed[N_ECALClusters]/F");
    _tree->Branch("clusterDEtaToSeed"     , _clusterDEtaToSeed.get()    , "clusterDEtaToSeed[N_ECALClusters]/F");  
    _tree->Branch("clusterDPhiToCentroid" , _clusterDPhiToCentroid.get(), "clusterDPhiToCentroid[N_ECALClusters]/F");
    _tree->Branch("clusterDEtaToCentroid" , _clusterDEtaToCentroid.get(), "clusterDEtaToCentroid[N_ECALClusters]/F");
    _tree->Branch("clusterInMustache"     , _clusterInMustache.get()    , "clusterInMustache[N_ECALClusters]/I");
    _tree->Branch("clusterInDynDPhi"      , _clusterInDynDPhi.get()     , "clusterInDynDPhi[N_ECALClusters]/I");
    _tree->Branch("clusterLeakage"        , _clusterLeakage.get()       , "clusterLeakage[N_ECALClusters]/I");

    _tree->Branch("clusterMaxDR"          , &_clusterMaxDR              , "clusterMaxDR/F");
    _tree->Branch("clusterMaxDRDPhi"      , &_clusterMaxDRDPhi          , "clusterMaxDRDPhi/F");
    _tree->Branch("clusterMaxDRDEta"      , &_clusterMaxDRDEta          , "clusterMaxDRDEta/F");
    _tree->Branch("clusterMaxDRRawEnergy" , &_clusterMaxDRRawEnergy     , "clusterMaxDRRawEnergy/F");

    _tree->Branch("clustersMeanRawEnergy" , &_clustersMeanRawEnergy     , "clustersMeanRawEnergy/F");
    _tree->Branch("clustersRMSRawEnergy"  , &_clustersRMSRawEnergy      , "clustersRMSRawEnergy/F");
    _tree->Branch("clustersMeanDRToSeed"  , &_clustersMeanDRToSeed      , "clustersMeanDRToSeed/F");
    _tree->Branch("clustersMeanDEtaToSeed", &_clustersMeanDEtaToSeed    , "clustersMeanDEtaToSeed/F");
    _tree->Branch("clustersMeanDPhiToSeed", &_clustersMeanDPhiToSeed    , "clustersMeanDPhiToSeed/F");

    // preshower information
    _psClusterRawEnergy.reset(new float[1], array_deleter<float>());
    _psClusterEta      .reset(new float[1], array_deleter<float>());
    _psClusterPhi      .reset(new float[1], array_deleter<float>());
    _tree->Branch("psClusterRawEnergy"    , _psClusterRawEnergy.get()   , "psClusterRawEnergy[N_PSClusters]/F");
    _tree->Branch("psClusterEta"          , _psClusterEta.get()         , "psClusterEta[N_PSClusters]/F");
    _tree->Branch("psClusterPhi"          , _psClusterPhi.get()         , "psClusterPhi[N_PSClusters]/F");


    _tree->Branch("isMatched"      , &_isMatched      , "isMatched/I");
    _tree->Branch("genPt"          , &_genPt          , "genPt/F");
    _tree->Branch("genEta"         , &_genEta         , "genEta/F");
    _tree->Branch("genPhi"         , &_genPhi         , "genPhi/F");
    _tree->Branch("genEnergy"      , &_genEnergy      , "genEnergy/F");
    _tree->Branch("genDEoE"        , &_genDEoE        , "genDEoE/F");
    _tree->Branch("genDRToCentroid", &_genDRToCentroid, "genDRToCentroid/F");
    _tree->Branch("genDRToSeed"    , &_genDRToSeed    , "genDRToSeed/F");

    _clusterDPhiToGen.reset(new float[1],array_deleter<float>());
    _clusterDEtaToGen.reset(new float[1],array_deleter<float>());
    _tree->Branch("clusterDPhiToGen", _clusterDPhiToGen.get(), "clusterDPhiToGen[N_ECALClusters]/F");
    _tree->Branch("clusterDEtaToGen", _clusterDEtaToGen.get(), "clusterDPhiToGen[N_ECALClusters]/F");


    _phInput                  = p.getParameter<edm::InputTag>("photonSrc");
    _ebReducedRecHitCollection = p.getParameter<edm::InputTag>("ebReducedRecHitCollection");
    _eeReducedRecHitCollection = p.getParameter<edm::InputTag>("eeReducedRecHitCollection");
    _PrimaryVertex             = p.getParameter<edm::InputTag>("VtxLabel"); 
    _geninput                  = p.getParameter<edm::InputTag>("genSrc");    
    _matchMinERatio            = p.getParameter<double>("matchMinERatio");
    _matchMaxDR                = p.getParameter<double>("matchMaxDR");
    _matchingType              = p.getParameter<int>("matchingType");// 0=DR; 1=DE/E, 2=DR+DE/E

}


void GedPhotonRegTreeMaker::analyze(const edm::Event& e, const edm::EventSetup& es) 
{

    _lazyTool = new EcalClusterLazyTools(e, es, _ebReducedRecHitCollection, _eeReducedRecHitCollection); 

    edm::Handle<reco::PhotonCollection> phs;
    edm::Handle<reco::VertexCollection> PV;
    e.getByLabel(_phInput, phs);  
    e.getByLabel(_PrimaryVertex, PV);
    _eventNumber = e.id().event();
    _nVtx = PV->size();
    _scIndex = 0;
    _matchedParticles.clear();
    _matchedPhotons.clear();

    //Get Rho
    edm::Handle<double> hRhoKt6PFJets;
    e.getByLabel(edm::InputTag("fixedGridRhoFastjetAll",""), hRhoKt6PFJets);
    _rho = (*hRhoKt6PFJets);

    //edm::InputTag theRhoTag = InputTag("fixedGridRhoFastjetAll","rho","RECO");
  //InputTag theRhoTag = InputTag("fixedGridRhoFastjetAll","rho","RECO");
    //edm::Handle<double> rhoHandle;
    //e.getByLabel("fixedGridRhoFastjetAll", rhoHandle);
    //double rho = *rhoHandle;
    //_rho = rho ;  

    edm::Handle<reco::GenParticleCollection> gens;
    e.getByLabel(_geninput,gens);
    if(phs.isValid() && gens.isValid())
    {
        matchPhotonToGen(*phs, *gens);
    }
    else 
    {
        throw cms::Exception("GedPhotonRegTreeMaker")
            << "Product ID was invalid!"
            << std::endl;
    }
    std::vector< std::pair<reco::GenParticle, int> >::const_iterator it = _matchedParticles.begin();
    std::vector< std::pair<reco::GenParticle, int> >::const_iterator itE = _matchedParticles.end();
    for(;it!=itE;++it)
    {

        if(it->second!=-1)
        {
            if(_matchedPhotons[it->second].superCluster().isNull())
            {
                throw cms::Exception("GedPhotonRegTreeMaker")
                    << "Invalid photon supercluster"
                    << std::endl;
            }
            processSuperClusterFillTree(e, es, *(_matchedPhotons[it->second].superCluster()), it->first);
	    //Fill H/E
	    _hoe=_matchedPhotons[it->second].hadronicOverEm();
        }
        else
        {
            processNoMatchFillTree(e, es, it->first);
        }
    }

    delete _lazyTool;
}

void GedPhotonRegTreeMaker::processSuperClusterFillTree(const edm::Event& e, const edm::EventSetup& es, const reco::SuperCluster& sc, const reco::GenParticle& gen) 
{
    const int N_ECAL = sc.clustersEnd() - sc.clustersBegin();
    const int N_PS   = sc.preshowerClustersEnd() - sc.preshowerClustersBegin();

    _N_ECALClusters = std::max(0,N_ECAL - 1); // minus 1 because of seed
    _N_PSClusters   = N_PS;


    // generated information
    _isMatched       = 1;
    _genEnergy       = gen.energy();
    _genPt           = gen.pt();
    _genEta          = gen.eta();
    _genPhi          = gen.phi();
    _genDEoE         = (gen.energy()-sc.energy())/gen.energy();
    _genDRToCentroid = reco::deltaR(gen,sc);
    _genDRToSeed     = reco::deltaR(gen,*(sc.seed()));


    // supercluster information
    _scIndex                += 1;
    _scRawEnergy             = sc.rawEnergy();
    _scCalibratedEnergy      = sc.energy();
    _scPreshowerEnergy       = sc.preshowerEnergy();
    _scPreshowerEnergyPlane1 = sc.preshowerEnergyPlane1();
    _scPreshowerEnergyPlane2 = sc.preshowerEnergyPlane2();
    _scIsEB                  = (sc.seed()->hitsAndFractions().at(0).first.subdetId()==EcalBarrel);
    _scEta                   = sc.position().Eta();
    _scPhi                   = sc.position().Phi();
    _scR                     = sc.position().R();
    _scPhiWidth              = sc.phiWidth();
    _scEtaWidth              = sc.etaWidth();

    // sc seed information
    edm::Ptr<reco::CaloCluster> theseed = sc.seed();
    _scSeedRawEnergy        = theseed->energy();
    _scSeedCalibratedEnergy = theseed->energy();
    _scSeedEta              = theseed->eta();
    _scSeedPhi              = theseed->phi();
    _scSeedSize             = theseed->hitsAndFractions().size();
    _scSeedE3x3              = _lazyTool->e3x3(*(theseed));
    _scSeedE5x5              = _lazyTool->e5x5(*(theseed));
    _scSeedR9               = _scSeedE3x3/_scRawEnergy;
    _scSeedEmax             = _lazyTool->eMax(*(theseed));
    _scSeedE2nd             = _lazyTool->e2nd(*(theseed));
    _scSeedE2x5max          = _lazyTool->e2x5Max(*(theseed));
    _scSeedELeft             = _lazyTool->eLeft(*(theseed));
    _scSeedERight            = _lazyTool->eRight(*(theseed));
    _scSeedETop              = _lazyTool->eTop(*(theseed));
    _scSeedEBottom           = _lazyTool->eBottom(*(theseed));
    _scSeedE2x5Left          = _lazyTool->e2x5Left(*(theseed));
    _scSeedE2x5Right         = _lazyTool->e2x5Right(*(theseed));
    _scSeedE2x5Top           = _lazyTool->e2x5Top(*(theseed));
    _scSeedE2x5Bottom        = _lazyTool->e2x5Bottom(*(theseed));
    _scSeedLeftRightAsym    = (_scSeedELeft+_scSeedERight!=0. ? (_scSeedELeft-_scSeedERight)/(_scSeedELeft+_scSeedERight) : 0.);
    _scSeedTopBottomAsym    = (_scSeedETop+_scSeedEBottom!=0. ? (_scSeedETop-_scSeedEBottom)/(_scSeedETop+_scSeedEBottom) : 0.);
    _scSeed2x5LeftRightAsym    = (_scSeedE2x5Left+_scSeedE2x5Right!=0. ? (_scSeedE2x5Left-_scSeedE2x5Right)/(_scSeedE2x5Left+_scSeedE2x5Right) : 0.);
    _scSeed2x5TopBottomAsym    = (_scSeedE2x5Top+_scSeedE2x5Bottom!=0. ? (_scSeedE2x5Top-_scSeedE2x5Bottom)/(_scSeedE2x5Top+_scSeedE2x5Bottom) : 0.);
    std::vector<float> vCov = _lazyTool->localCovariances(*(theseed));
    double see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
    double spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    double sep = 0.;
    if (see*spp > 0)
        sep = vCov[1] / (see * spp);
    else if (vCov[1] > 0)
        sep = 1.0;
    else
        sep = -1.0;
    _scSeedSigmaIetaIeta    = see;
    _scSeedSigmaIphiIphi    = spp;
    _scSeedSigmaIetaIphi    = sep;


    EcalClusterLocal ecalLocal;

    if(theseed->hitsAndFractions().at(0).first.subdetId()==EcalBarrel)
    {
        float cryPhi, cryEta, thetatilt, phitilt;
        int ieta, iphi;
        ecalLocal.localCoordsEB(*(theseed), es, cryEta, cryPhi, ieta, iphi, thetatilt, phitilt);
        _scSeedCryEta  = cryEta;
        _scSeedCryPhi  = cryPhi;
        _scSeedCryIeta = ieta;
        _scSeedCryIphi = iphi;
        _scSeedCryX  = 0;
        _scSeedCryY  = 0;
        _scSeedCryIx = 0;
        _scSeedCryIy = 0;
    }
    else
    {
        float cryX, cryY, thetatilt, phitilt;
        int ix, iy;
        ecalLocal.localCoordsEE(*(theseed), es, cryX, cryY, ix, iy, thetatilt, phitilt);
        _scSeedCryX  = cryX;
        _scSeedCryY  = cryY;
        _scSeedCryIx = ix;
        _scSeedCryIy = iy;
        _scSeedCryEta  = 0;
        _scSeedCryPhi  = 0;
        _scSeedCryIeta = 0;
        _scSeedCryIphi = 0;
    }

    //Take global coordinates of the seed

    DetId globseed = theseed->seed();
    if (theseed->hitsAndFractions().at(0).first.subdetId()==EcalBarrel) {
    EBDetId ebseed(globseed);
    _scSeedCryIeta_glob = ebseed.ieta();
    _scSeedCryIphi_glob = ebseed.iphi();
    }
    else {
    EEDetId eeseed(globseed);
    _scSeedCryIx_glob = eeseed.ix();
    _scSeedCryIy_glob = eeseed.iy();
    }


    // loop over all clusters that aren't the seed
    setTreeArraysForSize(_N_ECALClusters, _N_PSClusters);
    _clusterMaxDR     = 999.;
    _clusterMaxDRDPhi = 999.;
    _clusterMaxDRDEta = 999.;
    _clustersMeanDRToSeed = 999.;
    _clustersMeanDEtaToSeed = 999.;
    _clustersMeanDPhiToSeed = 999.;
    _clusterMaxDRRawEnergy      = 0.;
    _clustersMeanRawEnergy      = 0.;
    _clustersRMSRawEnergy       = 0.;
    float clustersMeanSquareRawEnergy = 0.;
    float subclustersRawEnergy = 0.;
    auto clusend = sc.clustersEnd();
    size_t iclus = 0;
    float maxDR = 0;
    edm::Ptr<reco::CaloCluster> pclus;
    for( auto clus = sc.clustersBegin(); clus != clusend; ++clus ) 
    {
        pclus = *clus;

        _clustersMeanRawEnergy += pclus->energy();
        clustersMeanSquareRawEnergy += (pclus->energy()*pclus->energy());

        if( theseed == pclus ) continue;
        _clusterRawEnergy.get()[iclus]      = pclus->energy();
        _clusterCalibEnergy.get()[iclus]    = pclus->energy();
        _clusterEta.get()[iclus]            = pclus->eta();
        _clusterPhi.get()[iclus]            = pclus->phi();
        _clusterDPhiToSeed.get()[iclus]     = TVector2::Phi_mpi_pi(pclus->phi() - theseed->phi());
        _clusterDEtaToSeed.get()[iclus]     = pclus->eta() - theseed->eta();
        _clusterDPhiToCentroid.get()[iclus] = TVector2::Phi_mpi_pi(pclus->phi() - sc.phi());
        _clusterDEtaToCentroid.get()[iclus] = pclus->eta() - sc.eta();
        // find cluster with max dR
        if(reco::deltaR(*pclus, *theseed) > maxDR)
        {
            maxDR = reco::deltaR(*pclus, *theseed);
            _clusterMaxDR = maxDR;
            _clusterMaxDRDPhi = _clusterDPhiToSeed.get()[iclus];
            _clusterMaxDRDEta = _clusterDEtaToSeed.get()[iclus];
            _clusterMaxDRRawEnergy = _clusterRawEnergy.get()[iclus];
        }

        subclustersRawEnergy += pclus->energy();
        _clustersMeanDRToSeed   = reco::deltaR(*pclus, *theseed)*pclus->energy();
        _clustersMeanDEtaToSeed = (pclus->eta() - theseed->eta())*pclus->energy();
        _clustersMeanDPhiToSeed = TVector2::Phi_mpi_pi(pclus->phi() - theseed->phi())*pclus->energy();

        _clusterDPhiToGen.get()[iclus] = TVector2::Phi_mpi_pi(pclus->phi() - gen.phi());
        _clusterDEtaToGen.get()[iclus] = pclus->eta() - gen.eta();
        _clusterInMustache.get()[iclus] = (int) MK::inMustache(theseed->eta(),
                theseed->phi(),
                pclus->energy(),
                pclus->eta(),
                pclus->phi());
        _clusterInDynDPhi.get()[iclus] = (int) MK::inDynamicDPhiWindow(
		    theseed->hitsAndFractions().at(0).first.subdetId()==EcalBarrel,
                    theseed->phi(),
                    pclus->energy(),
                    pclus->eta(),
                    pclus->phi());
        //if(IsLinkedByRecHit(theseed, pclus,50., 0.5))
        //    _clusterLeakage.get()[iclus] = 1;
        //else 
        //    _clusterLeakage.get()[iclus] = 0;



        ++iclus;
    }
    // compute mean and rms values
    _clustersMeanRawEnergy /= (double)(_N_ECALClusters+1);
    clustersMeanSquareRawEnergy /= (double)(_N_ECALClusters+1);
    _clustersRMSRawEnergy = sqrt(clustersMeanSquareRawEnergy - _clustersMeanRawEnergy*_clustersMeanRawEnergy);

    if(subclustersRawEnergy>0.)
    {
        _clustersMeanDRToSeed /= subclustersRawEnergy;
        _clustersMeanDEtaToSeed /= subclustersRawEnergy;
        _clustersMeanDPhiToSeed /= subclustersRawEnergy;
    }


    // loop over all preshower clusters 
    auto psclusend = sc.preshowerClustersEnd();
    size_t ipsclus = 0;
    edm::Ptr<reco::CaloCluster> ppsclus;
    for( auto psclus = sc.preshowerClustersBegin(); psclus != psclusend; ++psclus ) 
    {
        ppsclus = *psclus;
        _psClusterRawEnergy.get()[ipsclus] = ppsclus->energy();    
        _psClusterEta.get()[ipsclus]       = ppsclus->eta();    
        _psClusterPhi.get()[ipsclus]       = ppsclus->phi();
        ++ipsclus;
    }

    _tree->Fill();
}

void GedPhotonRegTreeMaker::processNoMatchFillTree(const edm::Event& e, const edm::EventSetup& es, const reco::GenParticle& gen) 
{

    _N_ECALClusters = 0; 
    _N_PSClusters   = 0;
    _hoe   = 0;

    // generated information
    _isMatched       = 0;
    _genEnergy       = gen.energy();
    _genPt           = gen.pt();
    _genEta          = gen.eta();
    _genPhi          = gen.phi();
    _genDEoE         = 999.;
    _genDRToCentroid = 999.;
    _genDRToSeed     = 999.;

    // supercluster information
    _scIndex                += 1;
    _scRawEnergy             = 0;
    _scCalibratedEnergy      = 0;
    _scPreshowerEnergy       = 0;
    _scPreshowerEnergyPlane1 = 0;
    _scPreshowerEnergyPlane2 = 0;
    _scIsEB                  = 0;
    _scEta                   = 0;
    _scPhi                   = 0;
    _scR                     = 0;
    _scPhiWidth              = 0;
    _scEtaWidth              = 0;


    // sc seed information
    _scSeedRawEnergy        = 0;
    _scSeedCalibratedEnergy = 0;
    _scSeedEta              = 0;
    _scSeedPhi              = 0;
    _scSeedSize             = 0;
    _scSeedR9               = 0;
    _scSeedEmax             = 0;
    _scSeedE2nd             = 0;
    _scSeedE2x5max          = 0;
    _scSeedLeftRightAsym    = 0;
    _scSeedTopBottomAsym    = 0;
    _scSeed2x5LeftRightAsym = 0;
    _scSeed2x5TopBottomAsym = 0;
    _scSeedSigmaIetaIeta    = 0;
    _scSeedSigmaIphiIphi    = 0;
    _scSeedSigmaIetaIphi    = 0;
    _scSeedE3x3          = 0;
    _scSeedE5x5          = 0;
    _scSeedELeft          = 0;
    _scSeedERight          = 0;
    _scSeedETop          = 0;
    _scSeedEBottom          = 0;
    _scSeedE2x5Left          = 0;
    _scSeedE2x5Right          = 0;
    _scSeedE2x5Top          = 0;
    _scSeedE2x5Bottom          = 0;

    _scSeedCryEta  = 0;
    _scSeedCryPhi  = 0;
    _scSeedCryIeta = 0;
    _scSeedCryIphi = 0;
    _scSeedCryX    = 0;
    _scSeedCryY    = 0;
    _scSeedCryIx   = 0;
    _scSeedCryIy   = 0;

    _scSeedCryIeta_glob = 0;
    _scSeedCryIphi_glob = 0;
    _scSeedCryIx_glob   = 0;
    _scSeedCryIy_glob   = 0;


    // loop over all clusters that aren't the seed
    setTreeArraysForSize(0, 0);
    _clusterMaxDR     = 999.;
    _clusterMaxDRDPhi = 999.;
    _clusterMaxDRDEta = 999.;
    _clustersMeanDRToSeed = 999.;
    _clustersMeanDEtaToSeed = 999.;
    _clustersMeanDPhiToSeed = 999.;
    _clusterMaxDRRawEnergy      = 0.;
    _clustersMeanRawEnergy      = 0.;
    _clustersRMSRawEnergy       = 0.;
    

    _tree->Fill();
}




void GedPhotonRegTreeMaker::setTreeArraysForSize(const size_t N_ECAL, const size_t N_PS) 
{
    float* cRE_new        = new float[N_ECAL];
    float* cCE_new        = new float[N_ECAL];
    float* cEta_new       = new float[N_ECAL];
    float* cPhi_new       = new float[N_ECAL];
    float* cDPhiSeed_new  = new float[N_ECAL];
    float* cDEtaSeed_new  = new float[N_ECAL];
    float* cDPhiCntr_new  = new float[N_ECAL];
    float* cDEtaCntr_new  = new float[N_ECAL];
    int*   cInMust_new    = new int[N_ECAL];
    int*   cInDynDPhi_new = new int[N_ECAL];
    float* psRE_new       = new float[N_PS];
    float* psEta_new      = new float[N_PS];
    float* psPhi_new      = new float[N_PS];

    _clusterRawEnergy     .reset(cRE_new       , array_deleter<float>());
    _clusterCalibEnergy   .reset(cCE_new       , array_deleter<float>());
    _clusterEta           .reset(cEta_new      , array_deleter<float>());
    _clusterPhi           .reset(cPhi_new      , array_deleter<float>());
    _clusterDPhiToSeed    .reset(cDPhiSeed_new , array_deleter<float>());
    _clusterDEtaToSeed    .reset(cDEtaSeed_new , array_deleter<float>());  
    _clusterDPhiToCentroid.reset(cDPhiCntr_new , array_deleter<float>());
    _clusterDEtaToCentroid.reset(cDEtaCntr_new , array_deleter<float>());
    _clusterInMustache    .reset(cInMust_new   , array_deleter<int>());
    _clusterInDynDPhi     .reset(cInDynDPhi_new, array_deleter<int>());
    _psClusterRawEnergy   .reset(psRE_new      , array_deleter<float>());
    _psClusterEta         .reset(psEta_new     , array_deleter<float>());
    _psClusterPhi         .reset(psPhi_new     , array_deleter<float>());

    _tree->GetBranch("clusterRawEnergy")     ->SetAddress(_clusterRawEnergy.get());
    _tree->GetBranch("clusterCalibEnergy")   ->SetAddress(_clusterCalibEnergy.get());
    _tree->GetBranch("clusterEta")           ->SetAddress(_clusterEta.get());
    _tree->GetBranch("clusterPhi")           ->SetAddress(_clusterPhi.get());
    _tree->GetBranch("clusterDPhiToSeed")    ->SetAddress(_clusterDPhiToSeed.get());
    _tree->GetBranch("clusterDEtaToSeed")    ->SetAddress(_clusterDEtaToSeed.get());
    _tree->GetBranch("clusterDPhiToCentroid")->SetAddress(_clusterDPhiToCentroid.get());
    _tree->GetBranch("clusterDEtaToCentroid")->SetAddress(_clusterDEtaToCentroid.get());
    _tree->GetBranch("clusterInMustache")    ->SetAddress(_clusterInMustache.get());
    _tree->GetBranch("clusterInDynDPhi")     ->SetAddress(_clusterInDynDPhi.get());
    _tree->GetBranch("psClusterRawEnergy")   ->SetAddress(_psClusterRawEnergy.get());
    _tree->GetBranch("psClusterEta")         ->SetAddress(_psClusterEta.get());
    _tree->GetBranch("psClusterPhi")         ->SetAddress(_psClusterPhi.get());

    float* cDPhiGen_new = new float[N_ECAL];
    float* cDEtaGen_new = new float[N_ECAL];
    _clusterDPhiToGen.reset(cDPhiGen_new, array_deleter<float>());
    _clusterDEtaToGen.reset(cDEtaGen_new, array_deleter<float>());
    _tree->GetBranch("clusterDPhiToGen") ->SetAddress(_clusterDPhiToGen.get());
    _tree->GetBranch("clusterDEtaToGen") ->SetAddress(_clusterDEtaToGen.get());
}


bool GedPhotonRegTreeMaker::IsLinkedByRecHit(const edm::Ptr<reco::PFCluster>& the_seed, const edm::Ptr<reco::PFCluster>& x, const double threshold, const double majority)
{

    if( the_seed->energy() < threshold ) 
        return false; 
    const auto& seedHitsAndFractions = the_seed->hitsAndFractions();
    const auto& xHitsAndFractions    = x->hitsAndFractions();      
    float x_rechits_tot              = xHitsAndFractions.size();
    float x_rechits_match            = 0.0;      
    for( const std::pair<DetId, float>& seedHit : seedHitsAndFractions ) 
    {
        for( const std::pair<DetId, float>& xHit : xHitsAndFractions ) 
        {
            if( seedHit.first == xHit.first ) 
            {	    
                x_rechits_match += 1.0;
            }
        }	
    }      
    return x_rechits_match/x_rechits_tot > majority;

}

void GedPhotonRegTreeMaker::matchPhotonToGen(const reco::PhotonCollection& phs, const reco::GenParticleCollection& gens)
{
    for( const auto& gen : gens ) 
    {
        if( abs(gen.pdgId()) == 22 &&  gen.status()==1) // looking at final state photons
        {
            double minDr = 1e6;
            double minDe = 1e6;
            double minDeDr = 1e6;
            double this_dr;
            double this_de;
            double this_dedr;
            const reco::Photon* bestmatch = NULL;
            // loop on EB SCs
            for( const auto& ph: phs) 
            {
                //math::XYZVector direction = ph.position() - gen.vertex();
                //double energy = ph.energy();
                //math::XYZVector momentum = direction.unit() * energy;
                //math::XYZTLorentzVector physicalSCP4(momentum.x(), momentum.y(), momentum.z(), energy );
                // TODO: is photon momentum already corrected for vertex?
                this_dr = reco::deltaR(gen,ph);
                this_de = fabs(gen.energy()-ph.energy())/gen.energy();
                this_dedr = sqrt(this_dr*this_dr + this_de*this_de);
                //if(this_dr < _matchMaxDR && sc.energy()/gen.energy()>_matchMinERatio && (this_dr < minDr || this_de<minDe))
                if(this_dr < _matchMaxDR && ph.energy()/gen.energy()>_matchMinERatio)
                {
                    if( (_matchingType==0 && this_dr<minDr) || 
                            (_matchingType==1 && this_de<minDe) ||
                            (_matchingType==2 && this_dedr<minDeDr))
                    {
                        minDr = this_dr;
                        minDe = this_de;
                        minDeDr = this_dedr;
                        bestmatch = &ph;
                    }
                }
            }
            if(bestmatch)
            {
                _matchedPhotons.push_back(*bestmatch);
                _matchedParticles.push_back( std::make_pair(gen, (int)_matchedPhotons.size()-1) );
            }
            else
            {
                _matchedParticles.push_back( std::make_pair(gen, -1) );
            }
        }
    }
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GedPhotonRegTreeMaker);
