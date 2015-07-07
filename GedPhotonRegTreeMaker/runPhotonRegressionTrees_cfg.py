import FWCore.ParameterSet.Config as cms

process = cms.Process("TUPLE")
process.load("Configuration.StandardSequences.Services_cff")
#process.load('Configuration.StandardSequences.Geometry_cff')
#process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

from Configuration.AlCa.autoCond import autoCond 
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('file:SamplesTest/photonRun2_bx25.root')
#fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/Fall14DR/DoublePhotonNoMaterial_FlatPt-0p01To100/GEN-SIM-RAW-RECO/Flat20to50bx25_NoMaterial_MCRUN2_72_V3-v1/00000/002E8190-A88B-E411-AF43-0025905A60E4.root') #run2 photon gun PU bx25
    )


# Event output
#process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('RegressionPhoton_ntuple.root')
    )

### ntuple makers
genParticleID = cms.int32(22)
matchMinERatio = cms.double(0.25)
matchMaxDR = cms.double(0.3)
matchingType = cms.int32(2)

#process.mustacheSCTree = cms.EDAnalyzer(
    #"PFSuperClusterRegTreeMaker",
    #doGen = cms.untracked.bool(True),
    #genSrc = cms.InputTag("genParticles"),
    #genParticleID = genParticleID,
    #matchMinERatio = matchMinERatio,
    #matchMaxDR = matchMaxDR,
    #matchingType = matchingType,
    ##superClusterSrcEB = cms.InputTag('particleFlowSuperClusterECALMustache', 'particleFlowSuperClusterECALBarrel'),
    ##superClusterSrcEE = cms.InputTag('particleFlowSuperClusterECALMustache', 'particleFlowSuperClusterECALEndcapWithPreshower'),
    #superClusterSrcEB = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALBarrel','reRECO'),
    #superClusterSrcEE = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALEndcapWithPreshower','reRECO'),
    #ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    #eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    #VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    #)

process.gedPhotonTree = cms.EDAnalyzer(
    "GedPhotonRegTreeMaker",
    doGen = cms.untracked.bool(True),
    genSrc = cms.InputTag("genParticles"),
    matchMinERatio = matchMinERatio,
    matchMaxDR = matchMaxDR,
    matchingType = matchingType,
    photonSrc = cms.InputTag('gedPhotons', ''),
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    )


#process.egSCTree = cms.EDAnalyzer(
    #"EGSuperClusterRegTreeMaker",
    #doGen = cms.untracked.bool(True),
    #genSrc = cms.InputTag("genParticles"),
    #genParticleID = genParticleID,
    #matchMinERatio = matchMinERatio,
    #matchMaxDR = matchMaxDR,
    #matchingType = matchingType,
    #superClusterSrcEB = cms.InputTag('correctedHybridSuperClusters', '','reRECO'),
    #superClusterSrcEE = cms.InputTag('correctedMulti5x5SuperClustersWithPreshower', '','reRECO'),
    #ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    #eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    #VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    #)

#process.regSeq = cms.Sequence( process.egSCTree +
#process.gedPhotonTree +
#process.mustacheSCTree )

process.regSeq = cms.Sequence( process.gedPhotonTree )

process.regPath = cms.Path(process.regSeq)

process.schedule = cms.Schedule(process.regPath)
