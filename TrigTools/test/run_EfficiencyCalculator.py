import FWCore.ParameterSet.Config as cms
import glob
import os

process = cms.Process("USER")

reRunL1 = False

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_dataRun3_Prompt_v4') 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


dir1Name = "/eos/cms/store/group/phys_egamma/ec/rsalvatico/MiniHoE/"#test/"
fileList1 = filter(os.path.isfile, glob.glob(dir1Name + "*.root"))
fList = []
file_counter = 0
for f in fileList1:
    fs = str(f).replace("/eos/","file:/eos/")
    fList.append(fs)
    file_counter += 1
    #if file_counter == 200: break

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(fList
                                    #'file:/eos/cms/store/group/phys_egamma/ec/rsalvatico/MiniHoE/stepMINI_7052362_0.root'
                                    #'file:/afs/cern.ch/user/r/rselvati/work/private/MyEffiPackage/CMSSW_14_0_14/src/stepMINIModified.root'
                )
                            )
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
#Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("Efficiency_Standard.root")
)


process.hltEgammaCandidates = cms.EDProducer( "EgammaHLTRecoEcalCandidateProducers",
    scHybridBarrelProducer = cms.InputTag( 'hltParticleFlowSuperClusterECALL1Seeded','hltParticleFlowSuperClusterECALBarrel' ),
    scIslandEndcapProducer = cms.InputTag( 'hltParticleFlowSuperClusterECALL1Seeded','hltParticleFlowSuperClusterECALEndcapWithPreshower' ),
    recoEcalCandidateCollection = cms.string( "" )
)

process.GlobalParametersRcdSource = cms.ESSource("EmptyESSource",
    recordName = cms.string('L1TGlobalParametersRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.GlobalParameters = cms.ESProducer("StableParametersTrivialProducer", 
    # bx in event
    #NumberBxInEvent = cms.int32(5),
                                    
    # trigger decision
                                            
    # number of physics trigger algorithms
    NumberPhysTriggers = cms.uint32(512),
                                              
                                              
    # trigger objects

    # muons
    NumberL1Muon = cms.uint32(8),
    
    # e/gamma and isolated e/gamma objects
    NumberL1EGamma = cms.uint32(12),

    #  jets
    NumberL1Jet = cms.uint32(12),

    # taus
    NumberL1Tau = cms.uint32(12),

    # hardware
                                            
    # number of maximum chips defined in the xml file
    NumberChips = cms.uint32(1),

    # number of pins on the GTL condition chips
    PinsOnChip = cms.uint32(512),

    # correspondence "condition chip - GTL algorithm word" in the hardware
    # e.g.: chip 2: 0 - 95;  chip 1: 96 - 128 (191)
    OrderOfChip = cms.vint32(1),
)

process.hltGtStage2Digis = cms.EDProducer( "L1TRawToDigi",
    CTP7 = cms.untracked.bool( False ),
    DmxFWId = cms.uint32( 0 ),
    FWId = cms.uint32( 0 ),
    FWOverride = cms.bool( False ),
    FedIds = cms.vint32( [ 1404 ] ),
    InputLabel = cms.InputTag( "rawDataCollector" ),
    MTF7 = cms.untracked.bool( False ),
    MinFeds = cms.uint32( 0 ),
    Setup = cms.string( "stage2::GTSetup" ),
    TMTCheck = cms.bool( True ),
    debug = cms.untracked.bool( False ),
    lenAMC13Header = cms.untracked.int32( 8 ),
    lenAMC13Trailer = cms.untracked.int32( 8 ),
    lenAMCHeader = cms.untracked.int32( 8 ),
    lenAMCTrailer = cms.untracked.int32( 0 ),
    lenSlinkHeader = cms.untracked.int32( 8 ),
    lenSlinkTrailer = cms.untracked.int32( 8 )
)

process.hltGtStage2ObjectMap = cms.EDProducer( "L1TGlobalProducer",
    AlgoBlkInputTag = cms.InputTag( "hltGtStage2Digis" ),
    AlgorithmTriggersUnmasked = cms.bool( True ),
    AlgorithmTriggersUnprescaled = cms.bool( True ),
    AlternativeNrBxBoardDaq = cms.uint32( 0 ),
    BstLengthBytes = cms.int32( -1 ),
    EGammaInputTag = cms.InputTag( "hltGtStage2Digis","EGamma" ),
    EmulateBxInEvent = cms.int32( 1 ),
    EtSumInputTag = cms.InputTag( "hltGtStage2Digis","EtSum" ),
    ExtInputTag = cms.InputTag( "hltGtStage2Digis" ),
    GetPrescaleColumnFromData = cms.bool( False ),
    JetInputTag = cms.InputTag( "hltGtStage2Digis","Jet" ),
    L1DataBxInEvent = cms.int32( 5 ),
    MuonInputTag = cms.InputTag( "hltGtStage2Digis","Muon" ),
    MuonShowerInputTag = cms.InputTag( "hltGtStage2Digis","MuonShower" ),
    PrescaleSet = cms.uint32( 1 ),
    PrintL1Menu = cms.untracked.bool( False ),
    ProduceL1GtDaqRecord = cms.bool( True ),
    ProduceL1GtObjectMapRecord = cms.bool( True ),
    RequireMenuToMatchAlgoBlkInput = cms.bool( True ),
    TauInputTag = cms.InputTag( "hltGtStage2Digis","Tau" ),
    TriggerMenuLuminosity = cms.string( "startup" ),
    Verbosity = cms.untracked.int32( 0 ),
    resetPSCountersEachLumiSec = cms.bool( True ),
    semiRandomInitialPSCounters = cms.bool( False ),
    useMuonShowers = cms.bool( True )
)

#process.hltL1sSingleEGor = cms.EDFilter( "HLTL1TSeed",
#    saveTags = cms.bool( True ),
#    L1SeedsLogicalExpression = cms.string( "L1_SingleEG38er2p5 OR L1_SingleEG40er2p5 OR L1_SingleEG42er2p5 OR L1_SingleEG45er2p5 OR L1_SingleEG60 OR L1_SingleEG34er2p5 OR L1_SingleEG36er2p5 OR L1_SingleIsoEG30er2p1 OR L1_SingleIsoEG32er2p1 OR L1_SingleIsoEG30er2p5 OR L1_SingleIsoEG32er2p5 OR L1_SingleIsoEG34er2p5"), #opt4
#    L1ObjectMapInputTag = cms.InputTag( "hltGtStage2ObjectMap" ),
#    L1GlobalInputTag = cms.InputTag( "hltGtStage2Digis" ),
#    L1MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
#    L1MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
#    L1EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
#    L1JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
#    L1TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
#    L1EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' )
#)

process.hltL1sSingleEGorModified = cms.EDFilter( "HLTL1TSeed",
    saveTags = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleEG60" ), #opt4
    L1ObjectMapInputTag = cms.InputTag( "hltGtStage2ObjectMap" ),
    L1GlobalInputTag = cms.InputTag( "hltGtStage2Digis" ),
    L1MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    L1MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
    L1EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
    L1JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
    L1TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
    L1EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' )
)

# process.hltL1sSingleAndDoubleEG = cms.EDFilter( "HLTL1TSeed",
#     saveTags = cms.bool( True ),
#     L1SeedsLogicalExpression = cms.string( "L1_SingleEG26er2p5 OR L1_SingleEG38er2p5 OR L1_SingleEG40er2p5 OR L1_SingleEG42er2p5 OR L1_SingleEG45er2p5 OR L1_SingleEG60 OR L1_SingleEG34er2p5 OR L1_SingleEG36er2p5 OR L1_SingleIsoEG30er2p1 OR L1_SingleIsoEG32er2p1 OR L1_SingleIsoEG30er2p5 OR L1_SingleIsoEG32er2p5 OR L1_SingleIsoEG34er2p5 OR L1_DoubleEG_22_10_er2p5 OR L1_DoubleEG_25_12_er2p5 OR L1_DoubleEG_25_14_er2p5 OR L1_DoubleEG_LooseIso22_12_er2p5" ),
#     L1ObjectMapInputTag = cms.InputTag( "hltGtStage2ObjectMap" ),
#     L1GlobalInputTag = cms.InputTag( "hltGtStage2Digis" ),
#     L1MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
#     L1MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
#     L1EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
#     L1JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
#     L1TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
#     L1EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' )
# )

# process.hltL1sSingleEG40to50 = cms.EDFilter( "HLTL1TSeed",
#     saveTags = cms.bool( True ),
#     L1SeedsLogicalExpression = cms.string( "L1_SingleEG40er2p5 OR L1_SingleEG42er2p5 OR L1_SingleEG45er2p5"),
#     L1ObjectMapInputTag = cms.InputTag( "hltGtStage2ObjectMap" ),
#     L1GlobalInputTag = cms.InputTag( "hltGtStage2Digis" ),
#     L1MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
#     L1MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
#     L1EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
#     L1JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
#     L1TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
#     L1EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' )
# )

# process.hltL1sSingleJet170IorSingleJet180IorSingleJet200 = cms.EDFilter( "HLTL1TSeed",
#     saveTags = cms.bool( True ),
#     L1SeedsLogicalExpression = cms.string( "L1_SingleJet180 OR L1_SingleJet200"),
#     L1ObjectMapInputTag = cms.InputTag( "hltGtStage2ObjectMap" ),
#     L1GlobalInputTag = cms.InputTag( "hltGtStage2Digis" ),
#     L1MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
#     L1MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
#     L1EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
#     L1JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
#     L1TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
#     L1EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' )
# )

# process.hltL1sSingleEGNonIsoOrWithJetAndTau = cms.EDFilter( "HLTL1TSeed",
#     saveTags = cms.bool( True ),
#     L1SeedsLogicalExpression = cms.string( "L1_SingleEG34er2p5 OR L1_SingleEG36er2p5 OR L1_SingleEG38er2p5 OR L1_SingleEG40er2p5 OR L1_SingleJet180 OR L1_SingleJet200 OR L1_SingleTau120er2p1 OR L1_SingleTau130er2p1 OR L1_SingleEG42er2p5 OR L1_SingleEG45er2p5 OR L1_SingleEG60"),
#     L1ObjectMapInputTag = cms.InputTag( "hltGtStage2ObjectMap" ),
#     L1GlobalInputTag = cms.InputTag( "hltGtStage2Digis" ),
#     L1MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
#     L1MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
#     L1EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
#     L1JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
#     L1TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
#     L1EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' )
# )

# process.hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau = cms.EDFilter( "HLTL1TSeed",
#     saveTags = cms.bool( True ),
#     L1SeedsLogicalExpression = cms.string( "L1_SingleEG34er2p5 OR L1_SingleEG36er2p5 OR L1_SingleEG38er2p5 OR L1_SingleEG40er2p5 OR L1_SingleEG42er2p5 OR L1_SingleEG45er2p5 OR L1_SingleEG60 OR L1_DoubleEG_22_10_er2p5 OR L1_DoubleEG_25_12_er2p5 OR L1_DoubleEG_25_14_er2p5 OR L1_SingleJet180 OR L1_SingleJet200 OR L1_SingleTau120er2p1 OR L1_SingleTau130er2p1"),
#     L1ObjectMapInputTag = cms.InputTag( "hltGtStage2ObjectMap" ),
#     L1GlobalInputTag = cms.InputTag( "hltGtStage2Digis" ),
#     L1MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
#     L1MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
#     L1EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
#     L1JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
#     L1TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
#     L1EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' )
# )

# # process.hltL1sSingleEGorFilter = cms.EDFilter( "HLTEgammaL1TMatchFilterRegional",
# #     saveTags = cms.bool( True ),
# #     candIsolatedTag = cms.InputTag( "hltEgammaCandidates" ),
# #     l1IsolatedTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
# #     candNonIsolatedTag = cms.InputTag( "" ),
# #     l1NonIsolatedTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
# #     L1SeedFilterTag = cms.InputTag( "hltL1sSingleEGor" ),
# #     l1CenJetsTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
# #     l1TausTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
# #     ncandcut = cms.int32( 1 ),
# #     doIsolated = cms.bool( False ),
# #     region_eta_size = cms.double( 0.522 ),
# #     region_eta_size_ecap = cms.double( 1.0 ),
# #     region_phi_size = cms.double( 1.044 ),
# #     barrel_end = cms.double( 1.4791 ),
# #     endcap_end = cms.double( 2.65 )
# # )


process.EfficiencyCalculator = cms.EDAnalyzer('EfficiencyCalculator',
                                              stageL1Trigger = cms.uint32(2)
)



if reRunL1:
    process.myTask = cms.Task(process.GlobalParametersRcdSource,process.GlobalParameters,process.hltGtStage2Digis,process.hltGtStage2ObjectMap)
    process.t = cms.Path(process.myTask)
    process.p = cms.Path(process.hltL1sSingleEGorModified + process.EfficiencyCalculator)
    #process.p = cms.Path(process.hltL1sSingleAndDoubleEG + process.EfficiencyCalculator)
    #process.p = cms.Path(process.hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau + process.EfficiencyCalculator)
    #process.p = cms.Path(process.hltL1sSingleEGNonIsoOrWithJetAndTau + process.EfficiencyCalculator)
    #process.p = cms.Path(process.hltL1sSingleEG40to50 + process.EfficiencyCalculator)
    process.schedule = cms.Schedule(process.t,process.p)
else:
    process.p = cms.Path(process.EfficiencyCalculator)
    process.schedule = cms.Schedule(process.p)
