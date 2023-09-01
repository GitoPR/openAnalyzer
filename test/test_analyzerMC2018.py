import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from Configuration.AlCa.GlobalTag import GlobalTag
process = cms.Process("Demo")


# intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# **********************************************************************
# set the maximum number of events to be processed                     *
#    this number (argument of int32) is to be modified by the user     *
#    according to need and wish                                        *
#    default is preset to -1 (all events)                              *
# **********************************************************************
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(30) )

# set the number of events to be skipped (if any) at end of file below

# define JSON file for 2012 data (not needed for MC)
#goodJSON = '/home/cms-opendata/CMSSW_5_3_32/src/Demo/DemoAnalyzer/datasets/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'

#myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

# ****************************************************************************
# define the input data set here by inserting the appropriate .txt file list *
# ****************************************************************************
import FWCore.Utilities.FileUtils as FileUtils

#print('Check1')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
#print('Check2')
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_design_v9', '')
#globaltag for 2018 MC
process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/102X_upgrade2018_design_v9.db')
#process.GlobalTag.globaltag = '102X_upgrade2018_design_v9'
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")

#print('Check3')                                                                                                                                                                               
#
# ********************************************************************
# load the data set                                                  * 
# this example uses the 2012 Higgs->4lepton simulated dataset        *
# ********************************************************************
#
# *** 2012 Higgs->4lepton simulated data set (299973 events) ***
#files2012data = FileUtils.loadListFromFile ('/home/cms-opendata/CMSSW_5_3_32/src/Demo/DemoAnalyzer/MCsets/CMS_MonteCarlo2012_Summer12_DR53X_SMHiggsToZZTo4L_M-125_8TeV-powheg15-JHUgenV3-pythia6_AODSIM_PU_S10_START53_V19-v1_10000_file_index.txt')
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(*files2012data    
#    )
#)
#
# to speed up, read only first file with 7499 events
##IO: LOAD ROOT FILES HERE
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "file:/hdfs/store/user/ekauffma/CICADA/ZeroBias/CICADA_2023RunC_ZB_21Jul2023/230721_084757/0000/output_1-366.root" 

    )
)

#print('Check1')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.StandardSequences.Services_cff')
#print('Check2')
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START42_V17B.db')
#process.GlobalTag.globaltag = 'START42_V17B::All'
#print('Check3')

# apply JSON file (not for MC)
#   (needs to be placed *after* the process.source input file definition!)
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)

# *************************************************
# number of events to be skipped (0 by default)   *
# *************************************************
process.source.skipEvents = cms.untracked.uint32(0)

##IO: RUN YOUR ANALYZER
process.demo = cms.EDAnalyzer('compAnalyzer', 
                              pfcands  = cms.InputTag("packedPFCandidates", "", "RECO"),
                              L1CaloRegion = cms.InputTag("simCaloStage2Layer1Digis" , "", "NTUPLIZE"), 

                              
    #                          UCTRegion = cms.InputTag('uct2016EmulatorDigis'),
   #                           folderName = cms.untracked.string("Stage3Regions"),
  #                            hcalDigis  = cms.InputTag( 'l1tCaloLayer1Digis'),
 #                             ecalDigis = cms.InputTag( 'l1tCaloLayer1Digis'),
#                              vertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
)


# ***********************************************************
# output file name                                          *
# default is openAnalyzer.root                              *
# ***********************************************************
##IO: OUTPUT DATAFILE GENERATED
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('/nfs_scratch/jorgeeh/compAnalyzer_MC2018.root')
                                   )

process.p = cms.Path(process.demo)
