import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os

from EFTGenReader.GenReader.DatasetHelper import DatasetHelper

options = VarParsing.VarParsing('analysis')

# Setup and register default options
options.maxEvents = 10
options.register("singleFile","",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,"name of a single root file")
options.register("dataset","central_ttH",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,"name of the dataset as it appears in the JSON file")
options.register("test",False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool, "changes the output name to a dummy value")
options.register("debug",False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool, "run in debug mode")
options.register("normType",1,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,"how to normalize the histograms; 0 - no norm, 1 - unit norm (default), 2 - xsec norm")
options.register("intgLumi",1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"intg. lumi to scale the histograms to (no effect for unit norm mode)")
options.register("fnSuffix","_output_tree",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,"string to append to the end of the output root file")
options.register("minPtJet",-1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max pt cut for genjets")
options.register("maxEtaJet",-1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max eta cut for genjets (-1 means no cut)")
options.register("minPtLep",-1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max pt cut for genleptons")
options.register("maxEtaLep",-1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max eta cut for genleptons (-1 means no cut)")

# Get and parse the command line arguments
options.parseArguments()

nd_redirect = "root://ndcms.crc.nd.edu/"
fnal_redirect = "root://cmsxrootd.fnal.gov/"
global_redirect = "root://cms-xrd-global.cern.ch/"

cmssw_base_dir = os.environ['CMSSW_BASE']
dataset_fpath = os.path.join(cmssw_base_dir,"src/EFTGenReader/GenReader/data/JSON/datasets.json")

ds_helper = DatasetHelper()
ds_helper.load(dataset_fpath)
ds_helper.root_redirect = nd_redirect

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Demo", eras.Run2_2017)

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents) # number of events
)

process.MessageLogger.cerr.FwkReport.reportEvery = 500

ds_name   = options.dataset
files     = ds_helper.getFiles(ds_name)
is_eft    = ds_helper.getData(ds_name,'is_eft')
xsec_norm = ds_helper.getData(ds_name,'central_xsec')
datatier  = ds_helper.getData(ds_name,'datatier')


file_name_in = options.singleFile
if file_name_in != "":
    out_fname = "EFTMaodHists_output_tree.root"
    files = ["file:"+file_name_in]
    is_eft = True
elif options.test:
    out_fname = "TEST_output_tree.root"
else:
    out_fname = "%s%s.root" % (ds_name,options.fnSuffix)

out_path = os.path.join("output",out_fname)

print "Using Sample: %s" % (ds_name)
print "Save output to: %s" % (out_path)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        *files
    )
)

process.load("EFTGenReader.GenReader.EFTMaodHists_cfi")

process.EFTMaodHists.debug     = options.debug
process.EFTMaodHists.norm_type = options.normType      # 0 - No norm, 1 - unit norm, 2 - xsec norm
process.EFTMaodHists.intg_lumi = options.intgLumi
process.EFTMaodHists.gp_events = 500

# Cut settings
process.EFTMaodHists.min_pt_jet = options.minPtJet
process.EFTMaodHists.min_pt_lep = options.minPtLep
process.EFTMaodHists.max_eta_jet = options.maxEtaJet
process.EFTMaodHists.max_eta_lep = options.maxEtaLep

process.EFTMaodHists.iseft = is_eft
process.EFTMaodHists.xsec_norm = xsec_norm

if datatier == "MINIAODSIM":
    process.EFTMaodHists.GenParticles = cms.InputTag("prunedGenParticles")
    process.EFTMaodHists.GenJets      = cms.InputTag("slimmedGenJets")
    process.EFTMaodHists.PatElectrons = cms.InputTag("slimmedElectrons")
    process.EFTMaodHists.PatMuons     = cms.InputTag("slimmedMuons")
    process.EFTMaodHists.PatJets      = cms.InputTag("slimmedJets")
else:
    print "[ERROR] Unknown datatier: {}".format(datatier)
    raise RuntimeError

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(out_path)
)

process.p = cms.Path(process.EFTMaodHists)

# summary
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)
