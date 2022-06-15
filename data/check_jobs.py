#! /usr/bin/env python

import os
import sys
import optparse
import datetime
import subprocess
import ROOT

from glob import glob

usage = "usage: To be run from RootTreeAnalyzer: python scripts/check_jobs.py -i /eos/cms/store/group/phys_exotica/lq-LQ-lq/rootReducedTrees/test_prod1 -l output/test_prod1 --merge --match QCD"

parser = optparse.OptionParser(usage)

parser.add_option("-i", "--input", dest="input",
                  help="batch directory to be checked", default = "/eos/user/c/cbasile/HLTapplyDATA")

parser.add_option("-o", "--output", dest="output",
                  help="directory to save the merged files", default = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/data")
#parser.add_option("-l", "--loginput", dest="loginput",
#                  help="the directory contains the logs to be checked")

parser.add_option("--merge", action="store_true", dest="doMerging")

parser.add_option("--match", dest="matchString")

(opt, args) = parser.parse_args()
#print(opt)

#if not opt.input:   
#    parser.error('input directory not provided')

#if not opt.loginput:   
#    parser.error('login directory not provided')

################################################

condordir = opt.input
mergedir = opt.output +"/"+"merged"

outTreeName = "B0sidebands"
mergeOutFileName = "BKG_MergData17.root"
mergeOutHistName = "HLTapply_MergData17.root"

builtList = True
HistoList = []

os.system("mkdir -p "+mergedir)

subDirList = next(os.walk(condordir))[2]
print subDirList[1]

TotEvents = 0
SampleEvs = 0
chain = ROOT.TChain("B0sidebands")
haddHistoString = "hadd " + mergedir + "/" + mergeOutHistName 

## Loop over datasets
for sample in subDirList:

	jobOutPath = condordir + "/" + sample
	jobOutFile = ROOT.TFile.Open(jobOutPath)	

	## --> TChain with output trees
	if (sample.endswith("_tree.root") ) : 

		dataTreePath = jobOutPath + "/" + outTreeName
		dataTree = jobOutFile.Get(outTreeName)
		TotEvents += dataTree.GetEntries() 
		print ("  pick " + sample + " --> EV " + str(dataTree.GetEntriesFast()))
		chain.Add(dataTreePath)
		print (" --> TOT EVENTS " + str(chain.GetEntries()))

		jobOutFile.Close()
	## --> add together all histograms
	if (sample.endswith("_histo.root") ):

		haddHistoString += " " + jobOutPath
		#print	haddHistoString

mergeOutFilePath = mergedir + "/" + mergeOutFileName

mergeOutFile = ROOT.TFile.Open(mergeOutFilePath, "RECREATE")
chain.Write()
mergeOutFile.Close()


os.system(haddHistoString)

