include("GeneratorUtils/StdAnalysisSetup.py")
theApp.EvtMax = 5000000

import os
import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = [
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00050.pool.root",
]


from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from Rivet_i.Rivet_iConf import Rivet_i

rivet = Rivet_i()
rivet.AnalysisPath = os.environ['PWD']
#rivet.Analyses += [ 'ATLAS_2017_I1512776' ]
rivet.Analyses += [ 'MC_tchan_particle' ]
rivet.RunName = ""
#rivet.HistoFile = "MG5_atNLO_Herwig"
rivet.HistoFile = "Acer"
#rivet.OutputLevel = DEBUG
#rivet.CrossSection = 26.587 #<- xSec from https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TopMC12SingleTopSamples#Singletop_Samples 
rivet.CrossSection = 28.4332127 #<- (Acer - 110069) with kFactor
#rivet.CrossSection = 27.337 #<- (Acer - 110069) without kFactor
job += rivet
theApp.EvtMax = 500000
