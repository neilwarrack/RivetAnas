
theApp.EvtMax = -1

import os
import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = [
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00001.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00002.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00003.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00004.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00005.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00006.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00007.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00008.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00009.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._000010.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._000011.pool.root",

# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00002.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00005.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00006.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00015.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00020.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00021.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00024.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00040.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00047.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00048.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00052.pool.root",
]


from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from Rivet_i.Rivet_iConf import Rivet_i

rivet = Rivet_i()
rivet.AnalysisPath = os.environ['PWD']
rivet.Analyses += [ 'ATLAS_2017_I1512776' ]
rivet.RunName = ""
rivet.HistoFile = "MG5_atNLO_Herwig"
rivet.OutputLevel = DEBUG
#rivet.CrossSection = 26.587 #<- xSec from https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TopMC12SingleTopSamples#Singletop_Samples 
rivet.CrossSection = 28.4348 #<- with kFactor
job += rivet
