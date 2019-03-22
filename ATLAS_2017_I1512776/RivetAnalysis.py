include("GeneratorUtils/StdAnalysisSetup.py")
theApp.EvtMax = 5000000

import os
import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = [
#"/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00002.pool.root",
#"/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00005.pool.root",
#"/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00006.pool.root",
#"/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00015.pool.root",
#"/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00020.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00021.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00024.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00040.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00047.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00048.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV.evgen/110121.Herwig_singletop_tchan2to3nlo_lept_CT10f4_8TeV._00052.pool.root",

"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00001.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00002.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00003.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00004.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00005.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00006.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00007.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00008.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00009.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00010.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00011.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00012.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00013.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00014.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00015.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00016.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00017.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00018.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00019.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00020.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00021.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00022.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00023.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00024.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00025.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00026.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00027.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00028.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00029.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00030.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00031.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00032.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00033.pool.root",
"/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00034.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00035.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00036.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00037.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00038.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00039.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00040.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00041.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00042.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00043.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00044.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00045.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00046.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00047.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00048.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00049.pool.root",
# "/nfs/atlas/nwarrack/user.dhirsch.AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV.evgen/AcerMC39.110069.singletop_tchan_lept_scale_115_CTEQ6L1_8TeV._00050.pool.root",
]


from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from Rivet_i.Rivet_iConf import Rivet_i

rivet = Rivet_i()
rivet.AnalysisPath = os.environ['PWD']
rivet.Analyses += [ 'ATLAS_2017_I1512776' ]
rivet.RunName = ""
#rivet.HistoFile = "MG5_atNLO_Herwig"
rivet.HistoFile = "Acer"
#rivet.OutputLevel = DEBUG
#rivet.CrossSection = 26.587 #<- xSec from https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TopMC12SingleTopSamples#Singletop_Samples 
rivet.CrossSection = 28.4348 #<- with kFactor
job += rivet
theApp.EvtMax = 500000
