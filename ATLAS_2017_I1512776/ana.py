# Neil's notes: 
# this file is from Andy Buckley
# to use it do:
# setupATLAS then do asetup 21.6.2,AthGeneration then, to run this, do:
# athena ana.py -c"FILE='/path/to/your/EVNT.12345._54321.pool.root.1'"

#import sys
#FILE = sys.argv[1] if len(sys.argv) > 1 else "foo.EVNT.pool.root"
#print "******", FILE, "*******"
assert "FILE" in dir()


import ROOT
from AthenaPython import PyAthena
import McParticleEvent.Pythonizations


class Ana(PyAthena.Alg):

    def __init__(self, name="Ana"):
        PyAthena.Alg.__init__(self, name=name)

    def initialize(self):
        self.sg = PyAthena.py_svc('StoreGateSvc')
        self.rootfile = ROOT.TFile("out.root", "RECREATE")
        self.h_pids = ROOT.TH1F("pids", "", 1201, -600.5, 600.5)
        self.h_pids12 = ROOT.TH1F("pids12", "", 1201, -600.5, 600.5)
        return PyAthena.StatusCode.Success

    def execute(self):
        mcevts = self.sg["GEN_EVENT"]
        for evt in mcevts:
            for p in evt.particles:
                pid = p.pdg_id()
                st = p.status()
                if abs(pid) == 5:
                    print "Neil - found b quark in state: ", st
                self.h_pids.Fill(pid)
                if st in [1,2]:
                    self.h_pids12.Fill(pid)
                    if abs(pid) >= 500 and abs(pid) < 600:
                        print "B-meson: {}, status={}".format(pid, st)


            break #< don't look at pile-up events
        return PyAthena.StatusCode.Success

    def finalize(self):
        self.rootfile.Write()
        return PyAthena.StatusCode.Success


###############

theApp.EvtMax = 1000
import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = [FILE]

from AthenaCommon.AlgSequence import AlgSequence
topSeq = AlgSequence()
topSeq += Ana()
