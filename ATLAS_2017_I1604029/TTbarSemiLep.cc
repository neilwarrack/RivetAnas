#include "Sacrifice/UserHooksFactory.hh"
#include "Pythia8/UserHooks.h"
#include "Pythia8/PhaseSpace.h"

namespace Pythia8 {


  struct TTbarSemiLep : public UserHooks {

    bool canVetoResonanceDecays() { return true; }

    bool doVetoResonanceDecays(Event& process) {
      int ntlep = 0, nthad = 0;
      for (int i = 1; i < process.size(); ++i) {
        Particle& p = process[i];
        if (p.idAbs() != 6) continue; 
        // cout << i << ": PID=" << p.id() << " status=" << p.status() << " d1=" << p.daughter1() << ", d2=" << p.daughter2() << endl;
        Particle* w = nullptr;
        if (p.daughter1() != 0 && process[p.daughter1()].idAbs() == 24) w = &process[p.daughter1()];
        if (p.daughter2() != 0 && process[p.daughter2()].idAbs() == 24) w = &process[p.daughter2()];
        if (w == nullptr) continue;
        Particle& wd1 = process[w->daughter1()];
        Particle& wd2 = process[w->daughter2()];
        // cout << wd1.id() << " " << wd2.id() << endl;
        if (wd1.isLepton() && wd2.isLepton()) ntlep += 1;
        if (wd1.isQuark() && wd2.isQuark()) nthad += 1;
      }
      if (ntlep + nthad != 2) return false;
      return ntlep != nthad;
    }

  };

  Sacrifice::UserHooksFactory::Creator<Pythia8::TTbarSemiLep> TTBARSEMILEP("TTbarSemiLep");

}
