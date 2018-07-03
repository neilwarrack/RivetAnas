#include "Sacrifice/UserHooksFactory.hh"
#include "Pythia8/UserHooks.h"
#include "Pythia8/PhaseSpace.h"

namespace Pythia8 {


  struct WW_to_Es_and_MUs : public UserHooks {

    bool canVetoResonanceDecays() {cout << "HOOK: CAN_VETO_RESONANCES!"; return true; }

    bool doVetoResonanceDecays(Event& process) {
     


      int nT = 0, nH = 0, nE_MU_TAU = 0, nEs_and_MUs = 0, ntE_MU = 0, ntbarE_MU=0;
      bool tFlag = false;
      bool tbarFlag = false;


      for (int i = 1; i < process.size(); ++i) {

        Particle& p = process[i];
	if (p.idAbs() !=6 ) continue; // find tops 
        if (p.id() == 6) { //top

	tFlag = true;

        Particle* tw = nullptr;
        if (p.daughter1() != 0 && process[p.daughter1()].idAbs() == 24) tw = &process[p.daughter1()];
     	if (p.daughter2() != 0 && process[p.daughter2()].idAbs() == 24) tw = &process[p.daughter2()];
     	if (tw == nullptr) continue;
        Particle& twd1 = process[tw->daughter1()];
        Particle& twd2 = process[tw->daughter2()];
	if (twd1.idAbs() == 11 || twd1.idAbs() == 13 || twd2.idAbs() == 11 || twd2.idAbs() == 13){
	  ntE_MU++;
	}
	} else { // antitop

	  tbarFlag = true;

	  Particle* tbarw = nullptr;
	  if (p.daughter1() != 0 && process[p.daughter1()].idAbs() == 24) tbarw = &process[p.daughter1()];
	  if (p.daughter2() != 0 && process[p.daughter2()].idAbs() == 24) tbarw = &process[p.daughter2()];
	  if (tbarw == nullptr) continue;
	  Particle& tbarwd1 = process[tbarw->daughter1()];
	  Particle& tbarwd2 = process[tbarw->daughter2()];
	  if (tbarwd1.idAbs() == 11 || tbarwd1.idAbs() == 13 || tbarwd2.idAbs() == 11 || tbarwd2.idAbs() == 13){
	    ntbarE_MU++;
	  }
	}
      }
      

      if (tFlag && tbarFlag){ // then a top and anti top are found.
	if ( !(ntE_MU == 1 && ntbarE_MU == 1) ) { // then there isn't two leptonic decays
	 
	  cout << "." ;
          return ( !(ntE_MU == 1 && ntbarE_MU == 1) ) ; // redecay
        } else {
          cout << "userhook success!" << endl;
	  return ( !(ntE_MU == 1 && ntbarE_MU == 1) ) ; // accept decay
        }
      } else {
	cout << "this hook requires only ttbar decays, no single tops! Something's not going to work here!" << endl;
	return false;
      }
    
 
	



    }
  };
  
  Sacrifice::UserHooksFactory::Creator<Pythia8::WW_to_Es_and_MUs> TTOE("WW_to_Es_and_MUs");
  
}

