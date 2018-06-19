#include "Sacrifice/UserHooksFactory.hh"
#include "Pythia8/UserHooks.h"
#include "Pythia8/PhaseSpace.h"

namespace Pythia8 {


  struct W_to_E_MU : public UserHooks {

    bool canVetoResonanceDecays() {cout << "HOOK: CAN_VETO_RESONANCES!"; return true; }

    bool doVetoResonanceDecays(Event& process) {
      //      cout << "HOOK: doVetoResonanceDecays(Event& process) init..." << endl;    


      int nT = 0, nH = 0, nE_MU_TAU = 0, nE_MU = 0;
      bool topFlag = false;


      for (int i = 1; i < process.size(); ++i) {
	//	cout << "_loop";
        Particle& p = process[i];
        if (p.idAbs() != 6) continue; // find tops 
	//	cout << "HOOK: topfound!";      
	topFlag = true;
	// cout << i << ": PID=" << p.id() << " status=" << p.status() << " d1=" << p.daughter1() << ", d2=" << p.daughter2() << endl;
        Particle* w = nullptr;
        if (p.daughter1() != 0 && process[p.daughter1()].idAbs() == 24) w = &process[p.daughter1()];
     	if (p.daughter2() != 0 && process[p.daughter2()].idAbs() == 24) w = &process[p.daughter2()];
     	if (w == nullptr) continue;
        Particle& wd1 = process[w->daughter1()];
        Particle& wd2 = process[w->daughter2()];

	//	       	cout << "HOOK: wd1.id()= " << wd1.id() << ", wd2.id()= " << wd2.id() << endl; 

	if ( wd1.idAbs() == 11 || wd1.idAbs() == 13 || wd2.idAbs() == 11 || wd2.idAbs() == 13) {nE_MU += 1;}// cout << "nE_MU+=1" << endl ;}
	
	else if ( wd1.isQuark()  && wd2.isQuark()  ) { 
	  nH += 1; 
	  //cout << "nH+=1" << endl ;
	}


      // USER INFO:
      /*
      if (!topFlag) {
	cout << "HOOK: no tops found!"<< endl;
	return 1;
      }
      */

      if (!topFlag) {
	cout << "No tops found?!" << endl;
	return false;
      }
      
      //cout << "HOOK: No. of tops found = " << nE_MU_TAU + nH << endl;
      //cout << "HOOK: nE_MU_TAU = " << nE_MU_TAU ;
      //cout << "; return (nE_MU_TAU == 0) = " << (nE_MU_TAU == 0) << endl ;
      /*
	if ( nE_MU_TAU == 0 ) {
	  //cout << "HOOK: VETO! retry..." << endl;
	  topFlag = false; // reset flag
	  return (nE_MU_TAU == 0) ; 
	} else {
	  //cout << "HOOK: Accept event..." << endl;
	  topFlag = false; // reset flag
	  return (nE_MU_TAU == 0) ;
	}
*/
	if ( nE_MU == 0 ) {
	  //	  cout << "HOOK: VETO! retry..." << endl;
	  topFlag = false; // reset flag
	  return (nE_MU == 0) ; 
	} else {
	  //cout << "HOOK: Accept event..." << endl;
	  topFlag = false; // reset flag
	  return (nE_MU == 0) ;
	}


      }
     

      /*  if ( nE_MU_TAU + nH != 2) { //cout << "single W in process -- return 'false'" << endl ;
	cout << "HOOK: No. of tops found = " << nE_MU_TAU + nH << endl;
	cout << "HOOK: nE_MU_TAU = " << nE_MU_TAU << endl;
	cout << "HOOK: return " << (nE_MU_TAU == 0) ;
	if ( nE_MU_TAU == 0 ) cout << " <- VETO! retry..." << endl;
	return nE_MU_TAU == 0 ; } // keeps single top events with leptons
      else {
	cout << "HOOK: No. of tops found = " << nE_MU_TAU + nH << endl;
	cout << "HOOK: nE_MU_TAU = " << nE_MU_TAU << endl;
	cout << "HOOK: return " << (nE_MU_TAU == 0) ;
	if ( nE_MU_TAU == 0 ) cout << " <- VETO! retry..." << endl;
      return (nE_MU_TAU == 0) ; // vetos fully hadronic ttbar events 
      */    
    }
    
  };
  
  Sacrifice::UserHooksFactory::Creator<Pythia8::W_to_E_MU> TTOE("W_to_E_MU");
  
}
