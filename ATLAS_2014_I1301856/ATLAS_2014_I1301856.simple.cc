// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Tools/ParticleUtils.hh"
#include "HepMC/IO_GenEvent.h"

// ============= Neil's Debug Toolkit ============== //
#define N1 do { std::cout << "1 "; } while(0)
#define N2 do { std::cout << "2 "; } while(0)
#define N3 do { std::cout << "3 "; } while(0)
#define N4 do { std::cout << "4 "; } while(0)
#define N5 do { std::cout << "5 "; } while(0)
#define N6 do { std::cout << "6 "; } while(0)
#define N7 do { std::cout << "7 "; } while(0)
#define N8 do { std::cout << "8 "; } while(0)
#define N9 do { std::cout << "9 "; } while(0)
#define WORRY   do { std::cout << "WORRY!"   << endl; } while(0)
#define SUCCESS do { std::cout << "SUCCESS!" << endl; } while(0)
#define BEGIN   do { std::cout << "Begin..." << endl; } while(0)
#define END     do { std::cout << "End..."   << endl; } while(0)
#define RET     do { std::cout << endl; } while(0)
// =================================================== //

namespace Rivet {


  /// @brief Parton-level top-pair cross-sections at 7TeV and 8TeV using dileptonic ($e \mu$) decays
  class ATLAS_2014_I1301856 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1301856);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      //if (fuzzyEquals(sqrtS(), 7*TeV)) {                                                          
      //} else {                                                                                    
      //}                                                                                           


      
      // Parton Level Top Quarks
      declare(PartonicTops(PartonicTops::ELECTRON), "ElectronPartonTops");
      declare(PartonicTops(PartonicTops::MUON), "MuonPartonTops");
      declare(PartonicTops(PartonicTops::ALL), "AllPartonTops");

      // Book counter
      _c_emuCtr = bookCounter("ttbarXSec");

     
      //_hepmcout.reset(new HepMC::IO_GenEvent("inspectevents.hepmc"));    
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      //BEGIN;

 
      const Particles partonicTops_el = apply<ParticleFinder>(event, "ElectronPartonTops").particlesByPt();
      const Particles partonicTops_mu = apply<ParticleFinder>(event, "MuonPartonTops").particlesByPt();
      const Particles partonicTops_all = apply<ParticleFinder>(event, "AllPartonTops").particlesByPt();

      // veto event if it is not ttbar event
      if (partonicTops_all.size() != 2) {N1; RET; vetoEvent;} //< debug mode!
      //if (partonicTops_all.size() != 2) vetoEvent

      // veto if not a opposite charged e mu pair
      if (partonicTops_el.size() == 1 && partonicTops_mu.size() == 1){
	if (partonicTops_el[0].charge() + partonicTops_mu[0].charge() == 0){
	  _c_emuCtr->fill(event.weight());
	  N9;RET;
	} else {N3;RET;} //< debug mode!  <<  not opposite charge 
      } else {N2;RET;} //< debug mode!    <<  not e mu pair


      
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double BR = 0.032; //< branching ratio: ttbar -> e mu nu nu b b

      double SF = crossSection()/picobarn/sumOfWeights()/BR;  //< scale factor	
      cout << "SF=" << SF << endl;
     
      scale (_c_emuCtr,SF);


      
      //_hepmcout->clear(); _hepmcout.reset();

    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_emuCtr;
    //@}

    std::unique_ptr<HepMC::IO_GenEvent> _hepmcout;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1301856);


}
