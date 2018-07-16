// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
//#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

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


  /// @brief The inclusive and fiducial ttbar production cross-sections are measured in the lepton+jets channel at 8TeV
  class ATLAS_2017_I1644099 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1644099);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      const FinalState particlefs;
	    
      // Define Cuts
      Cut lepCuts = ( Cuts::abseta < 2.5  ) & ( Cuts::pT > 25*GeV );
      Cut plJetLepCuts = ( Cuts::abseta < 2.5  ) & ( Cuts::pT > 25*GeV ); 
      

      // Initialise and Register Projections


      //      declare(PartonicTops(PartonicTops::E_MU), "LeptonicPartonTops");
      //      declare(FastJets(particlefs, FastJets::ANTIKT, 0.4), "particleAntiKt04Jets"); 

      //      FastJets jets(vfs, FastJets::ANTIKT, 0.4);


      // Projetions for (particle level) charged leptons
      // (dressed prompt electrons and muons, including form tau decays)
      
      IdentifiedFinalState lepfs(fs, {PID::ELECTRON, -PID::ELECTRON, PID::MUON, -PID::MUON});
      PromptFinalState promptLeptons(lepfs,true); // bool accepts tau decays
      IdentifiedFinalState photons(fs, {PID::PHOTON, -PID::PHOTON});
      DressedLeptons leptons(photons, promptLeptons, 0.1, lepCuts);
      declare(leptons, "Leptons");

      // Find neutrinos for jet clustering
      IdentifiedFinalState neutrinos;
      neutrinos.acceptNeutrinos();
      PromptFinalState promptNeutrinos(neutrinos);
      promptNeutrinos.acceptTauDecays();

      // Project jets
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(promptNeutrinos);
      vfs.addVetoOnThisFinalState(leptons);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      declareProjection(jets, "Jets");


      // Auxiliry projection for ttbar filter
      declare(PartonicTops(PartonicTops::ALL), "PartonicTops");
      


      ////////////////////////
      /*
      IdentifiedFinalState partelectronidfs(particlefs);
      partelectronidfs.acceptIdPair(PID::ELECTRON);
      declareProjection(partelectronidfs, "electrons");
      
      IdentifiedFinalState partmuonidfs(particlefs);
      partmuonidfs.acceptIdPair(PID::MUON);
      declareProjection(partmuonidfs, "muons");

      IdentifiedFinalState photonidfs(particlefs);
      photonidfs.acceptIdPair(PID::PHOTON);

      DressedLeptons dressedElectrons(photonidfs, partelectronidfs, 0.1, plJetLepCuts, true);
      declareProjection(dressedElectrons, "dressedElectrons");
      DressedLeptons dressedMuons(photonidfs, partmuonidfs, 0.1, plJetLepCuts, true);
      declareProjection(dressedMuons, "dressedMuons");
      */

      // Book Histograms
      _c_inclusive = bookCounter("Inclusive");
      _c_fiducial  = bookCounter("Fiducial");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      
      const double weight = event.weight();

      
      // Quick veto on non ttbar events
      if (applyProjection<PartonicTops>(event, "PartonicTops").particles().size() < 2 ) vetoEvent; 
      _c_inclusive->fill(weight);
      
      // Find particle level objects
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      const Particles leptons = applyProjection<DressedLeptons>(event, "Leptons").particlesByPt();

      // Veto any events without single prompt charged lepton
      if (leptons.size() != 1) vetoEvent;
      
      // jet multiplicity veto
      if (jets.size() < 3) vetoEvent;

      // jet lepton overlap veto
      for (const Jet& j : jets){
	if (deltaR(j,leptons[0]) < 0.4) vetoEvent; 
	}


    
    // fill fiducial counter
    _c_fiducial->fill(weight);
     
    }


    
    /// Normalise histograms etc., after the run
    void finalize() {

      double scaleFactor = crossSection()/picobarn/sumOfWeights();

      scale( _c_inclusive, scaleFactor);
      scale( _c_fiducial,  scaleFactor);
     
    }


    //@}
private:

    CounterPtr _c_fiducial, _c_inclusive;
};


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1644099);
}
