// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
//#include "Rivet/Projections/MissingMomentum.hh"
//#include "Rivet/Projections/PartonicTops.hh"
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


      // Define Cuts
      Cut plJetLepCuts = ( Cuts::abseta < 2.5  ) & ( Cuts::pT > 25*GeV ); 


      // Initialise and Register Projections

      const FinalState particlefs;
      //      declare(PartonicTops(PartonicTops::E_MU), "LeptonicPartonTops");
      //      declare(FastJets(particlefs, FastJets::ANTIKT, 0.4), "particleAntiKt04Jets"); 

      //      FastJets jets(vfs, FastJets::ANTIKT, 0.4);


      // Find and dress particle level electrons and muons
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

      // Find neutrinos for jet clustering
      IdentifiedFinalState neutrinoidfs;
      neutrinoidfs.acceptNeutrinos();
      PromptFinalState neutrinopfs(neutrinoidfs);
      neutrinopfs.acceptTauDecays(true);

      // Jet clustering.
      VetoedFinalState vfs;
      //vfs.vetoNeutrinos();
      vfs.addVetoOnThisFinalState(dressedElectrons);
      vfs.addVetoOnThisFinalState(dressedMuons);
      vfs.addVetoOnThisFinalState(neutrinopfs);
      const FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      //jets.useInvisibles(true); // ??? What does this do?
      declareProjection(jets, "Jets");


      /*
      // MET (prompt neutrinos)
      VetoedFinalState ivfs;
      ivfs.addVetoOnThisFinalState(VisibleFinalState());
      declare(PromptFinalState(ivfs), "MET");

      // Jets
      VetoedFinalState jet_fs;
      jet_fs.vetoNeutrinos();
      jet_fs.addVetoPairId(PID::MUON);
      const FastJets fastjets(jet_fs, FastJets::ANTIKT, 0.4);
      declare(fastjets, "Jets");
      */
      // Book Histograms
      //inclusive = bookHisto1D("Inclusive");
      fiducial  = bookHisto1D("Fiducial");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt();     

      const double weight = event.weight();

      // Define success/fail flags
      bool pl_fail = false;
    
      vector<DressedLepton> dressedElectrons = applyProjection<DressedLeptons>( event, "dressedElectrons").dressedLeptons();
      vector<DressedLepton> dressedMuons     = applyProjection<DressedLeptons>( event, "dressedMuons"    ).dressedLeptons();


      // Particle level jets
      //      const Jets pl_jets = apply<FastJets>(event, "particleAntiKt04Jets").jetsByPt();
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);    
      //const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(plJetLepCuts);
      //      const Jets& jets = applyProjection<FastJets>(event, "particleAntiKt04Jets").jetsByPt();     
      //     particleAntiKt04Jets"); 
      //      unsigned int jetMultiplicity = jets.size();
      //if ( jetMultiplicity <= 2 ) { N1; pl_fail = true; }
      bool single_electron = (dressedElectrons.size() == 1) && (dressedMuons.empty());
      bool single_muon     = (dressedMuons.size() == 1) && (dressedElectrons.empty());

      DressedLepton* lepton = nullptr;
      if (single_electron)   lepton = &dressedElectrons[0];
      else if (single_muon)  lepton = &dressedMuons[0];
      
      if(!single_electron && !single_muon) { N1; RET; pl_fail = true; vetoEvent;}

      //      // (Particle level) Signal Region 2 requires a single isolated lepton
      //bool pl_singleLepton = ( leptonicPartonTops.size() == 1 );
      //if ( pl_singleLepton ) {

      //Particle pl_selectedLepton;
      //	pl_selectedLepton = leptonicPartonTops[0];


	//double pl_leptonPhi = leptonicPartonTops[0].phi();
	//double pl_leptonPT =  leptonicPartonTops[0].pT();

	/// reject event if a selected lepton is at a distance R < 0.4 of selected jet.
      
	for (const Jet i : jets) {
	  if (deltaR(i, *lepton) <= 0.4) { N2; RET; pl_fail = true; vetoEvent;}
	}	
	
	/*
	if (pl_fail) { N8; RET; vetoEvent;} 
	else {
	*/	  
	N9;
	RET;
	
	fiducial->fill(8.0, weight);      
  //}

    }


      /*
      cout << endl;////////////////////////// P R I N T //////////////////////////////////
      cout << " dl_muons.size() = " <<  dl_muons.size() << endl;
      cout << " dl_electrons.size() = " << dl_electrons.size() << endl;
      cout << " leptonicPartonTops.size() = " << leptonicPartonTops.size() << endl;
      */
      //bool good2 = ( dl_muons.size() + dl_electrons.size() == leptonicPartonTops.size() );
      //      if (!good2) WORRY;
      /*
      cout << "__________";
      cout << endl;///////////////////////////////////////////////////////////////////////
      */

  


    /// Normalise histograms etc., after the run
    void finalize() {

      //      normalize(_h_YYYY); // normalize to unity
      //scale( inclusive, crossSection()/picobarn/sumOfWeights());
      scale( fiducial,  crossSection()/picobarn/sumOfWeights()); // norm to cross section
      cout<<"sum of weights = " << sumOfWeights() << endl; 
      cout<<"cross section = " << crossSection() << endl;
    }


    //@}


    /// @name Histograms
    //@{

    Histo1DPtr fiducial;

 //@}


};


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1644099);


}
