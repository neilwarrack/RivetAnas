// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PartonicTops.hh"

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
      Cut mu_detector_level_eta =   Cuts::abseta < 2.5 ;
      Cut el_detector_level_eta = ( Cuts::abseta < 1.37 ) & ( Cuts::absetaIn(1.52, 2.47) );
      Cut particle_level        = ( Cuts::abseta < 2.5  ) & ( Cuts::pT > 25*GeV ); 


      // Initialise and Register Projections
      // General
      const FinalState detectorfs;
      const FinalState particlefs(particle_level);
      addProjection(MissingMomentum(detectorfs), "MissingP") ;
      declare(PartonicTops(PartonicTops::E_MU), "LeptonicPartonTops");
      declare(FastJets(detectorfs, FastJets::ANTIKT, 0.4), "detectorAntiKt04Jets"); 
      declare(FastJets(particlefs, FastJets::ANTIKT, 0.4), "particleAntiKt04Jets"); 

      // Detector level
      IdentifiedFinalState detelectronfs( el_detector_level_eta );
      detelectronfs.acceptIdPair(PID::ELECTRON);
      declare(detelectronfs, "detectorElectrons");

      IdentifiedFinalState detmuonfs( mu_detector_level_eta );
      detmuonfs.acceptIdPair(PID::MUON);
      declare(detmuonfs, "detectorMuons"); 

      // Particle level
      IdentifiedFinalState partelectronfs(particlefs);
      partelectronfs.acceptIdPair(PID::ELECTRON);
      declare(partelectronfs, "particleElectrons");

      IdentifiedFinalState partmuonfs(particlefs);
      partmuonfs.acceptIdPair(PID::MUON);
      declare(partmuonfs, "particleMuons");

      IdentifiedFinalState photonfs(particlefs);   // for dressing particle-level leptons
      photonfs.acceptIdPair(PID::PHOTON);
      DressedLeptons dressedElectronfs(photonfs, partelectronfs, 0.1, particle_level, true, true);
      declare(dressedElectronfs, "particleDressedElectrons");
      DressedLeptons dressedMuonfs(photonfs, partmuonfs, 0.1, particle_level, true, true);
      declare(dressedMuonfs, "particleDressedMuons");


      // Book Histograms
      inclusive = bookHisto1D("Inclusive");
      fiducial  = bookHisto1D("Fiducial");
    
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      //BEGIN;

      // Define success/fail flags
      bool pl_fail = false;
      bool dl_fail = false;
      // Get objects from projections
      const Particles leptonicPartonTops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();

      const Particles dl_electrons = apply<ParticleFinder>(event, "detectorElectrons").particlesByPt();
      const Particles dl_muons     = apply<ParticleFinder>(event, "detectorMuons"    ).particlesByPt();

      pl_dressedElectrons = apply<DressedLeptons>(  event, "particleDressedElectrons").dressedLeptons();
      pl_dressedMuons     = apply<DressedLeptons>(  event, "particleDressedMuons"    ).dressedLeptons();


      // Particle level jets
      const Jets pl_jets = apply<FastJets>(event, "particleAntiKt04Jets").jetsByPt();
      unsigned int pl_jetMultiplicity = pl_jets.size();
      if ( pl_jetMultiplicity != 4 ) { N1; pl_fail = true; }
            unsigned int pl_bCtr = 0 ; 
      for (const Jet j : pl_jets) { 
	if ( j.bTagged() ) pl_bCtr++;
      }     
      if ( pl_bCtr != 2 ) { N2; pl_fail = true; }

      // Detector level jets
      const Jets dl_jets = apply<FastJets>(event, "detectorAntiKt04Jets").jetsByPt();
      unsigned int dl_jetMultiplicity = dl_jets.size();
      if ( dl_jetMultiplicity != 4 ) { N3; dl_fail = true;
}
      unsigned int dl_bCtr = 0 ; 
      for (const Jet k : dl_jets) { 
	if ( k.bTagged() ) dl_bCtr++;
      }     
      if ( dl_bCtr != 2 ) { N4; dl_fail = true; }
      
      /*
      // (Particle level) Signal Region 2 requires a single isolated lepton
      bool pl_singleLepton = ( pl_dressedMuons.size() + pl_dressedElectrons.size() == 1 );
      if ( pl_singleLepton ) {

	Particle pl_selectedLepton;
	if ( pl_dressedMuons.size() == 1 ) {
	  pl_selectedLepton = pl_dressedMuons[0];
	} else {
	  pl_selectedLepton = pl_dressedElectrons[0];
	}

	//double pl_leptonPhi = leptonicPartonTops[0].phi();
	//double pl_leptonPT =  leptonicPartonTops[0].pT();

	/// reject event if a selected lepton is at a distance R < 0.4 of selected jet.
	for (const Jet i : pl_jets) {
	  if (deltaR(i, pl_selectedLepton) <= 0.4) {N6; pl_fail = true; }
	}
	
	//SUCCESS;
	
      } else { 
	N5;
	pl_fail = true;
	//	vetoEvent; 
      }
*/


      // (Particle level) Signal Region 2 requires a single isolated lepton
      bool pl_singleLepton = ( leptonicPartonTops.size() == 1 );
      if ( pl_singleLepton ) {

	Particle pl_selectedLepton;
	pl_selectedLepton = leptonicPartonTops[0];


	//double pl_leptonPhi = leptonicPartonTops[0].phi();
	//double pl_leptonPT =  leptonicPartonTops[0].pT();

	/// reject event if a selected lepton is at a distance R < 0.4 of selected jet.
	for (const Jet i : pl_jets) {
	  if (deltaR(i, pl_selectedLepton) <= 0.4) {N6; pl_fail = true; }
	}
	
	//SUCCESS;
	
      } else { 
	N5;
	pl_fail = true;
	//	vetoEvent; 
      }






      // (Detector level) Signal Region 2 requires a single isolated lepton
      bool dl_singleLepton = ( dl_muons.size() + dl_electrons.size() == 1 );
      if ( dl_singleLepton ) {
	
	Particle dl_selectedLepton;
	if ( dl_muons.size() == 1 ) {
	  dl_selectedLepton = dl_muons[0];
	} else {
	  dl_selectedLepton = dl_electrons[0];
	}

	//      double dl_leptonPhi = leptonicPartonTops[0].phi();
	//      double dl_leptonPT =  leptonicPartonTops[0].pT();

      /// reject event if a selected lepton is at a distance R < 0.4 of selected jet.
      for (const Jet l : dl_jets) {
	if (deltaR(l, dl_selectedLepton) <= 0.4) { N7; dl_fail = true; }
      }
      
      //      SUCCESS;
      } else {// WORRY; 
	N8; //<- prints "1 " to stdout 
	dl_fail = true;
	//vetoEvent; 
      }

      /*
      cout << endl;////////////////////////// P R I N T //////////////////////////////////
      cout << " dl_muons.size() = " <<  dl_muons.size() << endl;
      cout << " dl_electrons.size() = " << dl_electrons.size() << endl;
      cout << " leptonicPartonTops.size() = " << leptonicPartonTops.size() << endl;
      */
      bool good2 = ( dl_muons.size() + dl_electrons.size() == leptonicPartonTops.size() );
      //      if (!good2) WORRY;
      /*
      cout << "__________";
      cout << endl;///////////////////////////////////////////////////////////////////////
      */

      // --------------------------------------------------------------------------------------------------------
      // --------------------------------------------------------------------------------------------------------
      // --------------------------------------------------------------------------------------------------------


      if (!pl_fail && !dl_fail) SUCCESS;
      if (!pl_fail) N9; 
      //      if (!dl_fail) {N9; SUCCESS;}      
      RET;
      vetoEvent; 
      // MET
      const MissingMomentum missP = applyProjection<MissingMomentum>(event, "MissingP");
      const FourMomentum missingFourMomentum = -missP.visibleMomentum();
      const double met = missingFourMomentum.pT();
      //if ( met <= 25*GeV ) RET; vetoEvent;

      // W transverse mass cut
      double metPhi = missingFourMomentum.phi() ;
      double wmt = sqrt( 2*leptonicPartonTops[0].pT() * met * (1-cos(leptonicPartonTops[0].phi() - metPhi))) ;
      const bool wMassBelowThreashold = ( wmt < 30*GeV ) ;
      //if ( wMassBelowThreashold ) RET; vetoEvent ;
     

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //      normalize(_h_YYYY); // normalize to unity
      scale( inclusive, crossSection()/picobarn/sumOfWeights());
      scale( fiducial,  crossSection()/picobarn/sumOfWeights()); // norm to cross section
 
    }


    //@}


    /// @name Histograms
    //@{

    Histo1DPtr inclusive, fiducial;
    vector<DressedLepton> pl_dressedElectrons, pl_dressedMuons;
    Jets pl_jets;   
 //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1644099);


}
