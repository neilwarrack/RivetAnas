// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/JetUtils.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/WFinder.hh"
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


  /// @brief Total and fiducial cross-section measurements with normalised differential cross-sections and absolute differential cross-sections for single top and anti-top production in pp collisions at 8*TeV.
  
  class ATLAS_2017_I1512776 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1512776);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;

      // Particle-level electron and muon cuts
      Cut e_muCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV ;
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      declare(photons, "Photons");
      
      IdentifiedFinalState electron_fs;
      electron_fs.acceptIdPair(PID::ELECTRON);
      PromptFinalState promptElectrons(electron_fs, true, false);
      DressedLeptons dressedElectrons(photons, promptElectrons, 0.1, e_muCuts, false);
      declare(dressedElectrons, "DressedElectrons");
      
      IdentifiedFinalState muon_fs;
      muon_fs.acceptIdPair(PID::MUON);
      PromptFinalState promptMuons(muon_fs, true, false);
      DressedLeptons dressedMuons(photons, promptMuons, 0.1, e_muCuts, false);
      declare(dressedMuons, "DressedMuons");
      
      WFinder w_electron(electron_fs, e_muCuts,PID::ELECTRON,35*GeV, 8000*GeV, 30*GeV, 0.1,  WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERNODECAY,WFinder::TRACK,  WFinder::TRANSMASS);
      declare(w_electron, "W_Electron");
      
      WFinder w_muon(muon_fs,  e_muCuts, PID::MUON, 35*GeV, 100*GeV, 30*GeV, 0.1, WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERNODECAY, WFinder::TRACK, WFinder::TRANSMASS);
      declare(w_muon, "W_Muon");
      
      VetoedFinalState jet_input(fs);
      jet_input.addVetoOnThisFinalState(dressedElectrons);
      jet_input.addVetoOnThisFinalState(dressedMuons);
      declare(FastJets(jet_input, FastJets::ANTIKT, 0.4), "Jets");
      
      declare(PartonicTops(PartonicTops::E_MU), "LeptonicTops") ;
      declare(PartonicTops(PartonicTops::ALL ), "AllPartonicTops") ;
      
      
      // book histograms
      // top quark differential cross sections                      
      _h_AbsPtclDiffXsecTPt   = bookHisto1D(1,1,1);
      _h_AbsPtclDiffXsecTY    = bookHisto1D(1,1,2);
      _h_NrmPtclDiffXsecTPt   = bookHisto1D(2,1,1);
      _h_NrmPtclDiffXsecTY    = bookHisto1D(2,1,2);
      
      // top antiquark differential cross sections                  
      _h_AbsPtclDiffXsecTbarPt   = bookHisto1D(1,1,3);
      _h_AbsPtclDiffXsecTbarY    = bookHisto1D(1,1,4);
      _h_NrmPtclDiffXsecTbarPt   = bookHisto1D(2,1,3);
      _h_NrmPtclDiffXsecTbarY    = bookHisto1D(2,1,4);


      
      // debug hepmc print-out
      _hepmcout.reset(new HepMC::IO_GenEvent("inspectevents.hepmc"));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //BEGIN;

      // Find leptons, bosons and jets
      const Particles& electrons = apply<DressedLeptons>(event, "DressedElectrons").particlesByPt();
      const Particles& muons =     apply<DressedLeptons>(event, "DressedMuons").particlesByPt();

      Jets jets = apply<FastJets>(event, "Jets").jetsByPt();
      const WFinder& w_el = apply<WFinder>(event, "W_Electron");
      const WFinder& w_mu = apply<WFinder>(event, "W_Muon");
      const Particles w_elP= w_el.bosons();
      const Particles w_muP= w_mu.bosons();
      const Particles leptonicTops = apply<ParticleFinder>(event, "LeptonicTops").particlesByPt();
      const Particles allPartonicTops = apply<ParticleFinder>(event, "AllPartonicTops").particlesByPt();

      /// Particle-level analysis
      
      bool particleVeto = false;

      // find jets
      Cut jetCuts = Cuts::abseta < 4.5 && Cuts::pT > 30*GeV ;
      cout << "sizeOfJets1=" << jets.size()<< " "; RET;
      ifilter_select(jets, jetCuts);
      cout << "sizeOfJets2=" << jets.size()<< " "; RET;

      // Find the leptons
      cout << "muons.size()=" << muons.size(); RET;
      cout << "electrons.size()=" << electrons.size(); RET;
      Particles leptons;
      for (const Particle& e : electrons){
        leptons.push_back(e);
      }
      for (const Particle& m : muons){
        leptons.push_back(m);
      }
      cout << "leptons.size()=" << leptons.size(); RET;

      // isolate jets and leptons
      bool tooClose= false;;

      if (!leptons.empty() && !jets.empty()){
        for (const Jet& j : jets){
          for (const Particle& l : leptons){
          }
          if (tooClose) particleVeto = true;
          tooClose = false;
        }
      }


      // find and count b-jets
      Jet bJet;
      int b=0;

      for (const Jet& j : jets){
        if ( j.bTagged(Cuts::pT > 5*GeV) && j.abseta() < 2.5) {
          b++ ;
          bJet = j;
        }
      }

      // veto at particle level
      if(leptons.size() != 1) {N2;RET; vetoEvent;}
      else if(b!=1) {N3;RET;vetoEvent;}
      else if(jets.size() != 2) {N4;RET;vetoEvent;}
      else if(particleVeto) {N5;RET;vetoEvent;} // lepton too close to a jet


      if ( b==1 && leptons.size() == 1 && jets.size() == 2 && !particleVeto ){

	// Then this implies there must be 2 jets of which only one is
        // b-tagged and there is exactly one isolate particle level
        // electron or muon in the event.
        // make cut on the mass of lepton + b-jet system
        FourMomentum lepton4p = leptons[0];
        FourMomentum bJet4p = bJet;
        FourMomentum lb4p = lepton4p + bJet4p;

	if ( lb4p.mass() < 160*GeV ) particleVeto = true;
        if(particleVeto) {N6;RET;vetoEvent;}
  
      FourMomentum w4p;
      // cout << "W-el=" << w_elP.size() << ", " ;
      // cout << "W-mu=" << w_muP.size() << " " ;
      // cout << w_el.constituentLepton().pid() << " ";
      // cout << w_mu.constituentLepton().pid() ;
      // if ( w_elP.size() != w_muP.size() ) WORRY;

      if ( w_elP.size() + w_muP.size() == 1){
	if   ( w_elP.size() == 1 )    w4p = w_elP[0] ;
	else if ( w_muP.size() == 1 ) w4p = w_muP[0] ;
	else {particleVeto = true;}
	
	if (!particleVeto){
	  
	  // define pseudoTop
	  const FourMomentum pseudoTop4p = bJet4p + w4p;
	  
	  if (leptons[0].charge() > 0){
	    N7;
	    _c_fid_t->fill(event.weight()) ;
	    
	    _h_AbsPtclDiffXsecTPt->fill(pseudoTop4p.pT(),     event.weight()) ;
	    _h_AbsPtclDiffXsecTPt->fill(pseudoTop4p.absrap(), event.weight()) ;
	    
	    _h_NrmPtclDiffXsecTPt->fill(pseudoTop4p.pT(),     event.weight()) ;
	    _h_NrmPtclDiffXsecTPt->fill(pseudoTop4p.absrap(), event.weight()) ;
	    
	    
	  } else {
	    N8;
	    _c_fid_tbar->fill(event.weight()) ;
	    
	    _h_AbsPtclDiffXsecTbarPt->fill(pseudoTop4p.pT(),     event.weight()) ;
	    _h_AbsPtclDiffXsecTbarPt->fill(pseudoTop4p.absrap(), event.weight()) ;
	    
	    _h_NrmPtclDiffXsecTbarPt->fill(pseudoTop4p.pT(),     event.weight()) ;
	    _h_NrmPtclDiffXsecTbarPt->fill(pseudoTop4p.absrap(), event.weight()) ;
	    
	  }
	}
      }
      else if ( w_elP.size() + w_muP.size() > 1) {N1;}
      }//< end of particle-level analysis
    cout << "particleVeto=" << particleVeto;RET;
    

    /// Parton-level analysis
    bool partonVeto = false;

    // find partonic tops                                                                              

    if ( !allPartonicTops.empty() ){
      cout << "-Number of Tops:" << allPartonicTops.size(); RET;
      if ( allPartonicTops.size() != 1 ) {
	partonVeto = true ;
      } else {
	
	
	if (allPartonicTops[0].pid() > 0){
	  
	  cout << "q-allPartonicTops[0].pid()=" << allPartonicTops[0].pid();RET;
	  
	  //      _h_AbsPtonDiffXsecTPt->fill(leptonicTops[0].pT(), event.weight());
	}
	if (allPartonicTops[0].pid() < 0){

	  cout << "qbar-allPartonicTops[0].pid()=" << allPartonicTops[0].pid();
	  
	}
      }
    }
    
    RET;


    } //< end of analyze()            

   


    /// Normalise histograms etc., after the run
    void finalize() {
      double SF = crossSection()/femtobarn/sumOfWeights();  // scale factor
      scale(_c_fid_t, SF);
      scale(_c_fid_tbar, SF);
      
      scale(_h_AbsPtclDiffXsecTPt, SF);
      scale(_h_AbsPtclDiffXsecTY, SF);
      scale(_h_AbsPtclDiffXsecTbarPt, SF);
      scale(_h_AbsPtclDiffXsecTbarY, SF);
      scale(_h_AbsPtonDiffXsecTPt, SF);
      

    }

    //@}


    /// @name Histograms
    //@{

  CounterPtr _c_fid_t, _c_fid_tbar;
  Histo1DPtr _h_AbsPtclDiffXsecTPt, _h_AbsPtclDiffXsecTY,
    _h_NrmPtclDiffXsecTPt, _h_NrmPtclDiffXsecTY,
    _h_AbsPtclDiffXsecTbarPt, _h_AbsPtclDiffXsecTbarY,
    _h_NrmPtclDiffXsecTbarPt,  _h_NrmPtclDiffXsecTbarY,
    _h_AbsPtonDiffXsecTPt;
  
    //@}

  std::unique_ptr<HepMC::IO_GenEvent> _hepmcout;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1512776);


}
