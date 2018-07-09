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
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "HepMC/IO_GenEvent.h"


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
      DressedLeptons dressedElectrons(photons, promptElectrons, 0.1, e_muCuts, true);
      //DressedLepton dressedElectrons(promptElectrons, 0.1, e_muCuts, true);
      declare(dressedElectrons, "DressedElectrons");
      
      IdentifiedFinalState muon_fs;
      muon_fs.acceptIdPair(PID::MUON);
      PromptFinalState promptMuons(muon_fs, true, false); //< "true, false" accepts tau decays but not muon decays
      DressedLeptons dressedMuons(photons, promptMuons, 0.1, e_muCuts, true); //< last arg includes non-prompt photons
      declare(dressedMuons, "DressedMuons");
      
      WFinder w_electron(electron_fs, e_muCuts,PID::ELECTRON,35*GeV, 8000*GeV, 30*GeV, 0.1,  WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERNODECAY,WFinder::TRACK,  WFinder::TRANSMASS);
      // notes: ET_miss>30GeV (see: https://arxiv.org/pdf/1502.05923.pdf)
      // notes: mass window > 35GeV (< transverse mass) (see: https://arxiv.org/pdf/1502.05923.pdf)
      declare(w_electron, "W_Electron");
      
      WFinder w_muon(muon_fs,  e_muCuts, PID::MUON, 35*GeV, 100*GeV, 30*GeV, 0.1, WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERNODECAY, WFinder::TRACK, WFinder::TRANSMASS);
      declare(w_muon, "W_Muon");

      //the selected lepton (with it's photons) is not to be included in jet clustering
      VetoedFinalState jet_input(fs);
      //declare(jet_input, "JetInput");
      jet_input.addVetoOnThisFinalState(dressedElectrons);
      jet_input.addVetoOnThisFinalState(dressedMuons);
      declare(FastJets(jet_input, FastJets::ANTIKT, 0.4), "Jets");

      // do parton level analysis seperatelty!!
      //declare(PartonicTops(PartonicTops::E_MU), "LeptonicTops") ;
      //declare(PartonicTops(PartonicTops::ALL ), "AllPartonicTops") ;
      
      // book counters
      _c_fid_t = bookCounter("fiducialXSec_top");
      _c_fid_tbar = bookCounter("fiducialXSec_anti-top");
      
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



      // Find leptons, bosons and jets
      const Particles& electrons = apply<DressedLeptons>(event, "DressedElectrons").particlesByPt();
      const Particles& muons = apply<DressedLeptons>(event, "DressedMuons").particlesByPt();
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(); //< keep non-const for later filtering
      const WFinder& wFinder_el = apply<WFinder>(event, "W_Electron");
      const WFinder& wFinder_mu = apply<WFinder>(event, "W_Muon");
      const Particles WEl= wFinder_el.bosons();
      const Particles WMu= wFinder_mu.bosons();
      //const Particles leptonicTops = apply<ParticleFinder>(event, "LeptonicTops").particlesByPt();
      //const Particles allPartonicTops = apply<ParticleFinder>(event, "AllPartonicTops").particlesByPt();

      /// Particle-level analysis
      
      //      bool particleVeto = false;

      // find jets
      Cut jetCuts = (Cuts::abseta < 4.5 && Cuts::pT > 30*GeV);
      ifilter_select(jets, jetCuts);
    

      // find and count b-jets
      Jet bJet;
      bool bJetFound = false;
      for (const Jet& j : jets){
        if ( j.bTagged(Cuts::pT > 5*GeV) && j.abseta() < 2.5) {
	  if (bJetFound) vetoEvent;
	  else bJetFound = true;
	}
      }


      if (jets.size() != 2 || !bJetFound) vetoEvent;

      
      // Find the leptons
      Particles leptons;
      for (const Particle& e : electrons){
        leptons.push_back(e);
      }
      for (const Particle& m : muons){
        leptons.push_back(m);
      }
      
      
      // discared event if selected jets and single lepton are not isolated
      bool tooClose= false;
      if (leptons.size()==1){
        for (const Jet& j : jets){
	  if (deltaR(j, leptons[0]) < 0.4) tooClose = true ;
	}
	if (tooClose) vetoEvent;
      } else vetoEvent;
      
    
      // Then this implies there must be 2 jets of which only one is
      // b-tagged and there is exactly one isolated particle level
      // electron or muon in the event.
      //*************************************************************************
      
      
      // make cut on the mass of lepton + b-jet system
      FourMomentum lepton4p = leptons[0];
      FourMomentum bJet4p = bJet;
      FourMomentum lb4p = add(lepton4p, bJet4p);
      if ( lb4p.mass() > 160*GeV ) vetoEvent;


      
      FourMomentum w4p;
      if ( WEl.size() + WMu.size() == 1){
	if   (WEl.size() == 1 ) w4p = WEl[0] ;
	else if ( WMu.size() == 1 ) w4p = WMu[0] ;

      } else vetoEvent; //< means w-boson was not reconstructed
	

      // define pseudoTop
      const FourMomentum pseudoTop4p = add(bJet4p, w4p);


      if (leptons[0].charge() > 0){ // fill top quark diff xsecs
	_c_fid_t->fill(event.weight()) ;
	
        _h_AbsPtclDiffXsecTPt->fill(pseudoTop4p.pT(),     event.weight()) ;

	cout << "pseudoTop4p.pT()=" << pseudoTop4p.pT() << endl;
	cout << "pseudoTop4p.absrap()=" << pseudoTop4p.absrap() << endl;
	cout << "weight=" << event.weight() << endl;
        _h_AbsPtclDiffXsecTY->fill( pseudoTop4p.absrap(), event.weight()) ;


	
        _h_NrmPtclDiffXsecTPt->fill(pseudoTop4p.pT(),     event.weight()) ;
        _h_NrmPtclDiffXsecTY->fill( pseudoTop4p.absrap(), event.weight()) ;
	

      } else { // fill top anti-quark diff xsecs
        _c_fid_tbar->fill(event.weight()) ;
	
        _h_AbsPtclDiffXsecTbarPt->fill(pseudoTop4p.pT(),     event.weight()) ;
        _h_AbsPtclDiffXsecTbarY->fill( pseudoTop4p.absrap(), event.weight()) ;
	
        _h_NrmPtclDiffXsecTbarPt->fill(pseudoTop4p.pT(),     event.weight()) ;
        _h_NrmPtclDiffXsecTbarY->fill( pseudoTop4p.absrap(), event.weight()) ;
      }
    }

     

   


    /// Normalise histograms etc., after the run
    void finalize() {
      double SF = crossSection()/picobarn/sumOfWeights();  // scale factor

      scale(_c_fid_t, SF);
      scale(_c_fid_tbar, SF);
      
      scale(_h_AbsPtclDiffXsecTPt, SF);
      scale(_h_AbsPtclDiffXsecTY, SF);
      scale(_h_AbsPtclDiffXsecTbarPt, SF);
      scale(_h_AbsPtclDiffXsecTbarY, SF);
      
      

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
