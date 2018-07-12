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
      _c_fid_t =    bookCounter("fiducialXSec_top");
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

      BEGIN;

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
      cout << "sizeOfJetsPreFilter=" << jets.size()<< " "; RET;
      ifilter_select(jets, jetCuts);
      cout << "sizeOfJetsPostFilter=" << jets.size()<< " "; RET;


      // find and count b-jets
      Jet bJet;
      bool bJetFound = false;
      for (const Jet& j : jets){
        if ( j.bTagged(Cuts::pT > 5*GeV) && j.abseta() < 2.5) {
	  if (bJetFound){N2;RET;vetoEvent;}//< debug mode!
	  //if (bJetFound){vetoEvent;}
	  else {
	    bJetFound = true;
	    bJet = j;
	  }
	}
      }

      if (jets.size() != 2){N3;RET;vetoEvent;} //< debug mode!
      if (!bJetFound)      {N4;RET;vetoEvent;} //< debug mode!
      //if (jets.size() != 2 || !bJetFound) vetoEvent;

      
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

      
      // discared event if selected jets and single lepton are not isolated
      bool tooClose= false;
      if (leptons.size()==1){
        for (const Jet& j : jets){
	  if (deltaR(j, leptons[0]) < 0.4) tooClose = true ;
	}
	if (tooClose){N5;RET;vetoEvent;}//< debug mode!
	//if (tooClose) vetoEvent;
      } else //vetoEvent;
	{N6;RET;vetoEvent;}//< debug mode!
      
    
      // Then this implies there must be 2 jets of which only one is
      // b-tagged and there is exactly one isolated particle level
      // electron or muon in the event.
      //*************************************************************************
      

      // make cut on the mass of lepton + b-jet system
      FourMomentum lepton4p = leptons[0];
      FourMomentum bJet4p = bJet;
      cout << "bJet4p.pT=" << bJet4p.pT() << endl;
      cout << "bJet4p.absrap=" << bJet4p.absrap() << endl;
      cout << "lepton4p.pT=" << lepton4p.pT() << endl;
      cout << "lepton4p.absrap=" << lepton4p.absrap() << endl;
FourMomentum lb4p = add(lepton4p, bJet4p);
      if ( lb4p.mass() > 160*GeV )//vetoEvent;
	{N7;RET;vetoEvent;}//< debug mode!

      
      FourMomentum w4p;
      //cout << "WEl.size()=" << WEl.size() << ", " ;
      //cout << "WMu.size()=" << WMu.size() << " " ;
      //cout << "WEl.constituentLepton().pid()=" << WEl.constituentLepton().pid();RET;
      //cout << "WMu.constituentLepton().pid()=" << WMu.constituentLepton().pid();RET;
      //cout << w_mu.constituentLepton().pid() ;
      //if ( WEl.size() + WMu.size() != 1 ) WORRY;

      if ( WEl.size() + WMu.size() == 1){
	if (WEl.size() == 1 ) {
	  cout << "WEl[0].pT=" << WEl[0].pt() << endl;
	  cout << "WEl[0].absrap=" << WEl[0].absrap() << endl;
	  cout << "WEl[0].E=" << WEl[0].E() << endl;
	  cout << "WEl[0].px=" << WEl[0].px() << endl;
	  cout << "WEl[0].py=" << WEl[0].py() << endl;
	  cout << "WEl[0].pz=" << WEl[0].pz() << endl;
	  w4p = WEl[0];
	  cout <<"W reconstructed from electron!" << endl;
	  cout << "w4p.pT=" << w4p.pT() << endl;
	  cout << "w4p.absrap=" << w4p.absrap() << endl;
	}
	else if ( WMu.size() == 1 ) {
	  cout << "WMu[0].pT=" << WMu[0].pt() << endl;
	  cout << "WMu[0].absrap=" << WMu[0].absrap() << endl;

	  w4p = WMu[0] ;
	  cout <<"W reconstructed from muon!" << endl;
	  cout << "w4p.pT=" << w4p.pT() << endl;
	  cout << "w4p.absrap=" << w4p.absrap() << endl;	  
	}
	  else {WORRY;
	  cout << "WEl.size()=" << WEl.size() << " " <<  "WMu.size()=" << WMu.size() << " " << endl; 
	}
      } else {N8;RET;vetoEvent;} //< means w-boson was not reconstructed
	

      // define pseudoTop
      const FourMomentum pseudoTop4p = add(bJet4p, w4p);
      cout << "pseudoTop4p.pT()=" << pseudoTop4p.pT() << endl;
      cout << "pseudoTop4p.absrap()=" << pseudoTop4p.absrap() << endl;
      //      cout << "weight=" << event.weight() << endl;
      

      bool flag = false;
      if (leptons[0].charge() > 0){ // fill top quark diff xsecs
	N9; RET;//< debug mode!
	flag = true;
	cout << "top!"<<endl;
	_c_fid_t->fill(event.weight()) ;

	_h_AbsPtclDiffXsecTPt->fill(pseudoTop4p.pT()/GeV, event.weight()) ;
	_h_AbsPtclDiffXsecTY->fill( pseudoTop4p.absrap(), event.weight()) ;

	_h_NrmPtclDiffXsecTPt->fill(pseudoTop4p.pT(),     event.weight()) ;
	_h_NrmPtclDiffXsecTY->fill( pseudoTop4p.absrap(), event.weight()) ;

      
	
      } else {
	flag = true;
	N9;RET;//< debug mode!
	cout << "anti-top!"<< endl;
	_c_fid_tbar->fill(event.weight()) ;

	_h_AbsPtclDiffXsecTbarPt->fill(pseudoTop4p.pT()/GeV, event.weight()) ;
	_h_AbsPtclDiffXsecTbarY->fill( pseudoTop4p.absrap(), event.weight()) ;

	_h_NrmPtclDiffXsecTbarPt->fill(pseudoTop4p.pT(),     event.weight()) ;
	_h_NrmPtclDiffXsecTbarY->fill( pseudoTop4p.absrap(), event.weight()) ;

      }
      
      if (!flag){
	N1;
	cout << "should have vetoed by now!"<< endl;
	if (leptons.empty()) cout << "no leptons!" << endl;
	
	cout << "leptons[0].charge()=" << leptons[0].charge();
	  }
      
      END;
      /*


     

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
      */

      
    } //< end of analyze()            

   


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h_YYYY); // normalize to unity                                                   
      //scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section 

      double SFfbGeV = crossSection()/GeV/femtobarn/sumOfWeights();
      double SFfb = crossSection()/femtobarn/sumOfWeights();
      double SFpb = crossSection()/picobarn/sumOfWeights();
      double Norm_pT = 1/picobarn/sumOfWeights();
      double Norm_Y = 1/picobarn/sumOfWeights()/1000;
      
      scale(_c_fid_t, SFpb);
      scale(_c_fid_tbar, SFpb);

      // normalize to cross section (1,1,X)
      scale({_h_AbsPtclDiffXsecTPt,_h_AbsPtclDiffXsecTY,
	    _h_AbsPtclDiffXsecTbarPt,_h_AbsPtclDiffXsecTbarY}, SFfb);

      // normalize to unity/1000 (2,1,X)
      normalize({_h_NrmPtclDiffXsecTPt, _h_NrmPtclDiffXsecTY,
	    _h_NrmPtclDiffXsecTbarPt, _h_NrmPtclDiffXsecTbarY}, 1/1000);
      
      
      
      

      //debug info
      RET;
      cout << "debug codes info:"<< endl;
      
      cout << "#1: a worrying error!!! " << endl;
      cout << "#2: >1 jet tagged as a b-jet " << endl;
      cout << "#3: jets.size() != 2" << endl;
      cout << "#4: no b-jets found " << endl;
      cout << "#5: lepton and jet are too close " << endl;
      cout << "#6: lepton.size()!=1 " << endl;
      cout << "#7: lepton-bJet mass > 160GeV" << endl;
      cout << "#8: W-boson not reconstructed (see WFinder params)" << endl;
      cout << "#9: pass :) " << endl << endl;;
      
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
