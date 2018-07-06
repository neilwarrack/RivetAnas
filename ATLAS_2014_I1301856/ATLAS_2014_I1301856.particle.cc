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


      FinalState fs;

      // define cuts for muons, electrons and jets
      Cut electronCuts   = (Cuts::abseta < 2.47) && ( (Cuts::abseta <= 1.37) || (Cuts::abseta >= 1.52) )
	&& ( Cuts::Et > 25*GeV );
      Cut muonCuts  = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;    

      // Project out prompt electrons (including those from intermediate tau decays)
      PromptFinalState el_pfs(electronCuts, true, false); //< accepttaudecays=true, acceptmudecays=false
      IdentifiedFinalState el(el_pfs);
      el.acceptIdPair(PID::ELECTRON);
      declare(el, "Electrons");

      // Project out muons
      IdentifiedFinalState mu_fs({PID::MUON, -PID::MUON},muonCuts);
      declare(mu_fs, "Muons");
      
      // Parton Level Top Quarks
      declare(PartonicTops(PartonicTops::ELECTRON), "ElectronPartonTops");
      declare(PartonicTops(PartonicTops::MUON), "MuonPartonTops");
      declare(PartonicTops(PartonicTops::ALL), "AllPartonTops");

      // Projection for jets
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");

      // Book histograms
      _c_bQuarks = bookCounter("assumed_b-jet_counter");
      _c_bJetsSelected = bookCounter("actual_b-jet_counter");
      _c_1btag = bookCounter("1BTags");
      _c_2btag = bookCounter("2BTags-ttbarXSec");
      _c_ttbarEMuPairsOppChar = bookCounter("ttbar_xsec");
      _c_emuEvents = bookCounter("Nemu");
      //_c_error_mu = bookCounter("error_mu");
     
      //_hepmcout.reset(new HepMC::IO_GenEvent("inspectevents.hepmc"));    
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      //BEGIN;

      const Particles& electrons = apply<IdentifiedFinalState>(event, "Electrons").particlesByPt();

      const Particles& muons = apply<IdentifiedFinalState>(event, "Muons").particlesByPt();

      const Particles partonicTops_el = apply<ParticleFinder>(event, "ElectronPartonTops").particlesByPt();

      const Particles partonicTops_mu = apply<ParticleFinder>(event, "MuonPartonTops").particlesByPt();

      const Particles partonicTops_all = apply<ParticleFinder>(event, "AllPartonTops").particlesByPt();

      Jets jets = apply<FastJets>(event, "Jets").jetsByPt();


      // veto event if it is not ttbar event
      if (partonicTops_all.size() != 2) {N1; RET; vetoEvent;}

      // else count the presumed number of b-quarks/jets (2)
      _c_bQuarks->fill(event.weight());



      // count the b-jets            
      Jets b_Jets;
      for (const Jet& j : jets){
	if ( j.bTagged(Cuts::pT>5*GeV) ) b_Jets.push_back(j);
      }

      
      /*
      if (b_Jets.size() == 1) {
	N5;
	_c_1btag->fill(event.weight());
      }	 
      
      if (b_Jets.size() == 2) {
	N6;
	_c_2btag->fill(event.weight());
      }
      */

      
      Jets selectedJets;
      for (const Jet& bj : b_Jets){
	if (bj.pT() > 25*GeV && bj.abseta() < 2.5) {
	  selectedJets.push_back(bj);
	  _c_bJetsSelected->fill(event.weight());
	}
      }
      //cout << "selectedJets.size()=" << selectedJets.size() << endl;
      //cout << "_c_bJetsSelected->numEntries=" << _c_bJetsSelected->numEntries()<< endl;


      // remove leptons that are too close to jets
      Particles leptons;
      for (const Particle& e : electrons) leptons.push_back(e);
      for (const Particle& m : muons)     leptons.push_back(m);
      
      Particles selectedLeptons;
      bool tooClose = false;

      for (const Particle& l : leptons){
	for (const Jet& sj : selectedJets){
	  if (deltaR(sj, l) < 0.4) tooClose = true ;
	} 
	if (!tooClose) {
	  selectedLeptons.push_back(l);
	} else tooClose = false;
      }


      // veto if there is no opposite charged e mu pair
      if (selectedLeptons.size() != 2)
	{N2; RET; vetoEvent;}
      if (selectedLeptons[0].charge() + selectedLeptons[1].charge() != 0.0)
	{N3; RET; vetoEvent;}
      if (!(selectedLeptons[0].abspid()==PID::ELECTRON && selectedLeptons[1].abspid()==PID::MUON)
	  &&
	  !(selectedLeptons[0].abspid()==PID::MUON  && selectedLeptons[1].abspid()==PID::ELECTRON))
	{N4; RET; vetoEvent;}

      _c_emuEvents->fill(event.weight());

      // veto if there are no selected b-tagged jets
      if (selectedJets.size() == 0) {
	N5;RET;
	vetoEvent;
      }	 
      
      // veto if there are >2 selected b-tagged jets
      if (selectedJets.size() > 2) {
	N6;RET;
	vetoEvent;
      }


      
      if (selectedJets.size() == 1) {
	N7;
	_c_1btag->fill(event.weight());
      }	 
      
      if (selectedJets.size() == 2) {
	N8;
	_c_2btag->fill(event.weight());
      }
      

      _c_ttbarEMuPairsOppChar->fill(event.weight());
      RET;



    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double BR = 0.032; // branching ratio: ttbar -> e mu nu nu b b

      if ( _c_bJetsSelected->numEntries() != 0){

	double ttbarEvents = _c_bQuarks->numEntries();
	cout << "ttbarEvents=" << ttbarEvents << endl; //< debug mode!
	double bCtr = 2*ttbarEvents;
	cout << "bCtr=" << bCtr << endl; //< debug mode!
	double NemuEvents=_c_emuEvents->numEntries();
	cout << "NemuEvents=" << NemuEvents << endl; //< debug mode!
	double bJetCtr = _c_bJetsSelected->numEntries();
	cout << "bJetCtr=" << bJetCtr << endl; //< debug mode!
	double bJetProb = bJetCtr/bCtr;
	cout << "bJetProb=" << bJetProb << endl; //< debug mode!
	double emuFraction = NemuEvents/ttbarEvents;
	cout << "emuFraction=" << emuFraction << endl; //< debug mode!

	
	//	double SF = crossSection()/picobarn/sumOfWeights()/BR/(bJetProb*bJetProb);  // scale factor
	double SF1 = crossSection()/picobarn/sumOfWeights()/emuFraction/(bJetProb*bJetProb);  // scale factor (requires all ttbar events!)	
	double SF2 = crossSection()/picobarn/sumOfWeights()/0.008/(bJetProb*bJetProb);  // scale factor (requires all ttbar events!)	

	cout << "SF=" << SF2 << endl;
	scale (_c_2btag,SF2);
	
	
      } else { // print warning!
	cout << "some '1/nan' in finalize()" << endl; //< debug mode!
	cout << "CROSS-SECTIONS NOT SCALED!!!!!" << endl; //< debug mode!
      }
      
      //double b1 = 0.0, b2 = 0.0, XsecSelection = 0.0, XsecSelection2 = 0.0;
      //b1 = _c_1btag->numEntries();
      //b2 = _c_2btag->numEntries();
      //XsecSelection = (b2 + b1/2.0)*(b2 + b2/2.0)*crossSection()*picobarn/sumOfWeights()/b2;
      //XsecSelection2 = (b2/0.008 + b1/2.0)*(b2/0.008 + b2/2.0)*0.008*0.008*crossSection()*picobarn/sumOfWeights()/b2;
      //cout << "Calculated Xsec from e mu selection and btagged jets = " << XsecSelection << "pb" << endl;
      //cout << "Calculated Xsec from e mu selection and btagged jets (including e mu efficiency factor of 0.8%) = " << XsecSelection2 << "pb" << endl;


      //_hepmcout->clear(); _hepmcout.reset();

    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_ttbar, _c_1btag, _c_2btag, _c_ttbarEMuPairsOppChar, _c_error_e, _c_error_mu, _c_bQuarks, _c_bJetsSelected, _c_bJets, _c_emuEvents;
    //@}

    std::unique_ptr<HepMC::IO_GenEvent> _hepmcout;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1301856);


}
