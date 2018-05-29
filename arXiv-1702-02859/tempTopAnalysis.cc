// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/JetUtils.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PartonicTops.hh"
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
  class tempTopAnalysis : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(tempTopAnalysis);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;

      Cut lepCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV ;
     
      // Cut muCuts  = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;    
      //Cut jetCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;


      //IdentifiedFinalState electron_fs(lepCuts, PID::ELECTRON);
      IdentifiedFinalState electron_fs;
      electron_fs.acceptIdPair(PID::ELECTRON);
      declare(electron_fs, "Electrons");

      //      IdentifiedFinalState muon_fs(muCuts, PID::MUON);
      IdentifiedFinalState muon_fs;
      muon_fs.acceptIdPair(PID::MUON);
      declare(muon_fs, "Muons");
      
      // Projection for Parton Level Top Quarks
      declare(PartonicTops(PartonicTops::E_MU, lepCuts), "leptonicTops") ;
      //      declare(PartonicTops(PartonicTops::MUON, lepCuts), "MuonPartonTops") ;
      
      // Projection for jets
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");

      // Book histograms
      _c_ttbar = bookCounter("ttbarXsec");
      _c_1btag = bookCounter("1BTag");
      _c_2btag = bookCounter("2BTag");
      _c_ttbarEMuPairsOppChar_partonic = bookCounter("emupairsOppChar_partonic");
      _c_error_e  = bookCounter("error_e");
      _c_error_mu = bookCounter("error_mu");
     
      _hepmcout.reset(new HepMC::IO_GenEvent("inspectevents.hepmc"));    
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      const Particles& electrons = apply<IdentifiedFinalState>(event, "Electrons").particlesByPt();
      const Particles& muons =     apply<IdentifiedFinalState>(event, "Muons").particlesByPt();

      const Particles leptonicTops = apply<ParticleFinder>(event, "leptonicTops").particlesByPt();
      //      const Particles leptons = apply<ParticleFinder>(event, "leptonTop") ;
      if ( leptons.size() != 1 ) vetoEvent;

      //const Particles muonpartontops     = apply<ParticleFinder>(event, "MuonPartonTops").particlesByPt();
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt();



      //Particles selectedElectrons, selectedMuons;
      //bool tooClose = false;
     
      
      // find jets in fiducial pT and eta ranges
      Cut jetCuts = Cuts::abseta < 4.5 && Cuts::pT > 30*GeV ;
      ifilter_select(jets, jetCuts ) ;
      if ( jets.size() != 2 ) vetoEvent ;
      
      // find and count b-jets
      Jet bJet;
      int b = 0 ;
      for (const Jet& j : jets){
	if ( j.bTagged() ) {
	  b++ ;
	  if (b == 2) vetoEvent;
	  bJet = j;
	}
      }

      
      //	if (j.pT() > 25*GeV && j.abseta() < 2.5) selectedJets.push_back(j);

      // find invariant mass of lepton-b-jet system
      const FourMomentum bJt4p = bJet ;
      const FourMomentum lep4p = leptons[0] ;
      const FourMomentum lepBJt4p = bJt4p + lep4p ;
      if ( lepBJt4p.mass() < 160*GeV ) vetoEvent ;
      




      /*
    //      cout << "re:" << electrons.size() << " rm:" << muons.size() << " ";
      //      cout << "pe:" << electronpartontops.size() << " pm:" << muonpartontops.size() << " ";
  

  
    for (const Particle& e : electrons){
	for (const Jet& sj : selectedJets){
	  if (deltaR(sj, e) < 0.4) tooClose = true ;
	} 
	if (!tooClose) selectedElectrons.push_back(e);
	tooClose = false;      
      }


      for (const Particle& m : muons){
	for (const Jet& sj : selectedJets){
	  if (deltaR(sj, m) < 0.4) tooClose = true ;
	} 
	if (!tooClose) selectedMuons.push_back(m);
	tooClose = false;
      }

      //      int N_1=0, N_2=0;

      Jets bTaggedJets;
      if ((selectedElectrons.size() == 1 ) && (selectedMuons.size() == 1)){
	//N6;
	for (const Jet& bj : selectedJets){
	  if (bj.bTagged()) {
	    bTaggedJets.push_back(bj);
	  }
	}
	
	if (bTaggedJets.size() == 1) {
	  //N7;
	  _c_1btag->fill(event.weight());
	}	 

	if (bTaggedJets.size() == 2) {
	  //N8;
	  _c_2btag->fill(event.weight());
	}
      }

  
      cout << "______________" << endl;
      cout << "final_st: fj:" << jets.size()         << " fe:" << electrons.size()         << " fm:"<< muons.size()         << " " << endl;
      cout << "selected: sj:" << selectedJets.size() << " se:" << selectedElectrons.size() << " sm:"<< selectedMuons.size() << " " << endl;
      cout << "partonic:      pe:" << electronpartontops.size() << " sm:"<< muonpartontops.size() << " " << endl;
  

  
      if (electronpartontops.size() > electrons.size()){cout <<         "                             pe>fe" << endl;
	_c_error_e->fill(event.weight()); //error counter
	// print event record to file
	const HepMC::GenEvent* ge = event.genEvent();
	if (ge == nullptr) {cout << "nullpointer found: GenEvent*"<<endl;
	} else { //if (ge != nullptr) {
	  ge->print(); // sends (formatted) HepMC event to screen
	  _hepmcout->write_event(ge); // writes to file
	}
      }
      if (muonpartontops.size() > muons.size()){cout <<                 "                             pm>fm" << endl;
	_c_error_mu->fill(event.weight()); //error counter
      }
      if (selectedElectrons.size() > electronpartontops.size()) cout << "                             se>pe" << endl;      
      if (selectedMuons.size() > muonpartontops.size()) cout <<         "                             se>pe" << endl;      
  

  
      if ((selectedElectrons.size() == 1 ) && (selectedMuons.size() == 1)) N6;
      if (bTaggedJets.size() == 1) N7;
      if (bTaggedJets.size() == 2) N8;
  

      // veto if not an e mu pair
      const bool eMuPair = ( electronpartontops.size() == 1 && muonpartontops.size() == 1 );
      if ( !eMuPair ){ N1; RET; vetoEvent; }
      else { 

	// Fill histos
	_c_ttbarEMuPairsOppChar_partonic->fill(event.weight());
	
	// veto if e and mu are not of opposite chare
	const bool oppositeCharge ( electronpartontops[0].charge() + muonpartontops[0].charge() == 0 );

	if ( oppositeCharge ) {
	  N9;
	  _c_ttbar->fill(event.weight());
	} else { N2; RET; vetoEvent; }

      }
*/
      N9;
      RET;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      /*
      double BR = 0.032; // branching ratio
      double SF = crossSection()*picobarn/sumOfWeights()/BR;  // scale factor
     
      scale(_c_ttbar, SF);
      scale(_c_ttbarEMuPairsOppChar_partonic,SF);

      double b1 = 0.0, b2 = 0.0, XsecSelection = 0.0, XsecSelection2 = 0.0;
      b1 = _c_1btag->numEntries();
      b2 = _c_2btag->numEntries();
      XsecSelection = (b2 + b1/2.0)*(b2 + b2/2.0)*crossSection()*picobarn/sumOfWeights()/b2;
      XsecSelection2 = (b2/0.008 + b1/2.0)*(b2/0.008 + b2/2.0)*0.008*0.008*crossSection()*picobarn/sumOfWeights()/b2;
      //cout << "Calculated Xsec from e mu selection and btagged jets = " << XsecSelection << "pb" << endl;
      //cout << "Calculated Xsec from e mu selection and btagged jets (including e mu efficiency factor of 0.8%) = " << XsecSelection2 << "pb" << endl;


      _hepmcout->clear(); _hepmcout.reset();
      */
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_ttbar, _c_1btag, _c_2btag, _c_ttbarEMuPairsOppChar_partonic, _c_error_e, _c_error_mu;
    //@}

    std::unique_ptr<HepMC::IO_GenEvent> _hepmcout;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(tempTopAnalysis);


}
