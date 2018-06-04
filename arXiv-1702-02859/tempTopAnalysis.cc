// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/JetUtils.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PartonicTops.hh"
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


  /// @brief What it is, what it is
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

      //      if (fuzzyEquals(sqrtS(), 7*TeV)) {
	//	...
      //      } else {
	//	...
      //      } 

      
      //IdentifiedFinalState electron_fs(lepCuts, PID::ELECTRON);
      IdentifiedFinalState electron_fs;
      electron_fs.acceptIdPair(PID::ELECTRON);
      declare(electron_fs, "Electrons");

      //      IdentifiedFinalState muon_fs(muCuts, PID::MUON);
      IdentifiedFinalState muon_fs;
      muon_fs.acceptIdPair(PID::MUON);
      declare(muon_fs, "Muons");
      
      WFinder w_electron(electron_fs, 
			 lepCuts, 
			 PID::ELECTRON, 
			 35*GeV, 100*GeV, // Mass Window??? Max??? NB: W transverse mass > 35GeV is stated in paper which defines pseudo-top (https://arxiv.org/pdf/1502.05923.pdf)
			 30*GeV, //E_T_miss 
			 0.1, // delta R of clustering to QED dress electrons
			 WFinder::PROMPTCHLEPTONS, // can remove if not changed
			 WFinder::CLUSTERNODECAY, // this will not include non prompt photons, sould I include them? if so use: CLUSTERALL. Can remove if not changed
			 WFinder::TRACK, //can remove if not changed
			 WFinder::TRANSMASS); // tells mass window (>35*GeV) to be read as transverse mass.
			       
      declare(w_electron, "W_Electron");



      WFinder w_muon(muon_fs,  lepCuts, PID::MUON, 35*GeV, 100*GeV, 30*GeV, 0.1, WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERNODECAY, WFinder::TRACK, WFinder::TRANSMASS); 
			       
      declare(w_muon, "W_Muon");





      // Projection for Parton Level Top Quarks
      declare(PartonicTops(PartonicTops::E_MU, lepCuts), "leptonicTops") ;
      //      declare(PartonicTops(PartonicTops::MUON, lepCuts), "MuonPartonTops") ;
      
      // Projection for jets
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");

      // Book histograms
      _c_fidt    = bookCounter("fidTotXsectq");
      _c_fidtbar = bookCounter("fidTotXsectbarq");
      //      _c_2btag = bookCounter("2BTag");
      //_c_ttbarEMuPairsOppChar_partonic = bookCounter("emupairsOppChar_partonic");
      //_c_error_e  = bookCounter("error_e");
      //_c_error_mu = bookCounter("error_mu");
     
      _hepmcout.reset(new HepMC::IO_GenEvent("inspectevents.hepmc"));    
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      const Particles& electrons = apply<IdentifiedFinalState>(event, "Electrons").particlesByPt();
      const Particles& muons =     apply<IdentifiedFinalState>(event, "Muons").particlesByPt();

      const Particles leptonicTops = apply<ParticleFinder>(event, "leptonicTops").particlesByPt();
      //      const Particles leptons = apply<ParticleFinder>(event, "leptonTop") ;


      //const Particles muonpartontops     = apply<ParticleFinder>(event, "MuonPartonTops").particlesByPt();
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt();
      
      const WFinder& w_el = apply<WFinder>(event, "W_Electron");
      const Particles w_elP= w_el.bosons();
      const WFinder& w_mu = apply<WFinder>(event, "W_Muon");
      const Particles w_muP= w_el.bosons();
 
      
	//const FourMomentum w_mu4p = w_mu.boson() ;
      //const Particles tempparts = w_el.bosons() ;
      // cout << "tempparts.size()" << tempparts.size() <<endl;
      // const FourMomentum w_els4p = w_mu.bosons() ;

      //      cout <<" w_els4p.size()=" <<  w_els4p.size() << endl;
      // cout <<" w_mus4p.size()=" <<  w_mus4p.size() << endl;
      
      //Particles selectedElectrons, selectedMuons;
      //bool tooClose = false;
     
      
      // find jets in fiducial pT and eta ranges
      Cut jetCuts = Cuts::abseta < 4.5 && Cuts::pT > 30*GeV ;
      ifilter_select(jets, jetCuts ) ;
      if ( jets.size() != 2 ) {
	N4;
	RET;
	vetoEvent ;
      }
      // find and count b-jets
      Jet bJet;
      int b = 0 ;
      for (const Jet& j : jets){
	if ( j.bTagged() ) {
	  b++ ;
	  if (b == 2) {N5; RET; vetoEvent;}
	  bJet = j;
	}
      }
      if ( b == 0 ) {RET; vetoEvent;} // this implies there must be 2 jets of which only one is b-tagged


      FourMomentum wp4;
      if ( w_elP.size() + w_muP.size() != 0){

	//	cout << "W-el=" << w_elP.size() << " " ;
	//cout << "W-mu=" << w_muP.size() << " " ;
	if       ( w_elP.size() == 1 ) const FourMomentum w4p = w_elP[0] ;
	else if  ( w_muP.size() == 1 ) const FourMomentum w4p = w_muP[0] ;
	else N6; //hmmm
	RET;
	vetoEvent;
      }
      
      //	if (j.pT() > 25*GeV && j.abseta() < 2.5) selectedJets.push_back(j);

      // find invariant mass of lepton-b-jet system
      const FourMomentum bJt4p = bJet ;
      const FourMomentum pseudoTopp4 = bJt4p + wp4;


      //const FourMomentum lep4p = leptonicTops[0] ;
      //const FourMomentum lepBJt4p = bJt4p + lep4p ;
      //if ( lepBJt4p.mass() < 160*GeV ){
      //	N3; RET; vetoEvent ;
      //}

            
      if ( leptonicTops.size() != 1 ){ N7; RET; vetoEvent;}
      else {
	const bool charge = leptonicTops[0].charge() > 0 ;

	  if (charge){
	    _c_fidt   ->fill(event.weight());
	  } else {
	    _c_fidtbar->fill(event.weight());
	  }
      }
      


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
      double SF = crossSection()*picobarn/sumOfWeights();  // scale factor
      scale(_c_fidt, SF);
      scale(_c_fidtbar, SF);
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
    CounterPtr _c_fidt, _c_fidtbar;
    //@}

    std::unique_ptr<HepMC::IO_GenEvent> _hepmcout;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(tempTopAnalysis);


}
