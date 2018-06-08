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


  /// @brief
  /// 1) Total and fiducial cross-section measurements for single top quark and
  /// single top antiquark production in pp collisions at 8*TeV.
  /// 2) Ratio of measured total cross-section for top quark to top antiquark
  /// production in the t-channel (??) 
  /// 3) Normalised differential cross-sections and absolute differential
  /// cross-sections for single top quark and top antiquark production at
  /// particle and parton level as a function of quark ransverse momentum p_T
  /// and absolute value of rapidity Y.
  
  class tempTopAnalysis : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(tempTopAnalysis);

   
    /// @name Analysis methods
    //@{
   
    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;

      // Particle-level electron and muon cuts
      Cut e_muCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV ;
      MissingMomentum missmom(fs);
      declare(missmom, "Missp4");
      //      Cut eCuts  = Cuts::abseta < 2.47 && (Cuts::abseta < 1.37 || Cuts::abseta > 1.52);    
      //Cut jetCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;

      //      if (fuzzyEquals(sqrtS(), 7*TeV)) {
	//	...
      //      } else {
	//	...
      //      } 

      

      
      /// Project (particle level) prompt final state electrons and muons*
      /// *(including from tau decays but not from muon decays -- see PromptFinalState bools) 

      // Projection to dress leptons for lepton and jet definitions
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      declare(photons, "Photons");

      /// Electrons
      IdentifiedFinalState electron_fs;
      electron_fs.acceptIdPair(PID::ELECTRON);
      //      declare(electron_fs, "Electrons");
      PromptFinalState promptElectrons(electron_fs, true, false);
      DressedLeptons dressedElectrons(photons, promptElectrons, 0.1, e_muCuts, false);
      //^^use decay photons*? final state?? eh?
      //*DressedLeptons::DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
      //                                double dRmax, const Cut& cut, bool useDecayPhotons)
      declare(dressedElectrons, "DressedElectrons");

      /// Muons
      IdentifiedFinalState muon_fs;
      muon_fs.acceptIdPair(PID::MUON);
      declare(muon_fs, "Muons");
      PromptFinalState promptMuons(muon_fs, true, false);
      
      DressedLeptons dressedMuons(photons, promptMuons, 0.1, e_muCuts, false);
      declare(dressedMuons, "DressedMuons");

      
      // Projections for (particle level) W boson to define (pseudo) top quark.
      WFinder w_electron(electron_fs, 
			 e_muCuts, 
			 PID::ELECTRON, 
			 35*GeV, 8000*GeV, // Mass Window??? Max??? NB: W transverse mass > 35GeV is stated in paper which defines pseudo-top (https://arxiv.org/pdf/1502.05923.pdf)
			 30*GeV, //E_T_miss 
			 0.1, // delta R of clustering to QED dress electrons
			 WFinder::PROMPTCHLEPTONS, // can remove if not changed
			 WFinder::CLUSTERNODECAY, // this will not include non prompt photons, sould I include them? if so use: CLUSTERALL. Can remove if not changed
			 WFinder::TRACK, //can remove if not changed
			 WFinder::TRANSMASS); // tells mass window (>35*GeV) to be read as transverse mass.
      declare(w_electron, "W_Electron");

      WFinder w_muon(muon_fs,  e_muCuts, PID::MUON, 35*GeV, 100*GeV, 30*GeV, 0.1, WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERNODECAY, WFinder::TRACK, WFinder::TRANSMASS); 
      declare(w_muon, "W_Muon");
      /*
      IdentifiedFinalState leptons(fs);
      leptonfs.acceptIdPairs({PID::ELECTRON, PID::MUON});
      declare(leptonfs, "Leptons");
      PromptFinalState promptLeptons(leptonfs, true, false);// including via tau decays
      IdentifiedFinalState neutrinos(fs);
      neutrinos.acceptNeutrinos(); // How do we know these come from the w boson?!?!?
      declare(neutrinos, "Neutrinos");
      */
      
      VetoedFinalState jet_input(fs);
      jet_input.addVetoOnThisFinalState(dressedElectrons);
      jet_input.addVetoOnThisFinalState(dressedMuons);
      //jet_input.addVetoOnThisFinalState(neutrinos);
      //jet_input.vetoNeutrinos();
      declare(jet_input, "Jet_Input");

      // projection for jets
      declare(FastJets(jet_input, FastJets::ANTIKT, 0.4), "Jets");

      // Projection for Parton Level Top Quarks
      //declare(PartonicTops(PartonicTops::E_MU, e_muCuts), "LeptonicTops") ;
      //declare(PartonicTops(PartonicTops::ALL, e_muCuts), "AllPartonicTops") ;
      declare(PartonicTops(PartonicTops::E_MU), "LeptonicTops") ;
      declare(PartonicTops(PartonicTops::ALL ), "AllPartonicTops") ;
      
      // Book counters
      _c_fid_t    = bookCounter("fidTotXsectq");
      _c_fid_tbar = bookCounter("fidTotXsectbarq");

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

      _h_AbsPtonDiffXsecTPt = bookHisto1D(3,1,1);
      
      // _h_diffXsecParticlePt_tbar = bookHisto1D("diffXsecParticlePt_tbar");
      //_h_diffXsecParticleY_t = bookHisto1D("diffXsecParticleY_tq");
      //_h_diffXsecParticleY_tbar = bookHisto1D("diffXsecParticleY_tbarq");

      
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
      //cout << "sizeOfJets1=" << jets.size()<< " "; RET;
      ifilter_select(jets, jetCuts);
      //cout << "sizeOfJets2=" << jets.size()<< " "; RET;

      // Find the leptons
      //cout << "muons.size()=" << muons.size(); RET;
      //cout << "electrons.size()=" << electrons.size(); RET;
      Particles leptons;
      for (const Particle& e : electrons){
	leptons.push_back(e);
      } 
      for (const Particle& m : muons){
	leptons.push_back(m);
      }
      //      cout << "leptons.size()=" << leptons.size(); RET;

      // isolate jets and leptons
      bool tooClose= false;;
      
      for (const Jet& j : jets){

      for (const Particle& l : leptons){
	if (deltaR(j, l) < 0.4) tooClose = true ;
      }
      if (tooClose) particleVeto = true;
      tooClose = false;      
      }
      
      // find and count b-jets
      Jet bJet;
      int b=0;
      
      for (const Jet& j : jets){
	if ( j.bTagged() && j.abseta() < 2.5) {
	  b++ ;
	  bJet = j;
	}
      }

      if(leptons.size() != 1) {
	if(b!=1) N2;
	if(jets.size() != 2) N4;
	if(particleVeto) N5;
      }
      
      if ( b==1 && leptons.size() == 1 && jets.size() == 2 && !particleVeto ){
	
	// Then this implies there must be 2 jets of which only one is
	// b-tagged and there is exactly one isolate particle level
	// electron or muon in the event.
      	
	// make cut on the mass of lepton + b-jet system    
	FourMomentum lepton4p = leptons[0];
	FourMomentum bJet4p = bJet;
	FourMomentum lb4p = lepton4p + bJet4p;
	
	if ( lb4p.mass() < 160*GeV ) particleVeto = true;
	
	FourMomentum w4p;
	
	
	// cout << "W-el=" << w_elP.size() << ", " ;
	// cout << "W-mu=" << w_muP.size() << " " ;
	//cout << w_el.constituentLepton().pid() << " ";
	//	  cout << w_mu.constituentLepton().pid() ;
	//if ( w_elP.size() != w_muP.size() ) WORRY;
	
	  if ( w_elP.size() + w_muP.size() == 1){

	    if   ( w_elP.size() == 1 )    w4p = w_elP[0] ;
	    else if ( w_muP.size() == 1 ) w4p = w_muP[0] ;
	    else {particleVeto = true;}

	    if (!particleVeto){
	     
	      // define pseudoTop
	      const FourMomentum pseudoTop4p = bJet4p + w4p;
	      
	      if (leptons[0].charge() > 0){ // fill top-quark histos
		N8;
		_c_fid_t->fill(event.weight()) ;
		
		_h_AbsPtclDiffXsecTPt->fill(pseudoTop4p.pT(),     event.weight()) ;
		_h_AbsPtclDiffXsecTPt->fill(pseudoTop4p.absrap(), event.weight()) ;
		
		_h_NrmPtclDiffXsecTPt->fill(pseudoTop4p.pT(),     event.weight()) ;
		_h_NrmPtclDiffXsecTPt->fill(pseudoTop4p.absrap(), event.weight()) ;
		
		
	      } else { // Fill anti top-quark histos
		N8;
		_c_fid_tbar->fill(event.weight()) ;
		
		_h_AbsPtclDiffXsecTbarPt->fill(pseudoTop4p.pT(),     event.weight()) ;
		_h_AbsPtclDiffXsecTbarPt->fill(pseudoTop4p.absrap(), event.weight()) ;
		
		_h_NrmPtclDiffXsecTbarPt->fill(pseudoTop4p.pT(),     event.weight()) ;
		_h_NrmPtclDiffXsecTbarPt->fill(pseudoTop4p.absrap(), event.weight()) ;
		
	      }
	    }
	  }
	  else if ( w_elP.size() + w_muP.size() > 1) N1;
	} //< end of particle-level analysis 
	//cout << "particleVeto=" << particleVeto;

	
	/// Parton-level analysis
	bool partonVeto = false;
	
	// find partonic tops
	
	if ( !allPartonicTops.empty() ){
	  //cout << "-Number of Tops:" << allPartonicTops.size(); RET;
	  if ( allPartonicTops.size() != 1 ) {
	    partonVeto = true ;
	  } else {
	   
	    if (allPartonicTops[0].pid() > 0){ // fill top-quark histos
	      //cout << "q-allPartonicTops[0].pid()=" << allPartonicTops[0].pid();
	      N9;
	      //_h_AbsPtonDiffXsecTPt->fill(leptonicTops[0].pT(), event.weight());
	    } else {
	      //cout << "qbar-allPartonicTops[0].pid()=" << allPartonicTops[0].pid();
	      N9;
	    }
	    
	  }

	  
	} else {N1;}
	RET;
    }
    
    


      
      ////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////







      
      
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
      }
*/
      
    
  

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
      
      //      normalize(_h_NrmPtclDiffXsecTPt);
      //(_h_NrmPtclDiffXsecTY, SF);
      // scale(_h_NrmPtclDiffXsecTbarPt, SF);
      //scale(_h_NrmPtclDiffXsecTbarY, SF);

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
  DECLARE_RIVET_PLUGIN(tempTopAnalysis);


}
