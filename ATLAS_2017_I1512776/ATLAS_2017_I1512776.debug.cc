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


  /// @brief Total and fiducial cross-section measurements with normalised differential cross-sections and absolute differential cross-sections for single top and anti-top production in pp collisions at 8*TeV. This analysis requires ONLY t-channel single top events as input. TODO: make general so analysis can take t-channel single top events + other stuff, as input.
  
  
  class ATLAS_2017_I1512776 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1512776);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;

      
      // Particle-level lepton cuts
      Cut e_muCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV ;

      
      // photons for dressing leptons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      declare(photons, "Photons");

      
      // particle-level electrons (assumed to come from a W decay)
      IdentifiedFinalState electron_fs;
      electron_fs.acceptIdPair(PID::ELECTRON);
      PromptFinalState promptElectrons(electron_fs, true); //< "true" accepts as 'prompt' particles from prompt taus
      DressedLeptons dressedElectrons(photons, promptElectrons, 0.1, e_muCuts, false);
      // NB: the above bool ("false") avoids dressing using 'decay photons'...
      // NB: 'Decay photons' are defined in depricated manner, using:
      // NB: fromDecay() const { return fromHadron() || fromPromptTau(); } 
      declare(dressedElectrons, "DressedElectrons");


      // particle-level muons (assumed to come from a W decay)
      IdentifiedFinalState muon_fs;
      muon_fs.acceptIdPair(PID::MUON);
      PromptFinalState promptMuons(muon_fs, true); 
      DressedLeptons dressedMuons(photons, promptMuons, 0.1, e_muCuts, false); 
      declare(dressedMuons, "DressedMuons");


      // Particle-level W (with electron decay)
      WFinder w_el(fs, e_muCuts,PID::ELECTRON, 35*GeV, 8000*GeV, 30*GeV, 0.1,
		   WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERNODECAY,
		   WFinder::TRACK,  WFinder::TRANSMASS);
      // notes: ET_miss>30GeV (see: https://arxiv.org/pdf/1502.05923.pdf)
      // notes: mass window > 35GeV (< transverse mass) (see: https://arxiv.org/pdf/1502.05923.pdf)
      declare(w_el, "W_Electron");


      // Particle-level W (with muon decay)
      WFinder w_mu(fs, e_muCuts, PID::MUON, 35*GeV, 8000*GeV, 30*GeV, 0.1,
		     WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERNODECAY,
		     WFinder::TRACK, WFinder::TRANSMASS);
      declare(w_mu, "W_Muon");
      
      
      //the selected lepton (with it's photons) is not to be included in jet clustering
      VetoedFinalState jet_input(fs);
      //declare(jet_input, "JetInput");
      jet_input.addVetoOnThisFinalState(dressedElectrons);
      jet_input.addVetoOnThisFinalState(dressedMuons);
      declare(FastJets(jet_input, FastJets::ANTIKT, 0.4), "Jets");

      // The missing Et requirement
      MissingMomentum mET;
      declare(mET, "MET");
      
      // for parton level analysis
      declare(PartonicTops(PartonicTops::E_MU), "LeptonicTops") ;//including from taus
      declare(PartonicTops(PartonicTops::ALL ), "AllPartonicTops") ;
      declare(PartonicTops(PartonicTops::HADRONIC, true, true), "HadTops") ; // include hadronic taus
      //declare(PartonicTops(PartonicTops::TAU), "TauTops") ;
      
      // book counters
      _c_N_tot =  bookCounter("mytest1");
      _c_N_tot_tq =  bookCounter("mytest11");
      _c_N_tot_tbarq =  bookCounter("mytest111");
      _c_N_fid =  bookCounter("mytest2");
      _c_N_fid_tq =  bookCounter("mytest22");
      _c_N_fid_tbarq =  bookCounter("mytest22");
      _c_fid_t =    bookCounter("fiducialXSec_top");
      _c_fid_tbar = bookCounter("fiducialXSec_anti-top");
      _c_Xsec_fid_tq = bookCounter("Xsec_fir_tq");
      _c_Xsec_fid_tbarq = bookCounter("Xsec_fir_tbarq");

      
      // book histograms
      // top quark differential cross sections at particle level
      _h_AbsPtclDiffXsecTPt   = bookHisto1D(1,1,1);
      _h_AbsPtclDiffXsecTY    = bookHisto1D(1,1,3);
      _h_NrmPtclDiffXsecTPt   = bookHisto1D(2,1,1);
      _h_NrmPtclDiffXsecTY    = bookHisto1D(2,1,3);
      
      // top antiquark differential cross sections at particle level    
      _h_AbsPtclDiffXsecTbarPt   = bookHisto1D(1,1,2);
      _h_AbsPtclDiffXsecTbarY    = bookHisto1D(1,1,4);
      _h_NrmPtclDiffXsecTbarPt   = bookHisto1D(2,1,2);
      _h_NrmPtclDiffXsecTbarY    = bookHisto1D(2,1,4);

      // top quark differential cross sections at parton level
      _h_AbsPrtnDiffXsecTPt   = bookHisto1D(3,1,1);
      _h_AbsPrtnDiffXsecTY    = bookHisto1D(3,1,3);
      _h_NrmPrtnDiffXsecTPt   = bookHisto1D(4,1,1);
      _h_NrmPrtnDiffXsecTY    = bookHisto1D(4,1,3);
      
      // top antiquark differential cross sections at parton level    
      _h_AbsPrtnDiffXsecTbarPt   = bookHisto1D(3,1,2);
      _h_AbsPrtnDiffXsecTbarY    = bookHisto1D(3,1,4);
      _h_NrmPrtnDiffXsecTbarPt   = bookHisto1D(4,1,2);
      _h_NrmPrtnDiffXsecTbarY    = bookHisto1D(4,1,4);
      

      
      // debug hepmc print-out
      _hepmcout.reset(new HepMC::IO_GenEvent("inspectevents.hepmc"));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //ignore:
      //      if (Analysis.isCompatible(1)) {cout << "compatible!"; RET;}
      BEGIN;
      N9;

     
      // count events and perform fast veto if not a single top event
      // (NB: do not give this analysis s-channel single top events as input!)
      const Particles allTops = apply<ParticleFinder>(event, "AllPartonicTops").particlesByPt();
      if ( allTops.size() == 1 ){
	_c_N_tot->fill(event.weight());
	if ( allTops[0].charge() > 0 ) _c_N_tot_tq->fill(    event.weight() );
	if ( allTops[0].charge() < 0 ) _c_N_tot_tbarq->fill( event.weight() );
      }
      else {
	cout << "not a single top event! VETOING!!";RET;
	vetoEvent;
      }

      
      // Find leptons, partonic tops, missing energy and jets in event
      const Particles& electrons = apply<DressedLeptons>(event, "DressedElectrons").particlesByPt();
      const Particles& muons = apply<DressedLeptons>(event, "DressedMuons").particlesByPt();
      const WFinder& w_el = apply<WFinder>(event, "W_Electron");
      const WFinder& w_mu = apply<WFinder>(event, "W_Muon");
      const Particles leptonicTops = apply<ParticleFinder>(event, "LeptonicTops").particlesByPt();
      const Particles hadTops = apply<ParticleFinder>(event, "HadTops").particlesByPt();
      const MissingMomentum misMom = apply<MissingMomentum>(event, "MET");
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(); //< non-const for filtering
      

      // START debug info //////////////////////////////////////////////////
      cout << "EVENT INFO:";RET;
      cout << "weight: ";RET;
      cout << event.weight();RET;
      //      if ( event.weight() < 0.0) {N1;RET;} else {N2;RET;}
      if ( allTops.size() != 1){
	cout << "WOW! partonicTops(ALL).size()=" << allTops.size(); RET;
      }
      else {
	//cout << "one partonic top found";RET;N5;RET;
      }
      if (leptonicTops.size() == 1){
	//cout << "single partonic top (e/mu) found";RET;N4;RET;
      }
      else {
	cout << "e/muTop.size() = " << leptonicTops.size();RET;
      }
      if (hadTops.size() == 1){
	//      	cout << "single hadronic top found";RET;N6;RET;
      }
      //      if (tauTops.size() == 1){
      //cout << "single tau top found";RET;N7;RET;
      // }
      // END debug info //////////////////////////////////////////////////

    

      
 
      /// require exactly one b-jet and at least on other for particle level results
      // apply particle-level jet cuts
      Cut jetCuts = (Cuts::abseta < 4.5 && Cuts::pT > 30*GeV);
      cout << "sizeOfJetsPreFilter=" << jets.size(); RET;
      ifilter_select(jets, jetCuts);
      cout << "sizeOfJetsPostFilter=" << jets.size(); RET;
      
      // find b-jets
      FourMomentum bJet;
      int bJetsFound = 0;
      for (const Jet& j : jets){
	//if (bJetFound == 2) {cout << "More than one analysis level b-jet";RET;break;}
	if ( j.bTagged()) { cout << "b-tagged jet found";RET;}
	if ( j.bTagged(Cuts::pT > 5*GeV)) {cout << "b-tagged jet found with b>5GeV found";RET;}
	if ( j.bTagged(Cuts::pT > 5*GeV) && j.abseta() < 2.5) {
	  bJetsFound++;
	  bJet = j;
	}
      }
      cout << "Number of b-jets:" << bJetsFound;RET;
    
      bool jetSelection = false;
      if ( jets.size() == 2 && bJetsFound == 1 ) {
	cout << "two jets found (exactly one of which is a b-jet)!"; RET;
	jetSelection = true;
      }


      // missing energy cut for parton and particle level analysis
      bool metSelection = false;
      if (misMom.met() < 30*GeV) {
	cout << "met failure: met (GeV) = " << misMom.met(); RET;
      }
      else {
	metSelection = true;
      }


      // Count isolated leptons
      cout << "muons.size()=" << muons.size(); RET;
      cout << "electrons.size()=" << electrons.size(); RET;
      Particles leptons;
      for (const Particle& e : electrons){
        leptons.push_back(e);
      }
      for (const Particle& m : muons){
        leptons.push_back(m);
      }
    
      bool leptonSelection = false;
      cout << "leptons.size()=" << leptons.size(); RET;
      if ( leptons.size() == 1 ){
	leptonSelection = true;
      }
      

      
      // discared event if selected jets and single lepton are not isolated
      bool jl_tooClose= false;
      if (leptonSelection){
        for (const Jet& j : jets){
	  if (deltaR(j, leptons[0]) < 0.4) jl_tooClose = true ;
	}
      }


      // calculate mass of lepton + b-jet system
      bool massljSelection = false;
      if (leptonSelection && jetSelection){ // guarantees exactly 1 lepton and exactly 1 bJet 
	FourMomentum lepton4p = leptons[0];
	FourMomentum bJet4p = bJet;
	cout << "bJet4p.pT=" << bJet4p.pT() << endl;
	cout << "bJet4p.absrap=" << bJet4p.absrap() << endl;
	cout << "lepton4p.pT=" << lepton4p.pT() << endl;
	cout << "lepton4p.absrap=" << lepton4p.absrap() << endl;
	FourMomentum lb4p = add(lepton4p, bJet4p);
	if (lb4p.mass() < 160*GeV){
	  massljSelection = true;
	}
      }


      // construct pseudo-W

      // TODO: what if more than one W reconstructed?!
      
      FourMomentum w4p;
      bool WbosonFound = false;
     
            
      if ( w_el.bosons().size() != 0 ){
	WbosonFound = true;
	const Particles W_el= w_el.bosons();

	//debug info /// START /////////////////////////////////////
	cout << "w_el.bosons().size()=" << w_el.bosons().size(); RET;
	cout << "W_el.size()=" << W_el.size() << ", " ;
	cout << "W_el.constituentLepton().pid()=" << w_el.constituentLepton().pid() << ", ";
	//debug info /// END ////////////////////////////////////////

	w4p = W_el[0];
	cout <<"W reconstructed from electron (naming it w4p)!" << endl;
	cout << "w4p.pT=" << w4p.pT() << endl;
	cout << "w4p.absrap=" << w4p.absrap() << endl;
	
      }

      if ( w_mu.bosons().size() != 0 ){
	if ( WbosonFound ) {N8;RET;vetoEvent;} // not sure which W to use in pseudo-top definition
	else WbosonFound = true;
	const Particles W_mu= w_mu.bosons();

	//debug info /// START /////////////////////////////////////
	cout << "w_mu.bosons().size()=" << w_mu.bosons().size(); RET;
	cout << "W_mu.size()=" << W_mu.size() << ", " ;
	cout << "W_mu.constituentLepton().pid()=" << w_mu.constituentLepton().pid() << ", ";
	//debug info /// END ////////////////////////////////////////

	w4p = W_mu[0];
	cout <<"W reconstructed from muon (naming it w4p)!" << endl;
	cout << "w4p.pT=" << w4p.pT() << endl;
	cout << "w4p.absrap=" << w4p.absrap() << endl;
	
      }
	/*
          
      if ( W_el.size() + W_mu.size() == 1){
	if (W_el.size() == 1 ) {
	  cout << "W_el[0].pT=" << W_el[0].pt() << endl;
	  cout << "W_el[0].absrap=" << W_el[0].absrap() << endl;
	  cout << "W_el[0].E=" << W_el[0].E() << endl;
	  cout << "W_el[0].px=" << W_el[0].px() << endl;
	  cout << "W_el[0].py=" << W_el[0].py() << endl;
	  cout << "W_el[0].pz=" << W_el[0].pz() << endl;
	}
	else if ( W_mu.size() == 1 ) {
	  cout << "W_mu[0].pT=" << W_mu[0].pt() << endl;
	  cout << "W_mu[0].absrap=" << W_mu[0].absrap() << endl;
	  
	  w4p = W_mu[0] ;
	  cout <<"W reconstructed from muon (naming it w4p)!" << endl;
	  cout << "w4p.pT=" << w4p.pT() << endl;
	  cout << "w4p.absrap=" << w4p.absrap() << endl;	  
	}
	else {
	  cout << "W_el.size()=" << W_el.size() << " " <<  "W_mu.size()=" << W_mu.size() << " " << endl; 
	}
      }
      */

    
      // construct pseudo-top
      bool pseudoTopFound = false;
      FourMomentum pseudoTop4p;
      
      if ( bJetsFound == 1 && WbosonFound )
	{
	  pseudoTop4p = add(bJet, w4p);
	  pseudoTopFound = true;
	  cout << "pseudoTop4p.pT()=" << pseudoTop4p.pT(); RET;
	  cout << "pseudoTop4p.absrap()=" << pseudoTop4p.absrap(); RET;
	}
      
      
      // debug
      bool makesSelection = false;
      if ( jetSelection && leptonSelection && massljSelection && !jl_tooClose ){
	makesSelection = true;
      } 

      if (makesSelection && !WbosonFound && !pseudoTopFound){
	N1;RET;//?
	if ( metSelection ) {N5;RET;} else {N6;RET;}
      }
      else if ( makesSelection && WbosonFound && !pseudoTopFound){
	N2;RET; //b-jet finding strangenes (finds W presumably from top but not the bjet)
      }
      else if ( makesSelection && !WbosonFound && pseudoTopFound){
	N3;RET; //W finding strangeness strangenes (finds the bjet from a top but not the W)
      }
      else if ( makesSelection && WbosonFound && pseudoTopFound){
	N4;RET;
      }

      

      // for fiducial selections fill fid cross section counters and particle-level plots
      if ( jetSelection && leptonSelection && massljSelection &&
	   !jl_tooClose && WbosonFound && pseudoTopFound ){
	
	
	_c_N_fid->fill(event.weight());

		
	if ( leptons[0].charge() > 0 ) { // top quark 

	  _c_N_fid_tq->fill( event.weight() ); // counter for normalising hostograms
	  
	  _h_AbsPtclDiffXsecTPt->fill( pseudoTop4p.pT(),     event.weight() ); //(1,1,1)  
	  _h_AbsPtclDiffXsecTY->fill(  pseudoTop4p.absrap(), event.weight() ); //(1,1,3)
	  _h_NrmPtclDiffXsecTPt->fill( pseudoTop4p.pT(),     event.weight() ); //(2,1,1)
	  _h_NrmPtclDiffXsecTY->fill(  pseudoTop4p.absrap(), event.weight() ); //(2,1,3)
	  
	}
	
	if ( leptons[0].charge() < 0 ) { // top antiquark
	  
	  _c_N_fid_tbarq->fill( event.weight() ); // counter for normalising hostograms
	  
	  _h_AbsPtclDiffXsecTbarPt->fill( pseudoTop4p.pT(),     event.weight() );  //(1,1,2)
	  _h_AbsPtclDiffXsecTbarY->fill(  pseudoTop4p.absrap(), event.weight() );  //(1,1,4)
	  _h_NrmPtclDiffXsecTbarPt->fill( pseudoTop4p.pT(),     event.weight() );  //(2,1,2)
	  _h_NrmPtclDiffXsecTbarY->fill(  pseudoTop4p.absrap(), event.weight() );  //(2,1,4)
	  
	}
	
      }
	
      
      /// Parton-level analysis  ///////////////////////////////////////////////
      
      
      if (leptonicTops.size() == 1){ // do partonic analysis
	
	
	if (leptonicTops[0].charge() > 0){
	  
	  
	  cout << "add to plots: parton-level t!";RET;

	  
	  //	  _c_Xsec_fid_tq->fill(event.weight()) ;
	  
	  _h_AbsPrtnDiffXsecTPt->fill( leptonicTops[0].pT()    , event.weight()) ;
	  _h_AbsPrtnDiffXsecTY->fill(  leptonicTops[0].absrap(), event.weight()) ;
	  
	  _h_NrmPrtnDiffXsecTPt->fill( leptonicTops[0].pT(),     event.weight()) ;
	  _h_NrmPrtnDiffXsecTY->fill(  leptonicTops[0].absrap(), event.weight()) ;
	  
	  
	}
	else {
	  
	  cout << "add to plots: parton-level t~";RET;
	  
	 	  //  _c_Xsec_fid_tbarq->fill(event.weight()) ;
	  
	  _h_AbsPrtnDiffXsecTbarPt->fill( leptonicTops[0].pT()    , event.weight()) ;
	  _h_AbsPrtnDiffXsecTbarY->fill(  leptonicTops[0].absrap(), event.weight()) ;
	  
	  _h_NrmPrtnDiffXsecTbarPt->fill( leptonicTops[0].pT(),     event.weight()) ;
	  _h_NrmPrtnDiffXsecTbarY->fill(  leptonicTops[0].absrap(), event.weight()) ;
	  
	}
      }

      
      // debug START /////////////////////////////////////////////
      if ( WbosonFound ){
	if ( makesSelection ){
	  cout << "report finds W and makes selection!";
	}
	else {
	  cout << "report finds W but fails selection!";
	}
      }
      else {
	if ( makesSelection ){
	  cout << "report doesn't find W but makes selection!";
	  N8;RET;
	}
	else {
	  cout << "report doesn't find W - doesn't makes selection!";
	}
      }
      // debug end ///////////////////////////////////////////////

   
    }   //< end of analyze()            
    
    


    /// Normalise histograms etc., after the run
    void finalize() {
      cout << "BEGIN FINALIZE"; RET;
      /// compute cross sections
      // tq+tbarq inclusive
      double Xsec_tqtbarq = crossSection()*_c_N_tot->val()/sumOfWeights();
      //info:
      cout << "crossSection()=" << crossSection();RET;
      cout << "_c_N_tot->val()=" << _c_N_tot->val();RET;
      cout << "sumOfWeights()=" << sumOfWeights();RET;
      cout << "Xsec_tqtbarq=" << Xsec_tqtbarq;RET;
      
      // tq inclusive
      scale( _c_N_tot_tq, Xsec_tqtbarq/_c_N_tot->val() );
      //info:
      cout << "_c_N_tot_tq->val()=" << _c_N_tot_tq->val();RET;
    
      // tbarq inclusive
      scale( _c_N_tot_tbarq, Xsec_tqtbarq/_c_N_tot->val() );
      //info:
      cout << "_c_N_tot_tbarq->val()=" << _c_N_tot_tbarq->val();RET;
          
      
      //normalize(_h_YYYY); // normalize to unity                                                   
      //scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section 

      double SFfbGeV = crossSection()/GeV/femtobarn/sumOfWeights();
      double SFfb = crossSection()/femtobarn/sumOfWeights();
      double SFpb = crossSection()/picobarn/sumOfWeights();
      double Norm_pT = 1/picobarn/sumOfWeights();
      double Norm_Y = 1/picobarn/sumOfWeights()/1000;
      
      scale(_c_fid_t, SFpb);
      scale(_c_fid_tbar, SFpb);
      scale(_c_Xsec_fid_tq, crossSection()/sumOfWeights() );
      scale(_c_Xsec_fid_tbarq, crossSection()/sumOfWeights() );
     
      // Abs cross Sections (normalized to cross section)
      scale(_h_AbsPrtnDiffXsecTPt,       _c_Xsec_fid_tq->val()/sumOfWeights() );
      scale(_h_AbsPrtnDiffXsecTbarPt, _c_Xsec_fid_tbarq->val()/sumOfWeights() );

      scale(_h_AbsPrtnDiffXsecTY,       _c_Xsec_fid_tq->val()/sumOfWeights() );
      scale(_h_AbsPrtnDiffXsecTbarY, _c_Xsec_fid_tbarq->val()/sumOfWeights() );

      // Nrm cross Sections (normalized to 1)
      scale(_h_NrmPrtnDiffXsecTPt,       1/sumOfWeights() );
      scale(_h_NrmPrtnDiffXsecTbarPt,    1/sumOfWeights() );
      
      scale(_h_NrmPrtnDiffXsecTY,    0.001/sumOfWeights() );
      scale(_h_NrmPrtnDiffXsecTbarY, 0.001/sumOfWeights() );      
    

      // end of anaysis!

      
      // Print info (just for me)
      ////////////////////////////////////////////////////
      RET;
      
      //cout << "picobarn = " << picobarn;RET;
      //cout << "femtobarn = " << femtobarn;RET;
      cout << "sum of weights = " << sumOfWeights(); RET;
      cout << "cross section = " << crossSection(); RET;
      cout << "_c_Xsec_fid_tq = " << _c_Xsec_fid_tq->val(); RET;
      cout << "_c_Xsec_fid_tbarq = " << _c_Xsec_fid_tbarq->val(); RET;
      /////////////////////////////////////////////////////
      RET;
           
    }
    
    //@}
    
    
    /// @name Histograms
    //@{
    
    CounterPtr _c_fid_t, _c_fid_tbar, _c_Xsec_fid_tq, _c_Xsec_fid_tbarq, _c_N_tot, _c_N_fid, _c_N_fid_tq, _c_N_fid_tbarq, _c_N_tot_tq, _c_N_tot_tbarq;
    Histo1DPtr _h_AbsPtclDiffXsecTPt, _h_AbsPtclDiffXsecTY, _h_NrmPtclDiffXsecTPt, _h_NrmPtclDiffXsecTY, _h_AbsPtclDiffXsecTbarPt, _h_AbsPtclDiffXsecTbarY, _h_NrmPtclDiffXsecTbarPt,  _h_NrmPtclDiffXsecTbarY, _h_AbsPrtnDiffXsecTPt, _h_AbsPrtnDiffXsecTY, _h_NrmPrtnDiffXsecTPt, _h_NrmPrtnDiffXsecTY, _h_AbsPrtnDiffXsecTbarPt, _h_AbsPrtnDiffXsecTbarY, _h_NrmPrtnDiffXsecTbarPt, _h_NrmPrtnDiffXsecTbarY;
    
    //@}
    
    std::unique_ptr<HepMC::IO_GenEvent> _hepmcout;
    
    };


// The hook for the plugin system
DECLARE_RIVET_PLUGIN(ATLAS_2017_I1512776);

}

