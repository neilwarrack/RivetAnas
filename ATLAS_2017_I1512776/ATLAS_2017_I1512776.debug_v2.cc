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
      
      DressedLeptons dressedElectrons(photons, promptElectrons, 0.1); // no cuts and last bool (see below) is false by default anyway)
      //DressedLeptons dressedElectrons(photons, promptElectrons, 0.1, e_muCuts, false);
      // NB: the above bool ("false") avoids dressing using 'decay photons'...
      // NB: 'Decay photons' are defined in depricated manner, using:
      // NB: fromDecay() const { return fromHadron() || fromPromptTau(); } 
      declare(dressedElectrons, "DressedElectrons");


      // particle-level muons (assumed to come from a W decay)
      IdentifiedFinalState muon_fs;
      muon_fs.acceptIdPair(PID::MUON);
      PromptFinalState promptMuons(muon_fs, true); 
      DressedLeptons dressedMuons(photons, promptMuons, 0.1); 
      //DressedLeptons dressedMuons(photons, promptMuons, 0.1, e_muCuts, false); 
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
      _c_sumw_singleTop =  bookCounter("mytest1");
      _c_N_tot_tq =  bookCounter("mytest11");
      _c_N_tot_tbarq =  bookCounter("mytest111");
      _c_sumw_fid =  bookCounter("mytest2");
      _c_sumw_fid_tq =  bookCounter("mytest22");
      _c_sumw_fid_tbarq =  bookCounter("mytest22");
      _c_sumw_prtn_tq =  bookCounter("sumw_parton_tq");
      _c_sumw_prtn_tbarq =  bookCounter("sumw_parton_tbarq");
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

      // 2nd jet differential cross sections at particle level in top quark channel 
      _h_AbsPtclDiffXsecTjPt  = bookHisto1D(5,1,1);
      _h_AbsPtclDiffXsecTjY   = bookHisto1D(5,1,3);
      _h_NrmPtclDiffXsecTjPt  = bookHisto1D(6,1,1);
      _h_NrmPtclDiffXsecTjY   = bookHisto1D(6,1,3);

      // top antiquark differential cross sections at particle level    
      _h_AbsPtclDiffXsecTbarPt   = bookHisto1D(1,1,2);
      _h_AbsPtclDiffXsecTbarY    = bookHisto1D(1,1,4);
      _h_NrmPtclDiffXsecTbarPt   = bookHisto1D(2,1,2);
      _h_NrmPtclDiffXsecTbarY    = bookHisto1D(2,1,4);

      // 2nd jet differential cross sections at particle level in top quark channel 
      _h_AbsPtclDiffXsecTbarjPt  = bookHisto1D(5,1,2);
      _h_AbsPtclDiffXsecTbarjY   = bookHisto1D(5,1,4);
      _h_NrmPtclDiffXsecTbarjPt  = bookHisto1D(6,1,2);
      _h_NrmPtclDiffXsecTbarjY   = bookHisto1D(6,1,4);

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
     
      // count events and perform fast veto if not a single top event
      // (NB: do not give this analysis s-channel single top events as input!)
      const Particles allTops = apply<ParticleFinder>(event, "AllPartonicTops").particlesByPt();
      if ( allTops.size() == 1 ){
	
	// count single top events
	if ( allTops[0].charge() > 0 ) _c_N_tot_tq->fill(    event.weight() );
	// else count single anti-top events
	else _c_N_tot_tbarq->fill( event.weight() );
      }

      else { // if allTops.size() != 1
	cout << "not a single top event! VETOING!!";RET;
	N8; RET;
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
	cout << "one partonic top found";RET;
      }
      
      if (leptonicTops.size() == 1){
	cout << "leptonic event! (e/mu partonic top found)";RET;
	N9;RET;
      }
      else {
	cout << "e/muTop.size() = " << leptonicTops.size();RET;
      }
      if (hadTops.size() == 1){
	cout << "single hadronic top found";RET;
      }
      // END debug info //////////////////////////////////////////////////

    

      
 
      /// require exactly one bJet and one non-bJet jet for particle level results
      // apply particle-level jet cuts
      Cut jetCuts = (Cuts::abseta < 4.5 && Cuts::pT > 30*GeV);
      cout << "sizeOfJetsPreFilter=" << jets.size(); RET;
      ifilter_select(jets, jetCuts);
      cout << "sizeOfJetsPostFilter=" << jets.size(); RET;
      
      // find and bJets
      Jet bJet;
      Jets noBtag;
      int bJetsFound = 0;
      for (const Jet& j : jets){
	
	if ( j.bTagged()) { cout << "b-tagged jet found";RET;}
	if ( j.bTagged(Cuts::pT > 5*GeV)) {cout << "b-tagged jet found with b>5GeV found";RET;}
	if ( j.bTagged(Cuts::pT > 5*GeV) && j.abseta() < 2.5) {
	  cout << "b-tagged jet found with b.pT>5GeV inside eta<2.5";RET;
	  bJetsFound++;
	  bJet = j;
	} else {
	  noBtag.push_back(j);
	}
      }

      cout << "Total number of particle-level b-jets:" << bJetsFound;RET;
    
      bool jetSelection = false;
      if ( noBtag.size() == 1 && bJetsFound == 1 ) {
	cout << "two jets found (exactly one of which is a b-jet)!"; RET;
	jetSelection = true;
      } else {
	if ( bJetsFound != 1) {N1;RET;}
	if ( noBtag.size() !=1) {N2;RET;}
      }
      if (!jetSelection) {cout << "jet selection failed";RET;N1;}
    else {cout << "jet selection passed";RET;}

    // not used!
    // missing energy cut for parton and particle level analysis
    bool metSelection = false;
    if (misMom.met() < 30*GeV) {
      cout << "NOT USED: metSelection failed (met (GeV) = " << misMom.met() << ")"; RET;
    }
    else {
      metSelection = true;
      cout << "NOT USED: metSelection pased (met (GeV) = " << misMom.met() << ")"; RET;
    }
    
    
    /// Select leptons
    bool leptonSelection = false;
    Particles leptons;
    
    cout << "muons.size()=" << muons.size(); RET;
    cout << "electrons.size()=" << electrons.size(); RET;

    // collect electrons and muons into single vector of leptons
    for (const Particle& e : electrons){
      leptons.push_back(e);
      cout << "pushed back one electron to 'leptons' vector";RET;
    }
    for (const Particle& m : muons){
      leptons.push_back(m);
      cout << "pushed back one muon to 'leptons' vector";RET;
    }
    
    // print out lepton eta and pT info for all leptons in vector
    cout << "lepton info for " << leptons.size() << " lepton(s):"; RET;
    for (const Particle& l : leptons){
      cout << "lep: eta=" << l.eta() << " pT=" << l.pT();RET;
    }
    
    //perform eta cuts on leptons:
    Cut e_mu_eta_Cut = Cuts::abseta < 2.5;
    ifilter_select(leptons, e_mu_eta_Cut);
    // print out lepton eta and pT info for all remaining leptons in vector
    cout << "post eta cut filter lepton info:";RET;
    for (const Particle& l : leptons){
      cout << "lep: eta=" << l.eta() << " pT=" << l.pT();RET;
    }
    
    //perform pT cuts on leptons
    Cut e_mu_pT_Cut = Cuts::pT > 25*GeV;
    ifilter_select(leptons, e_mu_pT_Cut);
    // print out lepton eta and pT info for all remaining leptons in vector
    cout << "post pT cut filter leptons.size()=" << leptons.size(); RET;
    for (const Particle& l : leptons){
      cout << "lep: eta=" << l.eta() << " pT=" << l.pT();RET;
    }
    
    if ( leptons.size() == 1 ){
      leptonSelection = true;
    }
    
    if (!leptonSelection) {cout << "lepton selection failed";RET;}
    else {cout << "lepon selection passed";RET;}
    
    
    
    // discared event if selected jets and single lepton are not isolated
    bool isolationSelection = false;
    if (leptonSelection && jets.size() > 0){
      bool jl_tooClose= false;
      for (const Jet& j : jets){
	if (deltaR(j, leptons[0]) < 0.4) {
	  jl_tooClose = true ;
	  cout << "jet found to be within deltaR<0.4 of lepton";RET;
	} else {
	  cout << "jet well separated from lepton";RET;
	}
      }
      isolationSelection = !jl_tooClose;
      if ( isolationSelection ) {cout << "isolation requirement met!";RET;}
      else {cout << "isolation requirement failed - a jet was too close to the single particle-level lepton!";RET;}
    }
    
      
    // calculate mass of lepton + b-jet system
    bool massljSelection = false;
    if (leptons.size() == 1 && bJetsFound == 1){
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
    if (!massljSelection) {cout << "jet_lepton mass selection failed";RET;}
    else {cout << "jet_lepton mass selection passed";RET;}
    
      
    // construct pseudo-W

    // TODO: what if more than one W reconstructed?!
      
    FourMomentum w4p;
    bool WbosonFound = false;
    //if (leptonSelection && jetSelection && isolationSelection && massljSelection){
            
    if ( w_el.bosons().size() != 0 ){
      if (w_el.bosons().size() != 1 ) {
	cout << "more than one electron W bosons found - veto";RET;
	N8;RET;vetoEvent;
      }
	  
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
      if ( WbosonFound ) { // not sure which W to use in pseudo-top definition
	cout << "single electron W bosons already found! Not sure which should be used - veto";RET;
	N8;RET;vetoEvent;
      } 
      else WbosonFound = true;
	  
      if (w_mu.bosons().size() != 1 ) {
	cout << "more than one muon W bosons found - veto";RET;
	N8;RET;vetoEvent;
      }
      else {	  
	const Particles W_mu= w_mu.bosons();
	
	//debug info /// START /////////////////////////////////////
	cout << "w_mu.bosons().size()=" << w_mu.bosons().size(); RET;
	cout << "W_mu.size()=" << W_mu.size() << ", " ;
	cout << "W_mu.constituentLepton().pid()=" << w_mu.constituentLepton().pid() << ", ";
	//debug info /// END ////////////////////////////////////////
	
	w4p = W_mu[0];
	cout << "W reconstructed from muon (naming it w4p)!" << endl;
	cout << "w4p.pT=" << w4p.pT() << endl;
	cout << "w4p.absrap=" << w4p.absrap() << endl;	
      }
    }

    
    // construct pseudo-top
    bool pseudoTopFound = false;
    FourMomentum pseudoTop4p;
    
    if ( bJetsFound == 1 && WbosonFound )
      {
	pseudoTop4p = add(bJet, w4p);
	pseudoTopFound = true;
	cout << "pseudo top created!";RET;
	cout << "pseudoTop4p.pT()=" << pseudoTop4p.pT(); RET;
	cout << "pseudoTop4p.absrap()=" << pseudoTop4p.absrap(); RET;
      }
    else {
      cout << "failed to create pseudo top due to ";
      
      if (WbosonFound && bJetsFound != 1){
	cout << "wrong number of b-jets (" << bJetsFound << ") despite W being reconstructed from a lepton and met";RET;
      }
      else if (!WbosonFound && bJetsFound == 1){
	cout << "correct number of b-jets (1) but no W reconstructed"; RET;
      }
      else if (!WbosonFound && bJetsFound != 1){
	cout << "wrong number of b-jets (" << bJetsFound << ") and no W reconstructed"; RET;
      } else {
	cout << "you screwed up - veto";RET;N8;RET; vetoEvent;
      }
    }
    
    
    bool makesSelection = false;
    if ( jetSelection && leptonSelection && massljSelection && isolationSelection ){
      makesSelection = true;
    }
    else { cout << "selection criteria failed";RET;}
    
    
    if (makesSelection && !WbosonFound && !pseudoTopFound){
      //if ( metSelection ) {N2;RET;} else {N3;RET;}
    }
    
    if ( makesSelection && WbosonFound && !pseudoTopFound){
      cout << "Rivet finds W presumed to be from the top but not the bjet!";RET;N8;RET;
    }
    
    if ( makesSelection && !WbosonFound && pseudoTopFound){
      cout << "This cannot happen" << endl;
      N8;RET;
    }
    
    
    
    
    // for fiducial selections fill fid cross section counters and particle-level plots
    if ( jetSelection && leptonSelection && massljSelection &&
	 isolationSelection && WbosonFound && pseudoTopFound ){
      
      
      
      _c_sumw_fid->fill(event.weight());
      
      if ( leptons[0].charge() > 0 ) { // top quark 
	
	N6;RET;  	  
	
	_c_sumw_fid_tq->fill( event.weight() ); // counter for normalising hostograms
	_c_Xsec_fid_tq->fill(event.weight()) ;

	_h_AbsPtclDiffXsecTPt->fill( pseudoTop4p.pT(),     event.weight() ); //(1,1,1)  
	_h_AbsPtclDiffXsecTY->fill(  pseudoTop4p.absrap(), event.weight() ); //(1,1,3)
	_h_NrmPtclDiffXsecTPt->fill( pseudoTop4p.pT(),     event.weight() ); //(2,1,1)
	_h_NrmPtclDiffXsecTY->fill(  pseudoTop4p.absrap(), event.weight() ); //(2,1,3)
	
	_h_AbsPtclDiffXsecTjPt->fill( noBtag[0].pT(),     event.weight() ); //(5,1,1);
	_h_AbsPtclDiffXsecTjY->fill(  noBtag[0].absrap(), event.weight() ); //(5,1,3);
	_h_NrmPtclDiffXsecTjPt->fill( noBtag[0].pT(),     event.weight() ); //(6,1,1);
	_h_NrmPtclDiffXsecTjY->fill(  noBtag[0].absrap(), event.weight() ); //(6,1,3);
	

      }
      
      else if ( leptons[0].charge() < 0 ) { // top antiquark
	
	N7;RET;
	
	_c_sumw_fid_tbarq->fill( event.weight() ); // counter for normalising hostograms
	_c_Xsec_fid_tbarq->fill(event.weight()) ;
	
	_h_AbsPtclDiffXsecTbarPt->fill( pseudoTop4p.pT(),     event.weight() );  //(1,1,2)
	_h_AbsPtclDiffXsecTbarY->fill(  pseudoTop4p.absrap(), event.weight() );  //(1,1,4)
	_h_NrmPtclDiffXsecTbarPt->fill( pseudoTop4p.pT(),     event.weight() );  //(2,1,2)
	_h_NrmPtclDiffXsecTbarY->fill(  pseudoTop4p.absrap(), event.weight() );  //(2,1,4)
	
	_h_AbsPtclDiffXsecTbarjPt->fill( noBtag[0].pT(),     event.weight() ); //(5,1,2);
	_h_AbsPtclDiffXsecTbarjY->fill(  noBtag[0].absrap(), event.weight() ); //(5,1,4);
	_h_NrmPtclDiffXsecTbarjPt->fill( noBtag[0].pT(),     event.weight() ); //(6,1,2);
	_h_NrmPtclDiffXsecTbarjY->fill(  noBtag[0].absrap(), event.weight() ); //(6,1,4);


      } else {cout << "well what is it now!?";RET;N8;RET;}
      
    }
    
    
    
    
    /// Parton-level analysis  ///////////////////////////////////////////////
    
    
    if (allTops.size() == 1){ // do partonic analysis
      
      
      if (allTops[0].charge() > 0){
	
	_c_sumw_prtn_tq->fill( event.weight() );
	
	cout << "add to plots: parton-level t!";RET;	  
	
	_h_AbsPrtnDiffXsecTPt->fill( allTops[0].pT()    , event.weight()) ;
	_h_AbsPrtnDiffXsecTY->fill(  allTops[0].absrap(), event.weight()) ;
	
	_h_NrmPrtnDiffXsecTPt->fill( allTops[0].pT(),     event.weight()) ;
	_h_NrmPrtnDiffXsecTY->fill(  allTops[0].absrap(), event.weight()) ;
	
      }
      else {
	
	_c_sumw_prtn_tbarq->fill( event.weight() );
	
	cout << "add to plots: parton-level t~";RET;
	
	_h_AbsPrtnDiffXsecTbarPt->fill( allTops[0].pT()    , event.weight()) ;
	_h_AbsPrtnDiffXsecTbarY->fill(  allTops[0].absrap(), event.weight()) ;
	
	_h_NrmPrtnDiffXsecTbarPt->fill( allTops[0].pT(),     event.weight()) ;
	_h_NrmPrtnDiffXsecTbarY->fill(  allTops[0].absrap(), event.weight()) ;
	
      }
    }
    
    
    
    // debug START /////////////////////////////////////////////
    if ( WbosonFound ){
      if ( makesSelection ){
	cout << "report: finds W and makes selection!";RET;
      }
      else {
	cout << "report: finds W but fails selection!";RET;
      }
    }
    else {
      if ( makesSelection ){
	cout << "report: doesn't find W but makes selection!";RET;
      }
      else {
	cout << "report: doesn't find W - doesn't makes selection!";RET;
      }
    }
    // debug end ///////////////////////////////////////////////
    
    RET;
    }
    
    
    /// Normalise histograms etc., after the run
    void finalize() {
      
      RET;RET;

      /// compute cross sections

      int nrm = 1000; // this is a corrective factor for a yet to be corrected
                      // mistake in the reference yoda file for ALL the normalised
                      // differential cross sections.

      
      // tq particle-level (fiducial)
      scale( _h_AbsPtclDiffXsecTPt, crossSection()/femtobarn/sumOfWeights() );
      scale( _h_AbsPtclDiffXsecTY,  crossSection()/sumOfWeights() );
      scale( _h_NrmPtclDiffXsecTPt, nrm / _c_Xsec_fid_tq->val() );
      scale( _h_NrmPtclDiffXsecTY,  nrm / _c_Xsec_fid_tq->val() );

      scale( _h_AbsPtclDiffXsecTjPt, crossSection()/femtobarn/sumOfWeights() );
      scale( _h_AbsPtclDiffXsecTjY,  crossSection()/sumOfWeights() );
      scale( _h_NrmPtclDiffXsecTjPt, nrm / _c_Xsec_fid_tq->val() );
      scale( _h_NrmPtclDiffXsecTjY,  nrm / _c_Xsec_fid_tq->val() );
            
          
      // tbarq particl-level (fiducial)
      scale( _h_AbsPtclDiffXsecTbarPt, crossSection()/femtobarn/sumOfWeights() );
      scale( _h_AbsPtclDiffXsecTbarY,  crossSection() / sumOfWeights() );
      scale( _h_NrmPtclDiffXsecTbarPt, nrm / _c_Xsec_fid_tbarq->val() );
      scale( _h_NrmPtclDiffXsecTbarY,  nrm / _c_Xsec_fid_tbarq->val() );

      scale( _h_AbsPtclDiffXsecTbarjPt, crossSection()/femtobarn/sumOfWeights() );
      scale( _h_AbsPtclDiffXsecTbarjY,  crossSection() / sumOfWeights() );
      scale( _h_NrmPtclDiffXsecTbarjPt, nrm / _c_Xsec_fid_tbarq->val() );
      scale( _h_NrmPtclDiffXsecTbarjY,  nrm / _c_Xsec_fid_tbarq->val() );


      // tq parton-level
      scale(_h_AbsPrtnDiffXsecTPt, crossSection()/femtobarn/sumOfWeights() );
      scale(_h_AbsPrtnDiffXsecTY,  crossSection() / sumOfWeights() );
      scale(_h_NrmPrtnDiffXsecTPt, nrm / _c_sumw_prtn_tq->val() );
      scale(_h_NrmPrtnDiffXsecTY,  nrm / _c_sumw_prtn_tq->val() );


      // tbarq parton-level
      scale(_h_AbsPrtnDiffXsecTbarPt, crossSection()/femtobarn/sumOfWeights() );
      scale(_h_AbsPrtnDiffXsecTbarY,  crossSection() / sumOfWeights() );
      scale(_h_NrmPrtnDiffXsecTbarPt, nrm / _c_sumw_prtn_tbarq->val() );
      scale(_h_NrmPrtnDiffXsecTbarY,  nrm / _c_sumw_prtn_tbarq->val() );
 
      
      
      
      // Print info (just for me)
      ////////////////////////////////////////////////////
      RET;
      cout << "sum of weights = " << sumOfWeights(); RET;
      cout << "cross section = " << crossSection(); RET;
      RET;RET;
      /////////////////////////////////////////////////////
      
           
    }
    
    //@}
    
    
    /// @name Histograms
    //@{
    
    CounterPtr _c_fid_t;
    CounterPtr _c_fid_tbar;
    CounterPtr _c_Xsec_fid_tq;
    CounterPtr _c_Xsec_fid_tbarq;
    CounterPtr _c_sumw_singleTop;
    CounterPtr _c_sumw_fid;
    CounterPtr _c_sumw_fid_tq;
    CounterPtr _c_sumw_fid_tbarq;
    CounterPtr _c_sumw_prtn_tq;
    CounterPtr _c_sumw_prtn_tbarq;
    CounterPtr _c_N_tot_tq;
    CounterPtr _c_N_tot_tbarq;
    
    Histo1DPtr _h_AbsPtclDiffXsecTPt;
    Histo1DPtr _h_AbsPtclDiffXsecTY;
    Histo1DPtr _h_NrmPtclDiffXsecTPt;
    Histo1DPtr _h_NrmPtclDiffXsecTY;
    Histo1DPtr _h_AbsPtclDiffXsecTbarPt;
    Histo1DPtr _h_AbsPtclDiffXsecTbarY;
    Histo1DPtr _h_NrmPtclDiffXsecTbarPt;
    Histo1DPtr _h_NrmPtclDiffXsecTbarY;
    Histo1DPtr _h_AbsPtclDiffXsecTjPt;
    Histo1DPtr _h_AbsPtclDiffXsecTjY;
    Histo1DPtr _h_NrmPtclDiffXsecTjPt;
    Histo1DPtr _h_NrmPtclDiffXsecTjY;
    Histo1DPtr _h_AbsPtclDiffXsecTbarjPt;
    Histo1DPtr _h_AbsPtclDiffXsecTbarjY;
    Histo1DPtr _h_NrmPtclDiffXsecTbarjPt;
    Histo1DPtr _h_NrmPtclDiffXsecTbarjY;
    Histo1DPtr _h_AbsPrtnDiffXsecTPt;
    Histo1DPtr _h_AbsPrtnDiffXsecTY;
    Histo1DPtr _h_NrmPrtnDiffXsecTPt;
    Histo1DPtr _h_NrmPrtnDiffXsecTY;
    Histo1DPtr _h_AbsPrtnDiffXsecTbarPt;
    Histo1DPtr _h_AbsPrtnDiffXsecTbarY;
    Histo1DPtr _h_NrmPrtnDiffXsecTbarPt;
    Histo1DPtr _h_NrmPrtnDiffXsecTbarY;
    
    //@}
    
    std::unique_ptr<HepMC::IO_GenEvent> _hepmcout;
    
    };


// The hook for the plugin system
DECLARE_RIVET_PLUGIN(ATLAS_2017_I1512776);

}

