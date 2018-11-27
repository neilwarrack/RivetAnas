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


// ============ Debug Toolkit ============= //
#define N1 do { std::cout << "1 "; } while(0)
#define N2 do { std::cout << "2 "; } while(0)
#define N3 do { std::cout << "3 "; } while(0)
#define N4 do { std::cout << "4 "; } while(0)
#define N5 do { std::cout << "5 "; } while(0)
#define N6 do { std::cout << "6 "; } while(0)
#define N7 do { std::cout << "7 "; } while(0)
#define N8 do { std::cout << "8 "; } while(0)
#define N9 do { std::cout << "9 "; } while(0)
#define NL do { std::cout << endl; } while(0)
// ========================================= // 


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
      PromptFinalState promptElectrons(electron_fs, true, true); //< "true, true" accepts as 'prompt' particles from intermediate taus and muons, respectively.
      DressedLeptons dressedElectrons(photons, promptElectrons, 0.1, e_muCuts, true);
      // NB: the above bool ("true") includes 'decay photons' in dressing...
      // NB: 'Decay photons' are defined in depricated manner, using:
      // NB: fromDecay() const { return fromHadron() || fromPromptTau(); } 
      declare(dressedElectrons, "DressedElectrons");


      // particle-level muons (similar to electrons above)
      IdentifiedFinalState muon_fs;
      muon_fs.acceptIdPair(PID::MUON);
      PromptFinalState promptMuons(muon_fs, true, true); 
      DressedLeptons dressedMuons(photons, promptMuons, 0.1, e_muCuts, true); 
      declare(dressedMuons, "DressedMuons");


      // Particle-level W (with electron decay)
      WFinder w_el(fs, e_muCuts,PID::ELECTRON, 35*GeV, 8000*GeV, 30*GeV, 0.1,
		   WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERALL,
		   WFinder::TRACK,  WFinder::TRANSMASS);
      // notes: ET_miss>30GeV (see: https://arxiv.org/pdf/1502.05923.pdf)
      // notes: mass window > 35GeV (< transverse mass) (see: https://arxiv.org/pdf/1502.05923.pdf)
      declare(w_el, "W_Electron");


      // Particle-level W (with muon decay)
      WFinder w_mu(fs, e_muCuts, PID::MUON, 35*GeV, 8000*GeV, 30*GeV, 0.1,
		     WFinder::PROMPTCHLEPTONS, WFinder::CLUSTERALL,
		     WFinder::TRACK, WFinder::TRANSMASS);
      declare(w_mu, "W_Muon");
      
      
      //the selected lepton (with it's photons) is not to be included in jet clustering
      VetoedFinalState jet_input(fs);
      jet_input.addVetoOnThisFinalState(dressedElectrons);
      jet_input.addVetoOnThisFinalState(dressedMuons);
      declare(FastJets(jet_input, FastJets::ANTIKT, 0.4), "Jets");

      
      // for parton level analysis
      declare(PartonicTops(PartonicTops::E_MU), "LeptonicTops") ;//includes from taus
      declare(PartonicTops(PartonicTops::ALL ), "AllPartonicTops") ;

      
      // book counters
      _c_sumw_fid =         bookCounter("sumw_fid_tq_tbarq");
      _c_sumw_fid_tq =      bookCounter("sumw_fid_tq");
      _c_sumw_fid_tbarq =   bookCounter("sumw_fid_tbarq");
      _c_sumw_prtn_tq =     bookCounter("sumw_parton_tq");
      _c_sumw_prtn_tbarq =  bookCounter("sumw_parton_tbarq");
      _c_Xsec_fid_tq =      bookCounter("Xsec_fid_tq");
      _c_Xsec_fid_tbarq =   bookCounter("Xsec_fid_tbarq");

      
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
      

    }
    


    
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      
      // perform fast veto if not a single top event
      // (NB: do not give this analysis s-channel single top events as input!)

      const Particles allTops = apply<ParticleFinder>(event, "AllPartonicTops").particlesByPt();
      
      if ( allTops.size() != 1 ){
	cout << "not a single top event! VETOING!!" << endl;
	N8; NL;
	vetoEvent;
      }

      
      // Find leptons, partonic tops, missing energy and jets in event
      const Particles& electrons = apply<DressedLeptons>(event, "DressedElectrons").particlesByPt();
      const Particles& muons = apply<DressedLeptons>(event, "DressedMuons").particlesByPt();
      const WFinder& w_el = apply<WFinder>(event, "W_Electron");
      const WFinder& w_mu = apply<WFinder>(event, "W_Muon");
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt();
      
 
      /// require exactly one bJet and one non-bJet jet for particle level results
      // apply particle-level jet cuts
      Cut jetCuts = (Cuts::abseta < 4.5 && Cuts::pT > 30*GeV);
      cout << "sizeOfJetsPreFilter=" << jets.size() << endl;
      ifilter_select(jets, jetCuts);
      cout << "sizeOfJetsPostFilter=" << jets.size() << endl;

      // find and bJets
      Jet bJet;
      Jets noBtag;
      int bJetsFound = 0;
      for (const Jet& j : jets){
	
	//if ( j.bTagged()) { cout << "b-tagged jet found" << endl;}
	//if ( j.bTagged(Cuts::pT > 5*GeV)) {cout << "b-tagged jet found with b>5GeV found" << endl;}
	if ( j.bTagged(Cuts::pT > 5*GeV) && j.abseta() < 2.5) {
	  cout << "b-tagged jet found with b.pT>5GeV inside eta<2.5" << endl;
	  bJetsFound++;
	  bJet = j;
	} else {
	  noBtag.push_back(j);
	}
      }

      cout << "Total number of particle-level b-jets:" << bJetsFound << endl;
    
      bool jetSelection = false;
      if ( noBtag.size() == 1 && bJetsFound == 1 ) {
	cout << "two jets found (exactly one of which is a b-jet)!" << endl;
	jetSelection = true;
      } else {
	if ( bJetsFound != 1) {N1;NL;}   // THIS PRINTS A COUNTABLE CODE (1)
	if ( noBtag.size() !=1) {N2;NL;} // THIS PRINTS A COUNTABLE CODE (2)
      }
      if (!jetSelection) {
	cout << "jet selection failed" << endl;
      }
      else {
	cout << "jet selection passed" << endl;
      }
    
    
    
      /// Select leptons
      bool leptonSelection = false;
      Particles leptons;
    
      cout << "muons.size()=" << muons.size() << endl;
      cout << "electrons.size()=" << electrons.size() << endl;

      // collect electrons and muons into single vector of leptons
      for (const Particle& e : electrons){
	leptons.push_back(e);
      }
      for (const Particle& m : muons){
	leptons.push_back(m);
      }
    
      // print out lepton eta and pT info for all leptons in vector
       for (const Particle& l : leptons){
	 cout << "lep: eta=" << l.eta() << " pT=" << l.pT() << endl;
       }
       
       if ( leptons.size() == 1 ){
	 leptonSelection = true;
       } else {
	 N3;NL; // THIS PRINTS A COUNTABLE CODE (3)
       }
       
       if (!leptonSelection) {cout << "lepton selection failed" << endl;}
       else {cout << "lepon selection passed" << endl;}
    
    
       // discared event if selected jets and single lepton are not isolated
       bool isolationSelection = false;
       if (leptonSelection && jets.size() > 0){
	 bool jl_tooClose= false;
	 for (const Jet& j : jets){
	   if (deltaR(j, leptons[0]) < 0.4) {
	     jl_tooClose = true ;
	     cout << "jet found to be within deltaR<0.4 of lepton" << endl;
	   } else {
	     cout << "jet well separated from lepton" << endl;
	   }
	 }
	 isolationSelection = !jl_tooClose;
	 if ( isolationSelection ) {cout << "isolation requirement met!" << endl;}
	 else {
	   cout << "isolation requirement failed - a jet was too close to the single "
		<< "particle-level lepton!" << endl;

	   N4;NL; // THIS PRINTS A COUNTABLE CODE (5)
	   
	 }
       }
       
       
       // calculate mass of lepton + b-jet system
       bool massljSelection = false;
       if (leptons.size() == 1 && bJetsFound == 1){
	 FourMomentum lepton4p = leptons[0];
	 FourMomentum bJet4p = bJet;
	 FourMomentum lb4p = add(lepton4p, bJet4p);
	 if (lb4p.mass() < 160*GeV){
	   massljSelection = true;
	 }
       }
       if (!massljSelection) {
	 cout << "jet_lepton mass selection failed" << endl;
	 N5;NL; // THIS PRINTS A COUNTABLE CODE (5)
       }
       else {cout << "jet_lepton mass selection passed" << endl;}
       
      
       // construct pseudo-W

       // TODO: what if more than one W reconstructed?!
      
       FourMomentum w4p;
       bool WbosonFound = false;
       
       if ( w_el.bosons().size() != 0 ){
	 if (w_el.bosons().size() != 1 ) {
	   cout << "more than one electron W bosons found - veto" << endl;
	   N8;NL; // THIS PRINTS OUT THE 'ERROR' CODE (8)
	   vetoEvent;
	 }
	 
	 WbosonFound = true;
	 const Particles W_el= w_el.bosons();
	 
	 w4p = W_el[0];
       }
       
       if ( w_mu.bosons().size() != 0 ){
	 if ( WbosonFound ) { // not sure which W to use in pseudo-top definition
	   cout << "single electron W bosons already found! Not sure which should "
		<< "be used - veto event!" << endl;
	   N8;NL; // THIS PRINTS OUT THE 'ERROR' CODE (8)
	   vetoEvent;
	 } 
	 else WbosonFound = true;
	  
	 if (w_mu.bosons().size() != 1 ) {
	   cout << "more than one muon W bosons found - veto" << endl;
	   N8;NL; // THIS PRINTS OUT THE 'ERROR' CODE (8)
	   vetoEvent;
	 }
	 else {	  
	   const Particles W_mu = w_mu.bosons();
	   w4p = W_mu[0];
	 }
       }
       
       if (!WbosonFound) {
	 N6;NL; // THIS PRINTS A COUNTABLE CODE (6)
       }
       
    
       // construct pseudo-top
       bool pseudoTopFound = false;
       FourMomentum pseudoTop4p;
       
       if ( bJetsFound == 1 && WbosonFound ) {
	 pseudoTop4p = add(bJet, w4p);
	 pseudoTopFound = true;
       }
       else {
	 cout << "failed to create pseudo top" << endl;
       }
       
       
       
       bool fidSelection = false;
       if ( jetSelection && leptonSelection && massljSelection && isolationSelection ){
	 fidSelection = true;
       }
       
       
    
    
       // for fiducial selections fill cross section counters and particle-level plots:
       if ( fidSelection && WbosonFound && pseudoTopFound ){
	 
	 _c_sumw_fid->fill(event.weight());
	 
	 if ( leptons[0].charge() > 0 ) { // top quark 
	   
	   _c_sumw_fid_tq->fill( event.weight() );
	   _c_Xsec_fid_tq->fill( event.weight() );
	   
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
	   
	   _c_sumw_fid_tbarq->fill( event.weight() );
	   _c_Xsec_fid_tbarq->fill( event.weight() );
	   
	   _h_AbsPtclDiffXsecTbarPt->fill( pseudoTop4p.pT(),     event.weight() ); //(1,1,2)
	   _h_AbsPtclDiffXsecTbarY->fill(  pseudoTop4p.absrap(), event.weight() ); //(1,1,4)
	   _h_NrmPtclDiffXsecTbarPt->fill( pseudoTop4p.pT(),     event.weight() ); //(2,1,2)
	   _h_NrmPtclDiffXsecTbarY->fill(  pseudoTop4p.absrap(), event.weight() ); //(2,1,4)
	   
	   _h_AbsPtclDiffXsecTbarjPt->fill( noBtag[0].pT(),     event.weight() ); //(5,1,2);
	   _h_AbsPtclDiffXsecTbarjY->fill(  noBtag[0].absrap(), event.weight() ); //(5,1,4);
	   _h_NrmPtclDiffXsecTbarjPt->fill( noBtag[0].pT(),     event.weight() ); //(6,1,2);
	   _h_NrmPtclDiffXsecTbarjY->fill(  noBtag[0].absrap(), event.weight() ); //(6,1,4);
	   
	   
	 }
	 
       }
       
       
       
       
       /// Parton-level ///
       
       
       if (allTops.size() == 1){ // fill partonic histos
	 
	 if (allTops[0].charge() > 0){ // top quark
	   
	   _c_sumw_prtn_tq->fill( event.weight() );
	   
	   _h_AbsPrtnDiffXsecTPt->fill( allTops[0].pT()    , event.weight()) ;
	   _h_AbsPrtnDiffXsecTY->fill(  allTops[0].absrap(), event.weight()) ;
	   
	   _h_NrmPrtnDiffXsecTPt->fill( allTops[0].pT(),     event.weight()) ;
	   _h_NrmPrtnDiffXsecTY->fill(  allTops[0].absrap(), event.weight()) ;
	   
	 }
	 else { // top antiquark
	   
	   _c_sumw_prtn_tbarq->fill( event.weight() );
	   
	   _h_AbsPrtnDiffXsecTbarPt->fill( allTops[0].pT()    , event.weight()) ;
	   _h_AbsPrtnDiffXsecTbarY->fill(  allTops[0].absrap(), event.weight()) ;
	   
	   _h_NrmPrtnDiffXsecTbarPt->fill( allTops[0].pT(),     event.weight()) ;
	   _h_NrmPrtnDiffXsecTbarY->fill(  allTops[0].absrap(), event.weight()) ;
	   
	 }
       }
    }
  
    
    /// Normalise histograms etc., after the run
    void finalize() {
      

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
 
      
      
      

      ////////////////////////////////////////////////////
      cout << endl;
      cout << "sum of weights = " << sumOfWeights() << endl;
      cout << "cross section = " << crossSection() << endl;
      cout << endl;
      /////////////////////////////////////////////////////
      
           
    }
    
    //@}
    
    
    /// @name Histograms
    //@{
    
    CounterPtr _c_Xsec_fid_tq;
    CounterPtr _c_Xsec_fid_tbarq;
    CounterPtr _c_sumw_fid;
    CounterPtr _c_sumw_fid_tq;
    CounterPtr _c_sumw_fid_tbarq;
    CounterPtr _c_sumw_prtn_tq;
    CounterPtr _c_sumw_prtn_tbarq;
    
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
    
        
    };


// The hook for the plugin system
DECLARE_RIVET_PLUGIN(ATLAS_2017_I1512776);

}

