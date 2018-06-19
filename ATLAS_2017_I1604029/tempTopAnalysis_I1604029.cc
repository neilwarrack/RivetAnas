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

  
  class tempTopAnalysis_I1604029 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(tempTopAnalysis_I1604029);

   
    /// @name Analysis methods
    //@{
   
    /// Book histograms and initialise projections before the run
    void init() {

      //      if (fuzzyEquals(sqrtS(), 7*TeV)) {
      //	...
      //      } else {
      //	...
      //      }
 

      FinalState fs;

      // Particle-level electron and muon cuts
      Cut e_muCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV ;
      Cut e_muForDressingCuts = Cuts::abseta < 2.7 && Cuts::pT > 10*GeV ;
      Cut hePhotonCuts = Cuts::abseta < 2.37 && Cuts::Et > 15*GeV ;

      // Photon projection to dress leptons for lepton & jet definitions
      //NBool: PromptFinalState(const Cut& c, bool accepttaudecays=false, bool acceptmudecays=false);
      PromptFinalState hardParticles(true, true); //selects particle with no hadrons in ancestor line
      IdentifiedFinalState photons(hardParticles);
      photons.acceptIdPair(PID::PHOTON);
      
      /// Leptons (Propper)
      IdentifiedFinalState leptons_fs({PID::ELECTRON, PID::MUON, -PID::ELECTRON, -PID::MUON}, e_muForDressingCuts);
      DressedLeptons dressedLeptons(photons, leptons_fs, 0.1, e_muCuts, false);//< for testing cuts (NB: not much impact))
      //DressedLeptons dressedLeptons(photons, leptons_fs, 0.1, e_muCuts, true);//< for testing cuts
      declare(dressedLeptons, "DressedLeptons");

      /// Leptons (Testing)
      IdentifiedFinalState allLeptons_fs(fs, {PID::ELECTRON, PID::MUON, -PID::ELECTRON, -PID::MUON});//< for testing cuts
      //IdentifiedFinalState allLeptons_fs(e_muForDressingCuts);//< for testing cuts
      //allLeptons_fs.acceptIdPair(PID::MUON);//< for testing cuts
      //allLeptons_fs.acceptIdPair(PID::ELECTRON);//< for testing cuts
      DressedLeptons dressedLeptonsNoCuts(photons, allLeptons_fs, 0.1, Cuts::open(), false);//< for testing cuts
      //*DressedLeptons::DressedLeptons(const FinalState& photons, const FinalState& bareleptons, double dRmax, const Cut& cut, bool useDecayPhotons)
      declare(dressedLeptonsNoCuts, "DressedLeptonsNoCuts");

      // projections for jets
      IdentifiedFinalState neutrinos(fs);
      neutrinos.acceptNeutrinos();
      IdentifiedFinalState muons(e_muForDressingCuts, {PID::MUON, -PID::MUON});
      //IdentifiedFinalState muon_fs(PID::MUON, e_muForDressingCuts); // needs anti-muons!
      DressedLeptons dressedMuons(photons, muons, 0.1, e_muCuts, false);
      //DressedLeptons dressedMuons(photons, muon_fs, 0.1, e_muCuts);
      VetoedFinalState jet_input(fs);
      jet_input.addVetoOnThisFinalState(muons);
      //jet_input.addVetoOnThisFinalState(dressedMuons);
      jet_input.addVetoOnThisFinalState(neutrinos);
      //declare(jet_input, "Jet_Input");
      declare(FastJets(jet_input, FastJets::ANTIKT, 0.4), "Jets");

      // Projection for high energy photons
      PromptFinalState heParticles(hePhotonCuts, true, true);
      IdentifiedFinalState hePhotons(heParticles);
      hePhotons.acceptIdPair(PID::PHOTON);
      declare(hePhotons, "highEnergyPhotons");


      
      // Projection for Parton Level Top Quarks
      //declare(PartonicTops(PartonicTops::E_MU, e_muCuts), "LeptonicTops") ;
      //declare(PartonicTops(PartonicTops::ALL, e_muCuts), "AllPartonicTops") ;
      //declare(PartonicTops(PartonicTops::E_MU), "LeptonicTops") ;
      //declare(PartonicTops(PartonicTops::ALL ), "AllPartonicTops") ;
      
      // Book counters
      _c_fidXsec = bookCounter("fidTotXsecttbargamma");
      //_c_fid_tbar = bookCounter("fidTotXsectbarq");

      // book histograms
 
      //_h_PhotonPt = bookHisto1D(1,1,1);
      //_h_PhotonAbsEta = bookHisto1D(1,1,2);

            
      // debug hepmc print-out
      _hepmcout.reset(new HepMC::IO_GenEvent("inspectevents.hepmc"));    

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      BEGIN;
          
      
      // Find leptons, jets and high energy photons
      const Particles& leptons = apply<DressedLeptons>(event, "DressedLeptons").particlesByPt();
      const Particles& leptonsNoCuts = apply<DressedLeptons>(event, "DressedLeptonsNoCuts").particlesByPt();
      const Particles& hePhotons = apply<IdentifiedFinalState>(event, "highEnergyPhotons").particlesByPt();
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt();

      
      /// object removal
      Particles electrons, selectedMuons, selectedElectrons;
      Jets selectedJetsEl, selectedJets;
      Jet bJet;
      bool tooClose = false ;
      double dR_min = 0.0 ;
      int rejectN=0, bJetN=0 ;
	
      // jet selection
      Cut jetCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV ;
      cout << "Jets=" << jets.size()<< " "; RET;
      ifilter_select(jets, jetCuts);
      cout << "JetsInFiducialPhaseSpace=" << jets.size()<< " "; RET;
      
      // veto if not enough jets
      if ( jets.size() < 4) {
	cout << "veto event: not enough jets to proceed";RET;
	cout << "VETOCODES:";RET;
	N3; RET; //NB: N3 is used twice!
	vetoEvent;
      }

      // testing lepton cuts:****************************************
      cout << "testing lepton cuts:****************************************"; RET;
      Particles electronsNoCuts;
      Particles muonsNoCuts;
      Particles freshCutLeptons;
      Particles selectedCutLeptons;
      if (leptonsNoCuts.empty() ) {
	cout << "no dressed (noCut)leptons"; RET; vetoEvent; 
      } else { //if there ARE leptons with no cuts
	cout << "searching through DressedLeptonsNoCuts...";RET;
	for (const Particle& l : leptonsNoCuts){
	  if (l.abspid() == PID::MUON){
	    cout << "found:muon";
	    if (l.abseta() > 2.5) {cout << "...failed eta cut.";RET;}
	    else if (l.pT() < 25*GeV ) {cout << "...failed pt cut.";RET;} 
	    muonsNoCuts.push_back(l);
	    RET;
	  }
	  else if (l.abspid() == PID::ELECTRON){
	    cout << "found:electron";
	    if (l.abseta() > 2.5) {cout << "...failed eta cut."; RET;}
	    else if (l.pT() < 25*GeV ) {cout << "...failed pt cut."; RET;}
	    electronsNoCuts.push_back(l);
	    RET;
	  }
	  else {cout<< "lepton wierdness!"; RET; N1; RET;}
	}
      
	cout << "electronsNoCuts.size()=" << electronsNoCuts.size(); RET;
	cout << "muonsNoCuts.size()=" << muonsNoCuts.size(); RET;
       
	for (const Particle& l : electronsNoCuts){freshCutLeptons.push_back(l);}
	for (const Particle& l : muonsNoCuts){freshCutLeptons.push_back(l);}
	int lepCtr = 0;
	for (const Particle& l : freshCutLeptons){
	  lepCtr++;
	  cout << "considering uncut dressed lepton No. " << lepCtr; RET;
	  if (l.pid == PID::ELECTRON || l.pid == -PID::ELECTRON) {
	    if (l.pT() < 25*Gev) { cout << "Electron: failed pT veto.";RET;}
	    if (l.abseta() > 2.5){ cout << "Electron: failed eta veto.";RET;}
	    if ( !(l.pT() < 25*Gev) && !(l.abseta() > 2.5) ) {
	      cout << "Electron->select!";RET;
	      selectedCutLeptons.push_back(l);
	    }
	  }
	  else if (l.pid == PID::MUON || l.pid == -PID::MUON) {
	    if (l.pT() < 25*Gev) { cout << "Muon: failed pT veto.";RET;}
	    if (l.abseta() > 2.5){ cout << "Muon: failed eta veto.";RET;}
	    if ( !(l.pT() < 25*Gev) && !(l.abseta() > 2.5) ) {
	      cout << "Muon->select!!";RET;
	      selectedCutLeptons.push_back(l);
	    }
	  }
	  else {cout<< "lepton wierdness(2)!"; RET; N1; RET;}
	}
      }
      cout << "end of lepton cuts test*************************************"; RET;
      // end of lepton cuts test*************************************
      
      // find muons that are seperated form jets.
      cout << "leptons.size()=" << leptons.size(); RET;
      if (!leptons.empty() ){
      cout << "searching through DressedLeptons...";RET;
      for (const Particle& l : leptons){
	if (l.abspid() == PID::MUON){
	  cout << "found:muon";RET;
	  for (const Jet& j : jets){
	    if (deltaR(j, l) < 0.4) tooClose = true;
	  }
	  
	  if (!tooClose) {
	    cout << "muon->select!";RET;
	    selectedMuons.push_back(l);
	  } else {
	    tooClose = false; //< reset flag
	    cout << "muon->reject!";RET;
	  }
	}
	
	else if (l.abspid() == PID::ELECTRON) {
	  cout << "found:electron";RET;
	  electrons.push_back(l);
	} else { WORRY;N1; cout << "NO LEPTONS!"<<endl;}
      }
      cout << "DressedLeptons...done!";RET;
      }

      
      // remove jet closest to electron if closer than deltaR = 0.2
     
      if (!electrons.empty()){
	cout << "removing jet closest to electron if closer than deltaR = 0.2..."; RET; 
	for (const Particle& e : electrons){
	  cout << "considering electron...";RET;
	  cout << "considering " << jets.size() << " jets";RET;
	  for (const Jet& j : jets){ // find minimum dR(e,j)
	   
	    if (dR_min == 0.0) { dR_min = deltaR(j, e);}
	    else if (deltaR(j,e) < dR_min) {dR_min = deltaR(j, e);}
	  }
	  cout << "dR_min=" << dR_min;RET; 
	// remove minimal dR jet if dR_min(e,j) < 0.2
	for (const Jet& j : jets){
	  
	  if ( dR_min < 0.2 && deltaR(j, e) > dR_min) {
	    cout << "jet->keep"; RET;
	    selectedJetsEl.push_back(j);
	  }
	  else if ( dR_min < 0.2 && deltaR(j, e) <= dR_min ){
	    cout << "jet->reject";RET;
	    rejectN++;
	  }
	  else if ( dR_min >= 0.2 ) {
	    cout << "jet->keep"; RET;
	    selectedJetsEl.push_back(j);
	  }
	}
	if (rejectN > 1) {WORRY; N1; cout<<"more than one jet under minimal j-e deltaR!"; RET;}
	if (rejectN == 0 && jets.size() != selectedJetsEl.size()) {WORRY; N1; cout<<"Something is wrong at ~L220";RET;}    
	rejectN = 0;
	dR_min = 0.0;
	jets.clear();
	for (const Jet& sj : selectedJetsEl){
	  jets.push_back(sj);
	}
	selectedJetsEl.clear();
	}

      }


      
      // remove jet closest to photon if closer than deltaR = 0.1
      if (!hePhotons.empty()){
	cout << "removing jet closest to photon if dR_min < 0.1..."; RET;
	cout << "this should be a zero->" << rejectN ; RET;
	for (const Particle& ph : hePhotons){
	cout<< "considering photon...";RET;
	cout << "considering " << jets.size() << " jets";RET;
	  for (const Jet& j : jets){ // find minimum dR
	    if (dR_min == 0.0) {dR_min = deltaR(j, ph);}
	    else if (deltaR(j, ph) < dR_min) {dR_min = deltaR(j, ph);}
	  }
	  cout << "dR_min=" << dR_min;RET; 
	  //remove minimal dR jet if dR_min < 0.1
	  for (const Jet& j : jets){
	      
	    if (dR_min < 0.1 && deltaR(j, ph) > dR_min) {
	      cout<<"jet->keep"; RET;
	      selectedJets.push_back(j);
	    }
	    else if ( dR_min < 0.1 && deltaR(j, ph) <= dR_min ) {
	      cout << "jet->reject"; RET;
	      rejectN++;
	    }
	    else if ( dR_min >= 0.1 ) {
	      cout << "jet->keep"; RET;
	      selectedJets.push_back(j);
	    }
	  }
	  if (rejectN > 1) {WORRY; N1; cout<<"more than one jet under minimal j-ph deltaR!"; RET;}
	  if (rejectN != 0 && jets.size() == selectedJets.size()) {WORRY; N1; cout<<"Something is wrong at ~L259";RET;}    
	  rejectN = 0;
	  dR_min = 0.0;
	  jets.clear();
	  for (const Jet& sj : selectedJets){
	    jets.push_back(sj);
	  }
	  selectedJets.clear();
	}
      }
    

      // For print out only - remove in final anaysis -
      // NB: must change selectedJets->jets in all the following code:	
      cout << "select kept jets...";RET;
      for (const Jet& j : jets){
	selectedJets.push_back(j);
	cout << "jet->select!";RET;
      }
      
    

      // remove electrons that are within delta R = 0.4 of a jet
      if (!electrons.empty()){
	cout << "removing any electrons that are within delta R = 0.4 of a jet"; RET;
	cout << "considering " << electrons.size() << " electrons";RET;
	for (const Particle& e : electrons){
	 	  for (const Jet& j : selectedJets){
	    if ( deltaR(e, j) < 0.4 ) tooClose = true;
	  }
	  if (!tooClose) {
	    cout<<"electron->select!"; RET;
	    selectedElectrons.push_back(e);
	  } else {
	    cout<<"electron->reject"; RET;
	  }
	  tooClose = false; //< reset flag

	}
      }

      
      // find and count b-jets
      cout << "find b-jets...";RET;
      for (const Jet& j : selectedJets){
	if ( j.bTagged(Cuts::pT>5*GeV) )  bJetN++ ;
      }
      if (bJetN!=0){cout << "found " << bJetN << " b-jets!"; RET;} else {cout << "no b-jets found"; RET;}

      
      /// event selection in fiducial phase space based on the above object collections.
    
      cout << "-Object:muons=" << selectedMuons.size(); RET;
      cout << "-Object:electrons=" << selectedElectrons.size(); RET;
      cout << "-Object:jets=" << selectedJets.size(); RET;
      cout << "-Object:b-jets=" << bJetN; RET;
      cout << "-Object:photons=" << hePhotons.size(); RET; 

      cout << "VetoCode:"; RET;
      ////////////////////////////////////////////////////////////////////////////////////
      if ( selectedElectrons.size() + selectedMuons.size() != 1 ) { RET; vetoEvent;}
      if ( selectedJets.size() < 4)                               {N3; RET; vetoEvent;}
      if ( bJetN == 0 )                                           {N4; RET; vetoEvent;}
      if ( hePhotons.size() != 1 )                                {N5; RET; vetoEvent;}
      if (selectedElectrons.size() == 1){
	if ( deltaR(selectedElectrons[0], hePhotons[0]) < 0.7) {N6; RET; vetoEvent;}
      }
      if (selectedMuons.size() ==1 ){
	if ( deltaR(selectedMuons[0], hePhotons[0])  < 0.7) {N7; RET; vetoEvent;}
      }
      for (Jet& j : selectedJets){
	if (deltaR(j, hePhotons[0]) < 0.5){N8; RET; vetoEvent;}
      }

      ////////////////////////////////////////////////////////////////////////////////////

      ////////////////////////////// VETOS w/o CODES /////////////////////////////////////
      /*
      if ( selectedElectrons.size() + selectedMuons.size() != 1 ) {vetoEvent;}
      if ( selectedJets.size() < 4)                               {vetoEvent;}
      if ( b == 0 )                                               {vetoEvent;}
      if ( hePhotons.size() != 1 )                                {vetoEvent;}
      if (selectedElectrons.size() == 1){
	if ( deltaR(selectedElectrons[0], hePhotons[0]) < 0.7) {vetoEvent;}
      }
      if (selectedMuons.size() ==1 ){
	if ( deltaR(selectedMuons[0], hePhotons[0])  < 0.7) {vetoEvent;}
      }
      for (Jet& j : selectedJets){
	if (deltaR(j, hePhotons[0]) < 0.5){vetoEvent;}
	}
      */
      ////////////////////////////////////////////////////////////////////////////////////

     

      
      // fill histos
      N9;RET;
      SUCCESS;
      //_h_PhotonAbsEta->fill(pseudoTop4p.abseta(), event.weight()) ;
      //_h_PhotonPt->fill(pseudoTop4p.pt(), event.weight()) ;
      cout << "END"; RET;
      _c_fidXsec->fill(event.weight()); //error counter                  //



      
      // DEBUGGING:
      
      //////////////////// HEPMC OUT (START) ////////////////////////////////
      /*                                                                   //
      if (electronpartontops.size() > electrons.size()){                   //
	cout << "   !!!" << endl;                                          //
	_c_error_e->fill(event.weight()); //error counter                  //
	                                                                   //
	// print event record to file                                      //
	const HepMC::GenEvent* ge = event.genEvent();                      //
	if (ge == nullptr) {cout << "nullpointer found: GenEvent*"<<endl;  //
	} else { //if (ge != nullptr)                                      //
	  ge->print(); // sends (formatted) HepMC event to screen          //
	  _hepmcout->write_event(ge); // writes to file                    //
	}                                                                  //
      }                                                                    //
      */                                                                   //
      //////////////////// HEPMC OUT (END) //////////////////////////////////


    }

    
    /// Normalise histograms etc., after the run
    void finalize() {

      double SF = crossSection()*femtobarn/sumOfWeights();  // scale factor

      scale(_c_fidXsec, SF);
      //scale(_c_fid_tbar, SF);

      //scale(_h_PhotonPt, SF);
      //scale(_h_PhotonAbsEta, SF);
 
      
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
    CounterPtr _c_fidXsec;
    Histo1DPtr _h_PhotonAbsEta, _h_PhotonPt;


    //@}

    std::unique_ptr<HepMC::IO_GenEvent> _hepmcout;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(tempTopAnalysis_I1604029);


}
