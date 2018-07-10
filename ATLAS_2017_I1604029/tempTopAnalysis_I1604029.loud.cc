// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/JetUtils.hh"
#include "Rivet/Tools/ParticleUtils.hh"
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

      // Particle-level electron, muon and high energy photon cuts
      Cut e_muCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV ;
      Cut e_muForDressingCuts = Cuts::abseta < 2.7 && Cuts::pT > 10*GeV ;
      Cut hePhotonCuts = Cuts::abseta < 2.37 && Cuts::Et > 15*GeV ;

      // Photon projection to dress leptons for lepton & jet definitions
      //NBool: PromptFinalState(const Cut& c, bool accepttaudecays=false, bool acceptmudecays=false); (bools in same order when no cuts supplied)
      PromptFinalState hardParticles(true, true); //selects particle with no hadrons in ancestor line.
      IdentifiedFinalState photons(hardParticles);
      photons.acceptIdPair(PID::PHOTON);

      
      /// Leptons (Propper)
      PromptFinalState hardLeptonCuts_pfs(e_muForDressingCuts, true, true);
      IdentifiedFinalState hardLeptons_ifs(hardLeptonCuts_pfs, {PID::ELECTRON, PID::MUON, -PID::ELECTRON, -PID::MUON});
      DressedLeptons dressedLeptons(photons, hardLeptons_ifs, 0.1, e_muCuts, false);//< for testing cuts (NB: not much impact)) NB: last bool ()accept decay photons) doesnt matter because the photons input is hard anyway!
      //DressedLeptons dressedLeptons(photons, leptons_fs, 0.1, e_muCuts, true);//< for testing cuts
      declare(dressedLeptons, "DressedLeptons");

      /// Leptons (Testing)
      PromptFinalState promptParticles_all(true, true); //accept muon and tau decays
      IdentifiedFinalState allLeptons_fs(promptParticles_all, {PID::ELECTRON, PID::MUON, -PID::ELECTRON, -PID::MUON});//< for testing cuts
      //IdentifiedFinalState allLeptons_fs(e_muForDressingCuts);//< for testing cuts
      //allLeptons_fs.acceptIdPair(PID::MUON);//< for testing cuts
      //allLeptons_fs.acceptIdPair(PID::ELECTRON);//< for testing cuts
      DressedLeptons dressedLeptonsNoCuts(photons, allLeptons_fs, 0.1, Cuts::open(), false);//< for testing cuts
      //*DressedLeptons::DressedLeptons(const FinalState& photons, const FinalState& bareleptons, double dRmax, const Cut& cut, bool useDecayPhotons)
      declare(dressedLeptonsNoCuts, "DressedLeptonsNoCuts");

      // projections for jets
      IdentifiedFinalState neutrinos(fs);
      neutrinos.acceptNeutrinos();
      IdentifiedFinalState muons(fs, {PID::MUON, -PID::MUON});
      //DressedLeptons dressedMuons(photons, muons, 0.1, e_muCuts, false); // <- I have chosen to veto non-dressed muons - paper is vague
      VetoedFinalState jet_input(fs);
      jet_input.addVetoOnThisFinalState(muons);
      //jet_input.addVetoOnThisFinalState(dressedMuons);
      jet_input.addVetoOnThisFinalState(neutrinos);
      //declare(jet_input, "Jet_Input");
      declare(FastJets(jet_input, FastJets::ANTIKT, 0.4), "Jets");

      // Projection for high energy photons (proper)
      PromptFinalState hePhotonCuts_pfs(hePhotonCuts, true, true);
      IdentifiedFinalState hePhotons_ifs(hePhotonCuts_pfs);
      hePhotons_ifs.acceptIdPair(PID::PHOTON);
      declare(hePhotons_ifs, "HighEnergyPhotons");

      // Projection for high energy photons (testing)
      PromptFinalState heParticlesNoCuts(true, true);
      IdentifiedFinalState hePhotonsNoCuts(heParticlesNoCuts);
      hePhotonsNoCuts.acceptIdPair(PID::PHOTON);
      declare(hePhotonsNoCuts, "highEnergyPhotonsNoCuts");
      IdentifiedFinalState allPhotons(fs);
      allPhotons.acceptIdPair(PID::PHOTON);
      declare(allPhotons, "AllPhotons");

      
      // Projection for Parton Level Top Quarks
      //declare(PartonicTops(PartonicTops::E_MU, e_muCuts), "LeptonicTops") ;
      //declare(PartonicTops(PartonicTops::ALL, e_muCuts), "AllPartonicTops") ;
      //declare(PartonicTops(PartonicTops::E_MU), "LeptonicTops") ;
      //declare(PartonicTops(PartonicTops::ALL ), "AllPartonicTops") ;
      
      // Book counters
      _c_fidXsec = bookCounter("fidTotXsecttbargamma");
      //_c_fid_tbar = bookCounter("fidTotXsectbarq");

      // book histograms
      _h_PhotonPt = bookHisto1D(1,1,1);
      _h_PhotonAbsEta = bookHisto1D(1,1,2);

            
      // debug hepmc print-out
      _hepmcout.reset(new HepMC::IO_GenEvent("inspectevents.hepmc"));    

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      BEGIN;
          
      
      // Find leptons, jets and high energy photons
      const Particles& leptons = apply<DressedLeptons>(event, "DressedLeptons").particlesByPt();
      const Particles& hePhotons = apply<IdentifiedFinalState>(event, "HighEnergyPhotons").particlesByPt();
      const Particles& photons = apply<IdentifiedFinalState>(event, "AllPhotons").particlesByPt();
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(); //<keep non const for filtering
      
      // for testing:
      const Particles& leptonsNoCuts = apply<DressedLeptons>(event, "DressedLeptonsNoCuts").particlesByPt();
      const Particles& hePhotonsNoCuts = apply<IdentifiedFinalState>(event, "highEnergyPhotonsNoCuts").particlesByPt();
      

      
      /// overlap object removal
      // From the paper:
      //    "Muons within a cone of ∆R = 0.4 around a jet are removed, as well as the closest jet within a cone of
      //    ∆R = 0.2 (0.1) around an electron (photon). Finally, electrons within a cone of ∆R = 0.4 around a jet are
      //    removed."
      
      Particles electrons, selectedMuons, selectedElectrons;
      Jets selectedJetsEl, selectedJets;
      Jet bJet;
      bool tooClose = false ;

      int rejectN=0, bJetN=0 ;
	
      // jet selection
      Cut jetCuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV ;
      cout << "Jets=" << jets.size(); RET;
      ifilter_select(jets, jetCuts);
      cout << "JetsInFiducialPhaseSpace=" << jets.size(); RET;
      
      // emediate (but non-essential) veto if not enough jets
      
      if ( jets.size() < 4) {
	cout << "not enough jets to proceed. veto!";RET;
	N3; RET; //NB: N3 is used twice! (see event selection section ~L460)
	vetoEvent;
      }
      
      
     
      // testing lepton cuts:****************************************
      
      cout << "START: Lepton Testing>>>>>>>>>>>>>>>>>>>>>>>>>>";RET;
      Particles electronsNoCuts;
      Particles muonsNoCuts;
      Particles freshCutLeptons;
      Particles selectedCutLeptons;
      if (leptonsNoCuts.empty() ) {
	cout << " no dressed (noCut)leptons";RET;
      //////////////////// HEPMC OUT (START) ////////////////////////////////
	/*
      //if ( selectedElectrons.size() + selectedMuons.size() != 1 )
      //{                   //
	//cout << "   !!!" << endl;                                          //
	//	_c_error_e->fill(event.weight()); //error counter                  //
	                                                                   //
	// print event record to file                                      //
	const HepMC::GenEvent* ge = event.genEvent();                      //
	if (ge == nullptr) {cout << "nullpointer found: GenEvent*"<<endl;  //
	} else { //if (ge != nullptr)                                      //
	  ge->print(); // sends (formatted) HepMC event to screen          //
	  _hepmcout->write_event(ge); // writes to file                    //
	}                                                                  //
	//}
	*/                                                                    //
      //////////////////// HEPMC OUT (END) //////////////////////////////////
      } else { //if there ARE leptons with no cuts
	cout << " searching through DressedLeptonsNoCuts...";RET;
	for (const Particle& l : leptonsNoCuts){
	  if (l.abspid() == PID::MUON){
	    cout << " found:muon";
	    if (l.abseta() > 2.5) {cout << "...failed eta cut";
	      if (l.pT() < 25*GeV ) {cout << " and pt cut";}
	      RET;
	    }
	    else if (l.pT() < 25*GeV ) {cout << "...failed pt cut";RET;} 
	    else {cout << "...passed cuts!";RET;}
	    muonsNoCuts.push_back(l);
	  }
	  else if (l.abspid() == PID::ELECTRON){
	    cout << " found:electron";
	    if (l.abseta() > 2.5) {cout << "...failed eta cut";
	      if (l.pT() < 25*GeV ) {cout << " and pt cut";}
	      RET;
	    }
	    else if (l.pT() < 25*GeV ) {cout << "...failed pt cut"; RET;}
	    else {cout << "...passed cuts!";RET;}
	    electronsNoCuts.push_back(l);
	  }
	  else {cout<< " lepton wierdness!"; RET; N1; RET;}
	}
      
	//cout << "electronsNoCuts.size()=" << electronsNoCuts.size(); RET;
	//cout << "muonsNoCuts.size()=" << muonsNoCuts.size(); RET;
       
	for (const Particle& l : electronsNoCuts){freshCutLeptons.push_back(l);}
	for (const Particle& l : muonsNoCuts){freshCutLeptons.push_back(l);}
	int lepCtr = 0;
	for (const Particle& l : freshCutLeptons){
	  lepCtr++;
	  cout << " considering uncut dressed lepton No. " << lepCtr; RET;
	  if (l.pid() == PID::ELECTRON || l.pid() == -PID::ELECTRON) {
	    //if (l.pT() < 25*GeV) { cout << "Electron: failed pT veto.";RET;}
	    //if (l.abseta() > 2.5){ cout << "Electron: failed eta veto.";RET;}
	    if ( !(l.pT() < 25*GeV) && !(l.abseta() > 2.5) ) {
	      cout << " Electron->SELECT!";RET;
	      selectedCutLeptons.push_back(l);
	    }
	  }
	  else if (l.pid() == PID::MUON || l.pid() == -PID::MUON) {
	    //if (l.pT() < 25*GeV) { cout << "Muon: failed pT veto.";RET;}
	    //if (l.abseta() > 2.5){ cout << "Muon: failed eta veto.";RET;}
	    if ( !(l.pT() < 25*GeV) && !(l.abseta() > 2.5) ) {
	      cout << " Muon->SELECT!";RET;
	      selectedCutLeptons.push_back(l);
	    }
	  }
	  else {cout<< " lepton wierdness(2)!"; RET; N1; RET;}
	}
      }
      cout << "END: Lepton Testing<<<<<<<<<<<<<<<<<<<<<<<<<<<<";RET;     
      // end of lepton cuts test*************************************
     
      
      // find muons that are seperated form jets.
      // From paper: "Muons within a cone of ∆R = 0.4 around a jet are removed"
      cout << "leptons.size()=" << leptons.size(); RET;
      if (leptons.empty()) {
	cout << "No leptons in event! Veto!" ; RET;
	N2; RET;
	vetoEvent;
      }
	  cout << "searching through DressedLeptons...";RET;
	for (const Particle& l : leptons){
	  if (l.abspid() == PID::MUON){
	    cout << "found:muon";RET;
	    for (const Jet& j : jets){
	      if (deltaR(j, l) < 0.4) tooClose = true;
	  }
	  
	  if (!tooClose) {
	    cout << "muon->SELECT!";RET;
	    selectedMuons.push_back(l);
	  } else {
	    tooClose = false; //< reset flag
	    cout << "muon->reject!";RET;
	  }
	}
	
	else if (l.abspid() == PID::ELECTRON) {
	  cout << "found:electron";RET;
	  electrons.push_back(l);
	}
      }
      

	
      // remove jet closest to electron if closer than deltaR = 0.2
      // From paper: "closest jet within a cone of ∆R = 0.2 around an electron [is removed]."      
      double el_dR_min = 0.2 ;
      if (!electrons.empty()){
	cout << "attempt to remove jet closest to electron(s) if closer than deltaR = 0.2...";RET; 
	for (const Particle& e : electrons){
	  cout << "considering electron...";RET;
	  cout << "considering " << jets.size() << " jets";RET;
	  for (const Jet& j : jets){ // find minimum dR(e,j)
	    if (deltaR(j,e) < el_dR_min) {el_dR_min = deltaR(j, e);}
	  }
	  cout << "el_dR_min=" << el_dR_min;RET;
	  cout<< "NB: if the above number is 0.2 then it probably means that no jet was found within 0.2 of an electron";RET;
	  // remove minimal dR jet if el_dR_min < 0.2
	  for (const Jet& j : jets){
	  
	    if ( el_dR_min < 0.2 && deltaR(j, e) > el_dR_min) {
	      cout << "jet->keep"; RET;
	      selectedJetsEl.push_back(j);
	    }
	    else if ( el_dR_min < 0.2 && deltaR(j, e) == el_dR_min ){
	    cout << "jet->reject";RET;
	      rejectN++;// counter for testing/sanity/debugging
	    }
	    else if ( el_dR_min == 0.2 ) {
	      cout << "jet->keep"; RET;
	      selectedJetsEl.push_back(j);
	    }
	  }
	  if (rejectN > 1) {WORRY; cout<<"more than one jet under minimal j-e deltaR!"; RET;}
	  if (rejectN == 0 && jets.size() != selectedJetsEl.size()) {WORRY; cout<<"Something is wrong at POINT1.0";RET;}    
	  rejectN = 0;
	  el_dR_min = 0.2;
	  jets.clear();
	  for (const Jet& sj : selectedJetsEl){
	    jets.push_back(sj);
	  }
	  selectedJetsEl.clear();
	}
      }
      
      // Test photons: START****************************************************************
      cout << "START: Photon Testing>>>>>>>>>>>>>>>>>>>>>>>>>>";RET;
      Particles freshCutPhotons, freshCutPhotonsNoHad;
      int hadronicPhotonsCtr=0;
      if (hePhotonsNoCuts.empty() ) {
	cout << " no prompt photons"; RET; 
      } else { //if there ARE photons with no cuts
	cout << " searching through hePhotonsNoCuts (prompt photons with no cuts)...";RET;
	for (const Particle& p : hePhotonsNoCuts){
	  if (p.abspid() == PID::PHOTON){
	    cout << " found:photon";
	    if (p.abseta() > 2.37) {cout << "...failed eta cut";
	      if (p.Et() < 15*GeV ) {cout << " and Et cut";}
	      RET;
	    }
	    else if (p.Et() < 15*GeV ) {cout << "...failed Et cut";RET;} 
	    else {freshCutPhotons.push_back(p); cout << "...passed cuts!"; RET;}
	  }
	  else {cout<< " photon wierdness!"; RET; N1; RET;}
	  cout << " photon.Et=" << p.Et();RET;
	}
      }

      if (photons.empty() ) {
	cout << " no photons at all!"; RET; 
      } else { //if there ARE photons with no cuts
	cout << " searching through photons (ie ALL final state photons)...";RET;
	for (const Particle& p : photons){
	  if (p.abspid() == PID::PHOTON){

	    if (p.fromHadron()) {hadronicPhotonsCtr++;continue;}
	    else {cout << " found:nonHadronicPhoton";
	      if (p.abseta() > 2.37) {cout << "...failed eta cut";
		if (p.Et() < 15*GeV ) {cout << " and pt cut";}
		RET;
	      }
	      else if (p.Et() < 15*GeV ) {cout << "...failed pt cut";RET;} 
	      else {freshCutPhotonsNoHad.push_back(p); cout << "...passed cuts!"; RET;}
	    }
	  }
	  else {cout<< " photon wierdness!"; RET; N1; RET;}
	}
      }

      //RET;
      cout << " (photons from hadronic decays) hadronicPhotonsCtr=" << hadronicPhotonsCtr;     RET;
      cout << " ('non-hadronic' photons passing eta and Et cuts) freshCutPhotonsNoHad.size()=" << freshCutPhotonsNoHad.size();    RET;
      cout << " (prompt photons) hePhotonsNoCuts.size()=" << hePhotonsNoCuts.size(); RET;
      cout << " (prompt photons with cuts) freshCutPhotons.size()=" << freshCutPhotons.size(); RET;
      cout << "END: Photons Testing<<<<<<<<<<<<<<<<<<<<<<<<";RET;
      // Test photons: END***********************************************************************************

      
      // remove jet closest to photon if closer than deltaR = 0.1
      // from paper: "closest jet within a cone of ∆R = 0.1 around a photon [is removed]"
      cout << "attempt to remove jet closest to photon if dR_min < 0.1..."; RET;
      if (hePhotons.empty()) {
	cout << "no high energy photons in event! Veto!"; RET;
	N5; RET;
	vetoEvent;
      }
      cout << "high energy photons found in event";RET;
      double ph_dR_min = 0.1 ;
      for (const Particle& ph : hePhotons){
	cout<< "considering photon...";RET;
	cout << "considering " << jets.size() << " jets";RET;
	for (const Jet& j : jets){ // find minimum dR
	  if (deltaR(j, ph) < ph_dR_min) {ph_dR_min = deltaR(j, ph);}
	}
	cout << "ph_dR_min=" << ph_dR_min;RET;
	cout << "NB: if the above number is 0.1 then it probably means no jet was found within 0.1 of a photon";RET;
	  //remove minimal dR jet if ph_dR_min < 0.1
	  for (const Jet& j : jets){
	      
	    if (ph_dR_min < 0.1 && deltaR(j, ph) > el_dR_min) {
	      cout<<"jet->keep"; RET;
	      selectedJets.push_back(j);
	    }
	    else if ( ph_dR_min < 0.1 && deltaR(j, ph) == ph_dR_min ) {
	      cout << "jet->reject"; RET;
	      rejectN++;
	    }
	    else if ( ph_dR_min == 0.1 ) {
	      cout << "jet->keep"; RET;
	      selectedJets.push_back(j);
	    }
	  }
	  if (rejectN > 1) {WORRY; cout<<"^^more than one jet under minimal j-ph deltaR!"; RET;}
	  if (rejectN != 0 && jets.size() == selectedJets.size()) {WORRY; cout<<"Something is wrong at POINT0.2";RET;}    
	  rejectN = 0;
	  ph_dR_min = 0.1;
	  jets.clear();
	  for (const Jet& sj : selectedJets){
	    jets.push_back(sj);
	  }
	  selectedJets.clear();
      }
      
    

      // For print out only - remove in final anaysis -
      // NB: must change selectedJets->jets in all the following code:	
      cout << "select kept jets...";RET;
      for (const Jet& j : jets){
	selectedJets.push_back(j);
	cout << "jet->SELECT!";RET;
      }

      
      // remove electrons that are within delta R = 0.4 of a jet
      if (!electrons.empty()){
	cout << "attempting to removie any electrons that are within deltaR = 0.4 of a jet"; RET;
	cout << "considering " << electrons.size() << " electrons";RET;
	for (const Particle& e : electrons){
	  for (const Jet& j : selectedJets){
	    if ( deltaR(e, j) < 0.4 ) tooClose = true;
	  }
	  if (!tooClose) {
	    cout<<"electron->SELECT!"; RET;
	    selectedElectrons.push_back(e);
	  } else {
	    cout<<"electron->reject"; RET;
	    tooClose = false; //< reset flag
	  }
	}
      }

      
      // find and count b-jets
      cout << "find b-jets...";RET;
      for (const Jet& j : selectedJets){
	if ( j.bTagged(Cuts::pT>5*GeV) )  bJetN++ ;
      }
      if (bJetN!=0){cout << "found " << bJetN << " b-jets!"; RET;} else {cout << "no b-jets found"; RET;}

      /*
      // Test photons: START****************************************************************
      cout << "START: Photon Testing>>>>>>>>>>>>>>>>>>>>>>>>>>";RET;
      Particles freshCutPhotons, freshCutPhotonsNoHad;
      int hadronicPhotonsCtr=0;
      if (hePhotonsNoCuts.empty() ) {
	cout << " no prompt photons"; RET; 
      } else { //if there ARE photons with no cuts
	cout << " searching through hePhotonsNoCuts (prompt photons with no cuts)...";RET;
	for (const Particle& p : hePhotonsNoCuts){
	  if (p.abspid() == PID::PHOTON){
	    cout << " found:photon";
	    if (p.abseta() > 2.37) {cout << "...failed eta cut";
	      if (p.Et() < 15*GeV ) {cout << " and pt cut";}
	      RET;
	    }
	    else if (p.Et() < 15*GeV ) {cout << "...failed pt cut";RET;} 
	    else {freshCutPhotons.push_back(p); cout << "...passed cuts!"; RET;}
	  }
	  else {cout<< " photon wierdness!"; RET; N1; RET;}
	}
      }

      if (photons.empty() ) {
	cout << " no photons at all!"; RET; 
      } else { //if there ARE photons with no cuts
	cout << " searching through photons (ie ALL final state photons)...";RET;
	for (const Particle& p : photons){
	  if (p.abspid() == PID::PHOTON){

	    if (p.fromHadron()) {hadronicPhotonsCtr++;continue;}
	    else {cout << " found:nonHadronicPhoton";
	      if (p.abseta() > 2.37) {cout << "...failed eta cut";
		if (p.Et() < 15*GeV ) {cout << " and pt cut";}
		RET;
	      }
	      else if (p.Et() < 15*GeV ) {cout << "...failed pt cut";RET;} 
	      else {freshCutPhotonsNoHad.push_back(p); cout << "...passed cuts!"; RET;}
	    }
	  }
	  else {cout<< " photon wierdness!"; RET; N1; RET;}
	}
      }

      //RET;
      cout << " (photons from hadronic decays) hadronicPhotonsCtr=" << hadronicPhotonsCtr;     RET;
      cout << " ('non-hadronic' photons passing eta and Et cuts) freshCutPhotonsNoHad.size()=" << freshCutPhotonsNoHad.size();    RET;
      cout << " (prompt photons) hePhotonsNoCuts.size()=" << hePhotonsNoCuts.size(); RET;
      cout << " (prompt photons with cuts) freshCutPhotons.size()=" << freshCutPhotons.size(); RET;
      cout << "END: Photons Testing<<<<<<<<<<<<<<<<<<<<<<<<";RET;
          
      // Test photons: END***********************************************************************************
      */
      
      // Print out Selected Object totals
      
      cout << "-Object:muons=" << selectedMuons.size(); RET;
      cout << "-Object:electrons=" << selectedElectrons.size(); RET;
      cout << "-Object:jets=" << selectedJets.size(); RET;
      cout << "-Object:b-jets=" << bJetN; RET;
      cout << "-Object:photons=" << hePhotons.size(); RET; 
      


      
      /// event selection in fiducial phase space based on the above object collections.
      //cout << "VetoCode:"; RET;
      ////////////////////////////////////////////////////////////////////////////////////
      if ( selectedElectrons.size() + selectedMuons.size() != 1 ) {
	cout << "no single lepton, veto!";RET;
	N2; RET; vetoEvent;} // N2 already used if leptons.size()==0
      if ( selectedJets.size() < 4) {
	cout << "not enough jets, veto!"; RET;
	N3; RET; vetoEvent;}
      if ( bJetN == 0 ) {
	cout << "no bjets, veto!";RET;
	N4; RET; vetoEvent;}
      if ( hePhotons.size() != 1 ) {
	cout << "no single photon, veto!";RET;
	N5; RET; vetoEvent;} // N5 already used if hePhotons.size()==0
      if (selectedElectrons.size() == 1){
	if ( deltaR(selectedElectrons[0], hePhotons[0]) < 0.7) {
	  cout << "chosen lepton (electron) and photon are too close, veto!"; RET;
	  N6; RET; vetoEvent;}
      }
      if (selectedMuons.size() ==1 ){
	if ( deltaR(selectedMuons[0], hePhotons[0])  < 0.7) {
	  cout << "chosen lepton (muon) and photon are too close, veto!";RET;
	  N7; RET; vetoEvent;}
      }
      for (Jet& j : selectedJets){
	if (deltaR(j, hePhotons[0]) < 0.5){
	  cout << "jet and photon are too close, veto!";RET;
	  N8; RET; vetoEvent;}
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
      _h_PhotonAbsEta->fill(hePhotons[0].abseta(), event.weight()) ;
      _h_PhotonPt->fill(hePhotons[0].pt()/GeV, event.weight()) ;
      _c_fidXsec->fill(event.weight());
      

      
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

      double SF = crossSection()/femtobarn/sumOfWeights();  // scale factor

      scale({_h_PhotonAbsEta, _h_PhotonPt}, SF);
      scale(_c_fidXsec, SF);

      // user info:_________________________________________________________________________________
      RET;
      cout << "VETO CODES:";RET;
      cout << "#1 test warning!";RET; 
      cout << "   --> NB: no veto for test warning!! => Total = vetos + N_{#1}";RET;
      cout << "#2 no single lepton";RET;
      cout << "#3 not enough jets"; RET;
      cout << "#4 no bjets";RET;
      cout << "#5 no single photon";RET;
      cout << "#6 chosen lepton (electron) and photon are too close"; RET;
      cout << "#7 chosen lepton (muon) and photon are too close";RET;
      cout << "#8 jet and photon are too close";RET;
      RET;
      

      cout << "Calculated inclusive Xsec from single lepton channel in fiducial region is "
	   << _c_fidXsec->numEntries()*SF << "fb" << endl;

      cout << "Generator Xsec is:" << crossSection() << "pb" << endl;

      RET;
      //_________________________________________________________________________________________


      _hepmcout->clear(); _hepmcout.reset();
      
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
