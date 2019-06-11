// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"

#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
// #include "Rivet/Cuts.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Math/Vector3.hh"
#include "Rivet/Math/VectorN.hh"

// Vetoed Jet inputs
#include "Rivet/Projections/VetoedFinalState.hh"

// event record specific
#include "Rivet/Particle.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "HepMC/GenEvent.h"

#include <math.h>
#include <string>
#include <sstream>
#include <iostream>

// #include <TMath.h>
#include <TLorentzVector.h>
#include <TMinuit.h>

#include "neutrino.h"

using namespace std;

namespace patch
{
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}

namespace Rivet {
  using namespace Cuts;

  class ATLAS_2017_I1512776 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2017_I1512776()
      : Analysis("ATLAS_2017_I1512776")
    {
      setNeedsCrossSection(true);
    }
  public:

    //=====================================================================
    // Set up projections and book histograms
    //=====================================================================
    void init() {
      /// Dressing of the leptons taken from MC_TTbar_TruthSel
      /// Code based on MC_SingleTop_Truth, MC_TTbar_TruthSel and ATLAS_2014_I1304688

      Cut eta_full = (Cuts::abseta < 5) && (Cuts::pT >= 1.0*MeV);
      Cut lep_cuts = (Cuts::abseta < 2.5) && (Cuts::pT >= 20*GeV);
      Cut lep_veto_cuts = (Cuts::abseta < 2.5) && (Cuts::pT >= 10*GeV);

      // All final state particles
      FinalState fs(Cuts::abseta < 5);

      // Get photons to dress leptons
      IdentifiedFinalState photons(Cuts::abseta < 5);
      photons.acceptIdPair(PID::PHOTON);

      // Projection to find the electrons
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);

      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);
      addProjection(electrons, "electrons");
      DressedLeptons dressedelectrons(photons, electrons, 0.1,lep_cuts,true,true );
      addProjection(dressedelectrons, "dressedelectrons");
      DressedLeptons vetodressedelectrons(photons, electrons, 0.1, lep_veto_cuts, true, true);
      addProjection(vetodressedelectrons, "vetodressedelectrons");
      DressedLeptons ewdressedelectrons(photons, electrons, 0.1,  Cuts::abseta < 3,true, true);
      addProjection(ewdressedelectrons, "ewdressedelectrons");

      // Projection to find the muons
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);

      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      addProjection(muons, "muons");
      DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts, true, true);
      addProjection(dressedmuons, "dressedmuons");
      DressedLeptons vetodressedmuons(photons, muons, 0.1, lep_veto_cuts, true, true);
      addProjection(vetodressedmuons, "vetodressedmuons");
      DressedLeptons ewdressedmuons(photons, muons, 0.1,  eta_full, true,true);
      addProjection(ewdressedmuons, "ewdressedmuons");

      // Projection to find prompt neutrinos (definition from ATLAS_2014_I1304688)
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);
      addProjection(neutrinos, "neutrinos");
      
      // Projection to get MET from all neutrinos
      //IdentifiedFinalState neutrinos_met(-4.5, 4.5, 0.0*GeV);       
      //neutrinos_met.acceptNeutrinos();
      //addProjection(neutrinos_met, "neutrinos_met");
      IdentifiedFinalState nu_id_met;
      nu_id_met.acceptNeutrinos();
      PromptFinalState neutrinos_met(nu_id);
      neutrinos_met.acceptTauDecays(true);// was false??
      addProjection(neutrinos_met, "neutrinos_met");

      // Jet clustering.
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(dressedelectrons);
      vfs.addVetoOnThisFinalState(dressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos_met);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      addProjection(jets, "jets");

      // Projection to find the taus
      IdentifiedFinalState taus;
      taus.acceptIdPair(PID::TAU);
      addProjection(taus, "taus");

      // ChargedFinalState (charged final particles)
      ChargedFinalState cfs(-4.5, 4.5, 5.0*GeV);
      addProjection(cfs, "CFS");

      // booking and initialising histograms
      BookHistograms();
      
      // setting the counter for electron and muon events to 0
      _el_weights = 0;
      _mu_weights = 0;
//      cutnames.at(0)="Init";
      eventcount=0;
      
      // setting the cutflow to zero
      for(int i = 0; i < 9; i++) {
        cutflow[i] = 0;
        cutflow_ele[i] = 0;
        cutflow_mu[i] = 0;
      }
      
      for (int i = 0; i < 1; i++) {
        other_checks[i] = 0;
      }

      // ----- NEW HISTOGRAMS ------ //
      for (size_t itop = 0; itop < 2; ++itop) {
	//
        _h_AbsPtclDiffXsecTPt[itop]   = bookHisto1D(1,1,1+itop); // top/antitop particle-level pT normalised to top/antitop Xsec
        _h_AbsPtclDiffXsecTY[itop]    = bookHisto1D(1,1,3+itop); // top/antitop particle-level Y  normalised to top/antitop Xsec
        _h_NrmPtclDiffXsecTPt[itop]   = bookHisto1D(2,1,1+itop); // top/antitop particle-level pT normalised to 1
        _h_NrmPtclDiffXsecTY[itop]    = bookHisto1D(2,1,3+itop); // top/antitop particle-level Y  normalised to 1
        //
        _h_AbsPtclDiffXsecJPt[itop]   = bookHisto1D(5,1,1+itop); // top/antitop particle-level jet pT normalised to top/antitop Xsec
        _h_AbsPtclDiffXsecJY[itop]    = bookHisto1D(5,1,3+itop); // top/antitop particle-level jet Y  normalised to top/antitop Xsec
        _h_NrmPtclDiffXsecJPt[itop]   = bookHisto1D(6,1,1+itop); // top/antitop particle-level jet pT normalised to 1
        _h_NrmPtclDiffXsecJY[itop]    = bookHisto1D(6,1,3+itop); // top/antitop particle-level jet Y  normalised to 1
        //
        // _h_AbsPrtnDiffXsecTPt[itop]	  = bookHisto1D(3,1,1+itop); // top/antitop parton-level pT normalised to top/antitop Xsec
	// _h_AbsPrtnDiffXsecTY[itop]	  = bookHisto1D(3,1,3+itop); // top/antitop parton-level Y  normalised to top/antitop Xsec
	// _h_NrmPrtnDiffXsecTPt[itop]	  = bookHisto1D(4,1,1+itop); // top/antitop parton-level pT normalised to 1
	// _h_NrmPrtnDiffXsecTY[itop]	  = bookHisto1D(4,1,3+itop); // top/antitop parton-level Y  normalised to 1
      }
      // --------------------------- //

    }

    //=====================================================================
    // Perform the per-event analysis
    //=====================================================================
    void analyze(const Event& event) {
      // bool _debug = false;
      _weight = event.weight();
      _h_weight->fill(_weight, _weight);
      MSG_DEBUG("mc weight: " << _weight);
      

      // get the needed objects using the projections
      _dressedelectrons     = sortByPt(applyProjection<DressedLeptons>(event, "dressedelectrons").dressedLeptons());
      _vetodressedelectrons = sortByPt(applyProjection<DressedLeptons>(event, "vetodressedelectrons").dressedLeptons());
      _ewdressedelectrons   = sortByPt(applyProjection<DressedLeptons>(event, "ewdressedelectrons").dressedLeptons());
      _dressedmuons         = sortByPt(applyProjection<DressedLeptons>(event, "dressedmuons").dressedLeptons());
      _vetodressedmuons     = sortByPt(applyProjection<DressedLeptons>(event, "vetodressedmuons").dressedLeptons());
      _ewdressedmuons       = sortByPt(applyProjection<DressedLeptons>(event, "ewdressedmuons").dressedLeptons());
//      _photons     = applyProjection<PromptFinalState>(event, "photons").particlesByPt();
      //_neutrinos_met        = applyProjection<IdentifiedFinalState>(event, "neutrinos_met").particlesByPt();
      const Jets alljets = applyProjection<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 4.5);
      Jets good_bjets, good_ljets, good_jets;

      // select only good electrons and muons
      vector<DressedLepton> good_electrons;
      SelectGoodElectrons(_dressedelectrons, good_electrons, 25.0*GeV, 2.5);
      vector<DressedLepton> good_muons;
      SelectGoodMuons(_dressedmuons, good_muons, 25.0*GeV, 2.5);

      vector<DressedLepton> good_leptons;
      good_leptons = good_electrons + good_muons;

      // distribution of veto and ew leptons for comparision purposes
      vector<DressedLepton> veto_electrons;
      SelectGoodElectrons(_vetodressedelectrons, veto_electrons, 15.0*GeV, 2.5);
      vector<DressedLepton> ew_electrons;
      SelectGoodElectrons(_ewdressedelectrons, ew_electrons, 0.0*GeV, 5.0);
      
      vector<DressedLepton> veto_muons;
      SelectGoodMuons(_vetodressedmuons, veto_muons, 15.0*GeV, 2.5);
      vector<DressedLepton> ew_muons;
      SelectGoodMuons(_ewdressedmuons, ew_muons, 0.0*GeV, 5.0);

      // selection of the jets, remove event if overlap detected
      cutflow[0]++;
      //std::cout << "\tcutflow[0]++ \n";
      if (SelectGoodJets(alljets, good_bjets, good_ljets, good_jets, 30.0*GeV, 4.5, good_electrons, good_muons))
              vetoEvent;
      // selection of the prompt neutrinos
      const FinalState& NeutrinoCollection = applyProjection<FinalState>(event, "neutrinos");
      vector<Particle> good_neutrinos_prompt = returnPropertiesFromCollection(NeutrinoCollection);
      
      // selection of all neutrinos
      const FinalState& NeutrinoCollection_met = applyProjection<FinalState>(event, "neutrinos_met");
      vector<Particle> good_neutrinos_met = returnPropertiesFromCollection(NeutrinoCollection_met);

      // get all charged particles in the event
      const ChargedFinalState& ChargedParticleCollection = applyProjection<ChargedFinalState>(event, "CFS");

      // -----------------------------------------------------------------------------------------------

      // calculate MET from NeutrinoCollection for MET
      FourMomentum _p_met(0., 0., 0., 0.);
      foreach(const Particle& _p, NeutrinoCollection_met.particlesByPt()) { _p_met = _p_met + _p.momentum(); }
      // MSG_INFO(_p_met.pT()/GeV << "\t" << _p_met.phi());

      // -----------------------------------------------------------------------------------------------

      // return B-hadrons collection (pt cut)
      vector<const HepMC::GenParticle*> B_hadrons = returnBHadrons(event, 2.5, 5.0);

      // -----------------------------------------------------------------------------------------------
      
      /// Up to this point all needed particles are stored in some vector.
      /// Cuts on jets pT, eta and overlap have been applied as well as on lepton pT and eta

      // veto event if there are neither leptons nor neutrinos
      // event has to match the t-channel production with leptonic W decay
      cutflow[1]++;
      //std::cout << "\tcutflow[1]++ \n";
      if (good_leptons.size() != 1) {vetoEvent; }
      // rejection if veto particles are present
      cutflow[2]++;

      if (good_electrons.size() == 1 && (veto_muons.size() != 0 || veto_electrons.size() > 1)) { vetoEvent; }
      if (good_muons.size() == 1 && (veto_electrons.size() != 0 || veto_muons.size() > 1)) { vetoEvent; }

      // separate all events into top and antitop channel
      bool isTop = false;
      if (good_leptons[0].pdgId() < 0) {
        isTop = true;   // If there is an anti-lepton, the decaying top is an anti top
      } else {
        isTop = false;  // If there is a lepton, the decaying top is a top
      }
      
      // veto event if there are no good b-jets
      // event has to match the t-channel production with top decay into W + b
      cutflow[3]++;

//      cout << "No. of b-jets: " << good_bjets.size() << endl;
//      cout << "No. of l-jets: " << good_ljets.size() << endl;
//      cout << "No. of jets: " << good_jets.size() << endl;
      if (good_jets.size() != 2) { vetoEvent; }
      cutflow[4]++;

      if (good_ljets.size() != 1) { vetoEvent; }
      // veto event if there are no good non b-jets
      // Event has to match the t-channel production with a light-quark jet
      cutflow[5]++;
      if (good_bjets.size() != 1) { vetoEvent; }

      // veto event if eTmiss cut not passed
      cutflow[6]++;
//      if ( _p_met.pT()/GeV < 30.0 ) { vetoEvent; }
      
      // veto event if mtw cut not passed
      cutflow[7]++;
      //std::cout << "\tcutflow[7]++ \n";

      Vector3 bjet_v3=good_bjets[0].momentum().vector3();
      FourMomentum bjet=good_bjets[0].momentum();
      // This is needed to get same acceptance as using ROOT
      bjet.setE(sqrt(bjet_v3.x()*bjet_v3.x()+bjet_v3.y()*bjet_v3.y()+bjet_v3.z()*bjet_v3.z()));

      FourMomentum lepton_b= good_leptons[0].momentum()+bjet;

      if((lepton_b.mass())>=160.0) {vetoEvent;}
      // -----------------------------------------------------------------------------------------------
      
      /// After this point no additional cuts are applied and all histograms (except for the total weight in the first lines) are filled
      cutflow[8]++;

      // select the spectator quark (only one remaining after cuts)
      _spectjet_from_top_quark = good_ljets[0];
      if(isTop) {
        for(unsigned int i=0; i < _h_ljet_pt_top.size(); i++) {
          _h_ljet_pt_top[i]->fill(_spectjet_from_top_quark.pT()*GeV, _weight);
          _h_ljet_pt_top_norm[i]->fill(_spectjet_from_top_quark.pT()*GeV, _weight);
        }
        for(unsigned int i=0; i < _h_ljet_rap_top.size(); i++) {
          _h_ljet_rap_top[i]->fill(_spectjet_from_top_quark.absrap(), _weight);
          _h_ljet_rap_top_norm[i]->fill(_spectjet_from_top_quark.absrap(), _weight);
        }
      } else {
        for(unsigned int i=0; i < _h_ljet_pt_antitop.size(); i++) {
          _h_ljet_pt_antitop[i]->fill(_spectjet_from_top_quark.pT()*GeV, _weight);
          _h_ljet_pt_antitop_norm[i]->fill(_spectjet_from_top_quark.pT()*GeV, _weight);
        }
        for(unsigned int i=0; i < _h_ljet_rap_antitop.size(); i++) {
          _h_ljet_rap_antitop[i]->fill(_spectjet_from_top_quark.absrap(), _weight);
          _h_ljet_rap_antitop_norm[i]->fill(_spectjet_from_top_quark.absrap(), _weight);
        }
      }

      // Reconstruction of the W boson
      //Particle _W_boson = returnWboson_byMET_woFit(good_leptons[0], _p_met);

      Particle _W_boson = returnWboson_byMET_wFit(good_leptons[0], _p_met.x(), _p_met.y(), _p_met.pT(), _p_met.phi());

      // fill top quark
      // Reconstruction of the top quark and top antiqark
      if(isTop) {
        Particle _top_quark = returnTopQuark(_W_boson, good_bjets);
        for (unsigned int i=0; i< _h_t_pt_top.size(); i++){
          _h_t_pt_top[i]->fill(_top_quark.momentum().pT()/GeV, _weight);
          _h_t_pt_top_norm[i]->fill(_top_quark.momentum().pT()/GeV, _weight);
        }
        for (unsigned int i=0; i< _h_t_rap_top.size(); i++){
          _h_t_rap_top[i]->fill(_top_quark.momentum().absrap(), _weight);
          _h_t_rap_top_norm[i]->fill(_top_quark.momentum().absrap(), _weight);
        }
	
	// Fill new top quark histos 
	_h_AbsPtclDiffXsecTPt[0]->fill(_top_quark.pT(),                   event.weight());
	_h_AbsPtclDiffXsecTY[0] ->fill(_top_quark.absrap(),               event.weight());
	_h_NrmPtclDiffXsecTPt[0]->fill(_top_quark.pT(),                   event.weight());
	_h_NrmPtclDiffXsecTY[0] ->fill(_top_quark.absrap(),               event.weight());
	_h_AbsPtclDiffXsecJPt[0]->fill(_spectjet_from_top_quark.pT(),     event.weight());
	_h_AbsPtclDiffXsecJY[0] ->fill(_spectjet_from_top_quark.absrap(), event.weight());
	_h_NrmPtclDiffXsecJPt[0]->fill(_spectjet_from_top_quark.pT(),     event.weight());
	_h_NrmPtclDiffXsecJY[0] ->fill(_spectjet_from_top_quark.absrap(), event.weight());

      } else {
        Particle _top_antiquark = returnTopQuark(_W_boson, good_bjets);
        for (unsigned int i=0; i< _h_t_pt_antitop.size(); i++){
          _h_t_pt_antitop[i]->fill(_top_antiquark.momentum().pT()/GeV, _weight);
          _h_t_pt_antitop_norm[i]->fill(_top_antiquark.momentum().pT()/GeV, _weight);
        }
        for (unsigned int i=0; i< _h_t_rap_antitop.size(); i++){
          _h_t_rap_antitop[i]->fill(_top_antiquark.momentum().absrap(), _weight);
          _h_t_rap_antitop_norm[i]->fill(_top_antiquark.momentum().absrap(), _weight);
        }

	// Fill new antitop quark histos 
	_h_AbsPtclDiffXsecTPt[1]->fill(_top_antiquark.pT(),                   event.weight());
	_h_AbsPtclDiffXsecTY[1] ->fill(_top_antiquark.absrap(),               event.weight());
	_h_NrmPtclDiffXsecTPt[1]->fill(_top_antiquark.pT(),                   event.weight());
	_h_NrmPtclDiffXsecTY[1] ->fill(_top_antiquark.absrap(),               event.weight());
	_h_AbsPtclDiffXsecJPt[1]->fill(_spectjet_from_top_quark.pT(),     event.weight());
	_h_AbsPtclDiffXsecJY[1] ->fill(_spectjet_from_top_quark.absrap(), event.weight());
	_h_NrmPtclDiffXsecJPt[1]->fill(_spectjet_from_top_quark.pT(),     event.weight());
	_h_NrmPtclDiffXsecJY[1] ->fill(_spectjet_from_top_quark.absrap(), event.weight());
      }

      //std::cout << "\tcutflow[8]++ \n";
      if(good_electrons.size() == 1) cutflow_ele[8]++;
      if(good_muons.size() == 1) cutflow_mu[8]++;
      eventcount++;
      
    } // end of analyze()
    
    // Normalise histograms etc., after the run
    void finalize() {
      double lum = 20239.3;
      double Br = 0.324;
      cout << "eventcount: " << eventcount << endl;        
      
      cout << "============================= tchan particle ================================" << endl;
      cout << "crossSection " << crossSection() << " sumOfWeights " << sumOfWeights() << endl;
      
      double scale_fac = crossSection()/Br/sumOfWeights();
      
      cout << "Scale fac.: " << scale_fac << endl;
      //cout << "No. of electron events: " << _el_weights * scale_fac << endl;
      //cout << "No. of muon events: " << _mu_weights * scale_fac << endl;
      cout << "=============================================================================" << endl;
      for (int i = 0; i < 9; i++) {
        cout << "Cutflow at step " << i << ": " << cutflow[i]<<" Afid: "<<(float)cutflow[i]/(float)cutflow[0] << endl;
      }
      for (int i = 0; i < 1; i++) {
        cout << "No. of events with M_T(W) above M(W): " << other_checks[i] << endl;
      }
      cout << "=============================================================================" << endl;
      // Scale all histograms by their cross section

      for (unsigned int i=0; i< _h_t_pt_top.size(); i++){ 
        scale(_h_t_pt_top[i], scale_fac);
        scale(_h_t_pt_top_norm[i], scale_fac);
        normalize(_h_t_pt_top_norm[i]);
      }
      for (unsigned int i=0; i< _h_t_pt_antitop.size(); i++){ 
        scale(_h_t_pt_antitop[i], scale_fac);
        scale(_h_t_pt_antitop_norm[i], scale_fac);
        normalize(_h_t_pt_antitop_norm[i]);
      }
      for (unsigned int i=0; i< _h_t_rap_top.size(); i++){ 
        scale(_h_t_rap_top[i], scale_fac);
        scale(_h_t_rap_top_norm[i], scale_fac);
        normalize(_h_t_rap_top_norm[i]);
      }
      for (unsigned int i=0; i< _h_t_rap_antitop.size(); i++){ 
        scale(_h_t_rap_antitop[i], scale_fac);
        scale(_h_t_rap_antitop_norm[i], scale_fac);
        normalize(_h_t_rap_antitop_norm[i]);
      }
      for (unsigned int i=0; i< _h_ljet_pt_top.size(); i++){ 
        scale(_h_ljet_pt_top[i], scale_fac);
        scale(_h_ljet_pt_top_norm[i], scale_fac);
        normalize(_h_ljet_pt_top_norm[i]);
      }
      for (unsigned int i=0; i< _h_ljet_pt_antitop.size(); i++){ 
        scale(_h_ljet_pt_antitop[i], scale_fac);
        scale(_h_ljet_pt_antitop_norm[i], scale_fac);
        normalize(_h_ljet_pt_antitop_norm[i]);
      }
      for (unsigned int i=0; i< _h_ljet_rap_top.size(); i++){ 
        scale(_h_ljet_rap_top[i], scale_fac);
        scale(_h_ljet_rap_top_norm[i], scale_fac);
        normalize(_h_ljet_rap_top_norm[i]);
      }
      for (unsigned int i=0; i< _h_ljet_rap_antitop.size(); i++){ 
        scale(_h_ljet_rap_antitop[i], scale_fac);
        scale(_h_ljet_rap_antitop_norm[i], scale_fac);
        normalize(_h_ljet_rap_antitop_norm[i]);
      } 

      // scale new histos
      for (size_t itop = 0; itop < 2; ++itop) {      
	//
	scale( _h_AbsPtclDiffXsecTPt[itop], crossSection()/femtobarn/sumOfWeights()/Br );
	scale( _h_AbsPtclDiffXsecTY[itop],  crossSection()/sumOfWeights()/Br );
	normalize(_h_NrmPtclDiffXsecTPt[itop]);
	normalize(_h_NrmPtclDiffXsecTY[itop]);
	//
	scale( _h_AbsPtclDiffXsecJPt[itop], crossSection()/femtobarn/sumOfWeights()/Br );
	scale( _h_AbsPtclDiffXsecJY[itop],  crossSection()/sumOfWeights()/Br );
	normalize(_h_NrmPtclDiffXsecJPt[itop]);
	normalize(_h_NrmPtclDiffXsecJY[itop]);
	//
	// scale(_h_AbsPrtnDiffXsecTPt[itop], crossSection()/femtobarn/sumOfWeights()/Br );
	// scale(_h_AbsPrtnDiffXsecTY[itop],  crossSection() / sumOfWeights()/Br );
	// normalize(_h_NrmPrtnDiffXsecTPt[itop]);
	// normalize(_h_NrmPrtnDiffXsecTY[itop]);
      }


    }

    //=====================================================================
    // fill properties from Collection - FinalState
    //=====================================================================
    /// created a function to get the vector of particles and to fill the histograms
    vector<Particle> returnPropertiesFromCollection(const FinalState& Collection, float etaCut=5.0, float ptCut=0.0)
    {
      vector<Particle> _partvect;
      
      foreach(const Particle& _p, Collection.particlesByPt()) {
        if (_p.momentum().pT() > ptCut && fabs(_p.momentum().eta()) < etaCut) {
          _partvect.push_back(_p);
        }
      }
      
      return _partvect;
    }

    //=====================================================================
    // fill properties from Collection - UnstableFinalState
    //=====================================================================
    /// created a function to get the vector of particles and to fill the histograms
    vector<Particle> returnPropertiesFromCollection(const IdentifiedFinalState& Collection, int pdfId, float etaCut=5.0, float ptCut=0.0) {
      vector<Particle> _partvect;
      
      foreach(const Particle& _p, Collection.particlesByPt()) {
        if(abs(_p.pdgId())!=pdfId) continue;
        if (_p.momentum().pT() > ptCut && fabs(_p.momentum().eta()) < etaCut) {
          _partvect.push_back(_p);
        }
      }
      
      return _partvect;
    }
    
    //=====================================================================
    // fill B-Hadrons
    //=====================================================================
    // created a function to get the vector of particles and to fill the histograms
    vector<const HepMC::GenParticle*> returnBHadrons(const Event& event, float etaCut=4.5, float ptCut=5.0) {
      vector<const HepMC::GenParticle*> B_hadrons;
      // Get the B-Hadrons with pT > ptCut GeV, to add to the final-state particles used for jet clustering.
      vector<const HepMC::GenParticle*> allParticles = particles(event.genEvent());
      
      for (size_t i = 0; i < allParticles.size(); i++) {
        const HepMC::GenParticle* _p = allParticles.at(i);

        // If the particle is a B-hadron and has pT > ptCut GeV, add it to the B-hadrons vector
        if(PID::isHadron(_p->pdg_id()) && PID::hasBottom(_p->pdg_id()) && fabs(_p->momentum().eta()) < etaCut && _p->momentum().perp() > ptCut) {
          Particle B_hadron = Particle(*_p);
          
          B_hadrons.push_back(_p);
        }
      }
      
      return B_hadrons;
    }
    
    //=====================================================================
    // fill W boson
    //=====================================================================
    /// Neutrino reconstruction from MET via W mass constraint without met fit
    Particle returnWboson_byMET_woFit(Particle _lepton, FourMomentum _p_met) {
      double pz = neutrinoPz(_lepton, _p_met);
      FourMomentum _neutrino_woFit( sqrt(sqr(_p_met.px()) + sqr(_p_met.py()) + sqr(pz)), _p_met.px(), _p_met.py(), pz);

      double pz_neu = _p_met.pz();

      FourMomentum _Wboson_mom = _lepton.momentum() + _neutrino_woFit;

      _lepton_from_W_boson = _lepton;
      Particle _Wboson(24, _Wboson_mom);

      // store mt(W)
      mtw = sqrt(2*_lepton.momentum().pT()*GeV*_p_met.pT()*GeV*(1-cos(_lepton.momentum().phi()-_p_met.phi())));
      cout << " pz_neu : " << pz_neu << " pz neutrinoPz : " << pz << " mtw : " << mtw <<  endl;

      return _Wboson;
    }

    //computing z component of neutrino momentum 
    double neutrinoPz(const Particle& lepton, FourMomentum& _p_met) const {
      double pz_nu;
      double m_W = 80.399; // in GeV
      double k = (( sqr( m_W ) - sqr( lepton.momentum().mass() ) ) / 2 ) + (lepton.momentum().px() * _p_met.px() + lepton.momentum().py() * _p_met.py());
      double a = sqr ( lepton.momentum().E() )- sqr ( lepton.momentum().pz() );
      double b = -2*k*lepton.momentum().pz();
      double c = sqr( lepton.momentum().E() ) * sqr( _p_met.pT() ) - sqr( k );
      double discriminant = sqr(b) - 4 * a * c;
      double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
      if (discriminant < 0)  pz_nu = - b / (2 * a); //if the discriminant is negative
      else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
        double absquad[2];
        for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
        if (absquad[0] < absquad[1])  pz_nu = quad[0];
        else                          pz_nu = quad[1];
      }
           
      return pz_nu;
    }

    /// Neutrino reconstruction from MET via W mass constraint with met fit
    Particle returnWboson_byMET_wFit(Particle _lepton, double MET_x, double MET_y, double MET, double MET_phi) {
      Vector3 neutrino2(MET_x, MET_y, 0.0);

      FourMomentum _neutrino = reconstructNeutrino(neutrino2, _lepton);

      FourMomentum _Wboson_mom = _lepton.momentum() + _neutrino;

      _lepton_from_W_boson = _lepton;
      Particle _Wboson(24, _Wboson_mom);

      mtw = sqrt(2*_lepton.momentum().pT()*GeV*MET*GeV*(1-cos(_lepton.momentum().phi()-MET_phi)));

      return _Wboson;
    }

    FourMomentum reconstructNeutrino(const Vector3& neutrino, const Particle& lepton) {
      TLorentzVector* met = new TLorentzVector(neutrino.x(), neutrino.y(), 0., sqrt(neutrino.x()*neutrino.x() + neutrino.y()*neutrino.y()));
      TLorentzVector* lep = new TLorentzVector(lepton.momentum().x(), lepton.momentum().y(), lepton.momentum().z(), lepton.momentum().E());

      double* pz  = new double[2];
      int nsol = pz_of_W(*lep, met, pz);

      double pz_sol = 0.0;

      if (nsol > 0) {
        pz_sol = (std::fabs(pz[0]) < std::fabs(pz[1]) ? pz[0] : pz[1] );
      } else
        cout << "ArachneRecoSingleTop::reconstruct: Problem with neutrino reco" << endl;

      double met_E = sqrt((*met).X()*(*met).X() + (*met).Y()*(*met).Y() + pz_sol*pz_sol);
      return FourMomentum(met_E, (*met).X(), (*met).Y(), pz_sol);
    }

    int pz_of_W(TLorentzVector lep, TLorentzVector* met, double* pz)
    {
      const double MW = 80.4;
      int nsol = 0;

      NeutrinoFit::FullReco_PTe = lep.Pt();
      NeutrinoFit::FullReco_Pxe = lep.Px();
      NeutrinoFit::FullReco_Pye = lep.Py();
      NeutrinoFit::FullReco_MET_X = met->Px();
      NeutrinoFit::FullReco_MET_Y = met->Py();

      double MisET2 = met->Px()*met->Px() + met->Py()*met->Py();
      double mWT = sqrt(2*(lep.Pt()*sqrt(MisET2) - lep.Px()*met->Px() - lep.Py()*met->Py()));

      double PxNu_1 = 0;
      double PxNu_2 = 0;
      double PyNu_1 = 0.;
      double PyNu_2 = 0.;
      double PxNu = 0.;
      double PyNu = 0.;
      bool isComplex = false;

      if(mWT>MW){
        other_checks[0]++;
        isComplex = true;

        //cout << " before met Px : " << met->Px() << " met Py : " << met->Py() << " met E : " << met->E() << endl;
        PxNu_1 = metfit(-1,1,MW); //(Printlevel, y-solution, MW)
        PxNu_2 = metfit(-1,2,MW);

        PyNu_1 = ((MW*MW*NeutrinoFit::FullReco_Pye + 2*NeutrinoFit::FullReco_Pxe*NeutrinoFit::FullReco_Pye*PxNu_1)
            -(MW*NeutrinoFit::FullReco_PTe)*(sqrt(MW*MW + 4*NeutrinoFit::FullReco_Pxe*PxNu_1)))
          /(2*NeutrinoFit::FullReco_Pxe*NeutrinoFit::FullReco_Pxe);
        PyNu_2 = ((MW*MW*NeutrinoFit::FullReco_Pye + 2*NeutrinoFit::FullReco_Pxe*NeutrinoFit::FullReco_Pye*PxNu_2)
            +(MW*NeutrinoFit::FullReco_PTe)*(sqrt(MW*MW + 4*NeutrinoFit::FullReco_Pxe*PxNu_2)))
          /(2*NeutrinoFit::FullReco_Pxe*NeutrinoFit::FullReco_Pxe);

        double delta1 =  sqrt((PxNu_1 - NeutrinoFit::FullReco_MET_X)*(PxNu_1 - NeutrinoFit::FullReco_MET_X) +(PyNu_1 - NeutrinoFit::FullReco_MET_Y)*(PyNu_1 - NeutrinoFit::FullReco_MET_Y));

        double delta2 =  sqrt((PxNu_2 - NeutrinoFit::FullReco_MET_X)*(PxNu_2 - NeutrinoFit::FullReco_MET_X) +(PyNu_2 - NeutrinoFit::FullReco_MET_Y)*(PyNu_2 - NeutrinoFit::FullReco_MET_Y));

        if(delta1<delta2) {
          PxNu = PxNu_1;
          PyNu = PyNu_1;
        } else {
          PxNu = PxNu_2;
          PyNu = PyNu_2;
        }
      }

      double pz1,pz2;
      if( !isComplex) {

        double mu = (MW*MW)/2 + met->Px()*lep.Px() + met->Py()*lep.Py();
        double a = (mu*lep.Pz()) / (lep.E()*lep.E() - lep.Pz()*lep.Pz());
        double a2 = TMath::Power(a,2);
        double b = (TMath::Power(lep.E(),2.)*MisET2 - TMath::Power(mu,2.))/
          (TMath::Power(lep.E(),2) - TMath::Power(lep.Pz(),2));


        if(a2-b < 0) {
          Warning("ArachneReco::pz_of_W",
                  "Complex result should not happen anymore!!!");
          pz1 = a;
          pz2 = a;
          nsol = 1;
        } else {
          double root = sqrt(a2-b);
          pz1 = a + root;
          pz2 = a - root;
          nsol = 2;
        }


      } else {

        double mu = (MW*MW)/2 + PxNu*lep.Px() + PyNu*lep.Py();
        double a = (mu*lep.Pz())/(lep.E()*lep.E() - lep.Pz()*lep.Pz());
        double a2 = TMath::Power(a,2);
        double b = (TMath::Power(lep.E(),2.)*(PxNu*PxNu+PyNu*PyNu) - TMath::Power(mu,2.))/
          (TMath::Power(lep.E(),2) - TMath::Power(lep.Pz(),2));

        pz1 = a;
        pz2 = a;
        nsol = 1;
      }

      if (fabs(pz1) <= fabs(pz2) ) {
        pz[0] = pz1;
        pz[1] = pz2;
      } else {
        pz[0] = pz2;
        pz[1] = pz1;
      }

      if (isComplex) {
        met->SetPx(PxNu);
        met->SetPy(PyNu);
        //cout << " after met Px : " << met->Px() << " met Py : " << met->Py() << " met E : " << met->E() << endl;
      }


      return nsol;
    }

    double metfit(double fitterPrintLevel, int ysol, double mW)
    {
      TMinuit* minu = new TMinuit(5);

      if (ysol == 1) minu->SetFCN(NeutrinoFit::delta1fcn);
      if (ysol == 2) minu->SetFCN(NeutrinoFit::delta2fcn);


      double arglist[20];
      int ierflg = 0;

      arglist[0] = fitterPrintLevel;
      minu->mnexcm("SET PRINT", arglist, 1, ierflg);
      arglist[0] =  2.0;
      minu->mnexcm("SET STRATEGY", arglist, 1, ierflg);
      arglist[0] =  1.0;
      minu->mnexcm("CALL FCN", arglist, 1, ierflg);

      Double_t upper = 0.0;
      Double_t lower = 0.0;
      Double_t start = 0.0;

      if(NeutrinoFit::FullReco_Pxe < 0) {
        upper = - mW*mW/(4*NeutrinoFit::FullReco_Pxe) -0.01;
        lower = -9999.;
      }

      if(NeutrinoFit::FullReco_Pxe == 0) {
        upper =  9999.;
        lower = -9999.;
      }

      if(NeutrinoFit::FullReco_Pxe > 0) {
        upper = 9999.;
        lower = - mW*mW/(4*NeutrinoFit::FullReco_Pxe) +0.01;
      }

     if(NeutrinoFit::FullReco_MET_X > upper) start = upper -1;
     else if(NeutrinoFit::FullReco_MET_X < lower) start = lower + 1;
     else start = NeutrinoFit::FullReco_MET_X;



      minu->mnparm(0, "Px", start, 0.001, lower, upper, ierflg);

      arglist[0] = .5;
      minu->mnexcm("SET ERR", arglist , 1, ierflg);

      arglist[0] = 0.0;
      minu->mnexcm("SET NOW", arglist , 1, ierflg);

      ierflg = 0;
      arglist[0] = 100;
      arglist[1] = 1.;
      minu->mnexcm("SIMPLEX", arglist, 2, ierflg);
      arglist[0] = 500;
      minu->mnexcm("MIGRAD", arglist, 1, ierflg);

      minu->mnmnos() ;

      double px_fit;
      double px_fit_error;

      int rtVal = minu->GetParameter(0, px_fit, px_fit_error);
      if (rtVal < 0)
        std::cerr << "Error with parameter." << std::endl;
      if (fitterPrintLevel > 0) {
        std::cout<<"*******************************"<<std::endl;
        std::cout << "Fit Results: Px(nu) = " << px_fit
            << " +- " << px_fit_error << std::endl;
        std::cout<<"*******************************"<<std::endl;
      }

      delete minu;
      return px_fit;
    }


    
    //=====================================================================
    // fill top quark
    //=====================================================================
    // created a function to get the vector of particles and to fill the histograms
    Particle returnTopQuark(Particle _Wboson,  Jets _bjets) {
      FourMomentum _bjet_selected(0, 0, 0, 0);
      FourMomentum _topquark_selected(0, 0, 0, 0);
      
      int index = 0;
      float topquark_mass = 172.5; // in GeV
      foreach(const Jet& _bjet, _bjets) {
        FourMomentum _bjet_candidate = _bjet.momentum();
        FourMomentum _topquark_candidate = _Wboson.momentum() + _bjet.momentum();

        if(index == 0){
          _topquark_selected = _topquark_candidate;
          _bjet_selected = _bjet_candidate;
        }
        if (fabs(_topquark_candidate.mass() - topquark_mass*GeV) < fabs(_topquark_selected.mass() - topquark_mass*GeV)) {
          _bjet_selected = _bjet_candidate;
          _topquark_selected = _topquark_candidate;   
        }
        index++;
      }

      _bjet_from_top_quark = _bjet_selected;
      Particle _topquark(6, _topquark_selected);

      return _topquark;
    }
    
    //=====================================================================
    // Functions from MC_TTbar_TruthSel
    //=====================================================================
    bool SelectGoodElectrons(const vector<DressedLepton> &ele, vector<DressedLepton> &good_ele, double _pt_el, double _eta_el){
      good_ele.clear();

      foreach(DressedLepton electron, ele) { 
        if(electron.momentum().pT() <= _pt_el)
          continue;
        if(electron.momentum().abseta() >= _eta_el)
          continue;
        good_ele.push_back(electron);
      }
      
      return true;
    }
    
    bool SelectGoodMuons(const vector<DressedLepton> &mu, vector<DressedLepton> &good_mu, double _pt_mu, double _eta_mu){
      good_mu.clear();
      foreach(DressedLepton muon, mu) {
        if(muon.momentum().pT() <= _pt_mu)
          continue;
        if(muon.momentum().abseta() >= _eta_mu)
          continue;
        good_mu.push_back(muon);
      }

      return true;
    }

    bool SelectGoodJets(const Jets &alljets, Jets &bjets, Jets &ljets, Jets &all_good, double _pt_j, double _eta_j, const vector<DressedLepton> &_dressedelectrons, const vector<DressedLepton> &_dressedmuons){
      bjets.clear();
      ljets.clear();
      all_good.clear();
      
      bool _overlap = false;
      Cut b_hadrons_cut = (Cuts::abseta < 5) && (Cuts::pT >= 5.0*GeV);

      for (unsigned int i = 0; i < alljets.size(); i++) {
        const Jet& jet = alljets[i];
        
        if(jet.momentum().pT() <= _pt_j || fabs(jet.momentum().eta()) >= _eta_j)
          continue; // minimum pT and maximum eta jet cuts   

        // If dR(el,jet) < 0.4 skip the event
        foreach (const Particle& el, _dressedelectrons) {
          if (deltaR(jet, el) < 0.4) { _overlap = true; break; }
        }

        // If dR(mu,jet) < 0.4 skip the event
        foreach (const Particle& mu, _dressedmuons) {
          if (deltaR(jet, mu) < 0.4) { _overlap = true; break; }
        }
        
        for (unsigned int j = i+1; j < alljets.size(); j++) {
          const Jet& jet2 = alljets[j];
          if (deltaR(jet, jet2) < 0.5) { _overlap = true; break; }
        }

        // Remove jets too close to an electron and create vectors for bjets and ljets
        if (!_overlap) {
          all_good.push_back(jet);
          if (jet.bTagged(b_hadrons_cut) && jet.momentum().abseta()<2.5) bjets.push_back(jet);
          else ljets.push_back(jet);
        }
      }

      return _overlap;
    }

    //=====================================================================
    // BookHistograms
    //=====================================================================
    void BookHistograms() {
	  // Booking of all histograms which are interesting. In most cases, are vector of histograms is created.
	  // Information on the plots stored in MC_tchan_particle.plot (e.g. title)

      std::vector <double> ptstudy1 {0,35,50,75,100,150,200,300};
      std::vector <double> ptstudy2 {0,35,50,75,100,150,300};
      std::vector <double> ptstudy3 {30,45,60,75,100,150,300};
      std::vector <double> rapstudy1 {0,0.15,0.30,0.45,0.70,1.00,1.30,2.20};
      std::vector <double> rapstudy2 {0,1.20,1.70,2.20,2.70,3.30,4.50};
      
      //_h_neu_pz.push_back(bookHisto1D("t_neu_pz", 50, 0, 500));
      _h_t_pt_top.push_back(bookHisto1D("t_pt_top", ptstudy1));
      _h_t_pt_top_norm.push_back(bookHisto1D("t_pt_top_norm", ptstudy1));
      _h_t_pt_antitop.push_back(bookHisto1D("t_pt_antitop", ptstudy2));
      _h_t_pt_antitop_norm.push_back(bookHisto1D("t_pt_antitop_norm", ptstudy2));
      _h_t_rap_top.push_back(bookHisto1D("t_rap_top", rapstudy1));
      _h_t_rap_top_norm.push_back(bookHisto1D("t_rap_top_norm", rapstudy1));
      _h_t_rap_antitop.push_back(bookHisto1D("t_rap_antitop", rapstudy1));
      _h_t_rap_antitop_norm.push_back(bookHisto1D("t_rap_antitop_norm", rapstudy1));
      _h_ljet_pt_top.push_back(bookHisto1D("ljet_pt_top", ptstudy3));
      _h_ljet_pt_top_norm.push_back(bookHisto1D("ljet_pt_top_norm", ptstudy3));
      _h_ljet_pt_antitop.push_back(bookHisto1D("ljet_pt_antitop", ptstudy3));
      _h_ljet_pt_antitop_norm.push_back(bookHisto1D("ljet_pt_antitop_norm", ptstudy3));
      _h_ljet_rap_top.push_back(bookHisto1D("ljet_rap_top", rapstudy2));
      _h_ljet_rap_top_norm.push_back(bookHisto1D("ljet_rap_top_norm", rapstudy2));
      _h_ljet_rap_antitop.push_back(bookHisto1D("ljet_rap_antitop", rapstudy2));
      _h_ljet_rap_antitop_norm.push_back(bookHisto1D("ljet_rap_antitop_norm", rapstudy2));

      _h_weight = bookHisto1D("weight", 20, -2, 2);
    }

    //@}

  private:

    double _weight;
    int cutflow[9];
    int cutflow_ele[9];
    int cutflow_mu[9];
    int eventcount;
//    std::string cutnames[9];
    int other_checks[1];
    
    double _el_weights, _mu_weights;

    FourMomentum _lepton_from_W_boson;
    FourMomentum _bjet_from_top_quark;
    FourMomentum _spectjet_from_top_quark;
    
    double mtw;
    
    vector<DressedLepton> _dressedelectrons;
    vector<DressedLepton> _vetodressedelectrons;
    vector<DressedLepton> _ewdressedelectrons;
    vector<DressedLepton> _dressedmuons;
    vector<DressedLepton> _vetodressedmuons;
    vector<DressedLepton> _ewdressedmuons;
    Particles _neutrinos;
    Particles _photons;

    /// @name Histograms
    //@{

    Histo1DPtr _h_weight;
    //vector<Histo1DPtr> _h_neu_pz;

    // top histograms
    vector<Histo1DPtr> _h_t_pt_top;
    vector<Histo1DPtr> _h_t_pt_top_norm;
    vector<Histo1DPtr> _h_t_rap_top;
    vector<Histo1DPtr> _h_t_rap_top_norm;
    vector<Histo1DPtr> _h_ljet_pt_top;
    vector<Histo1DPtr> _h_ljet_pt_top_norm;
    vector<Histo1DPtr> _h_ljet_rap_top;
    vector<Histo1DPtr> _h_ljet_rap_top_norm;
    
    // antitop histograms
    vector<Histo1DPtr> _h_t_pt_antitop;
    vector<Histo1DPtr> _h_t_pt_antitop_norm;
    vector<Histo1DPtr> _h_t_rap_antitop;
    vector<Histo1DPtr> _h_t_rap_antitop_norm;
    vector<Histo1DPtr> _h_ljet_pt_antitop;
    vector<Histo1DPtr> _h_ljet_pt_antitop_norm;
    vector<Histo1DPtr> _h_ljet_rap_antitop;
    vector<Histo1DPtr> _h_ljet_rap_antitop_norm;

    // new histo pointers
    Histo1DPtr _h_AbsPtclDiffXsecTPt[2], _h_AbsPtclDiffXsecTY[2], _h_NrmPtclDiffXsecTPt[2], _h_NrmPtclDiffXsecTY[2];
    Histo1DPtr _h_AbsPtclDiffXsecJPt[2], _h_AbsPtclDiffXsecJY[2], _h_NrmPtclDiffXsecJPt[2], _h_NrmPtclDiffXsecJY[2];
    //Histo1DPtr _h_AbsPrtnDiffXsecTPt[2], _h_AbsPrtnDiffXsecTY[2], _h_NrmPrtnDiffXsecTPt[2], _h_NrmPrtnDiffXsecTY[2];

    //@}
  };

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<ATLAS_2017_I1512776> plugin_ATLAS_2017_I1512776;

}
