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
#include "Rivet/Cuts.hh"
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

#include <TMath.h>
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

  class MC_tchan_particle : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_tchan_particle()
      : Analysis("MC_tchan_particle")
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
    }

    //=====================================================================
    // Perform the per-event analysis
    //=====================================================================
    void analyze(const Event& event) {
      // bool _debug = false;
      _weight = event.weight();
      _h_weight->fill(_weight, _weight);
      MSG_DEBUG("mc weight: " << _weight);
      
//      eventcount++;

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
      std::cout << "\tcutflow[0]++ \n";
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
      std::cout << "\tcutflow[1]++ \n";
      if (good_leptons.size() != 1) {vetoEvent; }
      // rejection if veto particles are present
      cutflow[2]++;

      if (good_electrons.size() == 1 && (veto_muons.size() != 0 || veto_electrons.size() > 1)) { vetoEvent; }
      if (good_muons.size() == 1 && (veto_electrons.size() != 0 || veto_muons.size() > 1)) { vetoEvent; }

      // separate all events into top and antitop channel
      bool isTop = true;
      if (good_leptons[0].pdgId() < 0) {
        isTop = true;   // If there is an anti-lepton, the decaying top is a top
      } else {
        isTop = false;  // If there is a lepton, the decaying top is an anti-top
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

      
      // fill W boson
      // Reconstruction of the W boson
//      Particle _W_boson = returnWboson_byMET(good_leptons[0], _p_met.x(), _p_met.y(), _p_met.pT(), _p_met.phi());
      
      // fill top quark
      // Reconstruction of the top quark
      // additionally stores the properties of the b-jet associated to the top
//      Particle _top_quark = returnTopQuark(_W_boson, good_bjets);
      
      // veto event if eTmiss cut not passed
      cutflow[6]++;
//      if ( _p_met.pT()/GeV < 30.0 ) { vetoEvent; }
      
      // veto event if mtw cut not passed
      cutflow[7]++;
      std::cout << "\tcutflow[7]++ \n";

      Vector3 bjet_v3=good_bjets[0].momentum().vector3();
      FourMomentum bjet=good_bjets[0].momentum();
      // This is needed to get same acceptance as using ROOT
      bjet.setE(sqrt(bjet_v3.x()*bjet_v3.x()+bjet_v3.y()*bjet_v3.y()+bjet_v3.z()*bjet_v3.z()));

      FourMomentum lepton_b= good_leptons[0].momentum()+bjet;

      if((lepton_b.mass())>=160.0) {vetoEvent;}
      // -----------------------------------------------------------------------------------------------
      
      /// After this point no additional cuts are applied and all histograms (except for the total weight in the first lines) are filled
      cutflow[8]++;
      std::cout << "\tcutflow[8]++ \n";
      if(good_electrons.size() == 1) cutflow_ele[8]++;
      if(good_muons.size() == 1) cutflow_mu[8]++;
      
//      std::ofstream outfile;
//      outfile.open("test.txt", std::ios_base::app);
//      outfile <<eventcount<<" "<<std::setprecision(5)<<good_leptons[0].momentum().pT()<<" "<<good_bjets[0].momentum().pT()<<"\n";
//      outfile.close();
      /*
      // filling of total weight in the channels and charged particles
      if(isTop) {
        _h_weight_top->fill(_weight, _weight);
        fillPropertiesFromCollection(ChargedParticleCollection, isTop);
      } else {
        _h_weight_antitop->fill(_weight, _weight);
        fillPropertiesFromCollection(ChargedParticleCollection, isTop);
      }
      
      // storing the total weight of both channels
      if(good_electrons.size() == 1) {
        _el_weights += _weight;
      } else {
        _mu_weights += _weight;
      }

      // fill lepton properties and MET
      if(isTop) {
        fillPropertiesFromCollection(good_electrons, _h_electrons_top);
        fillPropertiesFromCollection(veto_electrons, _h_veto_electrons_top);
        fillPropertiesFromCollection(ew_electrons, _h_ew_electrons_top);
        fillPropertiesFromCollection(good_muons, _h_muons_top);
        fillPropertiesFromCollection(veto_muons, _h_veto_muons_top);
        fillPropertiesFromCollection(ew_muons, _h_ew_muons_top);
        fillPropertiesFromCollection(good_neutrinos_prompt, _h_prompt_neutrinos_top);
        fillPropertiesFromCollection(good_neutrinos_met, _h_met_neutrinos_top);
        
        // fill leading lepton properties
        fillLeadingParticleProperties(good_electrons, _h_electron_1_top);
        fillLeadingParticleProperties(good_muons, _h_muon_1_top);
        fillLeadingParticleProperties(good_neutrinos_prompt, _h_prompt_neutrino_1_top);
        fillLeadingParticleProperties(good_neutrinos_met, _h_met_neutrino_1_top);
      } else {
        fillPropertiesFromCollection(good_electrons, _h_electrons_antitop);
        fillPropertiesFromCollection(veto_electrons, _h_veto_electrons_antitop);
        fillPropertiesFromCollection(ew_electrons, _h_ew_electrons_antitop);
        fillPropertiesFromCollection(good_muons, _h_muons_antitop);
        fillPropertiesFromCollection(veto_muons, _h_veto_muons_antitop);
        fillPropertiesFromCollection(ew_muons, _h_ew_muons_antitop);
        fillPropertiesFromCollection(good_neutrinos_prompt, _h_prompt_neutrinos_antitop);
        fillPropertiesFromCollection(good_neutrinos_met, _h_met_neutrinos_antitop);
        
        // fill leading lepton properties
        fillLeadingParticleProperties(good_electrons, _h_electron_1_antitop);
        fillLeadingParticleProperties(good_muons, _h_muon_1_antitop);
        fillLeadingParticleProperties(good_neutrinos_prompt, _h_prompt_neutrino_1_antitop);
        fillLeadingParticleProperties(good_neutrinos_met, _h_met_neutrino_1_antitop);
      }
      fillMET(_p_met, isTop);
      
      // fill jet properties
      fillPropertiesFromCollection(good_jets, isTop);
      fillJetsFlavours(good_bjets, good_ljets, isTop);
      fillBHadrons(B_hadrons, isTop, 2.5, 5.0);
      
      // fill W and top properties
      fillWboson(_W_boson, isTop);
      fillTopQuark(_W_boson, good_bjets, isTop);

      // select the spectator quark (only one remaining after cuts)
      _spectjet_from_top_quark = good_ljets[0];
      FillLJet(isTop);

      // get here the 4-momenta for the calc. of cos theta star
      // Some angular correlations between the interesting particles in the final state
      FillAngularDistributions(_top_quark.momentum(), _W_boson.momentum(), isTop);
      */
    } // end of analyze()
    
    // Normalise histograms etc., after the run
    void finalize() {
      double lum = 20276.9;
      
      cout << "============================= tchan particle ================================" << endl;
      cout << "crossSection " << crossSection() << " sumOfWeights " << sumOfWeights() << endl;
      
      double scale_fac = crossSection() * lum / sumOfWeights();
      
      cout << "No. of electron events: " << _el_weights * scale_fac << endl;
      cout << "No. of muon events: " << _mu_weights * scale_fac << endl;
      cout << "=============================================================================" << endl;
      for (int i = 0; i < 9; i++) {
        cout << "Cutflow at step " << i << ": " << cutflow[i]<<" Afid: "<<(float)cutflow[i]/(float)cutflow[0] << endl;
      }
      for (int i = 0; i < 1; i++) {
        cout << "No. of events with M_T(W) above M(W): " << other_checks[i] << endl;
      }
      cout << "=============================================================================" << endl;
      /*
      // Scale all histograms by their cross section
      scale(_h_weight_top, scale_fac);
      scale(_h_weight_antitop, scale_fac);
      
      for (unsigned int i=0; i< _h_charged_top.size(); i++) {
        scale(_h_charged_top[i], scale_fac);
        scale(_h_charged_antitop[i], scale_fac);
      }
      
      for (unsigned int i=0; i< _h_B_hadrons_top.size(); i++) {
        scale(_h_B_hadrons_top[i], scale_fac);
        scale(_h_B_hadrons_antitop[i], scale_fac);
      }
      
      for (unsigned int i=0; i< _h_jets_top.size(); i++) {
        scale(_h_jets_top[i], scale_fac);
        scale(_h_jets_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_ljets_top.size(); i++) {
        scale(_h_ljets_top[i], scale_fac);
        scale(_h_ljets_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_bjets_top.size(); i++) {
        scale(_h_bjets_top[i], scale_fac);
        scale(_h_bjets_antitop[i], scale_fac);
      }
      
      for (unsigned int i=0; i< _h_electrons_top.size(); i++) {
        scale(_h_electrons_top[i], scale_fac);
        scale(_h_electrons_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_muons_top.size(); i++) {
        scale(_h_muons_top[i], scale_fac);
        scale(_h_muons_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_taus_top.size(); i++) {
        scale(_h_taus_top[i], scale_fac);
        scale(_h_taus_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_prompt_neutrinos_top.size(); i++) {
        scale(_h_prompt_neutrinos_top[i], scale_fac);
        scale(_h_prompt_neutrinos_antitop[i], scale_fac);
        scale(_h_met_neutrinos_top[i], scale_fac);
        scale(_h_met_neutrinos_antitop[i], scale_fac);
      }
      
      for (unsigned int i=0; i< _h_veto_electrons_top.size(); i++) {
        scale(_h_veto_electrons_top[i], scale_fac);
        scale(_h_veto_electrons_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_veto_muons_top.size(); i++) {
        scale(_h_veto_muons_top[i], scale_fac);
        scale(_h_veto_muons_antitop[i], scale_fac);
      }
      
      for (unsigned int i=0; i< _h_ew_electrons_top.size(); i++) {
        scale(_h_ew_electrons_top[i], scale_fac);
        scale(_h_ew_electrons_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_ew_muons_top.size(); i++) {
        scale(_h_ew_muons_top[i], scale_fac);
        scale(_h_ew_muons_antitop[i], scale_fac);
      }
      
      for (unsigned int i=0; i< _h_electron_1_top.size(); i++) {
        scale(_h_electron_1_top[i], scale_fac);
        scale(_h_electron_1_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_muon_1_top.size(); i++) {
        scale(_h_muon_1_top[i], scale_fac);
        scale(_h_muon_1_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_prompt_neutrino_1_top.size(); i++) {
        scale(_h_prompt_neutrino_1_top[i], scale_fac);
        scale(_h_prompt_neutrino_1_antitop[i], scale_fac);
        scale(_h_met_neutrino_1_top[i], scale_fac);
        scale(_h_met_neutrino_1_antitop[i], scale_fac);
      }
      
      for (unsigned int i=0; i< _h_t_top.size(); i++) {
        scale(_h_t_top[i], scale_fac);
        scale(_h_t_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_W_top.size(); i++) {
        scale(_h_W_top[i], scale_fac);
        scale(_h_W_antitop[i], scale_fac);
      }
      
      for (unsigned int i=0; i< _h_MET_top.size(); i++) {
        scale(_h_MET_top[i], scale_fac);
        scale(_h_MET_antitop[i], scale_fac);
      }
      
      scale(_h_cosTheta_S_top, scale_fac);
      scale(_h_cosTheta_N_top, scale_fac);
      scale(_h_cosTheta_T_top, scale_fac);
      scale(_h_cosTheta_X_top, scale_fac);
      scale(_h_cosTheta_top, scale_fac);
      scale(_h_Phi_S_top, scale_fac);
      
      scale(_h_cosTheta_S_antitop, scale_fac);
      scale(_h_cosTheta_N_antitop, scale_fac);
      scale(_h_cosTheta_T_antitop, scale_fac);
      scale(_h_cosTheta_X_antitop, scale_fac);
      scale(_h_cosTheta_antitop, scale_fac);
      scale(_h_Phi_S_antitop, scale_fac);
      
      for (unsigned int i=0; i< _h_t_pt_top.size(); i++) {
        scale(_h_t_pt_top[i], scale_fac);
        scale(_h_t_pt_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_t_rap_top.size(); i++) {
        scale(_h_t_rap_top[i], scale_fac);
        scale(_h_t_rap_antitop[i], scale_fac);
      }
      
      for (unsigned int i=0; i< _h_ljet_pt_top.size(); i++) {
        scale(_h_ljet_pt_top[i], scale_fac);
        scale(_h_ljet_pt_antitop[i], scale_fac);
      }
      for (unsigned int i=0; i< _h_ljet_rap_top.size(); i++) {
        scale(_h_ljet_rap_top[i], scale_fac);
        scale(_h_ljet_rap_antitop[i], scale_fac);
      }
    */
    }


    //=====================================================================
    // FillAngularDistributions
    //=====================================================================
    void FillAngularDistributions(FourMomentum top, FourMomentum W, bool isTop) {
      // need to go into the W rest frame here

//       PrintFourMomentum(0, "spectator jet", _spectjet_from_top_quark);
//       PrintFourMomentum(0, "top quark", top);
//       PrintFourMomentum(0, "W boson",   W);
//       PrintFourMomentum(0, "lepton from W boson", _lepton_from_W_boson);
//       PrintFourMomentum(0, "b-jet from W boson", _bjet_from_top_quark);

      /*
      const Vector3 boostVec(W.px(), W.py(), W.pz());

      LorentzTransform boostTrafo; 
      boostTrafo.setBoost(-W.boostVector());

      const FourMomentum boostW   = boostTrafo.transform(W);
      //PrintFourMomentum(0,   "boosteW",   boostW);
      const FourMomentum boostTop = boostTrafo.transform(top);
      //PrintFourMomentum(0,  "BoostedTop", boostTop);
      const FourMomentum boostLep = boostTrafo.transform(lep); 
      //PrintFourMomentum(0, "BoostedLep",  boostLep);

      double cosThetaStar = cos(boostLep.angle(-boostTop));
      // MSG_INFO("cos(theta*): " << cosThetaStar);
      */

      // boost to top rest frame
      LorentzTransform boostTraf_Top;
      boostTraf_Top.setBoost(-top.boostVector());
      const FourMomentum SpectatorJetBoostedToTopCM = boostTraf_Top.transform(_spectjet_from_top_quark); // boost spectator jet to top rest frame
      const FourMomentum WBoostedToTopCM = boostTraf_Top.transform(W); // boost W boson to top rest frame
      const FourMomentum LeptonBoostedToTopCM = boostTraf_Top.transform(_lepton_from_W_boson); // boost lepton to top rest frame
      const FourMomentum bBoostedToTopCM = boostTraf_Top.transform(_bjet_from_top_quark); // boost b-jet to top rest frame

      // boost to W rest frame
      LorentzTransform boostTraf_W;
      boostTraf_W.setBoost(-WBoostedToTopCM.boostVector());
      const FourMomentum bJetBoostedToWCM = boostTraf_W.transform(bBoostedToTopCM); // boost b-jet from Top rest frame to W rest frame
      const FourMomentum LeptonBoostedToWCM = boostTraf_W.transform(LeptonBoostedToTopCM); // boost lepton from Top rest frame to W rest frame

      // angular distribution to extract the W boson helicity
      double cosTheta_S = cos(WBoostedToTopCM.angle(LeptonBoostedToWCM));
      // MSG_INFO("cos(theta*): " << cosTheta_S);
      if(isTop) {
        _h_cosTheta_S_top->fill(cosTheta_S, _weight);
      } else {
        _h_cosTheta_S_antitop->fill(cosTheta_S, _weight);
      }

      // define the Normal (N) and Transverse (T) direction
      // N = s x q (q ... W momentum direction)
      // T = q x N (s ... spectator jet mom direction)
      const Vector3 S = SpectatorJetBoostedToTopCM.vector3();
      const Vector3 q = WBoostedToTopCM.vector3();
      const Vector3 leptonInTopCM = LeptonBoostedToTopCM.vector3();
      const Vector3 l = LeptonBoostedToWCM.vector3();
      const Vector3 N = S.cross(q);
      const Vector3 T = q.cross(N);
      double cosTheta_N = cos(l.angle(N));
      double cosTheta_T = cos(l.angle(T));
      // MSG_INFO("cos(theta)_N: " << cosTheta_N);
      // MSG_INFO("cos(theta)_T: " << cosTheta_T);
      if(isTop) {
        _h_cosTheta_N_top->fill(cosTheta_N, _weight);
        _h_cosTheta_T_top->fill(cosTheta_T, _weight);
      } else {
        _h_cosTheta_N_antitop->fill(cosTheta_N, _weight);
        _h_cosTheta_T_antitop->fill(cosTheta_T, _weight);
      }
      
      // phi (from the two/three angle analyses)
      double Phi_S = atan2(cosTheta_N, -cosTheta_T);
      // MSG_INFO("Phi_S: " << Phi_S << " rad");
      while (Phi_S<0.) { Phi_S += 2*pi; }
      while (Phi_S>twopi) { Phi_S -= 2*pi; }
      Phi_S = Phi_S/pi;
      // MSG_INFO("Phi_S: " << Phi_S << " pi");
      if(isTop) {
        _h_Phi_S_top->fill(Phi_S, _weight);
      } else {
        _h_Phi_S_antitop->fill(Phi_S, _weight);
      }

      // cos(theta) (from the two/three angle analyses)
      double cosTheta = cos(S.angle(q));
      // MSG_INFO("cos(theta): " << cosTheta);
      if(isTop) {
        _h_cosTheta_top->fill(cosTheta, _weight);
      } else {
        _h_cosTheta_antitop->fill(cosTheta, _weight);
      }
      
      // angular distribution to extract the spin analyser times the top polarization, i.e. alpha_X\B7P.
      // If one considers the lepton in the top rest frame as the "spin analyser" then alpha_X=1
      double cosTheta_X = cos(S.angle(leptonInTopCM));
      // MSG_INFO("cos(theta)_X: " << cosTheta_X);
      if(isTop) {
        _h_cosTheta_X_top->fill(cosTheta_X, _weight);
      } else {
        _h_cosTheta_X_antitop->fill(cosTheta_X, _weight);
      }
    }
    
    // Filling of the light quark jet histograms
    void FillLJet(bool isTop) {
      if(isTop) {
        for(unsigned int i=0; i < _h_ljet_pt_top.size(); i++) {
          _h_ljet_pt_top[i]->fill(_spectjet_from_top_quark.pT()*GeV, _weight);
        }
        for(unsigned int i=0; i < _h_ljet_rap_top.size(); i++) {
          _h_ljet_rap_top[i]->fill(_spectjet_from_top_quark.absrap(), _weight);
        }
      } else {
        for(unsigned int i=0; i < _h_ljet_pt_antitop.size(); i++) {
          _h_ljet_pt_antitop[i]->fill(_spectjet_from_top_quark.pT()*GeV, _weight);
        }
        for(unsigned int i=0; i < _h_ljet_rap_antitop.size(); i++) {
          _h_ljet_rap_antitop[i]->fill(_spectjet_from_top_quark.absrap(), _weight);
        }
      }
    }

    //=====================================================================
    // PrintFourMomentum
    //=====================================================================
    void PrintFourMomentum(string jetname, Jets _jets) {
      int index = 0;
      
      foreach(const Jet& j, _jets) {
        PrintFourMomentum(index, jetname, j.momentum());
        index++;
      }
    }
    
    //=====================================================================
    // PrintFourMomentum
    //=====================================================================
    void PrintFourMomentum(vector<HepMC::GenParticle*> _partvect) {
      int index = 0;
      
      foreach (HepMC::GenParticle* _p, _partvect) {
        PrintFourMomentum(index, Particle(*_p));
        index++;
      }
    }

    //=====================================================================
    // PrintFourMomentum
    //=====================================================================
    void PrintFourMomentum(vector<Particle> _partvect) {
      int index = 0;
      
      foreach (Particle _p, _partvect) { 
        PrintFourMomentum(index, _p);
        index++;
      }
    }

    //=====================================================================
    // Mass
    //=====================================================================
    float Mass(Particle _p) { return Mass(_p.momentum()); }

    //=====================================================================
    // Mass
    //=====================================================================
    float Mass(FourMomentum FourMom) {
      float mass = sqrt(FourMom.mass2())/GeV;
      if (Rivet::isZero(mass, 1E-9)) { mass = 0.; }
      if (std::isnan(mass)) { mass = 0.; }
      return mass;
    }

    //=====================================================================
    // PrintFourMomentum
    //=====================================================================
    void PrintFourMomentum(int index, Particle _p) {
      string partname = "";
      if (fabs(_p.pdgId())==11) { partname = "electron"; }
      else if (fabs(_p.pdgId())==13) { partname = "muon"; }
      else if (fabs(_p.pdgId())==15) { partname = "tau"; }
      else if (fabs(_p.pdgId())==12 || fabs(_p.pdgId())==14 || fabs(_p.pdgId())==16) { partname = "neutrino"; }
      else if (fabs(_p.pdgId())==24) { partname = "W boson"; }
      else if (fabs(_p.pdgId())==5) { partname = "b quark"; }
      else if (fabs(_p.pdgId())==6) { partname = "top quark"; }
      else { partname = patch::to_string(_p.pdgId()); }

      // don't use _p.momentum().mass(), use this instead:      
      PrintFourMomentum(index, partname, _p.momentum());
    }

    //=====================================================================
    // PrintFourMomentum
    //=====================================================================
    void PrintFourMomentum(int index, string partname, FourMomentum FourMom) {
      MSG_INFO(index << ". (" << partname << ") " <<
 	       "pT = " << FourMom.pT()/GeV << " GeV, " <<
 	       "m = " << Mass(FourMom) << " GeV, " <<
	       "E = " << FourMom.E()/GeV << " GeV, " <<
 	       "ET = " << FourMom.Et()/GeV << " GeV, " <<
	       "p = [" << FourMom.px()/GeV << ", " << FourMom.py()/GeV << ", " << FourMom.pz()/GeV << "] GeV, " <<
 	       "\u03B7 = " << FourMom.eta() << ", "
 	       "\u0278 = " << FourMom.phi()<< " rad");
    }

    //=====================================================================
    // FillFourMomentum
    //=====================================================================
    void FillFourMomentum(Particle _p, vector<Histo1DPtr> _histo, int index=1, bool fillMass=false) {
      _histo.at(index)->fill(_p.momentum().E()/GeV, _weight);
      _histo.at(index+1)->fill(_p.momentum().Et()/GeV, _weight);
      _histo.at(index+2)->fill(_p.momentum().pT()/GeV, _weight);
      _histo.at(index+3)->fill(_p.momentum().pz()/GeV, _weight);
      _histo.at(index+4)->fill(_p.momentum().eta(), _weight);
      _histo.at(index+5)->fill(_p.momentum().phi(), _weight);
      if (fillMass) { _histo.at(index+6)->fill(Mass(_p), _weight); }
    }

    //=====================================================================
    // FillFourMomentum
    //=====================================================================
    void FillFourMomentum(const Jet& _j, vector<Histo1DPtr> _histo, int index=1, bool fillMass=false) {
      _histo.at(index)->fill(_j.momentum().E()/GeV, _weight);
      _histo.at(index+1)->fill(_j.momentum().Et()/GeV, _weight);
      _histo.at(index+2)->fill(_j.momentum().pT()/GeV, _weight);
      _histo.at(index+3)->fill(_j.momentum().pz()/GeV, _weight);
      _histo.at(index+4)->fill(_j.momentum().eta(), _weight);
      _histo.at(index+5)->fill(_j.momentum().phi(), _weight);
      if (fillMass) { _histo.at(index+6)->fill(Mass(_j.momentum()), _weight); }
    }

    //=====================================================================
    // fill properties from Collection - ChargedFinalState
    //=====================================================================
    // only used for the charged particles without interfering with the t-channel analysis, therefore no changes have been applied.
    vector<FourMomentum> fillPropertiesFromCollection(const ChargedFinalState& Collection, bool isTop, float etaCut=5.0, float ptCut=0.0) {
      vector<FourMomentum> _charged;
      
      foreach(const Particle& _p, Collection.particlesByPt()) {
        FourMomentum FourMom(0., 0., 0., 0.);
        FourMom = FourMom + _p.momentum();
        if (FourMom.pT() > ptCut && fabs(FourMom.eta()) < etaCut) {
          _charged.push_back(FourMom);
          if(isTop) {
            _h_charged_top.at(1)->fill(FourMom.E()/GeV, _weight);
            _h_charged_top.at(2)->fill(FourMom.Et()/GeV, _weight);
            _h_charged_top.at(3)->fill(FourMom.pT()/GeV, _weight);
            _h_charged_top.at(4)->fill(FourMom.pz()/GeV, _weight);
            _h_charged_top.at(5)->fill(FourMom.eta(), _weight);
            _h_charged_top.at(6)->fill(FourMom.phi(), _weight);
            _h_charged_top.at(7)->fill(Mass(FourMom), _weight);
            // _h_charged.at(8)->fill(sqrt(FourMom.mass2()+pow(FourMom.px(),2)+pow(FourMom.py(),2)), _weight);
          } else {
            _h_charged_antitop.at(1)->fill(FourMom.E()/GeV, _weight);
            _h_charged_antitop.at(2)->fill(FourMom.Et()/GeV, _weight);
            _h_charged_antitop.at(3)->fill(FourMom.pT()/GeV, _weight);
            _h_charged_antitop.at(4)->fill(FourMom.pz()/GeV, _weight);
            _h_charged_antitop.at(5)->fill(FourMom.eta(), _weight);
            _h_charged_antitop.at(6)->fill(FourMom.phi(), _weight);
            _h_charged_antitop.at(7)->fill(Mass(FourMom), _weight);
            // _h_charged.at(8)->fill(sqrt(FourMom.mass2()+pow(FourMom.px(),2)+pow(FourMom.py(),2)), _weight);
          }
        }
      }
      
      MSG_DEBUG("Number of charged particles = " << _charged.size());
      if(isTop) {
        _h_charged_top.at(0)->fill(Collection.size(),_weight);
      } else {
        _h_charged_antitop.at(0)->fill(Collection.size(),_weight);
      }
      // _h_charged.at(0)->fill(_charged.size(), _weight);
      return _charged;
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
    
    void fillPropertiesFromCollection(const vector<Particle> Collection, vector<Histo1DPtr> _histo, float etaCut=5.0, float ptCut=0.0) {
      vector<Particle> _partvect;
      Particle _p;
      
      for(unsigned int i = 0; i < Collection.size(); i++) {
        _p = Collection.at(i);
        if (_p.momentum().pT() > ptCut && fabs(_p.momentum().eta()) < etaCut) {
          _partvect.push_back(_p);
          FillFourMomentum(_p, _histo);
        }
      }
      
      _histo.at(0)->fill(_partvect.size(), _weight);
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
    
    void fillPropertiesFromCollection(const vector<Particle> Collection, vector<Histo1DPtr> _histo, int pdfId, float etaCut=5.0, float ptCut=0.0) {
      vector<Particle> _partvect;
      Particle _p;
      
      for(unsigned int i = 0; i < Collection.size(); ++i) {
        _p = Collection.at(i);
        if(abs(_p.pdgId())!=pdfId) continue;
        if (_p.momentum().pT() > ptCut && fabs(_p.momentum().eta()) < etaCut) {
          _partvect.push_back(_p);
          FillFourMomentum(_p, _histo);
        }
      }
      
      _histo.at(0)->fill(_partvect.size(), _weight);
    }

    //=====================================================================
    // fill properties of the leading particle - vector<Particle>
    //=====================================================================
    void fillLeadingParticleProperties(vector<Particle> _partvect, vector<Histo1DPtr> _histo) {
      if (_partvect.size()>0) {
        FillFourMomentum(_partvect[0], _histo, 0);
      }
    }
    
    //=====================================================================
    // fill MET values to histograms
    //=====================================================================
    void fillMET(FourMomentum _p_met, bool isTop) {
      if(isTop) {
        _h_MET_top.at(0)->fill(_p_met.pT()/GeV, _weight);
        _h_MET_top.at(1)->fill(_p_met.eta(), _weight);
        _h_MET_top.at(2)->fill(_p_met.phi(), _weight);
      } else {
        _h_MET_antitop.at(0)->fill(_p_met.pT()/GeV, _weight);
        _h_MET_antitop.at(1)->fill(_p_met.eta(), _weight);
        _h_MET_antitop.at(2)->fill(_p_met.phi(), _weight);
      }
    }

    //=====================================================================
    // fill properties from Collection - FastJets
    //=====================================================================    
    void fillPropertiesFromCollection(const Jets& Collection, bool isTop) {
      Jets _jets;
      for (unsigned int i = 0; i < Collection.size(); ++i) {
        Jet _jet = Collection.at(i);
        
        _jets.push_back(_jet);
        if(isTop) {
          FillFourMomentum(_jet, _h_jets_top, 1, true);
        } else {
          FillFourMomentum(_jet, _h_jets_antitop, 1, true);
        }
      }
      
      if(isTop) {
        _h_jets_top.at(0)->fill(_jets.size(), _weight);
      } else {
        _h_jets_antitop.at(0)->fill(_jets.size(), _weight);
      }
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
    
    void fillBHadrons(const vector<const HepMC::GenParticle*> B_hadrons, bool isTop, float etaCut=4.5, float ptCut=5.0) {
      for (size_t i = 0; i < B_hadrons.size(); i++) {
        const HepMC::GenParticle* _p = B_hadrons.at(i);

        // If the particle is a B-hadron and has pT > ptCut GeV, add it to the B-hadrons vector
        if(PID::isHadron(_p->pdg_id()) && PID::hasBottom(_p->pdg_id()) && fabs(_p->momentum().eta()) < etaCut && _p->momentum().perp() > ptCut) {
          Particle B_hadron = Particle(*_p);
          if(isTop) {
            FillFourMomentum(B_hadron, _h_B_hadrons_top, 1, true);
            _h_B_hadrons_top.at(8)->fill(sqrt(B_hadron.momentum().mass2()+pow(B_hadron.momentum().px(), 2)+pow(B_hadron.momentum().py(),2))/GeV, _weight);
          } else {
            FillFourMomentum(B_hadron, _h_B_hadrons_antitop, 1, true);
            _h_B_hadrons_antitop.at(8)->fill(sqrt(B_hadron.momentum().mass2()+pow(B_hadron.momentum().px(), 2)+pow(B_hadron.momentum().py(),2))/GeV, _weight);
          }
        } else {
          cout << "Different cuts in fillBHadrons applied? Should not be the case." << endl;
        }
      }
      
      if(isTop) {
        _h_B_hadrons_top.at(0)->fill(B_hadrons.size(), _weight);
      } else {
        _h_B_hadrons_antitop.at(0)->fill(B_hadrons.size(), _weight);
      }
    }

    //=====================================================================
    // fill jets flavours
    //=====================================================================
    // created a function to get the vector of particles and to fill the histograms
    void fillJetsFlavours(Jets& good_bjets, Jets& good_ljets, bool isTop) {
      if(isTop) {
        foreach(const Jet& _j, good_bjets) {
          FillFourMomentum(_j, _h_bjets_top, 1, true);
        }
        foreach(const Jet& _j, good_ljets) {
          FillFourMomentum(_j, _h_ljets_top, 1, true);
        }
        _h_ljets_top.at(0)->fill(good_ljets.size(), _weight);
        _h_bjets_top.at(0)->fill(good_bjets.size(), _weight);
      } else {
        foreach(const Jet& _j, good_bjets) {
          FillFourMomentum(_j, _h_bjets_antitop, 1, true);
        }
        foreach(const Jet& _j, good_ljets) {
          FillFourMomentum(_j, _h_ljets_antitop, 1, true);
        }
        _h_ljets_antitop.at(0)->fill(good_ljets.size(), _weight);
        _h_bjets_antitop.at(0)->fill(good_bjets.size(), _weight);
      }
    }

    //=====================================================================
    // fill W boson
    //=====================================================================
    /// Function not used! This is the old neutrino search code without reconstruction from MET
    Particle returnWboson(vector<Particle> _leptons,  vector<Particle> _neutrinos) {
      FourMomentum _lepton_selected(0, 0, 0, 0);
      FourMomentum _neutrino_selected(0, 0, 0, 0);
      FourMomentum _Wboson_selected(0, 0, 0, 0);

      int index = 0;
      float Wboson_mass = 80.399*GeV; // in GeV
      foreach (Particle _neutrino, _neutrinos) {
        foreach (Particle _lepton, _leptons) {
          FourMomentum _lepton_candidate   = _lepton.momentum();
          FourMomentum _neutrino_candidate = _neutrino.momentum();
          FourMomentum _Wboson_candidate   = _lepton.momentum() + _neutrino.momentum();

          if(index  == 0) {
            _lepton_selected   = _lepton_candidate;
            _neutrino_selected = _neutrino_candidate;
            _Wboson_selected   = _Wboson_candidate;
          } else if(fabs(_Wboson_candidate.mass() - Wboson_mass) < fabs(_Wboson_selected.mass() - Wboson_mass)) {	    
            _lepton_selected   = _lepton_candidate;
            _neutrino_selected = _neutrino_candidate;
            _Wboson_selected   = _Wboson_candidate;	    
          }
          index++;
        }
      }
      
      _lepton_from_W_boson = _lepton_selected; 
      Particle _Wboson(24, _Wboson_selected);

      // store mt(W)
      mtw = sqrt(2*_lepton_selected.pT()*GeV*_neutrino_selected.pT()*GeV*(1-cos(_lepton_selected.phi()-_neutrino_selected.phi())));

      return _Wboson;
    }
    
    /// Neutrino reconstruction from MET via W mass constraint
    Particle returnWboson_byMET(Particle _lepton, double MET_x, double MET_y, double MET, double MET_phi) {
      Vector3 neutrino2(MET_x, MET_y, 0.0);
      
      FourMomentum _neutrino = reconstructNeutrino(neutrino2, _lepton);
      //cout << "Fitting: pz=" << _neutrino.z() << endl;
      FourMomentum _neutrino_withoutFitting = reconstructNeutrino_withoutFitting(neutrino2, _lepton);
      //cout << "Without: pz=" << _neutrino_withoutFitting.z() << endl;
      
      FourMomentum _Wboson_mom = _lepton.momentum() + _neutrino;
      
      _lepton_from_W_boson = _lepton; 
      Particle _Wboson(24, _Wboson_mom);

      // store mt(W)
      mtw = sqrt(2*_lepton.momentum().pT()*GeV*MET*GeV*(1-cos(_lepton.momentum().phi()-MET_phi)));

      return _Wboson;
    }
    
    /* Neutrino code by Pienpen, adjusted to match the rivet vector classes */
    FourMomentum reconstructNeutrino(const Vector3& neutrino, const Particle& lepton) {
      TLorentzVector* met = new TLorentzVector(neutrino.x(), neutrino.y(), 0., sqrt(neutrino.x()*neutrino.x() + neutrino.y()*neutrino.y()));
      TLorentzVector* lep = new TLorentzVector(lepton.momentum().x(), lepton.momentum().y(), lepton.momentum().z(), lepton.momentum().E());
      
      double* pz  = new double[2];
      int nsol = pz_of_W(*lep, met, pz);
      
      double pz_sol = 0.0;

      if (nsol > 0) {
        pz_sol = (std::fabs(pz[0]) < std::fabs(pz[1]) ? pz[0] : pz[1] );
        //cout << pz[0] << "  " << pz[1] << "  -->  " << pz_sol:
      } else
        cout << "ArachneRecoSingleTop::reconstruct: Problem with neutrino reco" << endl;
      
      double met_E = sqrt((*met).X()*(*met).X() + (*met).Y()*(*met).Y() + pz_sol*pz_sol);
      return FourMomentum(met_E, (*met).X(), (*met).Y(), pz_sol);
    }
    
    int pz_of_W(TLorentzVector lep, TLorentzVector* met, double* pz)
    {
      // function for the reconstruction of the pz value of the W
      // using a W mass constrained
      // If the transverse mass of the W is bigger than the W mass itself
      // a fit changing the met is performed to ensure that the transverse
      // W mass is equal to the W mass

      const double MW = 80.4;
      int nsol = 0;
      //double pz[2]; // container for the two pz solutions;

      // filling the global variables neccessary for Minuit
      // lepton
      NeutrinoFit::FullReco_PTe = lep.Pt();
      NeutrinoFit::FullReco_Pxe = lep.Px();
      NeutrinoFit::FullReco_Pye = lep.Py();
      // met
      NeutrinoFit::FullReco_MET_X = met->Px();
      NeutrinoFit::FullReco_MET_Y = met->Py();

      double MisET2 = met->Px()*met->Px() + met->Py()*met->Py();
      //transverse W-mass
      double mWT = sqrt(2*(lep.Pt()*sqrt(MisET2) - lep.Px()*met->Px() - lep.Py()*met->Py()));


      double PxNu_1 = 0;
      double PxNu_2 = 0;
      double PyNu_1 = 0.;
      double PyNu_2 = 0.;
      double PxNu = 0.;
      double PyNu = 0.;


      bool isComplex = false;

      //check, whether transverse W-mass is larger than pole-mass mW or not.
      if(mWT>MW){
        other_checks[0]++;
        isComplex = true;

        // if mWT > mW correct PxNu and PyNu until mWT = MW
        PxNu_1 = metfit(-1,1,MW); //(Printlevel, y-solution, MW)
        PxNu_2 = metfit(-1,2,MW);


        PyNu_1 = ((MW*MW*NeutrinoFit::FullReco_Pye + 2*NeutrinoFit::FullReco_Pxe*NeutrinoFit::FullReco_Pye*PxNu_1)
            -(MW*NeutrinoFit::FullReco_PTe)*(sqrt(MW*MW + 4*NeutrinoFit::FullReco_Pxe*PxNu_1)))
          /(2*NeutrinoFit::FullReco_Pxe*NeutrinoFit::FullReco_Pxe);
        PyNu_2 = ((MW*MW*NeutrinoFit::FullReco_Pye + 2*NeutrinoFit::FullReco_Pxe*NeutrinoFit::FullReco_Pye*PxNu_2)
            +(MW*NeutrinoFit::FullReco_PTe)*(sqrt(MW*MW + 4*NeutrinoFit::FullReco_Pxe*PxNu_2)))
          /(2*NeutrinoFit::FullReco_Pxe*NeutrinoFit::FullReco_Pxe);


        //Calculate delta1 and delta2 from PxNu_1 and PxNu_2:

        double delta1 =  sqrt((PxNu_1 - NeutrinoFit::FullReco_MET_X)*(PxNu_1 - NeutrinoFit::FullReco_MET_X) +(PyNu_1 - NeutrinoFit::FullReco_MET_Y)*(PyNu_1 - NeutrinoFit::FullReco_MET_Y));

        double delta2 =  sqrt((PxNu_2 - NeutrinoFit::FullReco_MET_X)*(PxNu_2 - NeutrinoFit::FullReco_MET_X) +(PyNu_2 - NeutrinoFit::FullReco_MET_Y)*(PyNu_2 - NeutrinoFit::FullReco_MET_Y));


        //PxNu and PyNu(PxNu):
        if(delta1<delta2) {
          PxNu = PxNu_1;
          PyNu = PyNu_1;
        } else {
          PxNu = PxNu_2;
          PyNu = PyNu_2;
        }
      }

      double pz1,pz2;

      // z component of neutrino momentum ...
      if( !isComplex) {
        // ...for two real solutions (mWT < MW)

        double mu = (MW*MW)/2 + met->Px()*lep.Px() + met->Py()*lep.Py();
        double a = (mu*lep.Pz()) / (lep.E()*lep.E() - lep.Pz()*lep.Pz());
        double a2 = TMath::Power(a,2);
        double b = (TMath::Power(lep.E(),2.)*MisET2 - TMath::Power(mu,2.))/
          (TMath::Power(lep.E(),2) - TMath::Power(lep.Pz(),2));


        if(a2-b < 0) {
          //std::cout<<"Complex !!!"<<std::endl; // should not happen anymore!
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
        // ... for complex solutions (mWT > MW): PzNu_OLD and new correction (mWT = MW)

        double mu = (MW*MW)/2 + PxNu*lep.Px() + PyNu*lep.Py();
        double a = (mu*lep.Pz())/(lep.E()*lep.E() - lep.Pz()*lep.Pz());
        double a2 = TMath::Power(a,2);
        double b = (TMath::Power(lep.E(),2.)*(PxNu*PxNu+PyNu*PyNu) - TMath::Power(mu,2.))/
          (TMath::Power(lep.E(),2) - TMath::Power(lep.Pz(),2));

        if (a2-b >= 0) {
          //double root = sqrt(a2-b);
          //std::cout<<"a+sqrt(a2-b) = "<<(a+root)<<" a-sqrt(a2-b) = "<<(a-root)<<std::endl;
        }

        pz1 = a;
        pz2 = a;
        nsol = 1;
      }

      // smaller pz solution is written to entry 0
      if (fabs(pz1) <= fabs(pz2) ) {
        pz[0] = pz1;
        pz[1] = pz2;
      } else {
        pz[0] = pz2;
        pz[1] = pz1;
      }

      // change values of met in case of complex solution
      if (isComplex) {
        met->SetPx(PxNu);
        met->SetPy(PyNu);
      }


      return nsol;
    }

    double metfit(double fitterPrintLevel, int ysol, double mW)
    {

      //double Pxnu = 0.0;
      TMinuit* minu = new TMinuit(5);

      if (ysol == 1) minu->SetFCN(NeutrinoFit::delta1fcn);
      if (ysol == 2) minu->SetFCN(NeutrinoFit::delta2fcn);


      double arglist[20];
      int ierflg = 0;

      // Set print level.
      arglist[0] = fitterPrintLevel;
      minu->mnexcm("SET PRINT", arglist, 1, ierflg);
      // Set strategy. Possible values: 0, 1, 2
      arglist[0] =  2.0;
      minu->mnexcm("SET STRATEGY", arglist, 1, ierflg);
      arglist[0] =  1.0;
      minu->mnexcm("CALL FCN", arglist, 1, ierflg);

      // Calculate limits for the parameter:

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


      // Set parameters:

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


      //Pxnu  = px_fit;

      if (fitterPrintLevel > 0) {
        std::cout<<"*******************************"<<std::endl;
        std::cout << "Fit Results: Px(nu) = " << px_fit
            << " +- " << px_fit_error << std::endl;
        std::cout<<"*******************************"<<std::endl;
      }

      delete minu;
      return px_fit;
    }    
    /* End of neutrino code by Pienpen*/
    
    /* Neutrino code by Ozan, adjusted to match the rivet vector classes */
    FourMomentum reconstructNeutrino_withoutFitting(const Vector3& neutrino, const Particle& lepton) {
      float motherMass = 80.399*GeV; // in GeV
      Vector3 nuP(neutrino.x(), neutrino.y(), 0.0);

      vector<float> pZ = getNeutrinoPzSolutions(neutrino, lepton, motherMass);
      int nSolutions = pZ.size();
      if(nSolutions == 2) {
        nuP.setZ( std::fabs(pZ[0]) < std::fabs(pZ[1]) ? pZ[0] : pZ[1] );
      } else if(nSolutions == 0) {
        float mTSq = 2. * (nuP.mod()*lepton.momentum().pT()
                - nuP.x()*lepton.momentum().x()
                - nuP.y()*lepton.momentum().y());
        // nuP.SetPerp( motherMass*motherMass/mTSq*nuP.Pt() );
        nuP *= motherMass*motherMass/mTSq; // Scale down pt(nu) such that there is exactly 1 solution
        nuP.setZ(nuP.mod()/lepton.momentum().pT()*lepton.momentum().pz()); // p_nuT adjustment approach
      } else if(nSolutions == 1) {
        nuP.setZ(pZ[0]);
      } else {
        std::cout << "(ERROR)\tKinematics::reconstuctNeutrino():\tImpossible number of pZ solutions: " << nSolutions << ". Setting to NAN." << std::endl;
        nuP.setZ(NAN);
      }
      
      //cout << "Neutrino:   " << nuP.mod() << "   " << nuP.x() << "   " << nuP.y() << "   " << nuP.z() << endl;
      //cout << "Electron:   " << lepton.momentum().E() << "   " << lepton.momentum().x() << "   " << lepton.momentum().y() << "   " << lepton.momentum().z() << endl;
      //cout << "Result W:   " << sqrt(2.0 * (nuP.mod() * lepton.momentum().E() - nuP.x() * lepton.momentum().x() - nuP.y() * lepton.momentum().y() - nuP.z() * lepton.momentum().z()) + (lepton.momentum().E() * lepton.momentum().E() - lepton.momentum().x()*lepton.momentum().x() - lepton.momentum().y()*lepton.momentum().y() - lepton.momentum().z()*lepton.momentum().z())) << endl;

      return FourMomentum(nuP.mod(), nuP.x(), nuP.y(), nuP.z());
    }
    
    vector<float> getNeutrinoPzSolutions(const Vector3 & neutrino, const Particle& lepton, float motherMass) {
      vector<float> pz;

      float alpha = 0.5 * motherMass * motherMass + lepton.momentum().x() * neutrino.x() + lepton.momentum().y() * neutrino.y();
      float pT_lep2 = lepton.momentum().perp2();
      float discriminant = lepton.momentum().vector3().mod2() * (alpha * alpha - pT_lep2 * neutrino.mod2());
      if (discriminant < 0.)
        return pz;

      float pz_offset = alpha * lepton.momentum().z() / pT_lep2;

      float squareRoot = sqrt(discriminant);
      if(squareRoot / pT_lep2 < 1.e-6)
        pz.push_back(pz_offset);
      else {
        if(pz_offset > 0) {
          pz.push_back(pz_offset - squareRoot / pT_lep2);
          pz.push_back(pz_offset + squareRoot / pT_lep2);
        } else {
          pz.push_back(pz_offset + squareRoot / pT_lep2);
          pz.push_back(pz_offset - squareRoot / pT_lep2);
        }
      }

      return pz;
    }
    /* End of neutrino code by Ozan */
    
    void fillWboson(Particle _Wboson, bool isTop) {
      if(isTop) {
        FillFourMomentum(_Wboson, _h_W_top, 0, true);
        _h_W_top.at(7)->fill(mtw, _weight);
      } else {
        FillFourMomentum(_Wboson, _h_W_antitop, 0, true);
        _h_W_antitop.at(7)->fill(mtw, _weight);
      }
      //cout << Mass(_Wboson) << endl;
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
    
    void fillTopQuark(Particle _Wboson,  Jets _bjets, bool isTop) {
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
      
      if(isTop) {
        FillFourMomentum(_topquark, _h_t_top, 0, true);
        _h_t_top.at(7)->fill(sqrt(2*_Wboson.momentum().pT()/GeV*_bjet_selected.pT()/GeV*(1-cos(_Wboson.momentum().phi()-_bjet_selected.phi()))), _weight);
      } else {
        FillFourMomentum(_topquark, _h_t_antitop, 0, true);
        _h_t_antitop.at(7)->fill(sqrt(2*_Wboson.momentum().pT()/GeV*_bjet_selected.pT()/GeV*(1-cos(_Wboson.momentum().phi()-_bjet_selected.phi()))), _weight);
      }
      
      // Filling as in the tchan-parton routine
      if(isTop) {
        for (unsigned int i=0; i< _h_t_pt_top.size(); i++) 
          _h_t_pt_top[i]->fill(_topquark.momentum().pT()/GeV, _weight);
        for (unsigned int i=0; i< _h_t_rap_top.size(); i++) 
          _h_t_rap_top[i]->fill(_topquark.momentum().absrap(), _weight);
      } else {
        for (unsigned int i=0; i< _h_t_pt_antitop.size(); i++) 
          _h_t_pt_antitop[i]->fill(_topquark.momentum().pT()/GeV, _weight);
        for (unsigned int i=0; i< _h_t_rap_antitop.size(); i++) 
          _h_t_rap_antitop[i]->fill(_topquark.momentum().absrap(), _weight);
      }
    }

    //=====================================================================
    // fill properties of the leading particle - FastJets
    //=====================================================================
    void fillLeadingParticleProperties(const FastJets& Collection, vector<Histo1DPtr> _histo, int leadership, float ptCut, float etaCut) {
      if (Collection.size() > 0) {
        FourMomentum FourMom(0., 0., 0., 0.);
        FourMom = Collection.jetsByPt(ptCut*GeV)[leadership].momentum();
        if (fabs(FourMom.eta()) < etaCut) {
          _histo.at(0)->fill(FourMom.E()/GeV, _weight); 
          _histo.at(1)->fill(FourMom.Et()/GeV, _weight);  
          _histo.at(2)->fill(FourMom.pT()/GeV, _weight);  
          _histo.at(3)->fill(FourMom.pz()/GeV, _weight);  
          _histo.at(4)->fill(FourMom.eta(), _weight);
          _histo.at(5)->fill(FourMom.phi(), _weight);
        }
      }
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

      std::vector <double> ptstudy1 {0,50,100,150,200,500};
      std::vector <double> ptstudy2 {0,45,75,110,150,200,500};
      std::vector <double> ptstudy3 {0,40,75,110,145,175,210,500};
      std::vector <double> ptstudy4 {0,40,75,110,145,175,210,250,500};
      std::vector <double> rapstudy1 {0,0.3,0.7,1.3,3};
      std::vector <double> rapstudy2 {0,0.28,0.68,1.05,1.6,3};
      std::vector <double> rapstudy3 {0,0.28,0.68,1.05,1.6,2.3,3};
      
      /// Histograms for the tops
      
      _h_weight = bookHisto1D("weight", 20, -10.5, 9.5);
	  /*
      // mc weight
      _h_weight_top = bookHisto1D("weight_top", 20, -10.5, 9.5);

      // stable charged particles pT distribution
      _h_charged_top.push_back(bookHisto1D("charged_N_top", 25, 0, 25));
      _h_charged_top.push_back(bookHisto1D("charged_E_top", 50, 0, 250));
      _h_charged_top.push_back(bookHisto1D("charged_Et_top", 50, 0, 150));
      _h_charged_top.push_back(bookHisto1D("charged_pt_top", 50, 0, 150));
      _h_charged_top.push_back(bookHisto1D("charged_pz_top", 50, 0, 250));
      _h_charged_top.push_back(bookHisto1D("charged_eta_top", 40, -5.0, 5.0));
      _h_charged_top.push_back(bookHisto1D("charged_phi_top", 32, 0.0, twopi));
      _h_charged_top.push_back(bookHisto1D("charged_m_top", 50, 0, 50));

      // electrons
      _h_electrons_top.push_back(bookHisto1D("electrons_N_top", 5., 0., 5.));
      _h_electrons_top.push_back(bookHisto1D("electrons_E_top", 50, 0, 250));
      _h_electrons_top.push_back(bookHisto1D("electrons_Et_top", 50, 0, 150));
      _h_electrons_top.push_back(bookHisto1D("electrons_pt_top", 50, 0, 150));
      _h_electrons_top.push_back(bookHisto1D("electrons_pz_top", 50, 0, 250));
      _h_electrons_top.push_back(bookHisto1D("electrons_eta_top", 40, -5.0, 5.0));
      _h_electrons_top.push_back(bookHisto1D("electrons_phi_top", 32, 0.0, twopi));

      // muons
      _h_muons_top.push_back(bookHisto1D("muons_N_top", 5., 0., 5.));
      _h_muons_top.push_back(bookHisto1D("muons_E_top", 50, 0, 250));
      _h_muons_top.push_back(bookHisto1D("muons_Et_top", 50, 0, 150));
      _h_muons_top.push_back(bookHisto1D("muons_pt_top", 50, 0, 150));
      _h_muons_top.push_back(bookHisto1D("muons_pz_top", 50, 0, 250));
      _h_muons_top.push_back(bookHisto1D("muons_eta_top", 40, -5.0, 5.0));
      _h_muons_top.push_back(bookHisto1D("muons_phi_top", 32, 0.0, twopi));

      // electrons
      _h_veto_electrons_top.push_back(bookHisto1D("veto_electrons_N_top", 5., 0., 5.));
      _h_veto_electrons_top.push_back(bookHisto1D("veto_electrons_E_top", 50, 0, 250));
      _h_veto_electrons_top.push_back(bookHisto1D("veto_electrons_Et_top", 50, 0, 150));
      _h_veto_electrons_top.push_back(bookHisto1D("veto_electrons_pt_top", 50, 0, 150));
      _h_veto_electrons_top.push_back(bookHisto1D("veto_electrons_pz_top", 50, 0, 250));
      _h_veto_electrons_top.push_back(bookHisto1D("veto_electrons_eta_top", 40, -5.0, 5.0));
      _h_veto_electrons_top.push_back(bookHisto1D("veto_electrons_phi_top", 32, 0.0, twopi));

      // muons
      _h_veto_muons_top.push_back(bookHisto1D("veto_muons_N_top", 5., 0., 5.));
      _h_veto_muons_top.push_back(bookHisto1D("veto_muons_E_top", 50, 0, 250));
      _h_veto_muons_top.push_back(bookHisto1D("veto_muons_Et_top", 50, 0, 150));
      _h_veto_muons_top.push_back(bookHisto1D("veto_muons_pt_top", 50, 0, 150));
      _h_veto_muons_top.push_back(bookHisto1D("veto_muons_pz_top", 50, 0, 250));
      _h_veto_muons_top.push_back(bookHisto1D("veto_muons_eta_top", 40, -5.0, 5.0));
      _h_veto_muons_top.push_back(bookHisto1D("veto_muons_phi_top", 32, 0.0, twopi));

      // electrons
      _h_ew_electrons_top.push_back(bookHisto1D("ew_electrons_N_top", 10., 0., 10.));
      _h_ew_electrons_top.push_back(bookHisto1D("ew_electrons_E_top", 50, 0, 250));
      _h_ew_electrons_top.push_back(bookHisto1D("ew_electrons_Et_top", 50, 0, 150));
      _h_ew_electrons_top.push_back(bookHisto1D("ew_electrons_pt_top", 50, 0, 150));
      _h_ew_electrons_top.push_back(bookHisto1D("ew_electrons_pz_top", 50, 0, 250));
      _h_ew_electrons_top.push_back(bookHisto1D("ew_electrons_eta_top", 40, -5.0, 5.0));
      _h_ew_electrons_top.push_back(bookHisto1D("ew_electrons_phi_top", 32, 0.0, twopi));

      // muons
      _h_ew_muons_top.push_back(bookHisto1D("ew_muons_N_top", 10., 0., 10.));
      _h_ew_muons_top.push_back(bookHisto1D("ew_muons_E_top", 50, 0, 250));
      _h_ew_muons_top.push_back(bookHisto1D("ew_muons_Et_top", 50, 0, 150));
      _h_ew_muons_top.push_back(bookHisto1D("ew_muons_pt_top", 50, 0, 150));
      _h_ew_muons_top.push_back(bookHisto1D("ew_muons_pz_top", 50, 0, 250));
      _h_ew_muons_top.push_back(bookHisto1D("ew_muons_eta_top", 40, -5.0, 5.0));
      _h_ew_muons_top.push_back(bookHisto1D("ew_muons_phi_top", 32, 0.0, twopi));

      // neutrinos
      _h_prompt_neutrinos_top.push_back(bookHisto1D("neutrinos_N_prompt_top", 5., 0., 5.));
      _h_prompt_neutrinos_top.push_back(bookHisto1D("neutrinos_E_prompt_top", 50, 0, 300));
      _h_prompt_neutrinos_top.push_back(bookHisto1D("neutrinos_Et_prompt_top", 50, 0, 200));
      _h_prompt_neutrinos_top.push_back(bookHisto1D("neutrinos_pt_prompt_top", 50, 0, 200));
      _h_prompt_neutrinos_top.push_back(bookHisto1D("neutrinos_pz_prompt_top", 50, 0, 250));
      _h_prompt_neutrinos_top.push_back(bookHisto1D("neutrinos_eta_prompt_top", 40, -5.0, 5.0));
      _h_prompt_neutrinos_top.push_back(bookHisto1D("neutrinos_phi_prompt_top", 32, 0.0, twopi));

      // neutrinos
      _h_met_neutrinos_top.push_back(bookHisto1D("neutrinos_N_met_top", 10., 0., 10.));
      _h_met_neutrinos_top.push_back(bookHisto1D("neutrinos_E_met_top", 50, 0, 300));
      _h_met_neutrinos_top.push_back(bookHisto1D("neutrinos_Et_met_top", 50, 0, 200));
      _h_met_neutrinos_top.push_back(bookHisto1D("neutrinos_pt_met_top", 50, 0, 200));
      _h_met_neutrinos_top.push_back(bookHisto1D("neutrinos_pz_met_top", 50, 0, 250));
      _h_met_neutrinos_top.push_back(bookHisto1D("neutrinos_eta_met_top", 40, -5.0, 5.0));
      _h_met_neutrinos_top.push_back(bookHisto1D("neutrinos_phi_met_top", 32, 0.0, twopi));

      // leading electron
      _h_electron_1_top.push_back(bookHisto1D("electron_1_E_top", 50, 0, 250));
      _h_electron_1_top.push_back(bookHisto1D("electron_1_Et_top", 50, 0, 150));
      _h_electron_1_top.push_back(bookHisto1D("electron_1_pt_top", 50, 0, 150));
      _h_electron_1_top.push_back(bookHisto1D("electron_1_pz_top", 50, 0, 250));
      _h_electron_1_top.push_back(bookHisto1D("electron_1_eta_top", 40, -5.0, 5.0));
      _h_electron_1_top.push_back(bookHisto1D("electron_1_phi_top", 32, 0.0, twopi));

      // leading muon
      _h_muon_1_top.push_back(bookHisto1D("muon_1_E_top", 50, 0, 250));
      _h_muon_1_top.push_back(bookHisto1D("muon_1_Et_top", 50, 0, 150));
      _h_muon_1_top.push_back(bookHisto1D("muon_1_pt_top", 50, 0, 150));
      _h_muon_1_top.push_back(bookHisto1D("muon_1_pz_top", 50, 0, 250));
      _h_muon_1_top.push_back(bookHisto1D("muon_1_eta_top", 40, -5.0, 5.0));
      _h_muon_1_top.push_back(bookHisto1D("muon_1_phi_top", 32, 0.0, twopi));

      // leading neutrino
      _h_prompt_neutrino_1_top.push_back(bookHisto1D("neutrino_1_prompt_E_top", 50, 0, 250));
      _h_prompt_neutrino_1_top.push_back(bookHisto1D("neutrino_1_prompt_Et_top", 50, 0, 200));
      _h_prompt_neutrino_1_top.push_back(bookHisto1D("neutrino_1_prompt_pt_top", 50, 0, 200));
      _h_prompt_neutrino_1_top.push_back(bookHisto1D("neutrino_1_prompt_pz_top", 50, 0, 250));
      _h_prompt_neutrino_1_top.push_back(bookHisto1D("neutrino_1_prompt_eta_top", 40, -5.0, 5.0));
      _h_prompt_neutrino_1_top.push_back(bookHisto1D("neutrino_1_prompt_phi_top", 32, 0.0, twopi));

      // leading neutrino
      _h_met_neutrino_1_top.push_back(bookHisto1D("neutrino_1_met_E_top", 50, 0, 250));
      _h_met_neutrino_1_top.push_back(bookHisto1D("neutrino_1_met_Et_top", 50, 0, 200));
      _h_met_neutrino_1_top.push_back(bookHisto1D("neutrino_1_met_pt_top", 50, 0, 200));
      _h_met_neutrino_1_top.push_back(bookHisto1D("neutrino_1_met_pz_top", 50, 0, 250));
      _h_met_neutrino_1_top.push_back(bookHisto1D("neutrino_1_met_eta_top", 40, -5.0, 5.0));
      _h_met_neutrino_1_top.push_back(bookHisto1D("neutrino_1_met_phi_top", 32, 0.0, twopi));

      // MET
      _h_MET_top.push_back(bookHisto1D("MET_top", 50, 20, 200));
      _h_MET_top.push_back(bookHisto1D("MET_eta_top", 40, -5.0, 5.0));
      _h_MET_top.push_back(bookHisto1D("MET_phi_top", 32, 0.0, twopi));

      // "good" jets
      _h_jets_top.push_back(bookHisto1D("jets_N_top", 5, 0, 5));
      _h_jets_top.push_back(bookHisto1D("jets_E_top", 50, 0, 250));
      _h_jets_top.push_back(bookHisto1D("jets_Et_top", 50, 0, 250));
      _h_jets_top.push_back(bookHisto1D("jets_pt_top", 50, 0, 250));
      _h_jets_top.push_back(bookHisto1D("jets_pz_top", 50, 0, 400));
      _h_jets_top.push_back(bookHisto1D("jets_eta_top", 40, -5.0, 5.0));
      _h_jets_top.push_back(bookHisto1D("jets_phi_top", 32, 0.0, twopi));
      _h_jets_top.push_back(bookHisto1D("jets_m_top", 50, 0, 50));

      // B-hadrons
      _h_B_hadrons_top.push_back(bookHisto1D("B_hadrons_N_top", 10, 0, 10));
      _h_B_hadrons_top.push_back(bookHisto1D("B_hadrons_E_top", 50, 0, 300));
      _h_B_hadrons_top.push_back(bookHisto1D("B_hadrons_Et_top", 50, 0, 150));
      _h_B_hadrons_top.push_back(bookHisto1D("B_hadrons_pt_top", 50, 0, 150));
      _h_B_hadrons_top.push_back(bookHisto1D("B_hadrons_pz_top", 50, 0, 250));
      _h_B_hadrons_top.push_back(bookHisto1D("B_hadrons_eta_top", 40, -5.0, 5.0));
      _h_B_hadrons_top.push_back(bookHisto1D("B_hadrons_phi_top", 32, 0.0, twopi));
      _h_B_hadrons_top.push_back(bookHisto1D("B_hadrons_m_top", 50, 0, 25));
      _h_B_hadrons_top.push_back(bookHisto1D("B_hadrons_mt_top", 50, 0, 150));

      // light jets
      _h_ljets_top.push_back(bookHisto1D("ljets_N_top", 5, 0, 5));
      _h_ljets_top.push_back(bookHisto1D("ljets_E_top", 50, 0, 450));
      _h_ljets_top.push_back(bookHisto1D("ljets_Et_top", 50, 0, 250));
      _h_ljets_top.push_back(bookHisto1D("ljets_pt_top", 50, 0, 250));
      _h_ljets_top.push_back(bookHisto1D("ljets_pz_top", 50, 0, 500));
      _h_ljets_top.push_back(bookHisto1D("ljets_eta_top", 40, -5.0, 5.0));
      _h_ljets_top.push_back(bookHisto1D("ljets_phi_top", 32, 0.0, twopi));
      _h_ljets_top.push_back(bookHisto1D("ljets_m_top", 50, 0, 50));

      // b-jets
      _h_bjets_top.push_back(bookHisto1D("bjets_N_top", 5, 0, 5));
      _h_bjets_top.push_back(bookHisto1D("bjets_E_top", 50, 0, 400));
      _h_bjets_top.push_back(bookHisto1D("bjets_Et_top", 50, 0, 200));
      _h_bjets_top.push_back(bookHisto1D("bjets_pt_top", 50, 0, 200));
      _h_bjets_top.push_back(bookHisto1D("bjets_pz_top", 50, 0, 350));
      _h_bjets_top.push_back(bookHisto1D("bjets_eta_top", 40, -5.0, 5.0));
      _h_bjets_top.push_back(bookHisto1D("bjets_phi_top", 32, 0.0, twopi));
      _h_bjets_top.push_back(bookHisto1D("bjets_m_top", 50, 0, 50));

      // W boson
      _h_W_top.push_back(bookHisto1D("W_E_top", 50, 50, 500));
      _h_W_top.push_back(bookHisto1D("W_Et_top", 50, 0, 250));
      _h_W_top.push_back(bookHisto1D("W_pt_top", 50, 0, 200));
      _h_W_top.push_back(bookHisto1D("W_pz_top", 50, 0, 400));
      _h_W_top.push_back(bookHisto1D("W_eta_top", 40, -5.0, 5.0));
      _h_W_top.push_back(bookHisto1D("W_phi_top", 32, 0.0, twopi));
      _h_W_top.push_back(bookHisto1D("W_m_top", 50, 40, 120));
      _h_W_top.push_back(bookHisto1D("W_mt_top", 50, 40, 120));

      // top quark
      _h_t_top.push_back(bookHisto1D("t_E_top", 50, 100, 600));
      _h_t_top.push_back(bookHisto1D("t_Et_top", 50, 0, 250));
      _h_t_top.push_back(bookHisto1D("t_pt_top", 50, 0, 250));
      _h_t_top.push_back(bookHisto1D("t_pz_top", 50, 0, 250));
      _h_t_top.push_back(bookHisto1D("t_eta_top", 40, -5.0, 5.0));
      _h_t_top.push_back(bookHisto1D("t_phi_top", 32, 0.0, twopi));
      _h_t_top.push_back(bookHisto1D("t_m_top", 50, 100, 250));
      _h_t_top.push_back(bookHisto1D("t_mt_top", 50, 0, 250));

      // angular distributions
      _h_cosTheta_S_top = bookHisto1D("CosThetaS_top", 32, -1.0, 1.0);
      _h_cosTheta_N_top = bookHisto1D("CosThetaN_top", 32, -1.0, 1.0);
      _h_cosTheta_T_top = bookHisto1D("CosThetaT_top", 32, -1.0, 1.0);
      _h_cosTheta_X_top = bookHisto1D("CosThetaX_top", 32, -1.0, 1.0);
      _h_cosTheta_top = bookHisto1D("CosThetaW_top", 32, -1.0, 1.0);
      _h_Phi_S_top = bookHisto1D("PhiS_top", 32, 0.0, 2.0);

      // histograms from tchan_parton
      _h_t_pt_top.push_back(bookHisto1D("t_pt_1_top", ptstudy1));
      _h_t_pt_top.push_back(bookHisto1D("t_pt_2_top", ptstudy2));
      _h_t_pt_top.push_back(bookHisto1D("t_pt_3_top", ptstudy3));
      _h_t_pt_top.push_back(bookHisto1D("t_pt_4_top", ptstudy4));
      _h_t_pt_top.push_back(bookHisto1D("t_pt_5_top", 50, 0, 500));
      _h_t_pt_top.push_back(bookHisto1D("t_pt_6_top", 100, 0, 500));
      _h_t_rap_top.push_back(bookHisto1D("t_rap_1_top", rapstudy1));
      _h_t_rap_top.push_back(bookHisto1D("t_rap_2_top", rapstudy2));
      _h_t_rap_top.push_back(bookHisto1D("t_rap_3_top", rapstudy3));
      _h_t_rap_top.push_back(bookHisto1D("t_rap_4_top", 15, 0.0, 3.0));
      _h_t_rap_top.push_back(bookHisto1D("t_rap_5_top", 30, 0.0, 3.0));
      _h_t_rap_top.push_back(bookHisto1D("t_rap_6_top", 100, 0.0, 3.0));
      _h_t_rap_top.push_back(bookHisto1D("t_rap_7_top", 150, 0.0, 5.0));
      _h_ljet_pt_top.push_back(bookHisto1D("ljet_pt_1_top", 50, 0, 500));
      _h_ljet_pt_top.push_back(bookHisto1D("ljet_pt_2_top", 100, 0, 500));
      _h_ljet_rap_top.push_back(bookHisto1D("ljet_rap_1_top", 30, 0, 5.0));
      _h_ljet_rap_top.push_back(bookHisto1D("ljet_rap_2_top", 100, 0, 5.0));
      _h_ljet_rap_top.push_back(bookHisto1D("ljet_rap_3_top", 150, 0, 5.0));

      /// Histograms for the antitops

      // mc weight
      _h_weight_antitop = bookHisto1D("weight_antitop", 20, -10.5, 9.5);

      // stable charged particles pT distribution
      _h_charged_antitop.push_back(bookHisto1D("charged_N_antitop", 25, 0, 25));
      _h_charged_antitop.push_back(bookHisto1D("charged_E_antitop", 50, 0, 250));
      _h_charged_antitop.push_back(bookHisto1D("charged_Et_antitop", 50, 0, 150));
      _h_charged_antitop.push_back(bookHisto1D("charged_pt_antitop", 50, 0, 150));
      _h_charged_antitop.push_back(bookHisto1D("charged_pz_antitop", 50, 0, 250));
      _h_charged_antitop.push_back(bookHisto1D("charged_eta_antitop", 40, -5.0, 5.0));
      _h_charged_antitop.push_back(bookHisto1D("charged_phi_antitop", 32, 0.0, twopi));
      _h_charged_antitop.push_back(bookHisto1D("charged_m_antitop", 50, 0, 50));

      // electrons
      _h_electrons_antitop.push_back(bookHisto1D("electrons_N_antitop", 5., 0., 5.));
      _h_electrons_antitop.push_back(bookHisto1D("electrons_E_antitop", 50, 0, 250));
      _h_electrons_antitop.push_back(bookHisto1D("electrons_Et_antitop", 50, 0, 150));
      _h_electrons_antitop.push_back(bookHisto1D("electrons_pt_antitop", 50, 0, 150));
      _h_electrons_antitop.push_back(bookHisto1D("electrons_pz_antitop", 50, 0, 250));
      _h_electrons_antitop.push_back(bookHisto1D("electrons_eta_antitop", 40, -5.0, 5.0));
      _h_electrons_antitop.push_back(bookHisto1D("electrons_phi_antitop", 32, 0.0, twopi));

      // muons
      _h_muons_antitop.push_back(bookHisto1D("muons_N_antitop", 5., 0., 5.));
      _h_muons_antitop.push_back(bookHisto1D("muons_E_antitop", 50, 0, 250));
      _h_muons_antitop.push_back(bookHisto1D("muons_Et_antitop", 50, 0, 150));
      _h_muons_antitop.push_back(bookHisto1D("muons_pt_antitop", 50, 0, 150));
      _h_muons_antitop.push_back(bookHisto1D("muons_pz_antitop", 50, 0, 250));
      _h_muons_antitop.push_back(bookHisto1D("muons_eta_antitop", 40, -5.0, 5.0));
      _h_muons_antitop.push_back(bookHisto1D("muons_phi_antitop", 32, 0.0, twopi));

      // electrons
      _h_veto_electrons_antitop.push_back(bookHisto1D("veto_electrons_N_antitop", 5., 0., 5.));
      _h_veto_electrons_antitop.push_back(bookHisto1D("veto_electrons_E_antitop", 50, 0, 250));
      _h_veto_electrons_antitop.push_back(bookHisto1D("veto_electrons_Et_antitop", 50, 0, 150));
      _h_veto_electrons_antitop.push_back(bookHisto1D("veto_electrons_pt_antitop", 50, 0, 150));
      _h_veto_electrons_antitop.push_back(bookHisto1D("veto_electrons_pz_antitop", 50, 0, 250));
      _h_veto_electrons_antitop.push_back(bookHisto1D("veto_electrons_eta_antitop", 40, -5.0, 5.0));
      _h_veto_electrons_antitop.push_back(bookHisto1D("veto_electrons_phi_antitop", 32, 0.0, twopi));

      // muons
      _h_veto_muons_antitop.push_back(bookHisto1D("veto_muons_N_antitop", 5., 0., 5.));
      _h_veto_muons_antitop.push_back(bookHisto1D("veto_muons_E_antitop", 50, 0, 250));
      _h_veto_muons_antitop.push_back(bookHisto1D("veto_muons_Et_antitop", 50, 0, 150));
      _h_veto_muons_antitop.push_back(bookHisto1D("veto_muons_pt_antitop", 50, 0, 150));
      _h_veto_muons_antitop.push_back(bookHisto1D("veto_muons_pz_antitop", 50, 0, 250));
      _h_veto_muons_antitop.push_back(bookHisto1D("veto_muons_eta_antitop", 40, -5.0, 5.0));
      _h_veto_muons_antitop.push_back(bookHisto1D("veto_muons_phi_antitop", 32, 0.0, twopi));

      // electrons
      _h_ew_electrons_antitop.push_back(bookHisto1D("ew_electrons_N_antitop", 10., 0., 10.));
      _h_ew_electrons_antitop.push_back(bookHisto1D("ew_electrons_E_antitop", 50, 0, 250));
      _h_ew_electrons_antitop.push_back(bookHisto1D("ew_electrons_Et_antitop", 50, 0, 150));
      _h_ew_electrons_antitop.push_back(bookHisto1D("ew_electrons_pt_antitop", 50, 0, 150));
      _h_ew_electrons_antitop.push_back(bookHisto1D("ew_electrons_pz_antitop", 50, 0, 250));
      _h_ew_electrons_antitop.push_back(bookHisto1D("ew_electrons_eta_antitop", 40, -5.0, 5.0));
      _h_ew_electrons_antitop.push_back(bookHisto1D("ew_electrons_phi_antitop", 32, 0.0, twopi));

      // muons
      _h_ew_muons_antitop.push_back(bookHisto1D("ew_muons_N_antitop", 10., 0., 10.));
      _h_ew_muons_antitop.push_back(bookHisto1D("ew_muons_E_antitop", 50, 0, 250));
      _h_ew_muons_antitop.push_back(bookHisto1D("ew_muons_Et_antitop", 50, 0, 150));
      _h_ew_muons_antitop.push_back(bookHisto1D("ew_muons_pt_antitop", 50, 0, 150));
      _h_ew_muons_antitop.push_back(bookHisto1D("ew_muons_pz_antitop", 50, 0, 250));
      _h_ew_muons_antitop.push_back(bookHisto1D("ew_muons_eta_antitop", 40, -5.0, 5.0));
      _h_ew_muons_antitop.push_back(bookHisto1D("ew_muons_phi_antitop", 32, 0.0, twopi));

      // neutrinos
      _h_prompt_neutrinos_antitop.push_back(bookHisto1D("neutrinos_N_prompt_antitop", 5., 0., 5.));
      _h_prompt_neutrinos_antitop.push_back(bookHisto1D("neutrinos_E_prompt_antitop", 50, 0, 300));
      _h_prompt_neutrinos_antitop.push_back(bookHisto1D("neutrinos_Et_prompt_antitop", 50, 0, 200));
      _h_prompt_neutrinos_antitop.push_back(bookHisto1D("neutrinos_pt_prompt_antitop", 50, 0, 200));
      _h_prompt_neutrinos_antitop.push_back(bookHisto1D("neutrinos_pz_prompt_antitop", 50, 0, 250));
      _h_prompt_neutrinos_antitop.push_back(bookHisto1D("neutrinos_eta_prompt_antitop", 40, -5.0, 5.0));
      _h_prompt_neutrinos_antitop.push_back(bookHisto1D("neutrinos_phi_prompt_antitop", 32, 0.0, twopi));

      // neutrinos
      _h_met_neutrinos_antitop.push_back(bookHisto1D("neutrinos_N_met_antitop", 10., 0., 10.));
      _h_met_neutrinos_antitop.push_back(bookHisto1D("neutrinos_E_met_antitop", 50, 0, 300));
      _h_met_neutrinos_antitop.push_back(bookHisto1D("neutrinos_Et_met_antitop", 50, 0, 200));
      _h_met_neutrinos_antitop.push_back(bookHisto1D("neutrinos_pt_met_antitop", 50, 0, 200));
      _h_met_neutrinos_antitop.push_back(bookHisto1D("neutrinos_pz_met_antitop", 50, 0, 250));
      _h_met_neutrinos_antitop.push_back(bookHisto1D("neutrinos_eta_met_antitop", 40, -5.0, 5.0));
      _h_met_neutrinos_antitop.push_back(bookHisto1D("neutrinos_phi_met_antitop", 32, 0.0, twopi));

      // leading electron
      _h_electron_1_antitop.push_back(bookHisto1D("electron_1_E_antitop", 50, 0, 250));
      _h_electron_1_antitop.push_back(bookHisto1D("electron_1_Et_antitop", 50, 0, 150));
      _h_electron_1_antitop.push_back(bookHisto1D("electron_1_pt_antitop", 50, 0, 150));
      _h_electron_1_antitop.push_back(bookHisto1D("electron_1_pz_antitop", 50, 0, 250));
      _h_electron_1_antitop.push_back(bookHisto1D("electron_1_eta_antitop", 40, -5.0, 5.0));
      _h_electron_1_antitop.push_back(bookHisto1D("electron_1_phi_antitop", 32, 0.0, twopi));

      // leading muon
      _h_muon_1_antitop.push_back(bookHisto1D("muon_1_E_antitop", 50, 0, 250));
      _h_muon_1_antitop.push_back(bookHisto1D("muon_1_Et_antitop", 50, 0, 150));
      _h_muon_1_antitop.push_back(bookHisto1D("muon_1_pt_antitop", 50, 0, 150));
      _h_muon_1_antitop.push_back(bookHisto1D("muon_1_pz_antitop", 50, 0, 250));
      _h_muon_1_antitop.push_back(bookHisto1D("muon_1_eta_antitop", 40, -5.0, 5.0));
      _h_muon_1_antitop.push_back(bookHisto1D("muon_1_phi_antitop", 32, 0.0, twopi));

      // leading neutrino
      _h_prompt_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_prompt_E_antitop", 50, 0, 250));
      _h_prompt_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_prompt_Et_antitop", 50, 0, 200));
      _h_prompt_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_prompt_pt_antitop", 50, 0, 200));
      _h_prompt_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_prompt_pz_antitop", 50, 0, 250));
      _h_prompt_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_prompt_eta_antitop", 40, -5.0, 5.0));
      _h_prompt_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_prompt_phi_antitop", 32, 0.0, twopi));

      // leading neutrino
      _h_met_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_met_E_antitop", 50, 0, 250));
      _h_met_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_met_Et_antitop", 50, 0, 200));
      _h_met_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_met_pt_antitop", 50, 0, 200));
      _h_met_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_met_pz_antitop", 50, 0, 250));
      _h_met_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_met_eta_antitop", 40, -5.0, 5.0));
      _h_met_neutrino_1_antitop.push_back(bookHisto1D("neutrino_1_met_phi_antitop", 32, 0.0, twopi));

      // MET
      _h_MET_antitop.push_back(bookHisto1D("MET_antitop", 50, 20, 200));
      _h_MET_antitop.push_back(bookHisto1D("MET_eta_antitop", 40, -5.0, 5.0));
      _h_MET_antitop.push_back(bookHisto1D("MET_phi_antitop", 32, 0.0, twopi));

      // "good" jets
      _h_jets_antitop.push_back(bookHisto1D("jets_N_antitop", 5, 0, 5));
      _h_jets_antitop.push_back(bookHisto1D("jets_E_antitop", 50, 0, 250));
      _h_jets_antitop.push_back(bookHisto1D("jets_Et_antitop", 50, 0, 250));
      _h_jets_antitop.push_back(bookHisto1D("jets_pt_antitop", 50, 0, 250));
      _h_jets_antitop.push_back(bookHisto1D("jets_pz_antitop", 50, 0, 400));
      _h_jets_antitop.push_back(bookHisto1D("jets_eta_antitop", 40, -5.0, 5.0));
      _h_jets_antitop.push_back(bookHisto1D("jets_phi_antitop", 32, 0.0, twopi));
      _h_jets_antitop.push_back(bookHisto1D("jets_m_antitop", 50, 0, 50));

      // B-hadrons
      _h_B_hadrons_antitop.push_back(bookHisto1D("B_hadrons_N_antitop", 10, 0, 10));
      _h_B_hadrons_antitop.push_back(bookHisto1D("B_hadrons_E_antitop", 50, 0, 300));
      _h_B_hadrons_antitop.push_back(bookHisto1D("B_hadrons_Et_antitop", 50, 0, 150));
      _h_B_hadrons_antitop.push_back(bookHisto1D("B_hadrons_pt_antitop", 50, 0, 150));
      _h_B_hadrons_antitop.push_back(bookHisto1D("B_hadrons_pz_antitop", 50, 0, 250));
      _h_B_hadrons_antitop.push_back(bookHisto1D("B_hadrons_eta_antitop", 40, -5.0, 5.0));
      _h_B_hadrons_antitop.push_back(bookHisto1D("B_hadrons_phi_antitop", 32, 0.0, twopi));
      _h_B_hadrons_antitop.push_back(bookHisto1D("B_hadrons_m_antitop", 50, 0, 25));
      _h_B_hadrons_antitop.push_back(bookHisto1D("B_hadrons_mt_antitop", 50, 0, 150));

      // light jets
      _h_ljets_antitop.push_back(bookHisto1D("ljets_N_antitop", 5, 0, 5));
      _h_ljets_antitop.push_back(bookHisto1D("ljets_E_antitop", 50, 0, 450));
      _h_ljets_antitop.push_back(bookHisto1D("ljets_Et_antitop", 50, 0, 250));
      _h_ljets_antitop.push_back(bookHisto1D("ljets_pt_antitop", 50, 0, 250));
      _h_ljets_antitop.push_back(bookHisto1D("ljets_pz_antitop", 50, 0, 500));
      _h_ljets_antitop.push_back(bookHisto1D("ljets_eta_antitop", 40, -5.0, 5.0));
      _h_ljets_antitop.push_back(bookHisto1D("ljets_phi_antitop", 32, 0.0, twopi));
      _h_ljets_antitop.push_back(bookHisto1D("ljets_m_antitop", 50, 0, 50));

      // b-jets
      _h_bjets_antitop.push_back(bookHisto1D("bjets_N_antitop", 5, 0, 5));
      _h_bjets_antitop.push_back(bookHisto1D("bjets_E_antitop", 50, 0, 400));
      _h_bjets_antitop.push_back(bookHisto1D("bjets_Et_antitop", 50, 0, 200));
      _h_bjets_antitop.push_back(bookHisto1D("bjets_pt_antitop", 50, 0, 200));
      _h_bjets_antitop.push_back(bookHisto1D("bjets_pz_antitop", 50, 0, 350));
      _h_bjets_antitop.push_back(bookHisto1D("bjets_eta_antitop", 40, -5.0, 5.0));
      _h_bjets_antitop.push_back(bookHisto1D("bjets_phi_antitop", 32, 0.0, twopi));
      _h_bjets_antitop.push_back(bookHisto1D("bjets_m_antitop", 50, 0, 50));

      // W boson
      _h_W_antitop.push_back(bookHisto1D("W_E_antitop", 50, 50, 500));
      _h_W_antitop.push_back(bookHisto1D("W_Et_antitop", 50, 0, 250));
      _h_W_antitop.push_back(bookHisto1D("W_pt_antitop", 50, 0, 200));
      _h_W_antitop.push_back(bookHisto1D("W_pz_antitop", 50, 0, 400));
      _h_W_antitop.push_back(bookHisto1D("W_eta_antitop", 40, -5.0, 5.0));
      _h_W_antitop.push_back(bookHisto1D("W_phi_antitop", 32, 0.0, twopi));
      _h_W_antitop.push_back(bookHisto1D("W_m_antitop", 50, 40, 120));
      _h_W_antitop.push_back(bookHisto1D("W_mt_antitop", 50, 40, 120));

      // antitop quark
      _h_t_antitop.push_back(bookHisto1D("t_E_antitop", 50, 100, 600));
      _h_t_antitop.push_back(bookHisto1D("t_Et_antitop", 50, 0, 250));
      _h_t_antitop.push_back(bookHisto1D("t_pt_antitop", 50, 0, 250));
      _h_t_antitop.push_back(bookHisto1D("t_pz_antitop", 50, 0, 250));
      _h_t_antitop.push_back(bookHisto1D("t_eta_antitop", 40, -5.0, 5.0));
      _h_t_antitop.push_back(bookHisto1D("t_phi_antitop", 32, 0.0, twopi));
      _h_t_antitop.push_back(bookHisto1D("t_m_antitop", 50, 100, 250));
      _h_t_antitop.push_back(bookHisto1D("t_mt_antitop", 50, 0, 250));

      // angular distributions
      _h_cosTheta_S_antitop = bookHisto1D("CosThetaS_antitop", 32, -1.0, 1.0);
      _h_cosTheta_N_antitop = bookHisto1D("CosThetaN_antitop", 32, -1.0, 1.0);
      _h_cosTheta_T_antitop = bookHisto1D("CosThetaT_antitop", 32, -1.0, 1.0);
      _h_cosTheta_X_antitop = bookHisto1D("CosThetaX_antitop", 32, -1.0, 1.0);
      _h_cosTheta_antitop = bookHisto1D("CosThetaW_antitop", 32, -1.0, 1.0);
      _h_Phi_S_antitop = bookHisto1D("PhiS_antitop", 32, 0.0, 2.0);

      // histograms from tchan_parton
      _h_t_pt_antitop.push_back(bookHisto1D("t_pt_1_antitop", ptstudy1));
      _h_t_pt_antitop.push_back(bookHisto1D("t_pt_2_antitop", ptstudy2));
      _h_t_pt_antitop.push_back(bookHisto1D("t_pt_3_antitop", ptstudy3));
      _h_t_pt_antitop.push_back(bookHisto1D("t_pt_4_antitop", ptstudy4));
      _h_t_pt_antitop.push_back(bookHisto1D("t_pt_5_antitop", 50, 0, 500));
      _h_t_pt_antitop.push_back(bookHisto1D("t_pt_6_antitop", 100, 0, 500));
      _h_t_rap_antitop.push_back(bookHisto1D("t_rap_1_antitop", rapstudy1));
      _h_t_rap_antitop.push_back(bookHisto1D("t_rap_2_antitop", rapstudy2));
      _h_t_rap_antitop.push_back(bookHisto1D("t_rap_3_antitop", rapstudy3));
      _h_t_rap_antitop.push_back(bookHisto1D("t_rap_4_antitop", 15, 0.0, 3.0));
      _h_t_rap_antitop.push_back(bookHisto1D("t_rap_5_antitop", 30, 0.0, 3.0));
      _h_t_rap_antitop.push_back(bookHisto1D("t_rap_6_antitop", 100, 0.0, 3.0));
      _h_t_rap_antitop.push_back(bookHisto1D("t_rap_7_antitop", 150, 0.0, 5.0));
      _h_ljet_pt_antitop.push_back(bookHisto1D("ljet_pt_1_antitop", 50, 0, 500));
      _h_ljet_pt_antitop.push_back(bookHisto1D("ljet_pt_2_antitop", 100, 0, 500));
      _h_ljet_rap_antitop.push_back(bookHisto1D("ljet_rap_1_antitop", 30, 0, 5.0));
      _h_ljet_rap_antitop.push_back(bookHisto1D("ljet_rap_2_antitop", 100, 0, 5.0));
      _h_ljet_rap_antitop.push_back(bookHisto1D("ljet_rap_3_antitop", 150, 0, 5.0));
*/
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

    // top histograms
    Histo1DPtr _h_weight_top;
    vector<Histo1DPtr> _h_charged_top;
    vector<Histo1DPtr> _h_B_hadrons_top;
    vector<Histo1DPtr> _h_jets_top;
    vector<Histo1DPtr> _h_ljets_top;
    vector<Histo1DPtr> _h_bjets_top;
    vector<Histo1DPtr> _h_electrons_top;
    vector<Histo1DPtr> _h_muons_top;
    vector<Histo1DPtr> _h_taus_top;
    vector<Histo1DPtr> _h_prompt_neutrinos_top;
    vector<Histo1DPtr> _h_met_neutrinos_top;
    vector<Histo1DPtr> _h_veto_electrons_top;
    vector<Histo1DPtr> _h_veto_muons_top;
    vector<Histo1DPtr> _h_ew_electrons_top;
    vector<Histo1DPtr> _h_ew_muons_top;
    vector<Histo1DPtr> _h_electron_1_top;
    vector<Histo1DPtr> _h_muon_1_top;
    vector<Histo1DPtr> _h_prompt_neutrino_1_top;
    vector<Histo1DPtr> _h_met_neutrino_1_top;
    vector<Histo1DPtr> _h_t_top;
    vector<Histo1DPtr> _h_W_top;
    vector<Histo1DPtr> _h_MET_top;
    Histo1DPtr _h_cosTheta_S_top;
    Histo1DPtr _h_cosTheta_N_top;
    Histo1DPtr _h_cosTheta_T_top;
    Histo1DPtr _h_cosTheta_X_top;
    Histo1DPtr _h_cosTheta_top;
    Histo1DPtr _h_Phi_S_top;
    vector<Histo1DPtr> _h_t_pt_top;
    vector<Histo1DPtr> _h_t_rap_top;
    vector<Histo1DPtr> _h_ljet_pt_top;
    vector<Histo1DPtr> _h_ljet_rap_top;
    
    // antitop histograms
    Histo1DPtr _h_weight_antitop;
    vector<Histo1DPtr> _h_charged_antitop;
    vector<Histo1DPtr> _h_B_hadrons_antitop;
    vector<Histo1DPtr> _h_jets_antitop;
    vector<Histo1DPtr> _h_ljets_antitop;
    vector<Histo1DPtr> _h_bjets_antitop;
    vector<Histo1DPtr> _h_electrons_antitop;
    vector<Histo1DPtr> _h_muons_antitop;
    vector<Histo1DPtr> _h_taus_antitop;
    vector<Histo1DPtr> _h_prompt_neutrinos_antitop;
    vector<Histo1DPtr> _h_met_neutrinos_antitop;
    vector<Histo1DPtr> _h_veto_electrons_antitop;
    vector<Histo1DPtr> _h_veto_muons_antitop;
    vector<Histo1DPtr> _h_ew_electrons_antitop;
    vector<Histo1DPtr> _h_ew_muons_antitop;
    vector<Histo1DPtr> _h_electron_1_antitop;
    vector<Histo1DPtr> _h_muon_1_antitop;
    vector<Histo1DPtr> _h_prompt_neutrino_1_antitop;
    vector<Histo1DPtr> _h_met_neutrino_1_antitop;
    vector<Histo1DPtr> _h_t_antitop;
    vector<Histo1DPtr> _h_W_antitop;
    vector<Histo1DPtr> _h_MET_antitop;
    Histo1DPtr _h_cosTheta_S_antitop;
    Histo1DPtr _h_cosTheta_N_antitop;
    Histo1DPtr _h_cosTheta_T_antitop;
    Histo1DPtr _h_cosTheta_X_antitop;
    Histo1DPtr _h_cosTheta_antitop;
    Histo1DPtr _h_Phi_S_antitop;
    vector<Histo1DPtr> _h_t_pt_antitop;
    vector<Histo1DPtr> _h_t_rap_antitop;
    vector<Histo1DPtr> _h_ljet_pt_antitop;
    vector<Histo1DPtr> _h_ljet_rap_antitop;

    //@}
  };

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_tchan_particle> plugin_MC_tchan_particle;

}

