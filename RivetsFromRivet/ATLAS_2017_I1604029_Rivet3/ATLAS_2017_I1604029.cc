// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief ttbar + gamma at 8 TeV
  class ATLAS_2017_I1604029 : public Analysis {
  public:

    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1604029);

    // Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;

      // signal photons
      PromptFinalState prompt_ph(Cuts::abspid == PID::PHOTON && Cuts::pT > 15*GeV && Cuts::abseta < 2.37);
      declare(prompt_ph, "photons");

      // bare leptons
      Cut base_cuts = (Cuts::abseta < 2.7) && (Cuts::pT > 10*GeV);
      IdentifiedFinalState bare_leps(base_cuts);
      bare_leps.acceptIdPair(PID::MUON);
      bare_leps.acceptIdPair(PID::ELECTRON);
      declare(bare_leps, "bare_leptons");

      // dressed leptons
      Cut dressed_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 25*GeV);
      PromptFinalState prompt_mu(base_cuts && Cuts::abspid == PID::MUON);
      PromptFinalState prompt_el(base_cuts && Cuts::abspid == PID::ELECTRON);
      IdentifiedFinalState all_photons(fs, PID::PHOTON);
      DressedLeptons elecs(all_photons, prompt_el, 0.1, dressed_cuts);
      declare(elecs, "elecs");
      DressedLeptons muons(all_photons, prompt_mu, 0.1, dressed_cuts);
      declare(muons, "muons");

      // auxiliary projections for 'single-lepton ttbar filter'
      PromptFinalState prompt_lep(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      declare(prompt_lep, "prompt_leps");
      declare(UnstableParticles(), "ufs");

      // jets
      FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "jets");


      // BOOK HISTOGRAMS
      book(_h["pt"],  2,1,1);
      book(_h["eta"], 3,1,1);

    }


    // Perform the per-event analysis
    void analyze(const Event& event) {

      // analysis extrapolated to 1-lepton-plus-jets channel, where "lepton" cannot be a tau
      // (i.e. contribution from dileptonic ttbar where one of the leptons is outside
      // the detector acceptance has been subtracted as a background)
      if (applyProjection<PromptFinalState>(event, "prompt_leps").particles().size() != 1)  vetoEvent;
      for (const auto& p : apply<UnstableParticles>(event, "ufs").particles()) {
        if (p.fromPromptTau())  vetoEvent;
      }

      // photon selection
      Particles photons = applyProjection<PromptFinalState>(event, "photons").particlesByPt();
      Particles bare_leps  = apply<IdentifiedFinalState>(event, "bare_leptons").particles();
      for (const Particle& lep : bare_leps)
        ifilter_discard(photons, deltaRLess(lep, 0.1));
      if (photons.size() != 1)  vetoEvent;
      const Particle& photon = photons[0];

      // jet selection
      Jets jets = apply<JetAlg>(event, "jets").jetsByPt(Cuts::abseta < 2.5 && Cuts::pT > 25*GeV);

      // lepton selection
      const vector<DressedLepton>& elecs = apply<DressedLeptons>(event, "elecs").dressedLeptons();
      const vector<DressedLepton>& all_muons = apply<DressedLeptons>(event, "muons").dressedLeptons();

      // jet photon/electron overlap removal
      for (const DressedLepton& e : elecs)
        ifilter_discard(jets, deltaRLess(e, 0.2, RAPIDITY));
      for (const Particle& ph : photons)
        ifilter_discard(jets, deltaRLess(ph, 0.1, RAPIDITY));
	    if (jets.size() < 4)  vetoEvent;

      // photon-jet minimum deltaR
      double mindR_phjet = 999.;
      for (Jet jet : jets) {
        const double dR_phjet = deltaR(photon, jet);
        if (dR_phjet < mindR_phjet) mindR_phjet = dR_phjet;
      }
      if (mindR_phjet < 0.5)  vetoEvent;

      // muon jet overlap removal
      vector<DressedLepton> muons;
      for (const DressedLepton& mu : all_muons) {
        bool overlaps = false;
        for (const Jet& jet : jets) {
          if (deltaR(mu, jet) < 0.4) {
            overlaps = true;
            break;
          }
        }
        if (overlaps) continue;
        muons.push_back(mu);
      }

      // one electron XOR one muon
      bool isEl = elecs.size() == 1 && muons.size() == 0;
      bool isMu = muons.size() == 1 && elecs.size() == 0;
      if (!isEl && !isMu)  vetoEvent;

      // photon-lepton deltaR
      double mindR_phlep = deltaR(photon, isEl? elecs[0] : muons[0]);
      if (mindR_phlep < 0.7)  vetoEvent;

      // b-tagging
      Jets bjets;
      for (const Jet& jet : jets) {
        if (jet.bTagged(Cuts::pT > 5*GeV)) bjets +=jet;
      }
      if (bjets.empty())  vetoEvent;

      _h["pt"]->fill(photon.pT()/GeV);
      _h["eta"]->fill(photon.abseta());
    }


    // Normalise histograms etc., after the run
    void finalize() {
      const double normto(crossSection() / femtobarn / sumOfWeights());
      for (auto &hist : _h) {  scale(hist.second, normto);  }
    }


  private:

    map<string, Histo1DPtr> _h;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1604029);
}
