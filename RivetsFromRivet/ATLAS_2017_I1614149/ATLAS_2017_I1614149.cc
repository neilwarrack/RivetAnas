// -*- C++ -*
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "fastjet/tools/Filter.hh"

// substructure includes included in fjcontrib-1.021 (http://fastjet.hepforge.org/contrib/)
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"

namespace Rivet {
  class ATLAS_2017_I1614149 : public Analysis {
    public:

      /// Constructor
      ///@brief: Resolved and boosted ttbar l+jets cross sections at 13 TeV
      DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1614149);

      void init() {
        // Eta ranges
	      Cut eta_full = (Cuts::abseta < 5.0);
        Cut lep_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 25*GeV);

        // All final state particles
        FinalState fs(eta_full);

	      IdentifiedFinalState all_photons(fs);
        all_photons.acceptIdPair(PID::PHOTON);

	      // Get photons to dress leptons
        IdentifiedFinalState ph_id(fs);
        ph_id.acceptIdPair(PID::PHOTON);

        // Projection to find the electrons
        IdentifiedFinalState el_id(fs);
        el_id.acceptIdPair(PID::ELECTRON);

	      PromptFinalState photons(ph_id);
        photons.acceptTauDecays(true);
	      addProjection(photons, "photons");

        PromptFinalState electrons(el_id);
        electrons.acceptTauDecays(true);
        DressedLeptons dressedelectrons(photons, electrons, 0.1, lep_cuts);
        addProjection(dressedelectrons, "elecs");
	      DressedLeptons ewdressedelectrons(all_photons, electrons, 0.1, eta_full);

        // Projection to find the muons
        IdentifiedFinalState mu_id(fs);
        mu_id.acceptIdPair(PID::MUON);

        PromptFinalState muons(mu_id);
        muons.acceptTauDecays(true);
        DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts);
        addProjection(dressedmuons, "muons");
        DressedLeptons ewdressedmuons(all_photons, muons, 0.1, eta_full);

        // Projection to find MET
        declare(MissingMomentum(fs), "MET");

        // remove prompt neutrinos from jet clustering
        IdentifiedFinalState nu_id(fs);
        nu_id.acceptNeutrinos();
        PromptFinalState neutrinos(nu_id);
        neutrinos.acceptTauDecays(true);

        // Jet clustering.
        VetoedFinalState vfs(fs);
        vfs.addVetoOnThisFinalState(ewdressedelectrons);
        vfs.addVetoOnThisFinalState(ewdressedmuons);
        vfs.addVetoOnThisFinalState(neutrinos);
        FastJets jets(vfs, FastJets::ANTIKT, 0.4);
        jets.useInvisibles(true);
        addProjection(jets, "jets");

        // Addition of the large-R jets
        VetoedFinalState vfs1(fs);
        vfs1.addVetoOnThisFinalState(neutrinos);
        FastJets fjets(vfs1, FastJets::ANTIKT, 1.);
        fjets.useInvisibles(JetAlg::NO_INVISIBLES);
        fjets.useMuons(JetAlg::NO_MUONS);
        addProjection(fjets, "fjets");

	      bookHists("top_pt_res", 15);
        bookHists("top_absrap_res", 17);
        bookHists("ttbar_pt_res", 19);
        bookHists("ttbar_absrap_res", 21);
        bookHists("ttbar_m_res", 23);
        bookHists("top_pt_boost", 25);
        bookHists("top_absrap_boost", 27);

      }

      void analyze(const Event& event) {

        const double weight = event.weight();

        // Get the selected objects, using the projections.
        vector<DressedLepton> electrons = apply<DressedLeptons>(event, "elecs").dressedLeptons();
        vector<DressedLepton> muons     = apply<DressedLeptons>(event, "muons").dressedLeptons();
        const Jets& jets  = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
        const PseudoJets& all_fjets  = apply<FastJets>(event, "fjets").pseudoJetsByPt();


        // get MET
        const Vector3 met = apply<MissingMomentum>(event, "MET").vectorMPT();

        Jets bjets, lightjets;
	      for (const Jet& jet : jets) {
	        bool b_tagged = jet.bTags(Cuts::pT > 5*GeV).size();
	        if ( b_tagged && bjets.size() < 2)  bjets +=jet;
	        else lightjets += jet;
	      }

        // Implementing large-R jets definition
        // trim the jets
        PseudoJets trimmed_fatJets;
        float Rfilt = 0.2;
        float pt_fraction_min = 0.05;
        fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, Rfilt), fastjet::SelectorPtFractionMin(pt_fraction_min));
        for (PseudoJet pjet : all_fjets)  trimmed_fatJets += trimmer(pjet);
        trimmed_fatJets = fastjet::sorted_by_pt(trimmed_fatJets);
        PseudoJets trimmed_jets;
        for (unsigned int i = 0; i < trimmed_fatJets.size(); ++i) {
          FourMomentum tj_mom = momentum(trimmed_fatJets[i]);
          if (tj_mom.pt() <= 300*GeV)  continue;
          if (tj_mom.abseta() >= 2.0)  continue;
          trimmed_jets.push_back(trimmed_fatJets[i]);
        }

        bool single_electron = (electrons.size() == 1) && (muons.empty());
        bool single_muon = (muons.size() == 1) && (electrons.empty());

        DressedLepton *lepton = NULL;
        if (single_electron)   lepton = &electrons[0];
        else if (single_muon)  lepton = &muons[0];

        if (!single_electron && !single_muon) vetoEvent;

        bool pass_resolved = true;
	      bool num_b_tagged_jets = (bjets.size() == 2);
	      if (!num_b_tagged_jets)  pass_resolved = false;

	      if (jets.size() < 4) pass_resolved = false;

        bool pass_boosted = true;
        int fatJetIndex = -1;
        bool passTopTag = false;
        bool passDphi = false;
        bool passAddJet = false;
        bool goodLepJet = false;
        bool lepbtag = false;
        bool hadbtag=false;
        vector<int> lepJetIndex;
        vector<int> jet_farFromHadTopJetCandidate;
        if (met.mod() < 20*GeV)  pass_boosted = false;
        if (pass_boosted) {
          double transmass = _mT(lepton->momentum(), met);
          if (transmass + met.mod() < 60*GeV)  pass_boosted = false;
        }
        if (pass_boosted) {
          if (trimmed_jets.size() >= 1) {
            for (unsigned int j = 0; j<trimmed_jets.size(); ++j) {
              if (tau32( trimmed_jets.at(j), 1. ) < 0.75 &&
                  momentum(trimmed_jets.at(j)).mass() > 100*GeV &&
                  momentum(trimmed_jets.at(j)).pt() > 300*GeV &&
                  momentum(trimmed_jets.at(j)).pt() < 1500*GeV &&
                  fabs(momentum(trimmed_jets.at(j)).eta()) < 2.) {
                passTopTag = true;
                fatJetIndex = j;
                break;
              }
            }
          }
        }
        if(!passTopTag && fatJetIndex == -1)  pass_boosted = false;
        if (pass_boosted) {
          double dPhi_fatjet = deltaPhi(lepton->phi(), momentum(trimmed_jets.at(fatJetIndex)).phi());
          double dPhi_fatjet_lep_cut  = 1.0; //2.3
          if (dPhi_fatjet > dPhi_fatjet_lep_cut ) {
            passDphi = true;
          }
        }
        if (!passDphi)   pass_boosted = false;
        if (bjets.empty())   pass_boosted = false;
        if (pass_boosted) {
          for (unsigned int sj = 0; sj < jets.size(); ++sj) {
            double dR = deltaR(jets.at(sj).momentum(), momentum(trimmed_jets.at(fatJetIndex)));
            if(dR > 1.5) {
              passAddJet = true;
              jet_farFromHadTopJetCandidate.push_back(sj);
            }
          }
        }
        if (!passAddJet)  pass_boosted = false;
        if (pass_boosted) {
          for (int ltj : jet_farFromHadTopJetCandidate) {
            double dR_jet_lep = deltaR(jets.at(ltj).momentum(), lepton->momentum());
            double dR_jet_lep_cut = 2.0;//1.5
            if (dR_jet_lep < dR_jet_lep_cut) {
              lepJetIndex.push_back(ltj);
              goodLepJet = true;
            }
          }
        }
        if(!goodLepJet)  pass_boosted = false;
        if (pass_boosted) {
          for (int lepj : lepJetIndex) {
            lepbtag = jets.at(lepj).bTags(Cuts::pT > 5*GeV).size();
            if (lepbtag) break;
          }
        }
        double dR_fatBjet_cut = 1.0;
        if (pass_boosted) {
          for (const Jet& bjet : bjets) {
            hadbtag |= deltaR(momentum(trimmed_jets.at(fatJetIndex)), bjet) < dR_fatBjet_cut;
          }
        }

        if (!(lepbtag || hadbtag))  pass_boosted = false;

        FourMomentum pbjet1; //Momentum of bjet1
        FourMomentum pbjet2; //Momentum of bjet
        int Wj1index = -1, Wj2index = -1;

        if (pass_resolved) {

          if ( deltaR(bjets[0], *lepton) <= deltaR(bjets[1], *lepton) ) {
            pbjet1 = bjets[0].momentum();
            pbjet2 = bjets[1].momentum();
          } else {
            pbjet1 = bjets[1].momentum();
            pbjet2 = bjets[0].momentum();
          }

          double bestWmass = 1000.0*TeV;
          double mWPDG = 80.399*GeV;
          for (unsigned int i = 0; i < (lightjets.size() - 1); ++i) {
            for (unsigned int j = i + 1; j < lightjets.size(); ++j) {
              double wmass = (lightjets[i].momentum() + lightjets[j].momentum()).mass();
              if (fabs(wmass - mWPDG) < fabs(bestWmass - mWPDG)) {
                bestWmass = wmass;
                Wj1index = i;
                Wj2index = j;
              }
            }
          }

          FourMomentum pjet1 = lightjets[Wj1index].momentum();
          FourMomentum pjet2 = lightjets[Wj2index].momentum();

          // compute hadronic W boson
          FourMomentum pWhadron = pjet1 + pjet2;
          double pz = computeneutrinoz(lepton->momentum(), met);
          FourMomentum ppseudoneutrino( sqrt(sqr(met.x()) + sqr(met.y()) + sqr(pz)), met.x(), met.y(), pz);

          //compute leptonic, hadronic, combined pseudo-top
          FourMomentum ppseudotoplepton = lepton->momentum() + ppseudoneutrino + pbjet1;
          FourMomentum ppseudotophadron = pbjet2 + pWhadron;
          FourMomentum pttbar = ppseudotoplepton + ppseudotophadron;

          fillHists("top_pt_res", ppseudotophadron.pt()/GeV, weight);
          fillHists("top_absrap_res", ppseudotophadron.absrap(), weight);
          fillHists("ttbar_pt_res", pttbar.pt()/GeV, weight);
          fillHists("ttbar_absrap_res", pttbar.absrap(), weight);
          fillHists("ttbar_m_res", pttbar.mass()/GeV, weight);
        }
        
        if (pass_boosted) {// Boosted selection
          double hadtop_pt= momentum(trimmed_jets.at(fatJetIndex)).pt() / GeV;
          double hadtop_absrap= momentum(trimmed_jets.at(fatJetIndex)).absrap();
          fillHists("top_pt_boost", hadtop_pt, weight);
          fillHists("top_absrap_boost", hadtop_absrap, weight);
        }
      }

      void finalize() {
        // Normalize to cross-section
        const double sf = (crossSection() / sumOfWeights());
        for (HistoMap::value_type& hist : _h) {
          scale(hist.second, sf);
	        if (hist.first.find("_norm") != string::npos)  normalize(hist.second);
        }
      }

      void bookHists(std::string name, unsigned int index) {
         _h[name] = bookHisto1D(index, 1 ,1);
         _h[name + "_norm"] = bookHisto1D(index + 1, 1 ,1);
      }

      void fillHists(std::string name, double value, double weight) {
         _h[name]->fill(value, weight);
         _h[name + "_norm"]->fill(value, weight);
      }

      double _mT(const FourMomentum &l, const Vector3 &met) const {
        return  sqrt(2.0 * l.pT() * met.mod() * (1 - cos(deltaPhi(l, met))) );
      }

      double tau32(const fastjet::PseudoJet &jet, double jet_rad) const {
        double alpha = 1.0;
        fastjet::contrib::NormalizedCutoffMeasure normalized_measure(alpha, jet_rad, 1000000);
        // WTA definition
        // fastjet::contrib::OnePass_WTA_KT_Axes wta_kt_axes;
        // as in JetSubStructure recommendations 
        fastjet::contrib::KT_Axes kt_axes;

        /// NsubjettinessRatio uses the results from Nsubjettiness to calculate the ratio
        /// tau_N/tau_M, where N and M are specified by the user. The ratio of different tau values
        /// is often used in analyses, so this class is helpful to streamline code.
        fastjet::contrib::NsubjettinessRatio tau32_kt(3, 2, kt_axes, normalized_measure);

        double tau32 = tau32_kt.result(jet);
        return tau32;
      }

      double computeneutrinoz(const FourMomentum& lepton, const Vector3 &met) const {
        //computing z component of neutrino momentum given lepton and met
        double pzneutrino;
        double m_W = 80.399; // in GeV, given in the paper
        double k = (( sqr( m_W ) - sqr( lepton.mass() ) ) / 2 ) + (lepton.px() * met.x() + lepton.py() * met.y());
        double a = sqr ( lepton.E() )- sqr ( lepton.pz() );
        double b = -2*k*lepton.pz();
        double c = sqr( lepton.E() ) * sqr( met.mod() ) - sqr( k );
        double discriminant = sqr(b) - 4 * a * c;
        double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
        if (discriminant < 0)  pzneutrino = - b / (2 * a); //if the discriminant is negative
        else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
          double absquad[2];
          for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
          if (absquad[0] < absquad[1])  pzneutrino = quad[0];
          else                          pzneutrino = quad[1];
        }

        return pzneutrino;
      }

    private:

      /// @name Objects that are used by the event selection decisions
      typedef map<string, Histo1DPtr> HistoMap;
      HistoMap _h;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1614149);
}
