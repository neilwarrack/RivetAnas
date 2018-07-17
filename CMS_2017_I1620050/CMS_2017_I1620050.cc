// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Tools/JetUtils.hh"
//testing:
#include "Rivet/Projections/PartonicTops.hh"

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

  ///@brief: ttbar differential cross-sections in dilepton channel at 13 TeV
  // From paper: "events with leptons associated with tau lepton decays are treated as background"
  class CMS_2017_I1620050 : public Analysis {
  public:


    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2017_I1620050);

    // Book histograms and initialise projections before the run
    void init() {


      const FinalState fs;

      Cut lepCuts = ( Cuts::abseta < 2.4  ) & ( Cuts::pT > 20*GeV );
      IdentifiedFinalState lepfs(fs, {PID::ELECTRON, -PID::ELECTRON, PID::MUON, -PID::MUON});
      PromptFinalState promptLeptons(lepfs); // does not accepts tau decays
      IdentifiedFinalState photons(fs, {PID::PHOTON, -PID::PHOTON});
      DressedLeptons leptons(photons, promptLeptons, 0.1, lepCuts);
      declare(leptons, "Leptons");
      
      // Find neutrinos for jet clustering
      IdentifiedFinalState neutrinos;
      neutrinos.acceptNeutrinos();
      PromptFinalState promptNeutrinos(neutrinos);
      promptNeutrinos.acceptMuonDecays();
      declare(neutrinos, "Neutrinos");
      declare(promptNeutrinos, "PromptNeutrinos");
      
      // Project jets
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(promptNeutrinos);
      vfs.addVetoOnThisFinalState(leptons);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      declareProjection(jets, "Jets");


      // projections for reconstructed W bosons
      WFinder w_el(lepfs, lepCuts, PID::ELECTRON, 35*GeV, 8000*GeV, 10*GeV, 0.1);
      declare(w_el, "W_Electron");

      WFinder w_mu(lepfs,  lepCuts, PID::MUON, 35*GeV, 100*GeV, 10*GeV, 0.1);
      declare(w_mu, "W_Muon");

      // auxiliary projection for testing
      declare(PartonicTops(PartonicTops::ALL),"PartonicTops");

      // BOOK HISTOGRAMS
      _h["pt_lep"] = bookHisto1D(1,1,1);
      _h["pt_jet"] = bookHisto1D(1,1,2);

    }

    // Perform the per-event analysis
    void analyze(const Event& event) {
      
      const double weight = event.weight();

      // Select b-jets
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 2.4);
      Jets bJets;
      for (Jet j : jets){ if (j.bTagged()) bJets.push_back(j); }
      if (bJets.size() < 2) {N1;RET;vetoEvent;}//debug
      //if (bJets.size() < 2) vetoEvent;
      if (bJets.size() > 2) {cout<<"more than 2 b jets! how do we fill bjet pt now, hu?"<< endl;}//debug


      // Select leptons
      const vector<DressedLepton>& leptons = apply<DressedLeptons>(event, "Leptons").dressedLeptons();
      const Particles neutrinos = applyProjection<IdentifiedFinalState>(event, "Neutrinos").particlesByPt();
      const Particles promptNeutrinos = applyProjection<PromptFinalState>(event, "PromptNeutrinos").particlesByPt();

      
      // for testing! We actually need to reconstruct Ws from leptons and nutrinos
      /*
	const Particles w_el = applyProjection<WFinder>(event, "W_Electron").particlesByPt();
	const Particles w_mu = applyProjection<WFinder>(event, "W_Muon").particlesByPt();
	const Particles partonicTops = applyProjection<PartonicTops>(event, "PartonicTops").particlesByPt();
	if (partonicTops.size() != 2) vetoEvent;
	Particles wBosons;
	for (Particle w : w_mu) wBosons.push_back(w);
	for (Particle w : w_el) wBosons.push_back(w);
	cout << "W: " << wBosons.size() << endl;
      */


      
      // veto non dileptonic events
      if (leptons.size() != 2 ) {N2;RET;vetoEvent;}//debug
      //if (leptons.size() != 2 ) vetoEvent;
      cout << "No. of neutrinos(prompt):" << promptNeutrinos.size() << endl;
      if (promptNeutrinos.size() < 2) {N3;RET;vetoEvent;}//debug
      //if (promptNeutrinos.size() < 2) vetoEvent;

      
      // Construct pseudo-W-bosons

      cout << "COMPUTE Ws"<<endl;
      /*for (Particle n : neutrinos){ cout << "n" ;}
	cout << endl;
	for (Particle n : promptNeutrinos){ cout << "p" ;}
	cout << endl;
      */
      FourMomentum p4w1, p4w2, p4lep1=leptons[0], p4lep2=leptons[1], p4n1, p4n2, selectedW1, selectedW2;
      int ctr_n1=0, ctr_n2=0;
      double minMassSumW=8000*GeV;
      
      for (Particle n1 : promptNeutrinos){ // should I use prompt neutrinos ?? the paper just says 'neutrinos' 
	cout << "n1:";
	ctr_n1++;
	cout << ctr_n1 << endl;
	p4n1 = n1;
	p4w1 = add(p4lep1, p4n1);
	for (Particle n2 : promptNeutrinos) {
	  cout << "n2:";
	  ctr_n2++;
	  cout << ctr_n2 << endl;

	  if (ctr_n1 == ctr_n2){cout << "skipping like-neutrinos" << endl; continue; }
	  p4n2 = n2;
	  p4w2 = add(p4lep2, p4n2);

	  if ( abs(80*GeV - p4w1.mass()*GeV) + abs(80*GeV - p4w2.mass()*GeV) < minMassSumW){
	    minMassSumW = abs(80*GeV - p4w1.mass()*GeV) + abs(80*GeV - p4w2.mass()*GeV); // replace minMass
	    selectedW1 = p4w1;
	    selectedW2 = p4w2;
	    cout << "minMass=" << minMassSumW << endl;
	  } else {cout << "min remains unchaged"<<endl;}
	  
	}
	cout << "end of n2 loop" << endl;
	ctr_n2=0; // reset counter
      }
      cout << "end of n1 loop" << endl;
      ctr_n1=0; // reset counter
      cout << "minMass=" << minMassSumW << ". Done!" << endl;


      
      // Construct pseudo-Tops

      cout << "COMPUTE TOPs"<<endl;
      FourMomentum p4t1, p4t2, p4b1, p4b2, selectedTop1, selectedTop2;
      int ctr_b1=0, ctr_b2=0;
      double minMassSumTop=8000*GeV;
      Jet selectedBJet1, selectedBJet2;
      
      for (Jet b1 : bJets){ // should I use prompt neutrinos ?? the paper just says 'neutrinos' 
	cout << "b1:";
	ctr_b1++;
	cout << ctr_b1 << endl;
	p4b1 = b1;
	p4t1 = add(selectedW1, p4b1);
	for (Jet b2 : bJets) {
	  cout << "b2:";
	  ctr_b2++;
	  cout << ctr_b2 << endl;
	  
	  if (ctr_b1 == ctr_b2){cout << "skipping like-jets" << endl; continue; }
	  
	  p4b2 = b2;
	  p4t2 = add(selectedW2, p4b2);
	  
	  if ( abs(173*GeV - p4t1.mass()*GeV) + abs(173*GeV - p4t2.mass()*GeV) < minMassSumTop){
	    minMassSumTop = abs(173*GeV - p4t1.mass()*GeV) + abs(173*GeV - p4t2.mass()*GeV); // replace minMass
	    selectedTop1 = p4t1;
	    selectedTop2 = p4t2;
	    selectedBJet1 = b1;
	    selectedBJet2 = b2;
	    cout << "minMassSumTop=" << minMassSumTop << endl;
	  } else {cout << "min remains unchaged"<<endl;}
	  
	}
	cout << "end of b2 loop" << endl;
	ctr_b2=0; // reset counter
      }
      cout << "end of b1 loop" << endl;
      ctr_b1=0; // reset counter
      cout << "minMassSumTop=" << minMassSumTop << ". Done!" << endl;


      FourMomentum ttbar = add( selectedTop1, selectedTop2);
      double dPhi = deltaPhi(   selectedTop1, selectedTop2);
      
      // Fill single top histos

      _h["pt_lep"]->fill(  leptons[0].pT()/GeV,     weight);
      _h["pt_lep"]->fill(  leptons[1].pT()/GeV,     weight);

      _h["pt_jet"]->fill(  selectedBJet1.pT()/GeV, weight);
      _h["pt_jet"]->fill(  selectedBJet2.pT()/GeV, weight);
      /*
      _h["pt_top"]->fill(  selectedTop1.pT()/GeV,  weight);
      _h["pt_top"]->fill(  selectedTop2.pT()/GeV,  weight);
      _h["rap_top"]->fill( selectedTop1.rap(),     weight);
      _h["rap_top"]->fill( selectedTop2.rap(),     weight);


      // Fill ttbar histos

      _h["pt_ttbar"]->fill(   ttbar.pT()/GeV, weight);
      _h["rap_ttbar"]->fill(  ttbar.rap(),    weight);
      _h["mass_ttbar"]->fill( ttbar.mass(),   weight);
      _h["dPhi_ttbar"]->fill( dPhi,           weight);
      */

      

    }

    // Normalise histograms etc., after the run
    void finalize() {
      const double normto( 1.0 / sumOfWeights() );
      for (auto &hist : _h) {  scale(hist.second, normto);  }

      //      for (auto &hist : _h) {  normalize(hist.second);  }
    }

  private:

    map<string, Histo1DPtr> _h;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1620050);
}
