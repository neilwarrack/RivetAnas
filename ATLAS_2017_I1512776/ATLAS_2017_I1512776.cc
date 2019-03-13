// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Total and fiducial single-top and anti-top cross-section measurements at 8 TeV
  class ATLAS_2017_I1512776 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1512776);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// For parton level analysis
      declare(PartonicTops(PartonicTops::E_MU_TAU), "LeptonicTops");

      /// For particle level analysis
      // Leptons
      Cut lep_cuts = Cuts::abseta < 2.5 
	&& (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      PromptFinalState bareLeptons(lep_cuts,true,true); // true,true accepts tau and
                                                        // muon decays respectively
      IdentifiedFinalState photons(Cuts::abseta < 2.6 && Cuts::abspid == PID::PHOTON);
      DressedLeptons dressedLeptons(photons, bareLeptons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      declare(dressedLeptons, "Leptons");

      // Jets
      VetoedFinalState vfs;
      vfs.vetoNeutrinos();
      vfs.addVetoOnThisFinalState(dressedLeptons);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      declare(jets, "Jets");

      
      declare(MissingMomentum(FinalState(Cuts::abseta < 4.9)), "MET");


      // Book histograms and counters

      // Loop for equivalent t and tbar selections
      for (size_t itop = 0; itop < 2; ++itop) {
        // const string xtop = (itop == 0) ? "tq" : "tbarq";
        // _c_sumw_fid[itop] =      bookCounter("sumw_fid_"+xtop);
        // _c_Xsec_fid[itop] =      bookCounter("Xsec_fid_"+xtop);
        //
        _h_AbsPtclDiffXsecTPt[itop]   = bookHisto1D(1,1,1+itop);
        _h_AbsPtclDiffXsecTY[itop]    = bookHisto1D(1,1,3+itop);
        _h_NrmPtclDiffXsecTPt[itop]   = bookHisto1D(2,1,1+itop);
        _h_NrmPtclDiffXsecTY[itop]    = bookHisto1D(2,1,3+itop);
        //
        _h_AbsPtclDiffXsecJPt[itop]   = bookHisto1D(5,1,1+itop);
        _h_AbsPtclDiffXsecJY[itop]    = bookHisto1D(5,1,3+itop);
        _h_NrmPtclDiffXsecJPt[itop]   = bookHisto1D(6,1,1+itop);
        _h_NrmPtclDiffXsecJY[itop]    = bookHisto1D(6,1,3+itop);
        //
        _h_AbsPrtnDiffXsecTPt[itop]	  = bookHisto1D(3,1,1+itop);
	_h_AbsPrtnDiffXsecTY[itop]	  = bookHisto1D(3,1,3+itop);
	_h_NrmPrtnDiffXsecTPt[itop]	  = bookHisto1D(4,1,1+itop);
	_h_NrmPrtnDiffXsecTY[itop]	  = bookHisto1D(4,1,3+itop);
      }

      for(int i = 0; i < 9; i++) {
	cutflow[i] = 0;
	cutflow_ele[i] = 0;
	cutflow_mu[i] = 0;
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // PARTON-LEVEL ANALYSIS
      do { // trick to allow early exit from selection without full return

	//const Particles& allTops = apply<PartonicTops>(event, "AllPartonicTops").particles();
	const Particles& allTops = apply<PartonicTops>(event, "LeptonicTops").particles();
        if (allTops.size() != 1) break;
        const Particle& top = allTops[0];

        size_t itop = (top.charge() > 0) ? 0 : 1;
	   _h_AbsPrtnDiffXsecTPt[itop]->fill( allTops[0].pT(),      event.weight());
	   _h_AbsPrtnDiffXsecTY[itop] ->fill(  allTops[0].absrap(), event.weight());
	   _h_NrmPrtnDiffXsecTPt[itop]->fill( allTops[0].pT(),      event.weight());
	   _h_NrmPrtnDiffXsecTY[itop] ->fill(  allTops[0].absrap(), event.weight());

      } while (false);


      // PARTICLE-LEVEL ANALYSIS
      do { // trick to allow early exit from selection without full return
	cutflow[0]++;
        // Lepton selection
        const Particles& leps = apply<FinalState>(event, "Leptons").particlesByPt();
        MSG_DEBUG("  #leps = " << leps.size());
        if (leps.size() != 1) break;
        MSG_DEBUG("  Passed lepton selection");
        cutflow[1]++;
	const Particle& lep = leps[0];

        // Jet selection
        const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::abseta < 4.5 && Cuts::pT > 30*GeV);
	MSG_DEBUG("jets.size() = " << jets.size());
        // Filter into b-jets and light/non-b jets
        Jets bjets, ljets;
        for (const Jet& j : jets){
	  MSG_DEBUG("__new jet__");
	  if (j.abseta() < 2.5) MSG_DEBUG("jet abseta<2.5");
	  if (j.bTagged()) MSG_DEBUG("jet btagged");
	  if (j.bTagged(Cuts::pT > 5*GeV)) MSG_DEBUG("jet btagged with pT > 5GeV");
          if (j.abseta() < 2.5 && j.bTagged(Cuts::pT > 5*GeV)) bjets += j; else ljets += j;
       	}
	if (bjets.size() != 1) {MSG_DEBUG("N_bjets (w pT>5GeV) = " << bjets.size()); break;}
	cutflow[2]++;
	if (ljets.size() != 1) {MSG_DEBUG("N_ljets = " << ljets.size()); break;}
	cutflow[3]++;	
	//if (bjets.size() != 1 || ljets.size() != 1) break;
        MSG_DEBUG("  Passed jet selection");
        

        //const Jet& jet1 = jets[0];
        const Jet& bjet = bjets[0];
        const Jet& ljet = ljets[0];

        // Lepton-jet isolation cut
        if (any(jets, deltaRLess(lep, 0.4))) break;
        cutflow[4]++;

        // Mlb cut
        const FourMomentum plb = lep.mom() + bjet.mom();
        if (plb.mass() > 160*GeV) break;
        MSG_DEBUG("  Passed Mlb selection");
        MSG_DEBUG("    Passed fiducial selection");
        cutflow[6]++;

        // Pseudo-W and pseudo-top
        const MissingMomentum& mm = apply<MissingMomentum>(event, "MET");
        const FourMomentum pW = lep.mom() + mm.missingMom();
        const FourMomentum pTop = pW + bjet.mom();

        // Fill counters and histograms
        size_t itop = (lep.charge() > 0) ? 0 : 1;
        _h_AbsPtclDiffXsecTPt[itop]->fill(pTop.pT(),     event.weight());
	_h_AbsPtclDiffXsecTY[itop] ->fill(pTop.absrap(), event.weight());
	_h_NrmPtclDiffXsecTPt[itop]->fill(pTop.pT(),	 event.weight());
	_h_NrmPtclDiffXsecTY[itop] ->fill(pTop.absrap(), event.weight());
	_h_AbsPtclDiffXsecJPt[itop]->fill(ljet.pT(),	 event.weight());
	_h_AbsPtclDiffXsecJY[itop] ->fill(ljet.absrap(), event.weight());
	_h_NrmPtclDiffXsecJPt[itop]->fill(ljet.pT(),	 event.weight());
	_h_NrmPtclDiffXsecJY[itop] ->fill(ljet.absrap(), event.weight());
	
      } while (false); //< just one iteration

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // branching ratio of W->lepton
      double lepBR = 0.324;

      scale( _h_AbsPtclDiffXsecTPt, crossSection()/femtobarn/sumOfWeights() );
      scale( _h_AbsPtclDiffXsecTY,  crossSection()/sumOfWeights() );
      normalize(_h_NrmPtclDiffXsecTPt);
      normalize(_h_NrmPtclDiffXsecTY);
      //
      scale( _h_AbsPtclDiffXsecJPt, crossSection()/femtobarn/sumOfWeights());
      scale( _h_AbsPtclDiffXsecJY,  crossSection()/sumOfWeights() );
      normalize(_h_NrmPtclDiffXsecJPt);
      normalize(_h_NrmPtclDiffXsecJY);
      //
      scale(_h_AbsPrtnDiffXsecTPt, crossSection()/femtobarn/sumOfWeights()/lepBR );
      scale(_h_AbsPrtnDiffXsecTY,  crossSection() / sumOfWeights() / lepBR );
      normalize(_h_NrmPrtnDiffXsecTPt);
      normalize(_h_NrmPrtnDiffXsecTY);

      for (int i = 0; i < 9; i++) {
        cout << "Cutflow at step " << i << ": " << cutflow[i]<<" ratio: "<<(float)cutflow[i]/(float)cutflow[0] << endl;
      }


    }

    //@}


    /// @name Histograms
    //@{

    // CounterPtr _c_Xsec_fid_tq, _c_sumw_fid, _c_sumw_fid_tq, _c_sumw_prtn_tq;
    Histo1DPtr _h_AbsPtclDiffXsecTPt[2], _h_AbsPtclDiffXsecTY[2], _h_NrmPtclDiffXsecTPt[2], _h_NrmPtclDiffXsecTY[2];
    Histo1DPtr _h_AbsPtclDiffXsecJPt[2], _h_AbsPtclDiffXsecJY[2], _h_NrmPtclDiffXsecJPt[2], _h_NrmPtclDiffXsecJY[2];
    Histo1DPtr _h_AbsPrtnDiffXsecTPt[2], _h_AbsPrtnDiffXsecTY[2], _h_NrmPrtnDiffXsecTPt[2], _h_NrmPrtnDiffXsecTY[2];

    int cutflow[9];
    int cutflow_ele[9];
    int cutflow_mu[9];
    //@}


  };


  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1512776);

}
