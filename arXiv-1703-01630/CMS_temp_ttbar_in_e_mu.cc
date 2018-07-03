// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
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


  /// @brief Parton-level top-pair cross-sections at 7TeV and 8TeV using dileptonic ($e \mu$) decays
  class CMS_temp_ttbar_in_e_mu : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_temp_ttbar_in_e_mu);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      
      FinalState fs;

      Cut lepCuts = Cuts::abseta < 2.4 &&  Cuts::pT > 20*GeV ;
      Cut muCuts  = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;    

      IdentifiedFinalState electron_fs(PID::ELECTRON, lepCuts);
      declare(electron_fs, "Electrons");
      
      IdentifiedFinalState muon_fs(PID::MUON, lepCuts);
      declare(muon_fs, "Muons");
      
      // Parton Level Top Quarks
      declare(PartonicTops(PartonicTops::ELECTRON), "ElectronPartonTops") ;
      declare(PartonicTops(PartonicTops::MUON), "MuonPartonTops") ;

      // Projection for jets
      declare(FastJets(fs, FastJets::ANTIKT, 0.5), "Jets");

      // Book histograms
      /*
      _c_ttbar = bookCounter("ttbarXsec");
      _c_1btag = bookCounter("1BTag");
      _c_2btag = bookCounter("2BTag");
      _c_ttbarEMuPairsOppChar_partonic = bookCounter("emupairsOppChar_partonic");
      _c_error_e  = bookCounter("error_e");
      _c_error_mu = bookCounter("error_mu");
      */

      
      //DEBUGGING:
      _hepmcout.reset(new HepMC::IO_GenEvent("inspectevents.hepmc"));    
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      const Particles& electrons = apply<IdentifiedFinalState>(event, "Electrons").particlesByPt();
      const Particles& muons =     apply<IdentifiedFinalState>(event, "Muons").particlesByPt();

      const Particles elPartonTops = apply<ParticleFinder>(event, "ElectronPartonTops").particlesByPt();
      const Particles muPartonTops = apply<ParticleFinder>(event, "MuonPartonTops").particlesByPt();
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt();


      Jets selectedJets;
      //      Particles selectedElectrons, selectedMuons;
      bool eTopPos = false, eTopNeg=false, mTopPos=false,
	mTopNeg=false, opCharge=false;

      // make sure there is a muon and an electron
      if ( elPartonTops.empty() || muPartonTops.empty() ) {N1;}
      else {

      // make sure there are opposite charge e and mu.
      for (const Particle& eTop : elPartonTops){
	if (eTop.charge() > 0) eTopPos=true;
	if (eTop.charge() < 0) eTopNeg=true;
      }
      for (const Particle& mTop : muPartonTops){
	if (mTop.charge() > 0) mTopPos=true;
	if (mTop.charge() < 0) mTopNeg=true;
      }
      if (!( (eTopPos && mTopNeg) || (eTopNeg && mTopPos) )) {N2; RET; vetoEvent;}
      }
      

      
      // veto if not an e mu pair
      const bool eMuPair = ( elPartonTops.size() == 1 && muPartonTops.size() == 1 );
      if ( !eMuPair ){ N3;}

      else {
      
	// veto if e and mu are not of opposite chare
	const bool oppositeCharge = ( elPartonTops[0].charge() + muPartonTops[0].charge() == 0 );

	if ( oppositeCharge ) { N9;}
	else { N4; }
      }
    }

 

    /// Normalise histograms etc., after the run
    void finalize() {

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

    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_ttbar, _c_1btag, _c_2btag, _c_ttbarEMuPairsOppChar_partonic, _c_error_e, _c_error_mu;
    //@}

    std::unique_ptr<HepMC::IO_GenEvent> _hepmcout;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_temp_ttbar_in_e_mu);


}
