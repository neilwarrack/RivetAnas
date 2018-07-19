// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
//#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

// ============= Neil's Debug Toolkit ============== //                                                   
#define N1 do { std::cout << "1 "; } while(0)
#define V1 do { std::cout << "1 " << endl; vetoEvent; } while(0)
#define N2 do { std::cout << "2 "; } while(0)
#define V2 do { std::cout << "2 " << endl; vetoEvent; } while(0)
#define N3 do { std::cout << "3 "; } while(0)
#define V3 do { std::cout << "3 " << endl; vetoEvent; } while(0)
#define N4 do { std::cout << "4 "; } while(0)
#define V4 do { std::cout << "4 " << endl; vetoEvent; } while(0)
#define N5 do { std::cout << "5 "; } while(0)
#define V5 do { std::cout << "5 " << endl; vetoEvent; } while(0)
#define N6 do { std::cout << "6 "; } while(0)
#define V6 do { std::cout << "6 " << endl; vetoEvent; } while(0)
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


  /// @brief ttbar + gamma at 8 TeV fiducial cross section in dilepton channel
  class CMS_2017_I1607560 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2017_I1607560);


    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;


      // Project out particle level objects
      
      // Leptons
      Cut el_cuts = (Cuts::abseta < 2.5) && ((Cuts::abseta < 1.44) || (Cuts::abseta > 1.56))
	&& (Cuts::pT > 35*GeV) && (Cuts::abspid == PID::ELECTRON);
      Cut mu_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 26*GeV) && (Cuts::abspid == PID::MUON);

      PromptFinalState prompt_el(el_cuts); //< first bool accepts tau decays
      PromptFinalState prompt_mu(mu_cuts);

      declare(prompt_el, "Electrons");
      declare(prompt_mu, "Muons");

      // Photons
      PromptFinalState promptPhotons(Cuts::abspid == PID::PHOTON && Cuts::pT > 25*GeV && Cuts::abseta < 1.44);
      declare(promptPhotons, "Photons");
      
      // jets
      FastJets jets(fs, FastJets::ANTIKT, 0.5);
      declare(jets, "Jets");

      // Find neutrinos for pT_miss cut
      //IdentifiedFinalState neutrinos;
      //neutrinos.acceptNeutrinos();
      //declare(neutrinos, "Neutrinos");

      declare(MissingMomentum(fs), "MissingP");
      
      // Book histograms
      _c_fid_XSec = bookCounter("fiducial_XSec");
      //_h_XXXX = bookProfile1D(1, 1, 1);
      //_h_YYYY = bookHisto1D(2, 1, 1);
      //_h_ZZZZ = bookCounter(3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      Cut jetCuts = (Cuts::abseta < 2.4) && (Cuts::pT > 30*GeV);
      
      // Apply projections
      const Particles electrons = apply<PromptFinalState>(event, "Electrons").particlesByPt();
      const Particles muons =     apply<PromptFinalState>(event, "Muons").particlesByPt();
      const Particles photons =   apply<PromptFinalState>(event, "Photons").particlesByPt();
      const MissingMomentum missP = applyProjection<MissingMomentum>(event, "MissingP");
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(jetCuts);


      // Veto for jets
      if (jets.size() < 3) V1; //vetoEvent;
      
      // Veto if not exactly one particle-level lepton is present 
      if (electrons.size() + muons.size() != 1 ) V2;// vetoEvent;

      // Missing momentum veto
      /*
      const FourMomentum missingFourMomentum = -missP.visibleMomentum();
      if ( missingFourMomentum.pT() < 20*GeV ) V3;// vetoEvent;   
      cout << "MisingPt=" <<  missingFourMomentum.pT() << endl;   
      */
      if (missP.missingPt() < 20*GeV) V3;// vetoEvent;   
      cout << "MisingPt=" << missP.missingPt() << endl;   

      // Prompt photon veto
      if (photons.size() < 1) V4;
      
      // Fill counter
      _c_fid_XSec->fill(weight);
      N9; RET;
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_c_fid_XSec, crossSection()/femtobarn/sumOfWeights()); // norm to cross section

    }

    


    /// @name Histograms
    //@{
    //Profile1DPtr _h_XXXX;
    //Histo1DPtr _h_YYYY;
    CounterPtr _c_fid_XSec;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1607560);


}
