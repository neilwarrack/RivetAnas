// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
//#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  /// @brief Parton-level top-pair normalized differential distributions at 7 TeV.
  class ATLAS_2014_I1304289 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1304289);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(PartonicTops(PartonicTops::E_MU), "LeptonicPartonTops") ;
      declare(PartonicTops(PartonicTops::HADRONIC), "HadronicPartonTops") ;

      // Book histograms
      _h_hadronicTopPt   = bookHisto1D(1,1,1) ;
      _h_ttbarMass       = bookHisto1D(2,1,1) ;
      _h_ttbarPt          = bookHisto1D(3,1,1) ;
      _h_ttbarAbsY        = bookHisto1D(4,1,1) ;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Find partonic tops
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt() ;
      const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt() ;

      // veto if not semileptonic
      const bool semileptonic  = (leptonicpartontops.size() == 1 && hadronicpartontops.size() == 1 ) ;
      if ( !semileptonic) vetoEvent ;

      // top quark observables
      const FourMomentum t1P4 = leptonicpartontops[0] ;
      const FourMomentum t2P4 = hadronicpartontops[0] ;
      const FourMomentum ttbarP4 = add(t1P4, t2P4) ;

      // fill histograms with event weight
      const double weight = event.weight() ;

      _h_hadronicTopPt->fill( t2P4.pT(),        weight);
      _h_ttbarMass->fill(     ttbarP4.mass(),   weight);
      _h_ttbarPt->fill(       ttbarP4.pT(),     weight);
      _h_ttbarAbsY->fill(     ttbarP4.absrap(), weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //      normalize(_h_YYYY); // normalize to unity
      //scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section

      // Unit-normalized
      normalize({_h_hadronicTopPt, _h_ttbarMass, _h_ttbarPt}, 1000) ; //< Factor from x/y unit mismatch in HepData
      normalize(_h_ttbarAbsY) ;


    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr  _h_hadronicTopPt, _h_ttbarMass, _h_ttbarPt, _h_ttbarAbsY;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1304289);


}
