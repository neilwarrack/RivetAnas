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

#include "MendelMin.h"

using namespace std;

namespace Rivet {

  // Forward declaration of analysis specific functions for neutrino reconstruction
  double target1(const MendelMin::Params&, const MendelMin::Params& );
  double target2(const MendelMin::Params&, const MendelMin::Params& );
  double nuy1(const double, const double, const double , const double );
  double nuy2(const double, const double, const double , const double );
  double rand01();
  Vector3 reconstructNeutrino(Vector3&, FourMomentum&);
  Vector3 getNuXY(Vector3&, FourMomentum&, float, bool* is_complex );

  
  /// Total and fiducial single-top and anti-top cross-section measurements at 8 TeV
  class ATLAS_2017_I1512776 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1512776);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // For parton level analysis
      declare(PartonicTops(PartonicTops::DecayMode::ALL), "AllPartonicTops");
      declare(PartonicTops(PartonicTops::DecayMode::E_MU), "EmuLeptonicTops");
      declare(PartonicTops(PartonicTops::DecayMode::E_MU_TAU),"LeptonicTops");


      // Leptons, jets and MET
      //declare(DressedLeptons(PromptFinalState(Cuts::abseta < 2.6), 0.1, Cuts::abseta < 2.5 && Cuts::pT > 25*GeV), "Leptons");
      //DressedLeptons dressedLeps(PromptFinalState(Cuts::abseta < 2.6), 0.1, Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);

      ////////////////////////////////////////////////////////////////////////////////////
      //// What follows is a direct copy (with lines missing) from single top groups rivet
      //
      // Cut lep_cuts = (Cuts::abseta < 2.5) && (Cuts::pT >= 20*GeV);      
      // // All final state particles
      // FinalState fs(Cuts::abseta < 5);
      // // Get photons to dress leptons
      // IdentifiedFinalState photons(Cuts::abseta < 5);
      // photons.acceptIdPair(PID::PHOTON);
      // // Projection to find the electrons
      // IdentifiedFinalState el_id(fs);
      // el_id.acceptIdPair(PID::ELECTRON);
      // PromptFinalState electrons(el_id);
      // electrons.acceptTauDecays(true);
      // DressedLeptons dressedelectrons(photons, electrons, 0.1,lep_cuts,true,true );
      // // Projection to find the muons
      // IdentifiedFinalState mu_id(fs);
      // mu_id.acceptIdPair(PID::MUON);
      // PromptFinalState muons(mu_id);
      // muons.acceptTauDecays(true);
      // DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts, true, true);
      // IdentifiedFinalState nu;
      // nu.acceptNeutrinos();
      // PromptFinalState neutrinos(nu);
      // neutrinos.acceptTauDecays(true);
      // VetoedFinalState vfs;
      // vfs.addVetoOnThisFinalState(neutrinos);
      // vfs.addVetoOnThisFinalState(dressedmuons);
      // vfs.addVetoOnThisFinalState(dressedelectrons);
      // FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      // jets.useInvisibles();
      // declare(jets, "Jets");
      //////////////END///////////////////////////////////////////////////////////////////

      // MY REWORKING OF THE ABOVE CODE///////////////////////////////////////////////////
     

      // Leptons
      /*
      Cut lep_cuts = Cuts::abseta < 2.5 
	&& (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      PromptFinalState bareLeptons(lep_cuts,true); // 'true' accepts tau decays
      IdentifiedFinalState photons(Cuts::abseta < 2.6 && Cuts::abspid == PID::PHOTON);
      DressedLeptons dressedleptons(photons, bareLeptons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 25*GeV);
      declare(dressedleptons, "Leptons");
      */
      //DressedLeptons vetoDressedLeptons(photons, bareLeptons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 10*GeV);
      
      FinalState fs;
      Cut fs_w = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;
      // new attempt at dressed leptons
      // Lepton cuts
      //Cut FS_Wlept = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
      Cut FS_Wlept = Cuts::pT > 1*GeV;
      // Electrons and muons in Fiducial PS
      //PromptFinalState leptons(FinalState(fs_w && (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON)));
      PromptFinalState leptons(FinalState(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON));
      leptons.acceptTauDecays(true);
      // Get photons to dress leptons
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      DressedLeptons dressedleptons(photons, leptons, 0.1, FS_Wlept, true);
      declare(dressedleptons, "Leptons");


      
      // Jets
      VetoedFinalState vfs;
      vfs.vetoNeutrinos();
      vfs.addVetoOnThisFinalState(dressedleptons);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "Jets");
      //////// END ///////////////////////////////////////////////////////////////////////

      // for checking neutrino z
      IdentifiedFinalState nu;
      nu.acceptNeutrinos();
      PromptFinalState neutrinos(nu);
      //neutrinos.acceptTauDecays(true);
      
      declare(neutrinos, "Neutrinos");

      //vfs.addVetoOnThisFinalState(dressedLeps);
      
      //      declare(MissingMomentum(FinalState(Cuts::abseta < 4.5)), "MET");
      declare(MissingMomentum(FinalState()), "MET");
      //declare(MissingMomentum(lep_cuts), "MET");


      // Book histograms and counters
      // _c_sumw_fid =         bookCounter("sumw_fid_tot");
      // Loop for equivalent t and tbar selections
      for (size_t itop = 0; itop < 2; itop++) {
        // const string xtop = (itop == 0) ? "tq" : "tbarq";
        // _c_sumw_fid[itop] =      bookCounter("sumw_fid_"+xtop);
        // _c_Xsec_fid[itop] =      bookCounter("Xsec_fid_"+xtop);
        //
        //_h_AbsPtclDiffXsecTPt[itop]   = bookHisto1D(1,1,1+itop);
        //_h_AbsPtclDiffXsecTY[itop]    = bookHisto1D(1,1,3+itop);
	
	book(_h_NrmPtclDiffXsecTPt[itop], 2, 1, 1+itop);
	     //[itop]   = bookHisto1D(2,1,1+itop);
        book(_h_NrmPtclDiffXsecTY[itop],2,1,3+itop);
        //
        //_h_AbsPtclDiffXsecJPt[itop]   = bookHisto1D(5,1,1+itop);
        //_h_AbsPtclDiffXsecJY[itop]    = bookHisto1D(5,1,3+itop);
        book(_h_NrmPtclDiffXsecJPt[itop],6,1,1+itop);
        book(_h_NrmPtclDiffXsecJY[itop], 6,1,3+itop);
        //
        //_h_AbsPrtnDiffXsecTPt[itop]	  = bookHisto1D(3,1,1+itop);
	//_h_AbsPrtnDiffXsecTY[itop]	  = bookHisto1D(3,1,3+itop);
	book(_h_NrmPrtnDiffXsecTPt[itop],4,1,1+itop);
	book(_h_NrmPrtnDiffXsecTY[itop],4,1,3+itop);
      }
      
      for (int i = 0; i<9; i++) cutflow[i] = 0;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      _c_eventCtr += 1;
      //if (_c_eventCtr < 2434){vetoEvent;}


      // PARTON-LEVEL ANALYSIS
      do { // trick to allow early exit from selection without full return

	//const Particles& allTops = apply<PartonicTops>(event, "AllPartonicTops").particles();
	const Particles& allTops = apply<PartonicTops>(event, "EmuLeptonicTops").particles();
        if (allTops.size() != 1) break;
        const Particle& top = allTops[0];

        size_t itop = (top.charge() > 0) ? 0 : 1;
	if ((top.charge() > 0) && (itop == 1)) cout << "oops!" << endl ;
	//_h_AbsPrtnDiffXsecTPt[itop]->fill( top.pT(),     event.weight());
	//_h_AbsPrtnDiffXsecTY[itop] ->fill( top.absrap(), event.weight());
	_h_NrmPrtnDiffXsecTPt[itop]->fill( top.pT()/GeV);
	_h_NrmPrtnDiffXsecTY[itop] ->fill( top.absrap());

      } while (false);


      // PARTICLE-LEVEL ANALYSIS
      do { // trick to allow early exit from selection without full return
	cutflow[0]++;
        // Lepton selection
        const Particles& leps = apply<FinalState>(event, "Leptons").particlesByPt();


        MSG_DEBUG("  #leps = " << leps.size());
        if (leps.size() < 1) break;
        MSG_DEBUG("  Passed >1 lepton selection");
        cutflow[1]++;
	const Particle& lep = leps[0];

	//Cut FS_Wlept = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
	if (lep.abseta() > 2.5) break;
	cutflow[2]++;
	if (lep.pT() < 25*GeV) break;
	cutflow[3]++;

        // Jet selection
        const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::abseta < 4.5 && Cuts::pT > 30*GeV);
	MSG_DEBUG("jets.size() = " << jets.size());
        // Filter into b-jets and light/non-b jets
        Jets bjets, ljets;
        for (const Jet& j : jets){
	  if (j.bTagged(Cuts::pT > 5*GeV)) MSG_DEBUG("jet btagged with pT > 5GeV");
          if (j.abseta() < 2.5 && j.bTagged(Cuts::pT > 5*GeV)) bjets += j; else ljets += j;
       	}
	if (bjets.size() != 1) {MSG_DEBUG("N_bjets (w pT>5GeV) = " << bjets.size()); break;}
	cutflow[4]++;
	if (ljets.size() != 1) {MSG_DEBUG("N_ljets = " << ljets.size()); break;}
	cutflow[5]++;	
	//if (bjets.size() != 1 || ljets.size() != 1) break;
        MSG_DEBUG("  Passed jet selection");
        

        //const Jet& jet1 = jets[0];
        const Jet& bjet = bjets[0];
        const Jet& ljet = ljets[0];

        // Lepton-jet isolation cut
        if (any(jets, deltaRLess(lep, 0.4))) break;
        MSG_DEBUG("  Passed jet isolation cuts");
	cutflow[6]++;

        // Mlb cut
        const FourMomentum plb = lep.mom() + bjet.mom();
        if (plb.mass() > 160*GeV) break;
        MSG_DEBUG("  Passed Mlb selection");
        cutflow[7]++;

	// transverse mass cut
	const MissingMomentum& mm = apply<MissingMomentum>(event, "MET");
	
	double delta_phi_lep_met = lep.phi()-mm.missingMom().phi();
	//todo check this line:	
	double tmasslepmet = sqrt(2*lep.pT()*mm.met()*(1-cos(delta_phi_lep_met)));
	if (tmasslepmet < 50*GeV) break;
	MSG_DEBUG("  Passed transverse mass selection");
	cutflow[8]++;

	// lep/jet back-to-back cut
	double delta_phi_jet1_lep = jets[0].phi() - lep.phi();
	double back_to_back=40*GeV*(1-(3.14-abs(delta_phi_jet1_lep))/(3.14-1));
	if (back_to_back < 25*GeV){
	  if (lep.pT() < 25*GeV) break;
	} else {
	  if (lep.pT() < back_to_back) break;
	}
        MSG_DEBUG("  Passed lepton-jet back to back cut selection");
        MSG_DEBUG("  Passed fiducial selection");
	cutflow[9]++;


	// ---------------- top quark reconstruction --------------- //
        MSG_DEBUG("  #Passed fiducial selection, reconstructing neutrino, then W, then top:");

	
	FourMomentum lep_4p = lep.mom();
	Vector3 mmpt = mm.vectorMissingPt();
	
	// reconstruct neutrino
	Vector3 recoNeutrino = reconstructNeutrino(mmpt, lep_4p);
	
	// You have reconstructed a neutrino!
	// Print some stuff about the real (prompt) neutrinos (by pT) in event
	const Particles& neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
	for (const Particle &neutrino : neutrinos){
	  MSG_DEBUG("neutrino     : (" << neutrino.mom().x()<< ", " << neutrino.mom().y() << ", " << neutrino.mom().z() << ")");
	}
	// Print some stuff about the reconstructed neutrino in event
	MSG_DEBUG("reco-nu      : (" << recoNeutrino.x() << ", " << recoNeutrino.y() << ", " << recoNeutrino.z() << ")");


	// record some info to calculate the average differenence of reco-nu and nu z momenta
	_c_reco_nuCtr    += 1;
	_c_reco_nuZ_diff += abs(neutrinos[0].mom().z() - recoNeutrino.z());
	_c_abs_nuZ       += abs(neutrinos[0].mom().z());
	_c_abs_reco_nuZ  += abs(recoNeutrino.z());
	
	// Check that nu reco has worked
	if (recoNeutrino.dot(recoNeutrino) == 0.0) {
	  //MSG_ERROR("VETO DUE TO BAD NEUTRINO RECO");
	  vetoEvent;
	}

	// convert reconstructed neutrino momenta to a 4vector with zero mass.
	FourMomentum pNeutrino;
	pNeutrino.setPM(recoNeutrino.x(), recoNeutrino.y(), recoNeutrino.z(), 0.0);

	const FourMomentum pW = lep + pNeutrino;
	const FourMomentum pTop = pW + bjet.mom();
        

        // Fill particle-level fiducial differential distributions
        size_t itop = (lep.charge() > 0) ? 0 : 1;
        //_h_AbsPtclDiffXsecTPt[itop]->fill( pTop.pT(),     event.weight());
	//_h_AbsPtclDiffXsecTY[itop] ->fill( pTop.absrap(), event.weight());
	_h_NrmPtclDiffXsecTPt[itop]->fill( pTop.pT()/GeV);
	_h_NrmPtclDiffXsecTY[itop] ->fill( pTop.absrap());
	//_h_AbsPtclDiffXsecJPt[itop]->fill( ljet.pT(),	  event.weight());
	//_h_AbsPtclDiffXsecJY[itop] ->fill( ljet.absrap(), event.weight());
	_h_NrmPtclDiffXsecJPt[itop]->fill( ljet.pT()/GeV);
	_h_NrmPtclDiffXsecJY[itop] ->fill( ljet.absrap());
	
      } while (false); //< just one iteration

    }


    /// Normalise histograms and print info after the run
    void finalize() {

      // some debug print out
      cout << "Average difference between leading MC neutrino" << endl;
      cout << "z momentum and reconstructed neutrino z momentum: "
	   << _c_reco_nuZ_diff / _c_reco_nuCtr << endl;
      cout << "Average abs value of real leading MC neutrino z : "
	   << _c_abs_nuZ / _c_reco_nuCtr << endl;
      cout << "Average abs value of reconstructed neutrino z   : "
	   << _c_abs_reco_nuZ / _c_reco_nuCtr << endl;

      

      //double lepBR = 0.324; //branching ratio of W->lepton
      double lepBR = 0.29;  // branching ratio of W->e/mu (via taus also included)

      //double A_fid = 0.175; // Roughly. From paper.
      for (int i = 0; i < 2; i++) {

      // Fiducial particle-level distributions (tops)
      //scale( _h_AbsPtclDiffXsecTPt[i], crossSection()/femtobarn/sumOfWeights() / lepBR );
	//scale( _h_AbsPtclDiffXsecTPt[i], crossSection()/femtobarn/sumOfWeights() );
	//scale( _h_AbsPtclDiffXsecTY[i],  crossSection()/sumOfWeights() / lepBR );
	normalize(_h_NrmPtclDiffXsecTPt[i]);
	normalize(_h_NrmPtclDiffXsecTY[i]);
	// Fiducial particle-level distributions (jets)
	//scale( _h_AbsPtclDiffXsecJPt[i], crossSection()/femtobarn/sumOfWeights() / lepBR );
	//scale( _h_AbsPtclDiffXsecJY[i],  crossSection()/sumOfWeights() / lepBR );
	normalize(_h_NrmPtclDiffXsecJPt[i]);
	normalize(_h_NrmPtclDiffXsecJY[i]);
	// Full phase-space parton level results (tops)
	//scale(_h_AbsPrtnDiffXsecTPt[i], crossSection()/ femtobarn/ sumOfWeights()/ lepBR );
	//scale(_h_AbsPrtnDiffXsecTY[i],  crossSection() / sumOfWeights() / lepBR );
	normalize(_h_NrmPrtnDiffXsecTPt[i]);
	normalize(_h_NrmPrtnDiffXsecTY[i]);
      }

      // Print generic cutflow info
      for (int i = 0; i < 9; i++) {
        //cout << "Cutflow at step " << i << ": " << cutflow[i]<<" ratio: "<<(float)cutflow[i]/(float)cutflow[0] << endl;
      }
      
      // Print detailed cutflow details
      cout << "Cutflow: exactly 1 e or mu (dressed): " << cutflow[1] <<" ratio: "<<(float)cutflow[1]/(float)cutflow[0] << endl;
      cout << "Cutflow: lepton eta < 2.5 cut       : " << cutflow[2] <<" ratio: "<<(float)cutflow[2]/(float)cutflow[0] << endl;
      cout << "Cutflow: lepton pt < 25*GeV cut     : " << cutflow[3] <<" ratio: "<<(float)cutflow[3]/(float)cutflow[0] << endl;
      cout << "Cutflow: N_btaggedjets != 1         : " << cutflow[4] <<" ratio: "<<(float)cutflow[4]/(float)cutflow[0] << endl;
      cout << "Cutflow: N_light_jets  != 1         : " << cutflow[5] <<" ratio: "<<(float)cutflow[5]/(float)cutflow[0] << endl;
      cout << "Cutflow: jet isolation (deltaR<0.4) : " << cutflow[6] <<" ratio: "<<(float)cutflow[6]/(float)cutflow[0] << endl;
      cout << "Cutflow: lep bjet mass cut          : " << cutflow[7] <<" ratio: "<<(float)cutflow[7]/(float)cutflow[0] << endl;
      cout << "Cutflow: Tmass_lmet cut             : " << cutflow[8] <<" ratio: "<<(float)cutflow[8]/(float)cutflow[0] << endl;
      cout << "Cutflow: lepton jet back-to-bck cut : " << cutflow[9] <<" ratio: "<<(float)cutflow[9]/(float)cutflow[0] << endl;

      
    }

    //@}



      

    /// @name Histograms
    //@{
    int _c_eventCtr = 0;
    int _c_reco_nuCtr = 0;
    double _c_reco_nuZ_diff = 0.0;
    double _c_abs_nuZ = 0.0;
    double _c_abs_reco_nuZ = 0.0;


    // CounterPtr _c_Xsec_fid_tq, _c_sumw_fid, _c_sumw_fid_tq, _c_sumw_prtn_tq;
    //Histo1DPtr _h_AbsPtclDiffXsecTPt[2], _h_AbsPtclDiffXsecTY[2], _h_NrmPtclDiffXsecTPt[2], _h_NrmPtclDiffXsecTY[2];
    Histo1DPtr _h_NrmPtclDiffXsecTPt[2], _h_NrmPtclDiffXsecTY[2];
    //Histo1DPtr _h_AbsPtclDiffXsecJPt[2], _h_AbsPtclDiffXsecJY[2], _h_NrmPtclDiffXsecJPt[2], _h_NrmPtclDiffXsecJY[2];
    Histo1DPtr _h_NrmPtclDiffXsecJPt[2], _h_NrmPtclDiffXsecJY[2];
    //Histo1DPtr _h_AbsPrtnDiffXsecTPt[2], _h_AbsPrtnDiffXsecTY[2], _h_NrmPrtnDiffXsecTPt[2], _h_NrmPrtnDiffXsecTY[2];
    Histo1DPtr _h_NrmPrtnDiffXsecTPt[2], _h_NrmPrtnDiffXsecTY[2];

    int cutflow[9];

    //@}


  };


  Vector3 reconstructNeutrino(Vector3& mmpt, FourMomentum& lepton) {

      float wMass = 80.399*GeV;
      bool is_complex = false;

      // get x and y momentum values for the reconstructed neutrino
      Vector3 reco_nu = getNuXY(mmpt, lepton, wMass, &is_complex);

      // get z momentum value for recostructed neutrino
      double mu  = (wMass*wMass)/2 + reco_nu.x()*lepton.x() + reco_nu.y()*lepton.y();
      double mu2 = pow(mu,2);
      double nupt2 = reco_nu.x()*reco_nu.x() + reco_nu.y()*reco_nu.y();

      // define new neutrino z
      double nuz = mu*lepton.z()/pow(lepton.pt(),2);

      // if there are one or two non-complex solutions then over-write the
      // new nu z momentum value ('nuz') given above.
      if (!is_complex){
	double second_term = mu2*pow(lepton.z(),2)/pow(lepton.pt(),4) - (pow(lepton.E(),2)*nupt2 - mu2)/pow(lepton.pt(),2);
	
	if(second_term < 0) {
	  //MSG_ERROR("ERROR: COMPLEX SOLUTION - THIS SHOULDN'T HAPPEN!");
	} else {
	  
          double root = sqrt(second_term);

	  // take smallest val (from analysers original rivet code, not mentioned in int note)
	  double soln1 = nuz + root;
	  double soln2 = nuz - root;
	  if (abs(soln1) <= abs(soln2)) nuz = soln1; else nuz = soln2;
	  //MSG_DEBUG("solutions for nu z: " << soln1 << ", " << soln2 << " (choosing smaller abs val)");
	}
      }

      
      reco_nu.setZ(nuz);
      return reco_nu;

  }


  
  Vector3 getNuXY(Vector3 & neutrino, FourMomentum& lepton, float wMass, bool *is_complex) {

      Vector3 reco_nu;
      vector<float> pz;


      // Calculate transverse mass of W-boson
      double MisET2 = neutrino.x()*neutrino.x() + neutrino.y()*neutrino.y();
      double misET  = sqrt(MisET2);
      double mWT = sqrt(pow(lepton.pt() + misET,2) - pow(lepton.x() + neutrino.x(),2) - pow(lepton.y() + neutrino.y(),2) );


      if (mWT > wMass){	/// Redefine neutrino x and y
	*is_complex = true;
	// Define a null particle
	Vector3 nullP;
	nullP.setX(0.0) ; nullP.setY(0.0) ; nullP.setZ(0.0);
	
	cout << "Complex solution found (i.e wMass_T > wMass_pole), evaluating x and y..." << endl;

	double upper = 0.0;
	double lower = 0.0;;

	
	// Restrict range of neutrino x parameter
	if(lepton.x() < 0) {
	  upper = - wMass*wMass/(4*lepton.x()) -0.01;
	  lower = -9999.;
	}

	else if(lepton.x() == 0) {
	  upper =  9999.;
	  lower = -9999.;
	}

	else {
	  upper = 9999.;
	  lower = - wMass*wMass/(4*lepton.x()) +0.01;
	}

	//cout << "neutrino range: [" << lower << "," << upper << "]" << endl; 

	// Make an initial guess
	/*
	if(mmpt.x() > upper) initial = upper -1;
	else if(mmpt.x() < lower) initial = lower + 1;
	else initial = mmpt.x();
	*/

	bool x1_nan = false;
	bool x2_nan = false;

	
	// Define two possible solutions for neutrino x momentum
	double nux1 = 0.0;
	double nux2 = 0.0;


	// Define fixed parameters for mendelmin minimizer
	valarray<double> fixed_params{lepton.pt(),lepton.x(),lepton.y(), neutrino.x(), neutrino.y(), lower, upper};

	
	// shout if denominators in minamizable function are near zero	
	if (abs(lepton.x()) < 0.0001){
	  
	  cout << "lepton.x() near zero" << endl;
	  return nullP;
	}

	if ( abs((2*sqr(fixed_params[1])) - fixed_params[4]) < 0.0000001 ){
	  cout << "denominator near zero" << 2*sqr(fixed_params[1]) - fixed_params[4] << endl;
	  return nullP;
	}

	/// Compute delta1 and delta2 which will be minimized. The neutrino x mom used
	/// in the smaller delta will become the new neutrino x mom..

	// delta1
	MendelMin mm1(target1, fixed_params, rand01, 1);
	//const double best1 = mm1.evolve(200);
	mm1.evolve(200);
	valarray<double> fittest1 = mm1.fittest();
	double delta1 = target1(fittest1, fixed_params);
	nux1 = fittest1[0] * (upper - lower) + lower;
 	cout << "delta1 = " << delta1 << " (final x value: " << nux1 << ")" << endl;

	// delta2
	MendelMin mm2(target2, fixed_params, rand01, 1);
	//const double best2 = mm2.evolve(200);
	mm2.evolve(200);
	valarray<double> fittest2 = mm2.fittest();
	double delta2 = target2(fittest2, fixed_params);
	nux2 = fittest2[0] * (upper - lower) + lower;
	cout << "delta2 = " << delta2 << " (final x value: " << nux2 << ")" << endl;



	//cout << "nux1=" << nux1 << endl;
	//cout << "nux2=" << nux2 << endl;;

	
	double new_nux  = 0.0;
	double new_nuy  = 0.0;
	double new_nuy1 = 0.0;
	double new_nuy2 = 0.0;
	bool plus;
	
	if (delta1 <= delta2){
	  plus = true;
	  
	
	  new_nux = nux1; 
	  new_nuy1 = nuy1(nux1, lepton.x(), lepton.y(), lepton.pt());
	  new_nuy2 = nuy2(nux1, lepton.x(), lepton.y(), lepton.pt());
	  new_nuy = new_nuy1;

	  
	} else {
	  plus = false;

	  new_nux = nux2; 
	  new_nuy1 = nuy1(nux2, lepton.x(), lepton.y(), lepton.pt());
	  new_nuy2 = nuy2(nux2, lepton.x(), lepton.y(), lepton.pt());
	  new_nuy = new_nuy2;
	}

	cout << "recalculated y: " << new_nuy << endl;
	// set values of reconstructed neutrino x and y.
	reco_nu.setX(new_nux);
	reco_nu.setY(new_nuy);

	
	// print some info
	// cout << "(nux,nuy1)  = ("<< new_nux      <<","<< new_nuy1     <<")"<< endl;
	// cout << "(nux,nuy2)  = ("<< new_nux      <<","<< new_nuy2     <<")"<< endl;
	// cout << "(metx,mety) = ("<< neutrino.x() <<","<< neutrino.y() <<")"<< endl;
	// cout << "new (x,y,z)   = ("<< new_nux      <<", "<< new_nuy      <<", "<< reco_nu.z() << ")" << endl;

	// check for inconsistancies in execution of method
	if ( abs(neutrino.y() - new_nuy1) > abs(neutrino.y() - new_nuy2) && (plus) ) cout << "ERROR: this should not happen: delta1 mistakenly selected" << endl;

	if ( abs(neutrino.y() - new_nuy2) > abs(neutrino.y() - new_nuy1) && (!plus)) cout << "ERROR: this should not happen: delta2 mistakenly selected" << endl;
      } else {
	reco_nu.setX(neutrino.x());
	reco_nu.setY(neutrino.y());
      }
      return reco_nu;    
  }

  
  double target1(const MendelMin::Params& x_1, const MendelMin::Params& p) {
    
    double m_W = 80.4;
    
    // the function to minimize
    //return sqrt(pow( x_1[0]-p[3] ,2) + pow( (sqr(m_W)*p[2] + 2*p[1]*p[2]*x_1[0] + m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x_1[0]))/(2*sqr(p[1])) - p[4],2));
    
    // the function to minimize (with scaling s.t. x_1 -->  lower < x < upper
    double x = x_1[0]*(p[6]-p[5]) + p[5];
    
    return sqrt(pow( x-p[3] ,2) + pow( (sqr(m_W)*p[2] + 2*p[1]*p[2]*x + m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x))/(2*sqr(p[1])) - p[4],2));
    
  }
  
  
  double target2(const MendelMin::Params& x_2, const MendelMin::Params& p) {
    
    double m_W = 80.4; // TODO improve this accuracy
    
    // the function to minimize
    //return sqrt(pow( x_2[0]-p[3] ,2) + pow( (sqr(m_W)*p[2] + 2*p[1]*p[2]*x_2[0] - m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x_2[0]))/(2*sqr(p[1])) - p[4],2));
    
    
    // the function to minimize (with scaling s.t. x_2 -->  lower < x < upper
    double x = x_2[0]*(p[6]-p[5]) + p[5];
    
    return sqrt(pow( x-p[3] ,2) + pow( (sqr(m_W)*p[2] + 2*p[1]*p[2]*x - m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x))/(2*sqr(p[1])) - p[4],2));
    
  }
  
  double rand01() {
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
  }
  
  double nuy1(const double nux, const double lepx, const double lepy, const double leppt){
    
    double nuy = 0.0;
    double m_W = 80.4; // TODO improve this accuracy
    
    nuy = (sqr(m_W)*lepy + 2*lepx*lepy*nux + m_W*leppt*sqrt(sqr(m_W) + 4*lepx*nux))/(2*sqr(lepx));
    
    return nuy;
  }
  
  double nuy2(const double nux, const double lepx, const double lepy, const double leppt){
    
    double nuy = 0.0;
    double m_W = 80.4; // TODO improve this accuracy
    
    nuy = (sqr(m_W)*lepy + 2*lepx*lepy*nux - m_W*leppt*sqrt(sqr(m_W) + 4*lepx*nux))/(2*sqr(lepx));
    
    return nuy;
  }
  
  
  
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1512776);
  
}
