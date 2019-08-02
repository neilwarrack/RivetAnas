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

// #include "neutrino.h"
//using namespace std;
using namespace Rivet; // TODO what to do with this line??
//#include "NeilsMendelMin.h"
#include "MendelMin.h"
#include "multimin.h"
#include "minimizer.h"


inline double sqrd(double x) { return x*x; }

double AndysTarget(const MendelMin::Params& p, const MendelMin::Params& q) {
  const double x = p[0] - 0.3;
  const double y = p[1] - 0.5;
  
  //return (x*x + 2*y*y) + 0.1*sin(x)*sin(y);
  return (x*x + 2*y*y) + 0.1*sin(x)*sin(y) + q[0] - q[0];
}


double target1(const MendelMin::Params& x_1, const MendelMin::Params& p) {

  double m_W = 80.4;

  // the function to minimize
  //return sqrt(pow( x_1[0]-p[3] ,2) + pow( (sqrd(m_W)*p[2] + 2*p[1]*p[2]*x_1[0] + m_W*p[0]*sqrt(sqrd(m_W) + 4*p[1]*x_1[0]))/(2*sqrd(p[1])) - p[4],2));

  // the function to minimize (with scaling s.t. x_1 -->  lower < x < upper
  double x = x_1[0]*(p[6]-p[5]) + p[5];

  return sqrt(pow( x-p[3] ,2) + pow( (sqrd(m_W)*p[2] + 2*p[1]*p[2]*x + m_W*p[0]*sqrt(sqrd(m_W) + 4*p[1]*x))/(2*sqrd(p[1])) - p[4],2));

}


double target2(const MendelMin::Params& x_2, const MendelMin::Params& p) {

  double m_W = 80.4; // TODO improve this accuracy

  // the function to minimize
  //return sqrt(pow( x_2[0]-p[3] ,2) + pow( (sqrd(m_W)*p[2] + 2*p[1]*p[2]*x_2[0] - m_W*p[0]*sqrt(sqrd(m_W) + 4*p[1]*x_2[0]))/(2*sqrd(p[1])) - p[4],2));


  // the function to minimize (with scaling s.t. x_2 -->  lower < x < upper
  double x = x_2[0]*(p[6]-p[5]) + p[5];

  return sqrt(pow( x-p[3] ,2) + pow( (sqrd(m_W)*p[2] + 2*p[1]*p[2]*x - m_W*p[0]*sqrt(sqrd(m_W) + 4*p[1]*x))/(2*sqrd(p[1])) - p[4],2));

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

  nuy = (sqrd(m_W)*lepy + 2*lepx*lepy*nux + m_W*leppt*sqrt(sqrd(m_W) + 4*lepx*nux))/(2*sqrd(lepx));

  return nuy;
}

double nuy2(const double nux, const double lepx, const double lepy, const double leppt){

  double nuy = 0.0;
  double m_W = 80.4; // TODO improve this accuracy

  nuy = (sqrd(m_W)*lepy + 2*lepx*lepy*nux - m_W*leppt*sqrt(sqrd(m_W) + 4*lepx*nux))/(2*sqrd(lepx));

  return nuy;
}

namespace Rivet {
gsl_error_handler_t * gsl_set_error_handler_off();

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
      declare(PartonicTops(PartonicTops::ALL     ),  "AllPartonicTops");
      declare(PartonicTops(PartonicTops::E_MU    ),  "emuLeptonicTops");
      declare(PartonicTops(PartonicTops::E_MU_TAU),     "LeptonicTops");


      // Leptons, jets and MET
      //declare(DressedLeptons(PromptFinalState(Cuts::abseta < 2.6), 0.1, Cuts::abseta < 2.5 && Cuts::pT > 25*GeV), "Leptons");
      //DressedLeptons dressedLeps(PromptFinalState(Cuts::abseta < 2.6), 0.1, Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);


      // direct copy (with lines missing) from single top groups rivet///////////////////
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
      Cut lep_cuts = Cuts::abseta < 2.5 
	&& (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      PromptFinalState bareLeptons(lep_cuts,true); // 'true' accepts tau decays
      IdentifiedFinalState photons(Cuts::abseta < 2.6 && Cuts::abspid == PID::PHOTON);
      DressedLeptons dressedLeptons(photons, bareLeptons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 25*GeV);
      declare(dressedLeptons, "Leptons");

      DressedLeptons vetoDressedLeptons(photons, bareLeptons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 10*GeV);

      // Jets
      VetoedFinalState vfs;
      vfs.vetoNeutrinos();
      vfs.addVetoOnThisFinalState(dressedLeptons);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "Jets");
      //////// END ///////////////////////////////////////////////////////////////////////


      //vfs.addVetoOnThisFinalState(dressedLeps);
      
      //      declare(MissingMomentum(FinalState(Cuts::abseta < 4.5)), "MET");
      declare(MissingMomentum(FinalState()), "MET");


      // Book histograms and counters
      // _c_sumw_fid =         bookCounter("sumw_fid_tot");
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

      for(int i = 0; i < 9; i++) cutflow[i] = 0;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      _c_eventCtr += 1;
      //if (_c_eventCtr < 2434){vetoEvent;}
      // PARTON-LEVEL ANALYSIS
      do { // trick to allow early exit from selection without full return

	//const Particles& allTops = apply<PartonicTops>(event, "AllPartonicTops").particles();
	const Particles& allTops = apply<PartonicTops>(event, "emuLeptonicTops").particles();
        if (allTops.size() != 1) break;
        const Particle& top = allTops[0];

        size_t itop = (top.charge() > 0) ? 0 : 1;
	_h_AbsPrtnDiffXsecTPt[itop]->fill( top.pT(),     event.weight());
	_h_AbsPrtnDiffXsecTY[itop] ->fill( top.absrap(), event.weight());
	_h_NrmPrtnDiffXsecTPt[itop]->fill( top.pT(),     event.weight());
	_h_NrmPrtnDiffXsecTY[itop] ->fill( top.absrap(), event.weight());

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
        MSG_DEBUG("  Passed jet isolation cuts");
	cutflow[4]++;

        // Mlb cut
        const FourMomentum plb = lep.mom() + bjet.mom();
        if (plb.mass() > 160*GeV) break;
        MSG_DEBUG("  Passed Mlb selection");
        cutflow[5]++;

	// transverse mass cut
	const MissingMomentum& mm = apply<MissingMomentum>(event, "MET");
	
	double delta_phi_lep_met = lep.phi()-mm.missingMom().phi();
	double tmasslepmet = sqrt(2*lep.pT()*mm.met()*1-cos(delta_phi_lep_met));
	if (tmasslepmet < 50*GeV) break;
	MSG_DEBUG("  Passed transverse mass selection");
	cutflow[6]++;

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
	cutflow[7]++;


	// ---------------- top quark reconstruction --------------- //
	
	FourMomentum lep_4p = lep.mom();
	Vector3 mmpt = mm.vectorPt();
	double mmz = mm.missingMom().z();

	// reconstruct neutrino
	FourMomentum recoNeutrino = reconstructNeutrino(mmpt, lep_4p, mmz);

	const FourMomentum pW = lep + recoNeutrino;
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
      
      scale( _h_AbsPtclDiffXsecTPt, crossSection()/femtobarn/sumOfWeights() / lepBR );
      scale( _h_AbsPtclDiffXsecTY,  crossSection()/sumOfWeights() / lepBR );
      normalize(_h_NrmPtclDiffXsecTPt);
      normalize(_h_NrmPtclDiffXsecTY);
      //
      scale( _h_AbsPtclDiffXsecJPt, crossSection()/femtobarn/sumOfWeights() / lepBR );
      scale( _h_AbsPtclDiffXsecJY,  crossSection()/sumOfWeights() / lepBR );
      normalize(_h_NrmPtclDiffXsecJPt);
      normalize(_h_NrmPtclDiffXsecJY);
      //
      scale(_h_AbsPrtnDiffXsecTPt, crossSection()/femtobarn/sumOfWeights()/lepBR );
      scale(_h_AbsPrtnDiffXsecTY,  crossSection() / sumOfWeights() / lepBR );
      normalize(_h_NrmPrtnDiffXsecTPt);
      normalize(_h_NrmPrtnDiffXsecTY);

      for (int i = 0; i < 9; i++) {
        cout << "Cutflow at step " << i << ": " << cutflow[i]<<" Afid: "<<(float)cutflow[i]/(float)cutflow[0] << endl;
      }

      
    }

    //@}



    FourMomentum reconstructNeutrino(Vector3& mmpt, FourMomentum& lepton, double z_before) {

      float wMass = 80.399*GeV;
      FourMomentum return_neutrino;

      Vector3 neutrino;
      neutrino.setX(mmpt.x());
      neutrino.setY(mmpt.y());
      
      vector<float> pZ = getNeutrinoPzSolutions(mmpt, lepton, wMass);
      //cout << "pz=" << pZ << endl;      
      int nSolutions = pZ.size();
      if(nSolutions == 2) {
        neutrino.setZ( std::fabs(pZ[0]) < std::fabs(pZ[1]) ? pZ[0] : pZ[1] );
      } else if(nSolutions == 0) {
	cutflow[8]++;
	cout << "no neutrino z-solution!" << endl;
        // float mTSq = 2. * (neutrino.mod()*lepton.pT()
	// 		   - neutrino.x()*lepton.x()
	// 		   - neutrino.y()*lepton.y());
	// neutrino *= wMass*wMass/mTSq; // Scale down pt(nu) such that there is exactly 1 solution
        // neutrino.setZ(neutrino.mod()/lepton.pT()*lepton.pz()); // p_nuT adjustment approach


	

	
	//cout << "  z of missing momentum:" << z_before << endl;
	//cout << "scaled pT of nu: p_z = " << neutrino.z() << " (Before: " << z_before << ")" << endl;


	// minimizer...................................................
	
	double m_W = 80.4;
	double upper[1];
	double lower[1];
	double initial;
	
	if(lepton.x() < 0) {
	  upper[0] = - m_W*m_W/(4*lepton.x()) -0.01;
	  lower[0] = -9999.;
	}

	else if(lepton.x() == 0) {
	  upper[0] =  9999.;
	  lower[0] = -9999.;
	}

	else {
	  upper[0] = 9999.;
	  lower[0] = - m_W*m_W/(4*lepton.x()) +0.01;
	}

	if(mmpt.x() > upper[0]) initial = upper[0] -1;
	else if(mmpt.x() < lower[0]) initial = lower[0] + 1;
	else initial = mmpt.x();


	double delta_1, delta_2;
	bool x1_nan = false;
	bool x2_nan = false;


	/// define parameters
	double par[5] = { lepton.pt(), lepton.x(),lepton.y(), mmpt.x(), mmpt.y()};

	/// configure multimin
	// NB: 2nd to last integer argument n uses gradient functions
	// if (n<4 || n==5)... see multimin.c line 197
	struct multimin_params optim_par = {.1,1e-2,1000,1e-3,1e-5,7,0};
	

	/// shout if denominators are near zero	
	if (abs(par[1]) < 0.01){
	  cout << "lepton.x() near zero" << endl;
	}

	if ( abs((2*sqrd(par[1])) - par[4]) < 0.01 ){
	  cout << "denominator '2' are near zero" << 2*sqrd(par[1]) - par[4] << endl;
	}
	

	double x_1[1];
	double x_2[1];
	x_1[0] = initial;
	x_2[0] = initial;


	double p[5]={lepton.pt(),lepton.x(),lepton.y(), mmpt.x(), mmpt.y()};

	//                                                                                     |
	double fnVal_1 = sqrt(pow( x_1[0]-p[3] ,2) + pow( (sqrd(m_W)*p[2] + 2*p[1]*p[2]*x_1[0] + m_W*p[0]*sqrt(sqrd(m_W) + 4*p[1]*x_1[0]))/(2*sqrd(p[1])) - p[4],2));
	//                                                                                     |
	double fnVal_2 = sqrt(pow( x_2[0]-p[3] ,2) + pow( (sqrd(m_W)*p[2] + 2*p[1]*p[2]*x_2[0] - m_W*p[0]*sqrt(sqrd(m_W) + 4*p[1]*x_2[0]))/(2*sqrd(p[1])) - p[4],2));


	//_________________MendelMin START________________
	cout << "___MENDLMIN_START____" << endl;
	double nux1 = 0.0;
	double nux2 = 0.0;

	valarray<double> fixed_params{lepton.pt(),lepton.x(),lepton.y(), mmpt.x(), mmpt.y(), lower[0], upper[0]};
	//std::vector<MendelMin::Params> fixed_params = {lepton.pt(),lepton.x(),lepton.y(), mmpt.x(), mmpt.y()};;

	// DELTA 1
	MendelMin mm1(target1, fixed_params, rand01, 1);
	const double best1 = mm1.evolve(20);

	valarray<double> fittest1 = mm1.fittest();
	double delta1_mend = target1(fittest1, fixed_params);
 	cout << "NEIL:FOUND ANSWER 1: " << delta1_mend << " == " << best1 << "; ";
	for (double x : fittest1){
	  nux1 = x*(upper[0] - lower[0]) + lower[0];;
	  cout << "rescaled nu x (1) = " << nux1;
	}
	cout << endl;
	
	// DELTA 2
	MendelMin mm2(target2, fixed_params, rand01, 1);
	const double best2 = mm2.evolve(20);

	valarray<double> fittest2 = mm2.fittest();
	double delta2_mend = target2(fittest2, fixed_params);
	cout << "NEIL:FOUND ANSWER 2: " << delta2_mend << " == " << best2 << "; ";
	for (double x : fittest2){
	  nux2 = x*(upper[0] - lower[0]) + lower[0];
	  cout << "rescaled nu x (2) = " << nux2;
	}
	cout << endl;

	cout << "nux1=" << nux1 << endl;
	cout << "nux2=" << nux2 << endl;;

	double new_nux = 0.0;
	double new_nuy = 0.0;

	
	if (delta1_mend <= delta2_mend){
	  cout << "delta1_mend <= delta2_mend, choosing nu_x 1" << endl;
	  new_nux = nux1; 
	  //nuy = nuy(nux, lepx, lepy, leppT);
	  //lepton.pt(),lepton.x(),lepton.y(), mmpt.x(), mmpt.y(), lower[0], upper[0]};
	  new_nuy = nuy1(nux1, lepton.x(), lepton.y(), lepton.pt());

	} else {
	  cout << "delta1_mend > delta2_mend, choosing nu_x 2" << endl;
	  new_nuy = nuy2(nux2, lepton.x(), lepton.y(), lepton.pt());
  	  new_nux = nux2; 
	}

	cout << "(nux,nuy)   = (" << new_nux << ","  << new_nuy  << ")" << endl;
	cout << "(metx,mety) = (" << mmpt.x() << "," << mmpt.y() << ")" << endl;
	//_________________________
	// Andy's stuff (updated to take fixed_params and not use them)
	MendelMin Amm(AndysTarget, fixed_params, rand01, 2);
	const double Abest = Amm.evolve(20);

	valarray<double> Afittest = Amm.fittest();
	cout << "ANDY:FOUND ANSWER: " << AndysTarget(Afittest, fixed_params) << " == " << Abest << "; ";
	for (double x : Afittest)
	  cout << x << " ";
	cout << endl;

	valarray<double> Aright{0.3, 0.5};
	cout << "ANDY:RIGHT ANSWER: " << AndysTarget(Aright, fixed_params) << "; ";
	for (double x : Aright)
	  cout << x << " ";
	cout << endl;


	bool print_everything = true;

	//_________________MendelMin END__________________

	if (_c_eventCtr == 2435){
	  cout << "event No.:2435" << endl;
	  cout << lepton.pt() << ", " << lepton.x() << ", " << lepton.y() << ", " << mmpt.x() << ", " << mmpt.y() << endl;
	  cout << "fnVal_1 =" << fnVal_1 << endl;
	  cout << "fnVal_2 =" << fnVal_2 << endl;
	}
	
	if (_c_eventCtr > 2435){
	  cout << "event:" << _c_eventCtr << endl;
	}

	multimin(1,x_1,&delta_1,NULL,&lower[0],&upper[0],&f,&df,&fdf,(void *) par,optim_par);
	if (_c_eventCtr == 2435){cout << "you are here!"<< endl;}
	if (x_1[0] != x_1[0]) {
	  x1_nan = true;
	  cout << "  function value = " << fnVal_1  << endl;
	  cout << "  met.x()        = " << mmpt.x() << endl;
	  cout << "  minimum        = " << delta_1  << endl;
	  cout << "  initial        = " << initial  << endl;
	  cout << "  final          = " << x_1[0]   << endl;
	  cout << "  lower          = " << lower    << endl;
	  cout << "  upper          = " << upper    << endl;
	}
	
	multimin(1,x_2,&delta_2,NULL,&lower[0],&upper[0],&f2,&df2,&f2df2,(void *) par,optim_par);

	if (x_2[0] != x_2[0]) {
	  x2_nan = true;
	  cout << "  function value = " << fnVal_2  << endl;
	  cout << "  met.x()        = " << mmpt.x() << endl;
	  cout << "  minimum        = " << delta_2  << endl;
	  cout << "  initial        = " << initial  << endl;
	  cout << "  final          = " << x_2[0]   << endl;
	  cout << "  lower          = " << lower    << endl;
	  cout << "  upper          = " << upper    << endl;
	}
	
	if (x2_nan && x1_nan) {cout << "############ TWO NANS!!!!!!########" << endl;}
	if (!x2_nan && !x2_nan){
	  if (abs(x_1[0]) <= abs(x_2[0])){
	    cout << "x_1[0]=" << x_1[0] << " is smaller than x_2[0]=" << x_2[0] << endl;
	    //neutrino.setZ()
	  } else {
	    cout << "x_2[0]=" << x_2[0] << " is smaller than x_1[0]=" << x_1[0] << endl;
	  }
	}

	if (print_everything){
	  cout << "  functionValue1 = " << fnVal_1  << endl;
	  cout << "  functionValue2 = " << fnVal_2  << endl;
	  cout << "  met.x()        = " << mmpt.x() << endl;
	  cout << "  delta_1        = " << delta_1  << endl;
	  cout << "  delta_2        = " << delta_2  << endl;
	  cout << "  initial        = " << initial  << endl;
	  cout << "  final          = " << x_1[0]   << endl;
	  cout << "  lower          = " << lower    << endl;
	  cout << "  upper          = " << upper    << endl;
	}
	// ..........................................................
	
      } else if(nSolutions == 1) {
        neutrino.setZ(pZ[0]);
      }
      return_neutrino.setPM(neutrino.x(), neutrino.y(), neutrino.z(), 0.0);
      return return_neutrino;
    }
      

    vector<float> getNeutrinoPzSolutions(Vector3 & neutrino, FourMomentum& lepton, float wMass) {

      vector<float> pz;
      
      float alpha = 0.5 * wMass * wMass + lepton.x() * neutrino.x() + lepton.y() * neutrino.y();
      float pT_lep2 = lepton.perp2();
      float discriminant = lepton.vector3().mod2() * (alpha * alpha - pT_lep2 * neutrino.mod2());
      if (discriminant < 0.){

	return pz;
      }
      //cout << "non-complex solutions found" << endl;
      
      float pz_offset = alpha * lepton.z() / pT_lep2;
      
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


    



    /// @name Histograms
    //@{
    int _c_eventCtr = 0;
    // CounterPtr _c_Xsec_fid_tq, _c_sumw_fid, _c_sumw_fid_tq, _c_sumw_prtn_tq;
    Histo1DPtr _h_AbsPtclDiffXsecTPt[2], _h_AbsPtclDiffXsecTY[2], _h_NrmPtclDiffXsecTPt[2], _h_NrmPtclDiffXsecTY[2];
    Histo1DPtr _h_AbsPtclDiffXsecJPt[2], _h_AbsPtclDiffXsecJY[2], _h_NrmPtclDiffXsecJPt[2], _h_NrmPtclDiffXsecJY[2];
    Histo1DPtr _h_AbsPrtnDiffXsecTPt[2], _h_AbsPrtnDiffXsecTY[2], _h_NrmPrtnDiffXsecTPt[2], _h_NrmPrtnDiffXsecTY[2];

    int cutflow[9];

    //@}


  };


  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1512776);

}
