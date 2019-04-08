#if !defined(NEUTRINO)
#define NEUTRINO


namespace NeutrinoFit {

static Double_t FullReco_MET_X = 0.;
static Double_t FullReco_MET_Y = 0.;
static Double_t FullReco_Pxe = 0.;
static Double_t FullReco_Pye = 0.;
static Double_t FullReco_PTe = 0.;


static void delta1fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  double mW = 80.4;  // PDG 2006 value
  double delta = 0.0; 

  delta = sqrt((par[0]-FullReco_MET_X)*(par[0]-FullReco_MET_X) + TMath::Power((((mW*mW*FullReco_Pye + 2*FullReco_Pxe*FullReco_Pye*par[0])-(mW*FullReco_PTe)*(sqrt(mW*mW + 4*FullReco_Pxe*par[0])))/(2*FullReco_Pxe*FullReco_Pxe) - FullReco_MET_Y),2));

  f = delta;

}


static void delta2fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
 
  double mW = 80.4; // PDG 2006 value
  double delta = 0.0; 

  delta = sqrt((par[0]-FullReco_MET_X)*(par[0]-FullReco_MET_X) + TMath::Power((((mW*mW*FullReco_Pye + 2*FullReco_Pxe*FullReco_Pye*par[0])+(mW*FullReco_PTe)*(sqrt(mW*mW + 4*FullReco_Pxe*par[0])))/(2*FullReco_Pxe*FullReco_Pxe) - FullReco_MET_Y),2));
  

  f = delta;

}
}
#endif
