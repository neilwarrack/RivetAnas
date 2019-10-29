#include "multimin.h"
#include <cmath>

inline double sqr(double x) { return x*x; }

double g(double x, double y) {
  // const double r2 = sqr(x-3) + sqr(y-5);
  // const double r = sqrt(r2);
  // const double cos2r = sqr(cos(r));
  // return -cos2r/(r + 1e-6);
  return 10*sqr(x-1.5) + 20*sqr(y-2.5) + 30;
  
}


void f(const size_t n, const double *x,void *fparams,double *fval){
  double *p = (double *) fparams;
  // *fval=p[2]*(x[0]-p[0])*(x[0]-p[0]) + p[3]*(x[1]-p[1])*(x[1]-p[1]) + p[4];
  // *fval = g(x[0], x[1]);

  double m_W = 80.4;
  //form:
  //nuy = (sqr(m_W)*ly + 2*lx*ly*nux + m_W*lt*sqrt(sqr(m_W) + 4*lx*nux))/2*sqr(lx);

  *fval = (sqr(m_W)*x[0] + 2*p[1]*p[2]*x[0] + m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x[0]))/2*sqr(p[1]);
  
}

void df(const size_t n, const double *x,void *fparams,double *grad)
{
  double *p = (double *) fparams;
  // grad[0]=2*p[2]*(x[0]-p[0]);
  // grad[1]=2*p[3]*(x[1]-p[1]);
  /*
  const double delta = 1e-6;
  const double f0 = g(x[0], x[1]);
  const double f1x = g(x[0]+delta, x[1]);
  const double f1y = g(x[0], x[1]+delta);
  grad[0] = (f1x-f0)/delta;
  grad[1] = (f1y-f0)/delta;
  */
  double m_W = 80.4;
  grad[0] = sqr(m_W)*p[2]/(2*sqr(p[1])) + p[2]/p[1] + (m_W*p[0]*p[1])/(sqr(p[1])*sqrt(sqr(m_W) + 4*p[1]*x[0]));

}

void fdf(const size_t n, const double *x,void *fparams,double *fval,double *grad)
{
  // double *p = (double *) fparams;
  // *fval=p[2]*(x[0]-p[0])*(x[0]-p[0]) + p[3]*(x[1]-p[1])*(x[1]-p[1]) + p[4];
  // grad[0]=2*p[2]*(x[0]-p[0]);
  // grad[1]=2*p[3]*(x[1]-p[1]);
  f(n, x, fparams, fval);
  df(n, x, fparams, grad);
}


int main(){

  double par[3] = { 1.5, 2.5, 10.0};

  double x[1]={0};
  double minimum;

  double xmin[2], xmax[2];
  unsigned type[2];

  struct multimin_params optim_par = {.1,1e-2,100,1e-3,1e-5,2,0};

  /* unconstrained minimization */
  multimin(2,x,&minimum,NULL,NULL,NULL,&f,&df,&fdf,(void *) par,optim_par);

  printf("unconstrained minimum %e at x=%e\n",minimum,x[0]);

  // /* minimum constrained in [-1,2]x(-5,5] */
  // type[0]=3;
  // xmin[0]=-1;
  // xmax[0]=2;
  // x[0]=0;

  // type[1]=2;
  // xmin[1]=-5;
  // xmax[1]=5;
  // x[1]=0;

  // /* increase verbosity */
  // optim_par.verbosity=1;

  // multimin(2,x,&minimum,type,xmin,xmax,&f,&df,&fdf,par,optim_par);

  // printf("constrained minimum %e at x=%e y=%e\n",minimum,x[0],x[1]);


  return 0;
}
