
inline double sqr(double x) { return x*x; }

double g(double x, double y) {
  // const double r2 = sqr(x-3) + sqr(y-5);
  // const double r = sqrt(r2);
  // const double cos2r = sqr(cos(r));
  // return -cos2r/(r + 1e-6);
  return 10*sqr(x-1.5) + 20*sqr(y-2.5) + 30;
}


void f(const size_t n, const double *x,void *fparams,double *fval){
  // double *p = (double *) fparams;
  // *fval=p[2]*(x[0]-p[0])*(x[0]-p[0]) + p[3]*(x[1]-p[1])*(x[1]-p[1]) + p[4];
  *fval = g(x[0], x[1]);
}

void df(const size_t n, const double *x,void *fparams,double *grad)
{
  // double *p = (double *) fparams;
  // grad[0]=2*p[2]*(x[0]-p[0]);
  // grad[1]=2*p[3]*(x[1]-p[1]);
  const double delta = 1e-6;
  const double f0 = g(x[0], x[1]);
  const double f1x = g(x[0]+delta, x[1]);
  const double f1y = g(x[0], x[1]+delta);
  grad[0] = (f1x-f0)/delta;
  grad[1] = (f1y-f0)/delta;
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
