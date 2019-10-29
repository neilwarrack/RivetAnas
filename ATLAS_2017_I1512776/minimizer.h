/*
the function:
p_{y,nu}(p_{x,nu}) = (sqr(m_W)*p[2] + 2*p[1]*p[2]*x[0] +/- m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x[0]))/2*sqr(p[1])
its deriv:

deiv = (p[2]/p[1]) +/- (m_W*p[0]*p[1]/(sqr(p[1])*sqrt(sqr(m_W) + 4*p[1]*x[0])))
(sqr(m_W)*p[2]/(2*sqr(p[1])) + p[2]/p[1] + (m_W*p[0]*p[1])/(sqr(p[1])*sqrt(sqr(m_W) + 4*p[1]*x[0])))
*/

inline double sqr(double x) { return x*x; }

/* const double MenFn(){ */
/*   double m_W = 80.4; */
/*   double currentVal = sqrt(pow( x_1[0]-p[3] ,2) + pow( (sqr(m_W)*p[2] + 2*p[1]*p[2]*x_1[0] + m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x_1[0]))/(2*sqr(p[1])) - p[4],2)); */

/*   return currentVal; */
/* } */

const int MenFn_simple(int x){

  double currentVal = x*x;

  return currentVal;
}

const double randan(double s){
  double evalRand = rand()/double(RAND_MAX);
  return evalRand;
    }





void f(const size_t n, const double *x,void *fparams,double *fval){
  double *p = (double *) fparams;
  // *fval=p[2]*(x[0]-p[0])*(x[0]-p[0]) + p[3]*(x[1]-p[1])*(x[1]-p[1]) + p[4];
  // *fval = g(x[0], x[1]);
  
  double m_W = 80.4;
  //form:
  //nuy = (sqr(m_W)*ly + 2*lx*ly*nux + m_W*lt*sqrt(sqr(m_W) + 4*lx*nux))/2*sqr(lx);


  
  *fval = sqrt(pow( x[0]-p[3] ,2) + pow( (sqr(m_W)*p[2] + 2*p[1]*p[2]*x[0] + m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x[0]))/(2*sqr(p[1])) - p[4],2));
  
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
  grad[0] = 2*(x[0] - p[3])*x[0] + 2*((sqr(m_W)*p[2] + 2*p[1]*p[2]*x[0] + m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x[0]))/2*sqr(p[1]) - p[4]) * ((p[2]/p[1]) + (m_W*p[0]*p[1]/(sqr(p[1])*sqrt(sqr(m_W) + 4*p[1]*x[0])))) / (2*sqrt( pow(x[0]-p[3],2) + pow((sqr(m_W)*p[2] + 2*p[1]*p[2]*x[0] + m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x[0]))/2*sqr(p[1])- p[4],2)));

  
}

void fdf(const size_t n, const double *x,void *fparams,double *fval,double *grad)
{
  f(n, x, fparams, fval);
  df(n, x, fparams, grad);
}


void f2(const size_t n, const double *x,void *fparams,double *fval)
{
  double *p = (double *) fparams;
  // *fval=p[2]*(x[0]-p[0])*(x[0]-p[0]) + p[3]*(x[1]-p[1])*(x[1]-p[1]) + p[4];
  // *fval = g(x[0], x[1]);
  
  double m_W = 80.4;
  //form:
  //nuy = (sqr(m_W)*ly + 2*lx*ly*nux + m_W*lt*sqrt(sqr(m_W) + 4*lx*nux))/2*sqr(lx);


  
  *fval = sqrt(pow( x[0]-p[3] ,2) + pow( (sqr(m_W)*p[2] + 2*p[1]*p[2]*x[0] - m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x[0]))/(2*sqr(p[1])) - p[4],2));
  
}

void df2(const size_t n, const double *x,void *fparams,double *grad)
{
  double *p = (double *) fparams;
  double m_W = 80.4;
  grad[0] = 2*(x[0] - p[3])*x[0] + 2*((sqr(m_W)*p[2] + 2*p[1]*p[2]*x[0] - m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x[0]))/2*sqr(p[1]) - p[4]) * ((p[2]/p[1]) - (m_W*p[0]*p[1]/(sqr(p[1])*sqrt(sqr(m_W) + 4*p[1]*x[0])))) / (2*sqrt( pow(x[0]-p[3],2) + pow((sqr(m_W)*p[2] + 2*p[1]*p[2]*x[0] - m_W*p[0]*sqrt(sqr(m_W) + 4*p[1]*x[0]))/2*sqr(p[1])- p[4],2)));

  
}

void f2df2(const size_t n, const double *x,void *fparams,double *fval,double *grad)
{
  f2(n, x, fparams, fval);
  df2(n, x, fparams, grad);
}

