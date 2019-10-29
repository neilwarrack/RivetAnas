#include "MendelMin.h"
using namespace Rivet;
using namespace std;

double target(const MendelMin::Params& p) {
  const double x = p[0] - 0.3;
  const double y = p[1] - 0.5;
  return (x*x + 2*y*y) + 0.1*sin(x)*sin(y);
}

double rand01() {
  static random_device rd;
  static mt19937 gen(rd());
  static uniform_real_distribution<> dis(0.0, 1.0);
  return dis(gen);
}


int main() {

  MendelMin mm(target, rand01, 2);
  const double best = mm.evolve(20);

  valarray<double> fittest = mm.fittest();
  cout << "FOUND ANSWER: " << target(fittest) << " == " << best << "; ";
  for (double x : fittest)
    cout << x << " ";
  cout << endl;

  valarray<double> right{0.3, 0.5};
  cout << "RIGHT ANSWER: " << target(right) << "; ";
  for (double x : right)
    cout << x << " ";
  cout << endl;

  return 0;
}
