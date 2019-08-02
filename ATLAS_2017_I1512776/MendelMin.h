#ifndef MendelMin_H
#define MendelMin_H

#include <valarray>
#include <random>
#include <functional>
#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace Rivet {

/// MendelMin implements a home brewed genetic algorithm for finding
/// the minimum of a function defined on a unit hypercube returning a
/// mon-negative real number (eg. a Chi-squared value).
// template <typename FuncT, typename RndT = (*)()>
class MendelMin {
public:

  /// Typedef for a valaray of parameters to the function to be minimised.
  using Params = std::valarray<double>;
  /// Typedef for the function to be minimised
  using FuncT = std::function<double(const Params&, const Params&)>;
  /// Typedef for the [0,1] random number generator
  using RndT = std::function<double()>;


  /// Constructor taking as arguments: the function, @a fin, to be
  /// minimised; a random number generator, @rndin, returning double
  /// random numbers between 0 and 1; the dimension, @a ndim, of the
  /// unit hypercube for which @a fin is defined. Optional arguments
  /// are: the number, @a npop, of individuals in the population; and
  /// @a margin which determines how much ramndomness is involved when
  /// an individual is evolved twowards the fittest individual.
 MendelMin(const FuncT & fin, const Params & fixpar, const RndT & rndin, unsigned int ndim,
            unsigned int npop = 20, unsigned int ngen = 20,
            double margin = 0.1)
   : _f(fin), _q(fixpar), _rnd(rndin), _NDim(ndim), _margin(margin),
      _pop(npop), _fit(npop, -1.0), showTrace(false) {}


  /// Supply a best guess for the fittest parameter point to help
  /// things along.
  void guess(const Params & p) {
    _pop.push_back(p);
    limit01(_pop.back());
    _fit.push_back(-1.0);
  }


  /// Evolve the population a given number of generations and return
  /// the best fit value.
  double evolve(unsigned int NGen) {
    for ( unsigned n = 0; n < NGen; ++n ) {
      // Calculate the fitness.
      auto mm = minmax();
      // Always kill the fittest individual.
      if ( showTrace ) _debug();
      for ( unsigned int i = 1; i < _pop.size(); ++i ) {
        if ( _fit[i] > rnd()*(mm.second - mm.first) )
          // Kill all individuals that have low fitness or are just unlucky.
          _fit[i] = -1.0;
        else
          // Improve This individual to be more like the fittest.
          move(_pop[i],_pop[0]);
      }
    }
    return _fit[0];
  }

  /// Return the fittest parameter point found.
  Params fittest() const {
    return _pop[0];
  }

  /// Return the fittest value found.
  double fit() const {
    return _fit[0];
  }

  /// Simple wrapper around the random number generator.
  double rnd() const {
    return _rnd();
  }

  /// Return a random parameter point in the unit hypercube.
  Params rndParams() const {
    Params ret(_NDim);
    for ( unsigned int i = 0; i < _NDim; ++i ) ret[i] = rnd();
    return ret;
  }

  /// Limit a parameter point to inside the unit hypercube.
  void limit01(Params & p) const {
    for ( unsigned int i = 0; i < _NDim; ++i )
      p[i] = std::max(0.0, std::min(p[i], 1.0));
  }

  /// Move an @a bad parameter point towards a @a better one. The new
  /// point is picked randomly within the generalized hypercube where
  /// @bad and @better are at diagonally opposite corners, enlarged by
  /// a fraction _margin.
  void move(Params & bad, const Params & better) const {
      bad += (better - bad)*(rndParams()*(1.0 + 2.0*_margin) - _margin);
    limit01(bad);
  }

  /// Simple wrapper around the function to be minimised.
  double f(const Params & p, const Params & q) const {
    return _f(p, q);
  }


  /// Calculate the fitness values of all individuals and put the
  /// fittest one first. @return the best and worst fitness values.
  std::pair<double, double> minmax() {
    std::pair<double,double> mm(std::numeric_limits<double>::max(), 0.0);
    unsigned int iwin = 0;
    for ( unsigned int i = 0; i < _pop.size(); ++i ) {
      double & v = _fit[i];
      // negative fitness value means the individual is dead, so we
      // welocme a new immigrant.
      if ( v < 0.0 ) _pop[i] = rndParams();

      // The calculated fitness value cannot be negative.
      v = std::max(0.0, f(_pop[i],_q));

      // Compare to the best and worst fitness so far.
      if ( v < mm.first ) iwin = i;
      mm.first = std::min(v, mm.first);
      mm.second = std::max(mm.second, v);
    }

    // Move the winner to the top.
    if ( iwin != 0 ) {
      std::swap(_pop[0], _pop[iwin]);
      std::swap(_fit[0], _fit[iwin]);
    }
    return mm;
  }

  /// Inspect a population and the fitness of its individuals.
  void _debug() {
    std::cout << "GenAlgMax population status:" << std::endl;
    for ( unsigned int i = 0; i < _pop.size(); ++i ) {
      std::cout << std::setw(10) << _fit[i] << " (" << _pop[i][0];
      for ( unsigned int ip = 1; ip < _NDim; ++ip )
        std::cout << "," << _pop[i][ip];
      std::cout << ")" << std::endl;
    }
  }

private:

  /// The function to be minimised.
  const FuncT _f;

    /// The fixed parameter points.
  Params _q;
  
  /// The fixed parameters of the function.
  //  const double _q;

  /// The random number generator.
  const RndT _rnd;

  /// The dimension of the unit hypercube of allowed parameter points.
  unsigned int _NDim;

  /// WHen evolving less fit parameter point towards the fittest one
  /// we choose randomly in a generalised hypercube where these two
  /// parameter points are in diagonally opposite corners, expanded in
  /// each direction with This fraction.
  double _margin;

  /// The population of parameter points.
  std::vector<Params> _pop;  

  
  /// The fitness value of each individual in the population.
  std::vector<double> _fit;

public:

  /// Set true to get a verbose record of the evolution.
  bool showTrace;

};


/// Construct a MendelMin object taking as arguments: the function, @a
/// fin, to be minimised; a random number generator, @rndin, returning
/// double random numbers between 0 and 1; the dimension, @a ndim, of
/// the unit hypercube for which @a fin is defined; the number, @a
/// npop, of individuals in the population; and a @a margin which
/// determines how much ramndomness is involved when an individual is
/// evolved twowards the fittest individual.
// template <typename FuncT, typename RndT>
// MendelMin <FuncT, RndT>
// makeMendelMin(const FuncT & f, const RndT & rnd, unsigned int ndim,
//               unsigned int npop = 20, double margin = 0.1) {
//   return MendelMin<FuncT, RndT>(f, rnd, ndim, npop, margin);
// }

}

#endif // MendelMin_H
