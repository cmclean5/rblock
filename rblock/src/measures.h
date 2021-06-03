//define guards, so headers are declare only once.
#ifndef MEASURES_H
#define MEASURES_H

#include "Headers.h"
#include "network.h"

//Bernstein polynomial class
//#include "bernstein_polynomial.h"


class measures {

 public:
  measures();
  measures(network*);
  measures(network*, int);
 ~measures();

  void setOpenMP( int );
  int  getnCORES();
  
 int  getTwom();
 void bridgeness();
 void partitionDensity(double &);
 
 protected:
  network   *gg;//network vertices, see network.h for properties
  edgelist  *el;//list of edges
  double   **qz;//The edge probabilities, for the BlockModel

  int nCORES;
  int MAXEXPO;
  int N;        //Number of Vertices
  int M;        //Number of Edges
  int twom;     //Two times the number of edges 
  int K;        //Number of communities, for the BlockModel
  int El;       //Size of the edge-list - this may not be the same as M, since duplicate edges are removed in undirected networks.

  //GSL random number and seed
  unsigned long int seed;
  unsigned long int seedOffset;
  gsl_rng *g;
  void               testSeed(int=0);
  void               setSeed(bool=false, int=0);
  unsigned long int  getSeed();
  void               random_unity(int, double*);
  void               random_vec  (int, double*);
  double             entropy      (double);
  double             binaryEntropy(double);
  double             loggamma     (int);
  double             multico      (int, int);
  double             binco        (int, int);
  //void bernstein_polynomial(int, int, double, double &);
  void               min(int, double*, double &);
  void               max(int, double*, double &);
  
 private:
  void setup();
  void assignSpace();
  void freeSpace();
  
};

#endif
