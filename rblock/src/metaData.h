//define guards, so headers are declare only once.
#ifndef METADATA_H
#define METADATA_H

#include "measures.h"


class metaData : public measures {

 public:
  metaData();
  metaData( network * );
  metaData( network *, int );
  metaData( network *, int, int );
  ~metaData();

  int  getLb();
  void setLb( int );
  void CGRpoint( int, double & );
  void shiftedCGRpoint( int, double & );
  
 protected:

  int model;           // The type of Metadata model 0 = categoric, 1 = continous
  
  int MAXMETA;         // Maximum Metadata strings
  
  //--- meta-data
  int     *x;          // Metadata            [ Nx1 ]
  int     *nx;         // Number of nodes with each distinct metadata value [ 'x'x1 ] 
  double **nrx;        // Expected number in each group with value [ Kxnmlabels ]
  char **mlabel;       // Metadata strings
  int nmlabels;        // Number of distinct metadata strings  
  int SLOT;            // index in networks Nannos to use

  double **gmma;       // Prior parameters (spelled "gmma" rather than "gamma"
                       // [ Kxnmlabels ]

  double *xbern;       // continuous values for x [ Nx1 ] and in range [0,1]
  double **bernco;     // Bernstein coefficients [ KxD ]
  int      D;          // Degree of the Bernstein polynomials
  int     Lb;          // Length scale mapping for Bernstein Polynomials on semifinite Interval 
  double ***Qbern;     // [ NxDxK ]
  double **bvec;       // [ NxD ]
  double **Qsu;        // [ KxN ]
  
 private:
  string empty;

  void setup();
  void assignSpace();
  void modelCategoric();
  void modelContinuous();
  void freeSpace();
  void storeBernsteinPolValues();
  void unitInterval();
  void arbitraryInterval();
  void arbitraryInterval(double, double);
  void semifiniteInterval();
  
};

#endif
