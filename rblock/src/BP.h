//define guards, so headers are declare only once.
#ifndef BP_H
#define BP_H

#include "metaData.h"


class BP : public metaData {

 public:
  BP();
  BP( network * );
  BP( network *, int );
  BP( network *, int, int );  
  ~BP();

 protected:

  //int    doBeliefPropagation();
  int    bpSBM();      // message passing algorithm for standard Block Model 
  int    bpDCSBM();    // message passing algorithm for degree corrected Block Model  

  int BP_MAXSTEP;      // Maximum number of BP steps before aborting
  double BP_ACC;       // Required accuracy for BP to terminate
  double SMALL;
  
  double ***eta;       // Messages            [ NxDegreexK ]
  double **q;          // One-point marginals [ NxK ]

  double **omega;      // Mixing parameters   [ KxK ]

 private:
  void setup();
  void assignSpace();
  void freeSpace();
  
};

#endif
