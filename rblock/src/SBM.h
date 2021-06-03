//define guards, so headers are declare only once.
#ifndef SBM_H
#define SBM_H

#include "BP.h"

// Program based on Newman's 2014 code for a standard SBM
/* Program to perform the full two-group EM/BP calculation using the
 * standard SBM on an arbitrary network read from a file
 *
 * Written by Mark Newman  21 NOV 2014
 * Modified to include enumerative metadata stored in the "label" field of
 *   the network structure  23 NOV 2014
 */

//class SBM : public Measures {
class SBM : public BP {

 public:
 SBM ();
 SBM ( network *, int );
 SBM ( network *, int, int );
 SBM ( network *, int, int, int );
 ~SBM();

  void useSeed( int );
  
  int    runSBM();
  int    emSBM();
  double mdlSBM();
  
  void getVertexKProbs();
  void getEdgeKProbs();

  double getMaxLL();
  double getEMsteps();
  double getMaxDelta();
 
  
 private:
  void   initSBM();
  void   initOmegaSBM();
  double updateParmSBM();
  
 protected:
  int EM_MAXSTEP;      // Maximum number of EM steps before aborting
  double EM_ACC;       // Required accuracy for EM to terminate

  int    its;
  double LL;           // The Log-Likelihood 

  double KSMALL;
  
  bool NOCONVERGE;     // Set to abort if the solution doesn't converge
  bool calVprob;
  bool calEprob; 

  double STEPS;
  double MAXDELTA;
  
  double saved_LL;
  double saved_steps;
  double saved_maxDelta;
  
  double **saved_q;
  double **qz;    

  void   initialise();  
  void   initCategoric();
  void   initContinuous();

  void   updateCategoric();
  void   updateContinuous();

  void   updateBernco();
  void   updateQbern();

  void   energyCategoric ( int, double & );
  void   energyContinuous( int, double & );

  void   saveqCategoric();
  void   saveqContinuous();
  
  void   setup();
  void   assignSpace();
  void   freeSpace();

  void   checkValue( double );
  
};

#endif
