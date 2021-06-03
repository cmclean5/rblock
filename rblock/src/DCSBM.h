//define guards, so headers are declare only once.
#ifndef DCSBM_H
#define DCSBM_H

#include "SBM.h"

/*
  Based on M.E.J. Newman & A. Clauset. "Structure and inference in annotated networks", nature communitions, 7:11863, doi:1038/ncomms11863.
 */

//class DCSBM : public Measures {
class DCSBM : public SBM {

 public:
 DCSBM ();
 DCSBM ( network *, int );
 DCSBM ( network *, int, int );
 DCSBM ( network *, int, int, int );
 ~DCSBM();

  int    runDCSBM();
  int    emDCSBM();
  double mdlDCSBM();
  
 private: 

  void   initDCSBM();
  void   initOmegaDCSBM();
  double updateParmDCSBM();  
  
};

#endif
