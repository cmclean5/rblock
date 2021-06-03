//define guards, so headers are declare only once.
#ifndef VERTEX_H
#define VERTEX_H

#include "edge.h"

class vertex : edge {

 public:
  vertex();
  ~vertex();
  void copy( vertex* );
  void assignNannos  ( int );
  void assignAnno    ( int, string );
  void assignKprobs  ( int );
  void assignE       ( int );
  void assignIndeg   ( int );
  void assignIndeg   ( vector<int> &);
  void assignOutdeg  ( int );
  void assignOutdeg  ( vector<int> &);
  void freeNannos    ();
  void freeKprobs    ();  
  void freeE         ();
  void freeIndeg     ();
  void freeOutdeg    ();
  
  void freeSpace();
  void printV();
  void setPrint( bool );
  
  int id;            // GML ID number of vertex
  int degree;        // Degree of vertex (out-degree for directed nets)
  int inDeg;        
  int outDeg;
  int K;             // Concrete community;
  double bridge;     // Vertex bridgeness
  char *label;       // GML label of vertex.  NULL if no label specified
  string GeneName;   // GeneName
  string EntrezID;   // Gene Entrez ID
  edge   *E;         // Array of EDGE structs, one for each neighbor
  double *Kprobs;    // Community (K) probabilities
  string *Annos;     // N annotation types

  int *Indeg;
  int *Outdeg;
  
  int Nk;            // Size of Kprobs
  int Na;            // Size of Nannos;
 
 private:  
  bool PRINT;
  
};

#endif
