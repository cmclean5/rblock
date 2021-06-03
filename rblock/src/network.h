//define guards, so headers are declare only once.
#ifndef NETWORK_H
#define NETWORK_H

#include "Headers.h"
#include "vertex.h"

// struct to hold an edge
struct edgelist{
  int source;     // vertex source id
  int target;     // vertex target id  
  double weight;  // edge weight
 } ;

// Function to compare the IDs of two vertices
struct sortIDs{
    bool operator()(const vertex &v1p, const vertex &v2p) const {
      return v1p.id < v2p.id;
    }
  };


struct sortpairIntInt{
  bool operator()(const std::pair<int,int> &l, const std::pair<int,int> &r) const {

    if( l.first < r.first ) return true;
    if( l.first > r.first ) return false;
    return l.second < r.second;    
    //return l.first < r.first;

    }
  };



class network : vertex {

 public:
  network();
  ~network();
  void       removeVertices( int[], int, int );
  int        buildAdj( bool=true );
  void       emptyA ( int );
  void       buildA();
  void       checkA();
  void       buildEdgeList();
  void       countEdges();
  int        countDegree();
  void       freeV();
  void       freeA();
  void       freeModels();
  double *   getA();
  int        getN();
  int        getM();
  int        getM2();
  int        getSelf();
  edgelist*  getEdgeList();
  int        getelSize();
  int        getMinAnnosLength();
  int        getMaxAnnosLength();
  void       setN();
  void       setM();
  void       setM( int );
  void       setSelf( bool );
  void       setVDefaultAnno( string );
  void       findInOutDegrees();

  void       assignNmodels( int );
  void       assignModel( int, int );
  int        getModel( int );
  
  void printA();
  void printEdgelist();
  void printVertices();
  void printAnnoData();
  void printModels();
  int  getMaxK();
  int  getMinK();
  void reorderK();
  void offSetK( int );
  void setPrint( bool );
  
  int nvertices;     // Number of vertices in network
  int directed;      // 1 = directed network, 0 = undirected
  int sortK;         // 1 = we've reorder node's com. numbers, 0 = not
  vertex *V;         // Array of VERTEX structs, one for each vertex


  //Print _Nr x _Nc matrix _M
  static inline void printM( double *_M, int _Nr, int _Nc, const char* _Name ){

    int i,j,k,KK;
    
    cout << "Printing Matrix: " << _Name << endl;

    KK=_Nr * _Nc;
    for(k=0; k<KK; k++){
    i = floor(k/_Nr);
    j = k % _Nr;
    if( _M[(i*_Nr)+j] != 0 ) cout << " " << _M[(i*_Nr)+j] << " ";
    else                     cout << " . ";

    if( j == (_Nr-1) )       cout << "" << endl;
    }
    
 };

 
 private:
  edgelist *el;
  
  int   AKK;//Adjanceny matrix size
  double *A;//Adjacency matrix
  double *Awe;//Weighting Adjacency matrix, note some edges may have zero weightings.
  int     N;
  int     M;
  int     M2;
  int     elSize;
  int     loops;

  int     Nms;
  int     *models; //Metadata model type, 0 = categoric, 1 = continuous
  
  bool PRINT;
  bool self;

  void checkVertexDegree( double * );
  void checkVertexDegree( double *, double * );

 protected:
  typedef std::pair<int,int>    pairIntInt;
  
};

#endif

