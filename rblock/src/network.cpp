#include "network.h"

network::network() : vertex() {

  nvertices=0;
  directed =0;
  sortK    =0;
  V        =0;
  AKK      =0;
  A        =0;
  Awe      =0;
  N        =0;
  M        =0;
  M2       =0;
  loops    =0;
  Nms      =0;
  models   =0;
  el       =0;
  elSize   =0;
  self     = false;
  PRINT    = false;
}

network::~network(){

  //---free up memory of vertices
  freeV();
   
  //---free up memory of Adjacency matrix
  freeA();

  //---free up memory of models
  freeModels();
  
   //---delete a 1d array in C++
  if( el != 0 ){ delete[] el; }
  
}

void      network::freeV(){

  if( V!=0 && nvertices!=0 ){ delete[] V; }

  nvertices = 0;
  N         = 0;
  
}

void      network::freeA(){

  if( A!=0 && Awe!=0 && AKK != 0 ){ delete[] A; delete[] Awe; }

  AKK = 0;

}

void      network::freeModels(){

  if( models!=0 && Nms!=0 ){ delete[] models; }
  
}

void      network::assignNmodels( int slots ){
  
  int i;

  if( slots >= 0 ){

    freeModels();
    Nms = slots;
  
    models = new int[Nms];
  
    for(i=0; i<Nms; i++) models[i] = -1;

  }
    
}

void network::assignModel( int slot, int MODEL ){
  
  if( (slot >= 0) && (slot < Nms) ){
    models[slot] = MODEL;
  }
    
}

int network::getModel( int slot ){
  
  if( (slot >= 0) && (slot < Nms) ){
    return models[slot];
  }

  return -1;
}

void      network::emptyA( int _AKK ){

  int k,KK;

  N  = nvertices;
  KK = N * N;
  
  if( _AKK >=0 && _AKK == KK ){

    freeA();
    AKK = _AKK;

    A   = new double[AKK];
    Awe = new double[AKK];
    for(k=0; k<AKK; k++){ A[k] = 0.0; Awe[k] = 0.0; }
    
  }
  
}

void      network::buildA(){

  int i,j,I,J;

  N = nvertices;
  
  for(i=0; i<N; ++i){
    for(j=0; j<V[i].degree; j++){
      I = (int)V[i].id;
      J = (int)V[i].E[j].target;

      if( (!self) && (I == J) ){
	;
      } else {	
	A[(I*N)+J]   = 1.0;
	Awe[(I*N)+J] = (double)V[i].E[j].weight;
      }
    }
  }

 

}

void network::checkA(){

  int i,j,k,KK;
  double Eij, Eji;
  double wei, wej;
    
  N  = nvertices;
  
  //--- check A if network is undirected
  if( directed == 0 ){   

      for(i=0; i<N; i++){
	for(j=0; j<N; j++){

	  if( i <= j ){
      
	    Eij = (double)A[(i*N)+j];
	    Eji = (double)A[(j*N)+i];
      
	    wei = (double)Awe[(i*N)+j];
	    wej = (double)Awe[(j*N)+i];
		
	    if( Eij != Eji ){

	      if( Eij == 1 && Eji == 0 ){
		A[(j*N)+i] = Eij; Awe[(j*N)+i] = wei; }
	      if( Eji == 1 && Eij == 0 ){
		A[(i*N)+j] = Eji; Awe[(j*N)+i] = wej; }
	
	    }
	  }

	}
      }
  }


}

void network::buildEdgeList() {
  
  int i,j;
  
  N = nvertices;
  
  vector<pairIntInt> idpair; 

  for(i=0; i<N; ++i){
    for(j=0; j<V[i].degree; j++){

      if( (!self) && (i == j) ){
	;
      } else {
	idpair.push_back( pairIntInt(i,j) );
      }
    }
  }

  if( directed == 0 ){//undirected  
    //---unique list of source-target id pairs
    sort(idpair.begin(),idpair.end(),sortpairIntInt());
    idpair.erase(unique(idpair.begin(),idpair.end()),idpair.end());
  }

  elSize  = idpair.size();
  el = new edgelist[elSize];
    
  for( i=0; i<elSize; i++ ){
    int I = idpair[i].first;
    int J = idpair[i].second;
    el[i].source = (int)V[I].id;
    el[i].target = (int)V[I].E[J].target;
    el[i].weight = (double)V[I].E[J].weight;
  }

}

edgelist* network::getEdgeList() {

  if( el == 0 ) buildEdgeList();
  
  return el;

}

int network::getelSize(){ return elSize; }

void network::countEdges(){

  int i,j,k,KK;

  N = nvertices;

  KK = N * N ;
  
  //---test, counting number of edges in A
  for(k=0, M2=0, M=0, loops=0; k<KK; k++){
    if( A[k] == 1 ) M++;
    
    i = floor(k/N);
    j = k % N;

    if( j>=i && A[(i*N)+j] == 1 ){
      M2++;
      if(j==i){ loops++; }
    }

  }
       
  if( PRINT ){ cout << "TEST: M " << M << " M/2 " << M2 << " loops " << loops << " sum(D) " << countDegree() << endl; }
  
}

int network::countDegree(){

  int i,sum;

  N = nvertices;
  
  for(i=0,sum=0; i<N; ++i){ sum += (int)V[i].degree; }

  return sum;
  
}


int       network::getN()        { return nvertices; }

int       network::getM()        { return M; }

int       network::getM2()       {

  if( M2 == 0 ){ M2 = (directed == 0) ? (int)M/2 : M; } 

  return M2;

}

int       network::getSelf()     { return self; }

void      network::setN()        { if( V!=0 ) N = nvertices; }

void      network::setM( int Mm ){ if( V!=0 ) M = Mm; }

void      network::setSelf( bool useself ){ self = useself; }

void      network::setM(){

  int i;
  
  if( V!=0 ){
    
    N = nvertices;
    M = 0;  
    M = countDegree();

  }
    
}


double*  network::getA()        { return A; }


//build the network repersentations, Adjacency Matrix & edge list
int network::buildAdj( bool Check ){

  int i,j,m,k,I,J,KK;
  int status;
  double sum;
  
  if( V!= 0 ){

    N = nvertices;
    
    //linear indexing, size K is rows (N) x cols (N).
    KK=N*N;

    //---
    //---Create Adjacency matrix
    emptyA  ( KK );
    buildA ();

    //--- check A if network is undirected
    checkA();
    
    
    //---test, counting number of edges in A
    countEdges();

    //Check vertex degree from .gml file and our Adjacency matrix,
    //Always set to true if we're *not* reading in networks from .gml file. 
    if( Check == true ){

      checkVertexDegree( A, Awe );    

      //---test, counting number of edges in A
      countEdges();
    
    }
               

  } else { return -1; }

  return 0;

}

void network::checkVertexDegree( double *AA ){

  int i,j,k,K,sum;

  N = nvertices;
  
  if( PRINT ){cout << "checking vertex degrees..."; }
   
  for(i=0; i<N; ++i){
    sum=0;
    k=0;
    for(j=0; j<N; ++j){ if(AA[(i*N)+j] == 1) sum++; }
    V[i].assignE( sum );
    for(j=0; j<N; ++j){ 
      if(AA[(i*N)+j] == 1){
	V[i].E[k].target = j; 
	k++; 
      }//if
    }
  }

  if( PRINT ){cout << "done." << endl; }

 }


void network::checkVertexDegree( double *AA, double *AAwe ){

  int i,j,k,K,sum;

  N = nvertices;
  
  if( PRINT ){cout << "checking vertex degrees..."; }
   
  for(i=0; i<N; ++i){
    sum=0;
    k=0;
    for(j=0; j<N; ++j){ if(AA[(i*N)+j] == 1) sum++; }
    V[i].assignE( sum );
    for(j=0; j<N; ++j){ 
      if(AA[(i*N)+j] == 1){
	V[i].E[k].target = j; 
	V[i].E[k].weight = (double)AAwe[(i*N)+j];
	k++; 
      }//if
    }
  }

  if( PRINT ){cout << "done." << endl; }

 }

void network::removeVertices( int keys[], int length, int dummy ){

  int i,j,p,k,KK,KKin,Ng;

  N = nvertices;

  //find size of remaining vertices
  Ng = 0;
  for(i=0; i<N; i++){
    if( keys[i] !=  dummy ) Ng++;
  }
  
  if( (length == N) && (Ng != 0) && (Ng < N) ){

    if( PRINT ){cout << "Resizing V from " << N << " to " << Ng << "..." << endl; }
    
    //Set current size of the Adjaceny matrix
    KKin = N * N;

    //make sure A exists
    if( A==0 || AKK == 0 ){
      emptyA ( KKin );
      buildA();
    }

    //construct temp A
    KK = Ng * Ng;
    double tA[KK];
    
    for(k=0, p=0; k<KKin; k++){

      i = floor(k/N);
      j = k % N;

      if( keys[i] != dummy && keys[j] != dummy ){
	tA[p] = 0.0;
	tA[p] = A[(i*N)+j]; p++; }
      
    }
   
    //free A
    freeA();

    //copy V to temp V
    vertex *tV = new vertex[Ng];
    for(i=0, p=0; i<N; ++i){

      if( keys[i] != dummy ){
	tV[p].copy(&V[i]); p++;
      }

    }

    //free V
    freeV();

    //resize V to tV size Ng
    V = new vertex[Ng];
    nvertices = Ng;
    for(i=0; i<Ng; ++i){
      V[i].copy(&tV[i]);
      V[i].id = i;
    }

    //delete temp V
    if( tV!=0 && Ng!=0 ){ delete[] tV; Ng=0; }

    //Now assign edges to V.E using the Adjaceny matrix
    checkVertexDegree( tA  );    
    
    if( PRINT ){cout << "done." << endl; }

  }
  
  
}

void network::findInOutDegrees(){

  int i,j,u,k,KK;
  vector<int> inDeg;
  vector<int> outDeg;
      
  if( V!= 0 ){

    N = nvertices;
    
    //linear indexing, size K is rows (N) x cols (N).
    KK=N*N;

    //---
    //---Create Adjacency matrix
    emptyA  ( KK );
    buildA ();

    //--- check A if network is undirected
    checkA();

    //--- find each nodes in and out degrees
    for(i=0; i<N; i++){
      for(j=0; j<N; j++){
	if( A[(i*N)+j] == 1 ) outDeg.push_back(j);	  
	if( A[(j*N)+i] == 1 ) inDeg.push_back(i);	      
      }

      V[i].assignIndeg ( inDeg  );
      V[i].assignOutdeg( outDeg );
      inDeg.clear();
      outDeg.clear();      
    }

    freeA();
    
  }

}

void network::printA(){
  
  int i,j,k,K;

  if( A!=0 && AKK!=0 ){
  
    if( PRINT ){ cout << "Adjancency Matrix: " << endl; }
  N=nvertices;
  K=N*N;
  for(k=0; k<K; k++){
    i = floor(k/N);
    j = k % N;
    if( A[(i*N)+j] != 0 ) if( PRINT ){ cout << " " << A[(i*N)+j] << " "; }
    else                  if( PRINT ){ cout << " . "; }

    if( j == (N-1) )      if( PRINT ){ cout << "" << endl; }
  }

  if( PRINT ){ cout << "" << endl; }

  } else {
    if( PRINT ){ cout << "Empty container." << endl; }
  }
  
  
}

void network::printEdgelist(){
  
  int i;

  if( el!=0 && elSize!=0 ){
  
    if( PRINT ){ cout << "Edge list: " << endl; }

    for( i=0; i<elSize; i++ ){
      if( PRINT ){
	cout << i << ": " << el[i].source << "\t" << el[i].target << "\t" << el[i].weight << endl;  
      }
    }
  } else {
    if( PRINT ){ cout << "Empty container." << endl; }
  }
  
  
}


void network::printVertices(){

  
  int v;

  if( V!=0 && nvertices!=0 ){

    N    = nvertices;
    
    if( PRINT ){ cout << "Vertex Properties" << endl; }

    for(v=0; v<N; ++v){
      if( PRINT ){ cout << "id " << V[v].id << " label " << V[v].label
			<< " degree " << V[v].degree
			<< " in degree " << V[v].inDeg
			<< " out degree " << V[v].outDeg
			<< " K " << V[v].K << endl; }
    }

  } else {
    if( PRINT ){ cout << "Empty container." << endl; }
  }
  

}

void network::printAnnoData(){

  int v,j;

  if( V!=0 && nvertices!=0 ){

    N    = nvertices;
    
    if( PRINT ){ cout << "Vertex Annotation Data" << endl; }

  for(v=0; v<N; ++v){
    if( PRINT ){ cout << "id: " << V[v].id << " label: " << V[v].label << " = ["; }
    for(j=0; j<V[v].Na; j++){
      if( j == (V[v].Na-1) ){ if( PRINT ){ cout << V[v].Annos[j]; }
      } else                  if( PRINT ){ cout << V[v].Annos[j] << ","; }
    }
    if( PRINT ){ cout << "]" << endl; }
  }

  } else {
    if( PRINT ){ cout << "Empty container." << endl; }
  }
  
}

void network::printModels(){

  int i;

  if( models!=0 && Nms!=0 ){

    if( PRINT ){ cout << "Metadata model type(s)" << endl; }

    for(i=0; i<Nms; i++){
      if( PRINT ){
	if( models[i] == 0 ){
	  cout << "[" << i << "] = "
	       << "Metadata model type set to Categoric" << endl;
	}
	if( models[i] == 1 ){
	  cout << "[" << i << "] = "
	       << "Metadata model type set to Continuous" << endl;
	}	
      }
    }
  } else {
    if( PRINT ){ cout << "Empty container." << endl; }
  }
  
}


void network::setVDefaultAnno( string ANNO ){

  int v,j;

  if( V!=0 && nvertices!=0 ){

    N    = nvertices;

    for(v=0; v<N; v++){
      V[v].assignNannos( 1 );
      V[v].Annos[0] = ANNO;
    }

  }

}

int network::getMaxAnnosLength(){

  int v, max;

  if( V!=0 && nvertices!=0 ){

    N    = nvertices;
    max  = 0;
    
    for(v=0; v<N; v++){
      if( V[v].Na > max ) max = V[v].Na;
    }

  }

  return max;
  
}

int network::getMinAnnosLength(){

  int v, min;

  if( V!=0 && nvertices!=0 ){

    N    = nvertices;
    min  = getMaxAnnosLength();
    
    for(v=0; v<N; v++){
      if( V[v].Na < min ) min = V[v].Na;
    }

  }

  return min;
  
}

int network::getMaxK(){

  int v, Ncom;

  Ncom = -1;
  
  if( sortK == 0 ) reorderK();
  
  if( V!=0 && nvertices!=0 ){

    N    = nvertices;
    Ncom = 0;
    Ncom = V[0].K;
  
    for(v=0; v<N; ++v){
      if( V[v].K > Ncom ) Ncom = V[v].K; 
    }
  
  }  

  return Ncom;
  
}

int network::getMinK(){

  int v, Ncom;

  Ncom = -1;

  if( sortK == 0 ) reorderK();  
  
  if( V!=0 && nvertices!=0 ){

    N    = nvertices;
    Ncom = 0;
    Ncom = V[0].K;
  
    for(v=0; v<N; ++v){
      if( V[v].K < Ncom ) Ncom = V[v].K; 
    }
  
  }  

  return Ncom;
  
}

void network::reorderK(){

  int v,counter,Knew,Kmax;

  if( V!=0 && nvertices!=0 ){

    N = nvertices;
    
    int temp[N];
    counter = getMinK();
    Knew    = 1;
    Kmax    = getMaxK();

    while( counter <= Kmax ){

      bool found=false;
      
      for(v=0; v<N; v++){
	if( V[v].K == counter ){
	  temp[v] = Knew;
	  found=true;
	}
      }

      if(found) Knew++;
      
      counter++;
    }

    for(v=0; v<N; v++){ V[v].K = temp[v]; }
    
  }

  //save that we've reorder node community numbers
  sortK = 1;
  
}

void network::offSetK( int offset ){

  int v;

  if( V!=0 && nvertices!=0 && offset >= 0 ){

    N    = nvertices;
  
    for(v=0; v<N; ++v)
      V[v].K = (V[v].K - offset) + 1;

  }
    
}

void network::setPrint( bool status ){
  PRINT = status;
}
