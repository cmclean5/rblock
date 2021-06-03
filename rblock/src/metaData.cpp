#include "metaData.h"

metaData::metaData() : measures(){ setup(); }

metaData::metaData( network *gg ) : measures(gg) {
  setup();
  gg->setVDefaultAnno( empty );
  assignSpace();
}

metaData::metaData( network *gg, int K ) : measures(gg, K) {
  setup();
  gg->setVDefaultAnno( empty );
  assignSpace();
}

metaData::metaData( network *gg, int K, int annoSLOT ) : measures(gg, K) {
  setup();
  this->SLOT=annoSLOT;
  this->model=gg->getModel(SLOT);
  assignSpace();
}

void metaData::setup(){

  this->model=0;
  this->MAXMETA=300;
  this->SLOT=0;
  this->empty="";
  
  this->x=0;
  this->nx=0;
  this->nrx=0;
  this->mlabel=0;
  this->nmlabels=0;
  this->gmma=0;
  this->xbern=0;
  this->bernco=0;
  this->D=100;
  this->Lb=1;
  this->Qbern=0;
  this->bvec=0;
  this->Qsu=0;
}

void metaData::assignSpace(){

  if( model == 0 ) modelCategoric();
  else             modelContinuous();
 
}

void metaData::modelCategoric(){

  int u,i,r;

  //--- Make space for the metadata numbers and labels
  x = (int*)malloc(N*sizeof(int));
  mlabel = (char**)malloc(MAXMETA*sizeof(char*));
  for( i=0; i<MAXMETA; i++ ) mlabel[i] = (char*)malloc(MAXMETA*sizeof(char*));
  nmlabels = 0;

  //--- Go through the vertices
  for (u=0; u<N; u++) {

    //--- Check to see if this label is already in the list of labels
    for (i=0; i<nmlabels; i++) {
      if (strcmp(gg->V[u].Annos[SLOT].c_str(),mlabel[i])==0) break;
    }

    //--- If not, add it
    if (i==nmlabels) {
      strcpy(mlabel[nmlabels++],gg->V[u].Annos[SLOT].c_str());
    }

    //--- Record this as the metadata type for this vertex
    x[u] = i;

  }

  //--- Count how many nodes there are in each metadata group
  nx = (int*)calloc(nmlabels,sizeof(int));
  for (u=0; u<N; u++) nx[x[u]]++;

  cout << "Found " << nmlabels << " distinct metadata values:" << endl;
  for (i=0; i<nmlabels; i++) cout << " " << (i+1) << " " << mlabel[i] << endl;

  nrx = (double**)malloc(K*sizeof(double));
  for (r=0; r<K; r++) nrx[r] = (double*)malloc(nmlabels*sizeof(double));

  // Malloc space for parameters gmma and choose random initial values
  // gamma_sx, Prior probability for node assigned to 'K' communities, 
  // using annotation 'x' data. 
  //double ru[K];
  gmma = (double**)malloc(K*sizeof(double));
  for (r=0; r<K; r++) gmma[r] = (double*)malloc(nmlabels*sizeof(double));
  

}

void metaData::modelContinuous(){

  int u,i,r;
  
  // Make space for the Bernstein coefficients
  bernco = (double**)malloc(K*sizeof(double));
  for (r=0; r<K; r++) bernco[r] = (double*)malloc((D+1)*sizeof(double));  
  
  // Make space for the continuous data values
  // We'll need to read this in from SLOT
  double xx[N];
  xbern  = (double*)malloc(N*sizeof(double));

  //Store continuous values for this vertex at slot SLOT
  for(u=0; u<N; u++){
    xbern[u] = atof(gg->V[u].Annos[SLOT].c_str());
  }
    
  // Make space for Qbern [ NxDxK ] 
  Qbern = (double***)malloc(N*sizeof(double));
  for (u=0; u<N; u++) {
    Qbern[u] = (double**)malloc((D+1)*sizeof(double));
    for (i=0; i<=D; i++) {
      Qbern[u][i] = (double*)malloc(K*sizeof(double));
    }
  }

  // Make space for Qsu, Qbern summed over Bernstein's Degree [KxN], i.e.
  // Prior probability for node assigned to 'K' communities, 
  // using annotation 'x' data. 
  Qsu   = (double**)malloc(K*sizeof(double));
  for(r=0; r<K; r++){
    Qsu[r]   = (double*)malloc(N*sizeof(double));
  }
  
  // Make space for bernstein polynomial values [ NxD ] 
  bvec = (double**)malloc(N*sizeof(double));
  for (u=0; u<N; u++) {
    bvec[u] = (double*)malloc((D+1)*sizeof(double));
  }

  // calculate and store in bvec, the Bernstein Polynomial values for each x read in and each degree.
  storeBernsteinPolValues(); 
  
}

void metaData::freeSpace(){

  int i,u;
  
  if( x!= 0 )    { free(x); }

  if( nx!= 0 )   { free(nx); }

  if( nrx!=0 )   {for (i=0; i<K; i++) { free(nrx[i]);   } }

  if( gmma!=0 )  {for (i=0; i<K; i++) { free(gmma[i]);  } }

  if( xbern!=0 ) { free(xbern); }
  
  if( bernco!=0 ){for (i=0; i<K; i++) { free(bernco[i]);  } }

  if( bvec!=0 )  {for (i=0; i<N; i++) { free(bvec[i]); } }
  
  if( Qbern!=0 ){
    for (u=0; u<N; u++) {
      for (i=0; i<D; i++) free(Qbern[u][i]);
      free(Qbern[u]);
    }
    free(Qbern);
  }

  if( Qsu!=0 )   {for(i=0; i<K; i++){ free(Qsu[i]); } free(Qsu); }
  
  if( mlabel!=0 ){for (i=0; i<nmlabels; i++) { free(mlabel[i]); } }
  
}

metaData::~metaData(){ freeSpace(); }

//Check the interval of the continuous values, then
//calculate and store Bernstein Polynomial values in bvec.
void metaData::storeBernsteinPolValues(){

  double a,b;

  min(N,xbern,a);
  max(N,xbern,b);

  if( (a>=0) && (b<=1) ) unitInterval();
  else                   arbitraryInterval(a,b);
  
}

// calculate and store in bvec, the Bernstein Polynomial values for each x read in and degree, for x values in a unit interval, i.e. [0,1].
void metaData::unitInterval(){

  int u,j;
  double bico;
  
  for(j=0; j<=D; j++){
    bico = 0.0;
    bico = binco( D, j );
    for(u=0; u<N; u++){ 
      bvec[u][j] = bico * pow( xbern[u], j ) * pow( (1-xbern[u]), (D-j) );
    }
  }
  
  
}

// calculate and store in bvec, the Bernstein Polynomial values for each x read in and degree, for x values in an arbitrary interval, i.e. [a,b].
void metaData::arbitraryInterval(double a, double b){
  //Ref: https://www.sciencedirect.com/science/article/pii/S0377042715002459
  //This is effectively the max-min scaling, i.e. a linear transformation,
  //replacing x => (x-a)/(b-a) and 1-x => (b-x)/(b-a) 
  
  int u,j;
  double bico, num;

  for(j=0; j<=D; j++){
    bico = 0.0;
    bico = binco( D, j );
    for(u=0; u<N; u++){
      num  = 0.0;
      num  = pow( (xbern[u]-a), j) * pow( (b-xbern[u]),(D-j) );
      num /= pow( (b-a), D );
      bvec[u][j] = bico * num;
    }
  }
  
  
}

void metaData::arbitraryInterval(){
  //Ref: https://www.sciencedirect.com/science/article/pii/S0377042715002459
  
  int u,j;
  double bico, num, a, b;

  min(N,xbern,a);
  max(N,xbern,b);  
  
  for(j=0; j<=D; j++){
    bico = 0.0;
    bico = binco( D, j );
    for(u=0; u<N; u++){
      num  = 0.0;
      num  = pow( (xbern[u]-a), j) * pow( (b-xbern[u]),(D-j) );
      num /= pow( (b-a), D );
      bvec[u][j] = bico * num;
    }
  }
  
  
}

int metaData::getLb(){ return Lb; }

void metaData::setLb(int newLb ){ if( newLb > 0 ) Lb = newLb; } 

void metaData::CGRpoint( int j, double &xj ){
  //Ref: http://pu.edu.pk/images/journal/maths/PDF/Paper-7_49_1_17.pdf
  //Calculate the Chebyshev-Gauss_Radau point for D degree at base j  
  
  double num;

  if( j>=0 && j <= D ){  
    num  = 0;
    num  = 2*j*M_PI;
    num /= (2*D + 1);  
    xj = -1.0 * cos( num );
  }
  
}

void metaData::shiftedCGRpoint( int j, double &xj ){
  //Ref: http://pu.edu.pk/images/journal/maths/PDF/Paper-7_49_1_17.pdf
  //Calculate the shifted Chebyshev-Gauss_Radau point for D degree at base j  

  double num;
  
  if( j>=0 && j <= D ){
    xj = 0;
    CGRpoint(j,xj);  
    num  = 0;
    num  = (Lb*xj + Lb);
    num /= 2*(1-xj);  
    xj   = num;
  }
  
}


void metaData::semifiniteInterval(){
  //Ref: http://pu.edu.pk/images/journal/maths/PDF/Paper-7_49_1_17.pdf
  //See equations 3.11 and 3.12...
  //If x is defined on a semiinfinte interval, i.e. [0,inf), then
  //we can map back to a unit interval (i.e. [0,1]) via x' = x/(x+Lb)
 
  
  int u,j;
  double bico, num;
  
  for(j=0; j<=D; j++){
    bico = 0.0;
    bico = binco( D, j );
    for(u=0; u<N; u++){
      num  = 0.0;
      num  = pow( xbern[u], j) * pow( Lb,(D-j) );
      num /= pow( (xbern[u]+Lb), D );
      bvec[u][j] = bico * num;
    }
  }
  
  
}
