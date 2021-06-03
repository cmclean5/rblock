#include "measures.h"

measures::measures(){ setup(); }

measures::measures(network *gg){

  setup();
  this->gg   = gg;
  this->N    = gg->getN();
  this->M    = gg->getM();
  this->el   = gg->getEdgeList();
  this->El   = gg->getelSize();

  for( int i=0; i<N; i++ ) this->twom += gg->V[i].degree;
  
}

measures::measures(network *gg, int K){

  setup();
  this->gg   = gg;
  this->N    = gg->getN();
  this->M    = gg->getM();
  this->K    = (K > 0) ? K : 1;
  this->el   = gg->getEdgeList();
  this->El   = gg->getelSize();

  for( int i=0; i<N; i++ ) this->twom += gg->V[i].degree;
  
}

void measures::setup(){

  this->gg   =0;
  this->N    =0;
  this->M    =0;
  this->K    =1;

  this->el   =0;
  this->El   =0;

  this->twom =0;
  
  this->qz   =0;
  this->seed =0;
  this->g=0;

  this->MAXEXPO=700;
  this->nCORES=0;
  
  assignSpace();
  
}

void measures::setOpenMP( int ncores ){

  int max_threads = 0;
  int num_procs   = 0;
  max_threads     = omp_get_max_threads();
  num_procs       = omp_get_num_procs();  

  if( ncores > max_threads ){
    ncores = max_threads;
    omp_set_num_threads(ncores);
  }

  if( ncores <= 0 ){
    ncores = 1;
    omp_set_num_threads(1);
  }

  if( ncores > 0 && ncores <= max_threads ){
    omp_set_num_threads(ncores);
  }
  
  cout << "> OpenMP: cores set too: " << ncores << endl;


  nCORES = ncores;
  
}


int measures::getnCORES(){ return(nCORES); }

void measures::assignSpace(){

  //---set gsl random number generator, using taus2 generator
  //g = gsl_rng_alloc(gsl_rng_taus2);
  g = gsl_rng_alloc(gsl_rng_mt19937);  
  
}

void measures::freeSpace(){

  if( g!=0 ){ gsl_rng_free(g); }

}

measures::~measures(){ freeSpace(); }


//---Initialize random seed for testing:
void measures::testSeed (int newSeed ){

  seed = newSeed;
  gsl_rng_set(g, seed);
  
}

//---Initialize random seed:
void measures::setSeed (bool useSeed, int newSeed ){

  unsigned long int x,y,z,c,t;

  c = (unsigned long int) 6543217;

  x = (unsigned long int) time(NULL);

  seedOffset != 0 ? y = (unsigned long int) seedOffset : y = (unsigned long int) 987654321;

  z = (unsigned long int) getpid();

  cout << "" << endl;
  cout << "----" << endl;
  cout << "setting random number seed:" << endl;
  cout << "> x = " << x << endl;
  cout << "> y = " << y << endl;
  cout << "> z = " << z << endl;
  cout << "----" << endl;

  //--- JKISS RGN
  x  = 314527869 * x + 1234567;
  y ^= y << 5; y ^= y >> 7; y ^= y << 22;
  t  = 4294584393 * z + c; c = t >> 32; z = t; 

  seed = (unsigned long int) (x + y + z);
  
  cout << "> setSeed: " << seed << endl;
  cout << "----" << endl;

  gsl_rng_set(g , seed);
  
}



unsigned long int measures::getSeed(){ return seed; }


//---Function to generate d numbers at random that add up to unity
void measures::random_unity(int d, double *xx){
  int k;
  double sum=0.0;

  for (k=0; k<d; k++) {
    xx[k] = 0.0;
    xx[k] = gsl_rng_uniform(g);
    sum += xx[k];
  }
  for (k=0; k<d; k++) xx[k] /= sum;
}

//---Function to generate d random numbers
void measures::random_vec(int d, double *xx){
  int k;

  for (k=0; k<d; k++) {
    xx[k] = 0.0;
    xx[k] = gsl_rng_uniform(g);
  }

}

int measures::getTwom(){

  return twom;
  
}


void measures::bridgeness(){

  int i,k;

  for(i=0; i<N; i++){

    double b=0.0;

    for(k=0; k<K; k++){                                                       
      b  += (gg->V[i].Kprobs[k] - 1.0/(double)K) * (gg->V[i].Kprobs[k] - 1.0/(double)K);
}                     
                         
    gg->V[i].bridge  = 1.0 - sqrt((double)K/((double)K - 1.0) *  b);   
                        
  }

}

/*
Y. Ahn, J. Bagrow, S. Lehmann. Link communities reveal multiscale complexity in networks, Nature Letters, 15 August 2010.
*/
void measures::partitionDensity(double &PD){

  int v,k,m;
  
  double mc = 0.0;
  double nc = 0.0;
  PD        = 0.0;    
  for(k=0; k<K; k++){
      
    mc = 0.0;
    nc = 0.0;
    
    for(v=0; v<N; v++) nc += (double)gg->V[v].Kprobs[k];

    for(m=0; m<M; m++) mc += qz[m][k];      
    
    PD += ( mc - (nc - 1)/nc*(nc-1)/2-(nc-1) );
    
  }

}

// entropy function
// -1.0*x*log(x)
double measures::entropy( double x ){

  if( x == 0 ) return 0;
  else return -1.0*(x*log(x));
  
}

// binary Entropy function
// -1.0*(x*log(x) + (1-x)*log(1-x)), or
// -x*log(x) - (1-x)*log(1-x)
double measures::binaryEntropy( double x ){

  return (entropy(x) + entropy(1-x));
  
}


//return log gamma of n
double measures::loggamma( int Nn ){

  return gsl_sf_lnfact( (const unsigned int)Nn );
  
}

//multiset coefficient
double measures::multico( int Nn, int Kk ){

  return binco( (Nn+Kk-1), Kk );
  
}

//binomial coefficient
double measures::binco( int Nn, int Kk ){

  double result, temp, num, dem;
  
  result = 0.0; num = 0.0; dem = 0.0; temp = 0.0;

  num = gsl_sf_lnfact( (const unsigned int)Nn );
 
  dem = gsl_sf_lnfact( (const unsigned int)Kk ) +
        gsl_sf_lnfact( (const unsigned int)(Nn-Kk) );
 
  temp = num - dem;
  
  if( fabs(temp) > MAXEXPO ){
    temp > 0 ? result = exp( MAXEXPO ) : result = exp( -MAXEXPO );
  } else {  
    result = exp( temp ) ;
  }    
  
  return result;
  
}

void measures::min(int d, double *xx, double &min){

  int i;
  min = xx[0];
  for(i=0; i<d; i++){
    if( xx[i] < min ) min = xx[i];
  }
  
}

void measures::max(int d, double *xx, double &max){

  int i;
  max = xx[0];
  for(i=0; i<d; i++){
    if( xx[i] > max ) max = xx[i];
  }
  
}

