#include "BP.h"

BP::BP() : metaData(){ setup(); }

BP::BP( network *gg ) : metaData(gg) {
  setup();
  assignSpace();
}

BP::BP( network *gg, int K ) : metaData(gg, K) {
  setup();
  assignSpace();
}


BP::BP( network *gg, int K, int annoSLOT ) : metaData(gg, K, annoSLOT) {
  setup();
  assignSpace();
}

void BP::setup(){

  this->BP_MAXSTEP=20;
  this->BP_ACC=0.0001;
  this->SMALL=1.0e-100;

  this->eta=0;
  this->q=0;
  this->omega=0;
  
}

void BP::assignSpace(){

  int i,u,v,r,s;

  //---assign K communities to each vertices
  for(i=0; i<N; i++){ gg->V[i].assignKprobs(K); }

  
  //Make space for the marginals and initialize to random initial values
  //Randomly assign community probability to each node
  q       = (double**) malloc(N*sizeof(double));
  for (u=0; u<N; u++) {
    q[u]       = (double*)malloc(K*sizeof(double));
  }

  // Make space for the messages and initialize to the same values as the
  // marginals
  eta = (double***)malloc(N*sizeof(double));
  for (u=0; u<N; u++) {
    eta[u] = (double**)malloc(gg->V[u].degree*sizeof(double));
    for (i=0; i<gg->V[u].degree; i++) {
      eta[u][i] = (double*)malloc(K*sizeof(double));
    }
  }
  
  // Choose random values for the omegas, but with a bias toward
  // assortative choices (change if necessary for other networks)
  //double c[K][K];
  omega = (double**)malloc(K*sizeof(double));  
  for (r=0; r<K; r++) {
    omega[r] = (double*)malloc(K*sizeof(double));
  }

}

void BP::freeSpace(){

  int i,u;

  for(i=0; i<N; i++) { gg->V[i].freeKprobs(); }

  if( q!=0 )         {for (i=0; i<N; i++) { free(q[i]); } }

  if( eta!=0 ){
    for (u=0; u<N; u++) {
      for (i=0; i<gg->V[u].degree; i++) free(eta[u][i]);
      free(eta[u]);
    }
    free(eta);
  }  

  if( omega!=0 )    {for (i=0; i<K; i++) { free(omega[i]); } }

}

BP::~BP(){ freeSpace(); }


//---Do Belief Propagation
int BP::bpSBM(){

  int i,j,u,v,r,s;
  int    steps;
  double deltaeta,maxdelta;
  double logeta,neweta,sum,norm,largest, logmma;
  double nr[K];
  double logpre[K];
  double logqun[K];
  double ***logetaun; // [ NxDegreexK ] 

  
  // Make space for the new log-etas, which are called "logetaun" because
  // the are initially calculated in unnormalized form
  logetaun = (double***)malloc(N*sizeof(double));
  for (u=0; u<N; u++) {
    logetaun[u] = (double**)malloc(gg->V[u].degree*sizeof(double));
    for (i=0; i<gg->V[u].degree; i++) {
      logetaun[u][i] = (double*)malloc(K*sizeof(double));
    }
  }
 

  // Main BP loop
  steps    = 0;
  maxdelta = 1;
  while( (maxdelta > BP_ACC) && (steps++ <= BP_MAXSTEP) ){

    //Calculate the expected group size    
    for (r=0; r<K; r++) {
      nr[r] = 0.0;
      for (u=0; u<N; u++) nr[r] += q[u][r];
    }
    
    //---Calculate the log-prefactors... i.e. the auxiliary field for non-edges
    for (r=0; r<K; r++) {
      logpre[r] = 0.0;
      for (s=0; s<K; s++) logpre[r] -= omega[r][s]*nr[s];
    }

      
    //---Calculate new values for the one-vertex marginals
    for (u=0; u<N; u++) {
      for (r=0; r<K; r++) {

	if( model == 0 ){
	  // For categoric data we use gmma...
	  logqun[r] = log(gmma[r][x[u]]) + logpre[r];
	} else {
	  // For continuous data we'd use Qbern, i.e. Qsu...
	  logqun[r] = log(Qsu[r][u]) + logpre[r];
	}
	  
	for (i=0; i<gg->V[u].degree; i++) {
	  sum = 0.0;
	  for (s=0; s<K; s++) sum += eta[u][i][s]*omega[r][s];
	  if (sum<SMALL) sum = SMALL;
	  logqun[r] += log(sum);
	}
	if (r==0) largest = logqun[r];
	else if (logqun[r]>largest) largest = logqun[r];
      }
     
      //---Normalize
      norm = 0.0;
      for (r=0; r<K; r++) {
	logqun[r] -= largest;
	norm += exp(logqun[r]);
      }
      
      //new values for the one-vertex marginals, q
      for (r=0; r<K; r++) q[u][r] = exp(logqun[r])/norm;
      
    }
    
    //---Calculate (unnormalized) new values for the (log) messages
    for (u=0; u<N; u++) {
      for (i=0; i<gg->V[u].degree; i++) {
	v = gg->V[u].E[i].target;
	for (r=0; r<K; r++) {

	  if( model == 0 ){
	    logeta = log(gmma[r][x[v]]) + logpre[r];
	  } else {	    
	    // For coninuous data we'd use Qbern, i.e. Qsu, here...
	    logeta = log(Qsu[r][v]) + logpre[r];
	  }
	    
	  for (j=0; j<gg->V[v].degree; j++) {
	    if (gg->V[v].E[j].target!=u) {
	      sum = 0.0;
	      for (s=0; s<K; s++) sum += eta[v][j][s]*omega[r][s];
	      if (sum<SMALL) sum = SMALL;   // Prevent -Inf
	      logeta += log(sum);
	    }
	  }
	  logetaun[u][i][r] = logeta;
	}
      }
    }

    //---Normalize and calculate largest change 
    maxdelta = 0.0;
    for (u=0; u<N; u++) {
    for (i=0; i<gg->V[u].degree; i++) {
	norm = 0.0;
	largest = logetaun[u][i][0];
	for (r=1; r<K; r++) {
	  if (logetaun[u][i][r]>largest) largest = logetaun[u][i][r];
	}
	for (r=0; r<K; r++) {
	  logetaun[u][i][r] -= largest;
	  norm += exp(logetaun[u][i][r]);
	}	  
	for (r=0; r<K; r++) {
	  neweta = exp(logetaun[u][i][r])/norm;
	  deltaeta = fabs(neweta-eta[u][i][r]);
	  if (deltaeta>maxdelta) maxdelta = deltaeta;
	  eta[u][i][r] = neweta;
	}
      }
    }

    //cout << "done." << endl;
    
    //cout << "BP steps " << steps << " max change = " << maxdelta << endl;

  }//while 


  //---Free space
  for (u=0; u<N; u++) {
    for (i=0; i<gg->V[u].degree; i++) free(logetaun[u][i]);
    free(logetaun[u]);
  }
  free(logetaun);
  //---End Free space
  
  return steps;

}



//---Do Belief Propagation
//int BP::doBeliefPropagation(){
int BP::bpDCSBM(){

  int i,j,u,v,r,s;
  int    steps;
  double deltaeta;
  double logeta,neweta,sum,norm,largest, logmma;
  double d[K] = {};
  double logpre[K] = {};
  double logqun[K];
  double ***logetaun; // [ NxDegreexK ] 

			   
  //struct to store max delta
  struct RESULT { double param1; };

  //we'll need to define our own reduction clause
  #pragma omp declare reduction(storeMax : RESULT : omp_in.param1 > omp_out.param1 ? omp_out = omp_in : omp_out)
			   
			   
  // Make space for the new log-etas, which are called "logetaun" because
  // the are initially calculated in unnormalized form
  logetaun = (double***)malloc(N*sizeof(double));
  for (u=0; u<N; u++) {
    logetaun[u] = (double**)malloc(gg->V[u].degree*sizeof(double));
    for (i=0; i<gg->V[u].degree; i++) {
      logetaun[u][i] = (double*)malloc(K*sizeof(double));
    }
  }
 

  // Main BP loop
  steps    = 0;
  RESULT maxdelta{1};
   
  while( (maxdelta.param1 > BP_ACC) && (steps++ <= BP_MAXSTEP) ){
 
    for(r=0; r<K; r++){
      d[r]      = 0.0;
      logpre[r] = 0.0;
    }
    
    //Calculate the expected group degrees
    //part of the auxiliary field for non-edges...    
#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(K,N) private(r,u) shared(q,gg) reduction(+:d[:K])
    for (r=0; r<K; r++) {
      for (u=0; u<N; u++) d[r] += q[u][r]*gg->V[u].degree;
    }

    //---Calculate the log-prefactors (without the leading factor of d_i or
    //   the prior)
    #pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(K) private(r,s) shared(omega,d) reduction(+:logpre[:K])
    for (r=0; r<K; r++) {
      for (s=0; s<K; s++) logpre[r] -= omega[r][s]*d[s];
    }

      
    //---Calculate new values for the one-vertex marginals
    #pragma omp parallel for schedule(static) default(none) firstprivate(N,K,model,SMALL) private(u,r,i,s,sum,largest,norm,logqun) shared(gmma,x,gg,logpre,Qsu,eta,omega,q)
    for (u=0; u<N; u++) {
      for (r=0; r<K; r++) {

	if( model == 0 ){
	  // For categoric data we use gmma...
	  // First term log(..) is the prior, second term part of the auxiliary field for non-edges
	  logqun[r] = log(gmma[r][x[u]]) + gg->V[u].degree*logpre[r];
	} else {
	  // For continuous data we'd use Qbern, i.e. Qsu...
	  logqun[r] = log(Qsu[r][u]) + gg->V[u].degree*logpre[r];
	}
	  
	for (i=0; i<gg->V[u].degree; i++) {
	  sum = 0.0;
	  for (s=0; s<K; s++) sum += eta[u][i][s]*omega[r][s];
	  if (sum<SMALL) sum = SMALL;
	  logqun[r] += log(sum);
	}
	if (r==0) largest = logqun[r];
	else if (logqun[r]>largest) largest = logqun[r];
      }
     
      //---Normalize
      norm = 0.0;
      for (r=0; r<K; r++) {
	logqun[r] -= largest;
	norm += exp(logqun[r]);
      }
      
      //new values for the one-vertex marginals, q
      for (r=0; r<K; r++) q[u][r] = exp(logqun[r])/norm;
      
    }
    

    //---Calculate (unnormalized) new values for the (log) messages
    #pragma omp parallel for schedule(guided) default(none) firstprivate(N,K,model,SMALL) private(u,v,r,i,j,s,sum,logeta) shared(gmma,x,gg,logpre,Qsu,eta,omega,logetaun)
    for (u=0; u<N; u++) {
      for (i=0; i<gg->V[u].degree; i++) {
	v = gg->V[u].E[i].target;
	for (r=0; r<K; r++) {

	  if( model == 0 ){
	    logeta = log(gmma[r][x[v]]) + gg->V[v].degree*logpre[r];
	  } else {	    
	    // For coninuous data we'd use Qbern, i.e. Qsu, here...
	    logeta = log(Qsu[r][v]) + gg->V[v].degree*logpre[r];
	  }
	    
	  for (j=0; j<gg->V[v].degree; j++) {
	    if (gg->V[v].E[j].target!=u) {
	      sum = 0.0;
	      for (s=0; s<K; s++) sum += eta[v][j][s]*omega[r][s];
	      if (sum<SMALL) sum = SMALL;   // Prevent -Inf
	      logeta += log(sum);
	    }
	  }
	  logetaun[u][i][r] = logeta;
	}
      }
    }

    
    //---Normalize and calculate largest change   
    maxdelta.param1 = 0.0;
    #pragma omp parallel for schedule(static) default(none) firstprivate(N,K) private(u,r,i,norm,largest,neweta,deltaeta) shared(gg,logetaun,eta) reduction(storeMax:maxdelta)
    for (u=0; u<N; u++) {
    for (i=0; i<gg->V[u].degree; i++) {
	norm = 0.0;
	largest = logetaun[u][i][0];
	for (r=1; r<K; r++) {
	  if (logetaun[u][i][r]>largest) largest = logetaun[u][i][r];
	}
	for (r=0; r<K; r++) {
	  logetaun[u][i][r] -= largest;
	  norm += exp(logetaun[u][i][r]);
	}	  
	for (r=0; r<K; r++) {
	  neweta = exp(logetaun[u][i][r])/norm;
	  deltaeta = fabs(neweta-eta[u][i][r]);	 
	  if (deltaeta > maxdelta.param1) maxdelta = RESULT{deltaeta};
	  eta[u][i][r] = neweta;
	}
      }
    }

    //cout << "BP steps " << steps << " max change = " << maxdelta << endl;

  }//while 


  //---Free space
  for (u=0; u<N; u++) {
    for (i=0; i<gg->V[u].degree; i++) free(logetaun[u][i]);
    free(logetaun[u]);
  }
  free(logetaun);
  //---End Free space
  
  return steps;

}
