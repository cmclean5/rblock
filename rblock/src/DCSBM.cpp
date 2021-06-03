#include "DCSBM.h"

/*
  Based on M.E.J. Newman & A. Clauset. "Structure and inference in annotated networks", nature communitions, 7:11863, doi:1038/ncomms11863.
 */

DCSBM::DCSBM() : SBM() { }

DCSBM::DCSBM( network *gg, int K ) : SBM(gg, K){ }

DCSBM::DCSBM( network *gg, int K, int SLOT ) : SBM(gg, K, SLOT){ }

DCSBM::DCSBM( network *gg, int K, int SLOT, int iter ) : SBM(gg, K, SLOT, iter){ }


DCSBM::~DCSBM(){ }


void DCSBM::initDCSBM(){

  initialise();
  initOmegaDCSBM();  
 
}

void DCSBM::initOmegaDCSBM(){

  int r,s;
  
  // Choose random values for the omegas, but with a bias toward
  // assortative choices (change if necessary for other networks)
  double c[K][K];
  for (r=0; r<K; r++) {
    for (s=0; s<K; s++) {
                    c[r][s] = 0.0;
      if (r==s)     c[r][s] = 1 + gsl_rng_uniform(g);//<-- bias, assortative choice 
      else if (r<s) c[r][s] = gsl_rng_uniform(g);
      else          c[r][s] = c[s][r];
      omega[r][s] = c[r][s]/(double)twom;
    }
  }

  
}


int DCSBM::runDCSBM(){

  int l,v,k,status,statusEM;

  status=1;

  //its += 1;
  
  cout << "run algorithm " << its << " time(s)." << endl;

  for(l=0; l<its; l++){

    initDCSBM();
    statusEM   = 1;
    statusEM   = emDCSBM();    	 
    	 
    if( l == 0 ){
      saved_LL       = LL;
      saved_steps    = STEPS;
      saved_maxDelta = MAXDELTA;
      
      if( model == 0 ) saveqCategoric();
      else             saveqContinuous();
     
      continue;
    }     

    if( statusEM == 0 && gsl_finite(LL) ){

      if( fabs(MAXDELTA) < fabs(saved_maxDelta) ){
	saved_LL       = LL;
	saved_steps    = STEPS;
	saved_maxDelta = MAXDELTA;
	status         = 0;

	cout << "save run!" << endl;

	if( model == 0 ) saveqCategoric();
	else             saveqContinuous();

      }
    }
	 
      
  }//for
 
  
  cout << "Done." << endl;

  return status;

}


int DCSBM::emDCSBM(){

  int u,v,i,r,s,status;
  int steps,bpsteps;
  double norm,deltac,maxdelta;
  double c[K][K];
  double oldc[K][K];
  double ru[K];

  status=0;
  steps=0;
  bpsteps=0;
  maxdelta=1;
  LL=0;

  for (r=0; r<K; r++) {
           ru[r] = 0; 
    for (s=0; s<K; s++) {
         c[r][s] = 0;
      oldc[r][s] = 0;
    }
  }

  //---EM loop
  //cout << "Starting EM algorithm..." << endl;
  while( maxdelta > EM_ACC ){

    // Run BP to calculate the messages and one-vertex marginals
    bpsteps = bpDCSBM();//doBeliefPropagation();

    LL = 0.0;
    // Calculate the new values of the parameters
    LL = updateParmDCSBM();

    // Calculate the new values of the c variables
    for (r=0; r<K; r++) {
      for (s=0; s<K; s++) {
	oldc[r][s] = c[r][s];
	c[r][s] = omega[r][s]*(double)twom;
      }
    }

    // Find the largest change in any of the c's
    maxdelta = 0.0;
    for (r=0; r<K; r++) {
      for (s=0; s<K; s++) {
	deltac = fabs(c[r][s]-oldc[r][s]);
        if (deltac>maxdelta) maxdelta = deltac;
      }
    }    
  
    // Print out new values of the parameters
    //cout << "EM step " << steps << " , max change = " << maxdelta << endl;
   
    if( ++steps > EM_MAXSTEP) {
      cout << "Solution failed to converge in " << EM_MAXSTEP 
	   << " EM steps" << endl;
      break;
    }

  }//while
    
  // Print out new values of the parameters
  cout << "EM step " << steps << " , max change = " << maxdelta << endl;

  STEPS    = steps;
  MAXDELTA = maxdelta;
  
  if( bpsteps > BP_MAXSTEP ){
    cout << "BP failed converge on final EM step" << endl;
    status = 2;
  }

  cout << "Log-likelihood = " << LL << endl; 

  if( maxdelta == 0 ) status = 1;
  
  return status;

}




 //---Function to calculate new values of the parameters
 double DCSBM::updateParmDCSBM(){

    int u,v,i,j,r,s;
    double norm,esum,quvrs;
    double d[K] = {};
    double term[K][K];
    double sum[K][K];
    double tsum[(K*K)];
    double L;

    
    //---Calculate some basics

    
    //---Zero out the sum and d variables
    for (r=0; r<K; r++) {
      d[r]  = 0.0;
      for (s=0; s<K; s++){
	sum[r][s]     = 0.0;
	tsum[(r*K)+s] = 0.0;
      }
    }

    

    //#pragma omp parallel for schedule(static) default(none) firstprivate(K,N) private(r,u) shared(d,q,gg)
    #pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(K,N) private(r,u) shared(q,gg) reduction(+:d[:K])
    for (r=0; r<K; r++) {
      for (u=0; u<N; u++) {    
	d[r] += q[u][r]*gg->V[u].degree;
      }
    }
        
    
    if( model == 0 ) updateCategoric();
    else             updateContinuous();

  
  
  //---Calculate the new values of the omegas
 

  esum = 0.0;
  //---Perform the sums
  #pragma omp parallel for schedule(guided) default(none) firstprivate(N,K) private(u,v,r,s,i,j,norm,quvrs,term) shared(gg,eta,omega) reduction(+: esum) reduction(+:tsum[:(K*K)])
  for (u=0; u<N; u++) {
    for (i=0; i<gg->V[u].degree; i++) {
      v = gg->V[u].E[i].target;

      //---Find which edge leads back from v to u
      for (j=0; j<gg->V[v].degree; j++) {
	if (gg->V[v].E[j].target==u) break;
      }
      if (j==gg->V[v].degree) {
	//fprintf(stderr,"Error!\n");
	//exit(23);
	break;
      }

      //---Calculate the terms and the normalization factor
      // omega[r][s]*eta[u][i][r]*eta[v][j][s];
      // P(community mixing) * P(node u in coms r, if node v missing) 
      //                     * P(node v in coms s, if node u missing) 
      norm = 0.0;
      for (r=0; r<K; r++) {
	for (s=0; s<K; s++) {
	  term[r][s] = omega[r][s]*eta[u][i][r]*eta[v][j][s];
	  norm += term[r][s];
	}
      }

      //---Add to the running sums
      for (r=0; r<K; r++) {
	for (s=0; s<K; s++) {
          quvrs = term[r][s]/norm;
	  //sum[r][s] += quvrs;
	  tsum[(r*K)+s] += quvrs;
	  if( quvrs == 0.0 ) esum += 0;
	  else               esum += quvrs*log(quvrs);//eqn (18)
	}
      }
    }
  }

 
  //---Calculate the new values of the omega variables (after calculating
  //   the likelihood using the old omegas)
  for (r=0; r<K; r++) {
    for (s=0; s<K; s++){
      sum[r][s]   = tsum[(r*K)+s];
      omega[r][s] = sum[r][s]/(d[r]*d[s]);
    }
  }

  //---Calculate the expected log-likelihood

  //---Internal energy first
  L = 0.0;
  for (r=0; r<K; r++) {
    for (s=0; s<K; s++){
      if( omega[r][s] == 0.0 ) L += 0.5*sum[r][s]*log(omega[r][s]+SMALL);
      else                     L += 0.5*sum[r][s]*log(omega[r][s]);
    }

    if( model == 0 ) energyCategoric ( r, L );
    else             energyContinuous( r, L );
    
  
  }
  
  
  //---Now the entropy
  L -= 0.5*esum;
  for (u=0; u<N; u++) {
    for (r=0; r<K; r++) {      
      if (gg->V[u].degree>0) {	
	if( q[u][r] == 0.0 ) L += 0.0;
	else                 L += (gg->V[u].degree-1)*q[u][r]*log(q[u][r]);
      }
    }
  }

  return L;

  }


//--- Minimum Description Length for DCBM
//    MDL = DL_DCSBM + entropy_DCSBM
// ref: Tiago Peixoto, Physical Review X, 4, 011047 (2014) doi:10.1103/PhysRevX.4.011047
double DCSBM::mdlDCSBM(){

  int i,k,u,v,r,s,Kmax,ksum,ktot;  
  double ers,dqur,mdl=0;
  double qurmax;
  double nr[K];
  int    nK[N];
  double d[K];
 
  //--- Calculate some basics
  for (u=0, Kmax=0; u<N; u++)    
    if( gg->V[u].degree > Kmax ) Kmax = (int)gg->V[u].degree;

  //--- Calculate description length, DL_DCSBM  
  
  //--- emperical degree distribution of communities
  //--- some basics first
  /*
  for(u=0; u<N; u++){
    nK[u]  = 0;
    qurmax = q[u][0];
    for(r=1; r<K; r++){
      if( saved_q[u][r] > qurmax ){ nK[u] = r; qurmax = saved_q[u][r]; } 
    }
  }

  //--- now calculate the emperical degree distribution of communities
  int distK[Kmax];  
  for(k=0; k<Kmax; k++) distK[k] = 0;  
  for (r=0, ktot=0; r<K; r++) {
    nr[r] = 0.0;
    d[r]  = 0.0;    
    for (u=0; u<N; u++) {    
      nr[r] += saved_q[u][r];
      d[r]  += saved_q[u][r] * gg->V[u].degree;
      if( nK[u] == r ){
	for (i=0, ksum=0; i<gg->V[u].degree; i++){
	  v = gg->V[u].E[i].target;
	  if( nK[u] == nK[v] ) ksum++; 
	}
	distK[ksum]++;
	ktot += ksum;
      }
    }//n
    if( ktot > 0 ){
      for(k=0; k<Kmax; k++){
	mdl += (int)nr[r] * entropy( (double)(distK[k]*k)/(double)ktot );
	distK[k]=0;
      }
    } else distK[0]=0;
  }//r
  //--- done.
  */

  //--- the expected degree distribution of communities
  int distK[Kmax];  
  for(k=0; k<Kmax; k++) distK[k] = 0;  
  for (r=0; r<K; r++) {
    nr[r] = 0.0;
    d[r]  = 0.0;    
    for (u=0; u<N; u++) {    
      nr[r] += saved_q[u][r];
      dqur   = 0.0;
      dqur   = saved_q[u][r] * gg->V[u].degree;
      d[r]  += dqur;
      distK[(int)dqur]++;
    }
    if( d[r] > 0 ){
      for(k=0; k<Kmax; k++){
	mdl += (int)nr[r] * entropy( (double)(distK[k]*k)/(double)d[r] );
	distK[k]=0;
      } 
    } else distK[0] = 0;
  }
  //--- done.

  mdl += log( multico(0.5*K*(K+1),M) ) + log( multico(K,N) ) +
    loggamma(N);

  for(r=0; r<K; r++)
    mdl -= loggamma( (int)nr[r] );

  //--- Done with description length, DL_DCSBM 
  

  //--- Calculate entropy, entropy_DCSBM  

  //degree distribution of network
  for(k=0; k<Kmax; k++) distK[k] = 0;
  for(u=0; u<N; u++) distK[(int)gg->V[u].degree]++;
  for(k=0; k<Kmax; k++) mdl -= distK[k] * loggamma( (int)k );
  
  for( r=0, ers=0.0; r<K; r++ ){
    for( s=0; s<K; s++ ){
      ers +=  0.5 * omega[r][s] * d[r] * d[s];
      mdl -= (0.5 * omega[r][s] * d[r] * d[s]) * log( omega[r][s] );
    }
  }
  
  mdl -= ers;
  //--- Done with entropy, entropy_SBM
  
  return mdl;

}
    
