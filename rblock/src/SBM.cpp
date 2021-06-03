#include "SBM.h"

// Program based on Newman's 2014 code for a standard SBM
/* Program to perform the full two-group EM/BP calculation using the
 * standard SBM on an arbitrary network read from a file
 *
 * Written by Mark Newman  21 NOV 2014
 * Modified to include enumerative metadata stored in the "label" field of
 *   the network structure  23 NOV 2014
 */


SBM::SBM() : BP() { setup(); }

SBM::SBM( network *gg, int K ) : BP(gg, K){

  setup(); 
  assignSpace();
  
}


SBM::SBM( network *gg, int K, int SLOT ) : BP(gg, K, SLOT){

  setup();  
  assignSpace();
  
}

SBM::SBM( network *gg, int K, int SLOT, int iter ) : BP(gg, K, SLOT){

  setup();  
  this->its = ( iter <= 0 ) ? 1 : iter;
  assignSpace();
  
}


void SBM::setup(){

  this->its=10;

  this->EM_ACC=0.0001;
  this->EM_MAXSTEP=100;//1000;
  this->KSMALL=1.0e-3;
  this->NOCONVERGE=false;
  this->calVprob=false;
  this->calEprob=false;  

  this->saved_q=0;
  this->qz=0;
  
  this->LL=0;
  this->saved_LL=0;

  this->STEPS=100;
  this->MAXDELTA=1;

  this->saved_steps=100;
  this->saved_maxDelta=1;
  
  //--- Initialize gsl random number seed (in Measures)
  setSeed();
  //testSeed(1);//comment out... for testing only

  
}

void SBM::assignSpace(){

  int i;

  //--Set memory for qz, edge probabilities (mass action principle between two nodes) using network's edge-list.
  qz     = (double**)malloc(El*sizeof(double));
  for(i=0; i<El; i++){
    qz[i]  = (double*)malloc(K*sizeof(double));
  }

  //Make space for the marginals and initialize to random initial values
  //Randomly assign community probability to each node
  saved_q = (double**) malloc(N*sizeof(double));
  for (i=0; i<N; i++) {
    saved_q[i] = (double*)malloc(K*sizeof(double));
  }

}


SBM::~SBM(){ freeSpace(); }


void SBM::freeSpace(){

  int i;

  if(qz!=0)          { for (i=0; i<El; i++) { free(qz[i]); } }

  if( saved_q!=0 )   {for (i=0; i<N; i++) { free(saved_q[i]); } }

}


void SBM::useSeed( int newSeed ){ setSeed(true, newSeed); }

void SBM::initSBM(){

  initialise();
  initOmegaSBM();

}

void SBM::initialise(){

  int i,j,u,v,r;

  // Make space for the marginals and initialize to random initial values
  // Randomly assign community probability to each node
  // "marginals" are the non-conditional probabilities that node u belongs to 
  // community k. 
  for (u=0; u<N; u++) {
    random_unity(K,q[u]);
  }

  // Make space for the "messages" and initialize to the same values as the
  // "marginals"

  // messages are sent from one node to another along each edge (i,j).
  // "messages" are the eta's: they are the  probaility that node u would be in 
  // the group r if node v were absent from the network, i.e. eta[u][v][r].
  // each messages is initialised to the non-conditional probabilities that node
  // v belongs to community r, i.e. q[v][r].
  for (u=0; u<N; u++) {
    for (i=0; i<gg->V[u].degree; i++) {
      v = gg->V[u].E[i].target;
      for (r=0; r<K; r++){
	eta[u][i][r] = 0.0;
	eta[u][i][r] = q[v][r];
      }
    }
  }  

  // Initialise the Metadata model
  if( model == 0 ) initCategoric();
  else             initContinuous();
    

}

void SBM::initOmegaSBM(){

  int r,s;
  
  // Choose random values for the omegas, but with a bias toward
  // assortative choices (change if necessary for other networks)
  double c[K][K];
  double cnorm = (double)twom / (double)N;
  for (r=0; r<K; r++) {
    for (s=0; s<K; s++) {
      c[r][s] = 0.0;
      if (r==s)     c[r][s] = cnorm * (1 + gsl_rng_uniform(g));//<-- bias, assortative choice 
      else if (r<s) c[r][s] = cnorm * (gsl_rng_uniform(g));
      else          c[r][s] = c[s][r];
      omega[r][s] = c[r][s] / (double)N;
    }
  }
  
}

void SBM::initCategoric(){

  int i,r;
  double ru[K];
  
  // Malloc space for parameters gmma and choose random initial values
  // gamma_sx, Prior probability for node assigned to 'K' communities, 
  // using annotation 'x' data. 
  for (i=0; i<nmlabels; i++) {
    random_unity(K,ru);
    for (r=0; r<K; r++) gmma[r][i] = ru[r];
  }

  //For testing...
  //initContinuous();
  
}

void SBM::initContinuous(){

  int j,r;
  double ru[K];
  
  // Initialise Bernstein coefficients, subject to constain
  // Sum_s bernco_sj = 1
  for (j=0; j<=D; j++) {
    random_unity(K,ru);
    for (r=0; r<K; r++) bernco[r][j] = ru[r];
  }  
  
  //now initialise Qbern, and Qsu, using bernco_sj.
  updateQbern();
      
}

int SBM::runSBM(){

  int l,v,k,status,statusEM;

  status=1;

  //its += 1;
  
  cout << "run algorithm " << its << " time(s)." << endl;

  for(l=0; l<its; l++){

    initSBM();
    statusEM   = 1;
    statusEM   = emSBM();    	 
    	 
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
	//if( LL >= saved_LL ){
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


int SBM::emSBM(){

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
    bpsteps = bpSBM();//doBeliefPropagation();

    LL = 0.0;
    // Calculate the new values of the parameters
    LL = updateParmSBM();

    // Calculate the new values of the c variables
    for (r=0; r<K; r++) {
      for (s=0; s<K; s++) {
	oldc[r][s] = c[r][s];
	c[r][s] = omega[r][s] * (double)N;
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
 double SBM::updateParmSBM(){

    int u,v,i,j,r,s;
    double norm,esum,quvrs;
    double nr[K];
    double term[K][K];
    double sum[K][K];
    double L;

    //---Calculate some basics
    for (r=0; r<K; r++) {
      nr[r] = 0.0;
      for (u=0; u<N; u++) {    
	nr[r] += q[u][r];
      }
    }

    if( model == 0 ) updateCategoric();
    else             updateContinuous();

  
  
  //---Calculate the new values of the omegas

  //---Zero out the sum variables
  for (r=0; r<K; r++) {
    for (s=0; s<K; s++) sum[r][s] = 0.0;
  }

  esum = 0.0;
  //---Perform the sums
  for (u=0; u<N; u++) {
    for (i=0; i<gg->V[u].degree; i++) {
      v = gg->V[u].E[i].target;

      //---Find which edge leads back from v to u
      for (j=0; j<gg->V[v].degree; j++) {
	if (gg->V[v].E[j].target==u) break;
      }
      if (j==gg->V[v].degree) {
	fprintf(stderr,"Error!\n");
	exit(23);
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
	  sum[r][s] += quvrs;
	  if( quvrs == 0.0 ) esum += 0;
	  else               esum += quvrs*log(quvrs);//eqn (18)
	}
      }
    }
  }

  //---Calculate the new values of the omega variables (after calculating
  //   the likelihood using the old omegas)
  for (r=0; r<K; r++) {
    for (s=0; s<K; s++) omega[r][s] = sum[r][s]/(nr[r]*nr[s]);
  }

  //---Calculate the expected log-likelihood

  //---Internal energy first
  L = 0.0;
  for (r=0; r<K; r++) {
    for (s=0; s<K; s++){
      if( omega[r][s] == 0.0 || omega[r][s] == 1.0 ) L += 0.0;
      else 
	L += (nr[r]*nr[s])*(   omega[r][s]  * log(    omega[r][s]) +
			    (1-omega[r][s]) * log((1-omega[r][s])) );
      
    }
    
    if( model == 0 ) energyCategoric ( r, L );
    else             energyContinuous( r, L );       
  }
  
  
  //---Now the entropy
  ////L -= 0.5*esum;
  for (u=0; u<N; u++) {
    for (r=0; r<K; r++) {            	
	if( q[u][r] == 0.0 ) L -= 0.0;
	else                 L -= q[u][r]*log(q[u][r]);
      }
    }
 

  return L;

 }


void SBM::energyCategoric( int r, double &E ){

  int i;
  
  for (i=0; i<nmlabels; i++) {
    if( gmma[r][i] == 0.0 ) E += 0.0;
    else                    E += nx[i]*gmma[r][i]*log(gmma[r][i]);
  }
  
}

void SBM::energyContinuous( int r, double &E ){

  int u;
  
  for(u=0; u<N; u++){
    if( Qsu[r][u] == 0.0 ) E += 0.0;
    else                   E += Qsu[r][u]*log(Qsu[r][u]);
  }

}

void SBM::updateCategoric(){

  int i,r,u;
  
  for (r=0; r<K; r++) {
    for (i=0; i<nmlabels; i++) nrx[r][i] = 0.0;
  }

  for (u=0; u<N; u++) {
    for (r=0; r<K; r++) {
      nrx[r][x[u]] += q[u][r];
    }
  }

  //---Calculate new values of the gammas
  for (r=0; r<K; r++) {
    for (i=0; i<nmlabels; i++) gmma[r][i] = nrx[r][i]/nx[i];
  }


  //For testing...
  //updateContinuous();
  
}


void SBM::updateContinuous(){

  //---Calculate new bernco, eqn-29, using the new q values
  updateBernco();
  
  //---Update Qbern, eqn-28, and Qsu, after updating bernco.
  updateQbern();
  
}

//---Calculate new bernco, eqn-29, and Qsu here.
void SBM::updateBernco(){

  int u,r,j;  
  double bernco_norm[(D+1)];

  for(j=0; j<=D; j++){
    bernco_norm[j] = 0;
    for(r=0; r<K; r++){
      bernco[r][j] = 0;
      for(u=0; u<N; u++){
	bernco[r][j]   += q[u][r] * Qbern[u][j][r];
	bernco_norm[j] += q[u][r] * Qbern[u][j][r];
      }
    }
  }

  //cout << "bernco " << endl;
  for(j=0; j<=D; j++){
    for(r=0; r<K; r++){
      if( bernco[r][j] == 0 && bernco_norm[j] == 0 ){
	;
      } else {
	bernco[r][j] /= bernco_norm[j];
      }

    }
  }


}

// Update Qbern, eqn-28, and Qsu after updating bernco.
void SBM::updateQbern(){

  int u,r,j;
  double norm;
  
  for(u=0; u<N; u++){   
    for(r=0; r<K; r++){      
      norm = 0.0;
      for(j=0; j<=D; j++){norm += (bernco[r][j] * bvec[u][j]);}
      for(j=0; j<=D; j++){	
	Qbern[u][j][r] = 0.0;
	Qbern[u][j][r] = bernco[r][j] * bvec[u][j];
	if( Qbern[u][j][r] == 0 && norm == 0 ){
	  ;
	} else {
	  Qbern[u][j][r] /= norm;
	}	

      }    
    }
  }

  for(u=0; u<N; u++){
    for(r=0; r<K; r++){
      Qsu[r][u] = 0.0;
      for(j=0; j<=D; j++){
  	Qsu[r][u] +=  bernco[r][j] * bvec[u][j];
      }
    }
  } 

  
}

void SBM::saveqCategoric(){

  int u,r;

  for(u=0; u<N; u++){
    for(r=0; r<K; r++){
      saved_q[u][r] = 0.0;
      saved_q[u][r] = q[u][r];
    }
  }
  
}

void SBM::saveqContinuous(){

  int u,r;

  for(u=0; u<N; u++){
    for(r=0; r<K; r++){
      saved_q[u][r] = 0.0;
      saved_q[u][r] = Qsu[r][u];//q[u][r];//Qsu[r][u];
    }
  }
  
}


void SBM::getVertexKProbs(){

    int i,k;
    double max;
    double sum[N];

    for(i=0; i<N; i++){
      sum[i] = 0;
      for(k=0; k<K; k++){
	gg->V[i].Kprobs[k] = 0;
	if( q[i][k] > KSMALL )  gg->V[i].Kprobs[k] = saved_q[i][k];
	sum[i] += gg->V[i].Kprobs[k];
      }
    }

    for(i=0; i<N; i++){
      max = gg->V[i].Kprobs[0]/sum[i];
      gg->V[i].K=0;//1;
      for(k=0; k<K; k++){
	  gg->V[i].Kprobs[k] /= sum[i];
	  if( gg->V[i].Kprobs[k] > max ) { max=gg->V[i].Kprobs[k]; gg->V[i].K=k; }//(k+1); }
      }  
    }

    calVprob = true;

}



void SBM::getEdgeKProbs(){

  int m,k;

  //--- Calculate the probability for the edges
  for(m=0; m<El; m++){
      
    edgelist current = el[m];
      
    int ind_so = current.source;
    int ind_si = current.target;

    for(k=0; k<K; k++){

      double Pa   = gg->V[ind_so].Kprobs[k];
      double Pb   = gg->V[ind_si].Kprobs[k];
	
      qz[m][k] = 0.0;
      qz[m][k] = Pa*Pb;

    }
  }
  
  calEprob = true;

  }


//--- Minimum Description Length for SBM
//    MDL = DL_SBM + entropy_SBM
// ref: Tiago Peixoto, Physical Review X, 4, 011047 (2014) doi:10.1103/PhysRevX.4.011047
double SBM::mdlSBM(){

  int u,r,s;
  double mdl=0;
  double nr[K];
  double c[K][K];
 
  //--- Calculate description length, DL_SBM
  for (r=0; r<K; r++) {
    nr[r] = 0.0;
    for (u=0; u<N; u++) {    
      nr[r] += saved_q[u][r];
    }
  }
  
  mdl += log( multico(0.5*K*(K+1),M) ) + log( multico(K,N) ) +
    loggamma(N);

  for(r=0; r<K; r++)
    mdl -= loggamma( (int)nr[r] );

  //--- Done with description length, DL_SBM 


  //--- Calculate entropy, entropy_SBM
  for( r=0; r<K; r++ ){
    for( s=0; s<K; s++ ){
      mdl += 0.5 * nr[r] * nr[s] * binaryEntropy( omega[r][s] ); 
    }
  }

  //--- Done with entropy, entropy_SBM
  
  return mdl;

}
    
void SBM::checkValue( double x ){

  if( isnan(x) ){ cout << "isnan!"; }

  if( x == 0 ){ cout << "is zero!"; }  

  cout << "" << endl;
  
}

double SBM::getMaxLL(){
  return saved_LL;
}

double SBM::getEMsteps(){
  return saved_steps;
}

double SBM::getMaxDelta(){
  return saved_maxDelta;
}
