#include "Headers.h"
#include "readfile.h"
#include "network.h"
#include "DCSBM.h"

//----------------
//NOTE:
// When editing this file, i.e. adding or removing functions, need to carry out the following tasks:
// 1) edit file NAMESPACE, adding/removing the function name in export(...)
// Next run following in R, inside the package:
// 2) $> cd rblock/
// 3) $> R
// 4)  > library(Rcpp)
// 5)  > Rcpp::compileAttributes()
// running 5) will create new RcppExports.R in R/ and
// new RcppExports.cpp in src/
// If need to run (to build Makevars file in src):
// 1) $> cd rblock/
// 2) $> autoconfig
// 3) $> ./configure
// BUILD PACKAGE:
// 1) $> R CMD build rblock/
// INSTALL PACKAGE (note 'INSTALL' should be in capital letters):
// 1) $> R CMD INSTALL rblock_1.0.tar.gz
//----------------

// Global
network *gg=0;
DCSBM   *model=0;
int      verbose=0;

// [[Rcpp::export]]
void resetGraph(){ gg = 0; }

// [[Rcpp::export]]
void resetModel(){ model = 0; }

// [[Rcpp::export]]
void deleteGraph(){ if(gg != 0){ delete gg; }  }

// [[Rcpp::export]]
void deleteModel(){ if(model != 0){ delete model; }  }

// [[Rcpp::export]]
void print( Rcpp::IntegerVector PRINT=0 ){ if( PRINT[0] == 1 ) verbose = 1; }

// [[Rcpp::export]]
void loadGraph( Rcpp::DataFrame     DF,
		Rcpp::IntegerVector directed=0,
		Rcpp::IntegerVector names=1 ){

  int i,j,k,KK;

  bool useLoops    = false;
  bool checkM      = true;
  int alphaNumeric = 1;
  int Directed     = directed[0];
  
  //initialise network
  resetGraph();
  readfile *reader = 0;
  gg               = new network();
  gg->directed     = Directed;
  
  int ncols = DF.length();
  int nrows = DF.nrows();

  if( (ncols > 0) && (nrows > 0) ){

    //set size for our DATASET
    KK              = nrows*ncols;
    string *DATASET = new string[KK];
        
    
    if( verbose == 1 ) gg->setPrint( true );
    
    if( (names.length() == 1) ){
      if( names[0] == 0 ){
      	alphaNumeric = 0;
      }
    }
    

    if( ncols == 2 ){
      //unweighted networks
      
      Rcpp::StringVector       V1 = DF[0];
      Rcpp::StringVector       V2 = DF[1];
      
      for(k=0; k<KK; k++){
	i = floor(k/ncols);
	j = k % ncols;

	Rcpp::String v1(V1[i]);
	Rcpp::String v2(V2[i]);
	
	if( j == 0 ){ DATASET[(i*ncols)+j] = v1.get_cstring(); }
	if( j == 1 ){ DATASET[(i*ncols)+j] = v2.get_cstring(); }
	
      }    
      
    }

    if( ncols == 3 ){
      //weighted networks.
      
      Rcpp::StringVector       V1 = DF[0];
      Rcpp::StringVector       V2 = DF[1];
      Rcpp::StringVector       V3 = DF[2];
      
      for(k=0; k<KK; k++){
	i = floor(k/ncols);
	j = k % ncols;

	Rcpp::String v1(V1[i]);
	Rcpp::String v2(V2[i]);
	Rcpp::String v3(V3[i]);
	
	if( j == 0 ){ DATASET[(i*ncols)+j] = v1.get_cstring(); }
	if( j == 1 ){ DATASET[(i*ncols)+j] = v2.get_cstring(); }
	if( j == 2 ){ DATASET[(i*ncols)+j] = v3.get_cstring(); }
	
      }    
      
    }
    
    
    //load edgelist into network
    reader = new readfile( gg, DATASET, ncols, nrows, alphaNumeric );

    //build Adjaceny Matrix
    //gg->buildAdj( checkM );
    //---

    //gg->findInOutDegrees();
    //gg->printVertices();

  }
 
 
}//load network


// [[Rcpp::export]]
void printGraph(){

  if( gg != 0 ){

    cout << "  N : " << gg->getN()
	 << ", M : " << gg->getM() 
         << ", M2: " << gg->getM2() << endl;

  }

  //deleteGraph();
  
}

// [[Rcpp::export]]
void loadMetaData( Rcpp::DataFrame DF, Rcpp::IntegerVector MODELS ){//, Rcpp::IntegerVector PRINT=0 ){

  int i,j,u,N;
  
  int ncols   = DF.length();
  int nrows   = DF.nrows();
  int slots   = ncols-1;
  int nmodels = MODELS.length();
  //int print   = PRINT[0];
  
  if( (ncols > 0) && (nrows > 0) && (gg!=0) ){

    if( nmodels != slots ){ cout << "Number of Metadata columns in DF not equal to MODELS vector length."; return; }

    if( verbose == 1 ) gg->setPrint( true );
    
    gg->assignNmodels( slots );
    for(i=0; i<slots; i++) gg->assignModel(i, MODELS[i]);   

    N = gg->getN();
    for( u=0; u<N; u++ ){ gg->V[u].assignNannos( slots ); }     
    
    Rcpp::StringVector V1 = DF[0];

    for( j=0; j<slots; j++ ){

      Rcpp::StringVector Vj = DF[(j+1)];
       
      for( i=0; i<nrows; i++ ){
	Rcpp::String v1(V1[i]);
	Rcpp::String vj(Vj[i]);	
	for( u=0; u<N; u++ ){
	  if( strcmp( gg->V[u].label, v1.get_cstring() )==0 ){
	    //add annotation data to node
	    gg->V[u].Annos[j] = vj.get_cstring();
	    break;
	  }
	}       
      }//i
    }//j
     
     //print meta-data
    gg->printModels();
    gg->printAnnoData();
     
  }
  
}

// [[Rcpp::export]]
void runDCSBM( Rcpp::IntegerVector K=1,
	       Rcpp::IntegerVector SEED=0,
	       Rcpp::IntegerVector SLOT=0,
	       Rcpp::IntegerVector ITS=0,
	       Rcpp::IntegerVector NCORES=0){
  
  int k,min,slot,seed,its,ncores;
      
  k      = K[0];
  seed   = SEED[0];
  its    = ITS[0];
  ncores = NCORES[0];
  
  if( gg!=0 && k!=0 ){

    min  = gg->getMinAnnosLength();  
    slot = SLOT[0];
    if( slot >= 0 && slot < min ){
      ////slot = min-1;

      ////reset model
      resetModel();
    
      if( min > 0 ){ model = new DCSBM( gg, k, slot ); }
      else         { model = new DCSBM( gg, k );      }
    

      try{
	//run model
	model->useSeed(seed);
	model->setOpenMP(ncores);
	model->runDCSBM();
	  //model->ompTEST();
	model->getVertexKProbs();	  
      } catch (Rcpp::internal::InterruptedException& e)
	{
	  cout << "Caught an interrupt!" << endl;
	}	
    }
    
  }
  
}

// [[Rcpp::export]]
Rcpp::List fitParams(){

   Rcpp::NumericVector ll(1);
   ll[0] = 0.0;

   Rcpp::NumericVector EMsteps(1);
   EMsteps[0] = 0.0;

   Rcpp::NumericVector maxDelta(1);
   maxDelta[0] = 0.0;
   
  if( gg !=0 && model !=0 ){

    ll[0]       = model->getMaxLL();    
    EMsteps[0]  = model->getEMsteps();
    maxDelta[0] = model->getMaxDelta();
    
  }

  return Rcpp::List::create(Rcpp::Named("ll") = ll,
			    Rcpp::Named("EMsteps") = EMsteps,
			    Rcpp::Named("maxDelta") = maxDelta);
}


// [[Rcpp::export]]
Rcpp::List getVertexKProbs(){

  int i,k,N,K,A;

  if( gg!=0 ){
  
    N = gg->getN();

    model->getVertexKProbs();
    
    //--- output the node label and its cluster
    Rcpp::NumericVector ID    (N);
    Rcpp::StringVector  Label (N);
    Rcpp::NumericVector C     (N);

    for( i=0; i<N; i++ ){
      ID[i]    = gg->V[i].id;
      Label[i] = gg->V[i].label;
      C[i]     = gg->V[i].K;
    }

    // Create a named data.frame with the above quantities
    Rcpp::DataFrame DF = Rcpp::DataFrame::create(
				 Rcpp::Named("id")   = ID,
				 Rcpp::Named("name") = Label,
				 Rcpp::Named("K")    = C);

    
    //--- output node community probability values
    K = gg->V[0].Nk;

    Rcpp::NumericMatrix Kprobs ( N,(K+1) );    
    for( i=0; i<N; i++ ){
      Kprobs(i,0) = gg->V[i].id;
      for( k=0; k<K; k++ ){    
	Kprobs(i,(k+1)) = gg->V[i].Kprobs[k];
      }
    }       						 

    A = gg->V[0].Na;
    Rcpp::StringMatrix  ANNO  (N,A);
    for( i=0; i<N; i++ ){
      for( k=0; k<A; k++ ){
	ANNO(i,k) = gg->V[i].Annos[k];
      }
    }    

    
    return Rcpp::List::create(Rcpp::Named("df")   = DF,
			      Rcpp::Named("np")   = Kprobs,
			      Rcpp::Named("anno") = ANNO);
    
  } else {
    
    //--- output the node label and its cluster
    Rcpp::NumericVector ID    (0);
    Rcpp::StringVector  Label (0);
    Rcpp::NumericVector C     (0);
    Rcpp::NumericMatrix Kprobs(0);
    Rcpp::StringMatrix  ANNO(0);
    
    // Create a named data.frame with the above quantities
    Rcpp::DataFrame DF = Rcpp::DataFrame::create(
				 Rcpp::Named("id")   = ID,
				 Rcpp::Named("name") = Label,
				 Rcpp::Named("K")    = C);           
    
    return Rcpp::List::create(Rcpp::Named("df")   = DF,
			      Rcpp::Named("np")   = Kprobs,
			      Rcpp::Named("anno") = ANNO);
  }
    
}



