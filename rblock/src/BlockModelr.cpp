#include "Headers.h"
#include "readfile.h"
#include "network.h"
#include "DCSBM.h"

//----------------
//NOTE:
// When editing this file, i.e. adding or removing functions, need to carry out the following tasks:
// 1) edit file NAMESPACE, adding/removing the function name in export(...)
// Next run following in R, inside the package:
// 2) $> cd BlockModelr/
// 3) $> R
// 4)  > library(Rcpp)
// 5)  > Rcpp::compileAttributes()
// running 5) will create new RcppExports.R in R/ and
// new RcppExports.cpp in src/
//----------------

// Global
network *gg=0;
DCSBM   *model=0;
//unsigned long int seed;
//gsl_rng *g;

/*
// [[Rcpp::export]]
void setSeed(){
  g          = gsl_rng_alloc(gsl_rng_taus2);
  seed       = 0;
  gsl_rng_set(g , seed);  
}

// [[Rcpp::export]]
Rcpp::NumericVector rndTest(){
  Rcpp::NumericVector res(1);
  double val = gsl_rng_uniform(g);
  //cout << "  random no: " << val << endl;
  res[0] = val;
  return res;
}
*/

// [[Rcpp::export]]
void resetGraph(){ gg = 0; }

// [[Rcpp::export]]
void resetModel(){ model = 0; }

// [[Rcpp::export]]
void deleteGraph(){ if(gg != 0){ delete gg; }  }

// [[Rcpp::export]]
void deleteModel(){ if(model != 0){ delete model; }  }

// [[Rcpp::export]]
void loadGraph( Rcpp::DataFrame     DF,
		Rcpp::IntegerVector names=1){

  int i,j,k,KK;

  bool useLoops    = false;
  bool checkM      = true;
  int alphaNumeric = 1;

  //initialise network
  resetGraph();
  readfile *reader = 0;
  gg               = new network();
  
  int ncols = DF.length();
  int nrows = DF.nrows();

  if( (ncols > 0) && (nrows > 0) ){

    //set size for our DATASET
    KK              = nrows*ncols;
    string *DATASET = new string[KK];
        

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
void loadMetaData( Rcpp::DataFrame DF ){

  int i,j,u,N; 
  
  int ncols = DF.length();
  int nrows = DF.nrows();
  int slots = ncols-1;
  
  if( (ncols > 0) && (nrows > 0) && (gg!=0) ){

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
     gg->printAnnoData();
     
  }
  
}

// [[Rcpp::export]]
void runDCSBM( Rcpp::IntegerVector K=0 ){

  int k,min,slot;
      
  k = K[0];
  
  if( gg!=0 && k!=0 ){

    min  = gg->getMinAnnosLength();  
    slot = min-1;

    //reset model
    resetModel();
    
    if( min > 0 ){ model = new DCSBM( gg, k, slot ); }
    else         { model = new DCSBM( gg, k );      }
    

    model->run();
    model->getVertexKProbs();

    //gg->setPrint( true );
    //gg->printVertices();
    
  }
  
}

// [[Rcpp::export]]
Rcpp::List getVertexKProbs(){

  int i,N;

  if( gg!=0 ){
  
    N = gg->getN();

    model->getVertexKProbs();
    
    //--- output the node label and its cluster
    Rcpp::NumericVector ID    (N);
    Rcpp::StringVector  Label (N);
    Rcpp::NumericVector K     (N);

    for( i=0; i<N; i++ ){
      ID[i]    = gg->V[i].id;
      Label[i] = gg->V[i].label;
      K[i]     = gg->V[i].K;
    }

    // Create a named list with the above quantities
    return Rcpp::List::create(Rcpp::Named("id")   = ID,
			      Rcpp::Named("name") = Label,
			      Rcpp::Named("K")    = K);
    
  } else {
    
    //--- output the node label and its cluster
    Rcpp::NumericVector ID    (0);
    Rcpp::StringVector  Label (0);
    Rcpp::NumericVector K     (0);
    
    // Create a named list with the above quantities
    return Rcpp::List::create(Rcpp::Named("id")   = ID,
			      Rcpp::Named("name") = Label,
			      Rcpp::Named("K")    = K);
  }
    
}
