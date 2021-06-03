#include "vertex.h"

vertex::vertex() : edge(){

  id      = 0;          
  degree  = 0;      
  K       = 1;           
  bridge  = 0;   
  label   = 0;   //safer to set '0' than NULL.  
  GeneName= "";
  EntrezID= "";
  E       = 0;       
  Nk      = 0;
  Na      = 0;
  PRINT   = false;
  Kprobs  = 0;
  Annos   = 0;

  inDeg   = 0;
  Indeg   = 0;

  outDeg  = 0;
  Outdeg  = 0;
  
}

vertex::~vertex(){ freeSpace(); }


void vertex::copy( vertex *VV ){

  int i,length;
  
  id     = VV->id;
  degree = VV->degree;
  K      = VV->K;  
  bridge = VV->bridge;
  
  //---copy vertex label
  length = strlen(VV->label); 
  label  = new char[length];
  strcpy(label,VV->label);

  //---copy vertex GeneName
  GeneName.assign(VV->GeneName);

  //---copy vertex Entrez gene ID
  EntrezID.assign(VV->EntrezID);

  
  assignE       ( degree );
  assignKprobs  ( VV->Nk );
  assignNannos  ( VV->Na );
  assignIndeg   ( VV->inDeg );
  assignOutdeg  ( VV->outDeg );
    
  //fill
  if( Nk != 0 ){
    for(i=0; i<Nk; i++)
      Kprobs[i] = VV->Kprobs[i];
  }  

  if( Na != 0 ){
    for(i=0; i<Na; i++)
      Annos[i] = VV->Annos[i];
  }

   if( inDeg != 0 ){
    for(i=0; i<inDeg; i++)
      Indeg[i] = VV->Indeg[i];
  }

   if( outDeg != 0 ){
    for(i=0; i<outDeg; i++)
      Outdeg[i] = VV->Outdeg[i];
  }
   
}

void vertex::assignKprobs( int slots ){
  
  int i;

  if( slots >= 0 ){

    freeKprobs();
    Nk = slots;
  
    Kprobs = new double[Nk];
  
    for(i=0; i<Nk; i++) Kprobs[i] = 0.0;

  }
    
}

void vertex::assignNannos( int slots ){
  
  int i;

  if( slots >= 0 ){

    freeNannos();
    Na = slots;
  
    Annos = new string[Na];
  
    for(i=0; i<Na; i++) Annos[i] = "";

  }
    
}

void vertex::assignAnno( int slot, string X ){
  
  if( (slot >= 0) && (slot < Na) ){
    Annos[slot] = X;
  }
    
}



void vertex::freeKprobs(){

  //if(Kprobs!=0 && Nk!=0) free(Kprobs);
  if(Kprobs!=0 && Nk!=0){ delete[] Kprobs; }

  Nk = 0;

}


void vertex::freeNannos(){

  //if(Kprobs!=0 && Nk!=0) free(Kprobs);
  if(Annos!=0 && Na!=0){ delete[] Annos; }

  Na = 0;

}

void vertex::assignE( int _Ek ){
  
  int i;

  if( _Ek >= 0 ){

    freeE();
    degree = _Ek;
    
    //E = (edge*)malloc(degree*sizeof(edge));
    E = new edge[degree];

  }
    
}

void vertex::assignIndeg( int in ){

  int i;

  if( in >= 0 ){

    freeIndeg();
    inDeg = in;
  
    Indeg = new int[inDeg];
  
    for(i=0; i<inDeg; i++) Indeg[i] = 0;

  }
    
  
}

void vertex::assignIndeg( vector<int> &in ){

  int i;

  if( in.size() >= 0 ){

    freeIndeg();
    inDeg = in.size();
  
    Indeg = new int[inDeg];
  
    for(i=0; i<inDeg; i++) Indeg[i] = in[i];

  }
    
  
}

void vertex::assignOutdeg( int out ){

  int i;

  if( out >= 0 ){

    freeOutdeg();
    outDeg = out;
  
    Outdeg = new int[outDeg];
  
    for(i=0; i<outDeg; i++) Outdeg[i] = 0;

  }
    
  
}


void vertex::assignOutdeg( vector<int> &out ){

  int i;

  if( out.size() >= 0 ){

    freeOutdeg();
    outDeg = out.size();
  
    Outdeg = new int[outDeg];
  
    for(i=0; i<outDeg; i++) Outdeg[i] = out[i];

  }
    
  
}


void vertex::freeE(){   

  //if(E!=0 && degree!=0) free(E);     //as using malloc in readgml.C
  if(E!=0 && degree!=0){ delete[] E; }

  degree = 0;
}


void vertex::freeIndeg(){

  if(Indeg!=0 && inDeg!=0){ delete[] Indeg; }

  inDeg = 0;

}

void vertex::freeOutdeg(){

  if(Outdeg!=0 && outDeg!=0){ delete[] Outdeg; }

  outDeg = 0;

}

void vertex::freeSpace(){

  freeE();
  
  //if(label!=0)  free(label); //as using malloc in readgml.C  
  if(label!=0){ delete[] label; } 

  freeNannos();
  
  freeKprobs();

  freeIndeg();

  freeOutdeg();

}

void vertex::printV(){

  
  if( PRINT ){ cout << "Vertex properties: " << endl; }
  if( PRINT ){ cout << "id: " << id << ", degree " << degree << ", K " << K << ", bridge " << bridge << ", label " << label << ", E " << E << " [" << degree << "] , Kprobs " << Kprobs << " [" << Nk << "] , Annos " << Annos << " [" << Na << "]" << endl; }
  
}

void vertex::setPrint( bool status ){
  PRINT = status;
}
