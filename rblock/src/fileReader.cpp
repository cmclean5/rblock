#include "fileReader.h"

fileReader::fileReader(){

  this->stream  = 0;
  
  this->first   = 0;
  this->current = 0;
  
  this->Nlines  = 0;

  this->DATASET = 0;

  this->dels[0]='\t';
  //this->dels[1]=',';
  //this->dels[2]=' ';

  this->HEADdel[0]='#';
  
}

fileReader::fileReader(FILE *stream){

  this->stream = stream;

  if(this->stream == 0){
    perror("Error opeing file");
  }

  this->first   = 0;
  this->current = 0;

  this->Nlines  = 0;

  this->DATASET = 0;

  this->dels[0]='\t';
  //this->dels[1]=',';
  //this->dels[2]=' ';
  
  this->HEADdel[0]='#';
  
}

fileReader::~fileReader(){ free_buffer(); }


//Skip lines with header. Return 2 if header symbol found, Return 1
// if an EOF encountered, Return 0 if 
int fileReader::read_line(char line[LINELENGTH], char HEADER)
{
  if (fgets(line,LINELENGTH,stream)==NULL) return 1;

  //skip lines starting with '#' character
  if(line[0]==HEADER){ cout << "\n skipping: " << line << endl; return 2; }
  
  return 0;
  
}


// Function to read one line from a specified stream.  Return value is
// 1 if an EOF was encountered.  Otherwise 0.
int fileReader::read_line(char line[LINELENGTH])
{
  if (fgets(line,LINELENGTH,stream)==NULL) return 1;

  //skip lines starting with '#' character
  if(line[0]==HEADdel[0]){ cout << "\n skipping: " << line << endl; return 2; }
  
  line[strlen(line)-1] = '\0';   // Erase the terminating NEWLINE
  return 0;
}


// Function to read in the whole file into a linked-list buffer, so that we
// can do several passes on it, as required to read the GML format
// efficiently
int fileReader::fill_buffer()
{
 int length, fl;
  char line[LINELENGTH];
  LINE *previous;

  Nlines = 0;

  fl = read_line(line);          // first line
  
  if( fl==1 ){
    first = NULL;                // Indicates empty buffer
    return 1;
  }

  // Read any headers
  if(fl==2){
    fl = read_line(line);
    while( fl!=0 ){
      if( fl==1 ){
	first = NULL;
	return 1;
      } else {
	fl = read_line(line);
      }
    } 
  }

  length = strlen(line) + 1;
  first = (struct LINE*)malloc(sizeof(struct LINE));
  first->str = (char*)malloc(length*sizeof(char));
  strcpy(first->str,line);
  Nlines++;

  previous = first;
  while (read_line(line)==0) {
    length = strlen(line) + 1;
    previous->ptr = (struct LINE*)malloc(sizeof(struct LINE));
    previous = previous->ptr;
    previous->str = (char*)malloc(length*sizeof(char));
    strcpy(previous->str,line);
    Nlines++;
  }
  previous->ptr = NULL;          // Indicates last line

  return 0;
}


// Function to free up the buffer again
void fileReader::free_buffer()
{
  int i,j;
  LINE *thisptr;
  LINE *nextptr;

  thisptr = first;
  while (thisptr!=NULL) {
    nextptr = thisptr->ptr;
    free(thisptr->str);
    free(thisptr);
    thisptr = nextptr;
  }

  if(stream!=0){ fclose(stream); }

 
  if(DATASET!=0){ delete[] DATASET; }

}


// Function to reset to the start of the buffer again
void fileReader::reset_buffer()
{
  current = first;
}

// Function to get the next line in the buffer.  Returns 0 if there was
// a line or 1 if we've reached the end of the buffer.
int fileReader::next_line(char line[LINELENGTH])
{
  if (current==NULL) return 1;
  strcpy(line,current->str);
  current = current->ptr;
  return 0;
}


int fileReader::get_nlines(){
  return Nlines;
}

int fileReader::get_ncols(){
  return Ncols;
}


string* fileReader::readFile( int NCOLS ){

  int i,j,k,KK,start,end,length,row,col;
  vector<int>    del_pos;
  vector<string> tokens;
  bool           foundDEL;
  
  char line[LINELENGTH];
  char token[LINELENGTH];
  char *nonspace;
  
  if( Nlines == 0 ){ fill_buffer(); }
  
  if( NCOLS  == 0 ){//find number of delemiters

    //reset buffer
    reset_buffer();
    foundDEL = false;
    while( next_line(line) == 0 && foundDEL == false ){

      //skip lines starting with #
      if(line[0]==HEADdel[0]) continue;

      start  = 0;
      end    = 0;
      length = strlen(line);
      del_pos.clear();
      
      findDelimiterPos(line, length, dels, DELSIZE, del_pos );

      if( del_pos.size() != 0 ){ NCOLS = (del_pos.size() + 1); foundDEL = true; }

    }    
    
    Ncols = NCOLS;

  }

  if( NCOLS != 0 && Nlines != 0 ){

    KK      = Nlines*NCOLS;
    DATASET = new string[KK];

  
    //reset buffer 
    row=0; col=0; length=0; del_pos.clear();
    reset_buffer();

    //read each line in file buffer
    while( next_line(line) == 0 ){

      //skip lines starting with #
      if(line[0]==HEADdel[0]) continue;

      start  = 0;
      end    = 0;
      length = strlen(line);
    
      findDelimiterPos(line, length, dels, DELSIZE, del_pos );

      //if no delemiters skip
      if(del_pos.size() != 0){

	for(k=0; k<del_pos.size(); k++){

	  end=del_pos[k];
	
	  subString(line,token,(start+1),end);	  
	  for(nonspace=token; *nonspace==' '; nonspace++);
	  DATASET[(row*NCOLS)+col] = nonspace;
	  col++;
	  start=end+1;//+1 to miss out the delimiter
	}
      
	if( start < length ){
	  subString(line,token,(start+1),length);
	  for(nonspace=token; *nonspace==' '; nonspace++);
	  DATASET[(row*NCOLS)+col] = nonspace;
	}
      
      }   
    
      tokens.clear();
      del_pos.clear();
      col=0;
      row++;

    }

  } else {
    cout << "> could not find delemiters or records in file" << endl;
  }
  
    
  vector<int>().swap(del_pos);  //free space of vector
  vector<string>().swap(tokens);//free space of vector
  

  return DATASET;

}
