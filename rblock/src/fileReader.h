//define guards, so headers are declare only once.
#ifndef FILEREADER_H
#define FILEREADER_H

#include "Headers.h"

class fileReader {

 public:
  fileReader();
  fileReader(FILE *);
 ~fileReader();

 int  read_line   (char * , char );
 int  read_line   (char *);
 int  fill_buffer ();
 void free_buffer ();
 void reset_buffer();
 int  next_line   (char *);
 int  get_nlines  ();
 int  get_ncols   ();
 string* readFile ( int=0 );
 
 
static inline void subString(char s[], char sub[], int p, int l){
   int c    = 0;
   int diff = (l-p); 
 
   while (c <= diff) {
     sub[c] = s[p+c-1];
     c++;
   }
   sub[c] = '\0';
 };
 

static inline void findDelimiterPos( char Line[], int Length, char dels[], int Ndels, vector<int> &pos){ 

  int i,d;
  
  //loop over each character in Line
  for(i=0; i<Length; i++){
      
    //search position of delemiters in line
    for(d=0; d<Ndels; d++){
      if(Line[i]==dels[d]){ pos.push_back(i); }
    }

  }

};


 private:
 
 int Nlines;
 int Ncols;

 FILE *stream;

 //char ***DATASET;
 string *DATASET;
 
 // Constants
 #define LINELENGTH 1000

 // Types
 struct LINE {
   char *str;
   struct LINE *ptr;
 };

 // Globals
 LINE *first;
 LINE *current;

 //set file delineator(s)
 static const int DELSIZE = 1;
 char dels[DELSIZE];

 //set file header delineator(s)
 static const int HEADDELSIZE = 1;
 char HEADdel[HEADDELSIZE];
  
  
};

#endif
