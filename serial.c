/*
 Game of life
 initial pattern is a glider
 Compile:  gcc -o gameSerial serial.c
 Run:   ./gameSerial
*/
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#define k 27
#define position(row,col) printf("%c[%d;%dH",k,row,col)
#define clear() printf("%c[2J",k)
#define ROWS 28

char  *TOP[ROWS];
char  *BOTTOM[ROWS];
char  *FULL[ROWS] = {
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                     O                                            ",
  "                                      O                                           ",
  "                                    OOO                                           ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
        "                                                                                  "
};

void game( char**, char** );
void initialize( );
void display( char ** );


// --------------------------------------------------------------------------
// initialize
// initialize top and bottom of grid
void initialize( ) {
  int x;
  //allocate memory and copy input grid to top and bottom for processing
  for (x = 0; x< ROWS; x++ )  {
    TOP[x] = (char *) malloc( (strlen( FULL[0] ) + 1 ) * sizeof( char ) );
    strcpy( TOP[x], FULL[x] );
    BOTTOM[x] = (char *) malloc( (strlen( TOP[0] )+1) * sizeof( char )  );
    strcpy( BOTTOM[x], TOP[x] );
  }
}

// Display current iteration
void display( char* dish[] ) {
  int j;
  for (j=0; j<ROWS; j++ ) {
    position( j, 0 );
    printf( "%s \n", dish[j] );
  }

}

void  game( char** part, char** nextStep ) {
  /*
   * When given a portion of the grid, this function will calculate the 
   * next step and return the new portion of the grid
   */
  int m, n, row;
  int rowSize = strlen( part[0] );
  int partSize = ROWS;

  for (row = 0; row < ROWS; row++) {            // loop through the rows in grid
    for ( m = 0; m < rowSize; m++) {            // loop through each character in the row
      int o, n, side = 0;
      char curr = part[row][m];
      // search a 3 by 3 area in the grid for any live cells 'O'
      for ( o = row - 1; o <= row + 1; o++) {
	int realo = o;// traverse from bottom to top
	if (o == -1)
	  realo = partSize - 1;
	if (o == partSize)
	  realo = 0;
	for ( n = m - 1; n <= m + 1; n++) {
	  int realn = n; // traverse left to right
	  if (n == -1)
	    realn = rowSize - 1;
	  if (n == rowSize)
	    realn = 0;
	  if (o == row && n == m)
	    continue;
	  if (part[realo][realn] == 'O')
	    side++;
	}
      }
      if (curr == 'O') {
	if (side < 2 || side > 3)
	  nextStep[row][m] =  ' ';
	else
	  nextStep[row][m] = 'O';
      }
      if (curr == ' ') {
	if (side == 3)
	  nextStep[row][m] = 'O';
	else
	  nextStep[row][m] = ' ';
      }
    }
  }
}

//function to get wall time
double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


int main( int argc, char* argv[] ) {
  double start,end;
  start = get_wall_time();
  int steps = 1000;      // # of steps in game
  int u;
  char **part, **nextPart, **tmp;
  //empty screen for next step
  clear();
  initialize();
  part   = TOP;
  nextPart = BOTTOM;
  display( part );       // show initial step
  for ( u = 0; u < steps; u++) { // show rest of the steps
    // follow rules in game of life to current step and display next step
    game( part, nextPart );
    display( part );
    // This is just for user to be able to see whats going on. This must be removed for timing calculations
    sleep(1); 
    tmp = part; // obtain next step
    part = nextPart;
    nextPart = tmp;
  }
  display(part); // show result ie. final step
  
  end = get_wall_time();
  //printf("Total time elasped is %lf \n",(end-start));
}
