/*
Program Name: Parallel program that implements Conway's Game of Life using MPI
Created:  Mar 28, 2017
Modified: Apr 2, 20172:00 PM

Compile: mpicc mpi.c -o test.x -lm
Run    : mpirun -np 4 test.x
*/

#include "mpi.h"
#include <stdio.h>
#include <math.h>

#define GENERATIONS 1000
#define M 82

/* global variables */
int p, my_rank, sqr, adjust, my_grid_rank, grid_rank, m_cell, n_cell;
int dims[2], periods[2], coords[2], remain_dims[2];
int stateNum = 0, reorder = 1, ndims = 2, maxdims = 2;
MPI_Comm  m_comm, n_comm, grid_comm;

/* fucntion prototypes */
void init_grid(int grid[M][M]);
void init_glider(int grid[M][M]);
void life(int intial_state[adjust][adjust]);
void reconstruct_grid (int intial_state[adjust][adjust], int temp[adjust][adjust]);
void create_temp(int intial_state[adjust][adjust], int temp[adjust][adjust]);
void writeStep(int intial_state[adjust][adjust]);


int main(int argc, char* argv[]) {
  int i, j;
  int grid[M][M];
  double slavetime, totaltime, starttime;
  starttime = MPI_Wtime();
  
  MPI_Init(&argc, &argv);

  //  starttime = MPI_Wtime();

  /* Determine size and my rank in communicator */
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  sqr = (int) sqrt((double) p);
  dims[0] = sqr; /* dimension at 0 and 1 must equal sqrt */
  dims[1] = sqr;
  periods[0] = 1; /* enable coordinate wrapping */
  periods[1] = 1;

  /*  Create Cartesian topology  */
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &grid_comm);

  MPI_Comm_rank(grid_comm, &my_grid_rank);

  MPI_Cart_coords(grid_comm, my_grid_rank, maxdims, coords);
  MPI_Cart_rank(grid_comm, coords, &grid_rank);

  /*  split grid verticall   */
  remain_dims[0] = 0;
  remain_dims[1] = 1;

  MPI_Cart_sub(grid_comm, remain_dims, &m_comm);

  if (coords[1] == 0) {m_cell = coords[0];}
  else { m_cell = -1;}

  /* Broadcast the size to all processes */
  MPI_Bcast(&m_cell, 1, MPI_INT, 0, m_comm);

  /*  split grid horizontally   */
  remain_dims[0] = 1;
  remain_dims[1] = 0;

  MPI_Cart_sub(grid_comm, remain_dims, &n_comm);

  if (coords[0] == 0) {n_cell = coords[1];}
  else {n_cell = -1;}

  /* Broadcast the size to all processes */
  MPI_Bcast(&n_cell, 1, MPI_INT, 0, n_comm);

  adjust = M / sqr;
  int my_grid[adjust][adjust];

  /* Initialize graph and glider */
  init_grid(grid);
  init_glider(grid);

  /* Determine my part of the grid */
  for (i = 0; i < adjust; i++) {
    for (j = 0; j < adjust; j++) {
      my_grid[i][j] = grid[i + m_cell * adjust][j + n_cell * adjust];
    }
  }

  writeStep(my_grid);

  /* Run the game for given number of generations */
  for (i = 0; i < GENERATIONS; i++) {
    /* Call the life routine */
    life(my_grid);
    /* Call this routine to create output files */
    writeStep(my_grid);
  }

  /* Print the average time taken/processor */
  slavetime = MPI_Wtime() - starttime;
  MPI_Reduce (&slavetime, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, grid_comm);
  if (my_rank == 0)
    printf("Average time taken/processor: %lf s\n",totaltime/(double)p);

  MPI_Finalize();
  return 0;
}

/* The Life routine */
void life(int intial_state[adjust][adjust]) {
  int i, j;
  int ghost[adjust+2][adjust+2];    //ghost graph
  int temp[adjust][adjust];         //temp graph
  int left, right, above, below;           // neighbors
  int rankUpperLeft, upperRank, upperRightRank, lefRank, rigRank, lowerLeftRank, lowerRank, lowerRightRank;
  int upperLeft, upperRight, lowerLeft, lowerRight, firstRow[adjust],lastRow[adjust], firstColumn[adjust],
    lastColumn[adjust],shiftL[adjust],shiftR[adjust];

  /* Set neighbors */
  below = m_cell + 1;above = m_cell - 1;
  right = n_cell + 1; left = n_cell - 1;

  /* Initialize the boundaries of the life matrix */
  if (above < 0) above = sqr - 1; if (below == sqr) below = 0;
  if (left < 0)  left = sqr - 1; if (right == sqr) right = 0;

  /* Update local ranks */
  upperRank = above * sqr + n_cell;
  lowerLeftRank = below * sqr + left;
  lowerRank = below * sqr + n_cell;
  lowerRightRank = below * sqr + right;
  lefRank = m_cell * sqr + left;
  rankUpperLeft = above * sqr + left;
  upperRightRank = above * sqr + right;
  rigRank = m_cell * sqr + right;

  /* initialize ghost points from intial graph*/
  for (i = 0; i < adjust+2; i++) {
    for (j = 0; j < adjust+2; j++) {
      ghost[i][j] = 0;
    }
  }

  for (i = 0; i < adjust; i++) {
    for (j = 0; j < adjust; j++) {
      ghost[i+1][j+1] = intial_state[i][j];
    }
  }

  /* Swap the grids */
  for  (i = 0; i < adjust; i++){
    shiftL[i] = intial_state[i][0];
    shiftR[i] = intial_state[i][adjust-1];
  }

  /* Send and receive boundary information */
  MPI_Sendrecv(&intial_state[0][0], 1, MPI_INT, rankUpperLeft, 1, &lowerRight, 1, MPI_INT, lowerRightRank, 1, grid_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&intial_state[adjust-1][adjust-1], 1, MPI_INT, lowerRightRank, 1, &upperLeft, 1, MPI_INT, rankUpperLeft, 1, grid_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&intial_state[0][adjust-1], 1, MPI_INT, upperRightRank, 1, &lowerLeft, 1, MPI_INT, lowerLeftRank, 1, grid_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&intial_state[adjust-1][0], 1, MPI_INT, lowerLeftRank, 1, &upperRight, 1, MPI_INT, upperRightRank, 1, grid_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&intial_state[0], adjust, MPI_INT, upperRank, 1, &lastRow[0], adjust, MPI_INT, lowerRank, 1, grid_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&intial_state[adjust-1], adjust, MPI_INT, lowerRank, 1, &firstRow[0], adjust, MPI_INT, upperRank, 1, grid_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&shiftL[0], adjust, MPI_INT, lefRank, 1, &lastColumn[0], adjust, MPI_INT, rigRank, 1, grid_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&shiftR[0], adjust, MPI_INT, rigRank, 1, &firstColumn[0], adjust, MPI_INT, lefRank, 1, grid_comm, MPI_STATUS_IGNORE);

  /* ghost graph */
  ghost[0][0] = upperLeft;
  ghost[0][adjust+1] = upperRight;
  ghost[adjust+1][0] = lowerLeft;
  ghost[adjust+1][adjust+1] = lowerRight;

  for (i = 0; i < adjust; i++){
    ghost[0][i+1] = firstRow[i];
    ghost[adjust+1][i+1] = lastRow[i];
    ghost[i+1][0] = firstColumn[i];
    ghost[i+1][adjust+1] = lastColumn[i];
  }

  /* create temporary graph from initial graph  */
  create_temp(intial_state, temp);

  /* Update live/dead points */
  for (i = 0; i < adjust; i++){
    for(j = 0; j < adjust; j++) {
      int liveNeighbors = 0;

      /* check neighbors */
      if (ghost[i][j] == 1) liveNeighbors++;if (ghost[i][j+1] == 1) liveNeighbors++;
      if (ghost[i][j+2] == 1) liveNeighbors++; if (ghost[i+1][j] == 1) liveNeighbors++;

      /* check current point */
      if (ghost[i+1][j+2] == 1) liveNeighbors++;if (ghost[i+2][j] == 1) liveNeighbors++;
      if (ghost[i+2][j+1] == 1) liveNeighbors++;if (ghost[i+2][j+2] == 1) liveNeighbors++;

      /* check if the cell dies or life is born or no change */
      if (liveNeighbors == 3) {temp[i][j] = 1; } // comes to life
      else if (liveNeighbors == 2) { }// does not change change
      else if (liveNeighbors < 2)  { temp[i][j] = 0; }// dies a horrible death
      else { temp[i][j] = 0; } // dies a horrible death
    }
  }
  /* reconstruct the grid */
  reconstruct_grid (intial_state, temp);
}

/* writes the state of the whole distributed array to a single file */
void writeStep(int intial_state[adjust][adjust]){
  int i, j, k, l, z;
  int new_graph[M][M];
  char filename[100];
  FILE *fp;
  
  /* Open the output file */
  if (my_rank == 0) {
    sprintf(filename, "state%d.txt", stateNum);
    printf("%s",filename);
    printf("\n");
    fp = fopen(filename, "w+");
  }

  if (my_rank == 0){
    int result[p][adjust][adjust];
    for (i = 1; i < p; i++) {
      MPI_Recv(&result[i], adjust*adjust, MPI_INT, i, 1, grid_comm, MPI_STATUS_IGNORE);
    }

    for (i = 0; i < adjust; i++){
      for (j = 0; j < adjust; j++){
	result[0][i][j] = intial_state[i][j];
      }
    }

    for (i = 0; i < sqr; i++){
      for (j = 0; j < sqr; j++){
	for (k = 0; k < adjust; k++){
	  for (l = 0; l < adjust; l++){
	    new_graph[k + i * adjust][l + j * adjust] = result[i * sqr + j][k][l];
	  }
	}
      }
    }

  } else {
    MPI_Send(&intial_state[0], adjust*adjust, MPI_INT, 0, 1, grid_comm);
  }

  if (my_rank == 0) {
    for (i = 0; i < M; i++){
      for(j = 0; j < M; j++) {
	if (new_graph[i][j] == 1) {
	  fprintf(fp, "O");
	} else {
	  fprintf(fp," ");
	}
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
    stateNum++;
  }
}

/* initialize grid routine */
void init_grid(int grid[M][M]){
  int i, j;
  for (i = 0; i < M; i++){
    for(j = 0; j < M; j++) {
      grid[i][j] = 0;
    }
  }
}

/* intitialize glider routine */
void init_glider(int grid[M][M]) {
  grid[11][12] = 1;
  grid[12][13] = 1;
  grid[13][11] = 1;
  grid[13][12] = 1;
  grid[13][13] = 1;
}

/* reconstruct grid routine */
void reconstruct_grid (int intial_state[adjust][adjust], int temp[adjust][adjust]) {
  for (int i = 0; i < adjust; i++){
    for(int j = 0; j < adjust; j++) {
      intial_state[i][j] = temp[i][j];
    }
  }
}

void create_temp (int intial_state[adjust][adjust], int temp[adjust][adjust]) {
  for (int i = 0; i < adjust; i++){
    for(int j = 0; j < adjust; j++) {
      temp[i][j] = intial_state[i][j];
    }
  }
}
