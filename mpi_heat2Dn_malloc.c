#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NXPROB      20                 /* x dimension of problem grid */
#define NYPROB      20                 /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */
#define MAXWORKER   8                  /* maximum number of worker tasks */
#define MINWORKER   3                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define UTAG        4                  /* message tag */
#define DTAG        5                  /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */

struct Parms {
  float cx;
  float cy;
} parms = {0.1, 0.1};

int main (int argc, char *argv[]) {
  void inidat(),prtdat(),prtfdat(),update(),inidat2();
  float **u[2];                           /* array for grid */
  float **final_grid;                     /* the final grid that gets printed*/
  float *temp;                            /* for temporarily storing malloc-ed memory */
  int	taskid;                             /* this task's unique id */
	int numtasks;                           /* number of tasks */
	int ave_row,rows,extra_row;             /* for sending rows of data */
  int ave_column,columns,extra_column;    /* for sending columns of data*/
	int dest, source;                       /* to - from for message send-receive */
	int left,right,up,down;                 /* neighbor tasks */
	int msgtype;                            /* for message types */
	int rc;                                 /* misc */
	int i,j,z,ix,iy,iz,it;                  /* loop variables */
  MPI_Status status;
  MPI_Datatype MPI_ROW, MPI_COLUMN;       /* datatypes used for efficient data transfers between workers */
  MPI_Comm MPI_CART_COMM;                 /* Cartesian Communication World */
  int cart_ndims = 2;                     /* Number of dimensions of cartesian grid */
  int cart_dims[2] = {0, 0};              /* Size of each dimension in the cartesian grid */
  int cart_periods[2] = {0, 0};           /* Period of each dimension in the cartesian grid */
  int cart_reorder = 1;                   /* Node reorder option during cartesian grid construction */
  int coord[2];                           /* Process coordinates in the cartesian grid */
  char p_name[MPI_MAX_PROCESSOR_NAME];    /* Name of the processor the process is running on */
  int p_name_len;                         /* Processor name length */
  MPI_Request r_array[4];                 /* Handles for receiving information */
  MPI_Request s_array[4];                 /* Handles for sending information */


  /* First, find out my taskid and how many tasks are running */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  rows = 2;
  columns = 4;

  for(z=0; z<2; z++){
    u[z] = (float**) malloc((rows+2) * sizeof(float*));
    temp = (float*) malloc((rows+2) * (columns+2) * sizeof(float));
    for(i=0; i<rows+2; i++){
      u[z][i] = (float*) &(temp[i * (columns + 2)]);
    }
  }
  printf("Memory allocation finished in process #%d\n", taskid);

  /* Initialize everything - including the borders - to zero */
  for (iz=0; iz<2; iz++)
    for (ix=0; ix<rows+2; ix++)
      for (iy=0; iy<columns+2; iy++)
        u[iz][ix][iy] = 0.0;

  prtfdat(rows+2, columns+2, u[0]);

  MPI_Finalize(); //for debugging
  return 0;

  //inidat(rows+2, columns+2, u[0]);

  //prtfdat(rows+2, columns+2, u[0]);
}

/****************************** subroutine inidat *****************************/
void inidat(int nx, int ny, float *u) {
  int ix, iy;

  for (ix = 0; ix <= nx-1; ix++)
    for (iy = 0; iy <= ny-1; iy++)
      *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}
/****************************** subroutine prtfdat *****************************/
void prtfdat(int nx, int ny, float *u1) {
  int ix, iy;

  for (iy = 0; iy <= ny-1; iy++) {
    for (ix = 0; ix <= nx-1; ix++) {
      printf("%6.1f", *(u1+ix*ny+iy));
      if (ix != nx-1)
        printf(" ");
      else
        printf("\n");
    }
  }
}
