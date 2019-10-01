/*******************************************************************************
* FILE: mpi_heat2D.c
* DESCRIPTIONS:
*   HEAT2D Example - Parallelized C Version
*   This example is based on a simplified two-dimensional heat
*   equation domain decomposition.  The initial temperature is computed to be
*   high in the middle of the domain and zero at the boundaries.  The
*   boundaries are held at zero throughout the simulation.  During the
*   time-stepping, an array containing two domains is used; these domains
*   alternate between old data and new data.
*
*   In this parallelized version, the grid is decomposed by the master
*   process and then distributed by rows to the worker processes.  At each
*   time step, worker processes must exchange border data with neighbors,
*   because a grid point's current temperature depends upon it's previous
*   time step value plus the values of the neighboring grid points.  Upon
*   completion of all time steps, the worker processes return their results
*   to the master process.
*
*   Two data files are produced: an initial data set and a final data set.
* AUTHOR: Blaise Barney - adapted from D. Turner's serial C version. Converted
*   to MPI: George L. Gusciora (1/95)
* LAST REVISED: 04/02/05
*******************************************************************************/
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
  void inidat(), prtdat(), update();
  float  u[2][NXPROB][NYPROB];    /* array for grid */
  int	taskid;                     /* this task's unique id */
  int numworkers;                 /* number of worker processes */
	int numtasks;                   /* number of tasks */
	int ave_row,rows,offset_row,extra_row;   /* for sending rows of data */
  int ave_column,columns,
  offset_column,extra_column;     /* for sending columns of data*/
	int dest, source;               /* to - from for message send-receive */
	int left,right,up,down;         /* neighbor tasks */
  int row_start, row_end;         /* worker's row borders */
  int column_start, column_end;   /* worker's column borders */
	int msgtype;                    /* for message types */
	int rc;                         /* misc */
	int i,j,ix,iy,iz,it;            /* loop variables */
  double workers_root;            /* square root of workers to divide the grid*/
  MPI_Status status;
  MPI_Datatype MPI_ROW, MPI_COLUMN; /* datatypes used for efficient data transfers between workers */
  MPI_Comm MPI_CART_COMM;
  int cart_ndims = 2;
  int cart_dims[2] = {0, 0};
  int cart_periods[2] = {1, 1};
  int cart_reorder = 1;


  /* First, find out my taskid and how many tasks are running */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  numworkers = numtasks-1;

  /* Create a cartesian topology for the processes */
  MPI_Dims_create(numtasks, cart_ndims, cart_dims);
  if(taskid == 0){
    printf("Creating a cartesian topology of %d x %d dimensions", cart_dims[0], cart_dims[1]);
  }
  MPI_Cart_create(MPI_COMM_WORLD, cart_ndims, cart_dims, cart_periods, cart_reorder, &MPI_CART_COMM);
  MPI_Barrier(MPI_CART_COMM);

  /* Because of reordering the process id might have changed */
  MPI_Comm_rank(MPI_CART_COMM, &taskid);

  /* Get the name of the physical node the process is running in */
  char p_name[MPI_MAX_PROCESSOR_NAME];
  int p_name_len;
  MPI_Get_processor_name(p_name, &p_name_len);
  printf("Process #%d is running in processor: %s", taskid, p_name);

  MPI_Finalize(); //For debugging


  if (taskid == MASTER) {
    /******************************* master code ******************************/
    /* Check if numworkers is within range - quit if not */
    if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)) {
      printf("ERROR: the number of tasks must be between %d and %d.\n", MINWORKER+1, MAXWORKER+1);
      printf("Quitting...\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
      exit(1);
    }
    printf ("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

    /* Initialize grid */
    printf("Grid size: X= %d  Y= %d  Time steps= %d\n", NXPROB, NYPROB, STEPS);
    printf("Initializing grid and writing initial.dat file...\n");
    inidat(NXPROB, NYPROB, u);
    prtdat(NXPROB, NYPROB, u, "initial.dat");

    /* Distribute work to workers.  Must first figure out how many rows to
    *  send and what to do with extra rows. */

    workers_root = sqrt(numworkers);

    ave_row = NXPROB/(int)workers_root;
    extra_row = NXPROB%(int)workers_root;
    offset_row = 0;

    /*Same treatment for columns. Figure out how many columns to send
    * and what to do with extra columns. */
    ave_column = NYPROB/(int)numworkers;
    extra_column = NYPROB%(int)numworkers;
    offset_column = 0;

    for (i=1; i<=workers_root; i++) {
      rows = (i <= extra_row) ? ave_row+1 : ave_row; /*den eimai sigouros akoma an douleuei swsta twra auto*/
      for (j=1;j<=workers_root;j++)
      {
        /*The destination id is calculated using the number of full row sets, and column
        * sets we've assigned in this row. Given that, the worker's neighbors to the
        * left and right are +-1 respectively, and its up and down
        * neighbors are +-workers_root. */
        dest = (i-1)*workers_root + j;

        /* Tell each worker who its neighbors are, since they must exchange
        *  data with each other. */

          if (i == 1)
            up = NONE;
          else
            up = dest - workers_root;
          if (i == workers_root)
            down = NONE;
          else
            down = dest + workers_root;
          if (j == 1)
            left = NONE;
          else
            left = dest - 1;
          if (j == workers_root)
            right = NONE;
          else
            right = dest+1;

        /*apo edw kai katw ola idia - sigoura prepei na allaksei to offset apla den ithela na peiraksw mono ena pragma*/


        /* Now send startup information to each worker */
        MPI_Send(&offset_row, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&up, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&down, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&u[0][offset_row][0], rows*NYPROB, MPI_FLOAT, dest, BEGIN, MPI_COMM_WORLD);
        printf("Sent to task %d: rows= %d offset= %d ",dest, rows, offset_row);
        printf("up= %d down= %d\n", up, down);
        offset_row = offset_row + rows;
        offset_column = offset_column + rows;
      }
    }
    /* Now wait for results from all worker tasks */
    for (i=1; i<=numworkers; i++) {
      source = i;
      msgtype = DONE;
      // MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&u[0][offset_row][0], rows*NYPROB, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);
    }

    /* Write final output, call X graph and finalize MPI */
    printf("Writing final.dat file and generating graph...\n");
    prtdat(NXPROB, NYPROB, &u[0][0][0], "final.dat");
    printf("Click on MORE button to view initial/final states.\n");
    printf("Click on EXIT button to quit program.\n");

    MPI_Finalize();
  } /* End of master code */



  /******************************* workers code *******************************/
  if (taskid != MASTER) {
    /* Initialize everything - including the borders - to zero */
    for (iz=0; iz<2; iz++)
      for (ix=0; ix<NXPROB; ix++)
        for (iy=0; iy<NYPROB; iy++)
          u[iz][ix][iy] = 0.0;

    /* Receive my offsets, rows, columns, neighbors and grid partition from master */
    source = MASTER;
    msgtype = BEGIN;
    MPI_Recv(&offset_row, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&offset_column, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&columns, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&down, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&u[0][offset][0], rows*columns, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);

    /* Determine border elements.  Need to consider first and last columns.
    *  Obviously, row 0 can't exchange with row 0-1.  Likewise, the last
    *  row can't exchange with last+1. */
    if (offset_row == 0)
      row_start = 1;
    else
      row_start = offset_row;
    if ((offset_row+rows) == NXPROB)
      row_end = row_start+rows-2;
    else
      row_end = row_start+rows-1;

    /*  Τhe same goes for columns. */
    if (offset_column == 0)
      column_start = 1;
    else
      column_start = offset_column;
    if ((offset_column+columns) == NYPROB)
      column_end = column_start+columns-2;
    else
      column_end = column_start+columns-1;

    /* Create row and column datatypes. This way we can efficiently send a column
    * from the table without having to copy it to a buffer. More specifically:
    * A row has <column-size> blocks each of which contain 1 MPI_DOUBLE and there
    * is 1 block between the starts of each block in our table.
    * A column has <row-size> blocks each of which contain 1 MPI_DOUBLE and
    * there are <row-size> blocks between the start of each block in our table */
    MPI_Type_vector(columns, 1, 1, MPI_DOUBLE, &MPI_ROW);
    MPI_Type_vector(rows, 1, rows, MPI_DOUBLE , &MPI_COLUMN);
    MPI_Type_commit(&MPI_ROW);
    MPI_Type_commit(&MPI_COLUMN);

    /* Begin doing STEPS iterations.  Must communicate border rows with
    *  neighbors.  If I have the first or last grid row, then I only need
    *  to  communicate with one neighbor */
    printf("Task %d received work. Beginning time steps...\n",taskid);
    iz = 0;
    for (it = 1; it <= STEPS; it++) {
      if (left != NONE) {
        MPI_Send(&u[iz][offset_row][offset_column], 1, MPI_COLUMN, left, RTAG, MPI_COMM_WORLD);
        MPI_Recv(&u[iz][offset_row][offset_column-1], 1, MPI_COLUMN, left, LTAG, MPI_COMM_WORLD, &status);
      }
      if (up != NONE) {
        MPI_Send(&u[iz][offset_row][offset_column], 1, MPI_ROW, up, DTAG, MPI_COMM_WORLD);
        MPI_Recv(&u[iz][offset_row-1][offset_column], 1, MPI_ROW, up, UTAG, MPI_COMM_WORLD, &status);
      }
      if (right != NONE) {
        MPI_Send(&u[iz][offset_row][offset_column+columns-1], 1, MPI_COLUMN, right, LTAG, MPI_COMM_WORLD);
        MPI_Recv(&u[iz][offset_row][offset_column+columns], 1, MPI_COLUMN, right, RTAG, MPI_COMM_WORLD, &status);
      }
      if (down != NONE) {
        MPI_Send(&u[iz][offset_row+rows-1][offset_column], 1, MPI_ROW, down, UTAG, MPI_COMM_WORLD);
        MPI_Recv(&u[iz][offset_row+rows][offset_column], 1, MPI_ROW, down, DTAG, MPI_COMM_WORLD, &status);
      }
      /* Now call update to update the value of grid points */
      update(row_start, row_end, column_start, column_end, NYPROB, &u[iz][0][0],&u[1-iz][0][0]);
      iz = 1 - iz;
    }

    /* Finally, send my portion of final results back to master */
    MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
    MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
    MPI_Send(&u[iz][offset][0], rows*columns, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
    MPI_Finalize();
   }
}

/****************************** subroutine update *****************************/
/* u2[ix][iy] = u1[ix][iy]
                + cx * ( u1[ix+1][iy] + u1[ix-1][iy] - 2 * u1[ix][iy] )
                + cy * ( u1[ix][iy+1] + u1[ix][iy-1] - 2 * u1[xi][yi] ) */
void update(int x_start, int x_end, int y_start, int y_end,int ny, float *u1, float *u2) {
  int ix, iy;
  for (ix = x_start; ix <= x_end; ix++)
    for (iy = y_start; iy <= y_end; iy++)
      *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  +
                      parms.cx * (*(u1+(ix+1)*ny+iy) +
                      *(u1+(ix-1)*ny+iy) -
                      2.0 * *(u1+ix*ny+iy)) +
                      parms.cy * (*(u1+ix*ny+iy+1) +
                     *(u1+ix*ny+iy-1) -
                      2.0 * *(u1+ix*ny+iy));
}

/****************************** subroutine inidat *****************************/
void inidat(int nx, int ny, float *u) {
  int ix, iy;

  for (ix = 0; ix <= nx-1; ix++)
    for (iy = 0; iy <= ny-1; iy++)
      *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}

/****************************** subroutine prtdat *****************************/
void prtdat(int nx, int ny, float *u1, char *fnam) {
  int ix, iy;
  FILE *fp;

  fp = fopen(fnam, "w");
  for (iy = ny-1; iy >= 0; iy--) {
    for (ix = 0; ix <= nx-1; ix++) {
      fprintf(fp, "%6.1f", *(u1+ix*ny+iy));
      if (ix != nx-1)
        fprintf(fp, " ");
      else
        fprintf(fp, "\n");
    }
  }
  fclose(fp);
}
