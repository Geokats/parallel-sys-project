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

  if (taskid == 0) {
    /* Check if numtasks is within range - quit if not */
    if ((numtasks > MAXWORKER) || (numtasks < MINWORKER)) {
      printf("ERROR: the number of tasks must be between %d and %d.\n", MINWORKER+1, MAXWORKER+1);
      printf("Quitting...\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
      exit(1);
    }
    printf ("Starting mpi_heat2D with %d worker tasks.\n", numtasks);

    printf("Grid size: X= %d  Y= %d  Time steps= %d\n", NXPROB, NYPROB, STEPS);
    printf("Initializing grid and writing initial.dat file...\n");

  }

  /* Create a cartesian topology for the processes */
  MPI_Dims_create(numtasks, cart_ndims, cart_dims);
  if(taskid == 0){
    printf("Creating a cartesian topology of %d x %d dimensions\n", cart_dims[0], cart_dims[1]);
    printf("MPI_PROC_NULL == %d\n", MPI_PROC_NULL);
  }
  MPI_Cart_create(MPI_COMM_WORLD, cart_ndims, cart_dims, cart_periods, cart_reorder, &MPI_CART_COMM);
  MPI_Barrier(MPI_CART_COMM);

  /* Because of reordering the process id might have changed */
  MPI_Comm_rank(MPI_CART_COMM, &taskid);

  /* Find up, down, left and right neighbors in the grid */
  MPI_Cart_shift(MPI_CART_COMM, 0, 1, &left, &right);
  MPI_Cart_shift(MPI_CART_COMM, 1, 1, &down, &up);

  /* Get the name of the physical node the process is running in */
  MPI_Get_processor_name(p_name, &p_name_len);
  printf("Process #%d is running in processor: %s. up=%d, down=%d, left=%d, right=%d\n", taskid, p_name, up, down, left, right);
  /*Get the coordinates of the process in the cartesian process grid */

  MPI_Cart_coords(MPI_CART_COMM, taskid, cart_ndims, coord);
  printf("Process #%d is at coordinates (%d,%d) of the cartesian grid\n", taskid, coord[0], coord[1]);

  /* The Cartesian Topology provides an efficient way to split the grid based on
  * the topology grid. We only have to find the size of each worker's grid so
  * that all the wokrers' grid combined are equal to the initial grid's size. */

  ave_row = NXPROB/cart_dims[0];
  extra_row = NXPROB%cart_dims[0];
  rows = (coord[0] == cart_dims[0] - 1) ? ave_row + extra_row : ave_row;

  ave_column = NYPROB/cart_dims[1];
  extra_column = NYPROB%cart_dims[1];
  columns = (coord[0] == cart_dims[0] - 1) ? ave_column + extra_column : ave_column;

  printf("Process #%d gets a %d x %d grid (%d x %d including the halo)\n", taskid, rows, columns, rows+2, columns+2);

  /* to prwto print, prepei na doume pws tha ginetai
  *  prtdat(NXPROB, NYPROB, u, "initial.dat");*/

  /* Allocate grid memory for this process:
  * It is very important that the memory we allocate for this grid is consistent
  * (i.e. in consecutive memory addresses) otherwise our custom MPI_Datatypes
  * for rows and more specifically for columns will not work properly. To
  * achieve this we malloc all the memory needed at start and then create
  * pointers for each row that point to this memory. */

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
  //TODO: Do these need to be freed(?) somehow?

  /* Begin doing STEPS iterations.  Must communicate border rows with
  *  neighbors.  If I have the first  or last grid row, then I only need
  *  to  communicate with one neighbor. I get ready to receive first, then
  *  send my own information. While communication is ongoing, calculate
  *  the inside grid values, which require only already available information.
  *  After communication has been completed, update the outer values of the grid.*/

  /* We store the halo rows in rows 0 and rows+1 (first and last),
  *  and the columns respectively. Elements [0][0],[0][columns+1],
  *  [rows+1][0], [rows+1][columns+1], which are the corners
  *  of the extended grid, are never used*/

  iz = 0;
  for (it = 1; it <= STEPS; it++) {
    if (left != NONE) {
      MPI_Irecv(&u[iz][1][0], 1, MPI_COLUMN, left, LTAG, MPI_COMM_WORLD, &(r_array[0]));
      MPI_Isend(&u[iz][1][1], 1, MPI_COLUMN, left, RTAG, MPI_COMM_WORLD, &(s_array[0]));
    }
    if (up != NONE) {
      MPI_Irecv(&u[iz][0][1], 1, MPI_ROW, up, UTAG, MPI_COMM_WORLD, &(r_array[1]));
      MPI_Isend(&u[iz][1][1], 1, MPI_ROW, up, DTAG, MPI_COMM_WORLD, &(s_array[1]));
    }
    if (right != NONE) {
      MPI_Irecv(&u[iz][1][columns+1], 1, MPI_COLUMN, right, RTAG, MPI_COMM_WORLD, &(r_array[2]));
      MPI_Isend(&u[iz][1][columns], 1, MPI_COLUMN, right, LTAG, MPI_COMM_WORLD, &(s_array[2]));
    }
    if (down != NONE) {
      MPI_Irecv(&u[iz][rows+1][1], 1, MPI_ROW, down, DTAG, MPI_COMM_WORLD, &(r_array[3]));
      MPI_Isend(&u[iz][rows][1], 1, MPI_ROW, down, UTAG, MPI_COMM_WORLD, &(s_array[3]));
    }
    /* Now call update to update the value of inner grid points */
    update(2, rows-1, 2, columns-1, columns, &u[iz][0][0],&u[1-iz][0][0]);

    /* Wait for the receives to be over */
    MPI_Waitall(4,r_array,MPI_STATUSES_IGNORE);

    /* Update the outer values, based on the halos we have by now received*/
    update(1, 1, 1, columns, columns, &u[iz][0][0],&u[1-iz][0][0]);
    update(rows, rows, 1, columns, columns, &u[iz][0][0],&u[1-iz][0][0]);
    update(1, rows, 1, 1, columns, &u[iz][0][0],&u[1-iz][0][0]);
    update(1, rows, columns, columns, columns, &u[iz][0][0],&u[1-iz][0][0]);

    /* Wait for the sends to be over */
    MPI_Waitall(4,s_array,MPI_STATUSES_IGNORE);

    iz = 1 - iz;
  }

  /* Final data printing */
  if (taskid!=MASTER)
  {
    /* Finally, send my portion of final results back to master */
    // MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
    MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
    // MPI_Send(&u[iz][offset][0], rows*columns, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
    MPI_Finalize();

    /* Free the allocated grid memory*/
    free(u[0]);
    free(u[1]);
  }
  else{
      /*Allocate memory to gather all grids together*/
      for (i=0;i<NYPROB;i++)
      {
        final_grid = malloc((NXPROB*NYPROB)*sizeof(float*));
      }
      /* Now wait for results from all worker tasks */
      for (i=1; i<=numtasks; i++) {
        source = i;
        msgtype = DONE;
        // MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        // MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        // MPI_Recv(&final_grid[offset_row][0], rows*NYPROB, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);
      }
      /* Write final output, call X graph and finalize MPI */
      printf("Writing final.dat file and generating graph...\n");
      prtdat(NXPROB, NYPROB, &final_grid, "final.dat");
      printf("Click on MORE button to view initial/final states.\n");
      printf("Click on EXIT button to quit program.\n");

      free(final_grid);
      MPI_Finalize();
      /* End of master code */
  }

  /* Free allocated memory */
  for(z=0; z<2; z++){
    free(u[z][0]);
    free(u[z]);
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

void inidat2(int nx, int ny, int startx, int starty, float *u) {
  int ix, iy;
  int i,j;
  for (i=1;i<nx;i++)
  {
    *(u+i*ny) = 0.0;
    *(u+i*ny+ny+2) = 0.0;

  }
  for (i=1;i<ny;i++)
  {
    *(u+i) = 0.0;
    *(u+(nx+2)*ny+i) = 0.0;
  }

  for (ix = startx; ix <= nx; ix++)
    for (iy = starty; iy <= ny; iy++)
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
