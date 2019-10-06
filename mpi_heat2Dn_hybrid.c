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

#define NXPROB      80                 /* x dimension of problem grid */
#define NYPROB      64                 /* y dimension of problem grid */
#define STEPS       500                /* number of time steps */
#define MAXWORKER   160                /* maximum number of worker tasks */
#define MINWORKER   1                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define UTAG        4                  /* message tag */
#define DTAG        5                  /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */
#define CONV_PERIOD 10                 /* every how many time steps to check for convergence */

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
  MPI_Request r_array[2][4];              /* Handles for receiving information */
  MPI_Request s_array[2][4];              /* Handles for sending information */
  double t_start, t_end, t_run,           /* Count the time before and after calculations*/
         t_max, t_avg;                    /* Max and average time between all processes */
  int conv_local, conv_total;             /* If there are no changes in the local table, and in all processes */


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

  ave_row = NXPROB/cart_dims[1];
  extra_row = NXPROB%cart_dims[1];
  rows = (coord[1] == cart_dims[1] - 1) ? ave_row + extra_row : ave_row;

  ave_column = NYPROB/cart_dims[0];
  extra_column = NYPROB%cart_dims[0];
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

  /* Initialize everything - including the borders - to zero */
  for (iz=0; iz<2; iz++)
    for (ix=0; ix<rows+2; ix++)
      for (iy=0; iy<columns+2; iy++)
        u[iz][ix][iy] = 0.0;
  /* Initialize table values */
  inidat(rows+2, columns+2, u[0]);

  /* Create row and column datatypes. This way we can efficiently send a column
  * from the table without having to copy it to a buffer. More specifically:
  * A row has <column-size> blocks each of which contain 1 MPI_DOUBLE and there
  * is 1 block between the starts of each block in our table.
  * A column has <row-size> blocks each of which contain 1 MPI_DOUBLE and
  * there are <row-size> blocks between the start of each block in our table. */
  MPI_Type_vector(columns, 1, 1, MPI_FLOAT, &MPI_ROW);
  MPI_Type_vector(rows, 1, columns+2, MPI_FLOAT , &MPI_COLUMN);
  MPI_Type_commit(&MPI_ROW);
  MPI_Type_commit(&MPI_COLUMN);

  /* Synchronize all tasks by waiting until they all reach this point.
  *  Once they do, start counting time.*/
  MPI_Barrier(MPI_CART_COMM);
  t_start = MPI_Wtime();

  /* Begin doing STEPS iterations.  Must communicate border rows with
  *  neighbors.  If I have the first  or last grid row, then I only need
  *  to  communicate with one neighbor. I get ready to receive first, then
  *  send my own information. While communication is ongoing, calculate
  *  the inside grid values, which require only already available information.
  *  After communication has been completed, update the outer values of the grid.*/

  /* We store the halo rows in rows 0 and rows+1 (first and last),
  *  and the columns respectively. Elements [0][0],[0][columns+1],
  *  [rows+1][0], [rows+1][columns+1], which are the corners
  *  of the extended grid, are never used. */

  /* Create persistent communication requests for each neighbor */
  for (iz=0; iz<2; iz++){
    /* Creating persistent communication with left neighbor */
    MPI_Recv_init(&u[iz][1][0], 1, MPI_COLUMN, left, LTAG, MPI_CART_COMM, &(r_array[iz][0]));
    MPI_Send_init(&u[iz][1][1], 1, MPI_COLUMN, left, RTAG, MPI_CART_COMM, &(s_array[iz][0]));
    /* Creating persistent communication with up neighbor */
    MPI_Recv_init(&u[iz][0][1], 1, MPI_ROW, up, UTAG, MPI_CART_COMM, &(r_array[iz][1]));
    MPI_Send_init(&u[iz][1][1], 1, MPI_ROW, up, DTAG, MPI_CART_COMM, &(s_array[iz][1]));
    /* Creating persistent communication with right neighbor */
    MPI_Recv_init(&u[iz][1][columns+1], 1, MPI_COLUMN, right, RTAG, MPI_CART_COMM, &(r_array[iz][2]));
    MPI_Send_init(&u[iz][1][columns], 1, MPI_COLUMN, right, LTAG, MPI_CART_COMM, &(s_array[iz][2]));
    /* Creating persistent communication with down neighbor */
    MPI_Recv_init(&u[iz][rows+1][1], 1, MPI_ROW, down, DTAG, MPI_CART_COMM, &(r_array[iz][3]));
    MPI_Send_init(&u[iz][rows][1], 1, MPI_ROW, down, UTAG, MPI_CART_COMM, &(s_array[iz][3]));
  }

  MPI_Barrier(MPI_CART_COMM);

  #pragma omp parallel
  iz = 0;
  for (it = 1; it <= STEPS; it++) {
    conv_local = 1;
    /* Request and send data to left neighbor */

    MPI_Start(&r_array[iz][0]);
    MPI_Start(&s_array[iz][0]);
    /* Request and send data to up neighbor */

    MPI_Start(&r_array[iz][1]);
    MPI_Start(&s_array[iz][1]);
    /* Request and send data to right neighbor */

    MPI_Start(&r_array[iz][2]);
    MPI_Start(&s_array[iz][2]);
    /* Request and send data to down neighbor */

    MPI_Start(&r_array[iz][3]);
    MPI_Start(&s_array[iz][3]);

    /* Now call update to update the value of inner grid points */

    if(it % CONV_PERIOD == 0){
      // conv_local &= update_check_conv(2, rows-1, 2, columns-1, columns, u[iz], u[1-iz]);
      #pragma omp for schedule(static, 4)
      for (ix = 2; ix <= rows-1; ix++){
        #pragma omp for schedule(static, 4) reduction(&: conv)
        for (iy = 2; iy <= columns-1; iy++){
          u[1-iz][ix][iy] = u[iz][ix][iy]
                     + parms.cx * ( u[iz][ix+1][iy] + u[iz][ix-1][iy] - 2.0 * u[iz][ix][iy] )
                     + parms.cy * ( u[iz][ix][iy+1] + u[iz][ix][iy-1] - 2.0 * u[iz][ix][iy] );
          if(u[1-iz][ix][iy] != u[iz][ix][iy]){
            conv_local = 0;
          }
        }
      }
    }
    else{
      // update(2, rows-1, 2, columns-1, columns, u[iz], u[1-iz]);
    }

    /* Wait to receive data from left neighbor */

    MPI_Wait(&r_array[iz][0], MPI_STATUS_IGNORE);
    /* Wait to receive data from up neighbor */

    MPI_Wait(&r_array[iz][1], MPI_STATUS_IGNORE);
    /* Wait to receive data from right neighbor */

    MPI_Wait(&r_array[iz][2], MPI_STATUS_IGNORE);
    /* Wait to receive data from down neighbor */

    MPI_Wait(&r_array[iz][3], MPI_STATUS_IGNORE);

    /* Update the outer values, based on the halos we have by now received*/
    if(it % CONV_PERIOD == 0){
      // conv_local &= update_check_conv(1, 1, 1, columns, columns, u[iz], u[1-iz]);          //up
      // conv_local &= update_check_conv(rows, rows, 1, columns, columns, u[iz], u[1-iz]);    //down
      // conv_local &= update_check_conv(1, rows, 1, 1, columns, u[iz], u[1-iz]);             //left
      // conv_local &= update_check_conv(1, rows, columns, columns, columns, u[iz], u[1-iz]); //right
    }
    else{
      // update(1, 1, 1, columns, columns, u[iz], u[1-iz]);          //up
      // update(rows, rows, 1, columns, columns, u[iz], u[1-iz]);    //down
      // update(1, rows, 1, 1, columns, u[iz], u[1-iz]);             //left
      // update(1, rows, columns, columns, columns, u[iz], u[1-iz]); //right
    }

    /* Wait for data to be sent to left neighbor */
    MPI_Wait(&s_array[iz][0], MPI_STATUS_IGNORE);
    /* Wait for data to be sent to up neighbor */
    MPI_Wait(&s_array[iz][1], MPI_STATUS_IGNORE);
    /* Wait for data to be sent to right neighbor */
    MPI_Wait(&s_array[iz][2], MPI_STATUS_IGNORE);
    /* Wait for data to be sent to down neighbor */
    MPI_Wait(&s_array[iz][3], MPI_STATUS_IGNORE);

    /* Check for convergence after every <CONV_PERIOD> iterations */
    if(it % CONV_PERIOD == 0){
      MPI_Allreduce(&conv_local, &conv_total, 1, MPI_INT, MPI_LAND, MPI_CART_COMM);
      if(conv_total){
        /* In reality we would stop here but for this project we don't stop so
        * that all test have the same number of iterations. */
      }
    }

    iz = 1 - iz;
  }
  printf("[P%03d:%s] Finished time steps\n", taskid, p_name);

  /* Calculate the total time this process has run */
  t_end = MPI_Wtime();
  t_run = t_end - t_start;
  /* Calculate the max time a process has run, between all the processes. The
  * process with id == 0, takes care of gathering the result and printing it. */
  MPI_Reduce(&t_run, &t_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_CART_COMM);
  if(taskid == 0){
    printf("Max time = %f\n", t_max);
  }
  /* Calculate the average time a process has run, between all the processes */
  MPI_Reduce(&t_run, &t_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_CART_COMM);
  if(taskid == 0){
    t_avg = t_avg / (double) numtasks;
    printf("Average time = %f\n", t_avg);
  }

  /* Final data printing */
  if (taskid!=MASTER){
    /* Finally, send my portion of final results back to master */
    // MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
    //MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
    // MPI_Send(&u[iz][offset][0], rows*columns, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
  }
  else{
    /* Write final output, call X graph and finalize MPI */
    printf("Writing final.dat file and generating graph...\n");
    //prtdat(NXPROB, NYPROB, &final_grid, "final.dat");
    printf("Click on MORE button to view initial/final states.\n");
    printf("Click on EXIT button to quit program.\n");
  }

  /* Free the requests allocated*/
  for (iz=0;iz<2;iz++){
    MPI_Request_free(&(r_array[iz][0]));
    MPI_Request_free(&(s_array[iz][0]));
    MPI_Request_free(&(r_array[iz][1]));
    MPI_Request_free(&(s_array[iz][1]));
    MPI_Request_free(&(r_array[iz][2]));
    MPI_Request_free(&(s_array[iz][2]));
    MPI_Request_free(&(r_array[iz][3]));
    MPI_Request_free(&(s_array[iz][3]));
  }

  /* Free custom types we have created */
  MPI_Type_free(&MPI_ROW);
  MPI_Type_free(&MPI_COLUMN);

  /* Free allocated memory */
  for(z=0; z<2; z++){
    free(u[z][0]);
    free(u[z]);
  }
  MPI_Finalize();
  return 0;
}

/****************************** subroutine update *****************************/
void update(int x_start, int x_end, int y_start, int y_end,int ny, float **u1, float **u2) {
  int ix, iy;
  #pragma omp for schedule(static, 4)
  for (ix = x_start; ix <= x_end; ix++){
    #pragma omp for schedule(static, 4)
    for (iy = y_start; iy <= y_end; iy++){
      u2[ix][iy] = u1[ix][iy]
                 + parms.cx * ( u1[ix+1][iy] + u1[ix-1][iy] - 2.0 * u1[ix][iy] )
                 + parms.cy * ( u1[ix][iy+1] + u1[ix][iy-1] - 2.0 * u1[ix][iy] );
    }
  }
}

int update_check_conv(int x_start, int x_end, int y_start, int y_end,int ny, float **u1, float **u2) {
  int ix, iy;
  int conv = 1;

  #pragma omp for schedule(static, 4) reduction(&: conv) 
  for (ix = x_start; ix <= x_end; ix++){
    #pragma omp for schedule(static, 4)
    for (iy = y_start; iy <= y_end; iy++){
      u2[ix][iy] = u1[ix][iy]
                 + parms.cx * ( u1[ix+1][iy] + u1[ix-1][iy] - 2.0 * u1[ix][iy] )
                 + parms.cy * ( u1[ix][iy+1] + u1[ix][iy-1] - 2.0 * u1[ix][iy] );
      if(u2[ix][iy] != u1[ix][iy]){
        conv = 0;
      }
    }
  }
  return conv;
}

/****************************** subroutine inidat *****************************/
void inidat(int nx, int ny, float **u) {
  int ix, iy;

  for (ix = 0; ix <= nx-1; ix++) {
    for (iy = 0; iy <= ny-1; iy++) {
      u[ix][iy] = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
    }
  }
}

/****************************** subroutine prtdat *****************************/
void prtdat(int nx, int ny, float **u1, char *fnam) {
  int ix, iy;
  FILE *fp;

  fp = fopen(fnam, "w");
  for (iy = ny-1; iy >= 0; iy--) {
    for (ix = 0; ix <= nx-1; ix++) {
      fprintf(fp, "%6.1f", u1[ix][iy]);
      if (ix != nx-1)
        fprintf(fp, " ");
      else
        fprintf(fp, "\n");
    }
  }
  fclose(fp);
}

/****************************** subroutine prtfdat *****************************/
void prtfdat(int nx, int ny, float **u1) {
  int ix, iy;

  for (iy = 0; iy <= ny-1; iy++) {
    for (ix = 0; ix <= nx-1; ix++) {
      printf("%6.1f", u1[ix][iy]);
      if (ix != nx-1)
        printf(" ");
      else
        printf("\n");
    }
  }
}
