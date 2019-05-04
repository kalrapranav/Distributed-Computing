#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#define ROWS 9
#define COLS 10

int main(int argc, char *argv[]) {

   int size_mpi, rank_mpi, row_mpi, col_mpi;
   int i,j,p,ttlcols;
   int sizes[]= {2*ROWS,2*COLS};
   int subsizes[]= {ROWS,COLS};
   int starts[] = {0,0};
   int vals[ROWS][COLS];
   char hdr[] = "This is just a header.\n";
   MPI_Status stat_mpi;
   MPI_Datatype subarray;
   MPI_File fh;
   MPI_Offset offset, end_of_hdr;
   MPI_Info info_mpi;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&size_mpi);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank_mpi);

   ttlcols = 2*COLS;
   /* Where are we in the array of processes? */
   col_mpi = rank_mpi%2;
   row_mpi = rank_mpi/2;
   /* Populate the array */
   for (j=0; j<ROWS; j++){
      for (i=0; i<COLS; i++){
         vals[j][i] = ttlcols*(ROWS*row_mpi + j) +
                      COLS*col_mpi + i;
      }
   } 
   /* MPI derived datatype for setting a file view */    
   starts[0] = row_mpi*ROWS;
   starts[1] = col_mpi*COLS;
   MPI_Type_create_subarray(2, sizes, subsizes, starts,
                            MPI_ORDER_C, MPI_INT,
                            &subarray); 
   MPI_Type_commit(&subarray);
   /* open the file */    
   printf("opening file\n");
   MPI_File_open(MPI_COMM_WORLD, "arrdata.dat", 
                 MPI_MODE_WRONLY | MPI_MODE_CREATE,
                 MPI_INFO_NULL, &fh);
   printf("opened file\n");
   /* set the initial file view */    
   MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
   /* proc 0 writes first header */    
   if (rank_mpi == 0) {
      MPI_File_write(fh, (void*)hdr, strlen(hdr), MPI_CHAR, &stat_mpi);
      MPI_File_get_position(fh, &offset);
      MPI_File_get_byte_offset(fh, offset, &end_of_hdr); 
   }
   /* everybody has to know where proc 0 stopped writing */    
   MPI_Bcast((void*)&end_of_hdr, 1, MPI_INT, 0, MPI_COMM_WORLD);
   /* re-set file view for writing first array */    
   MPI_File_set_view(fh, end_of_hdr, MPI_INT,
                     subarray, "native",
                     MPI_INFO_NULL);
   /* and write the array */    
   MPI_File_write(fh, (void*)vals, ROWS*COLS, MPI_INT,
                  &stat_mpi);

   /* now go through the whole thing again to test */
   MPI_File_get_position(fh, &offset);
   MPI_File_get_byte_offset(fh, offset, &end_of_hdr); 
   MPI_File_set_view(fh, end_of_hdr, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
   if (rank_mpi == 0) {
      MPI_File_write(fh, (void*)hdr, strlen(hdr), MPI_CHAR, &stat_mpi);
      MPI_File_get_position(fh, &offset);
      MPI_File_get_byte_offset(fh, offset, &end_of_hdr); 
   }

   MPI_Bcast((void*)&end_of_hdr, 1, MPI_INT, 0, MPI_COMM_WORLD);

   MPI_File_set_view(fh, end_of_hdr, MPI_INT,
                     subarray, "native",
                     MPI_INFO_NULL);
   MPI_File_write(fh, (void*)vals, ROWS*COLS, MPI_INT,
                  &stat_mpi);
   MPI_File_close(&fh);

   MPI_Finalize();

   return 0;

}