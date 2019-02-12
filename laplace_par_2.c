/* laplace parrallel version

	author Ye Cai
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define MAX_SIZE 4096
#define EVEN_TURN 0 /* shall we calculate the 'red' or the 'black' elements */
#define ODD_TURN  1
#define FROM_MASTER 1
#define FROM_SLAVER 2
MPI_Status status;
MPI_Request request;

static double A[MAX_SIZE+2][MAX_SIZE+2]; /* (+2) - boundary elements */

struct globmem {
    int		N;		/* matrix size		*/
    int		maxnum;		/* max number of element*/
    char	*Init;		/* matrix init type	*/
    double	difflimit;	/* stop condition	*/
    double	w;		/* relaxation factor	*/
    int		PRINT;		/* print switch		*/
   // matrix	A;		/* matrix A		*/
} *glob;

/* forward declarations */
int work(int,char **);
void Init_Matrix();
void Print_Matrix();
void Init_Default();
int Read_Options(int, char **);
void Exchange_data(int rows,int rank,int size);


int 
main(int argc, char **argv)
{
    int i, timestart, timeend, iter;
 
    glob = (struct globmem *) malloc(sizeof(struct globmem));

    Init_Default();		/* Init default values	*/
    Read_Options(argc,argv);	/* Read arguments	*/
    Init_Matrix();		/* Init the matrix	*/
    iter = work(argc,argv);
    if (glob->PRINT == 1)
	Print_Matrix();
    printf("\nNumber of iterations = %d\n", iter);

    return 0;
}

int work(int argc,char **argv){
	double prevmax_even,prevmax_odd,odd,even,maxi,part_maxi,sum,w;
	int iteration=0;
	int finished=0;
	int turn = EVEN_TURN;

	int m,n,N,i;

	int rank,size;
	int rows;//rows of each slaver taking in charge

	prevmax_even = 0.0;
    prevmax_odd = 0.0;
   
    N = glob->N;
    w = glob->w;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double start_time,end_time;

	rows=N/size;

	if(rank==0){/*The master node init the matrix a and send rows to each else slaver node*/
		start_time=MPI_Wtime();
		//Init_Matrix();
		printf("program exe on %d nodes\n",size);
  		for(i=1;i<size;i++){
			MPI_Send(&A[i*rows][0],(rows+2)*(N+2),MPI_DOUBLE,i,FROM_MASTER,MPI_COMM_WORLD);
			MPI_Send(&rows,1,MPI_INT,i,FROM_MASTER,MPI_COMM_WORLD);
		}
	}

	if(rank!=0){/*The slaver nodes receive rows from master node*/
		MPI_Recv(&A[0][0],(rows+2)*(N+2),MPI_DOUBLE,0,FROM_MASTER,MPI_COMM_WORLD,&status);
		MPI_Recv(&rows,1,MPI_INT,0,FROM_MASTER,MPI_COMM_WORLD,&status);	
	}

	while (!finished) {
			
		//if(rank==0)
		iteration++;
			
		if(turn==EVEN_TURN){
			for (m = 1; m < rows+1; m++) {
				for (n = 1; n < N+1; n++) {
				    if (((m + n) % 2) == 0)
						A[m][n] = (1 - w) * A[m][n] + w * (A[m-1][n] + A[m+1][n] + A[m][n-1] + A[m][n+1]) / 4;
				}
			}

			part_maxi = -999999.0;
		    for (m = 1; m < rows+1; m++) {
				sum = 0.0;
				for (n = 1; n < N+1; n++)
			  	 	sum += A[m][n];	
				if (sum > part_maxi)			
			   		part_maxi= sum;
		    }
		    MPI_Barrier(MPI_COMM_WORLD);

		 
			if(rank!=0){
			 	MPI_Send(&part_maxi,1,MPI_DOUBLE,0,FROM_SLAVER,MPI_COMM_WORLD);
		    }
			if(rank==0){
			 	maxi=-999999.0;
				for(i=1;i<size;i++){
			 		MPI_Recv(&part_maxi,1,MPI_DOUBLE,i,FROM_SLAVER,MPI_COMM_WORLD,&status);
			 	    if(part_maxi>maxi)
		    			maxi=part_maxi;
			 	}
			 
			    if (fabs(maxi - prevmax_even) <= glob->difflimit)
					finished = 1;
				if ((iteration%100) == 0)
					printf("Iteration: %d, maxi = %f, prevmax_even = %f\n",iteration, maxi, prevmax_even);
				prevmax_even = maxi;
				
				for(i=1;i<size;i++){
					MPI_Send(&finished,1,MPI_INT,i,4,MPI_COMM_WORLD);
					//MPI_Send(&turn,1,MPI_INT,i,4,MPI_COMM_WORLD);
				}
			}
			else{
				MPI_Recv(&finished,1,MPI_INT,0,4,MPI_COMM_WORLD,&status);
				//MPI_Recv(&turn,1,MPI_INT,0,4,MPI_COMM_WORLD,&status);		
			}

			turn=ODD_TURN;
			Exchange_data(rows,rank,size);
			
			}
			else if(turn==ODD_TURN){
			for (m = 1; m < rows+1; m++) {
				for (n = 1; n < N+1; n++){
				    if (((m + n) % 2) == 1)
						A[m][n] = (1 - w) * A[m][n]  + w * (A[m-1][n] + A[m+1][n]  + A[m][n-1] + A[m][n+1]) / 4;
				}
			}

			part_maxi = -999999.0;
		    for (m = 1; m < rows+1; m++) {
				sum = 0.0;
				for (n = 1; n < N+1; n++)
			   		sum += A[m][n];	
				if (sum > part_maxi)			
			    	part_maxi = sum;
		    }

		   MPI_Barrier(MPI_COMM_WORLD);
		    
		    if(rank!=0){
		   		MPI_Send(&part_maxi,1,MPI_DOUBLE,0,FROM_SLAVER,MPI_COMM_WORLD);
		     }

		    if(rank==0){
		    	maxi=-999999.0;
		    	for(i=1;i<size;i++){
		    		MPI_Recv(&part_maxi,1,MPI_DOUBLE,i,FROM_SLAVER,MPI_COMM_WORLD,&status);
		    		if(part_maxi>maxi)
		    			maxi=part_maxi;
		    	}
		    	
			    if (fabs(maxi - prevmax_odd) <= glob->difflimit)
					finished = 1;
					//printf("rank %d finished even =1\n",rank);
				if ((iteration%100) == 0)
					printf("Iteration: %d, maxi = %f, prevmax_odd = %f\n",
					       iteration, maxi, prevmax_odd);
				prevmax_odd = maxi;
				    

				for(i=1;i<size;i++){
					MPI_Send(&finished,1,MPI_INT,i,3,MPI_COMM_WORLD);
				//	MPI_Send(&turn,1,MPI_INT,i,3,MPI_COMM_WORLD);
				}
			}
			else{
				MPI_Recv(&finished,1,MPI_INT,0,3,MPI_COMM_WORLD,&status);
			//	MPI_Recv(&turn,1,MPI_INT,0,3,MPI_COMM_WORLD,&status);
			}
			turn = EVEN_TURN;
		    Exchange_data(rows,rank,size);
		    if (iteration > 100000) {
		    /* exit if we don't converge fast enough */
				printf("Max number of iterations reached! Exit!\n");
			    finished = 1;
			}
		}
	}
	if(rank==0){
		end_time=MPI_Wtime();
		printf("execution time:%7.2fs \n",end_time-start_time);
	}
	MPI_Finalize();
	return iteration;
}





void Init_Matrix()
{
    int i, j, N, dmmy;
 
    N = glob->N;
    printf("\nsize      = %dx%d ",N,N);
    printf("\nmaxnum    = %d \n",glob->maxnum);
    printf("difflimit = %.7lf \n",glob->difflimit);
    printf("Init	  = %s \n",glob->Init);
    printf("w	  = %f \n\n",glob->w);
    printf("Initializing matrix...");
 
    /* Initialize all grid elements, including the boundary */
    for (i = 0; i < glob->N+2; i++) {
	for (j = 0; j < glob->N+2; j++) {
	    A[i][j] = 0.0;
	}
    }
    if (strcmp(glob->Init,"count") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
		A[i][j] = (double)i/2;
	    }
	}
    }
    if (strcmp(glob->Init,"rand") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
			A[i][j] = (rand() % glob->maxnum) + 1.0;
	    }
	}
    }
    if (strcmp(glob->Init,"fast") == 0) {
	for (i = 1; i < N+1; i++){
	    dmmy++;
	    for (j = 1; j < N+1; j++) {
		dmmy++;
		if ((dmmy%2) == 0)
		    A[i][j] = 1.0;
		else
		    A[i][j] = 5.0;
	    }
	}
    }

    /* Set the border to the same values as the outermost rows/columns */
    /* fix the corners */
    A[0][0] = A[1][1];
    A[0][N+1] = A[1][N];
    A[N+1][0] = A[N][1];
    A[N+1][N+1] = A[N][N];
    /* fix the top and bottom rows */
    for (i = 1; i < N+1; i++) {
		A[0][i] = A[1][i];
		A[N+1][i] = A[N][i];
    }
    /* fix the left and right columns */
    for (i = 1; i < N+1; i++) {
		A[i][0] = A[i][1];
		A[i][N+1] = A[i][N];
    }

    printf("done \n\n");
    if (glob->PRINT == 1)
		Print_Matrix();
}

void Print_Matrix()
{
    int i, j, N;
 	
    N = glob->N;
    for (i=0; i<N+2 ;i++){
    	//printf("rank %d ",rank);
	for (j=0; j<N+2 ;j++){
	    printf(" %7.2f",A[i][j]);
	}
	printf("\n");
    }
    printf("\n\n");
}

void 
Init_Default()
{
    glob->N = 2048;
    glob->difflimit = 0.00001*glob->N;
    glob->Init = "rand";
    glob->maxnum = 15.0;
    glob->w = 0.5;
    glob->PRINT = 0;
}
 
int
Read_Options(int argc, char **argv)
{
    char    *prog;
 
    prog = *argv;
    while (++argv, --argc > 0)
	if (**argv == '-')
	    switch ( *++*argv ) {
	    case 'n':
		--argc;
		glob->N = atoi(*++argv);
		glob->difflimit = 0.00001*glob->N;
		break;
	    case 'h':
		printf("\nHELP: try sor -u \n\n");
		exit(0);
		break;
	    case 'u':
		printf("\nUsage: sor [-n problemsize]\n");
		printf("           [-d difflimit] 0.1-0.000001 \n");
		printf("           [-D] show default values \n");
		printf("           [-h] help \n");
		printf("           [-I init_type] fast/rand/count \n");
		printf("           [-m maxnum] max random no \n");
		printf("           [-P print_switch] 0/1 \n");
		printf("           [-w relaxation_factor] 1.0-0.1 \n\n");
		exit(0);
		break;
	    case 'D':
		printf("\nDefault:  n         = %d ", glob->N);
		printf("\n          difflimit = 0.0001 ");
		printf("\n          Init      = rand" );
		printf("\n          maxnum    = 5 ");
		printf("\n          w         = 0.5 \n");
		printf("\n          P         = 0 \n\n");
		exit(0);
		break;
	    case 'I':
		--argc;
		glob->Init = *++argv;
		break;
	    case 'm':
		--argc;
		glob->maxnum = atoi(*++argv);
		break;
	    case 'd':
		--argc;
		glob->difflimit = atof(*++argv);
		break;
	    case 'w':
		--argc;
		glob->w = atof(*++argv);
		break;
	    case 'P':
		--argc;
		glob->PRINT = atoi(*++argv);
		break;
	    default:
		printf("%s: ignored option: -%s\n", prog, *argv);
		printf("HELP: try %s -u \n\n", prog);
		break;
	    } 
}

void Exchange_data(int rows,int rank,int size){
	int lastproc,nextproc;
	lastproc=(rank==0)? MPI_PROC_NULL:rank-1;
	nextproc=(rank==size-1)? MPI_PROC_NULL:rank+1;
	

	MPI_Sendrecv(&A[1][1], glob->N, MPI_DOUBLE, lastproc, 1000,&A[rows+1][1], glob->N, MPI_DOUBLE, nextproc, 1000,MPI_COMM_WORLD,&status);
	MPI_Sendrecv(&A[rows][1],glob->N, MPI_DOUBLE, nextproc,1000, &A[0][1], glob->N, MPI_DOUBLE, lastproc,1000, MPI_COMM_WORLD,&status);
  
}

 

