//assignment1 part 1 matrix mutiplication
//author Ye Cai
//divided the c matix in to m*2 processor



#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define FROM_MASTER 1
#define FROM_WORKER 2

#define SIZE 1024
#define DEBUG 0
static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
static double c[SIZE][SIZE];

MPI_Status status;

static void
init_matrix(void)
{
    int i, j;

    for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++) {
            /* Simple initialization, which enables us to easy check
             * the correct answer. Given SIZE size of the matrices, then
             * the output should be
             *     SIZE ... 2*SIZE ...
             *     ...
             *     2*SIZE ... 4*SIZE ...
             *     ...
             */
            a[i][j] = 1.0;
            if (i >= SIZE/2)

                a[i][j] = 2.0;
            b[i][j] = 1.0;
            if (j >= SIZE/2)
                b[i][j] = 2.0;
        }
}
static void Print(void){
        int i,j;
        for(i=0;i<SIZE;i++){
                for (j= 0; j< SIZE; ++j)
                        printf("%7.2f", c[i][j]);

                printf("\n");
        }
}

int main(int argc,char **argv)
{
       // printf("main");
        int myrank,mytype,nproc;
   		int i,j,k;
        int dest,src,offset_col,offset_row,rows,cols;
        double start_time,end_time;
        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nproc);



        if(myrank==0){//master node,responsible for send data and collect result
                printf("SIZE=%d,number of nodes=%d\n",SIZE,nproc);
                init_matrix();
                start_time=MPI_Wtime();
                
                if(nproc==1){
                    cols=SIZE;
                    rows=SIZE;
                }else{
                    rows=(SIZE/nproc)*2;
                    cols=SIZE/2;                   
                 }

                offset_row=0;
                offset_col=SIZE/2;
                for(dest=1;dest<nproc;dest++){
                        mytype=FROM_MASTER;
                        MPI_Send(&offset_col,1,MPI_INT,dest,mytype,MPI_COMM_WORLD);
                        MPI_Send(&offset_row,1,MPI_INT,dest,mytype,MPI_COMM_WORLD);
                        MPI_Send(&rows,1,MPI_INT,dest,mytype,MPI_COMM_WORLD);
                        MPI_Send(&a[offset_row][0],rows*SIZE,MPI_DOUBLE,dest,mytype,MPI_COMM_WORLD);
                        MPI_Send(&b[offset_col][0],(SIZE/2)*SIZE,MPI_DOUBLE,dest,mytype,MPI_COMM_WORLD);
                        if(dest%2==1){
                                offset_col=0;
                                offset_row+=rows;
                        }else offset_col=SIZE/2;
                }


                for(i=0;i<rows;i++){
                    for(j=0;j<cols;j++){
                            c[i][j]=0.0;
                            for(k=0;k<SIZE;k++)
                                    c[i][j]+=a[i][k]*b[k][j];
                            }
                }
                
                mytype=FROM_WORKER;
                for(src=1;src<nproc;src++){
                        MPI_Recv(&offset_row,1,MPI_INT,src,mytype,MPI_COMM_WORLD,&status);
                        MPI_Recv(&offset_col,1,MPI_INT,src,mytype,MPI_COMM_WORLD,&status);
                        MPI_Recv(&rows,1,MPI_INT,src,mytype,MPI_COMM_WORLD,&status);
                        for(i=offset_row;i<offset_row+rows;i++)
                        	MPI_Recv(&c[i][offset_col],SIZE/2,MPI_DOUBLE,src,mytype,MPI_COMM_WORLD,&status);

                        if(DEBUG)
                        printf("received %d rows c matrix from worker %d\n",rows,src);
                }
                end_time=MPI_Wtime();
                if(DEBUG)
                        Print();
                printf("execution time on %2d nodes: %f\n",nproc,end_time-start_time);
                
        }
        else{//worker node
                mytype=FROM_MASTER;
                MPI_Recv(&offset_col,1,MPI_INT,0,mytype,MPI_COMM_WORLD,&status);
                MPI_Recv(&offset_row,1,MPI_INT,0,mytype,MPI_COMM_WORLD,&status);
                MPI_Recv(&rows,1,MPI_INT,0,mytype,MPI_COMM_WORLD,&status);
                MPI_Recv(&a[offset_row][0],rows*SIZE,MPI_DOUBLE,0,mytype,MPI_COMM_WORLD,&status);
                MPI_Recv(&b[offset_col][0],SIZE*SIZE/2,MPI_DOUBLE,0,mytype,MPI_COMM_WORLD,&status);
                if(DEBUG)
                        printf("myrank=%d,offset_row=%d,offset_col=%d,rows=%d,a[offset_row][0]=%f,b[offset_col][0]=%f.\n",myrank,offset_row,offset_col,rows,a[offset_row][0],b[offset_col][0]);



                for(i=offset_row;i<offset_row+rows;i++){
                        for(j=offset_col;j<offset_col+SIZE/2;j++){
                                c[i][j]=0.0;
                                for(k=0;k<SIZE;k++)
                                        c[i][j]+=a[i][k]*b[k][j];
                        }
                }

                if(DEBUG)
                        printf("myrank=%d,offset_row=%d,offset_col=%d,rows=%d,c[offset_row][offset_col]=%f\n",myrank,offset_row,offset_col,rows,c[offset_row][offset_col]);
                mytype=FROM_WORKER;
                
                MPI_Send(&offset_row,1,MPI_INT,0,mytype,MPI_COMM_WORLD);
                MPI_Send(&offset_col,1,MPI_INT,0,mytype,MPI_COMM_WORLD);
                MPI_Send(&rows,1,MPI_INT,0,mytype,MPI_COMM_WORLD);
                for(i=offset_row;i<offset_row+rows;i++)
                     MPI_Send(&c[i][offset_col],SIZE/2,MPI_DOUBLE,0,mytype,MPI_COMM_WORLD);
        }

        MPI_Finalize();
        return 0;
}
