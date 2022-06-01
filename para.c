//gcc -O3 main.c `pkg-config --libs gsl`

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <string.h>
#include <mpi.h>

#define N 500
#define L 0.5
#define H 0.05
#define dt 0.001
#define TEND 10.0
#define STEPS 10000
#define NU 0.06


//#define DEBUG

#define RIGHT 0
#define LEFT 1

struct dgl
{
	double nu;
	double k;
	double h;
	int n;
	int lsize;
	double y_left, y_right;
};

int max(int x, int y)
{
	if(x > y)
		return x;
	else
		return y;
}

int min(int x, int y)
{
	if(x < y)
		return x;
	else
		return y;
}

//DGL for the 

int f(double t, const double y[], double dydt[], void *params)
{
	struct dgl *dgl;
	double nuk;
	
	int i;
	int world_rank, world_size;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	dgl = (struct dgl*) params;
	
	
	nuk = (dgl->nu * dgl->nu) / (dgl->k * dgl->k);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_size > 1)
	{
		//communicate the edges
		if(world_rank == 0)
		{
			MPI_Send(&y[dgl->lsize - 1], 1, MPI_DOUBLE, world_rank + 1, RIGHT, MPI_COMM_WORLD);
			MPI_Recv(&dgl->y_right, 1, MPI_DOUBLE, world_rank + 1, LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else if(world_rank == (world_size-1))
		{
			MPI_Send(&y[0], 1, MPI_DOUBLE, world_rank - 1, LEFT, MPI_COMM_WORLD);
			MPI_Recv(&dgl->y_left, 1, MPI_DOUBLE, world_rank - 1, RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Send(&y[dgl->lsize - 1], 1, MPI_DOUBLE, world_rank + 1, RIGHT, MPI_COMM_WORLD);
			MPI_Send(&y[0], 1, MPI_DOUBLE, world_rank - 1, LEFT, MPI_COMM_WORLD);
			MPI_Recv(&dgl->y_right, 1, MPI_DOUBLE, world_rank + 1, LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&dgl->y_left, 1, MPI_DOUBLE, world_rank - 1, RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		
		//calculate the right side of the DGL
		if(world_rank == 0)
		{		
			for(i=1; i<(dgl->lsize - 1); i++)
			{
				dydt[i]     = y[i + dgl->lsize];
				dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
			}
			dydt[i]     = y[i + dgl->lsize];
			dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * dgl->y_right;
		}
		else if(world_rank == (world_size - 1))
		{
			i=0;
			dydt[i]     = y[i + dgl->lsize];
			dydt[i + dgl->lsize] = nuk * dgl->y_left - 2 * nuk * y[i]  + nuk * y[i+1];	
			for(i=1; i<(dgl->lsize - 1); i++)
			{
				dydt[i]     = y[i + dgl->lsize];
				dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
			}		
		}
		else
		{
			i=0;
			dydt[i]     = y[i + dgl->lsize];
			dydt[i + dgl->lsize] = nuk * dgl->y_left - 2 * nuk * y[i]  + nuk * y[i+1];	
			for(i=1; i<(dgl->lsize - 1); i++)
			{
				dydt[i]     = y[i + dgl->lsize];
				dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
			}
			dydt[i]     = y[i + dgl->lsize];
			dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * dgl->y_right;
		}
	}
	else
	{
		for(i=1; i<(dgl->n); i++)
		{
			dydt[i]     = y[i + dgl->n];
			dydt[i + dgl->n] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
		}
	}
	
	return GSL_SUCCESS;
}



void dgl_y_init(double *y, struct dgl *dgl)
{
	int i, imin, imax;
	int world_rank, world_size;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	
	imin = (dgl->n + 1)  / world_size * world_rank;
	for(i=0; i<world_rank; i++)
	{
		if( ((dgl->n + 1) % world_size) > 0 && (world_size - ((dgl->n + 1) % world_size)) <= i )
			imin++;
	}
	
	
	imax = imin + dgl->lsize;

#ifdef DEBUG
	printf("%d imin: %d imax: %d\n", world_rank, imin, imax);
#endif

	for(i=max(0, imin); i<min(dgl->n / 2, imax); i++)
		y[i-imin] = 2.0 * dgl->h / (double)dgl->n * (double)i;
	for(i=max(dgl->n/2, imin); i<min(dgl->n+1, imax); i++)
		y[i-imin] = 2.0 * dgl->h * (1.0 - (double)i / (double)dgl->n);

}


//save function; save pos. and y in file
void save1(char *file, double *y, int n)
{
	FILE *fd;
	int i;
	
	fd=fopen(file,"w");
	
	for(i=0; i<n; i++)
	{
		fprintf(fd,"%E \t", (double)L/(double)N * i);
		fprintf(fd,"%E ",y[i]);
		
		fprintf(fd,"\n");
	}
	fclose(fd);

}

//save function; save t and y in file
void save2(char *file, double *y, int n)
{
	FILE *fd;
	int i;
	
	fd=fopen(file,"w");
	
	for(i=0; i<n; i++)
	{
		fprintf(fd,"%E \t", dt * i);
		fprintf(fd,"%E ",y[i]);
		
		fprintf(fd,"\n");
	}
	fclose(fd);

}

int main()
{
	int world_rank, world_size;
	int dim;
	int i, n, len, imin, imax;
	int steps;
	struct dgl dgl;
	double *gy, *y, *y0;
	double *uxt5, *uxt8, *ut3L4, *utL2;
	double *guxt5, *guxt8, *gut3L4, *gutL2;
	int pos_3L4, pos_L2;
	const gsl_odeiv2_step_type *type_ptr;
	gsl_odeiv2_step *step_ptr;
	gsl_odeiv2_evolve *evolve_ptr;
	gsl_odeiv2_system my_system;
	
	double t;
	
	//Initialisierung
	MPI_Init(NULL, NULL);
    	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	
	dgl.nu = NU;
	dgl.k = L / N;
	dgl.n = N;
	dgl.h = H;

	dgl.lsize = (dgl.n + 1)  / world_size;

	if( ((dgl.n + 1) % world_size) > 0 && (world_size - ((dgl.n + 1) % world_size)) <= world_rank )
		dgl.lsize++;
#ifdef DEBUG		
	printf("%3d: n=%d\n", world_rank, dgl.lsize);
#endif
	//vector of the left side
	dim = dgl.lsize * 2;
	y = calloc(sizeof(double), dim);
	
	//exitation at 5 and 8 seconds
	uxt5 = malloc(sizeof(double) * dim/2);
	uxt8 = malloc(sizeof(double) * dim/2);
	
	
	//global index imin and imax
	imin = (dgl.n + 1)  / world_size * world_rank;
	for(i=0; i<world_rank; i++)
	{
		if( ((dgl.n + 1) % world_size) > 0 && (world_size - ((dgl.n + 1) % world_size)) <= i )
			imin++;
	}
	imax = imin + dgl.lsize;
	
	//position of 3/4 L and 1/2 L
	pos_3L4 = 3 * dgl.n / 4 + 1;
	pos_L2 = dgl.n / 2 + 1;
	
	//time vector of position 3/4L and 1/2L
	ut3L4 = malloc(sizeof(double) * (STEPS+1));
	utL2 = malloc(sizeof(double) * (STEPS+1));
	
	//init the vector of the left side
	dgl_y_init(y, &dgl);


		
	
	type_ptr = gsl_odeiv2_step_rk4; //RK4 Solver
	step_ptr = gsl_odeiv2_step_alloc(type_ptr, dim); //stepping function
	evolve_ptr = gsl_odeiv2_evolve_alloc(dim); //evolution function
	

	
	
	//DGL-System (jacobian is not used)
	my_system.function = f;	/* the right-hand-side functions dy[i]/dt */
    	my_system.jacobian = NULL;//jacobi;	/* the Jacobian df[i]/dy[j] */
    	my_system.dimension = dim;	/* number of diffeq's */
    	my_system.params = &dgl;	/* parameters to pass to rhs and jacobian */
	
	
	//init time vectors
	t = 0.0;
	if(imin <= pos_3L4 && pos_3L4 <= imax)
		ut3L4[0] = y[pos_3L4-imin];
	if(imin <= pos_L2 && pos_L2 <= imax)
		utL2[0] = y[pos_L2-imin];
	
	//sim loop
	MPI_Barrier(MPI_COMM_WORLD);
	for(steps=0; steps<STEPS; steps++)
	{
	
		gsl_odeiv2_evolve_apply_fixed_step (evolve_ptr, NULL, step_ptr,
		                            &my_system, &t, dt, y);
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		//save exitation at t=5.0
		if(steps == (int)(5.0 / dt))
		{
			if(world_rank == 0)
			{
				guxt5 = malloc(sizeof(double) * (N+1));
				memcpy(guxt5, y, sizeof(double) * dgl.lsize);
				len = dgl.lsize;
				for(i=1; i<world_size; i++)
				{
					MPI_Recv(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&guxt5[len], n, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					len+=n;
				}
				
			}
			else
			{
				MPI_Send(&dgl.lsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(y, dgl.lsize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			}

		}
		//save exitation at t=8.0
		if(steps == (int)(8.0 / dt))
		{
			if(world_rank == 0)
			{
				guxt8 = malloc(sizeof(double) * (N+1));
				memcpy(guxt8, y, sizeof(double) * dgl.lsize);
				len = dgl.lsize;
				for(i=1; i<world_size; i++)
				{
					MPI_Recv(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&guxt8[len], n, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					len+=n;
				}
				
			}
			else
			{
				MPI_Send(&dgl.lsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(y, dgl.lsize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			}
			
		}
		//save exitation at 1/2L and 3/4L
		if(imin <= pos_3L4 && pos_3L4 <= imax)
			ut3L4[1 + steps] = y[pos_3L4-imin];
		
		if(imin <= pos_L2 && pos_L2 <= imax)
			utL2[1 + steps] = y[pos_L2-imin];
	}
	//get the results and save it
	if(world_rank == 0)
	{
		save1("uxt5.dat", guxt5, N+1);
		save1("uxt8.dat", guxt8, N+1);
		
		if(!(imin <= pos_3L4 && pos_3L4 <= imax))
			MPI_Recv(ut3L4, (STEPS+1), MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		save2("ut3L4.dat", ut3L4, STEPS + 1);
		
		if(!(imin <= pos_L2 && pos_L2 <= imax))
			MPI_Recv(utL2, (STEPS+1), MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		save2("utL2.dat", utL2, STEPS + 1);
		
	}
	else
	{
		if(imin <= pos_3L4 && pos_3L4 <= imax)
			MPI_Send(ut3L4, (STEPS+1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		if(imin <= pos_L2 && pos_L2 <= imax)
			MPI_Send(utL2, (STEPS+1), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}
	//final steps
	gsl_odeiv2_evolve_free (evolve_ptr);
    	gsl_odeiv2_step_free (step_ptr);
    	
    	MPI_Finalize();
    	
	return EXIT_SUCCESS;
}


