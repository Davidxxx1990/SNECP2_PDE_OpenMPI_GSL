//gcc -O3 main.c `pkg-config --libs gsl`

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <string.h>

#define N 500
#define L 0.5
#define H 0.05
#define dt 0.001
#define TEND 10.0
#define STEPS 10000
#define NU 0.06

struct dgl
{
	double nu;
	double k;
	double h;
	int n;
};


//DGL for the 
int f(double t, const double y[], double dydt[], void *params)
{
	struct dgl *dgl;
	double nuk;
	int i;
	
	dgl = (struct dgl*) params;
	
	
	nuk = (dgl->nu * dgl->nu) / (dgl->k * dgl->k);
	
	for(i=1; i<(dgl->n); i++)
	{
		dydt[i]     = y[i + dgl->n];
		dydt[i + dgl->n] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
	}
	
	return GSL_SUCCESS;
}

void dgl_y_init(double *y, struct dgl *dgl)
{
	int i;

	for(i=0; i<dgl->n / 2; i++)
		y[i] = 2.0 * dgl->h / (double)dgl->n * (double)i;
	for(i=dgl->n/2; i<dgl->n+1; i++)
		y[i] = 2.0 * dgl->h * (1.0 - (double)i / (double)dgl->n);

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
	int dim;
	int i;
	int steps;
	struct dgl dgl;
	double *y;
	double *uxt5, *uxt8, *ut3L4, *utL2;
	int pos_3L4, pos_L2;
	const gsl_odeiv2_step_type *type_ptr;
	gsl_odeiv2_step *step_ptr;
	gsl_odeiv2_evolve *evolve_ptr;
	gsl_odeiv2_system my_system;
	
	double t;
	
	
	dgl.nu = NU;
	dgl.k = L/N;
	dgl.n = N;
	dgl.h = H;
	
	//vector of the left side
	dim = 2 * (dgl.n + 1);
	y = calloc(sizeof(double), dim);
	
	//exitation at 5 and 8 seconds
	uxt5 = malloc(sizeof(double) * dim/2);
	uxt8 = malloc(sizeof(double) * dim/2);
	
	//time vector of position 3/4L and 1/2L
	ut3L4 = malloc(sizeof(double) * (STEPS+1));
	utL2 = malloc(sizeof(double) * (STEPS+1));
	
	dgl_y_init(y, &dgl);
	
	type_ptr = gsl_odeiv2_step_rk4; //RK4 Solver
	step_ptr = gsl_odeiv2_step_alloc(type_ptr, dim); //stepping function
	evolve_ptr = gsl_odeiv2_evolve_alloc(dim); //evolution function
	

	
	
	//DGL-System (jacobian is not used)
	my_system.function = f;	/* the right-hand-side functions dy[i]/dt */
    	my_system.jacobian = NULL;//jacobi;	/* the Jacobian df[i]/dy[j] */
    	my_system.dimension = dim;	/* number of diffeq's */
    	my_system.params = &dgl;	/* parameters to pass to rhs and jacobian */
	
	
	t = 0.0;
	//position of 3/4 L and 1/2 L
	pos_3L4 = 3 * dgl.n / 4 + 1;
	pos_L2 = dgl.n / 2 + 1;
	//init time vectors
	ut3L4[0] = y[pos_3L4];
	utL2[0] = y[pos_L2];
	
	//sim loop
	for(steps=0; steps<STEPS; steps++)
	{
		gsl_odeiv2_evolve_apply_fixed_step (evolve_ptr, NULL, step_ptr,
		                            &my_system, &t, dt, y);
		//save exitation at t=5.0
		if(steps == (int)(5.0 / dt))
			memcpy(uxt5, y, sizeof(double) * dim/2);
		//save exitation at t=8.0
		if(steps == (int)(8.0 / dt))
			memcpy(uxt8, y, sizeof(double) * dim/2);
	
		//save exitation at 1/2L and 3/4L
		ut3L4[1 + steps] = y[pos_3L4];
		utL2[1 + steps] = y[pos_L2];
	}
	//save the results
	save1("uxt5.dat", uxt5, dim/2);
	save1("uxt8.dat", uxt8, dim/2);
	save2("ut3L4.dat", ut3L4, STEPS + 1);
	save2("utL2.dat", utL2, STEPS + 1);
	
	
	//final steps
	gsl_odeiv2_evolve_free (evolve_ptr);
    	gsl_odeiv2_step_free (step_ptr);
    	
	return EXIT_SUCCESS;
}
