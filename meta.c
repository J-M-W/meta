#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mt19937ar.c"
#include "mt19937ar.h"
#include "potentials.h"

#define DIVIDER "--------------------\n"

typedef struct
{
	double* r;  /* Position */
	double E;   /* Energy at r */
	double h;   /* Height of Gaussian */
	double sig; /* Standard deviation of Gaussian */
}gauss;

typedef struct 
{
	/* PARAMETERS READ FROM INPUT */
	int nDim; /* Dimensionality of CV-space */
	int calcMethod[2]; /* 0th element: zeroth derivatives. 1st element: 2nd derivatives. 0: on-the-fly analytic; 1: on-the-fly numerical; 2: tabulated analytic; 3: tabulated numerical */
	int lscape; /* Landscape number */
	/* Dynamical parameters */
	double T;
	double dt;
	double g;
	double m;
	/* Gaussian parameters */
	double sigma;
	double w;
	double Vmax;
	/* Interpolation parameters */
	double n; /* Resolution of grid */
	double* bounds;
	double* deltas;
	double* eGrid; /* Energy landscape used for interpolation */
	double** fGrid; /* grad(E) used for interpolation. 2nd index corresponds to different dimensions */

	/* LANDSCAPE-SPECIFIC PARAMETERS */
	/* Harmonic parameters */
	double k;

	/* 'DERIVED' PARAMETERS */
	double c1;
	double c3;
}argstruct;



unsigned int get_input(char**);
int get_force(double* r, double* f, gauss* gaussians, int nGauss, argstruct* args);
void evolve_cvs(double* r, double* p, double* f, gauss* gaussians, int nGauss, argstruct* args);
double** malloc2d(int nx, int ny);
int initialise_metadynamics(int, char**, struct argstruct*);
double* grad(double* r, double h, struct argstruct*);
double* get_bounds(double**, int, int);
double** calc_landscape(double* bounds, int n, argstruct* args, gauss* gaussians, int nGauss, int tStep, FILE* f);
void recurse_landscape(int** indices, int* indCount, int* indTemp, int digit, int nDim, int n);
int append(double*** arrayPtr, int* rowPtr, int col, double* newRow, int rowIndex);
double rand_normal();
double eval_gaussian(gauss* gaussian, double* r, int nDim);
double eval_potential(double* r, argstruct* args, int deriv, int derivCoord);
double* populate_landscape(double* bounds, int n, argstruct* args);

int main(int argc, char* argv[])
{
    double dt=1E-2;
    int tStepsInit  = 1E1;        /* Initial number of timesteps - arrays will be extended beyond this if necessary */
	int tStepsLimit = 1E4;        /* Hard limit on the number of timesteps - simulation terminates at this point */ 
	int tStep = 0;                /* Number of tSteps simulated */
    double r[] = {0.2,0.5};    /* Position vector in configuration space, and its initial value */
	//double r[] = {4,1};    /* Position vector in configuration space, and its initial value */
    double p[] = {0,0};      /* Momentum vector in configuration space, and its initial value */
    int gaussInterval=10;
	int n=50; /* Dimensions on grid to evaluate the energy landscape on */

    double **visited; /* Array of visited positions */
	int visitedSize;  /* Number of rows in visited array */
    double **muGauss; /* Array of Gaussian means (a subset of 'visited') */
	int gaussArraySize;  /* Size of muGauss (at first) */
	argstruct args;   /* Struct of parsed args */
	int nGauss=0;     /* Number of Gaussians placed so far */
    int i, j, nDim;
    double* f;
	double* bounds;
	double* trajEnergy;     /* True energy at each visited point on the trajectory */
	double* metaTrajEnergy; /* Same as above, but with metadynamics Gaussians included */
    FILE* fOut;
	gauss* gaussians;

	/* Parse command line args; store them in the struct */
	//memset(&args, 0, sizeof(argstruct));
	initialise_metadynamics(argc, argv, &args);
	
	nDim = args.nDim;
	f = (double*)malloc(nDim*sizeof(double));

    /* Allocate visited position arrays */
	visitedSize = tStepsInit;
    visited = malloc2d(visitedSize, nDim);
	gaussArraySize = (visitedSize/gaussInterval)+1;
	printf("Initial array sizes: %d visited points, %d Gaussians\n", visitedSize, gaussArraySize);
    //muGauss = malloc2d(muGaussSize, nDim);
	gaussians = (gauss*)malloc(gaussArraySize*sizeof(gauss));
	/* Initialise force */
	get_force(r, f, gaussians, nGauss, &args);

    /* Begin simulation */
    for(tStep=0; tStep<tStepsLimit; tStep++)
    {
        /* Get force and iterate */
		//puts("Evolving CVs");
        evolve_cvs(r, p, f, gaussians, nGauss, &args);
        /* Append to visited */
		//puts("Appending to visited location data");
		append(&visited, &visitedSize, args.nDim, r, tStep);
		//printf("VISITED: (%lf, %lf)\n", visited[tStep][0], visited[tStep][1]);
		//puts("Done appending.");
		//puts("Adding Gaussian");
        if((tStep%gaussInterval)==0)
        {
			//append(&muGauss, &muGaussSize, nDim, r, nGauss);
			if(nGauss >= gaussArraySize)
			{
				gaussians = (gauss*)realloc(gaussians, 2*gaussArraySize*sizeof(gauss));
				gaussArraySize *= 2;
			}
			gaussians[nGauss].r = (double*)malloc(nDim*sizeof(double));
			// TODO - free these at the end. Memory leak.
			memcpy((gaussians[nGauss]).r,  r, nDim*sizeof(double));
			gaussians[nGauss].E = eval_potential(r, &args, 0, 0);
			gaussians[nGauss].h = (args.w)*exp(-gaussians[nGauss].E/args.Vmax);
			printf("E=%lf, h=%lf\n", gaussians[nGauss].E, gaussians[nGauss].h);
			gaussians[nGauss].sig = args.sigma;
            nGauss++;
        }
        //printf("Iteration %d\n", tStep);
		if(r[0] < -1000)
		{
			puts("DEBUG HALT");
		}
        //printf("r = (%lf, %lf)\n", r[0], r[1]);
        //printf("F = (%lf, %lf)\n", f[0], f[1]);
        //printf("--------------------------\n");
        //if(i>1E6){exit(EXIT_FAILURE);}
    }

    /* WRITE RESULTS */
    puts("Writing results to disk");
    fOut = fopen("output.txt", "w");
    if(fOut==NULL)
    {
        puts("I/O Error.");
        exit(EXIT_FAILURE);
    }
    /* Write parameters to file */
	fprintf(fOut, DIVIDER);
    fprintf(fOut, "PARAMETERS:\n");
	fprintf(fOut, "VERSION = 1.0\n");
	fprintf(fOut, "DIMENSIONS = %d\n", nDim);
    fprintf(fOut, "TIMESTEP = %g\n", dt);
    fprintf(fOut, "ITERATIONS = %d\n", tStep);
    fprintf(fOut, "INTERVAL BETWEEN GAUSSIANS = %d\n", gaussInterval);
	fprintf(fOut, "GAUSSIAN SIGMA = %g\n", args.sigma);
	fprintf(fOut, "GAUSSIAN W = %g\n", args.w);

	/* Landscape-specific parameters */
	if(args.lscape==0)
	{
		fprintf(fOut, "LANDSCAPE = MUELLER\n");
	}

    fprintf(fOut, DIVIDER);

    /* Write visited positions to file */
    fprintf(fOut, "TRAJECTORY (first column is timestep):\n");
    for(i=0; i<tStep; i++)
    {
		/* Print timestep number */
        fprintf(fOut, "%d", i);
		/* Iterate over dimensions */
		for(j=0; j<nDim; j++)
		{
			fprintf(fOut, "     %g", visited[i][j]);
		}
		fprintf(fOut, "\n");
    }
    fprintf(fOut, DIVIDER);

	/* Write energies to file */
	fprintf(fOut, "ENERGIES:\n");
	for(i=0; i<tStep; i++)
	{
		/* TODO: have the energies stored during the simulation time so this doesn't need recalculating */
		fprintf(fOut, "%g\n", eval_potential(visited[i], &args, 0, 0));
	}
	fprintf(fOut, DIVIDER);

	/* Write Gaussian positions to file */
	fprintf(fOut, "GAUSSIANS' LOCATIONS:\n");
	for(i=0; i<nGauss; i++)
	{
		/* Print timestep number */
		if(fprintf(fOut, "%d", i*gaussInterval) < 0)
		{
			puts("fprintf failed.");
			exit(EXIT_FAILURE);
		}
		/* Iterate over dimensions */
		for(j=0; j<nDim; j++)
		{
			fprintf(fOut, "     %g", gaussians[i].r[j]);
		}
		fprintf(fOut, "\n");
	}
    fprintf(fOut, DIVIDER);

	/* Evaluate energy landscape on a grid */
	puts("Getting bounds");
	bounds = get_bounds(visited, nDim, tStep);
	//bounds[0] = -1.5;
	//bounds[1] = 0.7;
	//bounds[2] = -0.2;
	//bounds[3] = 1.5;
	puts("Calculating energy landscape");
	calc_landscape(bounds, n, &args, gaussians, nGauss, tStep, fOut);

	//for(i=0; i<2*nDim; i++){printf("%g\n", bounds[i]);}
	/* Flush output and check for write errors. */
    puts("Done.");
	fflush(fOut);
	if(ferror(fOut) != 0)
	{
		puts("Output file error.");
		fclose(fOut);
		exit(EXIT_FAILURE);
	}
    fclose(fOut);
    return 0;
}

void evolve_cvs(double* r, double* p, double* f, gauss* gaussians, int nGauss, argstruct* args)
{
	/* Position r, momentum p, force f, mass m, timestep dt, temperature T, friction constant g, in nDim-dimensional space */
	double c1, c3;
	int i;
	double dt = args->dt;
	double m  = args->m;
    /* Euler integration :( */
	//int i;
	//for(i=0; i<nDim; i++)
	//{
	//	/* TODO: make friction a proper parameter */
	//	//f[i] += -0.0005*v[i];
	//	/* TODO: make +=? */
	//	r[i] = f[i]*dt;
	//}

	/* BAOAB integration scheme */
	c1 = exp(-(args->g)*dt);
	c3 = sqrt((args->T)*(1-c1*c1));
	for(i=0; i<(args->nDim); i++)
	{
		p[i] += dt*f[i]*0.5;
		r[i] += dt*p[i]/(2*m);
		p[i]  = c1*p[i] + c3*pow(m, 0.5)*rand_normal();
		r[i] += dt*p[i]/(2*m); 
	}
    get_force(r, f, gaussians, nGauss, args);
	for(i=0; i<(args->nDim); i++){p[i] += f[i]*dt*0.5;}
}

double eval_harmonic(double* r, argstruct* args)
{
	/* r must be of length ncv */
	int i;
	double energy=0;
	for(i=0; i<(args->nDim); i++)
	{
		energy += (args->k)*r[i]*r[i];
	}
	return energy;
}

double eval_double_harmonic(double* r, argstruct* args)
{
	/* r must be of length ncv */
	int i;
	double energy=0;
	for(i=0; i<(args->nDim); i++)
	{
		energy += (args->k)*r[i]*r[i]*r[i]*r[i];
		energy -= (args->k)*r[i]*r[i]*2;
	}
	return energy;
}

double eval_potential(double* r, argstruct* args, int deriv, int derivCoord)
{
	if(args->calcMethod[deriv]==0)
	{
		switch(args->lscape)
		{
			case 0:
				return mueller(r, deriv, derivCoord);
			default:
				printf("No potential with code %d\n", args->lscape);
				exit(EXIT_FAILURE);
		}
	}
	else
	{
		puts("Invalid arguments to eval_potential");
		exit(EXIT_FAILURE);
	}
}

double* grad(double* r, double h, argstruct* args)
{
	/* Numerically evaluates the gradient operator */
	double* rPlus;
	double* rMinus;
	double* output;
	int nDim = args->nDim;
	int i;

	/* Initialise the array which will contain gradient vector */
	output = (double*)malloc(nDim*sizeof(double));

	/* Arrays for the upper and lower finite difference points */
	rPlus  = (double*)malloc(nDim*sizeof(double));
	rMinus = (double*)malloc(nDim*sizeof(double));

	/* Iterate over components of grad */
	for(i=0; i<nDim; i++)
	{
		memcpy(rPlus,  r, nDim*sizeof(double));
		memcpy(rMinus, r, nDim*sizeof(double));
		rPlus[i]  += h;
		rMinus[i] -= h;
		/* Use f(x+h)-f(x-h)/2h to evaluate the partial derivative */
		output[i] = (eval_potential(rPlus, args, 0, 0) - eval_potential(rMinus, args, 0, 0))/(2*h);
		//printf("E_+ = %lf, E_- = %lf\n", args->energyPtr(rPlus, args), args->energyPtr(rMinus, args));
		//printf("r_%d = %lf, grad_%d = %lf\n", i, r[i], i, output[i]);
	}
	free(rPlus);
	free(rMinus);
	return output;
}

int get_force(double* r, double* f, gauss* gaussians, int nGauss, argstruct* args)
{
    /* f must already have been allocated! */
    int i, j;
    double w=(*args).w;
    double ds=(*args).sigma;
    double expFactor;
	double h=1E-8;
	double height;
	double* gradient;

    /* Numerically differentiate the energy to get the force */
    //f[0] = -r[0];
    //f[1] = -r[1];
	gradient = grad(r, h, args);
	for(i=0; i<(args->nDim); i++)
	{
		f[i] = -gradient[i];
		//printf("-grad U_%d = %lf\n", i, f[i]);
	}
	free(gradient);
	//printf("Force before Gaussians = %lf\n", f[0]);
    /* Metadynamics forces */
    for(i=0; i<nGauss; i++)
    {
		//puts("GAUSSIAN TIME OH YEAH");
		/* Construct the exponential factor, exp(-r^2/2ds^2) */
		expFactor=0;
		for(j=0; j<(args->nDim); j++)
		{
			/* Contribution to the argument of the exponent from the j-th dimensional coordinates */
			//expFactor += gaussians[i].r[j];
			expFactor += (r[j]-gaussians[i].r[j])*(r[j]-gaussians[i].r[j]);
			//printf("!!! %lf, %lf\n", r[j], (gaussians[i].r)[j]);
			//printf("%lf\n", expFactor);
		}
		//printf("Argument of exponential = %lf\n", expFactor);
		expFactor = exp(-expFactor/(2*ds));
		//printf("exponential = %lf\n", expFactor);
		/* Get the component along each axis */
		for(j=0; j<(args->nDim); j++)
		{
			f[j] += (gaussians[i].h)*(r[j]-gaussians[i].r[j])*expFactor/ds;
		}
		//printf("h=%lf\n", gaussians[i].h);
    }
	//printf("Force after Gaussians =  %lf\n", f[0]);
	//puts("--------------------------");
	return 0;
}

unsigned int get_input(char** strptr)
{
    /* Reads arbitrary-length standard input. Will return string length, with 0 for failure. */
    char tmp;
    unsigned int i=0; /* i counts the number of characters */

    *strptr = NULL;
    /* Read in the characters one-by-one */
    do
    {
        i++;
        tmp = getchar();
        /* Add space for one more character to the string */
        *strptr = (char*)realloc(*strptr, i*sizeof(char));
        if(*strptr == NULL)
		{
			puts("Could not allocate memory.");
			exit(EXIT_FAILURE);
		}
        /* Write tmp to strptr iff tmp is not EOL */
        (*strptr)[i-1] = tmp=='\n' ? '\0' : tmp;
    } while((tmp != '\n') && (tmp != '\0'));
    return i;
}

double** malloc2d(int nx, int ny)
{
    /* Allocates a 2d array of doubles with dimensions (nx, ny), exiting the program if malloc fails */
    double** array;
    int i;

    array = (double**)malloc(nx*sizeof(double));
    if(array==NULL)
    {
        puts("Malloc failed in function malloc2d.");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<nx; i++)
    {
        array[i] = (double*)malloc(ny*sizeof(double));
        if(array[i]==NULL)
        {
            puts("Malloc failed in function malloc2d.");
            exit(EXIT_FAILURE);
        }
    }
    return array;
}

int initialise_metadynamics(int argc, char* argv[], argstruct* structptr)
{
	double default_Vmax[] = {25};
	structptr->Vmax=-1.0;

	/* PARSE */

	/* TO BE IMPLEMENTED LATER, with input file; for now, just set the same options each time */
	structptr->nDim=2;
	structptr->lscape=0;
	structptr->k=1;
	structptr->sigma=1;
	structptr->w=1.0000;
	structptr->T=0.1;
	structptr->dt=1E-2;
	structptr->g=1E1;
	structptr->m=1;
	structptr->bounds = malloc(2*(structptr->nDim)*sizeof(double));
	structptr->bounds[0] = -5;
	structptr->bounds[1] = 5;
	structptr->bounds[2] = -5;
	structptr->bounds[3] = 5;
	structptr->calcMethod[0] = 0;
	structptr->calcMethod[1] = 1;
	/* Use default arguments if not manually overridden */
	if(structptr->Vmax<0){structptr->Vmax = default_Vmax[structptr->lscape];}

	/* INITIALISE */
	switch(structptr->calcMethod[0]) /* Calculation method for E(r) */
	{
		case 0:
			break;
		case 2:
			structptr->eGrid = populate_landscape(structptr->bounds, structptr->n, structptr);
			break;
		case 3:
			structptr->eGrid = populate_landscape(structptr->bounds, structptr->n, structptr);
			break;
		default:
			printf("No calculation method for E(r) with code %d\n", structptr->calcMethod[0]);
			exit(EXIT_FAILURE);
	}
	switch(structptr->calcMethod[1]) /* Calculation method for grad E(r) */
	{
		case 0:
			break;
		case 1:
			break;
		default:
			printf("No calculation method for grad E(r) with code %d\n", structptr->calcMethod[1]);
			exit(EXIT_FAILURE);
	}

	return 1;
}

int append(double*** arrayPtr, int* rowPtr, int col, double* newRow, int rowIndex)
{
	/* TODO: if memory use becomes an issue, ensure that the array is never expanded beyond tStepsMax */
	/* Writes newRow to the array, doubling the number of rows if necessary.
	rowIndex should be the index of the new row to be filled, and rowPtr a pointer to the number of rows total in the array.
	col should be the number of columns.
	array should be a pointer to a 2D array. rowIndex must be less than half the size of *rowPtr; THIS IS NOT CHECKED FOR. */
	int i;
	//printf("%d, %d\n", *rowPtr, rowIndex);
	/* Double the array size if necessary */
	if(rowIndex >= (*rowPtr))
	{
		/* Realloc the array */
		puts("REALLOCING");
		*arrayPtr = (double**)realloc(*arrayPtr, 2*(*rowPtr)*sizeof(double*));
		if((*arrayPtr)==NULL)
		{
			puts("Memory error - realloc failed in function 'append'.");
			exit(EXIT_FAILURE);
		}
		/* Malloc all the newly added rows */
		for(i=(*rowPtr); i<2*(*rowPtr); i++)
		{
			(*arrayPtr)[i] = (double*)malloc(col*sizeof(double));
			if((*arrayPtr)==NULL)
			{
				puts("Memory error - malloc failed in function 'append'.");
				exit(EXIT_FAILURE);
			}
		}

		/* Double the row count */
		(*rowPtr) *= 2;
	}

	/* Write the new row */
	//printf("nDim = %d\n", col);
	//for(i=0; i<col; i++)
	//{
	//	(*arrayPtr)[rowIndex][i] = newRow[i];
	//}
	memcpy((*arrayPtr)[rowIndex],  newRow, col*sizeof(double));
}

double* get_bounds(double** visited, int nDim, int tSteps)
{
	/* Calculates the limits of the area explored in the form (x1min, x1max, x2min, x2max, ...).
	'visited' must be a tSteps by nDim array */
	double* bounds;
	int i;
	int j;
	double extreme;

	bounds = (double*)malloc(2*nDim*sizeof(double));
	/* Iterate over collective variables */
	for(i=0; i<nDim; i++)
	{
		/* Go through all the values for this dimension yet visited, keeping track of the smallest */
		extreme=visited[0][i];
		for(j=1; j<tSteps; j++)
		{
			if(visited[j][i] < extreme){extreme = visited[j][i];}
		}
		bounds[2*i] = extreme;
		/* Repeat for the smallest values */
		extreme=visited[0][i];
		for(j=1; j<tSteps; j++)
		{
			if(visited[j][i] > extreme){extreme = visited[j][i];}		
		}
		bounds[2*i+1] = extreme;
		//printf("i=%d, min=%g, max=%g\n", i, bounds[2*i], bounds[2*i+1]);
	}

	return bounds;
}

double** calc_landscape(double* bounds, int n, argstruct* args, gauss* gaussians, int nGauss, int tStep, FILE* f)
{
	int nDim = args->nDim;
	int nPoints; /* Number of grid points in the n-dimensional grid */
	int** indices;
	int i, j, k, indCount;
	int* indTemp;
	double* r;
	double energy, metaEnergy, expFactor, expTotal;
    double w=(*args).w;
    double ds=(*args).sigma;

	/* Set nPoints=n^nDim */
	nPoints=1;
	for(i=0; i<nDim; i++){nPoints *= n;}

	r = (double*)malloc(nDim*sizeof(double));
	indTemp = (int*)malloc(nDim*sizeof(int));
	indices = (int**)malloc(nPoints*sizeof(int*));
	if(indices==NULL)
	{
		puts("Malloc failed in function calc_landscape");
		exit(EXIT_FAILURE);
	}
	for(i=0; i<nPoints; i++)
	{
		indices[i] = (int*)malloc(nDim*sizeof(int));
		if(indices[i]==NULL)
		{
			puts("Malloc failed in function calc_landscape");
			exit(EXIT_FAILURE);
		}
	}
	indCount=0;

	/* Get list of indices for grid points in the nDim-dimensional hypercuboid */
	puts("Recursively calculating energy landscape");
	recurse_landscape(indices, &indCount, indTemp, 0, nDim, n);
	fprintf(f, "ENERGY LANDSCAPE (format: (coords), (true energy), (cumulative energy from metadynamics gaussians)):\n");
	for(i=0; i<nPoints; i++)
	{
		/* Convert the indices to a position vector based on the calculated bounds and desired number of points */
		for(j=0; j<nDim; j++)
		{
			r[j] = (double)indices[i][j]*(bounds[2*j+1] - bounds[2*j])/(double)(n-1) + bounds[2*j];
			fprintf(f, "%g    ", r[j]);
			//printf("%g\n", r[j]);
		}

		/* Get true energy */
		energy = eval_potential(r, args, 0, 0);
		fprintf(f, "%g", energy);

		/* Get metadynamics energy (excluding the true energy) */
		expTotal=0;
		for(j=0; j<nGauss; j++)
		{
			expFactor=0;
			//for(k=0; k<nDim; k++)
			//{
			//	r[k] = (double)indices[i][k]*(bounds[2*k+1] - bounds[2*k])/(double)(n-1) + bounds[2*k];
			//	expFactor += (r[k]-muGauss[j][k])*(r[k]-muGauss[j][k]);
			//}
			expTotal += eval_gaussian(&gaussians[j], r, nDim);
			fprintf(f, "     %g", expTotal);
		}
		fprintf(f, "\n");
	}
	fprintf(f, DIVIDER);
	free(indTemp);
	free(indices);
	free(r);
}

void recurse_landscape(int** indices, int* indCount, int* indTemp, int digit, int nDim, int n)
{
	/* Used for recursive calculation of the energy landscape within an nDim-dimensional grid of dimensions n by n by ... by n
	indices is the table of indices to be filled, indTemp stores a single set of indices temporarily, indCount stores the extent to which indices has been filled */
	int i, j;
	//printf("Recursion level %d\n", digit);
	if(digit<(nDim-1))
	{
		for(i=0; i<n; i++)
		{
			//printf("indTemp[%d] = %d\n", digit, i);
			indTemp[digit] = i;
			recurse_landscape(indices, indCount, indTemp, digit+1, nDim, n);
		}
	}
	else
	{
		for(i=0; i<n; i++)
		{
			indTemp[digit] = i;
			//printf("WRITING TO INDICES at index %d\n", *indCount);
			for(j=0; j<nDim; j++)
			{
				//printf("indices[%d][%d] = %d\n", *indCount, j, indTemp[j]);
				indices[(*indCount)][j] = indTemp[j];
			}
			(*indCount)++;
		}
	}
}

double rand_normal()
{
	static int nCalls=0;
	double rand1, rSq, fac;
	static double rand2;

	/* Initialise Mersenne Twister */
	if(nCalls==0)
	{
		init_genrand(1);
	}
	/* Check if we have a random number left over; if so, return it, and if not, generate two new ones */
	if(nCalls%2==1)
	{
		//printf("RETURNING: %lf\n", rand2);
		return rand2;
	}
	else
	{
		/* Generate two rands in the interval [-1,1]; keep them if they represent a point within the unit circle */
		do
		{
		rand1 = 2*genrand_real1() - 1;
		rand2 = 2*genrand_real1() - 1;
		rSq = rand1*rand1 + rand2*rand2;
		} while(rSq > 1);
		fac = sqrt(-2*log(rSq)/rSq);
		
		/* Return the first rand and keep the second */
		rand2 = rand2*fac;
		//printf("RETURNING: %lf\n", rand1);
		return rand1*fac;
	}
	nCalls++;
}

double eval_gaussian(gauss* gaussian, double* r, int nDim)
{
	/* Evaluates contribution from gaussian to the energy at position r */
	double expFactor=0;
	int i;
	for(i=0; i<nDim; i++)
	{
		expFactor += (r[i]-gaussian->r[i])*(r[i]-gaussian->r[i]);
	}
	expFactor = exp(-expFactor/(2*(gaussian->sig)));
	return (gaussian->h)*expFactor;
}

double* populate_landscape(double* bounds, int n, argstruct* args)
{
	/* Returns an array of doubles constituting a flattened n-dimensional grid on whose points the energy is evaluated.
	bounds = array of bounds [x_min, x_max, y_min, y_max, ...]
	n  = [x-resolution, y-resolution, ...] 
	*/
	int nDim = args->nDim;
	int nPoints; /* Number of grid points in the n-dimensional grid */
	int** indices; /* nPoints-by-nDim array of index 'vectors' */
	int i, j, k, indCount;
	int* indTemp;
	double* delta;
	double* r;
	double* energies;
	double energy, metaEnergy, expFactor, expTotal;
    double w=(*args).w;
    double ds=(*args).sigma;

	/* Compute [dx, dy, ...] */
	delta = (double*)malloc(nDim*sizeof(double));
	for(i=0; i<nDim; i++)
	{
		delta[i] = (double)i*(bounds[2*i+1]-bounds[i])/n;
	}

	/* Set nPoints=n^nDim */
	nPoints=1;
	for(i=0; i<nDim; i++){nPoints *= n;}

	r = (double*)malloc(nDim*sizeof(double));
	indTemp = (int*)malloc(nDim*sizeof(int));
	indices = (int**)malloc(nPoints*sizeof(int*));
	if(indices==NULL)
	{
		puts("Malloc failed in function calc_landscape");
		exit(EXIT_FAILURE);
	}
	for(i=0; i<nPoints; i++)
	{
		indices[i] = (int*)malloc(nDim*sizeof(int));
		if(indices[i]==NULL)
		{
			puts("Malloc failed in function calc_landscape");
			exit(EXIT_FAILURE);
		}
	}
	indCount=0;

	/* Get list of indices for grid points in the nDim-dimensional hypercuboid */
	puts("Recursively calculating energy landscape");
	recurse_landscape(indices, &indCount, indTemp, 0, nDim, n);

	/* Temporarily change the calculation method from tabulated to on-the-fly so that the call to eval_potential doesn't try to interpolate it */


	/* Populate grid at the points specified by 'ind' */
	energies = (double*)malloc(nPoints*sizeof(double));
	for(i=0; i<nPoints; i++)
	{
		/* Convert the indices to a position vector based on the calculated bounds and desired number of points */
		for(j=0; j<nDim; j++)
		{
			r[j] = bounds[2*j] + (double)indices[i][j]*delta[j];
		}

		/* Get true energy at the current position and assign it to the array */
		energy = eval_potential(r, args, 0, 0);
		energies[i] = energy;
	}
	free(indTemp);
	free(indices);
	free(r);
	return energies;
}

int index2D(int indX, int indY, int nX)
{
    /* Converts 2D indices (indX, indY) into the equivalent index for the flattened array.
    2D array must be of x-length nX */
    return indX + nX*indY;
}

double interpolate(argstruct* args, double* r, int order)
{
    /* Interpolates the values at the point 'r'
    order=0 - nearest-neighbour
    order=1 - bilinear
    order=2 - bicubic */

    double energy=0;
    double ind[2]; /* surface-style indices of point 'r' (not constrained to integer values) */

	/* Check that the surface is 2D */
	if(args->nDim != 2)
	{
		puts("Error: 'interpolate' called with non-2D surface. Exiting.");
		exit(EXIT_FAILURE);
	}

    ind[0] = (r[0] - args->bounds[0])/args->deltas[0];
    ind[1] = (r[1] - args->bounds[1])/args->deltas[1];
    if(order==1)
    {
        double tl, tr, bl, br; /* Top-left, top-right, bottom-left, bottom-right.
                                  Bottom and right are in the direction of increasing x and y */
        double x1, x2, y1, y2;

        /* Calculate (x, y) coordinates of grid nodes */
        x1 = args->bounds[0] + (int)ind[0]*args->deltas[0];
        x2 = args->bounds[0] + ((int)ind[0]+1)*args->deltas[0];
        y1 = args->bounds[2] + (int)ind[0]*args->deltas[1];
        y2 = args->bounds[2] + ((int)ind[0]+1)*args->deltas[1];

        tl = args->eGrid[index2D((int)ind[0], (int)ind[1]+1, args->n)];
        tr = args->eGrid[index2D((int)ind[0]+1, (int)ind[1]+1, args->n)];
        bl = args->eGrid[index2D((int)ind[0], (int)ind[1], args->n)];
        br = args->eGrid[index2D((int)ind[0]+1, (int)ind[1], args->n)];

        energy=0;
        energy += bl*(x2-r[0])*(y2-r[1]);
        energy += br*(r[0]-x1)*(y2-r[1]);
        energy += tl*(x2-r[0])*(r[1]-y1);
        energy += tr*(r[0]-x1)*(r[1]-y1);
        energy /= (x2-x1)*(y2-y1);
        return energy;
    }
    else if(order==2)
    {
        puts("Bicubic interpolation not yet implemented.");
        exit(EXIT_FAILURE);
    }
    else if(order==0)
    {
        puts("Nearest-neighbour interpolation not yet implemented");
    }
    else
    {
        puts("Order most be 0, 1 or 2");
    }

    /* Get neighbours */

}