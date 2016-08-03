/*******************************************************************************/
/* This is a simple implementation of Genetic Learning PSO (GLPSO).            */
/* The codes are written in C.                                                 */
/* For any questions, please contact Y. Gong (gongyuejiaoATgmail.com).         */
/* GLPSO_1.0, Edited on March 6th, 2013.                                       */
/*******************************************************************************/


/* add your header files here */

#include<iostream.h>
#include<stdio.h>
#include<stdlib.h>
#include<fstream.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include"function.h"

/* Change any of these parameters to match your needs */

#define PARSIZE 50          /* particle population size */ 
#define FES 300000          /* max. number of function evaluations */ 
#define MAXGENS 3000        /* max. number of generations */
#define TIMES 30            /* number of runs */
#define dimensions 30       /* max. number of problem variables */
#define NSEL 10             /* tournament size in the selection*/    

const double w  = 0.7298;   /* inertia weight */
const double c1 = 1.49618;  /* accelerate coefficient */
const int    sg = 7;        /* stopping gap of generations */
const double pm = 0.01;     /* mutation probability */
double lbound, ubound;      /* the lower and upper bounds of variables */
const double RHO=0.2;
double vmax;                /* max velocity = RHO * (ubound - lbound) */


struct ParticleType         /* a particle in the population */      
{
	double velocity[dimensions];
	double position[dimensions];
	double pbest[dimensions];
	double pbestval;
	int	   stop;
};
ParticleType particle[PARSIZE];  
   
struct ExemplarType         /* the exemplar of a particle */
{
	double x[dimensions];
	double f;
	int stop;			
};
ExemplarType exemplar[PARSIZE];
ExemplarType newexemplar[PARSIZE];


/* control parameters used by the program */

int    generation;
int    fes;		
int    bestpar;
double gbestval;
double foundbest; 

#define betterthan <
#define worsethan >
#define worstvalue 1e300

/* declaration of functions used by this algorithm */

double (*function_name)(double pos[],int dim);
double randval(double low,double high);
double calfit(double pos[],int dim);
void   exemplar_crossover(int num);
void   exemplar_mutation(int num);
void   exemplar_selection(int num);
void   exemplar_tour_selection(int num);
void   exemplar_update(int num);
void   Initialize();
void   Evaluate(int num);
void   Update();
void   Process();


/***********************************************************/
/* Random value generator: generates a value within bounds */
/***********************************************************/

double randval(double low,double high)			
{
	return (double(rand())/RAND_MAX)*(high-low)+low;
}

/***********************************************************/
/* Fitness calculation: calculating the fitness of an      */
/* individual/position and fes ++.                         */
/* this takes a user defined function.                     */
/* "function_name" indicates the object function you are   */
/* testing, which should be initialized before used.       */                          
/***********************************************************/

double calfit(double pos[],int dim)		
{	
	fes++;
	return function_name(pos,dim);
}

/***********************************************************/
/* Exemplar crossover: each dimension of the new exemplar  */
/* is bred by either the uniform crossover of the pbest    */
/* and gbest or randomly copying another particle's pbest. */
/***********************************************************/

void exemplar_crossover(int num)
{
	double r;
	int i, n;
	for(i = 0; i < dimensions; i ++){	
		n = rand()%PARSIZE;
		if(particle[n].pbestval betterthan particle[num].pbestval || num == bestpar){  /* random copy */ 
			newexemplar[num].x[i] = particle[n].pbest[i];	
		}
		else {                                                                         /* uniform crossover */
			r = randval(0,1);   
			newexemplar[num].x[i] = r*particle[num].pbest[i]+(1-r)*particle[bestpar].pbest[i];  
		}		
	}	
}

/***********************************************************/
/* Exemplar mutation: reinitialize a dimension with a      */
/* mutation probability pm                                 */
/***********************************************************/

void exemplar_mutation(int num)
{
	for(int i = 0; i < dimensions; i ++){
		if(randval(0,1) < pm){
			newexemplar[num].x[i] = randval(lbound,ubound);
		}		
	}
}

/***********************************************************/
/* Exemplar selection: new exemplar replaces the old       */
/* one if it has a better fitness value.                   */
/***********************************************************/

void exemplar_selection(int num)
{
	newexemplar[num].f = calfit(newexemplar[num].x,dimensions);
	if(newexemplar[num].f betterthan exemplar[num].f){
		exemplar[num] = newexemplar[num];
		exemplar[num].stop = 0;		
		if(exemplar[num].f betterthan foundbest)
			foundbest=exemplar[num].f;
	}
	else
		exemplar[num].stop++;
}

/***********************************************************/
/* Exemplar tournament Selection: when the exemplar of a   */
/* particle is trapped, employ the 20%M-tournament         */
/* selection to update the particle's exemplar.             */     
/***********************************************************/

void exemplar_tour_selection(int num)
{
	int j;
	int child[NSEL], temp;
	double good;

	
	for(j=0; j<NSEL; j++){
		child[j]=rand()%PARSIZE;
	}
	good=exemplar[child[0]].f;
	temp=child[0];
	for(j=1; j<NSEL; j++){
		if(exemplar[child[j]].f betterthan good){
			good=exemplar[child[j]].f;
			temp=child[j];
		}
	}
	exemplar[num]=exemplar[temp];		
}


/*************************************************************/
/* Exemplar update: this function contains the overall       */
/* for update the exemplar of a particle.                    */
/*************************************************************/

void exemplar_update(int num)
{
	exemplar_crossover(num);
	exemplar_mutation(num);
	exemplar_selection(num);

	if(exemplar[num].stop > sg){
		exemplar[num].stop = 0;
		exemplar_tour_selection(num);		
	}	
}

/***************************************************************/
/* Initialization function: initializes the particles, pbests, */
/* gbest, and exempalrs                                        */
/***************************************************************/

void Initialize()
{
	int i,j;
	fes=0;
	generation=0;
	vmax=RHO*(ubound-lbound);

	/* initialize particles */
	for(i = 0; i < PARSIZE; i ++){                     
		for(j = 0; j < dimensions; j ++){
			particle[i].position[j] = randval(lbound, ubound);
			particle[i].velocity[j] = randval(lbound-particle[i].position[j], ubound-particle[i].position[j]);
			particle[i].pbest[j] = particle[i].position[j];			
		}
		particle[i].pbestval = calfit(particle[i].position,dimensions);
	}

	/* find the gbest */
	gbestval = particle[0].pbestval;
	bestpar = 0;
	for(i = 1; i < PARSIZE; i ++){
		if(particle[i].pbestval betterthan gbestval){  
			gbestval = particle[i].pbestval;
			bestpar = i;
		}
	}
	foundbest = gbestval;	

	/* initialize exemplars */
	for(i = 0; i  <PARSIZE; i ++) exemplar[i].f = worstvalue;		
	for(i = 0; i < PARSIZE; i ++) exemplar_update(i);
}

/***************************************************************/
/* Evaluation function: evaluate the fitness of particle i     */
/* and update the pbest and gbest if a better fitness value    */
/* is found. Otherwise, the stopping generations ++.           */
/***************************************************************/

void Evaluate(int num)
{
	int j;	
	double tmp;	
	
	tmp=calfit(particle[num].position,dimensions);
	if(tmp betterthan particle[num].pbestval){
		particle[num].stop=0;
		for(j=0;j<dimensions;j++)
			particle[num].pbest[j]=particle[num].position[j];
		particle[num].pbestval=tmp;
		if(tmp betterthan particle[bestpar].pbestval){
			bestpar=num;			
			if(tmp betterthan foundbest)
				foundbest=tmp;
		}
	}	
	else
		particle[num].stop++;
}

/***************************************************************/
/* Update function: for each particle, it first updates its    */
/* exemplar, then learns from the exemplar to update position  */
/* and velocity. Note that the velocity and position should be */
/* restricted in the specified range.                          */
/***************************************************************/

void Update()
{
	int i, j;
	double rnd;	
	for(i = 0; i < PARSIZE; i ++){
		if(fes >= FES) break; 
		exemplar_update(i);
		for(j = 0; j < dimensions; j ++){		
			rnd = randval(0,1);
			particle[i].velocity[j] = w * particle[i].velocity[j] + c1 * rnd * (exemplar[i].x[j] - particle[i].position[j]);				
			if(particle[i].velocity[j] < -vmax) particle[i].velocity[j] = -vmax;
			else if(particle[i].velocity[j] > vmax) particle[i].velocity[j] = vmax;
		
			particle[i].position[j] += particle[i].velocity[j];			           
			if(particle[i].position[j] < lbound){
				particle[i].position[j] = lbound;
				particle[i].velocity[j] *= -0.5;
			}
			else if(particle[i].position[j] > ubound){
				particle[i].position[j] = ubound;
				particle[i].velocity[j] *= -0.5;
			}
		}
		Evaluate(i);
	}
}

/*************************************************************/
/* Process: this function contains the overall procedure for */
/* running the genetic learning PSO.                         */
/*************************************************************/

void Process()
{
	
	Initialize();
	while(generation<MAXGENS && fes<FES){		
		Update();
		generation++;
	}
}

/*************************************************************/
/* The main function                                         */
/*************************************************************/
void main()
{

	int i;
	
    function_name = ff8;           /* set the objective function here */
    lbound = -500; ubound = 500;   /* set the variable range */
    for(i = 0; i < TIMES; i++){
        Process();                 /* run the GLPSO algorithm */
		
        /* you can output results using foundbest here, e.g., */
        printf("%g\n", foundbest);
    }

	printf("Successful!\n");
}	