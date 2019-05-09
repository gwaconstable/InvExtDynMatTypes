#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex>
#include <random>

//#include <libiomp/omp.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

// Gillespie Simulation of mating type dynamics described in Eqs. (1-2) of "Invasion and extinction dynamics of mating types under facultative sexual reproduction" 

int main(int argc, char * argv[]) {
    
	// Import parameters for simulation
    
	int output = atof(argv[1]);				// Define what data to output: [0,1,2] - 0=trajectory 1=fixation/extinction properties
	int paramSet = atof(argv[2]);				// Label the parameters for outputting data
	int run = atof(argv[3]);				// Label the run for outputting data 
	int initialCondSetting = atof(argv[4]);    		// Determine which initial conditions to use: Flat distribution (0)
   	double sampleRate = atof(argv[5])*atof(argv[9]);	// Timestep in between each output
    	double minTimeStep = atof(argv[6])*atof(argv[9]);	// Minimum timestep
    	double maxTimeStep = atof(argv[7])*atof(argv[9]);	// Maximum timestep
    
	int Mmax = atof(argv[8]);				// Maximum number of mating types allowed in memory 
	int N = atof(argv[9]);					// Initial number of individuals
	double m = atof(argv[10]);				// Mutation rate per "generation"
	double c = atof(argv[11]);				// Rate of clonal reproduction

	int M0 = atof(argv[12]);				//Initial number of mating types

	//double m = mg/(double)N;				// Mutation rate per timestep

	srand ( time(NULL) );
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	double r1;
	double r2;
	double sum;
	int z;

	int M = M0;

	// Initialise even distribution of number of mating types of i-th class, n[i] 
	std::vector<int> n;
	n.reserve(Mmax);
	n.push_back( N );
	for(int i=1; i<M0; i++)
	{
		n.push_back( round( N/(double)M0 ) );
		n[0] = n[0] - n[i];
	}

	std::vector<double> T;
	T.reserve(Mmax*Mmax+1);
	for(int i=0; i<M0*M0+1; i++)
	{
		T.push_back( 0 );
	}
	double TSUM;

	double t=0;
	double dT=0;
	double tau;

	int endProg = 0;

	while( endProg == 0 )
	{

		// Define all probability transition rates as per Eqs.(1-2) of paper
		TSUM=0;
		for(int i=0; i< M; i++)
		{
			for(int j=0; j<M; j++)
			{
				if( i!= j )
				{
					T[i*M+j] = ((1-c)/(double)2)*(n[i]/(double)N)*( 1 - (n[i]/(double)N) )*(n[j]/(double)N) + c*(n[i]/(double)N)*(n[j]/(double)N);
					TSUM = TSUM + T[i*M+j];
				}
				else
				{
					T[i*M+j] = 0;
					TSUM = TSUM + T[i*M+j];
				}
			}
		}
		T[M*M] = m;
		TSUM = TSUM + T[M*M];
		

		//generate random numbers for choosing time of next transition and which of the transitions it is
	        r1 = 0;
 	       	while(r1 == 0)
        	{
        	    r1 = dist(engine);
        	}
        	r2 = 0;
        	while(r2 == 0)
        	{
        	    r2 = dist(engine);
        	}

		// choose time to next reaction
		dT = (1/TSUM)*(log(1/r1));

		sum = 0;
		z = -1;

		// choose which reaction (z) occurs
		while(r2*TSUM >= sum)
		{
			z++;
		  	sum+=T[z];
		}

		// Note that the amount of times that we cycle through the data depends on dT (i.e. time until next mutation) but we dont want to output the updated state and time repeatedly

		while( (t+dT) > tau )
    		{
			if ( (t+dT) <= (maxTimeStep) &&  t >= minTimeStep   )
        		{
				if( output==1 )
				{
					cout << tau/(double)N <<" "<< M <<" ";
	            			for(int i=0; i< M; i++)
	            			{
	                			cout << n[i] <<" ";
	            			}
	            			cout <<" "<< endl;
				}
			}
	                else if( t > maxTimeStep )
	                {
				endProg = 1;
	                }
        
	               tau += sampleRate;
		}

        	// END: OUTPUT DATA -------------------------------------------------------------------------------------
		

		// Increment time and execute reaction
		t = t + dT;		
		if( z < M*M )
		{

			// Convert from z-th reaction to which type (i) increases by one and which type (j) decreases by one

			int zi=z/M;
			int zj=z%M;

			// Implement reaction	
			n[zi] = n[zi] + 1;
			n[zj] = n[zj] - 1;

			// If extinction occurs, remove jth class from n and all reactions involving type j
			if( n[zj] == 0 )
			{
				n.erase( n.begin() + zj );				// Erase appropriate element of n
				T.erase( T.begin() + 0, T.begin() + 2*M - 1 ); 		// Erase last column and row of Tij Matrix
				M = M - 1;
				if( output == 0 )					// If output = 0 output extinction time 
				{
					cout << t/(double)N << endl;
					//cout << "Working" << endl;
					return 0;
				}
			}
		}
		else if( z==M*M )
		{
			// If mutation happens, we need to pick type to mutate

			//generate random numbers
			r2 = 0;
			while(r2 == 0)
			{
			    r2 = dist(engine);
			}

			sum = 0;
			int zj = -1;

			while(r2*N >= sum)
			{
				zj++;
			  	sum+=n[zj];
			}

			n.push_back( 1 );
			n[zj] = n[zj] - 1;

			// If n[zj] is now zero, I need to delete n[zj], but as I've added a new type I don't need to edit M or T
			if( n[zj] == 0 )
			{
				n.erase( n.begin() + zj );				// Erase appropriate element of n
			}
			else
			{
				// if n[zj] is not zero then I need to add a bunch of transitions, but only for Tij, not for mutations
				for( int i = 0; i < 2*M + 1; i++ )
				{
					T.push_back( 0 );
				}
				M = M + 1;
			}
	
		}



		

	}
  	

	return 0;

}
