
#define MAXN WORDSIZE
#define ONE_WORD_SETS 1
#include "Configuration.h"
#include "Bank.h"
#include <iostream>
#include <fstream>

void readPoints(std::istream& f, Configuration& c){
	double d;
	for(int i=0; i<NUM_OF_SPHERES; i++){
		for(int j=0; j<3; j++){
			f>>d;
			c.p(3*i+j) = d;
		}
	}
	for(int i=0; i<NUM_OF_SPHERES-1; i++){
		for(int j=i+1; j<NUM_OF_SPHERES; j++){
			double dist = 0;
			for(int k=0; k<3; k++){
				dist+= pow(c.p(3*i+k)-c.p(3*j+k),2);
			}if(fabs(dist-1)<=1e-6){
				c.addEdge(i,j);
			}
		}
	}
	c.canonize();
}

int main(int argc, char** argv){
	
	Bank bank(true);
	std::cout.precision(16);
	Configuration n[11980];
	Configuration c[11994];
	std::ifstream f;
	f.open ("n12.txt");
	if (!f.is_open()){
		return 1;
	}
	
	
	std::ifstream g;
	g.open ("parsedoutput.txt");
	if (!g.is_open()){
		return 1;
	}
//	readPoints(f, c[0]);
//	readPoints(f, c[1]);

	for(int i=0; i<11980; i++){
		readPoints(f, n[i]);
	}
	for(int i=0; i<11980; i++){
		bank.add(n[i]);
	}
	
	for(int i=0; i<11994; i++){
		readPoints(g, c[i]);
	}
	for(int i=0; i<11994; i++){
		bank.add(c[i]);
	}
	
	
	bank.printDuplicates();
	std::cout<<bank.size()<<std::endl;
	
	return 0;
}













