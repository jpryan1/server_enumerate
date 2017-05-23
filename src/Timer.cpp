#include "Timer.h"
Timer::Timer(){
	memset(times, 0, 4*sizeof(int));
	memset(maxtimes, 0, 4*sizeof(int));
}
void Timer::start(){
	temp = time(0);
}
void Timer::end(int a){
	double t = time(0)-temp;
	times[a] += t;
	maxtimes[a] = std::max(maxtimes[a], t);
}
void Timer::display(){
	std::cout<<"Time spent finding dim - "<<times[0]<<std::endl;
	std::cout<<"Max time spent finding dim - "<<maxtimes[0]<<std::endl;
	
	std::cout<<"Time spent projecting - "<<times[1]<<std::endl;
	std::cout<<"Max time spent projecting - "<<maxtimes[1]<<std::endl;
	
}
