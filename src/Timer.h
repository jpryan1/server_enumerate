#ifndef  _TIMER_H_    /* only process this file once */
#define  _TIMER_H_

#include <iostream>
#include <ctime>
#include <string.h>
class Timer{
public:
	Timer();
	void start();
	void end(int a);
	void display();
private:
	double times[4];
	double maxtimes[4];
	double temp;
	
};




#endif




