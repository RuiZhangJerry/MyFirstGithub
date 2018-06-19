#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <iostream>
#include <math.h>

using namespace std;

class testfun{
public:
	testfun(){}
	double d(double x){
		//return 1.0/(1.0+25.0*x*x); // Test for Runge's function
		return sin(x);
	}
	double diffd(double x){
		return cos(x);
	}
};

class Interpolation{
	int nOfNode;
	int polydegree;
	double* intNodes;
	double* interval;
public:
	Interpolation(int n){
		this->polydegree = n;
		this->nOfNode = this->polydegree + 1;
		this->intNodes = new double[this->nOfNode];
		for(int i=0; i<this->nOfNode; i++)
			this->intNodes[i] = -1.0+(2.0/this->polydegree)*i;
	}

	Interpolation(int n, double* ncord){
		this->polydegree = n;
		this->nOfNode = this->polydegree + 1;
		this->intNodes = new double;
		this->intNodes = ncord;
	}

	~Interpolation(){
		//delete intNodes;
	}

	void setNodes(double* ncord){
		this->intNodes = new double;
		this->intNodes = ncord;
	}
	void setDegree(int n){
		this->polydegree = n;
	}
	void setInterval(double* intv){
		this->interval = new double[2];
		this->interval = intv;
	}

	double* getNodes(){return this->intNodes; }


	double Phi(int i, double x);
	double GradPhi(int i, double x);

	double mapping(double xi);
	double inversemapping(double x);
	double jacobi(double xi);
	
	double intplt(double* val, double x);
	double intpltGrad(double* val, double x);
};

#endif //INTERPOLATION_H

