#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <iostream>
#include <iomanip>
#include <math.h>

// GSL functions//
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_legendre.h>

using namespace std;


class integral_gauss{
	int degree;
	double* points;
	double* wts;
	const char* pointtype;

public:
	integral_gauss(int n, const char* p){
		this->degree = n-1;
		this->pointtype = p;
		this->points = new double[this->degree+1];
		this->wts = new double[this->degree+1];
	}
	~integral_gauss(){
		delete[] points;
		delete[] wts;
	}

	int getDegree(){return this->degree+1; }

	void ComputeIntegralPoints();
	void ComputeIntegralWeights();

	double getInitial(int i);
	double ComputeKernel(double x);
	void ComputeIntegralPoints_Newton();

	double* GetIntegralPoints(){
		return this->points;
	}
	double* GetIntegralWeights(){
		return this->wts;
	}

	void GetIntegralPW(double** PWs){
		for(int i=0; i<this->degree+1; i++){
			PWs[0][i] = this->points[i];
			PWs[1][i] = this->wts[i];
		}
	}


};

#endif //INTEGRAL_H