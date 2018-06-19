#ifndef FEM1D_H
#define FEM1D_H

#include "interpolation.h"
#include "Integral.h"

#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_blas.h>

#include <GL/glut.h>

class Coefficient{
public:
	Coefficient(){}
	virtual double d(double x){
		return 0.0;
	}
};

class RightHandSide:public Coefficient{
public:
	RightHandSide(){}
	virtual double d(double x){
		double pi = 4.0*atan(1.0);
		return pi*pi*sin(pi*x)+sin(pi*x);
	}
};

class DirichletBC:public Coefficient{
public:
	DirichletBC(){}
	virtual double d(double x){
		double pi = 4.0*atan(1.0);
		return sin(pi*x) ;
	}
};

class NuemannBC:public Coefficient{
public:
	NuemannBC(){}
	virtual double d(double x){
		double pi = 4.0*atan(1.0);
		return pi*cos(pi*x) ;
	}
};

class RobinBC:public Coefficient{
public:
	RobinBC(){}
	virtual double d(double x){
		double pi = 4.0*atan(1.0);
		if(x == 0.0)
			return -pi*cos(pi*x)+sin(pi*x) ;
		if(x == 1.0)
			return  pi*cos(pi*x)+sin(pi*x) ;
	}
};

// class RightHandSide:public Coefficient{
// public:
// 	RightHandSide(){}
// 	virtual double d(double x){
// 		double pi = 4.0*atan(1.0);
// 		return -pi*pi*exp(pi*x)+exp(pi*x);
// 	}
// };

// class DirichletBC:public Coefficient{
// public:
// 	DirichletBC(){}
// 	virtual double d(double x){
// 		double pi = 4.0*atan(1.0);
// 		return exp(pi*x) ;
// 	}
// };

// class NuemannBC:public Coefficient{
// public:
// 	NuemannBC(){}
// 	virtual double d(double x){
// 		double pi = 4.0*atan(1.0);
// 		return pi*exp(pi*x) ;
// 	}
// };

// class RightHandSide:public Coefficient{
// public:
// 	RightHandSide(){}
// 	virtual double d(double x){
// 		double pi = 4.0*atan(1.0);
// 		return 0.0;
// 	}
// };

// class DirichletBC:public Coefficient{
// public:
// 	DirichletBC(){}
// 	virtual double d(double x){
// 		double pi = 4.0*atan(1.0);
// 		return exp(x) ;
// 	}
// };

// class NuemannBC:public Coefficient{
// public:
// 	NuemannBC(){}
// 	virtual double d(double x){
// 		double pi = 4.0*atan(1.0);
// 		return exp(x) ;
// 	}
// };

// class RobinBC:public Coefficient{
// public:
// 	RobinBC(){}
// 	virtual double d(double x){
// 		double pi = 4.0*atan(1.0);
// 		if(x == 0.0)
// 			return 0.0 ;
// 		if(x == 1.0)
// 			return 2.0*exp(x);
// 	}
// };


class FiniteElementSpace{
	int nOfElement;
	int LocalDof;
	int GlobalDof;
	double* eleNodeCoord, *Vertical;
	int* INDEX, *DEINDEX, *MEINDEX;
	gsl_matrix* ConnectM;

public:
	FiniteElementSpace(int NE, int deg, double* interval){
		this->nOfElement = NE;
		this->LocalDof = deg+1;
		this->GlobalDof = (deg+1)*NE-(NE-1);
		this->ConnectM = gsl_matrix_calloc(NE, 2);
		for(int i=0; i<NE; i++){
			gsl_matrix_set(this->ConnectM,i,0,i+1);
			gsl_matrix_set(this->ConnectM,i,1,i+2);
		}
		this->eleNodeCoord = new double[NE+1];
		for(int i=0; i<NE+1; i++){
			this->eleNodeCoord[i] = interval[0]+(interval[1]-interval[0])/NE * i;
		}
		this->Vertical = new double[2];
		for(int i=0; i<2; i++){
			this->Vertical[i] = interval[i];
		}
		this->INDEX = new int[this->LocalDof];
		this->DEINDEX = new int[2];
		this->MEINDEX = new int[1];
	}

	FiniteElementSpace(int NE, int deg, double* interval, double* pts){
		this->nOfElement = NE;
		this->LocalDof = deg+1;
		this->GlobalDof = (deg+1)*NE-(NE-1);
		this->ConnectM = gsl_matrix_calloc(NE, 2);
		for(int i=0; i<NE; i++){
			gsl_matrix_set(this->ConnectM,i,0,i+1);
			gsl_matrix_set(this->ConnectM,i,1,i+2);
		}
		this->eleNodeCoord = new double[NE+1];
		for(int i=0; i<NE+1; i++){
			this->eleNodeCoord[i] = pts[i];
		}
		this->Vertical = new double[2];
		for(int i=0; i<2; i++){
			this->Vertical[i] = interval[i];
		}
		this->INDEX = new int[this->LocalDof];
		this->DEINDEX = new int[2];
		this->MEINDEX = new int[1];
	}
	~FiniteElementSpace(){
		gsl_matrix_free(this->ConnectM);
		delete[] this->eleNodeCoord;
		delete[] this->Vertical;
		delete[] this->INDEX;
		delete[] this->DEINDEX;
	}

	void setEleNodes(double* points){
		this->eleNodeCoord = points;
	}

	gsl_matrix* getConnectM(){return this->ConnectM; }
	double* geteleNodeCoord(){return this->eleNodeCoord; }
	int getLocalDof(){return this->LocalDof; }
	int getGlobalDof(){return this->GlobalDof; }
	int getnOfElement(){return this->nOfElement; }

	void IndexMapping(int i);
	int* getIndex(){return INDEX; }

	int getnOfBoundaryPoints(){return 2; }
	double* getBoundaryCoord(){return this->Vertical; };
	int* getDEINDEX();
	int* getMEINDEX(int mybool); // mybool == -1: left vertical; 
	//                              mybool == 1: right vertical !-)

	double mapping(int i, double xi);
	double jacobi(int i);
};

class FiniteElementMethodtoPoissonProblem{
	double coefa, coefb;
	Interpolation* itp_ptr;
	integral_gauss* itg_ptr;
	FiniteElementSpace* femspace;

	double* beta;

	gsl_matrix* eleStiffMatrix, *eleMassMatrix;
	gsl_spmatrix* gloStiffMatrix, *gloMassMatrix; 
	gsl_vector* RHS, *Solution;

	Coefficient** Coefficient_ptr;

public:
	FiniteElementMethodtoPoissonProblem(FiniteElementSpace* fem, 
										Interpolation* ip, integral_gauss* quad){
		this->itp_ptr = ip;
		this->itg_ptr = quad;
		this->femspace = fem;

		this->coefa = 1.0;
		this->coefb = 1.0;

		this->beta = new double [2];
		this->beta[0] = 1.0;
		this->beta[1] = 1.0;

		int n = fem->getLocalDof();
		int N = fem->getGlobalDof();


		this->gloStiffMatrix = gsl_spmatrix_alloc(N, N); 

		this->eleStiffMatrix = gsl_matrix_calloc(n, n);
		this->eleMassMatrix = gsl_matrix_calloc(n, n);

		this->RHS = gsl_vector_calloc(N);
		this->Solution = gsl_vector_calloc(N);

		this->Coefficient_ptr = new Coefficient* [4];
		this->Coefficient_ptr[0] = new  RightHandSide();
		this->Coefficient_ptr[1] = new  DirichletBC();
		this->Coefficient_ptr[2] = new  NuemannBC();
		this->Coefficient_ptr[3] = new  RobinBC();
	}
	~FiniteElementMethodtoPoissonProblem(){
		gsl_vector_free(this->RHS);
		gsl_vector_free(this->Solution);

		gsl_matrix_free(this->eleStiffMatrix);
		gsl_matrix_free(this->eleMassMatrix);

		gsl_spmatrix_free(this->gloStiffMatrix);
		
		delete[] this->Coefficient_ptr;
	}

	void ComputeRHS();

	void ComputeRefEleStiffMatrix();
	void ComputeRefEleMassMatrix();
	void ComputeLocalMatrix();

	void assemble();

	void imposeDirichletBC();
	void imposeMixedBC();
	void imposeRobinBC();
	void imposeBoundaryCondition(const char* BCtype);// Choose the type of 
													 // Boundary Condition:
		                                             // "_Dirichlet_", "_Mixed_" and
												     // "_Robin_" !-)
	void solve();

	double L2Error();

	gsl_vector* getSolution(){
		return this->Solution;
	}

	void PlotFem1d(int argc, char *argv[]);

};
void myDisplay();

typedef FiniteElementMethodtoPoissonProblem FEMPP;
typedef FiniteElementMethodtoPoissonProblem* pFEMPP;

#endif //FEM1D_H















