#include "FEM1D.h"

int main(int argc, char *argv[]){
	int Ne=1; int deg = 12; 
	double* intv = new double [2];
	intv[0] = -1.0, intv[1] = 1.0;

	double* coord;

	//integral_gauss* itg = new integral_gauss(deg+1, "_GaussLegendre_");
	integral_gauss* itg = new integral_gauss(deg+1, "_LegendreGaussLobato_");
	FiniteElementSpace* space = new FiniteElementSpace(Ne, deg, intv);
	itg->ComputeIntegralPoints();
	double* pts = itg->GetIntegralPoints();
	Interpolation* itp = new Interpolation(deg, pts);

	//space->setEleNodes(pts);

	cout<<"The DoF is: "<<space->getGlobalDof()<<endl;

	pFEMPP problem = new FEMPP(space, itp, itg);

	problem->solve();

	cout<<setprecision(15)<<"The L2 error is: "<<problem->L2Error()<<endl;

	problem->PlotFem1d(argc, argv);

	// gsl_vector* solution = problem->getSolution();
	// for(int i=0; i<space->getGlobalDof(); i++)
	// 	cout<<setprecision(15)<<gsl_vector_get(solution, i)<<endl;

	delete problem;
	//problem->ComputeRHS();

	// problem->ComputeRefEleStiffMatrix();

	// problem->ComputeRefEleMassMatrix();

	//problem->assemble();

	//problem->imposeDirichletBC();


	// space->IndexMapping(2);
	// int* INDEX = space->getIndex();
	// int eleDof = space->getLocalDof();
	// for(int i=0; i<eleDof; i++){
	// 	cout<<INDEX[i]<<endl;
	// }

	// gsl_matrix* M = space->getConnectM();
	// for(int i=0; i<Ne; i++){
	// 	for(int j=0; j<2; j++){
	// 		cout<<gsl_matrix_get(M, i, j)<<"  ";
	// 	}
	// 	cout<<endl;
	// }

	// coord = space->geteleNodeCoord();
	// for(int i=0; i<4; i++)
	// 	cout<<coord[i]<<endl;
}
