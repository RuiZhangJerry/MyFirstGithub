#include "FEM1D.h"

void FiniteElementSpace:: IndexMapping(int i){
	for(int j=0; j<this->LocalDof; j++){
		INDEX[j] = j+(this->LocalDof-1)*i;
	}
}

int* FiniteElementSpace:: getDEINDEX(){
	this->DEINDEX[0] = 0;
	this->DEINDEX[1] = this->GlobalDof-1;

	return this->DEINDEX;
}

int* FiniteElementSpace:: getMEINDEX(int mybool){
	if(mybool == -1){
		this->MEINDEX[0] = 0;
		return this->MEINDEX;
	}
	if(mybool == 1){
		this->MEINDEX[0] = 1;
		return this->MEINDEX;
	}
}

double FiniteElementSpace:: mapping(int i, double xi){
	double leftcoord = this->eleNodeCoord[i];
	double rightcoord = this->eleNodeCoord[i+1];

	double x = (rightcoord-leftcoord)/2.0 * xi + (rightcoord+leftcoord)/2.0;

	return x;
}

double FiniteElementSpace:: jacobi(int i){
	double leftcoord = this->eleNodeCoord[i];
	double rightcoord = this->eleNodeCoord[i+1];

	return (rightcoord-leftcoord)/2.0;
}

void FiniteElementMethodtoPoissonProblem:: ComputeRHS(){
	double* RHS_ptr = new double [this->femspace->getGlobalDof()];
	for(int i=0; i<this->femspace->getGlobalDof();i++){
		RHS_ptr[i] = 0.0;
	}
	double* eletem = new double [this->femspace->getLocalDof()];

	this->itg_ptr->ComputeIntegralPoints();
	this->itg_ptr->ComputeIntegralWeights();
	double* qp = this->itg_ptr->GetIntegralPoints();
	double* wts = this->itg_ptr->GetIntegralWeights();

	for(int i=0; i<this->femspace->getnOfElement(); i++){
		for(int iii=0; iii<this->femspace->getLocalDof();iii++){
			eletem[iii] = 0.0;
		}
		this->femspace->IndexMapping(i);
		int* INDEX = this->femspace->getIndex();
		//double* nodes = this->itp_ptr->getNodes();
		for(int ii=0; ii<this->femspace->getLocalDof(); ii++){
			double ja = this->femspace->jacobi(i);
			for(int jj=0; jj<this->itg_ptr->getDegree(); jj++){
				double x = this->femspace->mapping(i,qp[jj]);
				eletem[ii] += this->itp_ptr->Phi(ii, qp[jj])*\
								this->Coefficient_ptr[0]->d(x)*wts[jj]*fabs(ja);
			}
			//cout<<eletem[ii]<<endl;
			RHS_ptr[INDEX[ii]] += eletem[ii];
		}
		//cout<<i<<endl;
	}
	for(int i=0; i<this->femspace->getGlobalDof(); i++){
		gsl_vector_set(this->RHS, i, RHS_ptr[i]);
		//cout<<gsl_vector_get(this->RHS, i)<<endl;
	}

	delete[] RHS_ptr;
}

void FiniteElementMethodtoPoissonProblem:: ComputeRefEleStiffMatrix(){
	this->itg_ptr->ComputeIntegralPoints();
	this->itg_ptr->ComputeIntegralWeights();
	double* qp = this->itg_ptr->GetIntegralPoints();
	double* wts = this->itg_ptr->GetIntegralWeights();

	for(int i=0; i<this->femspace->getLocalDof(); i++){
		for(int j=0; j<this->femspace->getLocalDof(); j++){
			double temval = 0.0;
			for(int ii=0; ii<this->itg_ptr->getDegree(); ii++){
				temval += this->itp_ptr->GradPhi(i, qp[ii])*\
						  this->itp_ptr->GradPhi(j,qp[ii])*wts[ii];
			}
			gsl_matrix_set(this->eleStiffMatrix, i, j, temval);
		}
	}

	// for(int i=0; i<this->femspace->getLocalDof(); i++){
	// 	for(int j=0; j<this->femspace->getLocalDof(); j++){
	// 		cout<<gsl_matrix_get(this->eleStiffMatrix, i, j)<<"  ";
	// 	}
	// 	cout<<endl;
	// }
}

// Computation the Refelement Mass Matrices for FEM !-) //
// void FiniteElementMethodtoPoissonProblem:: ComputeRefEleMassMatrix(){
// 	this->itg_ptr->ComputeIntegralPoints();
// 	this->itg_ptr->ComputeIntegralWeights();
// 	double* qp = this->itg_ptr->GetIntegralPoints();
// 	double* wts = this->itg_ptr->GetIntegralWeights();

// 	for(int i=0; i<this->femspace->getLocalDof(); i++){
// 		for(int j=0; j<this->femspace->getLocalDof(); j++){
// 			double temval = 0.0;
// 			for(int ii=0; ii<this->itg_ptr->getDegree(); ii++){
// 				temval += this->itp_ptr->Phi(i, qp[ii])*\
// 						  this->itp_ptr->Phi(j,qp[ii])*wts[ii];
// 			}
// 			gsl_matrix_set(this->eleMassMatrix, i, j, temval);
// 		}
// 	}

// 	// for(int i=0; i<this->femspace->getLocalDof(); i++){
// 	// 	for(int j=0; j<this->femspace->getLocalDof(); j++){
// 	// 		cout<<gsl_matrix_get(this->eleMassMatrix, i, j)<<"  ";
// 	// 	}
// 	// 	cout<<endl;
// 	// }
// }

// Computation the Refelement Mass Matrices for  LGL Spectral Element Method !-) //
void FiniteElementMethodtoPoissonProblem:: ComputeRefEleMassMatrix(){
	this->itg_ptr->ComputeIntegralPoints();
	this->itg_ptr->ComputeIntegralWeights();
	double* qp = this->itg_ptr->GetIntegralPoints();
	double* wts = this->itg_ptr->GetIntegralWeights();

	for(int i=0; i<this->femspace->getLocalDof(); i++){
		double temval = wts[i];
		gsl_matrix_set(this->eleMassMatrix, i, i, temval);
	}

	// for(int i=0; i<this->femspace->getLocalDof(); i++){
	// 	for(int j=0; j<this->femspace->getLocalDof(); j++){
	// 		cout<<gsl_matrix_get(this->eleMassMatrix, i, j)<<"  ";
	// 	}
	// 	cout<<endl;
	// }
}

void FiniteElementMethodtoPoissonProblem:: ComputeLocalMatrix(){
	this->ComputeRefEleStiffMatrix();
	this->ComputeRefEleMassMatrix();
}

void FiniteElementMethodtoPoissonProblem:: assemble(){
	double** GM = new double* [this->femspace->getGlobalDof()];
	for(int i=0; i<this->femspace->getGlobalDof(); i++){
		GM[i] = new double [this->femspace->getGlobalDof()];
		for(int j=0; j<this->femspace->getGlobalDof(); j++){
			GM[i][j] = 0.0;
		}
	}

	this->ComputeLocalMatrix();

	for(int i=0; i<this->femspace->getnOfElement(); i++){
		this->femspace->IndexMapping(i);
		int* INDEX = this->femspace->getIndex();
		double ja = this->femspace->jacobi(i);
		for(int ii=0; ii<this->femspace->getLocalDof(); ii++){
			for(int jj=0; jj<this->femspace->getLocalDof(); jj++){
				GM[INDEX[ii]][INDEX[jj]] += gsl_matrix_get(this->eleStiffMatrix,
											ii, jj)*fabs(1.0/ja)*this->coefa;
				GM[INDEX[ii]][INDEX[jj]] += gsl_matrix_get(this->eleMassMatrix,
									        ii, jj)*fabs(ja)*this->coefb;
			}
		}
	}

	for(int i=0; i<this->femspace->getGlobalDof(); i++){
		for(int j=0; j<this->femspace->getGlobalDof(); j++){
			gsl_spmatrix_set(this->gloStiffMatrix, i, j, GM[i][j]);
		}
	}

	delete[] GM;
}

void FiniteElementMethodtoPoissonProblem:: imposeDirichletBC(){
	double* Vertical = this->femspace->getBoundaryCoord();
	gsl_vector* pp = gsl_vector_calloc(this->femspace->getGlobalDof());
	int* DEINDEX = this->femspace->getDEINDEX();
	for(int i=0; i<this->femspace->getnOfBoundaryPoints(); i++){
		double val = Coefficient_ptr[1]->d(Vertical[i]);
		gsl_vector_set(pp, DEINDEX[i], val);
	}
	gsl_matrix* A = gsl_matrix_calloc(this->femspace->getGlobalDof(),\
		                              this->femspace->getGlobalDof());
	
	gsl_spmatrix_sp2d(A, this->gloStiffMatrix);
	gsl_blas_dgemv(CblasNoTrans, -1.0, A, pp, 1.0, this->RHS);

	for(int i=0; i<this->femspace->getnOfBoundaryPoints(); i++){
		for(int ii=0; ii<this->femspace->getGlobalDof(); ii++){
			gsl_spmatrix_set(this->gloStiffMatrix, DEINDEX[i], ii, 0.0);
			gsl_spmatrix_set(this->gloStiffMatrix, ii, DEINDEX[i], 0.0);
		}
		gsl_spmatrix_set(this->gloStiffMatrix, DEINDEX[i], DEINDEX[i], 1.0);
		gsl_vector_set(this->RHS, DEINDEX[i], Coefficient_ptr[1]->d(Vertical[i]));
	}

	gsl_vector_free(pp); 
	gsl_matrix_free(A);
}

void FiniteElementMethodtoPoissonProblem:: imposeMixedBC(){
	double* Vertical = this->femspace->getBoundaryCoord();
	//int* MEINDEX = this->femspace->getMEINDEX(1);
	gsl_vector* pp = gsl_vector_calloc(this->femspace->getGlobalDof());
	int* DEINDEX = this->femspace->getDEINDEX();
	for(int i=0; i<1; i++){
		double val = Coefficient_ptr[1]->d(Vertical[i]);
		gsl_vector_set(pp, DEINDEX[i], val);
	}
	gsl_matrix* A = gsl_matrix_calloc(this->femspace->getGlobalDof(),\
		                              this->femspace->getGlobalDof());
	
	gsl_spmatrix_sp2d(A, this->gloStiffMatrix);
	gsl_blas_dgemv(CblasNoTrans, -1.0, A, pp, 1.0, this->RHS);

	for(int i=0; i<1; i++){
		for(int ii=0; ii<this->femspace->getGlobalDof(); ii++){
			gsl_spmatrix_set(this->gloStiffMatrix, DEINDEX[i], ii, 0.0);
			gsl_spmatrix_set(this->gloStiffMatrix, ii, DEINDEX[i], 0.0);
		}
		gsl_spmatrix_set(this->gloStiffMatrix, DEINDEX[i], DEINDEX[i], 1.0);
		gsl_vector_set(this->RHS, DEINDEX[i], Coefficient_ptr[1]->d(Vertical[i]));
	}
	gsl_vector_set(this->RHS, DEINDEX[1], 
		gsl_vector_get(RHS, DEINDEX[1])+Coefficient_ptr[2]->d(Vertical[1]));

	gsl_vector_free(pp); 
	gsl_matrix_free(A);
}

void FiniteElementMethodtoPoissonProblem:: imposeRobinBC(){
	double* Vertical = this->femspace->getBoundaryCoord();
	int* DEINDEX = this->femspace->getDEINDEX();
	
	gsl_spmatrix_set(this->gloStiffMatrix, DEINDEX[0], DEINDEX[0], 
		gsl_spmatrix_get(gloStiffMatrix, DEINDEX[0],DEINDEX[0])+this->beta[0]);

	gsl_spmatrix_set(this->gloStiffMatrix, DEINDEX[1], DEINDEX[1], 
		gsl_spmatrix_get(gloStiffMatrix, DEINDEX[1],DEINDEX[1])+this->beta[1]);

	gsl_vector_set(this->RHS, DEINDEX[0], \
		gsl_vector_get(RHS, DEINDEX[0])+Coefficient_ptr[3]->d(Vertical[0]));
	gsl_vector_set(this->RHS, DEINDEX[1], \
		gsl_vector_get(RHS, DEINDEX[1])+Coefficient_ptr[3]->d(Vertical[1]));
}

void FiniteElementMethodtoPoissonProblem:: imposeBoundaryCondition(const char* BCtype){
	if(BCtype == "_Dirichlet_"){
		this->imposeDirichletBC();
	}
	if(BCtype == "_Mixed_"){
		this->imposeMixedBC();
	}
	if(BCtype == "_Robin_"){
		this->imposeRobinBC();
	}
}

void FiniteElementMethodtoPoissonProblem:: solve(){
	this->ComputeRHS();
	this->ComputeLocalMatrix();
	this->assemble();
	//this->imposeBoundaryCondition("_Dirichlet_");
	this->imposeBoundaryCondition("_Mixed_");
	
	//gsl_gmres_solver setup//
	const double tol = 1.0e-14; /* solution relative tolerance */
	const size_t max_iter = 50000;
	size_t iter = 0;
	const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
	gsl_splinalg_itersolve *work =\
		gsl_splinalg_itersolve_alloc(T, this->femspace->getGlobalDof(), 0);
	int state;
	do{
		state = gsl_splinalg_itersolve_iterate(this->gloStiffMatrix, this->RHS, tol, \
									this->Solution, work);
	}
	while(state == GSL_CONTINUE && ++iter<max_iter );


	gsl_splinalg_itersolve_free(work);
}

double FiniteElementMethodtoPoissonProblem:: L2Error(){
	int refine = 50;
	double* refnodes = new double[refine+1];
	for(int i=0; i<refine+1; i++){
		refnodes[i] = -1.0+2/refine*i;
	}

	double* eleError = new double[this->femspace->getnOfElement()];
	for(int i=0; i<this->femspace->getnOfElement(); i++){
		eleError[i] = 0.0;
		for(int j=0; j<refine+1; j++){
			double sum = 0.0;
			double x = this->femspace->mapping(i,refnodes[j]);
			for(int jj=0; jj<this->femspace->getLocalDof(); jj++){
				sum += this->itp_ptr->Phi(jj, refnodes[j])*\
				gsl_vector_get(this->Solution, jj+(this->femspace->getLocalDof()-1)*i);
			}
			eleError[i] += pow((Coefficient_ptr[1]->d(x)-sum),2);
		}
		eleError[i] = sqrt(eleError[i]);
		//cout<<eleError[i]<<endl;
	}

	double error = 0.0;
	for(int i=0; i<this->femspace->getnOfElement(); i++){
		error += eleError[i]*eleError[i];
	}

	return sqrt(error);

	delete[] eleError;
}

void myDisplay(){

	int Ne=1; int deg = 12; int refine = 1000;
	double* intv = new double [2];
	intv[0] = -1.0, intv[1] = 1.0;
	double pi = 4.0*atan(1.0);

	glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_LINES);
    	glColor3d(1.0,1.0,1.0);
     	glVertex2d(intv[0]-0.1, 0.0);
     	glVertex2d(intv[1]+0.1, 0.0);
     	glVertex2d(0.0, intv[0]-0.1);
     	glVertex2d(0.0, 10.0);
    glEnd();

	double* coord;
	FiniteElementSpace* space = new FiniteElementSpace(Ne, deg, intv);

	//integral_gauss* itg = new integral_gauss(deg+1, "_GaussLegendre_");
	integral_gauss* itg = new integral_gauss(deg+1, "_LegendreGaussLobato_");
	
	itg->ComputeIntegralPoints();
	double* pts = itg->GetIntegralPoints();
	Interpolation* itp = new Interpolation(deg, pts);
	//space->setEleNodes(pts);

	pFEMPP problem = new FEMPP(space, itp, itg);

	problem->solve();

	gsl_vector* val = problem->getSolution();

	Coefficient* exsol = new DirichletBC();
	Coefficient* exsoldiff = new NuemannBC();

	GLdouble* Val = new GLdouble[space->getGlobalDof()];
	GLdouble* eleVal = new GLdouble[space->getLocalDof()];
	GLdouble* elementNodes = new GLdouble[Ne+1];
	GLdouble* x = new GLdouble[refine+1];
	GLint i, j;

	for(i=0; i<space->getGlobalDof(); i++){
		Val[i] = gsl_vector_get(val, i);
	}
	double* elenode = new double[Ne+1];
	double* eleinterval = new double[2];
	elenode = space->geteleNodeCoord();
	for(i=0; i<Ne+1; i++){
		elementNodes[i] = elenode[i];
	}
	for(i=0; i<refine+1; i++){
		x[i] = -1.0+2.0/refine*i*1.0;
	}

	for(i=0; i<space->getnOfElement(); i++){
		space->IndexMapping(i);
		int* INDEX = space->getIndex();
		for(j=0; j<space->getLocalDof(); j++){
			eleVal[j] = Val[INDEX[j]];
		}
		eleinterval[0] = elenode[i];
		eleinterval[1] = elenode[i+1];
		itp->setInterval(eleinterval);

		glLineWidth (1.5);
		glBegin(GL_LINE_STRIP);
     		glColor3d(1.0,0.0,0.0);
			for(j=0; j<refine+1; j++){
				GLdouble xj = space->mapping(i, x[j]);
				GLdouble inv = itp->intplt(eleVal, x[j]);
				//cout<<xj<<endl;
				glVertex2d(xj, inv);
			}
		glEnd();

	 	glLineWidth (1.0);
        glLineStipple (1, 0x0F0F);  
		glBegin(GL_LINES);
     		glColor3d(1.0,1.0,0.0);
			for(j=0; j<refine+1; j++){
				GLdouble xj = space->mapping(i, x[j]);
				GLdouble exv = exsol->d(xj);
				//cout<<xj<<endl;
				glVertex2d(xj, exv);
			}
		glEnd();

		glLineWidth (1.5);
		glBegin(GL_LINE_STRIP);
     		glColor3d(1.0,0.0,0.0);
			for(j=0; j<refine+1; j++){
				GLdouble xj = space->mapping(i, x[j]);
				GLdouble inv = itp->intpltGrad(eleVal, x[j]);
				glVertex2d(xj, inv/pi);
			}
		glEnd();

		glLineWidth (1.5);
		glBegin(GL_LINE_STRIP);
     		glColor3d(1.0,1.0,0.8);
			for(j=0; j<refine+1; j++){
				GLdouble xj = space->mapping(i, x[j]);
				GLdouble exv = exsoldiff->d(xj);
				glVertex2d(xj, exv/pi);
			}
		glEnd();

		glFlush();
	}

}

void FiniteElementMethodtoPoissonProblem:: PlotFem1d(int argc, char *argv[]){
	
	glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(200, 200);
    glutInitWindowSize(800, 1200);


    glutCreateWindow("My First OpenGL Plot for FEM1D !-)");
    glClearColor(0.5, 0.5, 0.5, 0.0);    
    //设置投影类型：正投影    
    glMatrixMode(GL_PROJECTION||GL_MODELVIEW);    
    //观察参数：x坐标值从0到200，y是从0到150    
    gluOrtho2D(-1.1, 1.1, -1.1, 1.1);  
    glutDisplayFunc(myDisplay);

    // glutCreateWindow("My First OpenGL Plot for FEM1D !-)");
    // glClearColor(0.5, 0.5, 0.5, 0.0);    
    // //设置投影类型：正投影    
    // glMatrixMode(GL_PROJECTION||GL_MODELVIEW);    
    // //观察参数：x坐标值从0到200，y是从0到150    
    // gluOrtho2D(-1.1, 1.1, -1.1, 1.1);  
    // glutDisplayFunc(myDisplay);

   
    glutMainLoop();
}








