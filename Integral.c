#include "Integral.h"


void integral_gauss:: ComputeIntegralPoints(){
	gsl_matrix* M;

	if ((this->pointtype) == "_GaussLegendre_" ){
		gsl_eigen_symmv_workspace * space = gsl_eigen_symmv_alloc(this->degree+1);
		gsl_vector * eval = gsl_vector_calloc(this->degree+1);
		gsl_matrix* evec = gsl_matrix_calloc(this->degree+1, this->degree+1);
		M=gsl_matrix_calloc(this->degree+1, this->degree+1);
		for(int i=0; i<this->degree; i++){
			double b = (i+1)*(i+1)/(4.0*(i+1)*(i+1)-1);
			gsl_matrix_set(M, i, i+1, sqrt(b));
			gsl_matrix_set(M, i+1, i, sqrt(b));
		}
		gsl_eigen_symmv (M, eval, evec, space);
		gsl_eigen_symmv_free (space);
		gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
		for(int i=0;i<this->degree+1;i++){
			this->points[i]=gsl_vector_get(eval, i);
		}
	}

	if ((this->pointtype) == "_LegendreGaussLobato_"){
		this->points[0]=-1.0;  this->points[this->degree]=1.0;
		gsl_eigen_symmv_workspace * space = gsl_eigen_symmv_alloc(this->degree-1);
		gsl_vector * eval = gsl_vector_calloc(this->degree-1);
		gsl_matrix* evec = gsl_matrix_calloc(this->degree-1, this->degree-1);
		M = gsl_matrix_calloc(this->degree-1, this->degree-1);
		for(int i=0; i<this->degree-2; i++){
			double b = (i+1)*(i+3)/((2.0*(i+1)+1)*(2.0*(i+1)+3));
			gsl_matrix_set(M, i, i+1, sqrt(b));
			gsl_matrix_set(M, i+1, i, sqrt(b));
		}
		gsl_eigen_symmv(M, eval, evec, space);
		gsl_eigen_symmv_free (space);
		gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
		for(int i=0;i<this->degree-1;i++){
			this->points[i+1]=gsl_vector_get(eval, i);
		}
	}
	
}

double integral_gauss:: getInitial(int i){
	if(this->pointtype == "_LegendreGaussLobato_"){
		double pi = 4.0*atan(1.0);
		double thetak =  (pi/(4.0*this->degree+2.0))*(4.0*i-1.0);
		double sigmak = 1.0-(this->degree-1)*1.0/(8.0*pow(this->degree*1.0,3));
		sigmak -= 1.0/(384.0*pow(this->degree*1.0,4))*(39.0-28.0/(pow(sin(thetak), 2)));
		sigmak = sigmak * cos(thetak);

		return sigmak;
	}
	if(this->pointtype == "_GaussLegendre_"){
		double pi = 4.0*atan(1.0);
		double thetak =  (pi/(4.0*(this->degree+1)+2.0))*(4.0*i-1.0);
		double sigmak = 1.0-((this->degree+1)-1)*1.0/(8.0*pow((this->degree*1.0+1.0),3));
		sigmak -= 1.0/(384.0*pow((this->degree*1.0+1.0)*1.0,4))*(39.0-28.0/(pow(sin(thetak), 2)));
		sigmak = sigmak * cos(thetak);

		return sigmak;
	}
}

double integral_gauss:: ComputeKernel(double x){
	 double val[this->degree+1], dval[this->degree+1];
	 gsl_sf_legendre_Pl_deriv_array(this->degree, x, val, dval);
	 double v = (1.0-x*x)*dval[this->degree];
	 v = v/(2.0*x*dval[this->degree]- 1.0*this->degree*(this->degree+1)*val[this->degree]);

	 return v;
}

void integral_gauss:: ComputeIntegralPoints_Newton(){
	double eps = 1.0e-10;
	double val[this->degree+1], dval[this->degree+1];

	if(this->pointtype == "_GaussLegendre_"){
		for(int i=0; i<this->degree+1; i++){
			int j = i+1;
			double x0 = this->getInitial(j);
			double error = fabs(gsl_sf_legendre_Pl(this->degree+1, x0));
			while(error > eps){
				gsl_sf_legendre_Pl_deriv_array(this->degree+1, x0, val, dval);
				x0 -= gsl_sf_legendre_Pl(this->degree+1, x0)/dval[this->degree+1];
				error = fabs(gsl_sf_legendre_Pl(this->degree+1, x0));
			}
			this->points[i] = -x0;
		}
	}

	if(this->pointtype == "_LegendreGaussLobato_"){
		this->points[0] = -1.0; this->points[this->degree] = 1.0;
		for(int i=1; i<=this->degree-1; i++){
			double x0 = 0.5*((this->getInitial(i))+this->getInitial(i+1));
			gsl_sf_legendre_Pl_deriv_array(this->degree, x0, val, dval);
			double error = fabs(dval[this->degree]);
			while (error > eps){
				x0 -= this->ComputeKernel(x0);
				gsl_sf_legendre_Pl_deriv_array(this->degree, x0, val, dval);
				error = fabs(dval[this->degree]);
			}
			this->points[i] = -x0;
		}
	}
}


void integral_gauss:: ComputeIntegralWeights(){
	if(this->pointtype == "_LegendreGaussLobato_"){
		for(int i=0; i<this->degree+1; i++){
			double xi = this->points[i];
			double Lx = gsl_sf_legendre_Pl(this->degree, xi);
			this->wts[i] = (1.0/(Lx*Lx))*2.0/(this->degree*(this->degree+1));
		}
	}

	if(this->pointtype == "_GaussLegendre_"){
		double val[this->degree+1], dval[this->degree+1];
		for(int i=0; i<this->degree+1; i++){
			gsl_sf_legendre_Pl_deriv_array(this->degree+1, this->points[i], val, dval);
			double v = (1.0-this->points[i]*this->points[i])*dval[this->degree+1]*dval[this->degree+1];
			this->wts[i] = 2.0/v;
		}
	}
}










