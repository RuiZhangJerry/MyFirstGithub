#include "interpolation.h"

//double Interpolation:: Phi(int i, double x){
//	double phival = 1.0;
//	for(int ii=0; ii<i; ii++){
//		phival = phival * (x-this->intNodes[ii])/(this->intNodes[i]-this->intNodes[ii]);
//	}
//	for(int ii=i+1; ii<this->nOfNode; ii++){
//		phival = phival * (x-this->intNodes[ii])/(this->intNodes[i]-this->intNodes[ii]);
//	}
//
//	return phival;
//}
//
//double Interpolation:: GradPhi(int i, double x){
//	double tem = 1.0;
//	for(int ii=0; ii<i; ii++){
//		tem = tem * (this->intNodes[i]-this->intNodes[ii]);
//	}
//	for(int ii=i+1; ii<this->nOfNode; ii++){
//		tem = tem * (this->intNodes[i]-this->intNodes[ii]);
//	}
//	tem = 1.0/tem;
//	
//	double sum = 0.0;
//	for(int k = 0; k<i; k++){
//		double tem2 = 1.0;
//		for(int ii=0; ii<i; ii++){	
//			if(ii == k) continue; 
//			tem2 = tem2 * (x-this->intNodes[ii]);
//		}
//		for(int ii=i+1; ii<this->nOfNode; ii++){
//			if(ii == k) continue;
//			tem2 = tem2 * (x-this->intNodes[ii]);	
//		}
//		sum += tem2;
//	}
//	for(int k = i+1; k<this->nOfNode; k++){
//		double tem2 = 1.0;
//		for(int ii=0; ii<i; ii++){
//			if(ii == k) continue;
//			tem2 = tem2 * (x-this->intNodes[ii]);
//		}
//		for(int ii=i+1; ii<this->nOfNode; ii++){
//			if(ii == k) continue;
//			tem2 = tem2 * (x-this->intNodes[ii]);	
//		}
//		sum += tem2;
//	}
//	
//	return tem*sum;
//} 

double Interpolation:: Phi(int i, double x){
	double phival = 1.0;
	for(int ii=0; ii<this->nOfNode; ii++){
		if(ii == i) continue;
		phival = phival * (x-this->intNodes[ii])/(this->intNodes[i]-this->intNodes[ii]);
	}

	return phival;
}

double Interpolation:: GradPhi(int i, double x){
	double tem = 1.0;
	for(int ii=0; ii<this->nOfNode; ii++){
		if(ii == i) continue;
		tem = tem * (this->intNodes[i]-this->intNodes[ii]);
	}

	tem = 1.0/tem;
	
	double sum = 0.0;
	for(int k = 0; k<this->nOfNode; k++){
		if(k == i) continue;
		double tem2 = 1.0;
		for(int ii=0; ii<this->nOfNode; ii++){	
			if(ii == k || ii == i) continue; 
			tem2 = tem2 * (x-this->intNodes[ii]);
		}
		sum += tem2;
	}
	
	return tem*sum;
} 

double Interpolation:: mapping(double xi){
	return (this->interval[1]-this->interval[0])/2.0 * xi + (this->interval[1]+this->interval[0])/2.0;
}

double Interpolation:: jacobi(double xi){
	return 2.0/(this->interval[1]-this->interval[0]);
}

double Interpolation:: intplt(double* val, double x){
	double sum = 0.0;
	for(int i=0; i<this->nOfNode; i++){
		sum += val[i]*this->Phi(i, x);
	}
	
	return sum;
}

double Interpolation::intpltGrad(double* val, double x){
	double sum = 0.0;
	for(int i=0; i<this->nOfNode; i++){
		sum += val[i]*this->GradPhi(i, x)*fabs(this->jacobi(x));
	}
	
	return sum;
}










