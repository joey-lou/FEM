/* Excercise on static beam strain FEM */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define E 1000
#define L 10

const int N = 100;
const double tx = 25;
const double bx = 10;
const double h = L/float(N-1);
const double M = 200;

double** alloc(double **arr, int m, int n);
void print2D(double **arr,int m, int n);
void gausselimi(double **x, double **A, double **b, int n);
void swaprow(double **arr, int i, int j, int n);
double shapecal(double x, int flag);

int main(int argc, char const *argv[])
{
	// allocate memories to matrices
	double** Kmat;
	Kmat = alloc(Kmat,N,N);
	double** Kele;
	Kele = alloc(Kele,2,2);
	Kele[0][0] = E/h; Kele[0][1] = -E/h; Kele[1][0] = -E/h; Kele[1][1] = E/h;

	int i,j;
	for (i=0;i<N-1;i++){
		for (j=i;j<i+2;j++){
			Kmat[i][j] += Kele[0][j-i];
			Kmat[i+1][j] += Kele[1][j-i];
		}
	}
	Kmat[0][0] = 1; Kmat[0][1] = 0; Kmat[1][0] = 0;

	double** fele;
	fele = alloc(fele,2,1);
	fele[0][0] = bx/2.0*h; fele[1][0] = bx/2.0*h;

	double** fvec;
	fvec = alloc(fvec,N,1);
	for (i=0;i<N/2;i++)
		for (j=0;j<2;j++)
			fvec[i+j][0]+=fele[j][0];
	for (i=N/2;i<N;i++){
		fvec[i][0] = 0;
	}
	fvec[0][0] = 0;
	fvec[1][0] -= Kele[1][0]*fvec[0][0];
	fvec[N-1][0] += tx;
	//print2D(fvec,N,1);
	//print2D(fele,2,1);
	
	double **u;
	u = alloc(u,N,1);
	gausselimi(u, Kmat, fvec, N);
	//print2D(u,N,1);

	// plot
	double **output, elex; int flag = 1, elenum = 1;
	output = alloc(output,M,3);
	//printf("h=%lf\n",h);
	// trouble here:
	for (i=0;i<M;i++){
		output[i][0] = L/(M-1)*i;
		if (output[i][0]-elenum*h>1e-6){
			elenum++;
		}
		elex = output[i][0] - (elenum-1)*h;
		//printf("elenum=%d,elex=%lf\n",elenum, elex);
		output[i][1] = shapecal(elex,1)*u[elenum-1][0];
		output[i][1] += shapecal(elex,2)*u[elenum][0];
		output[i][2] = -E/h*u[elenum-1][0]+E/h*u[elenum][0];
	}

	printf("x\t\t\tu\t\t\tsigma\n");
	print2D(output,M,3);
	FILE *fp1,*fp2;
	fp1 = fopen("plot1","w");
	fp2 = fopen("plot2","w");
	//fprintf(fp1, "x\t\tdisplacement\n");
	//fprintf(fp2, "x\t\tstress\n");
	for (i=0;i<M;i++){
		fprintf(fp1, "%f\t%f\n",output[i][0],output[i][1] );
		fprintf(fp2, "%f\t%f\n",output[i][0],output[i][2] );
	}
	fclose(fp1);
	fclose(fp2);
	return 0;
	
}

double shapecal(double x, int flag){
	double y;
	if (flag==1){
		y = 1-x/h;
	}
	else {
		y = x/h;
	}
	return y;
}

void gausselimi(double **x, double **A, double **b, int n){
	/* assumes left hand matrix is a square matrix */
	double **Gauss,vmax,f;
	int i,j,k,imax;
	Gauss = alloc(Gauss, n, n+1);
	for (i=0;i<n;i++){
		for (j=0;j<n;j++){
			Gauss[i][j] = A[i][j];
		}
		Gauss[i][n] = b[i][0];
	}

	for (k=0;k<n;k++){
		imax = k;
		vmax = Gauss[imax][k];
		for (i=k+1;i<n;i++){
			if (fabs(Gauss[i][k])>vmax){
				vmax = Gauss[i][k], imax = i;
			}
		}

		if(imax != k){
			swaprow(Gauss, k, imax, n+1);
		}

		for(i=k+1;i<n;i++){
			f = Gauss[i][k]/Gauss[k][k];
			for (j=k;j<n+1;j++){
				Gauss[i][j] -= Gauss[k][j]*f;
			}
		}
	}

	// back substitution
	for (i=n-1;i>-1;i--){
		x[i][0] = Gauss[i][n];
		for (j=i+1;j<n;j++){
			x[i][0] -= Gauss[i][j]*x[j][0];
		}
		x[i][0] = x[i][0]/Gauss[i][i];
	}
	
	return;
}

void swaprow(double **arr, int i, int j, int n){
	/* swap the rows in given matrix */
	int k; double temp;
	for (k=0;k<n;k++){
		temp = arr[i][k];
		arr[i][k] = arr[j][k];
		arr[j][k] = temp;
	}
	return;
}

double** alloc(double **arr, int m, int n){
	/* allocate memory for 2D arrays */
	arr = (double**)malloc(sizeof(double*)*m);
	for (int i=0;i<m;i++){
		arr[i] = (double*)malloc(sizeof(double)*n);
	}
	for (int i=0;i<m;i++){
		for (int j=0;j<n;j++){
			arr[i][j] = 0.0;
		}
	}
	return arr;
}

void print2D(double **arr,int m, int n){
	/* visualize matrices */
	int i,j;
	for (i=0;i<m;i++){
		for (j=0;j<n;j++){
			printf("%f\t", arr[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	return;
}
