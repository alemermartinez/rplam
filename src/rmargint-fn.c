// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>

#include <R.h>
#include "rmargint.h"

#define ZERO 1e-10
#define TOL_INVERSE ZERO
#define min(a,b) ((a)<(b)?(a):(b))



double psi_huber_w( double r, double k) {
double sal, aux1, aux2;
aux1 = fabs(r);
aux2 = k/aux1;
// sal = min(aux2, 1.0);
if( aux2 < 1 ) { sal = aux2; }
else { sal = 1;}
return( sal );
}


double psi_tukey_w( double r, double k) {
double sal=0.0, aux1;
aux1 = fabs(r/k);
if( aux1 < 1.0 ) {
	sal = (1+aux1) * (1+aux1) * (1-aux1) * (1-aux1);
};
return( sal );
}

//Nucleo de Epanechnikov

double kepan( double b ) {
double a, aux;
a = 0.75 * (1 - (b*b) );
aux = fabs(b);
if( aux < 1 ) { return(a); }
else { return(0); }
}

double kernel4( double b ) {
double a, aux;
a = 0.46875 * (1 - (b*b) ) * (3 - 7*(b*b) ); // 15/32
aux = fabs(b);
if( aux < 1 ) { return(a); }
else { return(0); }
}

double kernel6( double b ) {
double a, aux;
a = 0.41015625 * (1 - (b*b) ) * (5-30*b*b+33*b*b*b*b);
aux = fabs(b);
if( aux < 1 ) { return(a); }
else { return(0); }
}

double kernel8( double b ) {
double a, aux;
a = 0.076904296875 * (1 - (b*b) ) * (35-385*b*b+1001*b*b*b*b-715*b*b*b*b*b*b);
aux = fabs(b);
if( aux < 1 ) { return(a); }
else { return(0); }
}

double kernel10( double b ) {
double a, aux;
a = 0.75 * (1 - (b*b) ) * (2.4609375 - 3.28125*b*b + 0.984375*b*b*b*b-0.09374*b*b*b*b*b*b-b*b*b*b*b*b*b*b/384);
aux = fabs(b);
if( aux < 1 ) { return(a); }
else { return(0); }
}




double suma_vec(double *x, int n){
register int i;
double sal=0;
for(i=0;i<n;i++){
	sal = sal + x[i];
}
return(sal);
}

double prod_vec(double *x, int n){
register int i;
double sal=1;
for(i=0;i<n;i++){
	sal = sal*(x[i]);
}
return(sal);
}

double norma_dos(double *x, int n){
register int i;
double sal=0;
for(i=0;i<n;i++){
	sal = sal + (x[i])*(x[i]);
}
return(sqrt(sal));
}

double potencia(double x, int n){
register int i;
double sal=1;
for(i=0;i<n;i++){
	sal = sal*x;
}
return(sal);
}

double l2dist(double *x, double *y, int n){
register int i;
double sal=0;
for(i=0;i<n;i++){
	sal = sal + (x[i] - y[i]) * (x[i] - y[i]);
}
return(sqrt(sal));
}


void ini_mu_pos_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *ventanas, double *prob, double *salida) {

double kepan(double);
double prod_vec(double*, int);
void reset_vec(double*, int);
double median(double*, int);

register int i, j;
int n = *nrow, d= *ncol, dimredu=0;
double *pesos, *aux11, *yredu, **xx;

pesos = (double *) malloc( n * sizeof(double));
aux11 = (double *) malloc( d * sizeof(double));
yredu = (double *) malloc( n * sizeof(double));

xx = (double **) malloc ( n * sizeof(double *) );
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );

for(i=0;i<n;i++){
	for(j=0;j<d;j++) xx[i][j]= *(x + j*n + i);
};


for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
 	    aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
	}
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
}

reset_vec(yredu, n);

for(i=0;i<n;i++){
	if( pesos[i]!=0){
		yredu[dimredu] = y[i];
		dimredu = dimredu + 1;
	}
}

*salida = median(yredu,dimredu);

for(i=0; i<n; i++) free(xx[i]); free(xx);
free(pesos); free(aux11);
free(yredu);

return;

}


void kernel_cl_pos_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *ventanas, double *eps, double *prob, double *salida) {

double kepan(double);
double suma_vec(double*, int);
double prod_vec(double*, int);
register int i, j;
int n = *nrow, d= *ncol;
double muhat;
double *pesos, *aux1, *aux11, **xx;

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));

aux11 = (double *) malloc( d * sizeof(double));

xx = (double **) malloc ( n * sizeof(double *) );
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );


for(i=0;i<n;i++){
	for(j=0;j<d;j++) xx[i][j]= *(x + j*n + i);
};

for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
 	    aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
	}
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
	aux1[i] = pesos[i] * y[i];
};

	muhat = suma_vec(aux1,n) / suma_vec(pesos,n);

*salida = muhat;

for(i=0; i<n; i++) free(xx[i]); free(xx);
free(pesos); free(aux1); free(aux11);

return;
}



void kernel_huber_pos_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *muhat_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida) {

double kepan(double);
double psi_huber_w(double, double);
double suma_vec(double*, int);
double prod_vec(double*, int);
register int i, j;
int n = *nrow, d= *ncol, it=0;
double corte, muhat, muold, khh;
double *pesos, *aux1, *aux2, *res, *aux11, **xx;

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));

aux11 = (double *) malloc( d * sizeof(double));

xx = (double **) malloc ( n * sizeof(double *) );
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );


for(i=0;i<n;i++){
	for(j=0;j<d;j++) xx[i][j]= *(x + j*n + i);
};

for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
 	    aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
	}
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
}

corte = 10 * (*eps);
khh = *kh;
muhat = *muhat_initial;

while( (corte > (*eps)) && (it< (*maxit)) ){
	for(i=0;i<n;i++){
		res[i] = ( y[i] - muhat) / *sigmahat;
		aux1[i] = pesos[i] * psi_huber_w(res[i],khh) * y[i];
		aux2[i] = pesos[i] * psi_huber_w(res[i],khh);
	};
	muold = muhat;
	muhat = suma_vec(aux1,n) / suma_vec(aux2,n);
	corte = fabs(muold-muhat) / ( fabs(muold) + *eps );
	it = it + 1;
};
*salida = muhat;

for(i=0; i<n; i++) free(xx[i]); free(xx);
free(pesos); free(aux1); free(aux2); free(res); free(aux11);

return;
}




void kernel_tukey_pos_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *muhat_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida) {

double kepan(double);
double psi_huber_w(double, double);
double suma_vec(double*, int);
double prod_vec(double*, int);
register int i, j;
int n = *nrow, d= *ncol, it=0;
double corte, muhat, muold, ktt;
double *pesos, *aux1, *aux2, *res, *aux11, **xx;

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));

aux11 = (double *) malloc( d * sizeof(double));

xx = (double **) malloc ( n * sizeof(double *) );
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );


for(i=0;i<n;i++){
	for(j=0;j<d;j++) xx[i][j]= *(x + j*n + i);
};

for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
 	    aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
	}
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
}

corte = 10 * (*eps);
ktt = *kt;
muhat = *muhat_initial;

while( (corte > (*eps)) && (it< (*maxit)) ){
	for(i=0;i<n;i++){
		res[i] = ( y[i] - muhat) / *sigmahat;
		aux1[i] = pesos[i] * psi_tukey_w(res[i],ktt) * y[i];
		aux2[i] = pesos[i] * psi_tukey_w(res[i],ktt);
	};
	muold = muhat;
	muhat = suma_vec(aux1,n) / suma_vec(aux2,n);
	corte = fabs(muold-muhat) / ( fabs(muold) + *eps );
	it = it + 1;
};
*salida = muhat;

for(i=0; i<n; i++) free(xx[i]); free(xx);
free(pesos); free(aux1); free(aux2); free(res); free(aux11);

return;
}



void kernel_cl_lin_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *ventanas, double *eps, double *prob, double *salida) {


double kepan(double);
double suma_vec(double*, int);
double prod_vec(double*, int);
void reset_mat(double**, int, int);
void reset_vec(double*, int);
void sum_mat(double **a, double **b, double **c, int n, int m);
void vec_vecprime(double **a, double *v1, double *v2, int n);
void scalar_vec(double *a, double b, double *c, int n);
void sum_vec(double *a, double *b, double *c, int n);
double vecprime_vec(double *a, double *b, int n);
void mat_vec(double **a, double *b, double *c, int n, int m);
double norma_dos(double *x, int n);
double l2dist(double *x, double *y, int n);
int lu(double **a,int *P, double *x);
register int i, j;
int n = *nrow, d = *ncol, q= *ncol+1;
double *pesos, *aux1, *aux11, **zz, **xx, *beta;
double **tmp1, **tmp2, *tmp3, *tmp4;

zz = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) zz[i] = (double *) malloc ( (d+1) * sizeof(double) );

xx = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );


pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
beta = (double *) malloc( (d+1) * sizeof(double));
aux11 = (double *) malloc( d * sizeof(double));

tmp1 = (double **) malloc( (d+1) * sizeof(double *));
tmp2 = (double **) malloc( (d+1) * sizeof(double *));
tmp3 = (double *) malloc( (d+1) * sizeof(double));
tmp4 = (double *) malloc( (d+1) * sizeof(double));

for(j=0;j<(d+1);j++) {
	tmp1[j] = (double *) malloc( (d+1) * sizeof(double) );
	tmp2[j] = (double *) malloc( (d+2) * sizeof(double) );
}


for(i=0;i<n; i++) {
     for(j=0; j<d; j++) xx[i][j] = *(x + j*n + i);
};

for(i=0;i<n; i++) {
     zz[i][0] = 1;
     for(j=0; j<d; j++) zz[i][j+1] = xx[i][j]-punto[j];
};


for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
 	    aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
	}
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
}

reset_mat(tmp1, d+1, d+1);
reset_vec(tmp3, d+1);
reset_vec(beta, d+1);

reset_mat(tmp2, d+1, d+2);
reset_vec(tmp4, d+1);

for(i=0;i<n;i++){
	aux1[i] = pesos[i] * y[i];
	scalar_vec(zz[i], pesos[i], tmp3, d+1);
	vec_vecprime(tmp1, tmp3, zz[i], d+1);
	sum_mat(tmp1, tmp2, tmp2, d+1, d+1);
	scalar_vec(zz[i], aux1[i], tmp3, d+1);
	sum_vec(tmp3, tmp4, tmp4, d+1);
};

for(j=0;j<(d+1);j++) tmp2[j][d+1] = tmp4[j];

if( lu(tmp2, &q, beta) == 1) {
	for(j=0; j<(d+1); j++) beta[j] = NA_REAL;
};

for(j=0; j<(d+1); j++) salida[j] = beta[j];

free(pesos); free(aux1); free(aux11);
for(i=0; i<n; i++){
	free(zz[i]);
	free(xx[i]);
};
free(zz);
free(xx);
for(j=0;j<(d+1);j++) {
	free( tmp1[j] ); free( tmp2[j] );
}
free(tmp1); free(tmp2); free(tmp3); free(tmp4);
free(beta);
return;
}




void kernel_huber_lin_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *beta_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida) {


double kepan(double);
double psi_huber_w(double, double);
double suma_vec(double*, int);
double prod_vec(double*, int);
void reset_mat(double**, int, int);
void reset_vec(double*, int);
void sum_mat(double **a, double **b, double **c, int n, int m);
void vec_vecprime(double **a, double *v1, double *v2, int n);
void scalar_vec(double *a, double b, double *c, int n);
void sum_vec(double *a, double *b, double *c, int n);
double vecprime_vec(double *a, double *b, int n);
void mat_vec(double **a, double *b, double *c, int n, int m);
double norma_dos(double *x, int n);
double l2dist(double *x, double *y, int n);
int lu(double **a,int *P, double *x);
register int i, j;
int n = *nrow, it=0, d = *ncol, q= *ncol+1;
double corte, khh; // muhat, muold, khh;
double *pesos, *aux1, *aux11, *aux2, *res, **zz, **xx, *beta, *beta_old;
double **tmp1, **tmp2, *tmp3, *tmp4;

zz = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) zz[i] = (double *) malloc ( (d+1) * sizeof(double) );

xx = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );


pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
beta = (double *) malloc( (d+1) * sizeof(double));
beta_old = (double *) malloc( (d+1) * sizeof(double));
aux11 = (double *) malloc( d * sizeof(double));

tmp1 = (double **) malloc( (d+1) * sizeof(double *));
tmp2 = (double **) malloc( (d+1) * sizeof(double *));
tmp3 = (double *) malloc( (d+1) * sizeof(double));
tmp4 = (double *) malloc( (d+1) * sizeof(double));

for(j=0;j<(d+1);j++) {
	tmp1[j] = (double *) malloc( (d+1) * sizeof(double) );
	tmp2[j] = (double *) malloc( (d+2) * sizeof(double) );
}

corte = 10 * (*eps);
khh = *kh;

for(i=0;i<n; i++) {
     for(j=0; j<d; j++) xx[i][j] = *(x + j*n + i);
};

for(i=0;i<n; i++) {
     zz[i][0] = 1;
     for(j=0; j<d; j++) zz[i][j+1] = xx[i][j]-punto[j];
};


for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
 	    aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
	}
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
}

reset_mat(tmp1, d+1, d+1);
reset_vec(tmp3, d+1);
reset_vec(beta, d+1);
reset_vec(beta_old, d+1);
for(j=0; j<(d+1) ; j++) beta[j] = beta_initial[j];
while( (corte > (*eps)) && (it< (*maxit)) ){
	reset_mat(tmp2, d+1, d+2);
	reset_vec(tmp4, d+1);
	for(j=0; j<(d+1); j++) beta_old[j] = beta[j];
	for(i=0;i<n;i++){
		res[i] = ( y[i] - vecprime_vec(zz[i], beta, d+1) ) / *sigmahat;
		aux1[i] = pesos[i] * psi_huber_w(res[i],khh) * y[i];
		aux2[i] = pesos[i] * psi_huber_w(res[i],khh);
		scalar_vec(zz[i], aux2[i], tmp3, d+1);
		vec_vecprime(tmp1, tmp3, zz[i], d+1);
		sum_mat(tmp1, tmp2, tmp2, d+1, d+1);
		scalar_vec(zz[i], aux1[i], tmp3, d+1);
		sum_vec(tmp3, tmp4, tmp4, d+1);
	};
	for(j=0;j<(d+1);j++) tmp2[j][d+1] = tmp4[j];
	if( lu(tmp2, &q, beta) == 1) {
		it = *maxit; // System became singular
		for(j=0; j<(d+1); j++) beta[j] = NA_REAL;
	};
	corte = l2dist(beta, beta_old, d+1) / ( norma_dos(beta_old, d+1) + *eps);
	it = it + 1;
};
for(j=0; j<(d+1); j++) salida[j] = beta[j];
free(pesos); free(aux1); free(aux2); free(res); free(aux11);
for(i=0; i<n; i++){
	free(zz[i]);
	free(xx[i]);
};
free(zz);
free(xx);
for(j=0;j<(d+1);j++) {
	free( tmp1[j] ); free( tmp2[j] );
}
free(tmp1); free(tmp2); free(tmp3); free(tmp4);
free(beta); free(beta_old);
return;
}


void kernel_tukey_lin_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *beta_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida) {


double kepan(double);
double psi_tukey_w(double, double);
double suma_vec(double*, int);
double prod_vec(double*, int);
void reset_mat(double**, int, int);
void reset_vec(double*, int);
void sum_mat(double **a, double **b, double **c, int n, int m);
void vec_vecprime(double **a, double *v1, double *v2, int n);
void scalar_vec(double *a, double b, double *c, int n);
void sum_vec(double *a, double *b, double *c, int n);
double vecprime_vec(double *a, double *b, int n);
void mat_vec(double **a, double *b, double *c, int n, int m);
double norma_dos(double *x, int n);
double l2dist(double *x, double *y, int n);
int lu(double **a,int *P, double *x);
register int i, j;
int n = *nrow, it=0, d = *ncol, q= *ncol+1;
double corte, ktt;
double *pesos, *aux1, *aux11, *aux2, *res, **zz, **xx, *beta, *beta_old;
double **tmp1, **tmp2, *tmp3, *tmp4;

zz = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) zz[i] = (double *) malloc ( (d+1) * sizeof(double) );

xx = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );


pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
beta = (double *) malloc( (d+1) * sizeof(double));
beta_old = (double *) malloc( (d+1) * sizeof(double));
aux11 = (double *) malloc( d * sizeof(double));

tmp1 = (double **) malloc( (d+1) * sizeof(double *));
tmp2 = (double **) malloc( (d+1) * sizeof(double *));
tmp3 = (double *) malloc( (d+1) * sizeof(double));
tmp4 = (double *) malloc( (d+1) * sizeof(double));

for(j=0;j<(d+1);j++) {
	tmp1[j] = (double *) malloc( (d+1) * sizeof(double) );
	tmp2[j] = (double *) malloc( (d+2) * sizeof(double) );
}

corte = 10 * (*eps);
ktt = *kt;

for(i=0;i<n; i++) {
     for(j=0; j<d; j++) xx[i][j] = *(x + j*n + i);
};

for(i=0;i<n; i++) {
     zz[i][0] = 1;
     for(j=0; j<d; j++) zz[i][j+1] = xx[i][j]-punto[j];
};


for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
 	    aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
	}
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
}

reset_mat(tmp1, d+1, d+1);
reset_vec(tmp3, d+1);
reset_vec(beta, d+1);
reset_vec(beta_old, d+1);
for(j=0; j<(d+1) ; j++) beta[j] = beta_initial[j];
while( (corte > (*eps)) && (it< (*maxit)) ){
	reset_mat(tmp2, d+1, d+2);
	reset_vec(tmp4, d+1);
	for(j=0; j<(d+1); j++) beta_old[j] = beta[j];
	for(i=0;i<n;i++){
		res[i] = ( y[i] - vecprime_vec(zz[i], beta, d+1) ) / *sigmahat;
		aux1[i] = pesos[i] * psi_tukey_w(res[i],ktt) * y[i];
		aux2[i] = pesos[i] * psi_tukey_w(res[i],ktt);
		scalar_vec(zz[i], aux2[i], tmp3, d+1);
		vec_vecprime(tmp1, tmp3, zz[i], d+1);
		sum_mat(tmp1, tmp2, tmp2, d+1, d+1);
		scalar_vec(zz[i], aux1[i], tmp3, d+1);
		sum_vec(tmp3, tmp4, tmp4, d+1);
	};
	for(j=0;j<(d+1);j++) tmp2[j][d+1] = tmp4[j];
	if( lu(tmp2, &q, beta) == 1) {
		it = *maxit; // System became singular
		for(j=0; j<(d+1); j++) beta[j] = NA_REAL;
	};
	corte = l2dist(beta, beta_old, d+1) / ( norma_dos(beta_old, d+1) + *eps);
	it = it + 1;
};
for(j=0; j<(d+1); j++) salida[j] = beta[j];
free(pesos); free(aux1); free(aux2); free(res); free(aux11);
for(i=0; i<n; i++){
	free(zz[i]);
	free(xx[i]);
};
free(zz);
free(xx);
for(j=0;j<(d+1);j++) {
	free( tmp1[j] ); free( tmp2[j] );
}
free(tmp1); free(tmp2); free(tmp3); free(tmp4);
free(beta); free(beta_old);
return;
}


void kernel_cl_alpha_multi(double *punto, double *x, int *alpha, int *degree, int *ncol, int *nrow, int *ordenkernel , double *y, double *ventanas, double *eps, double *prob, double *salida) {

double potencia(double x, int n);
double kepan(double);
double kernel4(double);
double kernel6(double);
double kernel8(double);
double kernel10(double);
double suma_vec(double*, int);
double prod_vec(double*, int);
void reset_mat(double**, int, int);
void reset_vec(double*, int);
void sum_mat(double **a, double **b, double **c, int n, int m);
void vec_vecprime(double **a, double *v1, double *v2, int n);
void scalar_vec(double *a, double b, double *c, int n);
void sum_vec(double *a, double *b, double *c, int n);
double vecprime_vec(double *a, double *b, int n);
void mat_vec(double **a, double *b, double *c, int n, int m);
double norma_dos(double *x, int n);
double l2dist(double *x, double *y, int n);
int lu(double **a,int *P, double *x);
register int i, j;
int n = *nrow, alfa= *alpha, q= *degree, Q= *degree+1, d = *ncol, ordenk = *ordenkernel;
double *pesos, *aux1, *aux11, **zz, **xx, *beta;
double **tmp1, **tmp2, *tmp3, *tmp4;

zz = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) zz[i] = (double *) malloc ( (q+1) * sizeof(double) );

xx = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
beta = (double *) malloc( (q+1) * sizeof(double));
aux11 = (double *) malloc( d * sizeof(double));

tmp1 = (double **) malloc( (q+1) * sizeof(double *));
tmp2 = (double **) malloc( (q+1) * sizeof(double *));
tmp3 = (double *) malloc( (q+1) * sizeof(double));
tmp4 = (double *) malloc( (q+1) * sizeof(double));

for(j=0;j<(q+1);j++){
	tmp1[j] = (double *) malloc( (q+1) * sizeof(double) );
	tmp2[j] = (double *) malloc( (q+2) * sizeof(double) );
};


for(i=0;i<n; i++) {
     for(j=0; j<d; j++) xx[i][j] = *(x + j*n + i);
};

for(i=0;i<n; i++) {
     zz[i][0] = 1;
     for(j=0; j<q; j++) zz[i][j+1] = zz[i][j]*(xx[i][alfa-1]-punto[alfa-1]);
};


for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
		if(j==(alfa-1)){
			aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
		}else{
			if(ordenk==2){
				aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
			};
			if(ordenk==4){
				aux11[j] = kernel4((xx[i][j] - punto[j])/(ventanas[j]));
			};
			if(ordenk==6){
				aux11[j] = kernel6((xx[i][j] - punto[j])/(ventanas[j]));
			};
			if(ordenk==8){
				aux11[j] = kernel8((xx[i][j] - punto[j])/(ventanas[j]));
			};
			if(ordenk==10){
				aux11[j] = kernel10((xx[i][j] - punto[j])/(ventanas[j]));
			};
		};
	};
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
};

reset_mat(tmp1, q+1, q+1);
reset_vec(tmp3, q+1);
reset_vec(beta, q+1);

reset_mat(tmp2, q+1, q+2);
reset_vec(tmp4, q+1);

for(i=0;i<n;i++){
	aux1[i] = pesos[i] * y[i];
	scalar_vec(zz[i], pesos[i], tmp3, q+1);
	vec_vecprime(tmp1, tmp3, zz[i], q+1);
	sum_mat(tmp1, tmp2, tmp2, q+1, q+1);
	scalar_vec(zz[i], aux1[i], tmp3, q+1);
	sum_vec(tmp3, tmp4, tmp4, q+1);
};
for(j=0;j<(q+1);j++) tmp2[j][q+1] = tmp4[j];

if( lu(tmp2, &Q, beta) == 1) {
	for(j=0; j<(q+1); j++) beta[j] = NA_REAL;
};

for(j=0; j<(q+1); j++) salida[j] = beta[j];

free(pesos); free(aux1); free(aux11);
for(i=0; i<n; i++){
	free(zz[i]);
	free(xx[i]);
};
free(zz);
free(xx);
for(j=0;j<(q+1);j++) {
	free( tmp1[j] ); free( tmp2[j] );
}
free(tmp1); free(tmp2); free(tmp3); free(tmp4);
free(beta);
return;
}


void kernel_huber_alpha_multi(double *punto, double *x, int *alpha, int *degree, int *ncol, int *nrow, int *ordenkernel , double *y, double *beta_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida) {

double potencia(double x, int n);
double kepan(double);
double kernel4(double);
double kernel6(double);
double kernel8(double);
double kernel10(double);
double psi_huber_w(double, double);
double suma_vec(double*, int);
double prod_vec(double*, int);
void reset_mat(double**, int, int);
void reset_vec(double*, int);
void sum_mat(double **a, double **b, double **c, int n, int m);
void vec_vecprime(double **a, double *v1, double *v2, int n);
void scalar_vec(double *a, double b, double *c, int n);
void sum_vec(double *a, double *b, double *c, int n);
double vecprime_vec(double *a, double *b, int n);
void mat_vec(double **a, double *b, double *c, int n, int m);
double norma_dos(double *x, int n);
double l2dist(double *x, double *y, int n);
int lu(double **a,int *P, double *x);
register int i, j;
int n = *nrow, it=0, alfa= *alpha, q= *degree, Q= *degree+1, d = *ncol, ordenk = *ordenkernel;
double corte, khh;
double *pesos, *aux1, *aux11, *aux2, *res, **zz, **xx, *beta, *beta_old;
double **tmp1, **tmp2, *tmp3, *tmp4;

zz = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) zz[i] = (double *) malloc ( (q+1) * sizeof(double) );

xx = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
beta = (double *) malloc( (q+1) * sizeof(double));
beta_old = (double *) malloc( (q+1) * sizeof(double));
aux11 = (double *) malloc( d * sizeof(double));

tmp1 = (double **) malloc( (q+1) * sizeof(double *));
tmp2 = (double **) malloc( (q+1) * sizeof(double *));
tmp3 = (double *) malloc( (q+1) * sizeof(double));
tmp4 = (double *) malloc( (q+1) * sizeof(double));

for(j=0;j<(q+1);j++){
	tmp1[j] = (double *) malloc( (q+1) * sizeof(double) );
	tmp2[j] = (double *) malloc( (q+2) * sizeof(double) );
};

corte = 10 * (*eps);
khh = *kh;

for(i=0;i<n; i++) {
     for(j=0; j<d; j++) xx[i][j] = *(x + j*n + i);
};

for(i=0;i<n; i++) {
     zz[i][0] = 1;
     for(j=0; j<q; j++) zz[i][j+1] = zz[i][j]*(xx[i][alfa-1]-punto[alfa-1]);
};


for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
		if(j==(alfa-1)){
			aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
		}else{
			if(ordenk==2){
				aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
			};
			if(ordenk==4){
				aux11[j] = kernel4((xx[i][j] - punto[j])/(ventanas[j]));
			};
			if(ordenk==6){
				aux11[j] = kernel6((xx[i][j] - punto[j])/(ventanas[j]));
			};
			if(ordenk==8){
				aux11[j] = kernel8((xx[i][j] - punto[j])/(ventanas[j]));
			};
			if(ordenk==10){
				aux11[j] = kernel10((xx[i][j] - punto[j])/(ventanas[j]));
			};
		};
	};
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
};

reset_mat(tmp1, q+1, q+1);
reset_vec(tmp3, q+1);
reset_vec(beta, q+1);
reset_vec(beta_old, q+1);

for(j=0; j<(q+1) ; j++) beta[j] = beta_initial[j];
while( (corte > (*eps)) && (it< (*maxit)) ){
	reset_mat(tmp2, q+1, q+2);
	reset_vec(tmp4, q+1);
	for(j=0; j<(q+1); j++) beta_old[j] = beta[j];
	for(i=0;i<n;i++){
		res[i] = ( y[i] - vecprime_vec(zz[i], beta, q+1) ) / *sigmahat;
		aux1[i] = pesos[i] * psi_huber_w(res[i],khh) * y[i];
		aux2[i] = pesos[i] * psi_huber_w(res[i],khh);
		scalar_vec(zz[i], aux2[i], tmp3, q+1);
		vec_vecprime(tmp1, tmp3, zz[i], q+1);
		sum_mat(tmp1, tmp2, tmp2, q+1, q+1);
		scalar_vec(zz[i], aux1[i], tmp3, q+1);
		sum_vec(tmp3, tmp4, tmp4, q+1);
	};
	for(j=0;j<(q+1);j++) tmp2[j][q+1] = tmp4[j];
	if( lu(tmp2, &Q, beta) == 1) {
		it = *maxit; // System became singular
		for(j=0; j<(q+1); j++) beta[j] = NA_REAL;
	};
	corte = l2dist(beta, beta_old, q+1) / ( norma_dos(beta_old, q+1) + *eps);
	it = it + 1;
};
for(j=0; j<(q+1); j++) salida[j] = beta[j];
free(pesos); free(aux1); free(aux2); free(res); free(aux11);
for(i=0; i<n; i++){
	free(zz[i]);
	free(xx[i]);
};
free(zz);
free(xx);
for(j=0;j<(q+1);j++) {
	free( tmp1[j] ); free( tmp2[j] );
}
free(tmp1); free(tmp2); free(tmp3); free(tmp4);
free(beta); free(beta_old);
return;
}



void kernel_tukey_alpha_multi(double *punto, double *x, int *alpha, int *degree, int *ncol, int *nrow, int *ordenkernel , double *y, double *beta_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida) {

double potencia(double x, int n);
double kepan(double);
double kernel4(double);
double kernel6(double);
double kernel8(double);
double kernel10(double);
double psi_huber_w(double, double);
double suma_vec(double*, int);
double prod_vec(double*, int);
void reset_mat(double**, int, int);
void reset_vec(double*, int);
void sum_mat(double **a, double **b, double **c, int n, int m);
void vec_vecprime(double **a, double *v1, double *v2, int n);
void scalar_vec(double *a, double b, double *c, int n);
void sum_vec(double *a, double *b, double *c, int n);
double vecprime_vec(double *a, double *b, int n);
void mat_vec(double **a, double *b, double *c, int n, int m);
double norma_dos(double *x, int n);
double l2dist(double *x, double *y, int n);
int lu(double **a,int *P, double *x);
register int i, j;
int n = *nrow, it=0, alfa= *alpha, q= *degree, Q= *degree+1, d = *ncol, ordenk = *ordenkernel;
double corte, ktt;
double *pesos, *aux1, *aux11, *aux2, *res, **zz, **xx, *beta, *beta_old;
double **tmp1, **tmp2, *tmp3, *tmp4;

zz = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) zz[i] = (double *) malloc ( (q+1) * sizeof(double) );

xx = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) xx[i] = (double *) malloc ( d * sizeof(double) );

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
beta = (double *) malloc( (q+1) * sizeof(double));
beta_old = (double *) malloc( (q+1) * sizeof(double));
aux11 = (double *) malloc( d * sizeof(double));

tmp1 = (double **) malloc( (q+1) * sizeof(double *));
tmp2 = (double **) malloc( (q+1) * sizeof(double *));
tmp3 = (double *) malloc( (q+1) * sizeof(double));
tmp4 = (double *) malloc( (q+1) * sizeof(double));

for(j=0;j<(q+1);j++){
	tmp1[j] = (double *) malloc( (q+1) * sizeof(double) );
	tmp2[j] = (double *) malloc( (q+2) * sizeof(double) );
};

corte = 10 * (*eps);
ktt = *kt;

for(i=0;i<n; i++) {
     for(j=0; j<d; j++) xx[i][j] = *(x + j*n + i);
};

for(i=0;i<n; i++) {
     zz[i][0] = 1;
     for(j=0; j<q; j++) zz[i][j+1] = zz[i][j]*(xx[i][alfa-1]-punto[alfa-1]);
};


for(i=0; i<n; i++) {
	for(j=0; j<d; j++){
		if(j!=(alfa-1)){
			if(ordenk==2){
				aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
			}
			if(ordenk==4){
				aux11[j] = kernel4((xx[i][j] - punto[j])/(ventanas[j]));
			}
			if(ordenk==6){
				aux11[j] = kernel6((xx[i][j] - punto[j])/(ventanas[j]));
			}
			if(ordenk==8){
				aux11[j] = kernel8((xx[i][j] - punto[j])/(ventanas[j]));
			}
			if(ordenk==10){
				aux11[j] = kernel10((xx[i][j] - punto[j])/(ventanas[j]));
			}
		}else{
 	    		aux11[j] = kepan((xx[i][j] - punto[j])/(ventanas[j]));
		}
	}
	pesos[i] = prod_vec( aux11, d) / (prob[i]);
}

reset_mat(tmp1, q+1, q+1);
reset_vec(tmp3, q+1);
reset_vec(beta, q+1);
reset_vec(beta_old, q+1);

for(j=0; j<(q+1) ; j++) beta[j] = beta_initial[j];
while( (corte > (*eps)) && (it< (*maxit)) ){
	reset_mat(tmp2, q+1, q+2);
	reset_vec(tmp4, q+1);
	for(j=0; j<(q+1); j++) beta_old[j] = beta[j];
	for(i=0;i<n;i++){
		res[i] = ( y[i] - vecprime_vec(zz[i], beta, q+1) ) / *sigmahat;
		aux1[i] = pesos[i] * psi_tukey_w(res[i],ktt) * y[i];
		aux2[i] = pesos[i] * psi_tukey_w(res[i],ktt);
		scalar_vec(zz[i], aux2[i], tmp3, q+1);
		vec_vecprime(tmp1, tmp3, zz[i], q+1);
		sum_mat(tmp1, tmp2, tmp2, q+1, q+1);
		scalar_vec(zz[i], aux1[i], tmp3, q+1);
		sum_vec(tmp3, tmp4, tmp4, q+1);
	};
	for(j=0;j<(q+1);j++) tmp2[j][q+1] = tmp4[j];
	if( lu(tmp2, &Q, beta) == 1) {
		it = *maxit; // System became singular
		for(j=0; j<(q+1); j++) beta[j] = NA_REAL;
	};
	corte = l2dist(beta, beta_old, q+1) / ( norma_dos(beta_old, q+1) + *eps);
	it = it + 1;
};
for(j=0; j<(q+1); j++) salida[j] = beta[j];
free(pesos); free(aux1); free(aux2); free(res); free(aux11);
for(i=0; i<n; i++){
	free(zz[i]);
	free(xx[i]);
};
free(zz);
free(xx);
for(j=0;j<(q+1);j++) {
	free( tmp1[j] ); free( tmp2[j] );
}
free(tmp1); free(tmp2); free(tmp3); free(tmp4);
free(beta); free(beta_old);
return;
}


void lu_R(double *a, double *b, int *kk, double *y) {
int i,j, k = *kk;
double **m1;
int lu(double **a,int *P, double *x);
m1 = (double **) malloc( k * sizeof(double *) );
for(j=0;j<k;j++)
	m1[j] = (double *) malloc( (k+1) * sizeof(double) );
for(i=0;i<k;i++) {
	for(j=0;j<k;j++)
		m1[i][j] = *(a + j*k + i);
	m1[i][k] = b[i];
};
lu(m1, kk, y);
for(j=0;j<k;j++) free(m1[j]);
free(m1);
return;
}





int lu(double **a,int *P, double *x)
{
int *pp,p;
register int i,j,k;
double *kk,s;
p = *P;
if ((pp = (int *) malloc(p*sizeof(int)))==NULL)
	{ Rprintf("\nNot enough memory in LU\n");
	  return(-1); }
/* pp vector storing the permutations */
for(j=0;j<p;j++)   /* cols */
{ pp[j]=j;
  for(i=j;i<p;i++)   /* filas */
	if ( fabs( a[i][j] ) > fabs( a[pp[j]][j] ) )
		pp[j]=i;
  if ( pp[j] != j )       /* permuto las filas cambiando los punt */
	{ kk=a[j];
	  a[j]=a[pp[j]];
	  a[pp[j]]=kk;
	};
  /* salida si el sistema resulta singular (det=0)
   * se detecta si el pivote (j,j) es cero  */
/*  if ( a[j][j] == 0 ) {   free(pp);
				return(1);
				}; */
    if ( fabs(a[j][j]) < TOL_INVERSE ) {   free(pp);
				return(1);
				};
  for(k=(j+1);k<p;k++)
	a[k][j] = a[k][j] / a[j][j];
  for(k=(j+1);k<p;k++)
	for(i=(j+1);i<p;i++)
		a[k][i] = a[k][i] - a[k][j] * a[j][i];

};    /* cierra el for de j */
for(i=0;i<p;i++)
	{ s=0.0;
	  for(j=0;j<i;j++)
	    s += a[i][j] * x[j];
	    x[i] = a[i][p] - s;          /* y[i]=a[i][p] */
	};
for(i=(p-1);i>=0;i--)
	{ s=0;
	  for(j=(i+1);j<p;j++)
	    s += a[i][j] * x[j];
	  x[i] = (x[i] - s) / a[i][i];
	  };
free(pp);
return(0);
}


void sum_mat(double **a, double **b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		c[i][j] = a[i][j] + b[i][j];
}

void vec_vecprime(double **a, double *v1, double *v2, int n)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<n;j++)
		a[i][j] = v1[i] * v2[j];/* could take advantage of symmetry */
}

void scalar_mat(double **a, double b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
        for(j=0;j<m;j++)
	c[i][j]  = b * a[i][j];
}

void scalar_vec(double *a, double b, double *c, int n)
{
register int i;
for(i=0;i<n;i++)
	c[i]  = b * a[i];
}

double vecprime_vec(double *a, double *b, int n)
{
register int i;
double s = 0.0;
for(i=0;i<n;i++) s += a[i] * b[i];
return(s);
}

void sum_vec(double *a, double *b, double *c, int n)
{
register int i;
for(i=0;i<n;i++) c[i] = a[i] + b[i];
}

void dif_vec(double *a, double *b, double *c, int n)
{
register int i;
for(i=0;i<n;i++) c[i] = a[i] - b[i];
}

void dif_mat(double **a, double **b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++) c[i][j] = a[i][j] - b[i][j];
}

void mat_vec(double **a, double *b, double *c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(c[i]=0,j=0;j<m;j++) c[i] += a[i][j] * b[j];
}

void mat_mat(double **a, double **b, double **c, int n,
		int m, int l)
{
register int i,j,k;
for(i=0;i<n;i++)
	for(j=0;j<l;j++) {
	c[i][j] = 0;
	for(k=0;k<m;k++) c[i][j] += a[i][k] * b[k][j];
	};
}

/* void disp_vec(double *a, int n)
{
register int i;
Rprintf("\n");
for(i=0;i<n; i++) Rprintf("%lf ",a[i]);
Rprintf("\n");
}

void disp_mat(double **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++) {
Rprintf("\n");
for(j=0;j<m;j++) Rprintf("%10.8f ",a[i][j]);
};
Rprintf("\n");
}
*/

int inverse(double **a, double **b, int n)
{
int lu(double **, int *, double *);
void mat_vec(double **, double *, double *, int, int);
register int i,j,k;
double **c, *e;
c = (double **) malloc( n * sizeof(double *));
e = (double *) malloc( n * sizeof(double));
for(i=0;i<n;i++) c[i] = (double *) malloc ( (n+1) * sizeof(double) );
for(i=0;i<n;i++) {   /* i-th column */

for(j=0;j<n;j++)
	for(k=0;k<n;k++) c[j][k] = a[j][k];
for(j=0;j<i;j++) c[j][n] = 0.0;
c[i][n] = 1.0;
for(j=i+1;j<n;j++) c[j][n] = 0.0;
if( lu(c,&n,e) == 1) {
	for(i=0;i<n;i++) free(c[i]);
	free(c);free(e);
	return(1);
	};
for(j=0;j<n;j++) b[j][i] = e[j] ;
};
for(i=0;i<n;i++) free(c[i]);
free(c);free(e);
return(0);
}

void reset_mat(double **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		a[i][j] = 0.0;
}

void reset_vec(double *a, int n)
{
register int i;
for(i=0;i<n;i++) a[i] = 0.0;
}


double median(double *x, int n)
{
double kthplace(double *,int,int);
double *aux,t;
register int i;
if ( (aux = (double *) malloc (n*sizeof(double)) )==NULL)
{Rprintf("\nNot enought memory in median\n"); return(-99.0); };
for(i=0;i<n;i++) aux[i]=x[i];
if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2 ;
else	t = kthplace(aux,n, n/2+1 ) ;
free(aux);
return(t);
}

double median_abs(double *x, int n)
{
double kthplace(double *,int,int);
double *aux,t;
register int i;
if ( (aux = (double *) malloc (n*sizeof(double)) )==NULL )
{ Rprintf("\nNot enought memory in med_abs\n");return(-99.0);};
for(i=0;i<n;i++) aux[i]=fabs(x[i]);
if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2 ;
else 	t = kthplace(aux,n, n/2+1 ) ;
free(aux);
return(t);
}


double kthplace(double *a, int n, int k)
{
int jnc,j;
int l,lr;
double ax,w;
k--;
l=0;
lr=n-1;
while (l<lr)
	{ ax=a[k];
	  jnc=l;
	  j=lr;
	  while (jnc<=j)
		{ while (a[jnc] < ax) jnc++;
		  while (a[j] > ax) j--;
		  if (jnc <= j)
			{ w=a[jnc];
			  a[jnc]=a[j];
			  a[j]=w;
			  jnc++;
			  j--;
			};
		};
	  if (j<k) l=jnc;
	if (k<jnc) lr=j;
	};
return(a[k]);
}


void huber_pos(int *nrow, double *y, double *muhat_initial,
 double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida) {

double psi_huber_w(double, double);
double suma_vec(double*, int);
register int i;
int n = *nrow, it=0;
double corte, muhat, muold, khh;
double *aux1, *aux2, *res;

aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
corte = 10 * (*eps);
khh = *kh;
muhat = *muhat_initial;

while( (corte > (*eps)) && (it< (*maxit)) ){
	for(i=0;i<n;i++){
		res[i] = ( y[i] - muhat) / *sigmahat;
		aux1[i] = psi_huber_w(res[i],khh) * y[i] / prob[i];
		aux2[i] = psi_huber_w(res[i],khh) / prob[i];
	};
	muold = muhat;
	muhat = suma_vec(aux1,n) / suma_vec(aux2,n);
	corte = fabs(muold-muhat) / (fabs(muold) + *eps);
	it = it + 1;
};
*salida = muhat;
free(aux1); free(aux2); free(res);
return;
}

void tukey_pos(int *nrow, double *y, double *muhat_initial,
 double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida) {

double psi_tukey_w(double, double);
double suma_vec(double*, int);
register int i;
int n = *nrow, it=0;
double corte, muhat, muold, ktt;
double *aux1, *aux2, *res;

aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
corte = 10 * (*eps);
ktt = *kt;
muhat = *muhat_initial;

while( (corte > (*eps)) && (it< (*maxit)) ){
	for(i=0;i<n;i++){
		res[i] = ( y[i] - muhat) / *sigmahat;
		aux1[i] = psi_tukey_w(res[i],ktt) * y[i] / prob[i];
		aux2[i] = psi_tukey_w(res[i],ktt) / prob[i];
	};
	muold = muhat;
	muhat = suma_vec(aux1,n) / suma_vec(aux2,n);
	corte = fabs(muold-muhat) / (fabs(muold) + *eps);
	it = it + 1;
};
*salida = muhat;
free(aux1); free(aux2); free(res);
return;
}

