/* Copied from robustbase.h */

#include <R.h>
#include <Rinternals.h>
#include <complex.h>

/* For internationalized messages */
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#define dngettext(pkg, String, StringP, N) (N > 1 ? StringP : String)
#endif


double psi_huber_w( double r, double k);
double psi_tukey_w( double r, double k);
double kepan( double b );
double kernel4( double b );
double kernel6( double b );
double kernel8( double b );
double kernel10( double b );
double suma_vec(double *x, int n);
double prod_vec(double *x, int n);
double norma_dos(double *x, int n);
double potencia(double x, int n);
double l2dist(double *x, double *y, int n);
void lu_R(double *a, double *b, int *kk, double *y);
int lu(double **a,int *P, double *x);
void sum_mat(double **a, double **b, double **c, int n, int m);
void vec_vecprime(double **a, double *v1, double *v2, int n);
void scalar_mat(double **a, double b, double **c, int n, int m);
void scalar_vec(double *a, double b, double *c, int n);
double vecprime_vec(double *a, double *b, int n);
void sum_vec(double *a, double *b, double *c, int n);
void dif_vec(double *a, double *b, double *c, int n);
void dif_mat(double **a, double **b, double **c, int n, int m);
void mat_vec(double **a, double *b, double *c, int n, int m);
void mat_mat(double **a, double **b, double **c, int n, int m, int l);
int inverse(double **a, double **b, int n);
void reset_mat(double **a, int n, int m);
void reset_vec(double *a, int n);
double median(double *x, int n);
double median_abs(double *x, int n);
double kthplace(double *a, int n, int k);

void ini_mu_pos_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *ventanas, double *prob, double *salida);

void kernel_cl_pos_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *ventanas, double *eps, double *prob, double *salida);

void kernel_huber_pos_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *muhat_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida);

void kernel_tukey_pos_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *muhat_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida);

void kernel_cl_lin_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *ventanas, double *eps, double *prob, double *salida);

void kernel_huber_lin_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *beta_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida);

void kernel_tukey_lin_multi(double *punto, double *x, int *ncol, int *nrow, double *y, double *beta_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida);

void kernel_cl_alpha_multi(double *punto, double *x, int *alpha, int *degree, int *ncol, int *nrow, int *ordenkernel , double *y, double *ventanas, double *eps, double *prob, double *salida);

void kernel_huber_alpha_multi(double *punto, double *x, int *alpha, int *degree, int *ncol, int *nrow, int *ordenkernel , double *y, double *beta_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida);

void kernel_tukey_alpha_multi(double *punto, double *x, int *alpha, int *degree, int *ncol, int *nrow, int *ordenkernel , double *y, double *beta_initial, double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida);

void huber_pos(int *nrow, double *y, double *muhat_initial,
 double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida);

void tukey_pos(int *nrow, double *y, double *muhat_initial,
 double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida);

void matriz_salida(double *xx, int *n, int *p, double *salida);

