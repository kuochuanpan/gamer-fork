//#ifndef __NUCLEAREOS_H__
//#define __NUCLEAREOS_H__

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NTABLES 19
//#define DEBUG 1

// 6. Global variables for nuclear Eos

extern int nrho;
extern int ntemp;
extern int nye;

extern double *alltables;
extern double *logrho;
extern double *logtemp;
extern double *yes;
extern double energy_shift;
extern double dtemp, dtempi;
extern double drho, drhoi;
extern double dye, dyei;

// min and max values

extern double eos_rhomax, eos_rhomin;
extern double eos_tempmin, eos_tempmax;
extern double eos_yemin, eos_yemax;

// some vectors for selecting variables for more
// efficient interpolation
extern int ivs_short[19];

// table key
// 0 logpress
// 1 logenergy
// 2 entropy
// 3 munu
// 4 cs2
// 5 dedt
// 6 dpdrhoe
// 7 dpderho
// 8 muhat
// 9 mu_e
// 10 mu_p
// 11 mu_n
// 12 Xa
// 13 Xh
// 14 Xn
// 15 Xp
// 16 Abar
// 17 Zbar
// 18 Gamma

#if ( defined GRAVITY  &&  defined GREP )
   extern int    GREP_Center_Method;
   extern int    GREP_MaxIter;
   extern bool   GREP_LogBin;
   extern bool   GREP_RemoveEmptyBin;
   extern double GREP_LogBinRatio;
   extern double GREP_MaxRadius;
   extern double GREP_MinBinSize;
#endif


// frontend function declarations
void nuc_eos_C_short(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs, double* xent,
		     double *xcs2, double* xdedt, double* xdpderho,
		     double *xdpdrhoe, double* xmunu, int keytemp,
		     int *keyerr,double rfeps);

// core function declarations

void nuc_eos_C_ReadTable(char* nuceos_table_name);

void nuc_eos_C_linterp_many(double x, double y, double z,
			    double* f, double* ft,
			    int nx, int ny, int nz, int nvars,
			    double* xt,double*yt, double* zt);

void nuc_eos_C_linterp_some(double x, double y, double z,
			    double* f, double* ft,
			    int* ivs,
			    int nx, int ny, int nz, int nvars,
			    double* xt,double*yt, double* zt);

void nuc_eos_C_linterp_for_temp(double x, double y, double z,
				double* f, double* ft,
				int nx, int ny, int nz,
				double* xt, double*yt, double* zt,
				double* linterp_for_temp);

void nuc_eos_C_linterp_for_entr(double x, double y, double z,
				double* f, double* ft,
				int nx, int ny, int nz,
				double* xt, double*yt, double* zt,
				double* linterp_for_entr);

void nuc_eos_C_findtemp(double lr, double lt0, double ye,
			double leps, double prec, double *lt,
			int *keyerr);

void nuc_eos_C_findtemp_entropy(double lr, double lt0, double ye,
			double xs, double *lt, double prec, int *keyerr);

void nuc_eos_C_testing();

//#endif // __NUCLEAREOS_H__
