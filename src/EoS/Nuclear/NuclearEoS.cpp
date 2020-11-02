#include "NuclearEoS.h"

#if ( MODEL == HYDRO  &&  EOS == EOS_NUCLEAR )



#ifdef __CUDACC__

#include "cubinterp_some.cu"
#include "findenergy.cu"

#else

void nuc_eos_C_cubinterp_some( double x, double y, double z,
                               double *output_vars, const double *alltables,
                               int nx, int ny, int nz, int nvars,
                               const double *xt, const double *yt, const double *zt );;
void findenergy( double x, double y, double z, double *found_leps, const double *alltables_mode,
                 int nx, int ny, int nz, int neps, const double *xt, const double *yt, const double *zt,
                 const double *logeps, int keymode, int *keyerr );

#endif // #ifdef __CUDACC__ ... else ...




//-----------------------------------------------------------------------------------------------
// Function    :  nuc_eos_C_short
// Description :  Function to find thermodynamic varibles by searching
//                a pre-calculated nuclear equation of state table
//
// Note        :  1. It will strictly return values in cgs or MeV
//                2. Four modes are supported
//                   --> energy      mode (0)
//                       temperature mode (1)
//                       entropy     mode (2)
//                       pressure    mode (3)
//
// Parameter   :  xrho           : Input density (rho (g/cm^3))
//                xenr           : Specific internal energy (eps)
//                xye            : Electron fraction (Y_e)
//                xtemp          : Input (temperature mode) or ouput temperature in MeV
//                xent           : Input (entropy mode) or output specific entropy (e)
//                xprs           : Input (pressure mode) or output pressure
//                xcs2           : Output sound speed
//                xmunu          : Output chemical potential
//                energy_shift   : Energy shift
//                nrho           : Size of density array in the Nuclear EoS table
//                neps           : Size of energy  array in the Nuclear EoS table
//                nye            : Size of Y_e     array in the Nuclear EoS table
//                nmode          : Size of array for each mode in the Nuclear EoS table
//                                 --> log(T)  (1)
//                                     entropy (2)
//                                     log(P)  (3)
//                alltables      : Nuclear EoS table
//                alltables_mode : Auxiliary log(eps) arrays for different modes
//                logrho         : log(rho) array in the table
//                logeps         : log(eps) array in the table
//                yes            : Y_e      array in the table
//                logtemp_mode   : log(T)   array for the temperature mode
//                entr_mode      : entropy  array for the entropy mode
//                logprss_mode   : log(P)   array for the pressure mode
//                keymode        : Which mode we will use
//                                 --> 0 : energy mode      (input eps    )
//                                     1 : temperature mode (input T      )
//                                     2 : entropy mode     (input entropy)
//                                     3 : pressure mode    (input P      )
//                keyerr         : Output error
//                                 --> 667 : fail in finding energy (T, e, P modes)
//                                     101 : Y_e too high
//                                     102 : Y_e too low
//                                     103 : eps too high (if keymode = 0)
//                                     104 : eps too low  (if keymode = 0)
//                                     105 : rho too high
//                                     106 : rho too low
//                rfeps          : Tolerence for interpolations
//                                 --> Not used currently
//-----------------------------------------------------------------------------------------------
GPU_DEVICE
void nuc_eos_C_short( double xrho, double *xenr, double xye,
                      double *xtemp, double *xent, double *xprs,
                      double *xcs2, double *xmunu, double energy_shift,
                      int nrho, int neps, int nye, int nmode,
                      const double *alltables, const double *alltables_mode,
                      const double *logrho, const double *logeps, const double *yes,
                      const double *logtemp_mode, const double *entr_mode, const double *logprss_mode,
                      int keymode, int *keyerr, double rfeps )
{

   *keyerr = 0;

// check whether (rho, eps, Y_e) is within the table
   if ( log10(xrho) > logrho[nrho-1] )
   {
      *keyerr = 105;
      return;
   }

   if ( log10(xrho) < logrho[0] )
   {
      *keyerr = 106;
      return;
   }

   if ( xye > yes[nye-1] )
   {
      *keyerr = 101;
      return;
   }

   if ( xye < yes[0] )
   {
      *keyerr = 102;
      return;
   }

   if ( keymode == 0 )
   {
      if ( log10(*xenr + energy_shift) > logeps[neps-1] )
      {
         *keyerr = 103;
         return;
      }
      if ( log10(*xenr + energy_shift) < logeps[0] )
      {
         *keyerr = 104;
         return;
      }
   }

// set up local vars
   double lr = log10(xrho);
   double xeps = *xenr + energy_shift;
   double leps = log10(MAX(xeps, 1.0));

// find energy if needed
// ( keymode = 0: energy mode      )
// ( keymode = 1: temperature mode )
// ( keymode = 2: entropy mode     )
// ( keymode = 3  pressure mode    )
   if ( keymode == 1 )
   {
      double lt = log10(*xtemp);
      findenergy( lr, lt, xye, &leps, alltables_mode, nrho, nmode, nye, neps,
                  logrho, logtemp_mode, yes, logeps, keymode, keyerr );
      if ( *keyerr != 0 ) return;
   }

   else if ( keymode == 2 )
   {
      double entr = *xent;
      findenergy( lr, entr, xye, &leps, alltables_mode, nrho, nmode, nye, neps,
                  logrho, entr_mode, yes, logeps, keymode, keyerr );
      if ( *keyerr != 0 ) return;
   }

   else if ( keymode == 3 )
   {
      double lprs = log10(*xprs);
      findenergy( lr, lprs, xye, &leps, alltables_mode, nrho, nmode, nye, neps,
                  logrho, logprss_mode, yes, logeps, keymode, keyerr );
      if ( *keyerr != 0 ) return;
   }

   double res[5]; // result array

// linear interolation for other variables
// nuc_eos_C_linterp_some( lr, leps, xye, res, alltables,
//                         nrho, neps, nye, 5, logrho, logeps, yes );

// cubic interpolation for other variables
   nuc_eos_C_cubinterp_some( lr, leps, xye, res, alltables,
                             nrho, neps, nye, 5, logrho, logeps, yes );

   if ( keymode != 0 )  *xenr = pow(10.0, leps) - energy_shift;


// assign results
   *xprs  = pow(10.0, res[0]);
   *xtemp = pow(10.0, res[1]);
   *xent  = res[2];
   *xmunu = res[3];
   *xcs2  = res[4];


   return;

} // FUNCTION : nuc_eos_C_short



#endif // #if ( MODEL == HYDRO  &&  EOS == EOS_NUCLEAR )
