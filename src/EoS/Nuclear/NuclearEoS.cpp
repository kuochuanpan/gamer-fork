#include "NuclearEoS.h"

#if ( MODEL == HYDRO  &&  EOS == EOS_NUCLEAR )



#ifdef __CUDACC__

#include "cubinterp_some.cu"
#include "findenergy.cu"

#else

void nuc_eos_C_cubinterp_some( const real x, const real y, const real z,
                               real *output_vars, const real *alltables,
                               const int nx, const int ny, const int nz, const int nvars,
                               const real *xt, const real *yt, const real *zt );
void findenergy( const real x, const real y, const real z,
                 real *found_leps, const real *alltables_mode,
                 const int nx, const int ny, const int nz, const int neps,
                 const real *xt, const real *yt, const real *zt,
                 const real *logeps, const int keymode, int *keyerr );

#endif // #ifdef __CUDACC__ ... else ...




//-----------------------------------------------------------------------------------------------
// Function    :  nuc_eos_C_short
// Description :  Function to find thermodynamic varibles by searching
//                a pre-calculated nuclear equation of state table
//
// Note        :  1. It will strictly return values in cgs or MeV
//                2. Four modes are supported
//                   --> Energy      mode (0)
//                       Temperature mode (1)
//                       Entropy     mode (2)
//                       Pressure    mode (3)
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
//                                 --> 0 : Energy mode      (input eps    )
//                                     1 : Temperature mode (input T      )
//                                     2 : Entropy mode     (input entropy)
//                                     3 : Pressure mode    (input P      )
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
void nuc_eos_C_short( const real xrho, real *xenr, const real xye,
                      real *xtemp, real *xent, real *xprs,
                      real *xcs2, real *xmunu, const real energy_shift,
                      const int nrho, const int neps, const int nye, const int nmode,
                      const real *alltables, const real *alltables_mode,
                      const real *logrho, const real *logeps, const real *yes,
                      const real *logtemp_mode, const real *entr_mode, const real *logprss_mode,
                      const int keymode, int *keyerr, const real rfeps )
{

   real lr   = LOG10( xrho );
   real leps = ( keymode == NUC_MODE_ENGY ) ? LOG10( *xenr + energy_shift ) : NULL_REAL;


// check whether (rho, eps, Y_e) is within the table
   *keyerr = 0;

   if ( lr > logrho[nrho-1] )
   {
      *keyerr = 105;
      return;
   }

   if ( lr < logrho[0] )
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

   if ( keymode == NUC_MODE_ENGY )
   {
      if ( leps > logeps[neps-1] )
      {
         *keyerr = 103;
         return;
      }

      if ( leps < logeps[0] )
      {
         *keyerr = 104;
         return;
      }
   }


// find energy
   switch ( keymode )
   {
      case NUC_MODE_ENGY :
      {
         leps = MAX( leps, (real)0.0 );
      }
      break;

      case NUC_MODE_TEMP :
      {
         const real lt = LOG10( *xtemp );

         findenergy( lr, lt, xye, &leps, alltables_mode, nrho, nmode, nye, neps,
                     logrho, logtemp_mode, yes, logeps, keymode, keyerr );

         if ( *keyerr != 0 )  return;
      }
      break;

      case NUC_MODE_ENTR :
      {
         const real entr = *xent;

         findenergy( lr, entr, xye, &leps, alltables_mode, nrho, nmode, nye, neps,
                     logrho, entr_mode, yes, logeps, keymode, keyerr );

         if ( *keyerr != 0 )  return;
      }
      break;

      case NUC_MODE_PRES :
      {
         const real lprs = LOG10( *xprs );

         findenergy( lr, lprs, xye, &leps, alltables_mode, nrho, nmode, nye, neps,
                     logrho, logprss_mode, yes, logeps, keymode, keyerr );

         if ( *keyerr != 0 )  return;
      }
      break;
   } // switch ( keymode )


   real res[5]; // result array

// linear interolation for other variables
// nuc_eos_C_linterp_some( lr, leps, xye, res, alltables,
//                         nrho, neps, nye, 5, logrho, logeps, yes );

// cubic interpolation for other variables
   nuc_eos_C_cubinterp_some( lr, leps, xye, res, alltables, nrho, neps, nye, 5,
                             logrho, logeps, yes );

   if ( keymode != NUC_MODE_ENGY )  *xenr = POW( (real)10.0, leps ) - energy_shift;


// assign results
   *xprs  = POW( (real)10.0, res[0] );
   *xtemp = POW( (real)10.0, res[1] );
   *xent  = res[2];
   *xmunu = res[3];
   *xcs2  = res[4];


   return;

} // FUNCTION : nuc_eos_C_short



#endif // #if ( MODEL == HYDRO  &&  EOS == EOS_NUCLEAR )
