#include "GAMER.h"
#include "NuclearEos.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Src_Deleptonization( real fluid[], const double x, const double y, const double z, const double Time,
                      const int lv, double AuxArray[], const double dt, const double EngyB );

void (*Src_Deleptonization_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                      const int lv, double AuxArray[], const double dt, const double EngyB ) = Src_Deleptonization;

double YeOfRhoFunc(double);

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Deleptonization
// Description :  User-defined source terms
//
// Note        :  1. Invoked by Src_AdvanceDt() using the function pointer "Src_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Enabled by the runtime option "SRC_USER"
//
// Parameter   :  fluid    : Fluid array storing both the input and updated values
//                           --> Array size = NCOMP_TOTAL (so it includes both active and passive variables)
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//                dt       : Time interval to advance solution
//
// Return      :  fluid[]
//-------------------------------------------------------------------------------------------------------
void Src_Deleptonization( real fluid[], const double x, const double y, const double z, const double Time,
               const int lv, double AuxArray[], const double dt, const double EngyB )
{

#  if ( EOS == NUCLEAR )

// example
   /*
   const double CoolRate = 1.23; // set arbitrarily here
   double Ek, Eint;              // kinetic and internal energies

   Ek    = (real)0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   Eint  = fluid[ENGY] - Ek;
   Eint -= CoolRate*dt;
   fluid[ENGY] = Ek + Eint;
   */

   const double delep_minDens = 1.e6; // [g/cm^3]
   const double mev_to_kelvin = 1.1604447522806e10;
   const double Q = 1.293333;

   double xdens,dens, ener, entr, ye, xmom, ymom, zmom;
   double del_ye, del_entr;
   double xtmp, xenr, xprs, xent, xcs2, xdedt, xdpderho, xdpdrhoe, xmunu;
   int keyerr;
   const double rfeps = 1.0e-10;

   double debug1, debug2, debug3;
   double yout;
   const double Gamma_m1   = GAMMA - 1.0;

   if (EOS_POSTBOUNCE)
   {
      return;
   }

   // Deleptonization

   dens  = fluid[DENS]; // code units
   xdens = dens*UNIT_D; // [g/cm^3]
   entr  = fluid[ENTR]/dens;
   xmom  = fluid[MOMX]; // code unit
   ymom  = fluid[MOMY];
   zmom  = fluid[MOMZ];
   ye    = fluid[YE]/dens;
   ener  = fluid[ENGY];
   ener  = ener - 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/dens; // internal energy
#  ifdef MHD
   ener -= EngyB;
#  endif
   xenr  = (ener/dens*UNIT_V*UNIT_V) - energy_shift; // specific internal energy [need nuclear EoS]


   del_ye   = 0.0;
   del_entr = 0.0;

   if (xdens <= delep_minDens)
   {
     del_ye = 0.0;
   } else
   {
     yout = YeOfRhoFunc(xdens);
     del_ye = yout - ye;
     del_ye = MIN(0.0, del_ye); // Deleptonization cannot increase Ye
   }

   if (del_ye < 0.0)
   {

      xtmp = 10.0; // trial value

      nuc_eos_C_short(xdens,&xtmp,ye,&xenr, &xprs, &xent, &xcs2, &xdedt, &xdpderho,
          &xdpdrhoe, &xmunu, 0, &keyerr, rfeps); // energy mode


      xmunu += Q;

      if ((xmunu < DELEP_ENU) || (xdens >= 2.e12))
      {
        del_entr = 0.0;
      } else
      {
        del_entr = - del_ye * (xmunu - DELEP_ENU) / xtmp;
      }

      fluid[ENTR] = dens*(entr + del_entr);
      fluid[YE]   = dens*(ye + del_ye);

      xent = entr + del_entr;
      ye   = ye + del_ye;

      nuc_eos_C_short(xdens,&xtmp,ye,&xenr, &xprs, &xent, &xcs2, &xdedt, &xdpderho,
          &xdpdrhoe, &xmunu, 2, &keyerr, rfeps); // entropy mode

      fluid[ENGY] = (dens/(UNIT_V*UNIT_V))*(xenr + energy_shift) + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
#  ifdef MHD
      fluid[ENGY] += EngyB;
#  endif

      // if using Dual energy
#     ifdef DUAL_ENERGY
#     if (DUAL_ENERGY == DE_ENPY)
      xprs = xprs * UNIT_P;
      fluid[ENPY] = Hydro_DensPres2Entropy( dens, xprs, Gamma_m1 );
#     endif
#     endif

   }

#  endif // # if ( EOS == NUCLEAR )

} // FUNCTION : Src_User

double YeOfRhoFunc(double xdens)
{
# ifdef DELEPTIONIZATION

  double xofrho, ye;

  xofrho = 2.0*log10(xdens) - log10(DELEP_RHO2) - log10(DELEP_RHO1);
  xofrho = xofrho / (log10(DELEP_RHO2) - log10(DELEP_RHO1));
  xofrho = MAX(-1.0, MIN(1.0, xofrho));
  ye = 0.5*(DELEP_YE2 + DELEP_YE1) + 0.5*xofrho*(DELEP_YE2 - DELEP_YE1);
  ye = ye + DELEP_YEC*(1.0 - fabs(xofrho));
  ye = ye + DELEP_YEC*4.0*fabs(xofrho)*(fabs(xofrho)-0.5)*(fabs(xofrho) - 1.0);
  return ye;

# endif
}
