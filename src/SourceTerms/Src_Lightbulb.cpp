#include "GAMER.h"
#include "NuclearEos.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Src_LightBulb( real fluid[], const double x, const double y, const double z, const double Time,
                      const int lv, double AuxArray[], const double dt );

// this function pointer may be overwritten by various test problem initializers
void (*Src_LightBulb_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                      const int lv, double AuxArray[], const double dt ) = Src_LightBulb;


//-------------------------------------------------------------------------------------------------------
// Function    :  Src_LightBulb
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
void Src_LightBulb( real fluid[], const double x, const double y, const double z, const double Time,
               const int lv, double AuxArray[], const double dt )
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

   const double mev_to_kelvin = 1.1604447522806e10;

   double xdens,dens, ener, entr, ye;
   double xtmp, xenr, xprs, xent, xcs2, xdedt, xdpderho, xdpdrhoe, xmunu;
   double radius, xc, yc, zc;
   int keyerr;
   const double rfeps = 1.0e-10;

   double xXp, xXn;
   double dEneut, T6;

   const double  BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

   double logd, logt;
   double res[19];

   const double Gamma_m1   = GAMMA - 1.0;

   if (!EOS_POSTBOUNCE)
   {
        return;
   }

   dEneut = 0.0 ;

   dens  = fluid[DENS]; // code units
   xdens = dens*UNIT_D; // [g/cm^3]
   entr  = fluid[ENTR]/dens;
   ye    = fluid[YE]/dens;
   ener  = fluid[ENGY];
   ener  = ener - 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/dens; // internal energy
   xenr  = (ener/dens*UNIT_V*UNIT_V) - energy_shift; // specific internal energy [need nuclear EoS]

   xtmp = 10.0; // trial value [MeV]
   nuc_eos_C_short(xdens,&xtmp,ye,&xenr, &xprs, &xent, &xcs2, &xdedt, &xdpderho,
                     &xdpdrhoe, &xmunu, 0, &keyerr, rfeps); // energy mode

   logd = MIN(MAX(xdens, eos_rhomin), eos_rhomax);
   logt = MIN(MAX(xtmp,  eos_tempmin), eos_tempmax);

   logd = log10(logd);
   logt = log10(logt);


   // find xp xn
   nuc_eos_C_linterp_some(logd, logt, ye, res, alltables,
           ivs_short, nrho, ntemp, nye, 19,
           logrho, logtemp, yes);

   xXn = res[14];
   xXp = res[15];

   //printf("debug: xXp %13.7e  xXn %13.7e \n", xXp, xXn);

   xc = x - BoxCenter[0]; // [code unit]
   yc = y - BoxCenter[1];
   zc = z - BoxCenter[2];

   radius = sqrt(xc*xc + yc*yc + zc*zc);
   radius = radius * UNIT_L; // [cm]

   // calculate heating
   dEneut = 1.544e20 * (LB_LNU/1.e52) * SQR(1.e7 / radius) * SQR(LB_TNU / 4.);

   // now subtract cooling
   T6 = (0.5*xtmp)*(0.5*xtmp)*(0.5*xtmp)*(0.5*xtmp)*(0.5*xtmp)*(0.5*xtmp);
   dEneut = dEneut - 1.399e20 * T6;

   dEneut = dEneut * exp(-xdens*1.e-11);
   dEneut = dEneut * (xXp + xXn);  // [cgs]

   xenr = xenr + dEneut * (UNIT_T * dt);

   fluid[ENGY] = (dens/(UNIT_V*UNIT_V))*(xenr + energy_shift) + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

   nuc_eos_C_short(xdens,&xtmp,ye,&xenr, &xprs, &xent, &xcs2, &xdedt, &xdpderho,
                     &xdpdrhoe, &xmunu, 0, &keyerr, rfeps); // energy mode

   // update entropy using the new energy
   fluid[ENTR] = dens * xent ;
   fluid[YE]   = dens * ye;  // lb doesn't change ye

   // if using Dual energy
#  ifdef DUAL_ENERGY
#  if (DUAL_ENERGY == DE_ENPY)
   nuc_eos_C_short(xdens,&xtmp,ye,&xenr, &xprs, &xent, &xcs2, &xdedt, &xdpderho,
                     &xdpdrhoe, &xmunu, 0, &keyerr, rfeps); // energy mode
   xprs = xprs * UNIT_P;
   fluid[ENPY] = Hydro_DensPres2Entropy( dens, xprs, Gamma_m1 );
#  endif
#  endif

#  endif // #if ( EOS == NUCLEAR )

} // FUNCTION : Src_User
