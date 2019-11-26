#include "GAMER.h"
//#include "NuclearEos.h"


#if ( defined GRAVITY  &&  defined GREP )


static const double FourPI        = 4.0*M_PI;
static const double FourThirdPI   = FourPI/3.0;
static const double tolerance     = 1.e-5;


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_ComputeEffPot
// Description :  Construct the effective potential
//
// Note        :  1. Phi_eff is declared in NuclearEos.h and initialized here
//                2. Enabled if macro GRAVITY and GREP are set
//                3. The profile Phi_eff store the value of -Phi_NW(r) + Phi_TOV(r)
//                   at the left edge of bins
//-------------------------------------------------------------------------------------------------------
void CPU_ComputeEffPot( Profile_t *DensAve, Profile_t *EngyAve, Profile_t *VrAve, Profile_t *PresAve,
                        Profile_t *Phi_eff )
{

// iteratively construct the m_TOV and Gamma
   const double c2            = SQR( Const_c/UNIT_V );

   double  dMass;
   int     NIter              = GREP_MaxIter;
   int     NBin               = DensAve->NBin;
   double *Radius             = DensAve->Radius;
   double  Mass_TOV_USG[NBin] = { 0.0 };
   double  Mass_TOV    [NBin] = { 0.0 }; // store mass for \bar_Phi(r)     in Eq. (7) in Marek+ (2006)
   double  Mass_NW     [NBin] = { 0.0 }; // store mass for \bar_Phi(r)_TOV in Eq. (7) in Marek+ (2006)
   double  Gamma_TOV   [NBin] = { 1.0 }; // set initial guess to 1.0
   double  EdgeL       [NBin] = { 0.0 };
   double  EdgeLCubed  [NBin] = { 0.0 };
   double  RadiusCubed [NBin] = { 0.0 };

   for ( int i=1; i<NBin; i++ )   EdgeL[i] = ( GREP_LogBin ) ? sqrt( Radius[i - 1] * Radius[i] )
                                                             : 0.5*( Radius[i - 1] + Radius[i] );

   for ( int i=0; i<NBin; i++ )
   {
      EdgeLCubed [i] = CUBE( EdgeL [i] );
      RadiusCubed[i] = CUBE( Radius[i] );
   }


// construct Mass_NW
   Mass_NW[0] = FourThirdPI * RadiusCubed[0] * DensAve->Data[0];

   for ( int i=1; i<NBin-1; i++ )
   {
      dMass = FourThirdPI * ( EdgeLCubed [i] - RadiusCubed[i-1] ) * DensAve->Data[i-1]
            + FourThirdPI * ( RadiusCubed[i] - EdgeLCubed [i]   ) * DensAve->Data[i];

      Mass_NW[i] = Mass_NW[i-1] + dMass;
   }

   dMass = FourThirdPI * ( EdgeLCubed[NBin-1] - RadiusCubed[NBin-2] ) * DensAve->Data[NBin-2];
   Mass_NW[NBin-1] = Mass_NW[NBin-2] + dMass;


// construct Mass_TOV and Gamma_TOV
   while ( NIter )
   {
      NIter -= 1;

      if ( !NIter )   Aux_Error( ERROR_INFO, "Too many iterations in effective potentia\n" );

//    construct Mass_TOV
      Mass_TOV[0] = FourThirdPI * RadiusCubed[0] * Gamma_TOV[0]
                  * ( DensAve->Data[0] + EngyAve->Data[0] / c2 );

      for ( int i=1; i<NBin-1; i++ )
      {
         dMass = FourThirdPI * ( EdgeLCubed [i] - RadiusCubed[i-1] ) * Gamma_TOV[i-1]
               * ( DensAve->Data[i-1] + EngyAve->Data[i-1] / c2 )
               + FourThirdPI * ( RadiusCubed[i] - EdgeLCubed [i] ) * Gamma_TOV[i]
               * ( DensAve->Data[i]   + EngyAve->Data[i]   / c2 );

         Mass_TOV[i] = Mass_TOV[i-1] + dMass;
      }

      dMass = FourThirdPI * ( EdgeLCubed[NBin-1] - RadiusCubed[NBin-2] )
            * ( DensAve->Data[NBin-2] + EngyAve->Data[NBin-2] / c2 ) * Gamma_TOV[NBin-2];
      Mass_TOV[NBin-1] = Mass_TOV[NBin-2] + dMass;


//    compute Gamma_TOV
      for ( int i=0; i<NBin-1; i++ )
//       to avoid cases of negative value
         Gamma_TOV[i] = SQRT( MAX( 0.0,
                                   1.0 + ( SQR( VrAve->Data[i] ) - 2.0*NEWTON_G*Mass_TOV[i]/Radius[i] ) / c2 ) );

      Gamma_TOV[NBin-1] = SQRT( MAX( 0.0,
                                     1.0 +
                                     ( SQR( VrAve->Data[NBin-1] ) - 2.0*NEWTON_G*Mass_TOV[NBin-1]/EdgeL[NBin-1] ) / c2 ) );

//    check if tolerance is satisfied
      double error_max = 0.0;
      for ( int i=0; i<NBin; i++ )
      {
         double error = FABS( Mass_TOV_USG[i] - Mass_TOV[i] ) / Mass_TOV[i];
         if ( error_max <= error )   error_max = error;
      }

      if ( error_max <= tolerance )
         break;
      else
         for ( int i=0; i<NBin; i++ )   Mass_TOV_USG[i] = Mass_TOV[i];
   } // while ( NIter )


// copy bin information into Phi_eff
   Phi_eff->NBin        = DensAve->NBin;
   Phi_eff->LogBin      = DensAve->LogBin;
   Phi_eff->LogBinRatio = DensAve->LogBinRatio;
   Phi_eff->MaxRadius   = DensAve->MaxRadius;
   Phi_eff->AllocateMemory();
   for ( int d=0; d<3; d++ )              Phi_eff->Center[d] = DensAve->Center[d];
   for ( int b=0; b<Phi_eff->NBin; b++ )   Phi_eff->Radius[b] = DensAve->Radius[b];


// compute effective potential at the left edge
   Phi_eff->Data[NBin-1] = -NEWTON_G * (Mass_TOV[NBin-1] - Mass_NW[NBin-1]) / EdgeL[NBin-1];

   for ( int i=NBin-2; i>0; i-- )
   {
      double dr   = EdgeL[i+1] - EdgeL[i];
      double dPhi = -dr * NEWTON_G
                  * ( Mass_TOV[i] + FourPI * RadiusCubed[i] * PresAve->Data[i] / c2 )
                  / SQR( Radius[i] * Gamma_TOV[i] );

      dPhi *= 1.0 + ( EngyAve->Data[i] + PresAve->Data[i] ) / ( DensAve->Data[i] * c2 );
      dPhi -= -dr * NEWTON_G * Mass_NW[i] / SQR( Radius[i] );

      Phi_eff->Data[i] = Phi_eff->Data[i+1] + dPhi;
   }

   Phi_eff->Data[0] = Phi_eff->Data[1];


#ifdef GREP_DEBUG
   printf("\n# GREP_Center_Method: %d\n",                  GREP_Center_Method);
   printf("# Center              : %.15e\t%.15e\t%.15e\n", Phi_eff->Center[0], Phi_eff->Center[1], Phi_eff->Center[2]);
   printf("# MaxRadius           : %.15e\n",               Phi_eff->MaxRadius);
//   printf("# MinBinSize          : %.15e\n",               MinBinSize);
   printf("# LogBin              : %d\n",                  Phi_eff->LogBin);
   printf("# LogBinRatio         : %.15e\n",               Phi_eff->LogBinRatio);
   printf("# Num of Iteration    : %d\n",                  GREP_MaxIter - NIter);
   printf("# ============================================================");
   printf("# Profile info: NBin = %d\n", NBin);
   printf("# -- Bin -- NCell -- RADIUS -- DENS -- ENGY -- Vr -- Pressure -- Mass_NW -- Mass_TOV -- Gamma_TOV -- Eff_Pot --\n");
   for ( int i=0; i<NBin; i++ )
      printf("%5d\t%8ld\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
             i, DensAve->NCell[i], Radius[i], DensAve->Data[i], EngyAve->Data[i], VrAve->Data[i], PresAve->Data[i],
             Mass_NW[i], Mass_TOV[i], Gamma_TOV[i], Phi_eff->Data[i]);
#endif

} // FUNCTION : CPU_ComputeEffPot


#endif // #if ( defined GRAVITY  &&  defined GREP )
