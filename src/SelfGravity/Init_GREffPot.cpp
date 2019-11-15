#include "GAMER.h"

#if ( defined GRAVITY  &&  defined GREP )


Profile_t DensAve [NLEVEL+1];
Profile_t EngyAve [NLEVEL+1];
Profile_t VrAve   [NLEVEL+1];
Profile_t PresAve [NLEVEL+1];

// debug
Profile_t DensAve_t [NLEVEL+1];
Profile_t EngyAve_t [NLEVEL+1];
Profile_t VrAve_t   [NLEVEL+1];
Profile_t PresAve_t [NLEVEL+1];
// debug


Profile_t Phi_eff [2];

static const int    INTERNAL_ENGY    = 97;
static const int    VRAD             = 98;
static const int    PRESSURE         = 99;

static const int    NProf            = 4;
static const int    Quantity [NProf] = { DENS, INTERNAL_ENGY, VRAD, PRESSURE };

static       int    lv_old           = -1;
static       double Center   [3]     = { 0.0 };
static       double MaxRadius;
static       double MinBinSize;

void CombineProfile( Profile_t Prof[], const int Quant, const bool RemoveEmpty );



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_GREffPot
// Description :  Compute the spherical-averaged profile, and GR effective potential.
//                Then set up the CPU/GPU arrays for the GR potential correction
//
// Note        :  1. Invoked by Init_GAMER() and EvolveLevel()
//                2. Enabled by the macros GRAVITY and GREP
//                3. The total averaged profile is stored at QUANT[NLEVEL]
//-------------------------------------------------------------------------------------------------------
void Init_GREffPot( const int lv )
{

// Initialize the Center, MaxRadius, and MinBinSize at the first call;
   if ( lv == -1 )
   {
      switch ( GREP_Center_Method )
      {
         case 1:   for (int i=0; i<3; i++)   Center[i] = amr->BoxCenter[i];    break;
         default:  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GREP_Center_Method", GREP_Center_Method );
      }

//    Defaults to the distance between the center and the farthest box vertex
      MaxRadius  = ( GREP_MaxRadius > 0.0 )  ? GREP_MaxRadius
                                             : SQRT( SQR( MAX( amr->BoxSize[0] - Center[0], Center[0] ) )
                                             +       SQR( MAX( amr->BoxSize[1] - Center[1], Center[1] ) )
                                             +       SQR( MAX( amr->BoxSize[2] - Center[2], Center[2] ) ));

      MinBinSize = ( GREP_MinBinSize > 0.0 ) ? GREP_MinBinSize
                                             : amr->dh[MAX_LEVEL];
   }

   int lv_bak = lv;  // debug

// Compute the spherical-averaged profile
   if ( lv == -1 )
   {
//    compute the profile at all levels at the first call
      for (int lv_t=0; lv_t<NLEVEL; lv_t++)
      {
         Profile_t *Prof_all[NProf] = { &( DensAve[lv_t] ), &( EngyAve[lv_t] ), &( VrAve[lv_t] ), &( PresAve[lv_t] ) };
         Aux_ComputeProfile( Prof_all, Center, MaxRadius, MinBinSize, GREP_LogBin, GREP_LogBinRatio,
                             false, Quantity, NProf, lv_t );
      }
   }
   else
   {
//    update the profile at the current level
      Profile_t *Prof_all[NProf] = { &( DensAve[lv] ), &( EngyAve[lv] ), &( VrAve[lv] ), &( PresAve[lv] ) };
      Aux_ComputeProfile( Prof_all, Center, MaxRadius, MinBinSize, GREP_LogBin, GREP_LogBinRatio,
                          false, Quantity, NProf, lv );

//    update velocity at level lv_old when entering finer / coarser level
//    Step 8 in EvolveLevel only corrects the data in ghost zone, and does not affect the radial profile
      if ( lv_old != lv )
      {
         Profile_t *Prof_all[1] = { &( VrAve[lv_old] ) };
         int            quant[] = { VRAD };
         Aux_ComputeProfile( Prof_all, Center, MaxRadius, MinBinSize, GREP_LogBin, GREP_LogBinRatio,
                             false, quant, 1, lv_old );
      }


// debug
      for (int lv_t=0; lv_t<NLEVEL; lv_t++)
      {
         Profile_t *Prof_all_t[NProf] = { &( DensAve_t[lv_t] ), &( EngyAve_t[lv_t] ), &( VrAve_t[lv_t] ), &( PresAve_t[lv_t] ) };
         Aux_ComputeProfile( Prof_all_t, Center, MaxRadius, MinBinSize, GREP_LogBin, GREP_LogBinRatio,
                             false, Quantity, NProf, lv_t );

      for (int b=0; b < DensAve_t[lv_t].NBin; b++)
      if ( FABS(DensAve_t[lv_t].Data[b] - DensAve[lv_t].Data[b]) / DensAve_t[lv_t].Data[b] > 1.e-5 )
         printf("Dens differs at bin %d, current lv = %d and previous level = %d, Values are %.7e  %.7e\n", b, lv, lv_old, DensAve[lv_t].Data[b], DensAve_t[lv_t].Data[b]);

      for (int b=0; b < EngyAve_t[lv_t].NBin; b++)
      if ( FABS(EngyAve_t[lv_t].Data[b] - EngyAve[lv_t].Data[b]) / EngyAve_t[lv_t].Data[b] > 1.e-5 )
         printf("Engy differs at bin %d, current lv = %d and previous level = %d, Values are %.7e  %.7e\n", b, lv, lv_old, EngyAve[lv_t].Data[b], EngyAve_t[lv_t].Data[b]);

      for (int b=0; b < VrAve_t[lv_t].NBin; b++)
      if ( FABS(VrAve_t[lv_t].Data[b] - VrAve[lv_t].Data[b]) / VrAve_t[lv_t].Data[b] > 1.e-5 )
         printf("Vrad differs at bin %d, current lv = %d and previous level = %d, Values are %.7e  %.7e\n", b, lv, lv_old, VrAve[lv_t].Data[b], VrAve_t[lv_t].Data[b]);

      for (int b=0; b < PresAve_t[lv_t].NBin; b++)
      if ( FABS(PresAve_t[lv_t].Data[b] - PresAve[lv_t].Data[b]) / PresAve_t[lv_t].Data[b] > 1.e-5 )
         printf("Pres differs at bin %d, current lv = %d and previous level = %d, Values are %.7e  %.7e\n", b, lv, lv_old, PresAve[lv_t].Data[b], PresAve_t[lv_t].Data[b]);
      }
// debug

   }


// Combine the profile at each level
   CombineProfile( DensAve, DENS,          true );
   CombineProfile( EngyAve, INTERNAL_ENGY, true );
   CombineProfile( VrAve,   VRAD,          true );
   CombineProfile( PresAve, PRESSURE,      true );

//REVISE: copy Phi_eff[1] to Phi_eff[0] to support the feature Unsplit_Gravity


// compute the GR effective potential
   CPU_ComputeEffPot( &(DensAve[NLEVEL]), &(EngyAve[NLEVEL]), &(VrAve[NLEVEL]), &(PresAve[NLEVEL]), &(Phi_eff[1]) );


// initialize the auxiliary GPU arrays
#  ifdef GPU
   CUAPI_Init_GREffPot();
#  endif

// record the level
   lv_old = lv;

} // FUNCTION : Init_GREffPot



void CombineProfile( Profile_t Prof[], const int Quant, const bool RemoveEmpty )
{

// copy bin information into Prof[NLEVEL]
   Prof[NLEVEL].NBin        = Prof[0].NBin;
   Prof[NLEVEL].LogBin      = Prof[0].LogBin;
   Prof[NLEVEL].LogBinRatio = Prof[0].LogBinRatio;
   Prof[NLEVEL].MaxRadius   = Prof[0].MaxRadius;
   Prof[NLEVEL].AllocateMemory();
   for ( int d=0; d<3; d++ )               Prof[NLEVEL].Center[d] = Prof[0].Center[d];
   for ( int b=0; b<Prof[0].NBin; b++ )
   {
      Prof[NLEVEL].Radius[b] = Prof[0].Radius[b];
      Prof[NLEVEL].Data  [b] = 0.0;
      Prof[NLEVEL].Weight[b] = 0.0;
      Prof[NLEVEL].NCell [b] = 0;
   }


// combine the profile at each level
   for (int b=0; b<Prof[0].NBin; b++)
   {
      switch ( Quant )
      {
         case DENS         :
         case ENGY         :
         case MOMX         :
         case MOMY         :
         case MOMZ         :
         case PRESSURE     :
         case INTERNAL_ENGY:
            for (int lv=0; lv<NLEVEL; lv++)
            {
               if ( Prof[lv].NCell[b] == 0L )  continue;

               Prof[NLEVEL].Data  [b] += Prof[lv].Weight[b] * Prof[lv].Data[b];
               Prof[NLEVEL].Weight[b] += Prof[lv].Weight[b];
               Prof[NLEVEL].NCell [b] += Prof[lv].NCell [b];
            }

            if ( Prof[NLEVEL].NCell[b] > 0L )
               Prof[NLEVEL].Data[b] /= Prof[NLEVEL].Weight[b];
         break;

         case VRAD:
            for (int lv=0; lv<NLEVEL; lv++)
            {
               if ( Prof[lv].NCell[b] == 0L )  continue;

               Prof[NLEVEL].Data  [b] += Prof[lv].Weight[b] * Prof[lv].Data[b];
               Prof[NLEVEL].Weight[b] += Prof[lv].Weight[b];
               Prof[NLEVEL].NCell [b] += Prof[lv].NCell [b];
            }

            if ( Prof[NLEVEL].NCell[b] > 0L  &&  Prof[NLEVEL].Weight[b] > 0.0 )
               Prof[NLEVEL].Data[b] /= Prof[NLEVEL].Weight[b];
            break;
      } // switch ( Quantity )
   } // for (int b=0; b<Prof[lv].NBin; b++)


// Remove empty bins
   if ( RemoveEmpty )
   for (int b=0; b<Prof[NLEVEL].NBin; b++)
   {
      if ( Prof[NLEVEL].NCell[b] != 0L )   continue;

//    for cases of consecutive empty bins
      int b_up;
      for (b_up=b+1; b_up<Prof[NLEVEL].NBin; b_up++)
         if ( Prof[NLEVEL].NCell[b_up] != 0L )   break;

      const int stride = b_up - b;

      for (int b_up=b+stride; b_up<Prof[NLEVEL].NBin; b_up++)
      {
         const int b_up_ms = b_up - stride;

         Prof[NLEVEL].Radius[b_up_ms] = Prof[NLEVEL].Radius[b_up];
         Prof[NLEVEL].Data  [b_up_ms] = Prof[NLEVEL].Data  [b_up];
         Prof[NLEVEL].Weight[b_up_ms] = Prof[NLEVEL].Weight[b_up];
         Prof[NLEVEL].NCell [b_up_ms] = Prof[NLEVEL].NCell [b_up];
      }

//    reset the total number of bins
      Prof[NLEVEL].NBin -= stride;

//    reduce counter since all bins above b have been shifted downward
      b --;
   } // for (int b=0; b<Prof.NBin; b++)

// reset the maximum radius
   const int b = Prof[NLEVEL].NBin;

   Prof[NLEVEL].MaxRadius = ( Prof[NLEVEL].LogBin ) ? SQR ( Prof[NLEVEL].Radius[b - 1] ) / Prof[NLEVEL].Radius[b - 2]
                                                    : 2.0 * Prof[NLEVEL].Radius[b - 1]   - Prof[NLEVEL].Radius[b - 2];

}


#endif // #if ( defined GRAVITY  &&  defined GREP )
