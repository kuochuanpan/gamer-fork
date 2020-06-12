#include "GAMER.h"

#if ( defined GRAVITY  &&  defined GREP )


Profile_t *DensAve [NLEVEL+1];
Profile_t *EngyAve [NLEVEL+1];
Profile_t *VrAve   [NLEVEL+1];
Profile_t *PresAve [NLEVEL+1];

Profile_t *Phi_eff [2];

static int    level_old        = -1;
static int    NPatch_TopLv_old = -1;
static double Center   [3]     = { 0.0 };
static double MaxRadius;
static double MinBinSize;

void CombineProfile( Profile_t *Prof[], const bool RemoveEmpty );


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_GREffPot
// Description :  Compute the spherical-averaged profile, and GR effective potential.
//                Then set up the CPU/GPU arrays for the GR potential correction
//
// Note        :  1. Invoked by Init_GAMER() and EvolveLevel()
//                2. Enabled by the macros GRAVITY and GREP
//                3. The total averaged profile is stored at QUANT[NLEVEL]
//-------------------------------------------------------------------------------------------------------
void Init_GREffPot( const int level )
{
   int NPatch_TopLv_new;


// Initialize the Center, MaxRadius, and MinBinSize at the first call;
   if ( level == -1 )
   {
      switch ( GREP_CENTER_METHOD )
      {
         case 1:   for (int i=0; i<3; i++)   Center[i] = amr->BoxCenter[i];
                   break;
         default:  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GREP_CENTER_METHOD", GREP_CENTER_METHOD );
      }

//    Defaults to the distance between the center and the farthest box vertex
      MaxRadius  = ( GREP_MAXRADIUS > 0.0 )  ? GREP_MAXRADIUS
                                             : SQRT( SQR( MAX( amr->BoxSize[0] - Center[0], Center[0] ) )
                                             +       SQR( MAX( amr->BoxSize[1] - Center[1], Center[1] ) )
                                             +       SQR( MAX( amr->BoxSize[2] - Center[2], Center[2] ) ));

      MinBinSize = ( GREP_MINBINSIZE > 0.0 ) ? GREP_MINBINSIZE
                                             : amr->dh[MAX_LEVEL];
   }


// Compute the spherical-averaged profile
   if ( level == -1 )
   {
//    Initialize pointer array of Profile
      for (int lv=0; lv<NLEVEL+1; lv++)
      {
         DensAve[lv] = new Profile_t();
         EngyAve[lv] = new Profile_t();
         VrAve  [lv] = new Profile_t();
         PresAve[lv] = new Profile_t();
      }

      for (int PROFID=0; PROFID<2; PROFID++)
         Phi_eff[PROFID] = new Profile_t();

//    compute the profile at all levels at the first call
      for (int lv=0; lv<NLEVEL; lv++)
      {
         long       TVar [] = {       _DENS,     _VELR,       _PRES,   _EINT_DER };
         Profile_t *Prof [] = { DensAve[lv], VrAve[lv], PresAve[lv], EngyAve[lv] };

         Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO,
                             false, TVar, 4, lv );
      }
   }

   else
   {

/*
      if ( level == 0 )
      {
         for (int lv=0; lv<NLEVEL; lv++)
         {
            long       TVar [] = {       _DENS,     _VELR,       _PRES,   _EINT_DER };
            Profile_t *Prof [] = { DensAve[lv], VrAve[lv], PresAve[lv], EngyAve[lv] };

            Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO,
                                false, TVar, 4, lv );
         }
      }
*/

//    update the profile at the current level
      {
         int          lv    = level;
//         long       TVar [] = {       _DENS,     _VELR,       _PRES,   _EINT_DER };
//         Profile_t *Prof [] = { DensAve[lv], VrAve[lv], PresAve[lv], EngyAve[lv] };
         long       TVar [] = {       _DENS,   _EINT_DER,     _VELR,       _PRES };
         Profile_t *Prof [] = { DensAve[lv], EngyAve[lv], VrAve[lv], PresAve[lv] };

         Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO,
                             false, TVar, 4, lv );
      }

//    update profiles related to momentum and energy at last level when entering finer / coarser level
      if ( level_old != -1  &&  level_old != level )
      {
//       collect the number of Patches at TOP_LEVEL
         MPI_Allreduce( &amr->NPatchComma[TOP_LEVEL][1], &NPatch_TopLv_new, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

         if ( level_old == TOP_LEVEL  &&  NPatch_TopLv_new != NPatch_TopLv_old )
         {
//          also update the density if new Patches are created in TOP_LEVEL (Step 9 in EvolveLevel()).
//          For new patches created in other level, the profile will be updated in code below.
            int          lv    = TOP_LEVEL;
            long       TVar [] = {       _DENS,     _VELR,       _PRES,   _EINT_DER };
            Profile_t *Prof [] = { DensAve[lv], VrAve[lv], PresAve[lv], EngyAve[lv] };

            Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO,
                                false, TVar, 4, lv );
         }
         else
         {
            int          lv    = level_old;
            long       TVar [] = {     _VELR,       _PRES,   _EINT_DER };
            Profile_t *Prof [] = { VrAve[lv], PresAve[lv], EngyAve[lv] };

            Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO,
                                false, TVar, 3, lv );
         }
      }

//    update correction from finer level in Step 8 of EvolveLevel()
//           new patches                    Step 9
      if ( level_old > level )
      {
         for (int lv=level+1; lv<TOP_LEVEL; lv++)
         {
            long       TVar [] = {       _DENS,     _VELR,       _PRES,   _EINT_DER };
            Profile_t *Prof [] = { DensAve[lv], VrAve[lv], PresAve[lv], EngyAve[lv] };

            Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO,
                                false, TVar, 4, lv );
         }
      }

   }


// Combine the profile at each level
   CombineProfile( DensAve, true );
   CombineProfile( EngyAve, true );
   CombineProfile( VrAve,   true );
   CombineProfile( PresAve, true );

//REVISE: copy Phi_eff[1] to Phi_eff[0] to support the feature Unsplit_Gravity

// compute the GR effective potential
   CPU_ComputeEffPot( DensAve[NLEVEL], EngyAve[NLEVEL], VrAve[NLEVEL], PresAve[NLEVEL], Phi_eff[1] );

// initialize the auxiliary GPU arrays
#  ifdef GPU
   CUAPI_Init_GREffPot();
#  endif

// record the level and number of patches in finest level
   if ( level_old != level )   NPatch_TopLv_old = NPatch_TopLv_new;
   level_old = level;

} // FUNCTION : Init_GREffPot



void CombineProfile( Profile_t *Prof[], const bool RemoveEmpty )
{

// copy bin information into Prof[NLEVEL]
   Prof[NLEVEL]->NBin        = Prof[0]->NBin;
   Prof[NLEVEL]->LogBin      = Prof[0]->LogBin;
   Prof[NLEVEL]->LogBinRatio = Prof[0]->LogBinRatio;
   Prof[NLEVEL]->MaxRadius   = Prof[0]->MaxRadius;
   Prof[NLEVEL]->AllocateMemory();

   for ( int d=0; d<3; d++ )
      Prof[NLEVEL]->Center[d] = Prof[0]->Center[d];

   for ( int b=0; b<Prof[0]->NBin; b++ )
   {
      Prof[NLEVEL]->Radius[b] = Prof[0]->Radius[b];
      Prof[NLEVEL]->Data  [b] = 0.0;
      Prof[NLEVEL]->Weight[b] = 0.0;
      Prof[NLEVEL]->NCell [b] = 0;
   }


// combine the profile at each level
   for (int b=0; b<Prof[0]->NBin; b++)
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
         if ( Prof[lv]->NCell[b] == 0L )  continue;

         Prof[NLEVEL]->Data  [b] += Prof[lv]->Weight[b] * Prof[lv]->Data[b];
         Prof[NLEVEL]->Weight[b] += Prof[lv]->Weight[b];
         Prof[NLEVEL]->NCell [b] += Prof[lv]->NCell [b];
      }

      if ( Prof[NLEVEL]->NCell[b] > 0L  &&  Prof[NLEVEL]->Weight[b] > 0.0 )
         Prof[NLEVEL]->Data[b] /= Prof[NLEVEL]->Weight[b];
   } // for (int b=0; b<Prof[0]->NBin; b++)


// remove the empty bins
   if ( RemoveEmpty )
   {
      for (int b=0; b<Prof[NLEVEL]->NBin; b++)
      {
         if ( Prof[NLEVEL]->NCell[b] != 0L )   continue;

   //    for cases of consecutive empty bins
         int b_up;
         for (b_up=b+1; b_up<Prof[NLEVEL]->NBin; b_up++)
            if ( Prof[NLEVEL]->NCell[b_up] != 0L )   break;

         const int stride = b_up - b;

         for (int b_up=b+stride; b_up<Prof[NLEVEL]->NBin; b_up++)
         {
            const int b_up_ms = b_up - stride;

            Prof[NLEVEL]->Radius[b_up_ms] = Prof[NLEVEL]->Radius[b_up];
            Prof[NLEVEL]->Data  [b_up_ms] = Prof[NLEVEL]->Data  [b_up];
            Prof[NLEVEL]->Weight[b_up_ms] = Prof[NLEVEL]->Weight[b_up];
            Prof[NLEVEL]->NCell [b_up_ms] = Prof[NLEVEL]->NCell [b_up];
         }

   //    reset the total number of bins
         Prof[NLEVEL]->NBin -= stride;

   //    reduce counter since all bins above b have been shifted downward
         b --;
      } // for (int b=0; b<Prof[NLEVEL]->NBin; b++)

//    update the maximum radius since the last bin may have not been removed
      const int LastBin = Prof[NLEVEL]->NBin-1;

      Prof[NLEVEL]->MaxRadius = ( Prof[NLEVEL]->LogBin ) ? Prof[NLEVEL]->Radius[LastBin] * sqrt( Prof[NLEVEL]->LogBinRatio )
                                                         : Prof[NLEVEL]->Radius[LastBin] + 0.5*MinBinSize;
   }

} // FUNCTION : CombineProfile



#endif // #if ( defined GRAVITY  &&  defined GREP )
