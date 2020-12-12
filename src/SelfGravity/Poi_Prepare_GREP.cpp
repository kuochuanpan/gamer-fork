#include "GAMER.h"

#if ( defined GRAVITY  &&  defined GREP )


static void Update_GREP_Profile( const int lv, const int Sg, const double PrepTime );
static void Combine_GREP_Profile( Profile_t *Prof[][2], const int lv, const int Sg, const double PrepTime,
                                  const bool RemoveEmpty );
extern void SetTempIntPara( const int lv, const int Sg_Current, const double PrepTime, const double Time0, const double Time1,
                            bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT );


Profile_t *DensAve [NLEVEL+1][2];
Profile_t *EngyAve [NLEVEL+1][2];
Profile_t *VrAve   [NLEVEL+1][2];
Profile_t *PresAve [NLEVEL+1][2];
Profile_t *Phi_eff [NLEVEL  ][2];

int    GREP_LvUpdate;
double GREP_Prof_MaxRadius;
double GREP_Prof_MinBinSize;

int    GREPSg     [NLEVEL];
double GREPSgTime [NLEVEL][2];
double GREP_Prof_Center   [3];


// temporary switch for test different schemes for temporal interpolation:
//   True : apply temporal interpolation in Aux_ComputeProfile()
//   False: apply temporal interpolation when combining the stored profiles
static bool Do_TEMPINT_in_ComputeProfile = true;




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_Prepare_GREP
// Description :  Compute the spherical-averaged profiles and GR effective potential.
//
// Note        :  1. Enabled if macro GRAVITY and GREP are set
//                2. Invoked by Poi_UserWorkBeforePoisson_GREP()
//                3. The GREP evalutaed at each level is stored at Phi_eff[level]
//-------------------------------------------------------------------------------------------------------
void Poi_Prepare_GREP( const double Time, const int lv )
{

// compare the input Time with stored time to choose the suitable SaveSg
   int Sg;

   if      (  Mis_CompareRealValue( Time, GREPSgTime[lv][0], NULL, false )  )   Sg = 0;
   else if (  Mis_CompareRealValue( Time, GREPSgTime[lv][1], NULL, false )  )   Sg = 1;
   else                                                                         Sg = 1 - GREPSg[lv];


// update the level, Sg, and SgTime
   GREP_LvUpdate      = lv;
   GREPSg    [lv]     = Sg;
   GREPSgTime[lv][Sg] = Time;


// update the spherical-averaged profiles
   if ( Do_TEMPINT_in_ComputeProfile )   Update_GREP_Profile( lv, Sg, Time );
   else                                  Update_GREP_Profile( lv, Sg, -1.0 );


// combine the profile at each level
   Combine_GREP_Profile( DensAve, lv, Sg, Time, true );
   Combine_GREP_Profile( EngyAve, lv, Sg, Time, true );
   Combine_GREP_Profile( VrAve,   lv, Sg, Time, true );
   Combine_GREP_Profile( PresAve, lv, Sg, Time, true );


// compute the effective GR potential
   CPU_ComputeGREP( DensAve[NLEVEL][Sg], EngyAve[NLEVEL][Sg], VrAve[NLEVEL][Sg], PresAve[NLEVEL][Sg],
                    Phi_eff[lv]    [Sg] );

} // FUNCTION : Poi_Prepare_GREP



//-------------------------------------------------------------------------------------------------------
// Function    :  Update_GREP_Profile
// Description :  Update the spherical-averaged profiles
//
// Note        :  1. The contribution from     leaf patches on level = lv (<= lv) is stored QUANT[    lv]
//                                         non-leaf patches on level = lv         is stored QUANT[NLEVEL]
//                2. Apply temporal interpolation in Aux_ComputeProfile if PrepTime >= 0
//                3. To avoid inconsistent leaf- and non-leaf profiles during combining,
//                   we retain the empty bins in the profiles obtained here because
//                   they could have different empty bins.
//
//-------------------------------------------------------------------------------------------------------
static void Update_GREP_Profile( const int lv, const int Sg, const double PrepTime )
{

   long       TVar         [] = {               _DENS,             _VELR,               _PRES,           _EINT_DER };
   Profile_t *Prof_Leaf    [] = { DensAve[    lv][Sg], VrAve[    lv][Sg], PresAve[    lv][Sg], EngyAve[    lv][Sg] };
   Profile_t *Prof_NonLeaf [] = { DensAve[NLEVEL][Sg], VrAve[NLEVEL][Sg], PresAve[NLEVEL][Sg], EngyAve[NLEVEL][Sg] };


   if ( Do_TEMPINT_in_ComputeProfile )
//###CHECK: does leaf patch transit to non-leaf path at sub-cycling?
//    update the profile from leaf patches on level <= lv
      Aux_ComputeProfile   ( Prof_Leaf,     GREP_Prof_Center, GREP_Prof_MaxRadius, GREP_Prof_MinBinSize,
                             GREP_LOGBIN,   GREP_LOGBINRATIO, false, TVar, 4, -1, lv, PATCH_LEAF,    PrepTime );

   else
   {
//    update the profile from leaf patches on level = lv
      Aux_ComputeProfile   ( Prof_Leaf,     GREP_Prof_Center, GREP_Prof_MaxRadius, GREP_Prof_MinBinSize,
                             GREP_LOGBIN,   GREP_LOGBINRATIO, false, TVar, 4, lv, -1, PATCH_LEAF,    PrepTime );

//###CHECK: does the refinment correction affect leaf patch? If no, this part is not necesary.
//    update the USG profile from leaf patches on level = lv to account the correction from finer level
      if ( ( lv < TOP_LEVEL )  &&  ( GREPSgTime[lv][1 - Sg] >= 0.0 ) )
      {
         int               Sg_USG    = 1 - Sg;
         double          Time_USG    = GREPSgTime[lv][Sg_USG];
         Profile_t *Prof_Leaf_USG [] = { DensAve[lv][Sg_USG], VrAve[lv][Sg_USG], PresAve[lv][Sg_USG], EngyAve[lv][Sg_USG] };

         Aux_ComputeProfile( Prof_Leaf_USG, GREP_Prof_Center, GREP_Prof_MaxRadius, GREP_Prof_MinBinSize,
                             GREP_LOGBIN,   GREP_LOGBINRATIO, false, TVar, 4, lv, -1, PATCH_LEAF,    Time_USG );
      }
   }


// update the profile from the non-leaf patches on level = lv
   Aux_ComputeProfile      ( Prof_NonLeaf,  GREP_Prof_Center, GREP_Prof_MaxRadius, GREP_Prof_MinBinSize,
                             GREP_LOGBIN,   GREP_LOGBINRATIO, false, TVar, 4, lv, -1, PATCH_NONLEAF, PrepTime );

} // FUNCTION : Update_GREP_Profile



//-------------------------------------------------------------------------------------------------------
// Function    :  Combine_GREP_Profile
// Description :  Combine the separated spherical-averaged profiles
//                and handling the empty bins in the combined profile
//
// Note        :  1. The total averaged profile is stored at QUANT[NLEVEL]
//
// Parameter   :  Prof        : Profile_t object array to be combined
//                RemoveEmpty : true  --> remove empty bins from the data
//                              false --> these empty bins will still be in the profile arrays with
//                                        Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0
//-------------------------------------------------------------------------------------------------------
void Combine_GREP_Profile( Profile_t *Prof[][2], const int lv, const int Sg, const double PrepTime,
                           const bool RemoveEmpty )
{

   Profile_t *Prof_NonLeaf = Prof[NLEVEL][Sg];


// combine the contributions from leaf patches on level <= lv and non-leaf patches on level = lv,
// and store the sum in 'Prof_NonLeaf'
   if ( Do_TEMPINT_in_ComputeProfile )
   {
//    temporal interpolation is done in Aux_ComputeProfile()
      Profile_t *Prof_Leaf = Prof[lv][Sg];

      for (int b=0; b<Prof_Leaf->NBin; b++)
      {
         if ( Prof_Leaf->NCell[b] == 0L )  continue;

         Prof_NonLeaf->Data  [b]  = Prof_NonLeaf->Weight[b] * Prof_NonLeaf->Data[b]
                                  + Prof_Leaf   ->Weight[b] * Prof_Leaf   ->Data[b];
         Prof_NonLeaf->Weight[b] += Prof_Leaf   ->Weight[b];
         Prof_NonLeaf->NCell [b] += Prof_Leaf   ->NCell [b];

         Prof_NonLeaf->Data  [b] /= Prof_NonLeaf->Weight[b];
      }
   }

   else
   {
//    apply temporal interpolation to leaf patches on level <= lv
      for (int level=0; level<=lv; level++)
      {
//       temporal interpolation parameters
         bool FluIntTime;
         int  FluSg, FluSg_IntT;
         real FluWeighting, FluWeighting_IntT;

         SetTempIntPara( level, GREPSg[level], PrepTime, GREPSgTime[level][0], GREPSgTime[level][1],
                         FluIntTime, FluSg, FluSg_IntT, FluWeighting, FluWeighting_IntT );

         Profile_t *Prof_Leaf      = Prof[level][FluSg];
         Profile_t *Prof_Leaf_IntT = ( FluIntTime ) ? Prof[level][FluSg_IntT] : NULL;


//REVISE: If the profile center is allowed to change with time in future,
//        the number of bin and its location could be different between Prof_Leaf and Prof_Leaf_IntT
         for (int b=0; b<Prof_Leaf->NBin; b++)
         {
            if ( Prof_Leaf->NCell[b] == 0L )  continue;

            Prof_NonLeaf->Data  [b]  = ( FluIntTime )
                                     ?                       Prof_NonLeaf  ->Weight[b] * Prof_NonLeaf  ->Data[b]
                                       + FluWeighting      * Prof_Leaf     ->Weight[b] * Prof_Leaf     ->Data[b]
                                       + FluWeighting_IntT * Prof_Leaf_IntT->Weight[b] * Prof_Leaf_IntT->Data[b]
                                     :                       Prof_NonLeaf  ->Weight[b] * Prof_NonLeaf  ->Data[b]
                                       +                     Prof_Leaf     ->Weight[b] * Prof_Leaf     ->Data[b];

            Prof_NonLeaf->Weight[b] += ( FluIntTime )
                                     ?   FluWeighting      * Prof_Leaf     ->Weight[b]
                                       + FluWeighting_IntT * Prof_Leaf_IntT->Weight[b]
                                     :                       Prof_Leaf     ->Weight[b];

            Prof_NonLeaf->NCell [b] += Prof_Leaf   ->NCell [b];
            Prof_NonLeaf->Data  [b] /= Prof_NonLeaf->Weight[b];
         } // for (int b=0; b<Prof_Leaf->NBin; b++)
      } // for (int level=0; level<=lv; level++)
   } // if ( Do_TEMPINT_in_ComputeProfile )


// remove the empty bins in the combined profile stored in 'Prof_NonLeaf'
   if ( RemoveEmpty )
   {
      for (int b=0; b<Prof_NonLeaf->NBin; b++)
      {
         if ( Prof_NonLeaf->NCell[b] != 0L )   continue;

//       for cases of consecutive empty bins
         int b_up;
         for (b_up=b+1; b_up<Prof_NonLeaf->NBin; b_up++)
            if ( Prof_NonLeaf->NCell[b_up] != 0L )   break;

         const int stride = b_up - b;

         for (int b_up=b+stride; b_up<Prof_NonLeaf->NBin; b_up++)
         {
            const int b_up_ms = b_up - stride;

            Prof_NonLeaf->Radius[b_up_ms] = Prof_NonLeaf->Radius[b_up];
            Prof_NonLeaf->Data  [b_up_ms] = Prof_NonLeaf->Data  [b_up];
            Prof_NonLeaf->Weight[b_up_ms] = Prof_NonLeaf->Weight[b_up];
            Prof_NonLeaf->NCell [b_up_ms] = Prof_NonLeaf->NCell [b_up];
         }

//       reset the total number of bins
         Prof_NonLeaf->NBin -= stride;
      } // for (int b=0; b<Prof_NonLeaf->NBin; b++)

//    update the maximum radius since the last bin may have not been removed
      const int LastBin = Prof_NonLeaf->NBin-1;

      Prof_NonLeaf->MaxRadius = ( Prof_NonLeaf->LogBin )
                              ? Prof_NonLeaf->Radius[LastBin] * sqrt( Prof_NonLeaf->LogBinRatio )
                              : Prof_NonLeaf->Radius[LastBin] + 0.5*GREP_Prof_MinBinSize;
   } // if ( RemoveEmpty )

} // FUNCTION : Combine_GREP_Profile



#endif // #if ( defined GRAVITY  &&  defined GREP )
