#include "GAMER.h"

#if ( defined GRAVITY  &&  defined GREP )


extern void SetExtPotAuxArray_GREP( double AuxArray[] );


extern Profile_t *Phi_eff [NLEVEL][2];
extern int        GREPSg  [NLEVEL];
extern int        GREP_LvUpdate;

extern double *h_GREP_Lv_Data_New;
extern double *h_GREP_FaLv_Data_New;
extern double *h_GREP_FaLv_Data_Old;
extern double *h_GREP_Lv_Radius_New;
extern double *h_GREP_FaLv_Radius_New;
extern double *h_GREP_FaLv_Radius_Old;
extern int     h_GREP_Lv_NBin_New;
extern int     h_GREP_FaLv_NBin_New;
extern int     h_GREP_FaLv_NBin_Old;




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_UserWorkBeforePoisson_GREP
// Description :  Compute the GREP, transfer data to GPU device, and update CPU/GPU data pointer
//                before invoking the Poisson solver
//
// Note        :  1. Invoked by Gra_AdvanceDt() using the function pointer "Poi_UserWorkBeforePoisson_Ptr"
//
// Parameter   :  Time : Target physical time
//                lv   : Target refinement level
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Poi_UserWorkBeforePoisson_GREP( const double Time, const int lv )
{

// compute GREP
   Poi_Prepare_GREP( Time, lv );


// update the auxiliary arrays for GREP
   SetExtPotAuxArray_GREP( ExtPot_AuxArray );


// update the CPU pointer
   const int Lv   = GREP_LvUpdate;
   const int FaLv = ( Lv > 0 ) ? Lv - 1 : Lv;

   const int Sg_Lv   = GREPSg[Lv];
   const int Sg_FaLv = GREPSg[FaLv];

   Profile_t *Phi_Lv_New   = Phi_eff[ Lv   ][     Sg_Lv   ];
   Profile_t *Phi_FaLv_New = Phi_eff[ FaLv ][     Sg_FaLv ];
   Profile_t *Phi_FaLv_Old = Phi_eff[ FaLv ][ 1 - Sg_FaLv ];

   h_GREP_Lv_Data_New     = Phi_Lv_New  ->Data;
   h_GREP_FaLv_Data_New   = Phi_FaLv_New->Data;
   h_GREP_FaLv_Data_Old   = Phi_FaLv_Old->Data;
   h_GREP_Lv_Radius_New   = Phi_Lv_New  ->Radius;
   h_GREP_FaLv_Radius_New = Phi_FaLv_New->Radius;
   h_GREP_FaLv_Radius_Old = Phi_FaLv_Old->Radius;
   h_GREP_Lv_NBin_New     = Phi_Lv_New  ->NBin;
   h_GREP_FaLv_NBin_New   = Phi_FaLv_New->NBin;
   h_GREP_FaLv_NBin_Old   = Phi_FaLv_Old->NBin;


// update the auxiliary GPU arrays
#  ifdef GPU
   CUAPI_SetConstMemory_ExtAccPot();
   CUAPI_SetConstMemory_ExtPot_GREP();
#  endif

} // FUNCTION : Poi_UserWorkBeforePoisson_GREP



#endif // #if ( defined GRAVITY  &&  defined GREP )
