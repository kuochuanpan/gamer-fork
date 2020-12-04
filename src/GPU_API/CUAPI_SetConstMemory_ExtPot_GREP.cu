#include "Macro.h"
#include "Profile.h"
#include "CUAPI.h"
#include "CUDA_ConstMemory.h"

#if ( defined GPU  &&  defined GRAVITY  &&  defined GREP )


extern int        GREP_LvUpdate;
extern int        GREPSg  [NLEVEL];
extern Profile_t *Phi_eff [NLEVEL][2];




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SetConstMemory_ExtPot_GREP
// Description :  Set the constant memory variables on GPU used by GREP
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by Poi_UserWorkBeforePoisson_GREP()
//                3. EXT_POT_GREP_NAUX_MAX is defined in Macro.h (default = 4000)
//
// Return      :  c_GREP_Lv_Data_New[], c_GREP_Radius[], c_GREP_NBin
//-------------------------------------------------------------------------------------------------------
void CUAPI_SetConstMemory_ExtPot_GREP()
{

   Profile_t *Phi_Lv_New = Phi_eff[GREP_LvUpdate][ GREPSg[GREP_LvUpdate] ];


   if ( Phi_Lv_New->NBin > EXT_POT_GREP_NAUX_MAX )
      Aux_Error( ERROR_INFO, "Too many bins in GREP profiles %d !!\n", Phi_Lv_New->NBin );


   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_GREP_Lv_Data_New,    Phi_Lv_New->Data,   EXT_POT_GREP_NAUX_MAX*sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_GREP_Lv_Radius_New,  Phi_Lv_New->Radius, EXT_POT_GREP_NAUX_MAX*sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_GREP_Lv_NBin_New,   &Phi_Lv_New->NBin,                         sizeof(int   ) )  );

} // FUNCTION : CUAPI_SetConstMemory_ExtPot_GREP



#endif // #if ( defined GPU  &&  defined GRAVITY  &&  defined GREP )
