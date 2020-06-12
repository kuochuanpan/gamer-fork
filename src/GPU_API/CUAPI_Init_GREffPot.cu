#include "CUAPI.h"
#include "Profile.h"

#if ( defined GPU  &&  defined GRAVITY  &&  defined GREP )


#include "CUPOT.h"

extern Profile_t *Phi_eff[2];


// declare the GPU kernel requiring GREP_Data, GREP_EdgeL, GREP_Center, and r_max2
int CUPOT_SetConstMem_GREffPot( double h_GREP_Data[], double h_GREP_Radius[], double h_GREP_Center[],
                                int    h_GREP_NBin );



//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Init_GREffPot
// Description :  Set the auxiliary GPU constant-memory arrays for the GR effective potential

// Note        :  1. Invoked by Init_GREffPot()
//
// Parameter   :  None
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void CUAPI_Init_GREffPot()
{

// check
   if ( Phi_eff[1]->NBin > GR_POT_NAUX_MAX )
      Aux_Error( ERROR_INFO, "Too many bins in average radial Profile %d !!\n", Phi_eff[1]->NBin );


   int Exitcode = CUPOT_SetConstMem_GREffPot( Phi_eff[1]->Data, Phi_eff[1]->Radius, Phi_eff[1]->Center, Phi_eff[1]->NBin );
   if (  Exitcode != 0  )
      Aux_Error( ERROR_INFO, "CUPOT_SetConstMem_GREffPot failed... Exitcode %d...\n", Exitcode );

} // FUNCTION : CUAPI_Init_GREffPot



#endif // #if ( defined GPU  &&  defined GRAVITY  &&  defined GREP )
