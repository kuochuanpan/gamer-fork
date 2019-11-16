#include "CUAPI.h"
#include "Profile.h"

#if ( defined GPU  &&  defined GRAVITY  &&  defined GREP )


#include "CUPOT.h"

extern Profile_t *Phi_eff[2];


// declare the GPU kernel requiring GREP_Data, GREP_EdgeL, GREP_Center, and r_max2
int CUPOT_SetConstMem_GREffPot( double h_GREP_Data[], double h_GREP_Edge[], double h_GREP_Center[],
                                double h_r_max2, int h_GREP_NBin );



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

   int     NBin        = Phi_eff[1]->NBin;
   double *Radius      = Phi_eff[1]->Radius;
   double  r_max2      = SQR( Phi_eff[1]->MaxRadius );
   double  Edge[NBin+1];

// check
   if ( NBin > GR_POT_NAUX_MAX )
      Aux_Error( ERROR_INFO, "Too many bins in 1D Profile %d !!\n", NBin );

// compute the location of edge
   Edge[0] = 0.0;
   for ( int i=1; i<NBin; i++ )   Edge[i] = ( Phi_eff[1]->LogBin ) ? sqrt( Radius[i - 1] * Radius[i] )
                                                                   : 0.5*( Radius[i - 1] + Radius[i] );
   Edge[NBin] = ( Phi_eff[1]->LogBin ) ? SQR ( Edge[NBin - 1] ) / Edge[NBin - 2]
                                       : 2.0 * Edge[NBin - 1]   - Edge[NBin - 2];

   int Exitcode = CUPOT_SetConstMem_GREffPot( Phi_eff[1]->Data, Edge, Phi_eff[1]->Center, r_max2, NBin );
   if (  Exitcode != 0  )
      Aux_Error( ERROR_INFO, "CUPOT_SetConstMem_GREffPot failed... Exitcode %d...\n", Exitcode );

} // FUNCTION : CUAPI_Init_GREffPot



#endif // #if ( defined GPU  &&  defined GRAVITY  &&  defined GREP )
