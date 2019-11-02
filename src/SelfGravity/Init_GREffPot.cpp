#include "GAMER.h"

#if ( defined GRAVITY  &&  defined GREP )


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_GREffPot
// Description :  Set the auxiliary CPU/GPU arrays for the GR potential correction
//
// Note        :  1. Invoked by Init_GAMER() and EvolveLevel()
//                2. Enabled by the macros GRAVITY and GREP
//-------------------------------------------------------------------------------------------------------
void Init_GREffPot()
{

// compute the GR effective potential
   CPU_ComputeEffPot();

// initialize the auxiliary GPU arrays
#  ifdef GPU
   CUAPI_Init_GREffPot();
#  endif

} // FUNCTION : Init_GREffPot


#endif // #if ( defined GRAVITY  &&  defined GREP )
