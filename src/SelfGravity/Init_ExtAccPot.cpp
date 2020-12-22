#include "GAMER.h"

#ifdef GRAVITY


// prototypes of built-in ExtPot
void Init_ExtPot_Tabular();
void Init_ExtPot_GREP();

// these function pointers must be set by a test problem initializer
void (*Init_ExtAcc_Ptr)() = NULL;
void (*Init_ExtPot_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtAccPot
// Description :  Initialize external acceleration and potential
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Enabled by the runtime options "OPT__EXT_ACC" and "OPT__EXT_POT"
//                3. Function pointers Init_ExtAcc_Ptr and Init_ExtPot_Ptr must be set in advance by a
//                   test problem initializer
//                4. Must invoke either CUAPI_SetConstMemory() or CUAPI_SetConstMemory_ExtAccPot() afterward
//                   to set the GPU constant memory
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExtAccPot()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// external acceleration
   if ( OPT__EXT_ACC )
   {
      if ( Init_ExtAcc_Ptr != NULL )   Init_ExtAcc_Ptr();
      else                             Aux_Error( ERROR_INFO, "Init_ExtAcc_Ptr == NULL !!\n" );
   }

// external potential
   if ( OPT__EXT_POT )
   {
//    set the default function pointer of the built-in external potential
//    --> necessary only if it has not been overwritten by a test problem initializer
//    (1) tabular
      if ( OPT__EXT_POT == EXT_POT_TABLE  &&  Init_ExtPot_Ptr == NULL )    Init_ExtPot_Ptr = Init_ExtPot_Tabular;

//    (2) GREP
      if ( OPT__EXT_POT == EXT_POT_GREP  &&  Init_ExtPot_Ptr == NULL )     Init_ExtPot_Ptr = Init_ExtPot_GREP;

      if ( Init_ExtPot_Ptr != NULL )   Init_ExtPot_Ptr();
      else                             Aux_Error( ERROR_INFO, "Init_ExtPot_Ptr == NULL !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ExtAccPot



#endif // #ifdef GRAVITY
