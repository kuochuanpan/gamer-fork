#include "GAMER.h"

// function pointer for the user-defined source terms
extern void (*Src_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                             const int lv, double AuxArray[], const double dt);


//# ifdef DELEPTIONIZATION
extern void (*Src_Deleptonization_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                             const int lv, double AuxArray[], const double dt, const double EngyB );
//# endif

//# if ( NEUTRINO_SCHEME == LIGHTBULB )
extern void (*Src_LightBulb_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                             const int lv, double AuxArray[], const double dt, const double EngyB );
//# endif

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_AdvanceDt
// Description :  Add various local source terms
//
// Note        :  1. Invoked by Evolve()
//                2. Grackle library is treated separately
//
// Parameter   :  lv      : Target refinement level
//                TimeNew : Target physical time to reach
//                TimeOld : Physical time before update
//                          --> This function updates physical time from TimeOld to TimeNew
//                dt      : Time interval to advance solution
//                FluSg   : Target fluid sandglass (for both input and output)
//
// Return      : fluid[] in all patches
//-------------------------------------------------------------------------------------------------------
void Src_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt, const int FluSg )
{

// skip if there is no source term (check all source terms below)
#  if ( NEUTRINO_SCHEME == LIGHTBULB )
   const bool Src_Lightbulb = true;
#  else
   const bool Src_Lightbulb = false;
#  endif

   if ( !Src_Lightbulb  &&  !SRC_DELEPTONIZATION  &&  !SRC_USER )    return;


// check
   if ( SRC_USER  &&  Src_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "Src_User_Ptr == NULL !!\n" );

// TODO: add delepetonian and lightbulb


   const double dh = amr->dh[lv];
   real   fluid[NCOMP_TOTAL];
   double x, y, z, x0, y0, z0;


#  pragma omp parallel for private( fluid, x, y, z, x0, y0, z0 ) schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

      for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
      for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
      for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

//       compatible for MHD
#        ifdef MHD
         real B[3];

         MHD_GetCellCenteredBField( B,
                                    amr->patch[ amr->FluSg[lv] ][lv][PID]->magnetic[MAGX],
                                    amr->patch[ amr->FluSg[lv] ][lv][PID]->magnetic[MAGY],
                                    amr->patch[ amr->FluSg[lv] ][lv][PID]->magnetic[MAGZ],
                                    PS1, PS1, PS1, i, j, k );

         real EngyB = 0.5 * ( SQR( B[MAGX] ) + SQR( B[MAGY] ) + SQR( B[MAGZ] ) );
#        else
         real EngyB = NULL_REAL;
#        endif

//       get the input array
         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//###REVISE: should we make different source terms "commutative"?
//       add source terms

//       (1) lightbulb neutrino scheme
         if ( Src_Lightbulb )
//#        if ( NEUTRINO_SCHEME == LIGHTBULB )
            Src_LightBulb_Ptr ( fluid, x, y, z, TimeNew, lv, NULL, dt, EngyB );
//#        endif


//       (2) deleptonization
#        ifdef DELEPTIONIZATION
         if ( SRC_DELEPTONIZATION )
            Src_Deleptonization_Ptr( fluid, x, y, z, TimeNew, lv, NULL, dt, EngyB );
#        endif

//       (3) user-defined
         if ( SRC_USER )
            Src_User_Ptr       ( fluid, x, y, z, TimeNew, lv, NULL, dt );

//       store the updated results
         for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];
      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : Src_AdvanceDt
