#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Region
// Description :  Check if the element (i,j,k) of the input patch is within the regions allowed to be refined
//
// Note        :  To use this functionality, please turn on the option "OPT__FLAG_REGION" and then specify the
//                target regions in this file
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[0][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//
// Return      :  "true/false"  if the input cell "is/is not" within the region allowed for refinement
//-------------------------------------------------------------------------------------------------------
bool Flag_Region( const int i, const int j, const int k, const int lv, const int PID )
{

   const double dh     = amr->dh[lv];                                         // cell size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,     // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   bool Within = true;


// put the target region below
// ##########################################################################################################
/*
// Example : sphere
   const double Center[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double dR[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R         = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );
   const double MaxR      = 1.0;

   Within = R <= MaxR;
*/
// ##########################################################################################################

   const double Center[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double dR[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R         = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

   const double baseR     = 32.0e5/UNIT_L; // km

   double fmaxR;

   // For Core Collapse Supernova Problem
   if (TESTPROB_ID == 17) {

        if  (lv >= NLEVEL -1) {
            fmaxR = baseR;
        } else if (lv >= NLEVEL-2) {
            fmaxR = baseR*2.0 ;
        } else if (lv >= NLEVEL-3) {
            fmaxR = baseR*4.0 ;
        } else if (lv >= NLEVEL-4) {
            fmaxR = baseR*8.0 ;
        } else if (lv >= NLEVEL-5) {
            fmaxR = baseR*16.0 ;
        } else if (lv >= NLEVEL-6) {
            fmaxR = baseR*32.0 ;
        } else if (lv >= NLEVEL-7) {
            fmaxR = baseR*64.0 ;
        } else {
            fmaxR = baseR*128.0 ;
        }
        Within = R <= fmaxR;
        return Within;
   }else {
        return Within;
   }

} // FUNCTION : Flag_Region

