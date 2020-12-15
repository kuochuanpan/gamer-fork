#include "GAMER.h"

#ifdef GREP


extern Profile_t *DensAve [NLEVEL+1][2];
extern Profile_t *EngyAve [NLEVEL+1][2];
extern Profile_t *VrAve   [NLEVEL+1][2];
extern Profile_t *PresAve [NLEVEL+1][2];
extern Profile_t *Phi_eff [NLEVEL  ][2];




//-------------------------------------------------------------------------------------------------------
// Function    :  End_GREP
// Description :  Free memory previously allocated by Init_GREP()
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_GREP()
{

   for (int Sg=0; Sg<2; Sg++)
   for (int lv=0; lv<=NLEVEL; lv++)
   {
      DensAve [lv][Sg]->FreeMemory();
      EngyAve [lv][Sg]->FreeMemory();
      VrAve   [lv][Sg]->FreeMemory();
      PresAve [lv][Sg]->FreeMemory();

      if ( lv < NLEVEL )
      Phi_eff [lv][Sg]->FreeMemory();
   }

} // FUNCTION : End_GREP



#endif // #ifdef GREP
