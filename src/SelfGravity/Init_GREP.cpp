#include "GAMER.h"

#if ( defined GRAVITY  &&  defined GREP )


extern Profile_t *DensAve    [NLEVEL+1][2];
extern Profile_t *EngyAve    [NLEVEL+1][2];
extern Profile_t *VrAve      [NLEVEL+1][2];
extern Profile_t *PresAve    [NLEVEL+1][2];
extern Profile_t *Phi_eff    [NLEVEL  ][2];

extern int    GREPSg     [NLEVEL];
extern double GREPSgTime [NLEVEL][2];
extern double GREP_Prof_Center   [3];

extern double GREP_Prof_MaxRadius;
extern double GREP_Prof_MinBinSize;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_GREP
// Description :  Initialize the GREP Profiles_t objects and parameters
//-------------------------------------------------------------------------------------------------------
void Init_GREP()
{

// (1) intialize the GREP profiles
   for (int Sg=0; Sg<2; Sg++)
   for (int lv=0; lv<=NLEVEL; lv++)
   {
      DensAve [lv][Sg] = new Profile_t();
      EngyAve [lv][Sg] = new Profile_t();
      VrAve   [lv][Sg] = new Profile_t();
      PresAve [lv][Sg] = new Profile_t();

      if ( lv < NLEVEL )
      Phi_eff [lv][Sg] = new Profile_t();
   }


// (2) initialize GREP Sg and SgTime
   for (int lv=0; lv<NLEVEL; lv++)
   {
      GREPSg[lv] = 0;
      for (int Sg=0; Sg<2; Sg++)   GREPSgTime[lv][Sg] = -__FLT_MAX__;
   }


// (3) initialize the GREP parameters
   switch ( GREP_CENTER_METHOD )
   {
      case 1:
         for (int i=0; i<3; i++)   GREP_Prof_Center[i] = amr->BoxCenter[i];
      break;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GREP_CENTER_METHOD", GREP_CENTER_METHOD );
   }

   GREP_Prof_MinBinSize = ( GREP_MINBINSIZE > 0.0 ) ? GREP_MINBINSIZE : amr->dh[MAX_LEVEL];

   GREP_Prof_MaxRadius  = ( GREP_MAXRADIUS > 0.0 )
                        ? GREP_MAXRADIUS
                        : SQRT( SQR( MAX( amr->BoxSize[0] - GREP_Prof_Center[0], GREP_Prof_Center[0] ) )
                        +       SQR( MAX( amr->BoxSize[1] - GREP_Prof_Center[1], GREP_Prof_Center[1] ) )
                        +       SQR( MAX( amr->BoxSize[2] - GREP_Prof_Center[2], GREP_Prof_Center[2] ) ) );

} // FUNCTION : Init_GREP



#endif // #if ( defined GRAVITY  &&  defined GREP )
