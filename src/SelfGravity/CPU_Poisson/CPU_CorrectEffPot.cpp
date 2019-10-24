#include "GAMER.h"
//#include "NuclearEos.h"


#if ( defined GRAVITY  &&  defined GREP )

#define LinearInterp( x, xa, xb, ya, yb )   ( ( ((x) - (xa)) * (yb) + ((xb) - (x)) * (ya) ) / ((xb) - (xa)) )

extern Profile_t Phi_eff;


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_CorrectEffPot
// Description :  Apply/Undo the GR potential correction
//
// Note        :  1. CPU only
//
// Parameter   :
//
// TODO        :  1. clean up and add document
//                2. Combine the two loops in lv = 0 case into one
//                3. divide routine for lv > 0 into another file, and add GPU version
//-------------------------------------------------------------------------------------------------------
void CPU_CorrectEffPot(       real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
                              real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
                        const double g_Corner_Array [][3],
                        const int    NPatchGroup,
                        const real dh, const bool Undo, const bool USG)
{

// Profile information
   const double *Data                = Phi_eff.Data;
         double *Center              = Phi_eff.Center;
         double  r_max2              = SQR( Phi_eff.MaxRadius );
         double  dr_min              = ( Phi_eff.LogBin ) ? Phi_eff.Radius[0] / pow( Phi_eff.LogBinRatio, -0.5 )
                                                          : Phi_eff.Radius[0] * 2.0;
         double  EdgeL[Phi_eff.NBin] = { 0.0 };

   for ( int i=0; i<Phi_eff.NBin; i++ )
      EdgeL[i] = ( Phi_eff.LogBin ) ? dr_min*pow( Phi_eff.LogBinRatio, (real)(i - 1) )
                                    : (real)i*dr_min;

// declare index for loop
#  ifdef UNSPLIT_GRAVITY
   const int IDX    = ( USG ) ? USG_NXT_G      : GRA_NXT;
   const int IDX_GZ = ( USG ) ? USG_GHOST_SIZE : GRA_GHOST_SIZE;
#  else
   const int IDX    = GRA_NXT;
   const int IDX_GZ = GRA_GHOST_SIZE;
#  endif

   const int IDX_sqr = SQR (IDX);


#  pragma omp parallel for schedule( runtime )
   for (int P=0; P<NPatchGroup*8; P++)
   {
#     ifdef UNSPLIT_GRAVITY
      real *pot_new = ( USG ) ? g_Pot_Array_USG[P]      : g_Pot_Array_New[P];
#     else
      real *pot_new = g_Pot_Array_New[P];
#     endif

      double phi;

//    correct potential (including ghost zone)
      for (int idx_g0=0; idx_g0<CUBE(IDX); idx_g0++)
      {
         const int i_g0    = idx_g0 % IDX;
         const int j_g0    = idx_g0 % IDX_sqr / IDX;
         const int k_g0    = idx_g0 / IDX_sqr;

         const double dx = g_Corner_Array[P][0] + (double)((i_g0-IDX_GZ)*dh) - Center[0];
         const double dy = g_Corner_Array[P][1] + (double)((j_g0-IDX_GZ)*dh) - Center[1];
         const double dz = g_Corner_Array[P][2] + (double)((k_g0-IDX_GZ)*dh) - Center[2];

         const double r2 = SQR(dx) + SQR(dy) + SQR(dz);

         if ( r2 < r_max2 )
         {
            const double r   = SQRT( r2 );
            const int    bin = ( Phi_eff.LogBin ) ? (  (r<dr_min) ? 0 : int( log(r/dr_min)/log(Phi_eff.LogBinRatio) ) + 1  )
                                                  : int( r/dr_min );
//          prevent from round-off errors
            if ( bin >= Phi_eff.NBin )   continue;

//          check
#           ifdef GAMER_DEBUG
            if ( bin < 0 )    Aux_Error( ERROR_INFO, "bin (%d) < 0 !!\n", bin );
#           endif

            phi = ( bin == Phi_eff.NBin-1 ) ? Data[bin]
                                            : LinearInterp( r, EdgeL[bin], EdgeL[bin+1], Data[bin], Data[bin+1] );

            pot_new[idx_g0] += ( Undo ) ? -(real)phi : (real)phi;
         } // if ( r2 < r_max2 )
      } // for (int idx_g0=0; idx_g0<CUBE(IDX); idx_g0++)
   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU_CorrectEffPot


#endif // #if ( defined GRAVITY  &&  defined GREP )
