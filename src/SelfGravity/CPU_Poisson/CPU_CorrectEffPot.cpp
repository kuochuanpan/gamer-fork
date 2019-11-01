#include "GAMER.h"
//#include "NuclearEos.h"


#if ( defined GRAVITY  &&  defined GREP )

#define LinearInterp( x, xa, xb, ya, yb )   ( ( ((x) - (xa)) * (yb) + ((xb) - (x)) * (ya) ) / ((xb) - (xa)) )

extern Profile_t Phi_eff;


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_CorrectEffPot
// Description :  Do/Undo the GR potential correction
//
// Note        :  1. Support CPU only in current version
//                2. The potential correction calculated at the current step is applied to
//                   both g_Pot_Array_New and g_Pot_Array_USG in current version
//
// Parameter   :  g_Pot_Array_New   : Array storing the input potential (at the current step)
//                                    --> _New: to be distinguishable from g_Pot_Array_USG[], which is defined at the previous step
//                g_Pot_Array_USG   : Array storing the input potential for UNSPLIT_GRAVITY (at the previous step)
//                g_Corner_Array    : Array storing the physical corner coordinates of each patch
//                NPatchGroup       : Number of input patch groups (for CPU only)
//                dh                : Cell size
//                Undo              : Add (true) or subtract (false) potential correction to the input potential
//                USG               : Flag to indicate which potential is input
//
//-------------------------------------------------------------------------------------------------------
void CPU_CorrectEffPot(       real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
                              real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
                        const double g_Corner_Array [][3],
                        const int    NPatchGroup,
                        const real dh, const bool Undo, const bool USG)
{

// Profile information
   const double *Data   = Phi_eff.Data;
   const double *Radius = Phi_eff.Radius;
         double *Center = Phi_eff.Center;
         double  r_max2 = SQR( Phi_eff.MaxRadius );
         double  EdgeL[Phi_eff.NBin] = { 0.0 };

   for ( int i=1; i<Phi_eff.NBin; i++ )   EdgeL[i] = ( Phi_eff.LogBin ) ? sqrt( Radius[i - 1] * Radius[i] )
                                                                        : 0.5*( Radius[i - 1] + Radius[i] );

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
      real *pot_new = ( USG ) ? g_Pot_Array_USG[P] : g_Pot_Array_New[P];
#     else
      real *pot_new = g_Pot_Array_New[P];
#     endif

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

//          use binary search algorithm to find the index of bin
            int bin;
            for ( int i=0, j=Phi_eff.NBin - 1; j - i != 1; bin = (i + j) / 2 )
            {
               int mid = (i + j) / 2;
               if ( r > EdgeL[mid] )   i = mid;
               else                    j = mid;
            }

//          prevent from round-off errors
            if ( bin >= Phi_eff.NBin )   continue;

//          check
#           ifdef GAMER_DEBUG
            if ( bin < 0 )    Aux_Error( ERROR_INFO, "bin (%d) < 0 !!\n", bin );
#           endif

            double phi = ( bin == Phi_eff.NBin-1 ) ? Data[bin]
                                                   : LinearInterp( r, EdgeL[bin], EdgeL[bin+1], Data[bin], Data[bin+1] );

            if ( (r < EdgeL[bin])  || ( r > EdgeL[bin + 1]) )
            printf("r = %.6e\tEdgeL = %.6e\tEdgeR = %.6e\n", r, EdgeL[bin], EdgeL[bin+1]);

            pot_new[idx_g0] += ( Undo ) ? -(real)phi : (real)phi;
         } // if ( r2 < r_max2 )
      } // for (int idx_g0=0; idx_g0<CUBE(IDX); idx_g0++)
   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU_CorrectEffPot


#endif // #if ( defined GRAVITY  &&  defined GREP )
