#include "GAMER.h"
//#include "NuclearEos.h"


#if ( defined GRAVITY  &&  defined GREP )

#define LinearInterp( x, xa, xb, ya, yb )   ( ( ((x) - (xa)) * (yb) + ((xb) - (x)) * (ya) ) / ((xb) - (xa)) )

extern Profile_t *Phi_eff[2];


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

// REVISE: support USG feature
      Profile_t *Phi    = Phi_eff[1];

// Profile information
   const    int  NBin   = Phi->NBin;
   const double *Data   = Phi->Data;
   const double *Radius = Phi->Radius;
   const double *Center = Phi->Center;
   const double  r_max2 = SQR( Phi->MaxRadius );
         double  Edge[NBin+1] = { 0.0 };

   for ( int i=1; i<NBin; i++ )   Edge[i] = ( Phi->LogBin ) ? sqrt( Radius[i - 1] * Radius[i] )
                                                            : 0.5*( Radius[i - 1] + Radius[i] );

   Edge[NBin] = ( Phi->LogBin ) ? SQR ( Edge[NBin - 1] ) / Edge[NBin - 2]
                                : 2.0 * Edge[NBin - 1]   - Edge[NBin - 2];


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
         const int i_g0  = idx_g0 % IDX;
         const int j_g0  = idx_g0 % IDX_sqr / IDX;
         const int k_g0  = idx_g0 / IDX_sqr;

         const double dx = g_Corner_Array[P][0] + (double)((i_g0-IDX_GZ)*dh) - Center[0];
         const double dy = g_Corner_Array[P][1] + (double)((j_g0-IDX_GZ)*dh) - Center[1];
         const double dz = g_Corner_Array[P][2] + (double)((k_g0-IDX_GZ)*dh) - Center[2];

         const double r2 = SQR(dx) + SQR(dy) + SQR(dz);

         if ( r2 < r_max2 )
         {
            const double r = SQRT( r2 );
            int bin;

//          if empty bins are removed, the separations between bins are not equal in linear/logarithmic scale
//          use binary search algorithm to find the index of bin
            for ( int i=0, j=NBin; j-i != 1; bin = (i+j)/2 )
            {
               int mid = (i+j)/2;
               if ( r > Edge[mid] )   i = mid;
               else                   j = mid;
            }

            double phi = ( bin == NBin-1 ) ? Data[bin]
                                           : LinearInterp( r, Edge[bin], Edge[bin+1], Data[bin], Data[bin+1] );

//          for debug, clean later
            if ( (r < Edge[bin])  || ( r > Edge[bin + 1]) )
            {
            printf("Bin = %d  NBin = %d\n", bin, NBin);
            printf("rmax = %.6e  \n", Phi->MaxRadius);
            printf("Incorrect r: r = %.6e\tEdgeL = %.6e\tEdgeR = %.6e  Radius[i-1] = %.6e  Radius[i] = %.6e\n", r, Edge[bin], Edge[bin+1], Radius[bin-1], Radius[bin]);
            }

            pot_new[idx_g0] += ( Undo ) ? -(real)phi : (real)phi;
         } // if ( r2 < r_max2 )
      } // for (int idx_g0=0; idx_g0<CUBE(IDX); idx_g0++)
   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU_CorrectEffPot


#endif // #if ( defined GRAVITY  &&  defined GREP )
