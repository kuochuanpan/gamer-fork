#include "CUPOT.h"
#include <stdio.h>

#if ( defined GRAVITY  &&  defined GREP )


// external functions and GPU-related set-up
#ifdef __CUDACC__

// variables reside in constant memory
__constant__ double c_GREP_Data[GR_POT_NAUX_MAX];
__constant__ double c_GREP_Edge[GR_POT_NAUX_MAX];
__constant__ double c_GREP_Center[3];
__constant__ double c_r_max2;
__constant__ int    c_GREP_NBin;


//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_SetConstMem_GREffPot
// Description :  Set the constant memory used by CUPOT_CorrectEffPot()
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by CUAPI_Init_GREffPot()
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
__host__
int CUPOT_SetConstMem_GREffPot( double h_GREP_Data[], double h_GREP_Edge[], double h_GREP_Center[],
                                double h_r_max2, int h_GREP_NBin )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( c_GREP_Data,   h_GREP_Data,   GR_POT_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   if (  cudaSuccess != cudaMemcpyToSymbol( c_GREP_Edge,   h_GREP_Edge,   GR_POT_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -2;

   if (  cudaSuccess != cudaMemcpyToSymbol( c_GREP_Center, h_GREP_Center,               3*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -3;

   if (  cudaSuccess != cudaMemcpyToSymbol( c_r_max2,     &h_r_max2,                      sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -4;

   if (  cudaSuccess != cudaMemcpyToSymbol( c_GREP_NBin,  &h_GREP_NBin,                   sizeof(int),
                                            0, cudaMemcpyHostToDevice)  )
      return -5;

   return 0;

} // FUNCTION : CUPOT_SetConstMem_GREffPot

#endif // ifdef __CUDACC__


#define LinearInterp( x, xa, xb, ya, yb )   ( ( ((x) - (xa)) * (yb) + ((xb) - (x)) * (ya) ) / ((xb) - (xa)) )




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_CorrectEffPot
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
//                dh                : Cell size
//                Undo              : Add (true) or subtract (false) potential correction to the input potential
//                USG               : Flag to indicate which potential is input
//
//-------------------------------------------------------------------------------------------------------
__global__
void CUPOT_CorrectEffPot(       real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
                                real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
                          const double g_Corner_Array [][3],
                          const real dh, const bool Undo, const bool USG )
{

// declare index for loop
#  ifdef UNSPLIT_GRAVITY
   const int IDX    = ( USG ) ? USG_NXT_G      : GRA_NXT;
   const int IDX_GZ = ( USG ) ? USG_GHOST_SIZE : GRA_GHOST_SIZE;
#  else
   const int IDX    = GRA_NXT;
   const int IDX_GZ = GRA_GHOST_SIZE;
#  endif

   const int IDX_sqr = SQR (IDX);

   const int P = blockIdx.x;
   {
//    loop over all cells of the target patch
//    _g0: indices for the arrays without any ghost zone
      CGPU_LOOP( idx_g0, CUBE(IDX) )
      {

         const int i_g0 = idx_g0 % IDX;
         const int j_g0 = idx_g0 % IDX_sqr / IDX;
         const int k_g0 = idx_g0 / IDX_sqr;

         const double dx = g_Corner_Array[P][0] + (double)((i_g0-IDX_GZ)*dh) - c_GREP_Center[0];
         const double dy = g_Corner_Array[P][1] + (double)((j_g0-IDX_GZ)*dh) - c_GREP_Center[1];
         const double dz = g_Corner_Array[P][2] + (double)((k_g0-IDX_GZ)*dh) - c_GREP_Center[2];

         const double r2 = SQR(dx) + SQR(dy) + SQR(dz);


         if ( r2 < c_r_max2 )
         {
            const double r = SQRT( r2 );

//          use binary search algorithm to find the index of bin
            int bin;
            for ( int i=0, j=c_GREP_NBin; j-i != 1; bin = (i+j)/2 )
            {
               int mid = (i+j)/2;
               if ( r > c_GREP_Edge[mid] )   i = mid;
               else                          j = mid;
            }

            double phi = ( bin == c_GREP_NBin-1 ) ? c_GREP_Data[bin]
                                                  : LinearInterp( r, c_GREP_Edge[bin], c_GREP_Edge[bin+1],
                                                                     c_GREP_Data[bin], c_GREP_Data[bin+1] );

            if ( Undo )   phi = -phi;

#           ifdef UNSPLIT_GRAVITY
            if ( USG )
               g_Pot_Array_USG[P][idx_g0] += (real)phi;
            else
               g_Pot_Array_New[P][idx_g0] += (real)phi;
#           else
               g_Pot_Array_New[P][idx_g0] += (real)phi;
#           endif
         } // if ( r2 < r_max2 )

      } // CGPU_LOOP( idx_g0, CUBE(PS1) )
   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU_CorrectEffPot


#endif // #if ( defined GRAVITY  &&  defined GREP )
