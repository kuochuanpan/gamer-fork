#include "NuclearEoS.h"

#if ( MODEL == HYDRO  &&  EOS == EOS_NUCLEAR )



#ifdef __CUDACC__

#include "linterp_some.cu"

GPU_DEVICE static
void nuc_eos_C_cubinterp_some( double x, double y, double z,
                               double *output_vars, const double *alltables,
                               int nx, int ny, int nz, int nvars,
                               const double *xt, const double *yt, const double *zt );
#else

void nuc_eos_C_linterp_some( double x, double y, double z,
                             double *output_vars, const double *alltables,
                             int nx, int ny, int nz, int nvars,
                             const double *xt, const double *yt, const double *zt ) ;

#endif // #ifdef __CUDACC__ ... else ...




//-------------------------------------------------------------------------------------
// Function    :  nuc_eos_C_cubinterp_some
// Description :  Find thermodynamic variables using 3D Catmull-Rom cubic interpolation
//                formula by searching the tabulated nuclear EoS
//
// Note        :  1. Invoked by nuc_eos_C_short()
//
// Parameter   :  x           : Input vector of first  variable (rho)
//                y           : Input vector of second variable (eps)
//                z           : Input vector of third  variable (Y_e)
//                output_vars : Output variables of interpolated function values
//                alltables   : 3D array of tabulated variables
//                nx          : X-dimension of table
//                ny          : Y-dimension of table
//                nz          : Z-dimension of table
//                nvars       : Number of variables we will find
//                xt          : Vector of x-coordinates of table
//                yt          : Vector of y-coordinates of table
//                zt          : Vector of z-coordinates of table
//
// Return      :  output_vars
//-------------------------------------------------------------------------------------
GPU_DEVICE
void nuc_eos_C_cubinterp_some( double x, double y, double z,
                               double *output_vars, const double *alltables,
                               int nx, int ny, int nz, int nvars,
                               const double *xt, const double *yt, const double *zt )
{

   const double *pv = NULL;

   double dx, dy, dz, dxi, dyi, dzi;
   double delx, dely, delz;
   double u[4], v[4], w[4];
   double r[4], q[4];
   double vox;
   int    ix, iy, iz;
   int    nxy;

   nxy = nx*ny;

// determine spacing parameters of equidistant (!!!) table
#  if 1
   dx = ( xt[nx-1] - xt[0] ) / ( 1.0*(nx-1) );
   dy = ( yt[ny-1] - yt[0] ) / ( 1.0*(ny-1) );
   dz = ( zt[nz-1] - zt[0] ) / ( 1.0*(nz-1) );

   dxi = 1.0/dx;
   dyi = 1.0/dy;
   dzi = 1.0/dz;
#  endif

#  if 0
   dx = drho;
   dy = deps;
   dz = dye;

   dxi = drhoi;
   dyi = depsi;
   dzi = dyei;
#  endif


// determine location in table
   ix = (int)( ( x - xt[0] + 1.0e-10 ) * dxi );
   iy = (int)( ( y - yt[0] + 1.0e-10 ) * dyi );
   iz = (int)( ( z - zt[0] + 1.0e-10 ) * dzi );

   if ( ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz )
      return;

// linear interpolation at boundaries
   if ( ix == 0 || iy == 0 || iz == 0 ||
        ix == nx-2 || iy == ny-2 || iz == nz-2 )
   {
      nuc_eos_C_linterp_some( x, y, z, output_vars, alltables,
                              nx, ny, nz, nvars, xt, yt, zt );
      return;
   }


// difference
   delx = ( x - xt[ix] ) * dxi;
   dely = ( y - yt[iy] ) * dyi;
   delz = ( z - zt[iz] ) * dzi;

// factors for Catmull-Rom interpolation
   u[0] = -0.5*CUBE(delx) +     SQR(delx) - 0.5*delx;
   u[1] =  1.5*CUBE(delx) - 2.5*SQR(delx) + 1.0;
   u[2] = -1.5*CUBE(delx) + 2.0*SQR(delx) + 0.5*delx;
   u[3] =  0.5*CUBE(delx) - 0.5*SQR(delx);

   v[0] = -0.5*CUBE(dely) +     SQR(dely) - 0.5*dely;
   v[1] =  1.5*CUBE(dely) - 2.5*SQR(dely) + 1.0;
   v[2] = -1.5*CUBE(dely) + 2.0*SQR(dely) + 0.5*dely;
   v[3] =  0.5*CUBE(dely) - 0.5*SQR(dely);

   w[0] = -0.5*CUBE(delz) +     SQR(delz) - 0.5*delz;
   w[1] =  1.5*CUBE(delz) - 2.5*SQR(delz) + 1.0;
   w[2] = -1.5*CUBE(delz) + 2.0*SQR(delz) + 0.5*delz;
   w[3] =  0.5*CUBE(delz) - 0.5*SQR(delz);

   for (int iv=0; iv<nvars; iv++)
   {
      vox = 0.0;
      pv = alltables + iv + NUC_TABLE_NVAR*((ix-1) + (iy-1)*nx + (iz-1)*nxy);

      for (int k=0; k<4; k++)
      {
         q[k] = 0.0;

         for (int j=0; j<4; j++)
         {
            r[j] = 0.0;

            for (int i=0; i<4; i++)
            {
               r[j] += u[i]* *pv;
               pv   += NUC_TABLE_NVAR;
            }

            q[k] += v[j]*r[j];
            pv   += NUC_TABLE_NVAR*nx - 4*NUC_TABLE_NVAR;
         }

         vox += w[k]*q[k];
         pv  += nxy*NUC_TABLE_NVAR - 4*NUC_TABLE_NVAR*nx;
      } // for (int k=0; k<4; k++)

      output_vars[iv] = vox;
   } // for (int iv=0; iv<nvars; iv++)

// linear interpolation at boundaries
   if ( isnan(*output_vars) )
   {
      nuc_eos_C_linterp_some( x, y, z, output_vars, alltables,
                              nx, ny, nz, nvars, xt, yt, zt );
      return;
   }


   return;

} // FUNCTION : nuc_eos_C_cubinterp_some



#endif // #if ( MODEL == HYDRO  &&  EOS == EOS_NUCLEAR )
