#include "NuclearEoS.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#include "CUFLU_Shared_FluUtility.cu"
#endif

#if ( MODEL == HYDRO  &&  EOS == EOS_NUCLEAR )



#ifdef __CUDACC__

#include "NuclearEoS.cu"

#else

// global variables
int    g_nrho;
int    g_neps;
int    g_nye;
int    g_nmode;
double g_energy_shift;

double *g_alltables      = NULL;
double *g_alltables_mode = NULL;
double *g_logrho         = NULL;
double *g_logeps         = NULL;
double *g_yes            = NULL;
double *g_logtemp_mode   = NULL;
double *g_entr_mode      = NULL;
double *g_logprss_mode   = NULL;

// prototypes
void nuc_eos_C_short( double xrho, double *xenr, double xye,
                      double *xtemp, double *xent, double *xprs,
                      double *xcs2, double *xmunu, double energy_shift,
                      int nrho, int neps, int nye, int nmode,
                      const double *alltables, const double *alltables_mode,
                      const double *logrho, const double *logeps, const double *yes,
                      const double *logtemp_mode, const double *entr_mode, const double *logprss_mode,
                      int keymode, int *keyerr, double rfeps );
void nuc_eos_C_ReadTable( char *nuceos_table_name );
void CUAPI_PassNuclearEoSTable2GPU();
#endif // #ifdef __CUDACC__ ... else ...




/********************************************************
1. Nuclear EoS (EOS_NUCLEAR)

2. This file is shared by both CPU and GPU

   GPU_EoS_Nuclear.cu -> CPU_EoS_Nuclear.cpp

3. Three steps are required to implement an EoS

   I.   Set an EoS auxiliary array
   II.  Implement EoS conversion functions
   III. Set EoS initialization functions

4. All EoS conversion functions must be thread-safe and
   not use any global variable

5. When an EoS conversion function fails, it is recommended
   to return NAN in order to trigger auto-correction such as
   "OPT__1ST_FLUX_CORR" and "AUTO_REDUCE_DT"
********************************************************/



// =============================================
// I. Set an EoS auxiliary array
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_SetAuxArray_Nuclear
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by EoS_Init_Nuclear()
//                2. AuxArray_Flt/Int[] have the size of EOS_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_Nuclear( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[NUC_AUX_ESHIFT] = g_energy_shift;

   AuxArray_Int[NUC_AUX_NRHO  ] = g_nrho;
   AuxArray_Int[NUC_AUX_NEPS  ] = g_neps;
   AuxArray_Int[NUC_AUX_NYE   ] = g_nye;
   AuxArray_Int[NUC_AUX_NMODE ] = g_nmode;

} // FUNCTION : EoS_SetAuxArray_Nuclear
#endif // #ifndef __CUDACC__



// =============================================
// II. Implement EoS conversion functions
//     (1) EoS_DensEint2Pres_*
//     (2) EoS_DensPres2Eint_*
//     (3) EoS_DensPres2CSqr_*
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_Nuclear
// Description :  Convert gas mass density and internal energy density to gas pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_Nuclear() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Pres_Nuclear( const real Dens, const real Eint, const real Passive[], const double AuxArray_Flt[],
                                       const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  if ( NCOMP_PASSIVE > 0 )
   if ( Passive == NULL )  printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
#  endif
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );

   if ( Hydro_CheckNegative(Dens) )
      printf( "ERROR : invalid input density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Dens, __FILE__, __LINE__, __FUNCTION__ );
#  endif // GAMER_DEBUG


   const double EnergyShift = AuxArray_Flt[NUC_AUX_ESHIFT];
   const int    NRho        = AuxArray_Int[NUC_AUX_NRHO  ];
   const int    NEps        = AuxArray_Int[NUC_AUX_NEPS  ];
   const int    NYe         = AuxArray_Int[NUC_AUX_NYE   ];
   const int    NMode       = AuxArray_Int[NUC_AUX_NMODE ];

   int  Mode    = NUC_MODE_ENGY;
   real sEint   = Eint / Dens;
   real Ye      = Passive[ YE - NCOMP_FLUID ] / Dens;
   real Pres    = NULL_REAL;
   real Useless = NULL_REAL;
   int  Err     = NULL_INT;

   nuc_eos_C_short( Dens, &sEint, Ye, &Useless, &Useless, &Pres, &Useless, &Useless,
                    EnergyShift, NRho, NEps, NYe, NMode,
                    Table[NUC_TAB_ALL], Table[NUC_TAB_ALL_MODE], Table[NUC_TAB_RHO], Table[NUC_TAB_EPS],
                    Table[NUC_TAB_YE], Table[NUC_TAB_TEMP_MODE], Table[NUC_TAB_ENTR_MODE], Table[NUC_TAB_PRES_MODE],
                    Mode, &Err, NULL_REAL );

   if ( Err )
   {
      Pres = NAN;

#     ifdef GAMER_DEBUG
      printf( "ERROR : EoS failed with error %d !!\n", Err );
#     endif
   }


// check
#  ifdef GAMER_DEBUG
   if ( Hydro_CheckNegative(Pres) )
   {
      printf( "ERROR : invalid output pressure (%13.7e) in %s() !!\n", Pres, __FUNCTION__ );
      printf( "        Dens=%13.7e, Eint=%13.7e\n", Dens, Eint );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "        Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Pres;

} // FUNCTION : EoS_DensEint2Pres_Nuclear



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2Eint_Nuclear
// Description :  Convert gas mass density and pressure to gas internal energy density
//
// Note        :  1. See EoS_DensEint2Pres_Nuclear()
//
// Parameter   :  Dens       : Gas mass density
//                PresIn     : Gas pressure
//                Passive    : Passive scalars
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas internal energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2Eint_Nuclear( const real Dens, const real PresIn, const real Passive[], const double AuxArray_Flt[],
                                       const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  if ( NCOMP_PASSIVE > 0 )
   if ( Passive == NULL )  printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
#  endif
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );

   if ( Hydro_CheckNegative(Dens) )
      printf( "ERROR : invalid input density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Dens, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(Pres) )
      printf( "ERROR : invalid input pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Pres, __FILE__, __LINE__, __FUNCTION__ );
#  endif // GAMER_DEBUG



   const double EnergyShift = AuxArray_Flt[NUC_AUX_ESHIFT];
   const int    NRho        = AuxArray_Int[NUC_AUX_NRHO  ];
   const int    NEps        = AuxArray_Int[NUC_AUX_NEPS  ];
   const int    NYe         = AuxArray_Int[NUC_AUX_NYE   ];
   const int    NMode       = AuxArray_Int[NUC_AUX_NMODE ];

   int  Mode    = NUC_MODE_PRES;
   real Eint    = NULL_REAL;
   real sEint   = NULL_REAL;
   real Ye      = Passive[ YE - NCOMP_FLUID ] / Dens;
   real Pres    = PresIn;
   real Useless = NULL_REAL;
   int  Err     = NULL_INT;

   nuc_eos_C_short( Dens, &sEint, Ye, &Useless, &Useless, &Pres, &Useless, &Useless,
                    EnergyShift, NRho, NEps, NYe, NMode,
                    Table[NUC_TAB_ALL], Table[NUC_TAB_ALL_MODE], Table[NUC_TAB_RHO], Table[NUC_TAB_EPS],
                    Table[NUC_TAB_YE], Table[NUC_TAB_TEMP_MODE], Table[NUC_TAB_ENTR_MODE], Table[NUC_TAB_PRES_MODE],
                    Mode, &Err, NULL_REAL );

   if ( Err )
   {
      sEint = NAN;

#     ifdef GAMER_DEBUG
      printf( "ERROR : EoS failed with error %d !!\n", Err );
#     endif
   }

   Eint = sEint*Dens;


   return Eint;

} // FUNCTION : EoS_DensPres2Eint_Nuclear



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_Nuclear
// Description :  Convert gas mass density and pressure to sound speed squared
//
// Note        :  1. See EoS_DensEint2Pres_Nuclear()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Sound speed square
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_Nuclear( const real Dens, const real PresIn, const real Passive[], const double AuxArray_Flt[],
                                       const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  if ( NCOMP_PASSIVE > 0 )
   if ( Passive == NULL )  printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
#  endif
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );

   if ( Hydro_CheckNegative(Dens) )
      printf( "ERROR : invalid input density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Dens, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(Pres) )
      printf( "ERROR : invalid input pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Pres, __FILE__, __LINE__, __FUNCTION__ );
#  endif // GAMER_DEBUG


   const double EnergyShift = AuxArray_Flt[NUC_AUX_ESHIFT];
   const int    NRho        = AuxArray_Int[NUC_AUX_NRHO  ];
   const int    NEps        = AuxArray_Int[NUC_AUX_NEPS  ];
   const int    NYe         = AuxArray_Int[NUC_AUX_NYE   ];
   const int    NMode       = AuxArray_Int[NUC_AUX_NMODE ];

   int  Mode    = NUC_MODE_PRES;
   real Ye      = Passive[ YE - NCOMP_FLUID ] / Dens;
   real Pres    = PresIn;
   real Cs2     = NULL_REAL;
   real Useless = NULL_REAL;
   int  Err     = NULL_INT;

   nuc_eos_C_short( Dens, &Useless, Ye, &Useless, &Useless, &Pres, &Cs2, &Useless,
                    EnergyShift, NRho, NEps, NYe, NMode,
                    Table[NUC_TAB_ALL], Table[NUC_TAB_ALL_MODE], Table[NUC_TAB_RHO], Table[NUC_TAB_EPS],
                    Table[NUC_TAB_YE], Table[NUC_TAB_TEMP_MODE], Table[NUC_TAB_ENTR_MODE], Table[NUC_TAB_PRES_MODE],
                    Mode, &Err, NULL_REAL );

   if ( Err )
   {
      Cs2 = NAN;

#     ifdef GAMER_DEBUG
      printf( "ERROR : EoS failed with error %d !!\n", Err );
#     endif
   }


// check
#  ifdef GAMER_DEBUG
   if ( Hydro_CheckNegative(Cs2) )
   {
      printf( "ERROR : invalid output sound speed squared (%13.7e) in %s() !!\n", Cs2, __FUNCTION__ );
      printf( "        Dens=%13.7e, Pres=%13.7e\n", Dens, Pres );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "        Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Cs2;

} // FUNCTION : EoS_DensPres2CSqr_Nuclear



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_Nuclear;
FUNC_SPACE EoS_DP2E_t EoS_DensPres2Eint_Ptr = EoS_DensPres2Eint_Nuclear;
FUNC_SPACE EoS_DP2C_t EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_Nuclear;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_SetCPU/GPUFunc_Nuclear
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_Nuclear()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_SetFunc_Nuclear( CPU_FuncPtr, GPU_FuncPtr );
//
//                3. Call-by-reference
//
// Parameter   :  EoS_DensEint2Pres_CPU/GPUPtr : CPU/GPU function pointers to be set
//                EoS_DensPres2Eint_CPU/GPUPtr : ...
//                EoS_DensPres2CSqr_CPU/GPUPtr : ...
//
// Return      :  EoS_DensEint2Pres_CPU, EoS_DensPres2Eint_CPU/GPUPtr, EoS_DensPres2CSqr_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void EoS_SetGPUFunc_Nuclear( EoS_DE2P_t &EoS_DensEint2Pres_GPUPtr,
                             EoS_DP2E_t &EoS_DensPres2Eint_GPUPtr,
                             EoS_DP2C_t &EoS_DensPres2CSqr_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Pres_GPUPtr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2Eint_GPUPtr, EoS_DensPres2Eint_Ptr, sizeof(EoS_DP2E_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_Nuclear( EoS_DE2P_t &EoS_DensEint2Pres_CPUPtr,
                             EoS_DP2E_t &EoS_DensPres2Eint_CPUPtr,
                             EoS_DP2C_t &EoS_DensPres2CSqr_CPUPtr )
{
   EoS_DensEint2Pres_CPUPtr = EoS_DensEint2Pres_Ptr;
   EoS_DensPres2Eint_CPUPtr = EoS_DensPres2Eint_Ptr;
   EoS_DensPres2CSqr_CPUPtr = EoS_DensPres2CSqr_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void EoS_SetAuxArray_Nuclear( double [], int [] );
void EoS_SetCPUFunc_Nuclear( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t & );
#ifdef GPU
void EoS_SetGPUFunc_Nuclear( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_Nuclear
// Description :  Initialize EoS
//
// Note        :  1. Set an auxiliary array by invoking EoS_SetAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU EoS routines by invoking EoS_SetCPU/GPUFunc_*()
//                3. Invoked by EoS_Init()
//                   --> Enable it by linking to the function pointer "EoS_Init_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void EoS_Init_Nuclear()
{

// check if the default maximum table size is large enough
   if ( EOS_NTABLE_MAX < NUC_TABLE_NPTR )
      Aux_Error( ERROR_INFO, "EOS_NTABLE_MAX (%d) < NUC_TABLE_NPTR (%d) for the nuclear EoS !!\n",
                 EOS_NTABLE_MAX, NUC_TABLE_NPTR );


   nuc_eos_C_ReadTable( NUC_TABLE );

   EoS_SetAuxArray_Nuclear( EoS_AuxArray_Flt, EoS_AuxArray_Int );
   EoS_SetCPUFunc_Nuclear( EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr, EoS_DensPres2CSqr_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_Nuclear( EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr, EoS_DensPres2CSqr_GPUPtr );
#  endif

#  ifdef GPU
   CUAPI_PassNuclearEoSTable2GPU();
#  endif

} // FUNCTION : EoS_Init_Nuclear

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO  &&  EOS == EOS_NUCLEAR )
