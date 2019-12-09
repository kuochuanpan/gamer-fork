#include "GAMER.h"
#include "TestProb.h"
#include "NuclearEos.h"


// problem-specific global variables
// =======================================================================================
int    GREP_Center_Method;    // center of radial profile                                    [1]
                              // (0=max density, 1=box center, 3={ 0.0 })
int    GREP_MaxIter;          // number of iteration for constructing mass_TOV and Gamma_TOV [1000]
bool   GREP_LogBin;           // scale of bin in radila profile                              [true]
                              // (true=logarithmic, false=linear)
double GREP_LogBinRatio;      // ratio in logarithmic bin                                    [1.25]
double GREP_MaxRadius;        // maximum radius in the radial profile                        [-1.0]
                              // (<0=spearation between vertex farthest from the center)
double GREP_MinBinSize;       // minimum bin size                                            [-1.0]
                              // (<0=use amr->dh[MAX_LEVEL])

double GREP_Bfield_Ab;        // magnetic field strength                                     [1e15]
double GREP_Bfield_np;        // dependence on the density                                   [0.0]
double GREP_Bfield_Pcut;      // times of ambient pressure                                   [1.455e25]
// =======================================================================================


// for initial condition
static double *NeutronStar_Prof = NULL;        // radial progenitor model
static char   NeutronStar_ICFile[MAX_STRING];  // Filename for initial condition
static int    NeutronStar_NBin;                // number of radial bins in the progenitor model

void LoadICTable();
void Record_MigrationTest();


//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// examples
/*
// errors
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC  ||  OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC for this test !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test !!\n" );
#     endif

      if ( FLAG_BUFFER_SIZE < 5 )
         Aux_Message( stderr, "WARNING : it's recommended to set FLAG_BUFFER_SIZE >= 5 for this test !!\n" );
   } // if ( MPI_Rank == 0 )
*/


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
//   ReadPara->Add( "var_bool",          &var_bool,              true,          Useless_bool,     Useless_bool      );
//   ReadPara->Add( "var_double",        &var_double,            1.0,           Eps_double,       NoMax_double      );
//   ReadPara->Add( "var_int",           &var_int,               2,             0,                5                 );
//   ReadPara->Add( "var_str",            var_str,               Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "GREP_Center_Method",  &GREP_Center_Method,    1,             0,                3                 );
   ReadPara->Add( "GREP_MaxIter",        &GREP_MaxIter,          1000,          100,              NoMax_int         );
   ReadPara->Add( "GREP_LogBin",         &GREP_LogBin,           true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "GREP_LogBinRatio",    &GREP_LogBinRatio,      1.25,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "GREP_MaxRadius",      &GREP_MaxRadius,        -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "GREP_MinBinSize",     &GREP_MinBinSize,       -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "NeutronStar_ICFile",  NeutronStar_ICFile,     Useless_str,   Useless_str,      Useless_str       );
#  ifdef MHD
   ReadPara->Add( "GREP_Bfield_Ab",      &GREP_Bfield_Ab,        1.0e15,        0.0,              NoMax_double      );
   ReadPara->Add( "GREP_Bfield_np",      &GREP_Bfield_np,        0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "GREP_Bfield_Pcut",    &GREP_Bfield_Pcut,      1.455e25,      0.0,              NoMax_double      );
#  endif

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",      TESTPROB_ID );
      Aux_Message( stdout, "  GREP_Center_Method        = %d\n",      GREP_Center_Method );
      Aux_Message( stdout, "  GREP_MaxIter              = %d\n",      GREP_MaxIter );
      Aux_Message( stdout, "  GREP_LogBin               = %d\n",      GREP_LogBin );
      Aux_Message( stdout, "  GREP_LogBinRatio          = %13.7e\n",  GREP_LogBinRatio );
      Aux_Message( stdout, "  GREP_MaxRadius            = %13.7e\n",  GREP_MaxRadius );
      Aux_Message( stdout, "  GREP_MinBinSize           = %13.7e\n",  GREP_MinBinSize );
      Aux_Message( stdout, "  NeutronStar_ICFile        = %s\n",      NeutronStar_ICFile );
#     ifdef MHD
      Aux_Message( stdout, "  GREP_Bfield_Ab            = %13.7e\n",  GREP_Bfield_Ab );
      Aux_Message( stdout, "  GREP_Bfield_np            = %13.7e\n",  GREP_Bfield_np );
      Aux_Message( stdout, "  GREP_Bfield_Pcut          = %13.7e\n",  GREP_Bfield_Pcut );
#     endif
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double  BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double *Table_R      = NeutronStar_Prof + 0*NeutronStar_NBin;
   const double *Table_Velr   = NeutronStar_Prof + 1*NeutronStar_NBin;
   const double *Table_Dens   = NeutronStar_Prof + 2*NeutronStar_NBin;
   const double *Table_Pres   = NeutronStar_Prof + 3*NeutronStar_NBin;

   double dens, velr, pres;

   const double x0 = x - BoxCenter[0];
   const double y0 = y - BoxCenter[1];
   const double z0 = z - BoxCenter[2];
   const double r = SQRT( SQR( x0 ) + SQR( y0 ) + SQR( z0 ) );

   dens = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r);
   velr = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Velr, r);
   pres = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r);

   fluid[DENS] = dens;
   fluid[MOMX] = dens*velr*x0/r;
   fluid[MOMY] = dens*velr*y0/r;
   fluid[MOMZ] = dens*velr*z0/r;
   fluid[ENGY] = pres / ( GAMMA - 1.0 )
               + 0.5*( SQR( fluid[MOMX] ) + SQR( fluid[MOMY] ) + SQR( fluid[MOMZ] ) ) / dens;

} // FUNCTION : SetGridIC


#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                2. Use vector field for generating poloidal B field, defined in Liu+ 2008, Phys. Rev. D78, 024012
//                     A_phi = Ab * \bar\omega^2 * (1 - rho / rho_max)^np * max(P - Pcut, 0)
//                   where
//                     \omega^2 = (x - x_center)^2 + y^2
//                   And
//                     A_x = -(y / \bar\omega^2) * A_phi;  A_y = (x / \bar\omega^2) * A_phi;  A_z = 0
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

   const double  BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double *Table_R      = NeutronStar_Prof + 0*NeutronStar_NBin;
   const double *Table_Dens   = NeutronStar_Prof + 2*NeutronStar_NBin;
   const double *Table_Pres   = NeutronStar_Prof + 3*NeutronStar_NBin;

   const double x0 = x - BoxCenter[0];
   const double y0 = y - BoxCenter[1];
   const double z0 = z - BoxCenter[2];

   const double dens_c = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, 0.0);
   const double Pcut   = GREP_Bfield_Pcut / UNIT_P;
   const double Ab     = GREP_Bfield_Ab / UNIT_B;

// Use finite difference to compute the B field
   double diff = amr->dh[TOP_LEVEL];
   double r_xm, r_xp, dens_xm, dens_xp, pres_xm, pres_xp;
   double r_ym, r_yp, dens_ym, dens_yp, pres_ym, pres_yp;
   double r_zm, r_zp, dens_zm, dens_zp, pres_zm, pres_zp;

   r_xm    = SQRT( SQR( y0 ) + SQR( z0 ) + SQR( x0 - diff ) );
   r_xp    = SQRT( SQR( y0 ) + SQR( z0 ) + SQR( x0 + diff ) );
   r_ym    = SQRT( SQR( z0 ) + SQR( x0 ) + SQR( y0 - diff ) );
   r_yp    = SQRT( SQR( z0 ) + SQR( x0 ) + SQR( y0 + diff ) );
   r_zm    = SQRT( SQR( x0 ) + SQR( y0 ) + SQR( z0 - diff ) );
   r_zp    = SQRT( SQR( x0 ) + SQR( y0 ) + SQR( z0 + diff ) );

   dens_xm = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r_xm);
   dens_xp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r_xp);
   dens_ym = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r_ym);
   dens_yp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r_yp);
   dens_zm = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r_zm);
   dens_zp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r_zp);

   pres_xm = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r_xm);
   pres_xp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r_xp);
   pres_ym = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r_ym);
   pres_yp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r_yp);
   pres_zm = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r_zm);
   pres_zp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r_zp);

/*
   double dAy_dx = ( ( x0 + diff ) * POW( 1.0 - dens_xp / dens_c, GREP_Bfield_np ) * MAX( pres_xp - Pcut, 0.0 )   \
                 -   ( x0 - diff ) * POW( 1.0 - dens_xm / dens_c, GREP_Bfield_np ) * MAX( pres_xm - Pcut, 0.0 ) ) \
                 / ( 2.0  * diff );

   double dAx_dy = ( -( y0 + diff ) * POW( 1.0 - dens_yp / dens_c, GREP_Bfield_np ) * MAX( pres_yp - Pcut, 0.0 )   \
                 -   -( y0 - diff ) * POW( 1.0 - dens_ym / dens_c, GREP_Bfield_np ) * MAX( pres_ym - Pcut, 0.0 ) ) \
                 / ( 2.0  * diff );

   double dAphi_dz = ( POW( 1.0 - dens_zp / dens_c, GREP_Bfield_np ) * MAX( pres_zp - Pcut, 0.0 )   \
                   -   POW( 1.0 - dens_zm / dens_c, GREP_Bfield_np ) * MAX( pres_zm - Pcut, 0.0 ) ) \
                   / ( 2.0  * diff );
*/
   double dAy_dx = ( ( x0 + diff ) * POW( 1.0 - dens_xp / dens_c, GREP_Bfield_np ) * ( (pres_xp - Pcut > 0.0)? 1.0: 0.0 )  \
                 -   ( x0 - diff ) * POW( 1.0 - dens_xm / dens_c, GREP_Bfield_np ) * ( (pres_xm - Pcut > 0.0)? 1.0: 0.0 ) ) \
                 / ( 2.0  * diff );

   double dAx_dy = ( -( y0 + diff ) * POW( 1.0 - dens_yp / dens_c, GREP_Bfield_np ) * ( (pres_yp - Pcut > 0.0)? 1.0: 0.0 )   \
                 -   -( y0 - diff ) * POW( 1.0 - dens_ym / dens_c, GREP_Bfield_np ) * ( (pres_ym - Pcut > 0.0)? 1.0: 0.0 ) ) \
                 / ( 2.0  * diff );

   double dAphi_dz = ( POW( 1.0 - dens_zp / dens_c, GREP_Bfield_np ) * ( (pres_zp - Pcut > 0.0)? 1.0: 0.0 )   \
                   -   POW( 1.0 - dens_zm / dens_c, GREP_Bfield_np ) * ( (pres_zm - Pcut > 0.0)? 1.0: 0.0 ) ) \
                   / ( 2.0  * diff );


   magnetic[MAGX] = -x0 * Ab * dAphi_dz;
   magnetic[MAGY] = -y0 * Ab * dAphi_dz;
   magnetic[MAGZ] =       Ab * ( dAy_dx - dAx_dy );

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_NeutronStar_MigrationTest
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_NeutronStar_MigrationTest()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();

// Load IC Table
   if ( OPT__INIT != INIT_BY_RESTART )   LoadICTable();

// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr         = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr  = SetBFieldIC;
#  endif
   Init_Field_User_Ptr            = NULL; // set NCOMP_PASSIVE_USER;          example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
   Flag_User_Ptr                  = NULL; // option: OPT__FLAG_USER;          example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr       = NULL; // option: OPT__DT_USER;            example: Miscellaneous/Mis_GetTimeStep_User.cpp
   BC_User_Ptr                    = NULL; // option: OPT__BC_FLU_*=4;         example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr       = NULL; // option: OPT__RESET_FLUID;        example: Fluid/Flu_ResetByUser.cpp
   Output_User_Ptr                = NULL; // option: OPT__OUTPUT_USER;        example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr            = Record_MigrationTest; // option: OPT__RECORD_USER;        example: Auxiliary/Aux_Record_User.cpp
   Init_User_Ptr                  = NULL; // option: none;                    example: none
   End_User_Ptr                   = NULL; // option: none;                    example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
   Src_User_Ptr                   = NULL; // option: SRC_USER
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr           = NULL; // option: OPT__GRAVITY_TYPE=2/3;   example: SelfGravity/Init_ExternalAcc.cpp
   Init_ExternalPot_Ptr           = NULL; // option: OPT__EXTERNAL_POT;       example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> Init_ExtPot()
   Poi_AddExtraMassForGravity_Ptr = NULL; // option: OPT__GRAVITY_EXTRA_MASS; example: none
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr        = NULL; // option: PAR_INIT=1;              example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr    = NULL; // set PAR_NATT_USER;               example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_NeutronStar_MigrationTest



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadICTable
// Description :  Load inpu table file for initial condition
//-------------------------------------------------------------------------------------------------------
void LoadICTable()
{

   const bool RowMajor_No  = false;           // load data into the column major OPT__RECORD_USER
   const bool AllocMem_Yes = true;            // allocate memort for NeutronStar_Prof
   const int  NCol         = 4;               // total number of columns to load
   const int  TargetCols[NCol] = { 0,1,2,3 }; // target columns: {radius, vr, density, pressure}

   double *Table_R, *Table_Dens, *Table_Pres, *Table_Velr;

   NeutronStar_NBin = Aux_LoadTable( NeutronStar_Prof, NeutronStar_ICFile, NCol, TargetCols, RowMajor_No, AllocMem_Yes );

   Table_R    = NeutronStar_Prof + 0*NeutronStar_NBin;
   Table_Velr = NeutronStar_Prof + 1*NeutronStar_NBin;
   Table_Dens = NeutronStar_Prof + 2*NeutronStar_NBin;
   Table_Pres = NeutronStar_Prof + 3*NeutronStar_NBin;

   // convert to code units (assuming progentior model is in cgs)
   for (int b=0; b<NeutronStar_NBin; b++)
   {
      Table_R[b]    /= UNIT_L;
      Table_Velr[b] /= UNIT_V;
      Table_Dens[b] /= UNIT_D;
      Table_Pres[b] /= UNIT_P;
   }

} // FUNCTION : LoadICTable()



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_MigrationTest
// Description :  Record the central density
//-------------------------------------------------------------------------------------------------------
void Record_MigrationTest()
{

   const char   filename_central_dens[] = "Record__CentralDens";
   const double BoxCenter[3]            = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
// make sure that the central region is always resolved by the finest level
   const double r_max2                  = SQR( amr->dh[NLEVEL - 1] );

// allocate memory for per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   double DataCoord[4] = { -__DBL_MAX__ }, **OMP_DataCoord=NULL;
   Aux_AllocateArray2D( OMP_DataCoord, NT, 4 );


#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize arrays
      OMP_DataCoord[TID][0] = -__DBL_MAX__;
      for (int b=1; b<4; b++)   OMP_DataCoord[TID][b] = 0.0;

      for (int lv=0; lv<NLEVEL; lv++)
      {
         const double dh = amr->dh[lv];

#        pragma omp for schedule( runtime )
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            if ( amr->patch[0][lv][PID]->son != -1 )  continue;

            for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh;
            for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh;
            for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh;

               const double dx = x - BoxCenter[0];
               const double dy = y - BoxCenter[1];
               const double dz = z - BoxCenter[2];
               const double r2 = SQR(dx) + SQR(dy) + SQR(dz);

               if ( r2 < r_max2 )
               {
                  const double dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i] * UNIT_D;

                  if ( dens > OMP_DataCoord[TID][0] )
                  {
                     OMP_DataCoord[TID][0] = dens;
                     OMP_DataCoord[TID][1] = x;
                     OMP_DataCoord[TID][2] = y;
                     OMP_DataCoord[TID][3] = z;
                  }

               } // if ( r2 < r_max2 )
            }}} // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // for (int lv=0; lv<NLEVEL; lv++)
   } // OpenMP parallel region


// find the maximum over all OpenMP threads
   for (int TID=0; TID<NT; TID++)
   {
      if ( OMP_DataCoord[TID][0] > DataCoord[0] )
      {
         for (int b=0; b<4; b++)   DataCoord[b] = OMP_DataCoord[TID][b];
      }
   }

// free per-thread arrays
   Aux_DeallocateArray2D( OMP_DataCoord );


// collect data from all ranks
# ifndef SERIAL
   {
      double DataCoord_All[4 * MPI_NRank] = { 0.0 };

      MPI_Allgather( &DataCoord, 4, MPI_DOUBLE, &DataCoord_All, 4, MPI_DOUBLE, MPI_COMM_WORLD );

      for (int i=0; i<MPI_NRank; i++)
      {
         if ( DataCoord_All[4 * i] > DataCoord[0] )
            for (int b=0; b<4; b++)   DataCoord[b] = DataCoord_All[4 * i + b];
      }
   }
# endif // ifndef SERIAL


// output to file
   if ( MPI_Rank == 0 )
   {

        static bool FirstTime = true;

        if ( FirstTime )
        {
//          write header before the first output
            if ( Aux_CheckFileExist(filename_central_dens) )
            {
                Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_central_dens );
            }
            else
            {
                FILE *file_max_dens = fopen( filename_central_dens, "w" );
                fprintf( file_max_dens, "#%19s   %10s   %14s   %14s   %14s   %14s\n", "Time", "Step", "Dens", "Posx", "Posy", "Posz" );
                fclose( file_max_dens );
            }

            FirstTime = false;
        }

     FILE *file_max_dens = fopen( filename_central_dens, "a" );
     fprintf( file_max_dens,
              "%20.14e   %10ld   %14.7e   %14.7e   %14.7e   %14.7e\n",
              Time[0], Step, DataCoord[0], DataCoord[1], DataCoord[2], DataCoord[3] );
     fclose( file_max_dens );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Record_MigrationTest()

