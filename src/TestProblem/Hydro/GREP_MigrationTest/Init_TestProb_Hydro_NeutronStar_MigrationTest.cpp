#include "GAMER.h"
#include "TestProb.h"


// problem-specific global variables
// =======================================================================================
// Parameters for the toroidal B field
#ifdef MHD
static double  Bfield_Ab;                       // magnetic field strength                            [1e15]
static double  Bfield_np;                       // dependence on the density                          [0.0]
#endif

// Parameters for GW emissions output
static int     GW_OUTPUT_OPT;                   // output GW signal (0=off)                           [0]
static int     GW_OUTPUT_VERBOSE;               // output spherically summed GW signal (0=off)        [0]
static double  GW_OUTPUT_DT;                    // output GW signals every GW_OUTPUT_DT time interval [0.0]

// Parameters for initial condition
static double *NeutronStar_Prof = NULL;         // radial progenitor model
static int     NeutronStar_NBin;                // number of radial bins in the progenitor model
static char    NeutronStar_ICFile[MAX_STRING];  // Filename for initial condition
// =======================================================================================

static double Mis_InterpolateFromTable_Ext( Profile_t *Phi, const double r );
static void   LoadICTable();
static void   Record_MigrationTest();
static void   Record_CentralDens();
static void   Record_GWSignal_1st();
static void   Record_GWSignal_2nd();

#if ( defined GRAVITY  &&  defined GREP )
extern void   Init_ExtPot_GREP();
extern void   Poi_UserWorkBeforePoisson_GREP( const double Time, const int lv );
#endif




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
   ReadPara->Add( "NeutronStar_ICFile",   NeutronStar_ICFile,    Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "GW_OUTPUT_OPT",       &GW_OUTPUT_OPT,         0,             0,                NoMax_int         );
   ReadPara->Add( "GW_OUTPUT_VERBOSE",   &GW_OUTPUT_VERBOSE,     0,             0,                NoMax_int         );
   ReadPara->Add( "GW_OUTPUT_DT",        &GW_OUTPUT_DT,          0.0,           0.0,              NoMax_double      );
#  ifdef MHD
   ReadPara->Add( "Bfield_Ab",           &Bfield_Ab,             1.0e15,        0.0,              NoMax_double      );
   ReadPara->Add( "Bfield_np",           &Bfield_np,             0.0,           NoMin_double,     NoMax_double      );
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
      Aux_Message( stdout, "  NeutronStar_ICFile        = %s\n",      NeutronStar_ICFile );
      Aux_Message( stdout, "  GW_OUTPUT_OPT             = %d\n",      GW_OUTPUT_OPT );
      Aux_Message( stdout, "  GW_OUTPUT_VERBOSE         = %d\n",      GW_OUTPUT_VERBOSE );
      Aux_Message( stdout, "  GW_OUTPUT_DT              = %13.7e\n",  GW_OUTPUT_DT );
#     ifdef MHD
      Aux_Message( stdout, "  Bfield_Ab                 = %13.7e\n",  Bfield_Ab );
      Aux_Message( stdout, "  Bfield_np                 = %13.7e\n",  Bfield_np );
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
   double momx, momy, momz, eint, etot;

   const double x0 = x - BoxCenter[0];
   const double y0 = y - BoxCenter[1];
   const double z0 = z - BoxCenter[2];
   const double r = SQRT( SQR( x0 ) + SQR( y0 ) + SQR( z0 ) );

   dens = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r);
   velr = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Velr, r);
   pres = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r);

   momx = dens*velr*x0/r;
   momy = dens*velr*y0/r;
   momz = dens*velr*z0/r;

   eint = EoS_DensPres2Eint_CPUPtr( dens, pres, NULL, EoS_AuxArray );   // assuming EoS requires no passive scalars
   etot = Hydro_ConEint2Etot( dens, momx, momy, momz, eint, 0.0 );      // do NOT include magnetic energy here

   fluid[DENS] = dens;
   fluid[MOMX] = momx;
   fluid[MOMY] = momy;
   fluid[MOMZ] = momz;
   fluid[ENGY] = etot;

} // FUNCTION : SetGridIC


#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                2. Generate poloidal B field from vector potential in a form similar
//                   to that defined in Liu+ 2008, Phys. Rev. D78, 024012
//                     A_phi = Ab * \bar\omega^2 * (1 - rho / rho_max)^np * (P / P_max)
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

   // Use the data at first row as the central density and pressure
   // to avoid incorrect values extrapolated from the input IC table
   const double dens_c = Table_Dens[0];
   const double pres_c = Table_Pres[0];
   const double Ab     = Bfield_Ab / UNIT_B;

   // Use finite difference to compute the B field
   double diff = amr->dh[TOP_LEVEL];
   double r,    dens,    pres;
   double r_xp, dens_xp, pres_xp;
   double r_yp, dens_yp, pres_yp;
   double r_zp, dens_zp, pres_zp;

   r       = SQRT( SQR( y0 ) + SQR( z0 ) + SQR( x0        ) );
   r_xp    = SQRT( SQR( y0 ) + SQR( z0 ) + SQR( x0 + diff ) );
   r_yp    = SQRT( SQR( z0 ) + SQR( x0 ) + SQR( y0 + diff ) );
   r_zp    = SQRT( SQR( x0 ) + SQR( y0 ) + SQR( z0 + diff ) );

   dens    = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r   );
   dens_xp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r_xp);
   dens_yp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r_yp);
   dens_zp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Dens, r_zp);

   pres    = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r   );
   pres_xp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r_xp);
   pres_yp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r_yp);
   pres_zp = Mis_InterpolateFromTable(NeutronStar_NBin, Table_R, Table_Pres, r_zp);

   double dAy_dx = ( ( x0 + diff )*POW( 1.0 - dens_xp/dens_c, Bfield_np )*( pres_xp / pres_c )   \
                 -   ( x0        )*POW( 1.0 - dens   /dens_c, Bfield_np )*( pres    / pres_c ) ) \
                 / diff;

   double dAx_dy = ( -( y0 + diff )*POW( 1.0 - dens_yp/dens_c, Bfield_np )*( pres_yp / pres_c )   \
                 -   -( y0        )*POW( 1.0 - dens   /dens_c, Bfield_np )*( pres    / pres_c ) ) \
                 / diff;

   double dAphi_dz = ( POW( 1.0 - dens_zp/dens_c, Bfield_np )*( pres_zp / pres_c )   \
                   -   POW( 1.0 - dens   /dens_c, Bfield_np )*( pres    / pres_c ) ) \
                   / diff;


   magnetic[MAGX] = -x0 * Ab * dAphi_dz;
   magnetic[MAGY] = -y0 * Ab * dAphi_dz;
   magnetic[MAGZ] =       Ab * ( dAy_dx - dAx_dy );

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_GREP_MigrationTest
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_GREP_MigrationTest()
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
   Poi_AddExtraMassForGravity_Ptr = NULL; // option: OPT__GRAVITY_EXTRA_MASS; example: none
#if ( defined GRAVITY  &&  defined GREP )
   Init_ExtPot_Ptr                = Init_ExtPot_GREP;
   Poi_UserWorkBeforePoisson_Ptr  = Poi_UserWorkBeforePoisson_GREP;
#endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr        = NULL; // option: PAR_INIT=1;              example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr    = NULL; // set PAR_NATT_USER;               example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_GREP_MigrationTest



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
// Description :  Interface for calling multiple record functions
//-------------------------------------------------------------------------------------------------------
void Record_MigrationTest()
{

// the maximum density near the domain center
   Record_CentralDens();

// GW Signal
   if ( GW_OUTPUT_OPT )
   {
      double GW_DumpTime = ( GW_OUTPUT_DT > 0.0 ) ? int( Time[0]/GW_OUTPUT_DT )*GW_OUTPUT_DT : Time[0];

//    output data when the elapsed time ~ GW_OUTPUT_DT or dt > GW_OUTPUT_DT
      if ( Time[0] == 0.0  ||  fabs( (Time[0]-GW_DumpTime)/Time[0] ) <= 1.0e-8  ||  dTime_Base >= GW_OUTPUT_DT )
         Record_GWSignal_2nd();
   }

} // FUNCTION : Record_MigrationTest()



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_CentralDens
// Description :  Record the maximum central density
//-------------------------------------------------------------------------------------------------------
void Record_CentralDens()
{

   const char   filename_central_dens[] = "Record__CentralDens";
   const double BoxCenter[3]            = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
// the farthest distance of cells to be considered
   const double r_max2                  = SQR( amr->dh[TOP_LEVEL] );

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
                  const double dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];

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
         for (int b=0; b<4; b++)   DataCoord[b] = OMP_DataCoord[TID][b];
   }

// free per-thread arrays
   Aux_DeallocateArray2D( OMP_DataCoord );


// collect data from all ranks
# ifndef SERIAL
   {
      double DataCoord_All[4 * MPI_NRank];

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

//    output file header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(filename_central_dens) )
         {
             Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_central_dens );
         }

         else
         {
             FILE *file_max_dens = fopen( filename_central_dens, "w" );
             fprintf( file_max_dens, "#%14s %12s %16s %16s %16s %16s\n",
                                     "Time", "Step", "Dens", "PosX", "PosY", "PosZ" );
             fclose( file_max_dens );
         }

         FirstTime = false;
      }

      FILE *file_max_dens = fopen( filename_central_dens, "a" );
      fprintf( file_max_dens, "%15.7e %12ld %16.7e %16.7e %16.7e %16.7e\n",
               Time[0]*UNIT_T, Step, DataCoord[0]*UNIT_D, DataCoord[1]*UNIT_L, DataCoord[2]*UNIT_L, DataCoord[3]*UNIT_L );
      fclose( file_max_dens );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Record_CentralDens()



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_InterpolateFromTable_Ext
// Description :  Extended Mis_InterpolateFromTable() for extrapolation
//-------------------------------------------------------------------------------------------------------
static double Mis_InterpolateFromTable_Ext( Profile_t *Phi, const double r )
{

   const int    NBin = Phi->NBin;
   const double rmin = Phi->Radius[0];
   const double rmax = Phi->Radius[NBin - 1];

   double Phi_interp;

   if      ( r < rmin )   Phi_interp = Phi->Data[0];
   else if ( r < rmax )   Phi_interp = Mis_InterpolateFromTable( NBin, Phi->Radius, Phi->Data, r );
   else                   Phi_interp = Phi->Data[NBin-1];

   return Phi_interp;

}



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_GWSignal_1st
// Description :  Record the first-order time derivative of mass quadrupole moments
//                see Nakamura & Oohara (1989), Oohara et al. (1997)
//-------------------------------------------------------------------------------------------------------
void Record_GWSignal_1st()
{

   const char   filename_QuadMom_1st[] = "Record__QuadMom_1st";
   const double BoxCenter[3]           = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

// allocate memory for per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

// in order of xx, xy, xz, yy, yz, zz
   const int NData = 6;

   double QuadMom_1st[NData] = { 0.0 };
   double **OMP_QuadMom_1st  = NULL;
   Aux_AllocateArray2D( OMP_QuadMom_1st, NT, NData );


#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize arrays
      for (int b=0; b<NData; b++)   OMP_QuadMom_1st[TID][b] = 0.0;

      for (int lv=0; lv<NLEVEL; lv++)
      {
         const double dh = amr->dh[lv];
         const double dv = CUBE( dh );

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

               const double momx = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i];
               const double momy = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i];
               const double momz = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i];

               const double trace = dx * momx + dy * momy + dz * momz;

               OMP_QuadMom_1st[TID][0] += dv * ( 2.0 * dx * momx - (2.0 / 3.0) * trace     );  // xx
               OMP_QuadMom_1st[TID][1] += dv * (       dx * momy +               dy * momx );  // xy
               OMP_QuadMom_1st[TID][2] += dv * (       dx * momz +               dz * momx );  // xz
               OMP_QuadMom_1st[TID][3] += dv * ( 2.0 * dy * momy - (2.0 / 3.0) * trace     );  // yy
               OMP_QuadMom_1st[TID][4] += dv * (       dy * momz +               dz * momy );  // yz
               OMP_QuadMom_1st[TID][5] += dv * ( 2.0 * dz * momz - (2.0 / 3.0) * trace     );  // zz

            }}} // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // for (int lv=0; lv<NLEVEL; lv++)
   } // OpenMP parallel region


// sum over all OpenMP threads
   for (int b=0; b<NData; b++)
   for (int t=0; t<NT; t++)
   {
      QuadMom_1st[b] += OMP_QuadMom_1st[t][b];
   }

// free per-thread arrays
   Aux_DeallocateArray2D( OMP_QuadMom_1st );


// collect data from all ranks (in-place reduction)
#  ifndef SERIAL
   if ( MPI_Rank == 0 )   MPI_Reduce( MPI_IN_PLACE, QuadMom_1st, NData, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
   else                   MPI_Reduce( QuadMom_1st,  NULL,        NData, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#  endif // ifndef SERIAL


// multiply the coefficient and unit
   const double UNIT_QuadMom_1st = UNIT_M * UNIT_L * UNIT_V;
   const double coe = 2.0 * Const_NewtonG / pow( Const_c, 4.0 );

   for (int b=0; b<NData; b++)   QuadMom_1st[b] *= coe * UNIT_QuadMom_1st;


// output to file
   if ( MPI_Rank == 0 )
   {

      static bool FirstTime = true;

//    output file header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(filename_QuadMom_1st) )
         {
             Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_QuadMom_1st );
         }
         else
         {
             FILE *file_QuadMom_1st = fopen( filename_QuadMom_1st, "w" );
             fprintf( file_QuadMom_1st, "#%14s %7s %16s %16s %16s %16s %16s %16s\n",
                                        "Time", "Step", "xx", "xy", "xz", "yy", "yz", "zz" );
             fclose( file_QuadMom_1st );
         }

         FirstTime = false;
      }

      FILE *file_QuadMom_1st = fopen( filename_QuadMom_1st, "a" );

                                    fprintf( file_QuadMom_1st, "%15.7e %7ld", Time[0] * UNIT_T, Step );
      for (int b=0; b<NData; b++)   fprintf( file_QuadMom_1st, "%17.7e", QuadMom_1st[b] );
                                    fprintf( file_QuadMom_1st, "\n" );

      fclose( file_QuadMom_1st );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Record_GWSignal_1st()



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_GWSignal_2nd
// Description :  Record the second-order time derivative of mass quadrupole moments
//                see Nakamura & Oohara (1989), Oohara et al. (1997)
//-------------------------------------------------------------------------------------------------------
void Record_GWSignal_2nd()
{

#  if ( defined GRAVITY  &&  defined GREP )

   const char   filename_QuadMom_2nd[ ] = "Record__QuadMom_2nd";
   const double BoxCenter           [3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

// allocate memory for per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

// in order of xx, xy, xz, yy, yz, zz
   const int NData   = 6;
   const int ArrayID = 0;
   const int NPG_Max = POT_GPU_NPGROUP;

   double QuadMom_2nd[NData] = { 0.0 };
   double **OMP_QuadMom_2nd  = NULL;
   Aux_AllocateArray2D( OMP_QuadMom_2nd, NT, NData );


// compute and output the spherically summed GW signal if GW_OUTPUT_VERBOSE is enabled
   Profile_t *QuadMom_Prof[NData];
   double dr_min, r_max, Center[3];
   double ***OMP_Data =NULL;
   long   ***OMP_NCell=NULL;


   if ( GW_OUTPUT_VERBOSE )
   {
//    use GREP parameters to set up the GW signal profiles
      switch ( GREP_CENTER_METHOD )
      {
         case 1:
            for (int i=0; i<3; i++)   Center[i] = amr->BoxCenter[i];
         break;

         default:
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GREP_CENTER_METHOD", GREP_CENTER_METHOD );
      }

      dr_min = ( GREP_MINBINSIZE > 0.0 ) ? GREP_MINBINSIZE : amr->dh[MAX_LEVEL];
      r_max  = ( GREP_MAXRADIUS  > 0.0 ) ? GREP_MAXRADIUS
                                         : SQRT( SQR( MAX( amr->BoxSize[0] - Center[0], Center[0] ) )
                                         +       SQR( MAX( amr->BoxSize[1] - Center[1], Center[1] ) )
                                         +       SQR( MAX( amr->BoxSize[2] - Center[2], Center[2] ) ) );

//    initialize the profile objects
      for (int p=0; p<NData; p++)
      {

         QuadMom_Prof[p] = new Profile_t();

//       get the total number of radial bins and the corresponding maximum radius
         if ( GREP_LOGBIN )
         {
            QuadMom_Prof[p]->NBin      = int( log(r_max/dr_min)/log(GREP_LOGBINRATIO) ) + 2;
            QuadMom_Prof[p]->MaxRadius = dr_min*pow( GREP_LOGBINRATIO, QuadMom_Prof[p]->NBin-1 );
         }

         else // linear bin
         {
            QuadMom_Prof[p]->NBin      = (int)ceil( r_max / dr_min );
            QuadMom_Prof[p]->MaxRadius = dr_min*QuadMom_Prof[p]->NBin;
         }

//       record profile parameters
         for (int d=0; d<3; d++)   QuadMom_Prof[p]->Center[d] = Center[d];

         QuadMom_Prof[p]->LogBin = GREP_LOGBIN;

         if ( GREP_LOGBIN )   QuadMom_Prof[p]->LogBinRatio = GREP_LOGBINRATIO;

         QuadMom_Prof[p]->AllocateMemory();

//       record radial coordinates
         if ( GREP_LOGBIN )
            for (int b=0; b<QuadMom_Prof[0]->NBin; b++)   QuadMom_Prof[p]->Radius[b] = dr_min*pow( GREP_LOGBINRATIO, b-0.5 );
         else
            for (int b=0; b<QuadMom_Prof[0]->NBin; b++)   QuadMom_Prof[p]->Radius[b] = (b+0.5)*dr_min;

      } // for (int p=0; p<NData; p++)


      Aux_AllocateArray3D( OMP_Data,  NData, NT, QuadMom_Prof[0]->NBin );
      Aux_AllocateArray3D( OMP_NCell, NData, NT, QuadMom_Prof[0]->NBin );
   }


#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize arrays
      for (int b=0; b<NData; b++)   OMP_QuadMom_2nd[TID][b] = 0.0;

      if ( GW_OUTPUT_VERBOSE )
      {
         for (int p=0; p<NData; p++)
         for (int b=0; b<QuadMom_Prof[0]->NBin; b++)
         {
            OMP_Data  [p][TID][b] = 0.0;
            OMP_NCell [p][TID][b] = 0;
         }
      }


      for (int lv=0; lv<NLEVEL; lv++)
      {

         const double dh = amr->dh[lv];
         const double dv = CUBE( dh );

         const double TimeNew   = Time[lv];
         const int    NTotal    = amr->NPatchComma[lv][1] / 8;
               int   *PID0_List = new int [NTotal];

         for (int t=0; t<NTotal; t++)   PID0_List[t] = 8*t;


         for (int Disp=0; Disp<NTotal; Disp+=NPG_Max)
         {

            int NPG = ( NPG_Max < NTotal-Disp ) ? NPG_Max : NTotal-Disp;

            Prepare_PatchData( lv, TimeNew, &h_Pot_Array_P_Out[ArrayID][0][0][0][0], NULL,
                               GRA_GHOST_SIZE, NPG, PID0_List+Disp, _POTE, _NONE,
                               OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, (GRA_GHOST_SIZE==0)?NSIDE_00:NSIDE_06, false,
                               OPT__BC_FLU, OPT__BC_POT, -1.0, -1.0, false );

#        pragma omp for schedule( runtime )
         for (int PID_IDX=0; PID_IDX<8*NPG; PID_IDX++)
         {
            int PID = 8*Disp + PID_IDX;

            if ( amr->patch[0][lv][PID]->son != -1 )  continue;

            for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh; const int kk = k + GRA_GHOST_SIZE;
            for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh; const int jj = j + GRA_GHOST_SIZE;
            for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh; const int ii = i + GRA_GHOST_SIZE;

               const double dx = x - BoxCenter[0];
               const double dy = y - BoxCenter[1];
               const double dz = z - BoxCenter[2];
               const double r = SQRT( SQR(dx) + SQR(dy) + SQR(dz) );

               const double dens  = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
               const double momx  = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i];
               const double momy  = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i];
               const double momz  = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i];
               const double _dens = 1.0 / dens;

               const real (*PrepPotPtr)[GRA_NXT][GRA_NXT] = h_Pot_Array_P_Out[ArrayID][PID_IDX];

               const double dPhi_dx = ( PrepPotPtr[kk  ][jj  ][ii+1] - PrepPotPtr[kk  ][jj  ][ii-1] ) / (2.0 * dh);
               const double dPhi_dy = ( PrepPotPtr[kk  ][jj+1][ii  ] - PrepPotPtr[kk  ][jj-1][ii  ] ) / (2.0 * dh);
               const double dPhi_dz = ( PrepPotPtr[kk+1][jj  ][ii  ] - PrepPotPtr[kk-1][jj  ][ii  ] ) / (2.0 * dh);

               const double trace = _dens * ( SQR(momx) + SQR(momy) + SQR(momz) )
                                  -  dens * ( dx * dPhi_dx + dy * dPhi_dy + dz * dPhi_dz );

               const double QuadMom_xx = dv * ( 2.0 * _dens * momx * momx - (2.0 / 3.0) * trace
                                              - 2.0 *  dens * dx * dPhi_dx                      );  // xx
               const double QuadMom_xy = dv * ( 2.0 * _dens * momx * momy
                                              -        dens * ( dx * dPhi_dy + dy * dPhi_dx )   );  // xy
               const double QuadMom_xz = dv * ( 2.0 * _dens * momx * momz
                                              -        dens * ( dx * dPhi_dz + dz * dPhi_dx )   );  // xz
               const double QuadMom_yy = dv * ( 2.0 * _dens * momy * momy - (2.0 / 3.0) * trace
                                              - 2.0 *  dens * dy * dPhi_dy                      );  // yy
               const double QuadMom_yz = dv * ( 2.0 * _dens * momy * momz
                                              -        dens * ( dy * dPhi_dz + dz * dPhi_dy )   );  // yz
               const double QuadMom_zz = dv * ( 2.0 * _dens * momz * momz - (2.0 / 3.0) * trace
                                              - 2.0 *  dens * dz * dPhi_dz                      );  // zz

               OMP_QuadMom_2nd[TID][0] += QuadMom_xx;
               OMP_QuadMom_2nd[TID][1] += QuadMom_xy;
               OMP_QuadMom_2nd[TID][2] += QuadMom_xz;
               OMP_QuadMom_2nd[TID][3] += QuadMom_yy;
               OMP_QuadMom_2nd[TID][4] += QuadMom_yz;
               OMP_QuadMom_2nd[TID][5] += QuadMom_zz;


//             store the spherically summed profiles
               if ( GW_OUTPUT_VERBOSE )
               {
                  const double dx = x - Center[0];
                  const double dy = y - Center[1];
                  const double dz = z - Center[2];
                  const double r = SQRT( SQR(dx) + SQR(dy) + SQR(dz) );

                  if ( r <= QuadMom_Prof[0]->MaxRadius )
                  {
                     const int bin = ( GREP_LOGBIN ) ? (  (r<dr_min) ? 0 : int( log(r/dr_min)/log(GREP_LOGBINRATIO) ) + 1  )
                                                     : int( r/dr_min );
//                   prevent from round-off errors
                     if ( bin >= QuadMom_Prof[0]->NBin )   continue;

//                   store the data only
                     OMP_Data[0][TID][bin] += QuadMom_xx;
                     OMP_Data[1][TID][bin] += QuadMom_xy;
                     OMP_Data[2][TID][bin] += QuadMom_xz;
                     OMP_Data[3][TID][bin] += QuadMom_yy;
                     OMP_Data[4][TID][bin] += QuadMom_yz;
                     OMP_Data[5][TID][bin] += QuadMom_zz;

                     for (int p=0; p<NData; p++)   OMP_NCell [p][TID][bin] ++;
                  }
               }

            }}} // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
#        pragma omp barrier
         } // for (int Disp=0; Disp<NTotal; Disp+=NPG_Max)
      } // for (int lv=0; lv<NLEVEL; lv++)
   } // OpenMP parallel region


// sum over all OpenMP threads
   for (int b=0; b<NData; b++)
   for (int t=0; t<NT; t++)
   {
      QuadMom_2nd[b] += OMP_QuadMom_2nd[t][b];
   }

// free per-thread arrays
   Aux_DeallocateArray2D( OMP_QuadMom_2nd );


// collect data from all ranks (in-place reduction)
#  ifndef SERIAL
   if ( MPI_Rank == 0 )   MPI_Reduce( MPI_IN_PLACE, QuadMom_2nd, NData, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
   else                   MPI_Reduce( QuadMom_2nd,  NULL,        NData, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#  endif // ifndef SERIAL


// multiply the coefficient and unit
   const double UNIT_QuadMom_2nd = UNIT_M * UNIT_V * UNIT_V;
   const double coe = 2.0 * Const_NewtonG / pow( Const_c, 4.0 );

   for (int b=0; b<NData; b++)   QuadMom_2nd[b] *= coe * UNIT_QuadMom_2nd;


// output to file
   if ( MPI_Rank == 0 )
   {

      static bool FirstTime = true;

//    output file header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(filename_QuadMom_2nd) )
         {
             Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_QuadMom_2nd );
         }
         else
         {
             FILE *file_QuadMom_2nd = fopen( filename_QuadMom_2nd, "w" );
             fprintf( file_QuadMom_2nd, "#%14s %7s %16s %16s %16s %16s %16s %16s\n",
                                        "Time", "Step", "xx", "xy", "xz", "yy", "yz", "zz" );
             fclose( file_QuadMom_2nd );
         }

         FirstTime = false;
      }

      FILE *file_QuadMom_2nd = fopen( filename_QuadMom_2nd, "a" );

                                    fprintf( file_QuadMom_2nd, "%15.7e %7ld", Time[0] * UNIT_T, Step );
      for (int b=0; b<NData; b++)   fprintf( file_QuadMom_2nd, "%17.7e", QuadMom_2nd[b] );
                                    fprintf( file_QuadMom_2nd, "\n" );

      fclose( file_QuadMom_2nd );

   } // if ( MPI_Rank == 0 )


// postprocess and output the spherically summed profiles
   if ( GW_OUTPUT_VERBOSE )
   {
//    sum over all OpenMP threads
      for (int p=0; p<NData; p++)
      {
         for (int b=0; b<QuadMom_Prof[0]->NBin; b++)
         {
            QuadMom_Prof[p]->Data  [b]  = OMP_Data  [p][0][b];
            QuadMom_Prof[p]->NCell [b]  = OMP_NCell [p][0][b];
         }

         for (int t=1; t<NT; t++)
         for (int b=0; b<QuadMom_Prof[0]->NBin; b++)
         {
            QuadMom_Prof[p]->Data  [b] += OMP_Data  [p][t][b];
            QuadMom_Prof[p]->NCell [b] += OMP_NCell [p][t][b];
         }
      }

      Aux_DeallocateArray3D( OMP_Data );

//    collect data from all ranks (in-place reduction)
#     ifndef SERIAL
      for (int p=0; p<NData; p++)
      {
         const int NBin = QuadMom_Prof[p]->NBin;

         if ( MPI_Rank == 0 )
         {
            MPI_Reduce( MPI_IN_PLACE,           QuadMom_Prof[p]->Data,  NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce( MPI_IN_PLACE,           QuadMom_Prof[p]->NCell, NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
         }

         else
         {
            MPI_Reduce( QuadMom_Prof[p]->Data,  NULL,                   NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce( QuadMom_Prof[p]->NCell, NULL,                   NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
         }
      }
#     endif


// remove the empty bins
      for (int b=0; b<QuadMom_Prof[0]->NBin; b++)
      {
         if ( QuadMom_Prof[0]->NCell[b] != 0L )   continue;

//       remove consecutive empty bins at the same time for better performance
         int b_up;
         for (b_up=b+1; b_up<QuadMom_Prof[0]->NBin; b_up++)
            if ( QuadMom_Prof[0]->NCell[b_up] != 0L )   break;

         const int stride = b_up - b;

         for (b_up=b+stride; b_up<QuadMom_Prof[0]->NBin; b_up++)
         {
            const int b_up_ms = b_up - stride;

            for (int p=0; p<NData; p++)
            {
               QuadMom_Prof[p]->Radius[b_up_ms] = QuadMom_Prof[p]->Radius[b_up];
               QuadMom_Prof[p]->Data  [b_up_ms] = QuadMom_Prof[p]->Data  [b_up];
               QuadMom_Prof[p]->NCell [b_up_ms] = QuadMom_Prof[p]->NCell [b_up];
            }
         }

//       reset the total number of bins
         for (int p=0; p<NData; p++)   QuadMom_Prof[p]->NBin -= stride;
      } // for (int b=0; b<QuadMom_Prof->NBin; b++)


//    output to file
      if ( MPI_Rank == 0 )
      {

         char filename_QuadMom_Prof[50];
         sprintf( filename_QuadMom_Prof, "Record__QuadMom_Prof_%06ld", Step );

//       file header
         if ( Aux_CheckFileExist(filename_QuadMom_Prof) )
         {
             Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_QuadMom_Prof );
         }
         else
         {
             FILE *file_QuadMom_Prof = fopen( filename_QuadMom_Prof, "w" );
             fprintf( file_QuadMom_Prof, "# Time : %13.7e,  Step : %8ld\n", Time[0] * UNIT_T, Step );
             fprintf( file_QuadMom_Prof, "#%4s %8s %15s %15s %15s %15s %15s %15s %15s\n",
                                         "Bin", "NCell", "Radius", "xx", "xy", "xz", "yy", "yz", "zz" );
             fclose( file_QuadMom_Prof );
         }

         FILE *file_QuadMom_Prof = fopen( filename_QuadMom_Prof, "a" );

         for (int b=0; b<QuadMom_Prof[0]->NBin; b++)
         {
                                          fprintf( file_QuadMom_Prof, "%5d %8ld %15.7e",
                                                   b, QuadMom_Prof[0]->NCell[b], QuadMom_Prof[0]->Radius[b] * UNIT_L );
            for (int p=0; p<NData; p++)   fprintf( file_QuadMom_Prof, "%16.7e",
                                                   QuadMom_Prof[p]->Data[b] * coe * UNIT_QuadMom_2nd );
                                          fprintf( file_QuadMom_Prof, "\n" );
         }

         fclose( file_QuadMom_Prof );

      } // if ( MPI_Rank == 0 )

//    free memory
      for (int p=0; p<NData; p++)   QuadMom_Prof[p]->FreeMemory();
   } // if ( GW_OUTPUT_VERBOSE )

#  endif // if ( defined GRAVITY  &&  defined GREP )

} // FUNCTION : Record_GWSignal_2nd()

