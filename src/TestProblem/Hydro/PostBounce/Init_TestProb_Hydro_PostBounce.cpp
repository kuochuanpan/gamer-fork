#include "GAMER.h"
#include "TestProb.h"
#include "NuclearEos.h"


// problem-specific global variables
// =======================================================================================
static double  Bfield_Ab;                       // magnetic field strength                          [1e15]
static double  Bfield_np;                       // dependence on the density                        [0.0]

// Parameters for initial condition
static double *NeutronStar_Prof = NULL;         // radial progenitor model
static char    NeutronStar_ICFile[MAX_STRING];  // Filename for initial condition
static int     NeutronStar_NBin;                // number of radial bins in the progenitor model
// =======================================================================================


static double Mis_InterpolateFromTable_Ext( Profile_t *Phi, const double r );

static void LoadICTable();
static void Record_PostBounce();
static void Record_GWSignal_1st();
static void Record_GWSignal_2nd();
static void Record_CentralDens();
static bool Flag_User_PostBounce( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );


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
      Aux_Message( stdout, "  test problem ID           = %d\n",      TESTPROB_ID          );
      Aux_Message( stdout, "  NeutronStar_ICFile        = %s\n",      NeutronStar_ICFile   );
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
   const double *Table_R       = NeutronStar_Prof + 0*NeutronStar_NBin;
   const double *Table_Dens    = NeutronStar_Prof + 1*NeutronStar_NBin;
   const double *Table_Ye      = NeutronStar_Prof + 2*NeutronStar_NBin;
   const double *Table_Velr    = NeutronStar_Prof + 3*NeutronStar_NBin;
   const double *Table_Temp    = NeutronStar_Prof + 4*NeutronStar_NBin;
   const double *Table_Pres    = NeutronStar_Prof + 5*NeutronStar_NBin;
   const double *Table_Entropy = NeutronStar_Prof + 6*NeutronStar_NBin;

   const double x0 = x - BoxCenter[0];
   const double y0 = y - BoxCenter[1];
   const double z0 = z - BoxCenter[2];
   const double r  = SQRT( SQR( x0 ) + SQR( y0 ) + SQR( z0 ) );

   double dens, ye, velr, temp, pres, entropy;

   dens    = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Dens,    r );
   ye      = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Ye,      r );
   velr    = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Velr,    r );
   temp    = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Temp,    r );  // in Kelvin
   pres    = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Pres,    r );
   entropy = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Entropy, r );

   fluid[DENS] = dens;
   fluid[MOMX] = dens*velr*x0/r;
   fluid[MOMY] = dens*velr*y0/r;
   fluid[MOMZ] = dens*velr*z0/r;
   fluid[ENGY] = pres / ( GAMMA - 1.0 )
               + 0.5*( SQR( fluid[MOMX] ) + SQR( fluid[MOMY] ) + SQR( fluid[MOMZ] ) ) / dens;
#  if ( EOS == NUCLEAR )
   fluid[YE]   = ye*dens;       // electron fraction [dens]
   fluid[ENTR] = entropy*dens;  // entropy [kB/baryon * dens]
#  endif


// reset the energy and entropy from nuclear EoS table
#  if ( EOS == NUCLEAR )
// get the energy and entropy from Nuclear EoS table (use temperature mode 1)
   double xtmp, xenr, xprs, xent, xcs2, xdedt, xdpderho, xdpdrhoe, xmunu;
   int keyerr;

   const double rfeps = 1.0e-10;
   const double kelvin_to_mev = Const_kB_eV / 1.0e6;
   const double mev_to_kelvin = 1.1604447522806e10;

// compute temperature using ideal EoS: P = \rho*K*T / ( mu*m_H )
//   xtmp = (pres*UNIT_P) / (dens*UNIT_D) * MOLECULAR_WEIGHT * Const_mH / Const_kB;  // in Kelvin
   xtmp = temp / mev_to_kelvin;  // to MeV

   nuc_eos_C_short((dens*UNIT_D), &xtmp, ye, &xenr, &xprs, &xent, &xcs2,
                   &xdedt, &xdpderho, &xdpdrhoe, &xmunu,
                   1, &keyerr, rfeps);

   if (keyerr != 0)   printf("debug: keyerr not zero %d\n",keyerr);

   fluid[YE]   = ye*dens;    // electron fraction [dens]
   fluid[ENTR] = xent*dens;  // entropy [kB/baryon * dens]
   fluid[ENGY] = ( dens/(UNIT_V*UNIT_V) )*( xenr + energy_shift )
               + 0.5*( SQR( fluid[MOMX] ) + SQR( fluid[MOMY] ) + SQR( fluid[MOMZ] ) ) / dens;
#  endif

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
   const double *Table_Dens   = NeutronStar_Prof + 1*NeutronStar_NBin;
   const double *Table_Pres   = NeutronStar_Prof + 5*NeutronStar_NBin;

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

   dens    = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Dens, r    );
   dens_xp = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Dens, r_xp );
   dens_yp = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Dens, r_yp );
   dens_zp = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Dens, r_zp );

   pres    = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Pres, r    );
   pres_xp = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Pres, r_xp );
   pres_yp = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Pres, r_yp );
   pres_zp = Mis_InterpolateFromTable( NeutronStar_NBin, Table_R, Table_Pres, r_zp );

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
// Function    :  Init_TestProb_Hydro_PostBounce
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_PostBounce()
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
   Flag_User_Ptr                  = Flag_User_PostBounce; // option: OPT__FLAG_USER;          example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr       = NULL; // option: OPT__DT_USER;            example: Miscellaneous/Mis_GetTimeStep_User.cpp
   BC_User_Ptr                    = NULL; // option: OPT__BC_FLU_*=4;         example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr       = NULL; // option: OPT__RESET_FLUID;        example: Fluid/Flu_ResetByUser.cpp
   Output_User_Ptr                = NULL; // option: OPT__OUTPUT_USER;        example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr            = Record_PostBounce; // option: OPT__RECORD_USER;        example: Auxiliary/Aux_Record_User.cpp
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

} // FUNCTION : Init_TestProb_Hydro_PostBounce



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadICTable
// Description :  Load inpu table file for initial condition
//                Temperature unit is UNIT_P / UNIT_D / Const_kB
//-------------------------------------------------------------------------------------------------------
void LoadICTable()
{

   const bool RowMajor_No  = false;           // load data into the column major OPT__RECORD_USER
   const bool AllocMem_Yes = true;            // allocate memort for NeutronStar_Prof
   const int  NCol         = 7;               // total number of columns to load

   // target columns: {radius, density, ye, radial velocity, temperature, pressure, entropy }
   const int  TargetCols[NCol] = { 0, 2, 3, 4, 5, 6, 7 };

   double *Table_R, *Table_Dens, *Table_Ye, *Table_Velr, *Table_Temp, *Table_Pres, *Table_Entropy;

   NeutronStar_NBin = Aux_LoadTable( NeutronStar_Prof, NeutronStar_ICFile, NCol, TargetCols, RowMajor_No, AllocMem_Yes );

   Table_R       = NeutronStar_Prof + 0*NeutronStar_NBin;
   Table_Dens    = NeutronStar_Prof + 1*NeutronStar_NBin;
   Table_Ye      = NeutronStar_Prof + 2*NeutronStar_NBin;
   Table_Velr    = NeutronStar_Prof + 3*NeutronStar_NBin;
   Table_Temp    = NeutronStar_Prof + 4*NeutronStar_NBin;
   Table_Pres    = NeutronStar_Prof + 5*NeutronStar_NBin;
   Table_Entropy = NeutronStar_Prof + 6*NeutronStar_NBin;

   // convert to code units (assuming progentior model is in cgs)
   for (int b=0; b<NeutronStar_NBin; b++)
   {
      Table_R[b]    /= UNIT_L;
      Table_Dens[b] /= UNIT_D;
      Table_Velr[b] /= UNIT_V;
      Table_Pres[b] /= UNIT_P;
   }

} // FUNCTION : LoadICTable()



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_PostBounce
// Description :  Interface for calling multiple record functions
//-------------------------------------------------------------------------------------------------------
void Record_PostBounce()
{

   Record_CentralDens();
   Record_GWSignal_1st();
   Record_GWSignal_2nd();

} // FUNCTION : Record_PostBounce()



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_CentralDens
// Description :  Record the central density
//-------------------------------------------------------------------------------------------------------
void Record_CentralDens()
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
//       write header before the first output
         if ( Aux_CheckFileExist(filename_central_dens) )
         {
             Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_central_dens );
         }
         else
         {
             FILE *file_max_dens = fopen( filename_central_dens, "w" );
             fprintf( file_max_dens, "#%19s %12s %16s %16s %16s %16s\n",
                                     "Time", "Step", "Dens", "Posx", "Posy", "Posz" );
             fclose( file_max_dens );
         }

         FirstTime = false;
      }

      FILE *file_max_dens = fopen( filename_central_dens, "a" );
      fprintf( file_max_dens,
               "%20.14e   %10ld   %14.7e   %14.7e   %14.7e   %14.7e\n",
              Time[0] * UNIT_T, Step, DataCoord[0], DataCoord[1], DataCoord[2], DataCoord[3] );
      fclose( file_max_dens );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Record_CentralDens()




//-------------------------------------------------------------------------------------------------------
// Function    :  Record_GWSignal_1st
// Description :  Record the first-order derivative of mass quadrupole moments
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
   double **OMP_QuadMom_1st=NULL;
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
   for (int b=0; b<NData; b++) {
   for (int t=0; t<NT; t++)    {
      QuadMom_1st[b] += OMP_QuadMom_1st[t][b];
   }}

// free per-thread arrays
   Aux_DeallocateArray2D( OMP_QuadMom_1st );



// collect data from all ranks (in-place reduction)
#  ifndef SERIAL
   if ( MPI_Rank == 0 )   MPI_Reduce( MPI_IN_PLACE, QuadMom_1st, NData, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
   else                   MPI_Reduce( QuadMom_1st,  NULL,        NData, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#  endif // ifndef SERIAL


// multiply the coefficient and unit
   const double UNIT_QuadMom_1st = UNIT_M * UNIT_L * UNIT_V;
   const double coe = 4.0 * Const_NewtonG / pow( Const_c, 4.0 );

   for (int b=0; b<NData; b++)   QuadMom_1st[b] *= coe * UNIT_QuadMom_1st;


// output to file
   if ( MPI_Rank == 0 )
   {

      static bool FirstTime = true;

      if ( FirstTime )
      {
//       write header before the first output
         if ( Aux_CheckFileExist(filename_QuadMom_1st) )
         {
             Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_QuadMom_1st );
         }
         else
         {
             FILE *file_QuadMom_1st = fopen( filename_QuadMom_1st, "w" );
             fprintf( file_QuadMom_1st, "#%19s %12s %16s %16s %16s %16s %16s %16s\n",
                                        "Time", "Step", "xx", "xy", "xz", "yy", "yz", "zz" );
             fclose( file_QuadMom_1st );
         }

         FirstTime = false;
      }

      FILE *file_QuadMom_1st = fopen( filename_QuadMom_1st, "a" );

                                    fprintf( file_QuadMom_1st, "%20.14e %12ld", Time[0] * UNIT_T, Step  );
      for (int b=0; b<NData; b++)   fprintf( file_QuadMom_1st, "%17.7e",        QuadMom_1st[b] );
                                    fprintf( file_QuadMom_1st, "\n"                            );

      fclose( file_QuadMom_1st );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Record_GWSignal_1st()



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

   if ( r < rmin )
      Phi_interp = Phi->Data[0];
   else if ( r < rmax )
      Phi_interp = Mis_InterpolateFromTable( NBin, Phi->Radius, Phi->Data, r );
   else
      Phi_interp = Phi->Data[NBin-1];

   return Phi_interp;

}



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_GWSignal_2nd
// Description :  Record the second-order derivative of mass quadrupole moments
//                see Nakamura & Oohara (1989), Oohara et al. (1997), eq. (6) in Andresen et al. (2017)
//-------------------------------------------------------------------------------------------------------
void Record_GWSignal_2nd()
{

#  if ( defined GRAVITY  &&  defined GREP )

   const char   filename_QuadMom_2nd[ ] = "Record__QuadMom_2nd";
   const double BoxCenter           [3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
         double Center              [3] = { 0.0 };
   const long   TVar                [ ] = { _DENS, _EINT_DER, _VELR, _PRES };
         double MaxRadius, MinBinSize;
   Profile_t *Phi_eff, *ProfAve[4];


// compute the GR effective potential
   switch ( GREP_CENTER_METHOD )
   {
      case 1:   for (int i=0; i<3; i++)   Center[i] = amr->BoxCenter[i];
                break;
      default:  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GREP_CENTER_METHOD", GREP_CENTER_METHOD );
   }

//    Defaults to the distance between the center and the farthest box vertex
   MaxRadius  = ( GREP_MAXRADIUS > 0.0 )  ? GREP_MAXRADIUS
                                          : SQRT( SQR( MAX( amr->BoxSize[0] - Center[0], Center[0] ) )
                                          +       SQR( MAX( amr->BoxSize[1] - Center[1], Center[1] ) )
                                          +       SQR( MAX( amr->BoxSize[2] - Center[2], Center[2] ) ));

   MinBinSize = ( GREP_MINBINSIZE > 0.0 ) ? GREP_MINBINSIZE
                                          : amr->dh[MAX_LEVEL];

                                         Phi_eff        = new Profile_t();
   for (int NProf=0; NProf<4; NProf++)   ProfAve[NProf] = new Profile_t();

   Aux_ComputeProfile( ProfAve, Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO, true, TVar, 4, -1 );
   CPU_ComputeEffPot ( ProfAve[0], ProfAve[1], ProfAve[2], ProfAve[3], Phi_eff );


   const int    NBin = Phi_eff->NBin;
   const double rmin = Phi_eff->Radius[0];
   const double rmax = Phi_eff->Radius[NBin-1];


// allocate memory for per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

// in order of xx, xy, xz, yy, yz, zz
   const int NData = 6;

   double QuadMom_2nd[NData] = { 0.0 };
   double **OMP_QuadMom_2nd=NULL;
   Aux_AllocateArray2D( OMP_QuadMom_2nd, NT, NData );


#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize arrays
      for (int b=0; b<NData; b++)   OMP_QuadMom_2nd[TID][b] = 0.0;

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

               const double dens  = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
               const double momx  = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i];
               const double momy  = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i];
               const double momz  = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i];
               const double _dens = 1.0 / dens;

               const double r = SQRT( SQR(dx) + SQR(dy) + SQR(dz) );
               const double Phi_eff_r = Mis_InterpolateFromTable_Ext( Phi_eff, r );

               double dPhi_dx, dPhi_dy, dPhi_dz;

//             dPhi_dx
               if ( i < PS1 -1 )
               {
                  const double rpx = SQRT( SQR(dx + dh) + SQR(dy) + SQR(dz) );
                  const double Phi_eff_rpx = Mis_InterpolateFromTable_Ext( Phi_eff, rpx );

                  dPhi_dx = ( amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j][i+1] + Phi_eff_rpx
                            - amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j][i]   - Phi_eff_r   ) / dh;
               }
               else
               {
                  const double rpx = SQRT( SQR(dx - dh) + SQR(dy) + SQR(dz) );
                  const double Phi_eff_rpx = Mis_InterpolateFromTable_Ext( Phi_eff, rpx );

                  dPhi_dx = ( amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j][i]   + Phi_eff_r
                            - amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j][i-1] - Phi_eff_rpx  ) / dh;
               }

//             dPhi_dy
               if ( j < PS1 -1 )
               {
                  const double rpy = SQRT( SQR(dx) + SQR(dy + dh) + SQR(dz) );
                  const double Phi_eff_rpy = Mis_InterpolateFromTable_Ext( Phi_eff, rpy );

                  dPhi_dy = ( amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j+1][i] + Phi_eff_rpy
                            - amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j][i]   - Phi_eff_r   ) / dh;
               }
               else
               {
                  const double rpy = SQRT( SQR(dx) + SQR(dy - dh) + SQR(dz) );
                  const double Phi_eff_rpy = Mis_InterpolateFromTable_Ext( Phi_eff, rpy );

                  dPhi_dy = ( amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j][i]   + Phi_eff_r
                            - amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j-1][i] - Phi_eff_rpy  ) / dh;
               }

//             dPhi_dz
               if ( k < PS1 -1 )
               {
                  const double rpz = SQRT( SQR(dx) + SQR(dy) + SQR(dz + dh) );
                  const double Phi_eff_rpz = Mis_InterpolateFromTable_Ext( Phi_eff, rpz );

                  dPhi_dz = ( amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k+1][j][i] + Phi_eff_rpz
                            - amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j][i]   - Phi_eff_r   ) / dh;
               }
               else
               {
                  const double rpz = SQRT( SQR(dx) + SQR(dy) + SQR(dz - dh) );
                  const double Phi_eff_rpz = Mis_InterpolateFromTable_Ext( Phi_eff, rpz );

                  dPhi_dz = ( amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k][j][i]   + Phi_eff_r
                            - amr->patch[ amr->FluSg[lv] ][lv][PID]->pot[k-1][j][i] - Phi_eff_rpz  ) / dh;
               }

               const double trace = _dens * ( SQR(momx) + SQR(momy) + SQR(momz) )
                                  -  dens * ( dx * dPhi_dx + dy * dPhi_dy + dz * dPhi_dz );

               OMP_QuadMom_2nd[TID][0] += dv * (       _dens * momx * momx - (1.0 / 3.0) * trace
                                               -        dens * dx * dPhi_dx                      );  // xx
               OMP_QuadMom_2nd[TID][1] += dv * (       _dens * momx * momy
                                               - 0.5 *  dens * ( dx * dPhi_dy + dy * dPhi_dx)    );  // xy
               OMP_QuadMom_2nd[TID][2] += dv * (       _dens * momx * momz
                                               - 0.5 *  dens * ( dx * dPhi_dz + dz * dPhi_dx)    );  // xz
               OMP_QuadMom_2nd[TID][3] += dv * (       _dens * momy * momy - (1.0 / 3.0) * trace
                                               -        dens * dy * dPhi_dy                      );  // yy
               OMP_QuadMom_2nd[TID][4] += dv * (       _dens * momy * momz
                                               - 0.5 *  dens * ( dy * dPhi_dz + dz * dPhi_dy)    );  // yz
               OMP_QuadMom_2nd[TID][5] += dv * (       _dens * momz * momz - (1.0 / 3.0) * trace
                                               -        dens * dz * dPhi_dz                      );  // zz

            }}} // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // for (int lv=0; lv<NLEVEL; lv++)
   } // OpenMP parallel region



// sum over all OpenMP threads
   for (int b=0; b<NData; b++) {
   for (int t=0; t<NT; t++)    {
      QuadMom_2nd[b] += OMP_QuadMom_2nd[t][b];
   }}

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

      if ( FirstTime )
      {
//       write header before the first output
         if ( Aux_CheckFileExist(filename_QuadMom_2nd) )
         {
             Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_QuadMom_2nd );
         }
         else
         {
             FILE *file_QuadMom_2nd = fopen( filename_QuadMom_2nd, "w" );
             fprintf( file_QuadMom_2nd, "#%19s %12s %16s %16s %16s %16s %16s %16s\n",
                                        "Time", "Step", "xx", "xy", "xz", "yy", "yz", "zz" );
             fclose( file_QuadMom_2nd );
         }

         FirstTime = false;
      }

      FILE *file_QuadMom_2nd = fopen( filename_QuadMom_2nd, "a" );

                                    fprintf( file_QuadMom_2nd, "%20.14e %12ld", Time[0] * UNIT_T, Step );
      for (int b=0; b<NData; b++)   fprintf( file_QuadMom_2nd, "%17.7e",        QuadMom_2nd[b]         );
                                    fprintf( file_QuadMom_2nd, "\n"                                    );

      fclose( file_QuadMom_2nd );

   } // if ( MPI_Rank == 0 )


// free memory
                                         Phi_eff       ->FreeMemory();
   for (int NProf=0; NProf<4; NProf++)   ProfAve[NProf]->FreeMemory();

#  endif // if ( defined GRAVITY  &&  defined GREP )

} // FUNCTION : Record_GWSignal_2nd()



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_User_PostBounce
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  1. Invoked by "Flag_Check" using the function pointer "Flag_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//                Threshold   : User-provided threshold for the flag operation, which is loaded from the
//                              file "Input__Flag_User"
//                              In order of radius_min, radius_max, threshold_dens
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_User_PostBounce( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   bool Flag = false;

   const double dh        = amr->dh[lv];                                                  // grid size
   const double Center[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double Pos   [3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                              amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double r = SQRT( SQR( Center[0] - Pos[0] )
                        + SQR( Center[1] - Pos[1] )
                        + SQR( Center[2] - Pos[2] ) );

   const real (*Rho )[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];  // density
   const real dens = Rho[k][j][i];

   if ( ( r > Threshold[0] )  &&  ( r < Threshold[1])  &&  ( dens > Threshold[2] ) )
      Flag = true;

   return Flag;

} // FUNCTION : Flag_User_PostBounce

