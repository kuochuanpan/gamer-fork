#include "GAMER.h"
#include "TestProb.h"
#include "NuclearEos.h"


// problem-specific global variables
// =======================================================================================
static bool   var_bool;
static double var_double;
static int    var_int;
static char   progenitor_file[MAX_STRING]; // The supernova progenitor file
static bool   bounce;                      // Core bounce
static double bounce_time;                 // Core bounce time


// Parameters for the progenitor model
static double *Progenitor_Prof = NULL; // radial progenitor model
static int    Progenitor_NBin ;        // number of radial bins in the progenitor model
// =======================================================================================




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

// errors
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif


// TODO: must enable nuclear EoS


/*
#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif
*/

#  ifdef GRAVITY
   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC  ||  OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC for this test !!\n" );
#  endif

/*
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
   ReadPara->Add( "var_bool",          &var_bool,              true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "var_double",        &var_double,            1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "var_int",           &var_int,               2,             0,                5                 );
   ReadPara->Add( "progenitor_file",   progenitor_file,        Useless_str,   Useless_str,      Useless_str       );

   ReadPara->Read( FileName );

   delete ReadPara;

// Load the progenitor model file

   if ( OPT__INIT != INIT_BY_RESTART )
   {
     const bool RowMajor_No  = false; // load data into the column major OPT__RECORD_USER
     const bool AllocMem_Yes = true; // allocate memort for Progenitor_Prof
     const int  NCol         = 6;    // total number of columns to load   TODO: use runtime paprameter
     const int  TargetCols[NCol] = {0,2,3,4,5,7}; // target columns: {radius, density, temp, pressure, velr, ye}

     double *Table_R, *Table_Dens, *Table_Temp, *Table_Pres, *Table_Velr, *Table_Ye;

     Progenitor_NBin = Aux_LoadTable(Progenitor_Prof, progenitor_file, NCol, TargetCols, RowMajor_No, AllocMem_Yes);

     Table_R    = Progenitor_Prof + 0*Progenitor_NBin;
     Table_Dens = Progenitor_Prof + 1*Progenitor_NBin;
     Table_Temp = Progenitor_Prof + 2*Progenitor_NBin;
     Table_Pres = Progenitor_Prof + 3*Progenitor_NBin;
     Table_Velr = Progenitor_Prof + 4*Progenitor_NBin;
     Table_Ye   = Progenitor_Prof + 5*Progenitor_NBin;

     // convert to code units (assuming progentior model is in cgs)
     for (int b=0; b<Progenitor_NBin; b++)
     {
       Table_R[b]    /= UNIT_L;
       Table_Dens[b] /= UNIT_D;
       Table_Pres[b] /= UNIT_P;
       Table_Velr[b] /= UNIT_V;
     }

   }

// load the EOS table
//nuc_eos_C_ReadTable("LS220.h5");
//printf("nrho: %d\n",nrho);
//printf("ntemp: %d\n",ntemp);
//printf("nye: %d\n",nye);

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
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      //Aux_Message( stdout, "  var_bool                  = %d\n",     var_bool );
      //Aux_Message( stdout, "  var_double                = %13.7e\n", var_double );
      //Aux_Message( stdout, "  var_int                   = %d\n",     var_int );
      Aux_Message( stdout, "  progenitor_file           = %s\n",     progenitor_file );
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
   const double *Table_R    = Progenitor_Prof + 0*Progenitor_NBin;
   const double *Table_Dens = Progenitor_Prof + 1*Progenitor_NBin;
   const double *Table_Temp = Progenitor_Prof + 2*Progenitor_NBin;
   const double *Table_Pres = Progenitor_Prof + 3*Progenitor_NBin;
   const double *Table_Velr = Progenitor_Prof + 4*Progenitor_NBin;
   const double *Table_Ye   = Progenitor_Prof + 5*Progenitor_NBin;

   double xc, yc, zc;
   double r, dens, temp, pres, velr, velx, vely, velz, ye, r_xy, v_xy, angle, sign;

   const double temp_mev_to_kelvin = 1.1604447522806e10;
   double xtmp, xenr, xprs, xent, xcs2, xdedt, xdpderho, xdpdrhoe, xmunu;
   int keyerr;
   const double rfeps = 1.0e-10;

   xc = x - BoxCenter[0];
   yc = y - BoxCenter[1];
   zc = z - BoxCenter[2];

   r = sqrt(SQR(xc) + SQR(yc) + SQR(zc));

   dens = Mis_InterpolateFromTable(Progenitor_NBin, Table_R, Table_Dens, r); // code unit
   temp = Mis_InterpolateFromTable(Progenitor_NBin, Table_R, Table_Temp, r); // [K]
   pres = Mis_InterpolateFromTable(Progenitor_NBin, Table_R, Table_Pres, r);
   velr = Mis_InterpolateFromTable(Progenitor_NBin, Table_R, Table_Velr, r); // code unit
   ye   = Mis_InterpolateFromTable(Progenitor_NBin, Table_R, Table_Ye, r);

   xtmp = temp/temp_mev_to_kelvin; // to MeV

   r_xy = sqrt( SQR(xc) + SQR(yc));

   if (r_xy == 0.0)
   {
     angle = M_PI/2.;
   } else
   {
     angle = atan(zc/r_xy);
   }
   velz = velr*sin(angle);

   v_xy = velr*cos(angle);
   if (xc == 0.0)
   {
     angle = M_PI/2.0;
   } else {
     angle = atan(yc/xc);
   }
   sign = xc/ abs(xc);
   velx = sign*v_xy*cos(angle); // code unit
   vely = sign*v_xy*sin(angle);

   // call EOS to get other variables
   // use temperature mode (1)

   nuc_eos_C_short((dens*UNIT_D),&xtmp,ye,&xenr, &xprs, &xent, &xcs2, &xdedt, &xdpderho,
        &xdpdrhoe, &xmunu, 1, &keyerr, rfeps);

   if (keyerr != 0)
   {
     printf("debug: keyerr not zero %d\n",keyerr);
   }

   //printf("debug: gamc %15.6E\n",(xcs2*(dens*UNIT_D)/xprs));
   //printf("debug: game %15.6E\n",(xprs/((dens*UNIT_D)*xenr)+1.0));

   fluid[DENS] = dens; // code unit
   fluid[MOMX] = dens*velx;
   fluid[MOMY] = dens*vely;
   fluid[MOMZ] = dens*velz;
   fluid[YE]   = ye*dens;    // electron fraction [dens]
   fluid[ENTR] = xent*dens;  // entropy [kB/baryon * dens]
   //fluid[ENGY] = pres / ( GAMMA - 1.0 ) + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   fluid[ENGY] = (dens/(UNIT_V*UNIT_V))*(xenr + energy_shift) + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

} // FUNCTION : SetGridIC



void End_CoreCollapse()
{
  delete [] Progenitor_Prof;
} // FUNCTION End_CoreCollapse

#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Template
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CoreCollapse()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr      = SetGridIC;
   Init_Field_User_Ptr         = NULL;    // set NCOMP_PASSIVE_USER;        example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
   Flag_User_Ptr               = NULL;    // option: OPT__FLAG_USER;        example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr    = NULL;    // option: OPT__DT_USER;          example: Miscellaneous/Mis_GetTimeStep_User.cpp
   BC_User_Ptr                 = NULL;    // option: OPT__BC_FLU_*=4;       example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr    = NULL;    // option: OPT__RESET_FLUID;      example: Fluid/Flu_ResetByUser.cpp
   Output_User_Ptr             = NULL;    // option: OPT__OUTPUT_USER;      example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr         = NULL;    // option: OPT__RECORD_USER;      example: Auxiliary/Aux_Record_User.cpp
   Src_User_Ptr                = NULL;       //
   End_User_Ptr                = End_CoreCollapse;    // option: none;                  example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr        = NULL;    // option: OPT__GRAVITY_TYPE=2/3; example: SelfGravity/Init_ExternalAcc.cpp
   Init_ExternalPot_Ptr        = NULL;    // option: OPT__EXTERNAL_POT;     example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> Init_ExtPot()
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr     = NULL;    // option: PAR_INIT=1;            example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr = NULL;    // set PAR_NATT_USER;             example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Template
