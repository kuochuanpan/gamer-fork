#include "GAMER.h"
#include "TestProb.h"
#include "NuclearEos.h"


// problem-specific global variables
// =======================================================================================
static bool   var_bool;
static double var_double;
static int    var_int;
static char   progenitor_file[MAX_STRING]; // The supernova progenitor file
static double rot_omega;                   // Initial rotational speed                  [rad/sec]
static double rot_A;                       // Scale factor of the rotation              [km]
static double Bfield_Ab;                  // magnetic field strength                          [1e15]
static double Bfield_np;                  // dependence on the density                        [0.0]

// Parameters for the progenitor model
static double *Progenitor_Prof = NULL; // radial progenitor model
static int    Progenitor_NBin ;        // number of radial bins in the progenitor model
// =======================================================================================

static bool Flag_CoreCollapse(const int i, const int j, const int k, const int lv, const int PID, const double *Threshold);


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

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "   Setting runtime parameters ...\n" );


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
   ReadPara->Add( "rot_omega",        &rot_omega,              0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "rot_A",            &rot_A,                  1000.0,        0.0,              NoMax_double      );
#  ifdef MHD
   ReadPara->Add( "Bfield_Ab",         &Bfield_Ab,             1.0e15,        0.0,              NoMax_double      );
   ReadPara->Add( "Bfield_np",         &Bfield_np,             0.0,           NoMin_double,     NoMax_double      );
#  endif 
   //ReadPara->Add( "bounce",            &bounce,                false,         Useless_bool,     Useless_bool      );
   //ReadPara->Add( "bounce_time",       &bounce_time,           0.0,           Eps_double,       NoMax_double      );

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
      Aux_Message( stdout, "  rot_omega                 = %13.7e\n",  rot_omega );
      Aux_Message( stdout, "  rot_A                     = %13.7e\n",  rot_A ); 
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

   // Add rotation
   const double vphi = rot_omega/(1.0+SQR(r/rot_A))*r_xy;
   velx = velx - sign*sin(angle)*vphi;
   vely = vely + sign*cos(angle)*vphi;


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
#  if ( EOS == NUCLEAR )
   fluid[YE]   = ye*dens;    // electron fraction [dens]
   fluid[ENTR] = xent*dens;  // entropy [kB/baryon * dens]
#  endif
   //fluid[ENGY] = pres / ( GAMMA - 1.0 ) + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   fluid[ENGY] = (dens/(UNIT_V*UNIT_V))*(xenr + energy_shift) + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

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
   const double *Table_R      = Progenitor_Prof + 0*Progenitor_NBin;
   const double *Table_Dens   = Progenitor_Prof + 1*Progenitor_NBin;
   const double *Table_Pres   = Progenitor_Prof + 3*Progenitor_NBin;

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

   dens    = Mis_InterpolateFromTable( Progenitor_NBin, Table_R, Table_Dens, r    ); // code unit
   dens_xp = Mis_InterpolateFromTable( Progenitor_NBin, Table_R, Table_Dens, r_xp );
   dens_yp = Mis_InterpolateFromTable( Progenitor_NBin, Table_R, Table_Dens, r_yp );
   dens_zp = Mis_InterpolateFromTable( Progenitor_NBin, Table_R, Table_Dens, r_zp );

   pres    = Mis_InterpolateFromTable( Progenitor_NBin, Table_R, Table_Pres, r    ); // code unit
   pres_xp = Mis_InterpolateFromTable( Progenitor_NBin, Table_R, Table_Pres, r_xp );
   pres_yp = Mis_InterpolateFromTable( Progenitor_NBin, Table_R, Table_Pres, r_yp );
   pres_zp = Mis_InterpolateFromTable( Progenitor_NBin, Table_R, Table_Pres, r_zp );

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


void Record_CoreCollapse()
{
    //printf("debug: bounce %d \n",EOS_POSTBOUNCE);
    // Detect Core Boubce

#   if ( EOS == NUCLEAR )
    if (!EOS_POSTBOUNCE) {

        const int CountMPI = 2;
        const double bounceDens = 2.0e14; // cgs
        const double shockEntr  = 3.0;    // [kB/by]

        double dens, max_dens_loc=-__DBL_MAX__;
        double entr, max_entr_loc=-__DBL_MAX__;
        double send[CountMPI], recv[CountMPI];
        double max_dens, max_entr;
        double radius, xc, yc, zc;

        const double  BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

        // collect local data
        for (int lv=0; lv<NLEVEL; lv++)
        for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
        {
            // skip non-leaf patches
            if ( amr->patch[0][lv][PID]->son != -1 )  continue;


            for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
            for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
            for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

                xc = x - BoxCenter[0];
                yc = y - BoxCenter[1];
                zc = z - BoxCenter[2];

                radius = sqrt(xc*xc + yc*yc + zc*zc);
                radius = radius * UNIT_L; // [cm]

                dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
                entr = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENTR][k][j][i];
                entr = entr/dens;     // [kb/by]
                dens = dens * UNIT_D; // [g/cm^3]

                if ( radius <= 3.e6)
                {
                    if ( dens > max_dens_loc)
                    {
                        max_dens_loc  = dens;
                    }
                    if ( entr > max_entr_loc)
                    {
                        max_entr_loc  = entr;
                    }
                }
            }}}
        }

        // gather data

        send[0] = max_dens_loc;
        send[1] = max_entr_loc;
        MPI_Allreduce( send, recv, CountMPI, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        max_dens = recv[0];
        max_entr = recv[1];

        //printf("debug: max dens %13.7e, max entr %13.7e, time %13.7e \n",max_dens, max_entr,Time[0]);
        if (max_dens > bounceDens && max_entr > shockEntr)
        {
            if (MPI_Rank ==0)
            {
                printf("Bounce! time = %13.7e \n", Time[0]);
                EOS_POSTBOUNCE = true;
                EOS_BOUNCETIME = Time[0];
            }
        }
    }

#   endif // if ( EOS == NUCLEAR )

    // =====================================================================================
    // other stuff here
    // Record time dependent quanitites here
    //
    // Ex.  central density, shock radius, pns radius, neutirno luminosity, gw signals ....
    //
    const int CountMPI = 4;

    const char filename_central_dens[] = "Record__CentralDens";

    double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3];
    double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];

    const double  BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
    double radius, xc, yc, zc;

    for (int lv=0; lv<NLEVEL; lv++)
    for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
    {
        // skip non-leaf patches
        if ( amr->patch[0][lv][PID]->son != -1 )  continue;


        for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
        for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
        for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

            xc = x - BoxCenter[0]; // x-distance to the center
            yc = y - BoxCenter[1];
            zc = z - BoxCenter[2];

            dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
            dens = dens * UNIT_D; // [g/cm^3]

            if ( dens > max_dens_loc)
            {
                max_dens_loc  = dens;
                max_dens_pos_loc[0] = xc;
                max_dens_pos_loc[1] = yc;
                max_dens_pos_loc[2] = zc;
            }
        }}}
    }

    // gather data

    send[0] = max_dens_loc;
    send[1] = max_dens_pos_loc[0];
    send[2] = max_dens_pos_loc[1];
    send[3] = max_dens_pos_loc[2];

    MPI_Gather( send, CountMPI, MPI_DOUBLE, recv[0], CountMPI, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    double max_dens      = -__DBL_MAX__;
    int    max_dens_rank = -1;

    if ( MPI_Rank == 0 )
    {
        for (int r=0; r<MPI_NRank; r++)
        {
            if ( recv[r][0] > max_dens )
            {
                max_dens      = recv[r][0];
                max_dens_rank = r;
            }
        }

        if ( max_dens_rank < 0  ||  max_dens_rank >= MPI_NRank )
        {
            Aux_Error( ERROR_INFO, "incorrect max_dens_rank (%d) !!\n", max_dens_rank );
        }

        static bool FirstTime = true;

        if ( FirstTime )
        {
            // Central Density
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
        fprintf( file_max_dens, "%20.14e   %10ld   %14.7e   %14.7e   %14.7e   %14.7e\n",
                               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][1], recv[max_dens_rank][2], recv[max_dens_rank][3] );
        fclose( file_max_dens );

    } // if ( MPI_Rank == 0 )


}  // FUNCTION Record CoreCollapse




//---------------------------------------------------------------------------------------------------------
// Function   : Flag_CoreCollapse
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
//---------------------------------------------------------------------------------------------------------
bool Flag_CoreCollapse( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{
    bool Flag = false;

    const double dh        = amr->dh[lv];
    const double Center[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
    const double Pos   [3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                               amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                               amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh };

    const double r = SQRT( SQR(Center[0] - Pos[0])
                          +SQR(Center[1] - Pos[1])
                          +SQR(Center[2] - Pos[2]));

    const real (*Rho)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS]; // Density
    const real (*Entropy)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENTR]; // Entropy

    real dens = Rho[k][j][i];     // code unit
    real entr = Entropy[k][j][i]; // code unit

    entr = entr/dens;     // [kb/by]
    dens = dens * UNIT_D; // [g/cm^3]

    //printf("KC debug: user flags %14.7e  %14.7e %14.7e \n",  Threshold[0],Threshold[1],Threshold[2]);
    //printf("KC debug: dens entr [cgs] %14.7e  %14.7e\n", dens, entr);
 
    if (!EOS_POSTBOUNCE) {
        // collapse
        if ( (r > Threshold[0]) &&  ( r < Threshold[1] ) && (dens > Threshold[2]) )
           Flag = true;

    } else {
        // postbounce
        if ( (r > Threshold[0]) &&  ( r < Threshold[1] ) && (dens > Threshold[2]) )
           Flag = true;
        // TODO: extra refinement in shocked region
    }
    return Flag;

} // FUNCTION : Flag_CoreCollapse


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
   Init_Function_User_Ptr         = SetGridIC;
# ifdef MHD
   Init_Function_BField_User_Ptr  = SetBFieldIC;
# endif
   Init_Field_User_Ptr         = NULL;    // set NCOMP_PASSIVE_USER;        example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
   Flag_User_Ptr               = Flag_CoreCollapse;    // option: OPT__FLAG_USER;        example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr    = NULL;    // option: OPT__DT_USER;          example: Miscellaneous/Mis_GetTimeStep_User.cpp
   BC_User_Ptr                 = NULL;    // option: OPT__BC_FLU_*=4;       example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr    = NULL;    // option: OPT__RESET_FLUID;      example: Fluid/Flu_ResetByUser.cpp
   Output_User_Ptr             = NULL;    // option: OPT__OUTPUT_USER;      example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr         = Record_CoreCollapse; // option: OPT__RECORD_USER;      example: Auxiliary/Aux_Record_User.cpp
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
