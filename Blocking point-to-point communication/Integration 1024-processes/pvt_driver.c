/*=============================================================================
  |  SYM.8 : pvt_driver.c
  |----------------------------------------------------------------------------
  |  Copyright 2004-2011 (c) Sienna Geodynamics & Consulting, Inc.
  |  author: Anthony Park, ajpark@sienna-geo.com
  |--------------------------------------------------------------------------*/
#include "h_parameters.h"
#include "h_define.h"
#include "h_global.h"

#include "h_balance.h"
#include "h_boundary.h"
#include "h_grid.h"
#include "h_numeric.h"
#include "h_recover.h"
#include "h_mpi.h"
#include "FEM.h"

extern double **allocate_double_2d ();
extern double ***allocate_double_3d ();

void         allocate_fixed_pv ();
enum status  PV_driver ();
void         PV_fdsolver_coeffs_1 ();
void         PV_fdsolver_coeffs_2 ();
void         PV_fdsolver_coeffs_3 ();
enum logical PV_perm_change ();
enum status  T_driver ();
void         T_constant_initial ();
double       T_geotherm_cell ();
void         T_geotherm_basin ();
void         T_geotherm_reservoir ();

/*===========================================================================*/
/*  primary driver for T solvers.  */

enum status T_driver (double dt)

{
  int ix, iy, iz;

  /*---------------------------------------------------------------------------
    |  solve or setup following properties in the order of
    |    domain configuration, and solver type.
    |
    |  IMPORTANT : for BCM simulations T is set in the
    |    boundary conditions handling function in boundary_set.c,
    |    update_boundary_conditions_bcm ().
    |------------------------------------------------------------------------*/

  /*  -----------------------------  invoke a temperature setting function.  */
  if (OPT_SIMULATOR == BCM)
    {
      //  bcm_bc_T is set using the database queried T setting.
      for_all_xyz
        cell[each]->pvt->TC = bcm_bc_T[iz];
    }
  else if (OPT_T == 0)  //  constant T from input file.
    T_constant_initial ();
  //  else if (OPT_SIMTYPE == RESERVOIR1D)
  else if (OPT_T == 1)  //  set T according to geothermal gradient.
    T_geotherm_reservoir ();
  else if (OPT_T == 2)  //  set T according to HKF model
    ;
  else if (OPT_SIMTYPE == BASIN)
    T_geotherm_basin ();
  else
    T_geotherm_reservoir ();

  //  ------------------------------------------------  set kelvin temperature.
  if (ntimesteps <= 1){
    for_all_xyz{
      cell[each]->pvt->TK = cell[each]->pvt->TC + 273.15;
      cell[each]->pvt->HKF_TK = cell[each]->pvt->TK;

      if(isnan(cell[each]->pvt->HKF_TK)) {
        printf("T_driver: nan: HKF_TK cell %d\n",ix);
        exit(1);
      }

    }

  }
  return (none);
}

/*===========================================================================*/
/*  use fixed initial background temperature.  */

void T_constant_initial ()

{
  int ix, iy, iz;
  double TC;
  /*---------------------------------------------------------------------------
    |  set temperature field using the initial background temperature.
    |------------------------------------------------------------------------*/
#pragma omp parallel for private (ix, iy, iz)
  for_all_xyz
    cell[each]->pvt->TC = temperature_initial;;

}
/*===========================================================================*/
/*  temperature solver reservoir option.  */

void T_geotherm_reservoir ()

{
  int ix, iy, iz;
  double TC;
  /*---------------------------------------------------------------------------
    |  set temperature field according to geothermal gradient data.
    |  the z coordinate is the depth from the surface.
    |
    |  09.07.07.  compute the temperature for the bottom-most unit, and
    |  set the rest of the reservoir to the same temperature.
    |--------------------------------------------------------------------------
    |  11.03.23.  set T according to geothermal gradient.
    |------------------------------------------------------------------------*/
  /*------------------------  this segment as of 09.07.07.  replaced, 11.03.23.
    TC = T_geotherm_cell (temperature_surface, -cell[1][1][1]->z,
    subs_use[index_subs_now]->geotherm);

    for_all_xyz
    cell[each]->pvt->TC = TC;
    ---------------------------------------------------------------------------*/

  //  following segment as of 11.03.23.
#pragma omp parallel for private (ix, iy, iz)
  for_all_xyz
    cell[each]->pvt->TC =
    T_geotherm_cell (temperature_surface, -cell[each]->z,
                     subs_use[index_subs_now]->geotherm);

}
/*===========================================================================*/
/*  temperature solver option for basin configuration.  */

void T_geotherm_basin ()

{
  int ix, iy, iz;
  double ztop = cell[1][1][nz]->z01;
  double z;
  /*---------------------------------------------------------------------------
    |  set temperature field according to geothermal gradient data.
    |  compute z, depth, as the difference between the top cell boundary
    |  and the center of the cell.
    |------------------------------------------------------------------------*/
#pragma omp parallel for private (ix, iy, iz, z)
  for_all_xyz
    {
      z = ztop - cell[each]->z;
      cell[each]->pvt->TC = T_geotherm_cell
        (temperature_surface, z, subs_use[index_subs_now]->geotherm);
    }  /*  end, for xyz.  */
}
/*===========================================================================*/
/*  temperature in celsius; geothermal gradient.  */

double T_geotherm_cell (double Tsurface, double z, double geotherm)

{
  double TC;

  /*---------------------------------------------------------------------------
    |  geotherm is in deg/cm, positive value for increasing T with z;
    |  z is in cm, positive with increasing depth.
    |------------------------------------------------------------------------*/
  TC = Tsurface + (z * geotherm);

  return (TC);
}
/*===========================================================================*/
/*  primary driver for P and V solvers.  */

enum status PV_driver (double dt, enum status stat_rxn)

{
  int ix=1, iy=1, iz=1;
  double var, v, weight, t1, t2;
  enum logical pvt_flag = true;

  /*---------------------------------------------------------------------------
    |  solve or setup following properties in the order shown:
    |    stress; and
    |    pressure fluid flow velocities.
    |--------------------------------------------------------------------------
    |  stress is not used, and not calculated.
    |
    |  take inlet and outlet boundary fluid pressure conditions from the
    |  burial history data set. take a linear average between two data
    |  intervals.
    |------------------------------------------------------------------------*/
  if (index_subs_now < nsubs)
    {
      //  take linearly interpolated pf values using time and subs segments.
      const_bc_p_inlet = subs_use[index_subs_now]->pf_inlet;
      const_bc_p_outlet = subs_use[index_subs_now]->pf_outlet;

      t1 = subs_use[index_subs_now]->time_seconds;
      t2 = subs_use[index_subs_now + 1]->time_seconds;
      weight = (time_seconds - t1) / (t2 - t1);

      const_bc_p_inlet = subs_use[index_subs_now]->pf_inlet +
        weight * (subs_use[index_subs_now + 1]->pf_inlet -
                  subs_use[index_subs_now]->pf_inlet);
      const_bc_p_outlet = subs_use[index_subs_now]->pf_outlet +
        weight * (subs_use[index_subs_now + 1]->pf_outlet -
                  subs_use[index_subs_now]->pf_outlet);

    }
  else
    {
      const_bc_p_inlet = subs_use[index_subs_now]->pf_inlet;
      const_bc_p_outlet = subs_use[index_subs_now]->pf_outlet;
    }
  /*
  if ( OPT_FEM == 1 )
    {
      if (OPT_SIMTYPE==RESERVOIR1D)
        {
          long ii;
          // Update inlet nodes, and populate dirichlet BCs for FEM_Grid
          for (ii=0; ii<FEM_Grid->n_bc_d; ii++)
            if (FEM_Grid->bc_d_identifier[ii]=='i')
              {
                FEM_Grid->bc_d_values[ii] = const_bc_p_inlet;
              }
            else if (FEM_Grid->bc_d_identifier[ii] == 'o')
              FEM_Grid->bc_d_values[ii] = const_bc_p_outlet;
        }
      else if (OPT_SIMTYPE==RESERVOIR2D)
        {
          // Handle boundary conditions in libMesh code.
        }
    }
  */
  /*---------------------------------------------------------------------------
    |  stat_rxn shortcut (2011.10.20) :
    |
    |    stat_rxn is set to error when the reaction solver fails and
    |    forces reiteration. when that happens, reset Pf to previously
    |    calculated values rather than re-solving. previous results
    |    are saved to PF_new.
    |
    |    velocities vx, vy, and vz are also saved to .._new.
    |------------------------------------------------------------------------*/
  if (stat_rxn == error)
    {
      for_all_xyz
        {
          cell[each]->pvt->Pf = cell[each]->pvt->Pf_new;
          cell[each]->pvt->vx = cell[each]->pvt->vx_new;
          cell[each]->pvt->vy = cell[each]->pvt->vy_new;
          cell[each]->pvt->vz = cell[each]->pvt->vz_new;
          /* fprintf(stdout,"error: reset to (p,v)=(%.2e,<%2.e,%.2e,%.2e>)\n", */
          /*      cell[each]->pvt->Pf, */
          /*      cell[each]->pvt->vx, */
          /*      cell[each]->pvt->vy, */
          /*      cell[each]->pvt->vz); */

        }  //  end, for xyz.
    }
  else if (OPT_SIMTYPE == RESERVOIR3D)
    {
      if (OPT_SIMULATOR == SYM8)
        {
#ifdef _MPI_
          mpi_broadcast_sym8_timestamp();
          mpi_send_sym8_cell_data();
          /*
          if( mpi_my_rank == MASTER )
            {
              int message = 0;
              int datasize = nx * ny * nz;
              double data[datasize];
              int k = 0;
              for (ix=1; ix<=nx; ix++)
                for (iy=1; iy<=ny; iy++)
                  for (iz=1; iz<=nz; iz++)
                    {
                      data[k++] = cell[ix][iy][iz]->sed->porosity_initial;
                    }

              // works:
              // MPI_Send(&data, datasize, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

              int worker;
              for( worker = 1; worker < mpi_np; worker++ )
                {
                  fprintf(stderr,"Rank %d process sending porosity field (%d cells) to worker %d\n",
                          mpi_my_rank,datasize,worker);
                  MPI_Send(&data, datasize, MPI_DOUBLE, worker, 1, MPI_COMM_WORLD);
                }
              */
              // MPI_Bcast does not work
              /* if( MPI_Bcast(&data, datasize, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS ) */
              /*        { */
              /*          fprintf(stderr,"Rank %d process failed sending porosity field to all workers\n", */
              /*                  mpi_my_rank); */
              /*        } */
              /* else */
              /*        { */
              /*          fprintf(stderr,"Rank %d process finished broadcaasting porosity field to all workers\n", */
              /*                  mpi_my_rank); */
              /*        } */
          //} // if master
#endif
          /* In Sym8, use libMesh for pressure and flux/velocity calculations */
        if( debug )
          fprintf(stderr, "Rank %d, starting solve (t=%.2fs, dt=%.4gs, ns=%d)\n",
                  mpi_my_rank, time_seconds, timestep, ntimesteps);
        libMesh_solve(time_seconds, timestep, ntimesteps, send_data);
        if( debug )
          fprintf(stderr, "Rank %d, completed solve\n", mpi_my_rank);
#ifdef _MPI_
          int error_code;

          //
          /* MPI_Gather(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count,  */
          if( debug )
            fprintf(stderr,"Rank %d process gather blocks of %d cells from all workers, %d total cells\n",
                    mpi_my_rank,send_count,recv_count);

          // recv_count is the count of elements received per process, not the total summation of
          // counts from all processes.
          //
          error_code = MPI_Gather(send_data, send_count, mpi_cell_data, recv_data, send_count, mpi_cell_data, 0, MPI_COMM_WORLD);

          if (error_code != MPI_SUCCESS)
            {
              char error_string[BUFSIZ];
              int length_of_error_string;

              MPI_Error_string(error_code, error_string, &length_of_error_string);
              fprintf(stderr, "%3d: %s\n", mpi_my_rank, error_string);
            }
          if( debug )
            fprintf(stderr,"Rank %d process end gather.  Sent %d structs totalling %d bytes. (%d=%d?)\n",
                    mpi_my_rank,
                    send_count,
                    send_count * sizeof(mpi_cell_data),
                    recv_count * sizeof(cell_data),
                    recv_count * sizeof(mpi_cell_data));

          int cellDataIndex;
          pt_cell *acell;
          int offset;
          int process;

          for(process = 0; process < mpi_np; process++)
            {
              offset = process * send_count;

              if( debug )
                fprintf(stderr,"Rank %d process spans recv_data[%d - %d], recv_data length = %d\n",
                        process,
                        offset,
                        (process + 1) * send_count * sizeof(cell_data) - 1,
                        recv_count * sizeof(cell_data));

              for(cellDataIndex=0; cellDataIndex<send_count; cellDataIndex++)
                {

                  if( recv_data[offset + cellDataIndex].ix == 0 || recv_data[offset + cellDataIndex].iy == 0 || recv_data[offset + cellDataIndex].iz == 0 ) {
                    fprintf(stderr,"pvt_driver.c: zero cell index (%d,%d,%d)\n",
                            recv_data[offset + cellDataIndex].ix,
                            recv_data[offset + cellDataIndex].iy,
                            recv_data[offset + cellDataIndex].iz);
                  }

                  acell = cell[recv_data[offset + cellDataIndex].ix][recv_data[offset + cellDataIndex].iy][recv_data[offset + cellDataIndex].iz];

                  // MPI process rank
                  //
                  acell->rank = recv_data[offset + cellDataIndex].rank;

                  /*
                  fprintf(stderr,"*** Update from rank %d, cell %d (%d,%d,%d), address=%0x, offset=%d, after MPI_Gather\n",
                          acell->rank,
                          cellDataIndex,
                          recv_data[offset + cellDataIndex].ix,
                          recv_data[offset + cellDataIndex].iy,
                          recv_data[offset + cellDataIndex].iz,
                          &recv_data[offset + cellDataIndex],
                          offset);
                  */

                  // pressure and temperature
                  //
                  acell->pvt->Pf = recv_data[offset + cellDataIndex].Pf;
                  acell->pvt->Pf_new = recv_data[offset + cellDataIndex].Pf;

                  // update cell temperature from libMesh solution
                  //
                  acell->pvt->TC_old = acell->pvt->TC;
                  acell->pvt->TC = recv_data[offset + cellDataIndex].TC;
                  acell->pvt->TK = acell->pvt->TC + 273.15;

                  // fluid velocity
                  //
                  acell->pvt->vx_old = acell->pvt->vx;
                  acell->pvt->vy_old = acell->pvt->vy;
                  acell->pvt->vz_old = acell->pvt->vz;
                  acell->pvt->vx = recv_data[offset + cellDataIndex].vx;
                  acell->pvt->vy = recv_data[offset + cellDataIndex].vy;
                  acell->pvt->vz = recv_data[offset + cellDataIndex].vz;
                  acell->pvt->vx_new = acell->pvt->vx;
                  acell->pvt->vy_new = acell->pvt->vy;
                  acell->pvt->vz_new = acell->pvt->vz;

                  // strain, stress, and displacement
                  //
                  acell->strain->xx = recv_data[offset + cellDataIndex].strain_xx;
                  acell->strain->xy = recv_data[offset + cellDataIndex].strain_xy;
                  acell->strain->xz = recv_data[offset + cellDataIndex].strain_xz;
                  acell->strain->yy = recv_data[offset + cellDataIndex].strain_yy;
                  acell->strain->yz = recv_data[offset + cellDataIndex].strain_yz;
                  acell->strain->zz = recv_data[offset + cellDataIndex].strain_zz;

                  acell->stress->xx = recv_data[offset + cellDataIndex].stress_xx;
                  acell->stress->xy = recv_data[offset + cellDataIndex].stress_xy;
                  acell->stress->xz = recv_data[offset + cellDataIndex].stress_xz;
                  acell->stress->yy = recv_data[offset + cellDataIndex].stress_yy;
                  acell->stress->yz = recv_data[offset + cellDataIndex].stress_yz;
                  acell->stress->zz = recv_data[offset + cellDataIndex].stress_zz;

                  acell->displacement_x = recv_data[offset + cellDataIndex].displacement_x;
                  acell->displacement_y = recv_data[offset + cellDataIndex].displacement_y;
                  acell->displacement_z = recv_data[offset + cellDataIndex].displacement_z;

                  // porosity, tortuosity, and permeability
                  //
                  acell->sed->porosity_old = acell->sed->porosity;
                  acell->volume->cell_old = acell->volume->cell;
                  acell->volume->pore_old = acell->volume->pore;
                  acell->volume->matrix_old = acell->volume->matrix;
                  acell->volume->cell = recv_data[offset + cellDataIndex].cell;
                  acell->volume->pore = recv_data[offset + cellDataIndex].pore;

                  update_volume_fractions(acell);
                  sediment_porosity(acell);
                  sediment_tortuosity(acell);
                  sediment_permeability(acell);

                  // fprintf(stderr,"*** cell %d updated\n", cellDataIndex);
                }
              if( debug )
                fprintf(stderr,"Rank %d process finished updating %d cells after MPI_Gather\n",
                        process,cellDataIndex);
            }
          if( debug )
            fprintf(stderr,"Rank %d process finished updating all %d cells after MPI_Gather\n",
                    mpi_my_rank,send_count*mpi_np);
#endif
        }
      else {
        /*  -----------------  configure inlet and outlet boundary P b.c. values.
            this is set only once at the beginning of the simulation.  */
        if ((ntimesteps <= 1) || (recover_status == true))
          {
            PV_P3D_inlet_boundary ();
            PV_P3D_outlet_boundary ();
          }  //  end, if ntimesteps.

        if (PV_perm_change () == true)
          {
            PV_fdsolver3D_linbcg (dt);

            PV_P3D_velocities ();
          }  //  end, if PV_perm_change.
      }
    }
  else if (OPT_SIMTYPE == RESERVOIR2D)
    {
      if (OPT_SIMULATOR == SYM8)
        {
          /* In Sym8, use libMesh for pressure and flux/velocity calculations */

#ifdef _MPI_
          int message = 0;
          int datasize = nx * ny * nz + 3;
          double data[datasize];
          int k = 0;
          for (ix=1; ix<=nx; ix++)
            for (iy=1; iy<=ny; iy++)
              for (iz=1; iz<=nz; iz++)
                {
                  /*
                    printf("MPI_Send: cell(%d,%d,%d) = data(%d) phi = %.2f %d bytes\n",
                    ix,iy,iz,k,
                    cell[ix][iy][iz]->sed->porosity_initial,
                    sizeof(cell[ix][iy][iz]));
                  */
                  data[k++] = cell[ix][iy][iz]->sed->porosity_initial;
                }
          data[k++] = time_seconds;
          data[k++] = timestep;
          data[k++] = (double)ntimesteps;
          fprintf(stderr,"Rank %d process sending porosity field, time_seconds, and timestep to all workers\n",mpi_my_rank);
          MPI_Send(&data, datasize, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
#endif
          libMesh_solve(time_seconds, timestep, ntimesteps, send_data);
        }
  else
  /*  -----------------  configure inlet and outlet boundary P b.c. values.
          this is set only once at the beginning of the simulation.  */
  {

    if ((ntimesteps <= 1) || (recover_status == true))
        {
          PV_P2D_inlet_boundary ();
          PV_P2D_outlet_boundary ();
        }  //  end, if ntimesteps.

      if (PV_perm_change () == true)
        {
          PV_fdsolver2D_linbcg (dt);

          PV_P2D_velocities ();
        }  //  end, if PV_perm_change.
  }
    }
  else if (OPT_SIMTYPE == RESERVOIR1D)
    {
      /*  -------  set water flow velocity according to constant-flux data.  */
      if (PV_perm_change () == true)
        {
          if ((OPT_P1D == 0) && (OPT_FEM==0))
            PV_flux_velocities ();  //  ----------------  constant flux option.
          else if ( (OPT_FEM==1) && (OPT_ORIENTATION == VERTICAL) )
            {
              poisson_solver_2D( FEM_Grid, dt ); //DISABLE
              for_all_z
                {
                  double    **cell_wt   = FEM_Grid->cell_weight;
                  long      nnodes      = FEM_Grid->n_nodes;

                  cell[each]->pvt->Pf_new = cell[each]->pvt->Pf = cell_avg(FEM_Grid->fluid_pressure,cell_wt,iz,nnodes);
                  cell[each]->pvt->Pover  = cell_avg(FEM_Grid->excess_pressure,cell_wt,iz,nnodes);
                  cell[each]->pvt->Phstat = cell_avg(FEM_Grid->hydrostatic_pressure,cell_wt,iz,nnodes);
                  cell[each]->pvt->breakdown_pressure_HW = cell_avg(FEM_Grid->breakdown_pressure_HW,cell_wt,iz,nnodes);
                  cell[each]->pvt->breakdown_pressure_HF = cell_avg(FEM_Grid->breakdown_pressure_HF,cell_wt,iz,nnodes);
                  *cell[each]->stress = cell_tensor_avg(FEM_Grid->stress,cell_wt,iz,nnodes);
                  *cell[each]->strain = cell_tensor_avg(FEM_Grid->strain,cell_wt,iz,nnodes);
                }
              PV_P1D_velocities_vertical();
            }
          else if ((OPT_P1D == 1) && (OPT_ORIENTATION == VERTICAL))
            {
              PV_fdsolver1D_linbcg_vertical (dt);

              PV_P1D_velocities_vertical ();
            }
          else if ((OPT_FEM==1) && (OPT_ORIENTATION == HORIZONTAL))
            {
              libMesh_solve(time_seconds, timestep, ntimesteps, send_data);

              /* poisson_solver_2D(FEM_Grid, dt); */
              /* for_all_x */
              /*        { */
              /*          double    **cell_wt   = FEM_Grid->cell_weight; */
              /*          long      nnodes      = FEM_Grid->n_nodes; */

              /*          cell[each]->pvt->Pf_new = cell[each]->pvt->Pf = cell_avg(FEM_Grid->fluid_pressure, cell_wt, ix, nnodes); */
              /*          cell[each]->pvt->Pover  = cell_avg(FEM_Grid->excess_pressure,cell_wt,ix,nnodes); */
              /*          cell[each]->pvt->Phstat = cell_avg(FEM_Grid->hydrostatic_pressure,cell_wt,ix,nnodes); */
              /*          cell[each]->pvt->breakdown_pressure_HW = cell_avg(FEM_Grid->breakdown_pressure_HW,cell_wt,ix,nnodes); */
              /*          cell[each]->pvt->breakdown_pressure_HF = cell_avg(FEM_Grid->breakdown_pressure_HF,cell_wt,ix,nnodes); */
              /*          *cell[each]->stress = cell_tensor_avg(FEM_Grid->stress,cell_wt,ix,nnodes);                  */
              /*          *cell[each]->strain = cell_tensor_avg(FEM_Grid->strain,cell_wt,ix,nnodes); */
              /*        } */
              /* PV_P1D_velocities_horizontal(); */
            }
          else if ((OPT_P1D == 1) && (OPT_ORIENTATION == HORIZONTAL))
            {
              /* Poroelastic pressure solver Jon Matthews */
              Poro_bksolver1D_horizontal(dt);

              PV_P1D_velocities_horizontal ();
            }
        }  //  end, if PV_perm_change.
    }  //  end, if else.
  else
    {
      printf ("\t#### ERR. UNKNOWN OPT_SIMTYPE, PV DRIVER.\n");
      simulation_terminate ("UNKNOWN OPT_SIMTYPE, PV DRIVER.");
    }  //  end, if OPT_SIMTYPE.

  /*  ----------------------------------------  set absolute flow velocity.  */
#pragma omp parallel for private (ix, iy, iz)
  for_all_xyz
    {
      v = (cell[each]->pvt->vx * cell[each]->pvt->vx) +
        (cell[each]->pvt->vz * cell[each]->pvt->vz);
      cell[each]->pvt->v = pow (v, 0.5);
    }  //  end, for xyz.

  // Update inlet and outlet velocity
  if (OPT_ORIENTATION==HORIZONTAL)
  {
    cell[0][1][1]->pvt->v = cell[1][1][1]->pvt->v;
    cell[ix][1][1]->pvt->v = cell[ix-1][1][1]->pvt->v;
  }
  else if (OPT_ORIENTATION==VERTICAL)
  {
    cell[1][1][0]->pvt->v = cell[1][1][1]->pvt->v;
    cell[1][1][iz]->pvt->v = cell[1][1][iz-1]->pvt->v;
  }

  /*  TEMPORARY : SPECIAL CASE, AVID 120207.
      set vx, vy, and vz with imposed pattern recorded in the grid file.  */
  if (OPT_SPCL_AVID == 1)
    for_all_xyz
      {
        if (idx_grid_vx != 0)
          cell[each]->pvt->vx = grid_var[idx_grid_vx]->var[each];

        if (idx_grid_vy != 0)
          cell[each]->pvt->vy = grid_var[idx_grid_vy]->var[each];

        if (idx_grid_vz != 0)
          cell[each]->pvt->vz = grid_var[idx_grid_vz]->var[each];
      }  //  end, for_all_xyz.

  /*  ----------------------------------------  set flow distance (v * dt).  */
  var = yrtosec (1.0);
#pragma omp parallel for private (ix, iy, iz)
  for_all_xyz
    {
      cell[each]->pvt->vx_distance = cell[each]->pvt->vx * var;
      cell[each]->pvt->vy_distance = cell[each]->pvt->vy * var;
      cell[each]->pvt->vz_distance = cell[each]->pvt->vz * var;
    }  //  end, for xyz.

  return (none);
}
/*===========================================================================*/
/*  evaluate perm change since last PV computation.  */

enum logical PV_perm_change ()

{
  enum logical flag = true;
  int ix, iy, iz;
  double diff_perm = 0.0;
  double diff = 0.0;

  if ((limit_pv_perm > 0.0) || (recover_status == true))
    {
      for_all_xyz
        {
          diff = 1.0 - (cell[each]->sed->perm_pv / cell[each]->sed->perm);
          if (diff > diff_perm) diff_perm = diff;
        }  //  end, for_all_xyz.

      if (diff_perm < limit_pv_perm) flag = false;  //  evaluate.

      if (flag == true)  //  save current perm to perm_pv.
        for_all_xyz
          cell[each]->sed->perm_pv = cell[each]->sed->perm;
    }  //  end, if limit_pv_perm.

  return (flag);
}
/*===========================================================================*/
/*  set the P equation coefficients.  */

void PV_fdsolver_coeffs_1 (pt_cell *acell, double dt)

{
  int ix = acell->ix;
  int iy = acell->iy;
  int iz = acell->iz;
  double var, perm;
  double water_density_molar =
    acell->water->density_mass / water_mol_wt;
  double water_density_molar_old =
    acell->water->density_mass_old / water_mol_wt;

  /*---------------------------------------------------------------------------
    |  P equation, incompressible porous medium.
    |
    |  coefficients are for a pressure pde of the form
    |
    |            rho(water) * perm
    |    grad [ ------------------- grad P ] -
    |                viscosity
    |
    |       [ phi * rho(water) - phi(t-dt) * rho(water, t-dt) ] / dt = 0.
    |
    |  this eqn. is derived from water conservation of mass and darcy
    |  flow equation, and by neglecting P compressibility and thermal
    |  expansivity. the gravity term is also removed.
    |------------------------------------------------------------------------*/

  pa[each] = pb[each] = pc[each] = 0.0;

  pa[each] += -acell->sed->perm * water_density_molar /
    acell->water->viscosity_kinematic;

  //  pb[each] += 0.0;

  pc[each] += ((water_density_molar * acell->sed->porosity) -
               (water_density_molar_old * acell->sed->porosity_old)) / dt;
}
/*===========================================================================*/
/*  set the P equation coefficients.  */

void PV_fdsolver_coeffs_2 (pt_cell *acell, double dt)

{
  int ix = acell->ix;
  int iy = acell->iy;
  int iz = acell->iz;
  double var, perm;
  double water_density_molar, water_density_molar_old;

  /*---------------------------------------------------------------------------
    |  P equation, incompressible porous medium.
    |  derived from equation of state for water, conservation of water,
    |  and darcy flow equation.
    |------------------------------------------------------------------------*/
  water_density_molar = acell->water->density_mass / water_mol_wt;
  water_density_molar_old = acell->water->density_mass_old / water_mol_wt;

  pa[each] = pb[each] = pc[each] = 0.0;

  pa[each] += -acell->sed->perm * water_density_molar /
    (acell->water->viscosity_kinematic * water_P_compressibility);

  pb[each] += acell->sed->porosity * water_density_mass_ref /
    (dt * water_mol_wt);

  var = water_density_mass_ref *
    (1.0 - water_T_expansivity * (acell->pvt->TC - reference_TC) -
     (water_P_compressibility * reference_P));

  pc[each] += acell->sed->porosity * var /
    (dt * water_mol_wt * water_P_compressibility);

  pc[each] += -acell->sed->porosity_old *
    water_density_molar_old / (dt * water_P_compressibility);

  var = water_density_mass_ref *
    (1.0 - water_T_expansivity * (acell->pvt->TC - reference_TC));

  pc[each] += acell->sed->porosity * var / (water_mol_wt * dt);

  pc[each] += -acell->sed->porosity_old * acell->water->density_mass_old /
    (water_mol_wt * dt);

  pa[each] *= timestep;
  pb[each] *= timestep;
  pc[each] *= timestep;

}
/*===========================================================================*/
/*  set the P equation coefficients.  */

void PV_fdsolver_coeffs_3 (pt_cell *acell, double dt)

{
  int ix = acell->ix;
  int iy = acell->iy;
  int iz = acell->iz;
  double var, perm;
  double water_density_molar, water_density_molar_old;

  /*---------------------------------------------------------------------------
    |  P equation, incompressible porous medium.
    |  derived from equation of state for water, conservation of water,
    |  and darcy flow equation (10.01.08).
    |
    |
    |------------------------------------------------------------------------*/
  water_density_molar = acell->water->density_mass / water_mol_wt;
  water_density_molar_old = acell->water->density_mass_old / water_mol_wt;

  pa[each] = pb[each] = pc[each] = 0.0;

  pa[each] += -acell->sed->perm * water_density_molar /
    (acell->water->viscosity_kinematic * water_P_compressibility);

  pb[each] += acell->sed->porosity * water_density_mass_ref /
    (dt * water_mol_wt);

  var = water_density_mass_ref *
    (1.0 - water_T_expansivity * (acell->pvt->TC - reference_TC) -
     (water_P_compressibility * reference_P));

  pc[each] += acell->sed->porosity * var /
    (dt * water_mol_wt * water_P_compressibility);

  pc[each] += -acell->sed->porosity_old *
    water_density_molar_old / (dt * water_P_compressibility);

  var = water_density_mass_ref *
    (1.0 - water_T_expansivity * (acell->pvt->TC - reference_TC));

  pc[each] += acell->sed->porosity * var / (water_mol_wt * dt);

  pc[each] += -acell->sed->porosity_old * acell->water->density_mass_old /
    (water_mol_wt * dt);

  pa[each] *= timestep;
  pb[each] *= timestep;
  pc[each] *= timestep;

}
/*===========================================================================*/
/*  allocate fixed space required for the pv solver modules.  */

void allocate_fixed_pv ()

{
  int ix, iy, iz;
  int nxyz = nx * ny * nz;

  /*---------------------------------------------------------------------------
    |  p solver variables, discritization data handling.
    |
    |  allocation of space is simulation dimension-specific.
    |------------------------------------------------------------------------*/
  //  psa = (double *) calloc ((nxyz * nxyz), sizeof (double));
  //pija = (int *) calloc ((nxyz * nxyz), sizeof (int));

  //pa = allocate_double_3d (nx, ny, nz);
  //pb = allocate_double_3d (nx, ny, nz);
  //pc = allocate_double_3d (nx, ny, nz);

  //  for single-variable global matrices.
  //  A = allocate_double_2d (nxyz, nxyz);
  //  Am = allocate_double_2d (nxyz, nxyz);

  //  for tri-diagonal global matrices. Added by Jonathan Matthews 08/16/2012
  D = (double *) calloc ((nxyz), sizeof (double));    //main diagonal
  DU = (double *) calloc ((nxyz-1), sizeof (double)); //upper diagonal
  DL = (double *) calloc ((nxyz-1), sizeof (double)); //lower diagonal

  order = (int *) calloc ((nxyz+1), sizeof (int));
  vv = (double *) calloc ((nxyz+1), sizeof (double));

  x = (double *) calloc ((nxyz+1), sizeof (double));
  xm = (double *) calloc ((nxyz+1), sizeof (double));
  b = (double *) calloc ((nxyz+1), sizeof (double));
  bm = (double *) calloc ((nxyz+1), sizeof (double));

  //  large global jacobian matrix, nx * ny * nz * (nsolutes + 1)
  //  GA = allocate_double_2d ((nxyz * nsolutes + 1), (nxyz * nsolutes));
  //  Gb = (double *) calloc ((nxyz * nsolutes + 1), sizeof (double));

  /*---------------------------------------------------------------------------
    |  p solver variables, bc setup handling.
    |  assume horizontal configuration.
    |------------------------------------------------------------------------*/
  bc_p_inlet_pct = (double **) calloc ((ny * nz + 1), sizeof (double*));
  bc_p_outlet_pct = (double **) calloc ((nx * ny + 1), sizeof (double*));
  for_all_y
    for_all_z
  {
    bc_p_inlet_pct[iy] = (double *) calloc ((nx + 1), sizeof (double));
    bc_p_outlet_pct[iy] = (double *) calloc ((nx + 1), sizeof (double));
  }  //  end, for iy, iz.

}
/*===========================================================================*/
/*  END OF FILE                                                              */
/*===========================================================================*/
