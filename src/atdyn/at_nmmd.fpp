!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_nmmd
!> @brief   Perform noraml mode molecular dynamics  (NMMD) simulation with velocity verlet algorithm
!! @authors Rémi Vuillemot
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_nmmd_mod

    use at_output_mod
    use at_dynvars_mod
    use at_ensemble_mod
    use at_constraints_mod
    use at_boundary_mod
    use at_pairlist_mod
    use at_energy_mod
    use at_energy_str_mod
    use at_output_str_mod
    use at_dynamics_str_mod
    use at_dynvars_str_mod
    use at_ensemble_str_mod
    use at_remd_str_mod
    use at_rpath_str_mod
    use at_constraints_str_mod
    use at_boundary_str_mod
    use at_pairlist_str_mod
    use at_enefunc_str_mod
    use at_enefunc_gbsa_mod
    use at_enefunc_gamd_mod
    use molecules_str_mod
    use fileio_rst_mod
    use fileio_pdb_mod
    use fileio_mod
    use timers_mod
    use random_mod
    use math_libs_mod
    use messages_mod
    use constants_mod
    use mpi_parallel_mod 
    use string_mod

    implicit none
    private

    ! subroutines
    public  :: nmmd_dynamics
    private :: initial_nmmd
    private :: langevin_thermostat_vv1
    private :: langevin_thermostat_vv2
    !private :: compute_nma
    private :: read_nma
    private :: compute_energy_nmmd
    !private :: output_nmmd
    !private :: check_nma
    !private :: include_id_to_nm_prefix

contains

    !======1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    nmmd_dynamics
    !> @brief        NMMD integrator
    !! @authors      RV
    !! @param[inout] output      : output information
    !! @param[inout] molecule    : molecule information
    !! @param[inout] enefunc     : potential energy functions information
    !! @param[inout] dynvars     : dynamic variables information
    !! @param[inout] dynamics    : dynamics information
    !! @param[inout] pairlist    : non-bond pair list information
    !! @param[inout] boundary    : boundary conditions information
    !! @param[inout] constraints : bond constraints information
    !! @param[inout] ensemble    : ensemble information
    !
    !======1=========2=========3=========4=========5=========6=========7=========8

    subroutine nmmd_dynamics(output, molecule, enefunc, dynvars, dynamics, &
                                pairlist, boundary, constraints, ensemble)

        ! formal arguments
        type(s_output), intent(inout) :: output
        type(s_molecule), target, intent(inout) :: molecule
        type(s_enefunc), intent(inout) :: enefunc
        type(s_dynvars), target, intent(inout) :: dynvars
        type(s_dynamics), target, intent(inout) :: dynamics
        type(s_pairlist), intent(inout) :: pairlist
        type(s_boundary), intent(inout) :: boundary
        type(s_constraints), intent(inout) :: constraints
        type(s_ensemble), intent(inout) :: ensemble

        ! local variables
        real(wp)                  :: simtim, dt
        integer                   :: i, j, natom
        integer                   :: nsteps, istart, iend
        character                 :: num*9
        logical                   :: savefile
        real(wp), pointer :: coord(:, :), coord_pbc(:,:), vel(:,:)
        real(wp)                  :: shrink_ratio_x, shrink_ratio_y, shrink_ratio_z 
        integer                   :: count_qm
        integer                   :: icount
        integer :: nmodes, exitstatus
        logical :: file_exists
        real(wp),         pointer :: force_omp(:,:,:)

        coord => dynvars%coord
        vel       => dynvars%velocity
        coord_pbc => dynvars%coord_pbc
        force_omp => dynvars%force_omp
        natom = molecule%num_atoms
        istart = dynamics%istart_step
        iend = dynamics%iend_step
        nsteps = dynamics%nsteps
        simtim = dynamics%initial_time
        nmodes = dynamics%nm_number
        dt = dynamics%timestep/AKMA_PS


        ! first-step MD
        !
        if (.not. dynamics%restart) then
            call initial_nmmd(output, molecule, enefunc, dynamics, pairlist, &
                                    boundary, ensemble, constraints, dynvars)
        else
            ! Remove velocity of fixed atoms
            ! 
            if (constraints%num_fixatm > 0) &
            call clear_fixatm_component(constraints, natom, dynvars%velocity)
        end if


        !
        ! stop NPT NVT and not langevin
        !
        if (ensemble%ensemble == EnsembleNPAT .or. &
            ensemble%ensemble == EnsembleNPT .or. &
            ensemble%ensemble == EnsembleNPgT.or. &
            ensemble%ensemble == EnsembleNVE) &
            call error_msg('NMMMD_dynamics> NVT only is allowed in NMMD')
        if (ensemble%tpcontrol /= TpcontrolLangevin) &
            call error_msg('NMMMD_dynamics> Langevin thermostat only is allowed in NMMD')

        ! Main MD loop

        do i = istart, iend

            call timer(TimerIntegrator, TimerOn)

            simtim = simtim + dt*AKMA_PS
            dynvars%time = simtim
            dynvars%step = i

            ! Compute velocities(t + 1/2 dt) and coordinates(t + dt)
            !
            call langevin_thermostat_vv1(molecule, dynamics, ensemble, &
                                            boundary, constraints, dynvars)


            call timer(TimerIntegrator, TimerOff)

            ! calculate potential energy(t + dt), force(t + dt), and virial(t + dt)
            !
            call compute_energy_nmmd(molecule, dynamics, enefunc, &
                    boundary, dynvars, pairlist, constraints, i)

            call timer(TimerIntegrator, TimerOn)

            call langevin_thermostat_vv2(molecule, dynamics, ensemble, &
                                                boundary, constraints, dynvars)

            call timer(TimerIntegrator, TimerOff)

            ! Remove translational and rotational motion about COM(t + dt)
            !
            if (dynamics%stoptr_period > 0) then

                if (mod(i, dynamics%stoptr_period) == 0) then
                    call stop_trans_rotation(molecule%num_atoms, molecule%mass, &
                                             dynamics%stop_com_translation, &
                                             dynamics%stop_com_rotation, &
                                             dynvars%coord, &
                                             constraints%fixatm, &
                                             dynvars%velocity)
                end if

            end if

            ! Update nonbond pairlist for coordinates(t + 2dt)
            !
            if (dynamics%nbupdate_period > 0) then
                if (mod(i, dynamics%nbupdate_period) == 0 .and. real_calc) then

                    ! Update boundary if pressure is controlled
                    !
                    if (ensemble%use_barostat) then
                        if (enefunc%forcefield == ForcefieldRESIDCG) then
                        call update_boundary_cg(enefunc%cg_pairlistdist_ele, &
                            enefunc%cg_pairlistdist_126,                     &
                            enefunc%cg_pairlistdist_PWMcos,                  &
                            enefunc%cg_pairlistdist_DNAbp,                   &
                            enefunc%cg_pairlistdist_exv,                     &
                            boundary)
                        else
                        call update_boundary(enefunc%table%table,               &
                            enefunc%pairlistdist,                               &
                            boundary)
                        end if
                    end if

                    ! Shrink box
                    !
                    if (dynamics%shrink_box) then
                        if (.not. ensemble%use_barostat .and. &
                            mod(i, dynamics%shrink_period) == 0) then
                        boundary%box_size_x     = boundary%box_size_x - dynamics%dbox_x
                        boundary%box_size_y     = boundary%box_size_y - dynamics%dbox_y
                        boundary%box_size_z     = boundary%box_size_z - dynamics%dbox_z
                        boundary%box_size_x_ref = boundary%box_size_x
                        boundary%box_size_y_ref = boundary%box_size_y
                        boundary%box_size_z_ref = boundary%box_size_z

                        shrink_ratio_x = boundary%box_size_x &
                                        / (boundary%box_size_x + dynamics%dbox_x)
                        shrink_ratio_y = boundary%box_size_y &
                                        / (boundary%box_size_y + dynamics%dbox_y)
                        shrink_ratio_z = boundary%box_size_z &
                                        / (boundary%box_size_z + dynamics%dbox_z)

                        do j = 1, natom
                            coord(1,j) = coord(1,j) * shrink_ratio_x
                            coord(2,j) = coord(2,j) * shrink_ratio_y
                            coord(3,j) = coord(3,j) * shrink_ratio_z
                            vel(1,j)   = vel(1,j) * shrink_ratio_x
                            vel(2,j)   = vel(2,j) * shrink_ratio_y
                            vel(3,j)   = vel(3,j) * shrink_ratio_z
                        end do

                        if (enefunc%forcefield == ForcefieldRESIDCG) then
                            call update_boundary_cg(enefunc%cg_pairlistdist_ele, &
                                enefunc%cg_pairlistdist_126,                     &
                                enefunc%cg_pairlistdist_PWMcos,                  &
                                enefunc%cg_pairlistdist_DNAbp,                   &
                                enefunc%cg_pairlistdist_exv,                     &
                                boundary)
                        else
                            call update_boundary(enefunc%table%table,               &
                                enefunc%pairlistdist,                               &
                                boundary)
                        end if
                        end if
                    end if

                    call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                                        coord_pbc, pairlist)


                end if
            end if

            ! Output energy(t + dt) and dynamical variables(t + dt)
            !
            call output_md(output, molecule, enefunc, dynamics, boundary, &
                           ensemble, dynvars)

        end do

        !call output_nmmd(dynvars, dynamics, natom, nmodes)

        return

    end subroutine nmmd_dynamics

    !======1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    initial_nmmd
    !> @brief        compute the first step (0+dt)
    !! @authors      RV
    !! @param[in]    output      : output information
    !! @param[in]    molecule    : molecule information
    !! @param[in]    enefunc     : potential energy functions information
    !! @param[in]    dynamics    : dynamics information
    !! @param[in]    pairlist    : non-bond pair list information
    !! @param[in]    boundary    : boundary conditions information
    !! @param[in]    ensemble    : ensemble information
    !! @param[inout] constraints : bond constraints information
    !! @param[inout] dynvars     : dynamic variables information
    !! @param[in]    nmodes      : number of modes
    !
    !======1=========2=========3=========4=========5=========6=========7=========8

    subroutine initial_nmmd(output, molecule, enefunc, dynamics, pairlist, &
                               boundary, ensemble, constraints, dynvars)

        ! formal arguments
        type(s_output), intent(in)    :: output
        type(s_molecule), target, intent(inout) :: molecule
        type(s_enefunc), intent(inout) :: enefunc
        type(s_dynamics), target, intent(inout) :: dynamics
        type(s_pairlist), intent(inout)    :: pairlist
        type(s_boundary), intent(in)    :: boundary
        type(s_ensemble), intent(in)    :: ensemble
        type(s_constraints), intent(inout) :: constraints
        type(s_dynvars), target, intent(inout) :: dynvars

        ! local variables
        real(wp)                  :: dt, simtim
        integer                   :: i,j, natom, nmodes, nm_init_size

        real(wp), pointer :: coord(:, :), coord_ref(:, :)
        real(wp), pointer :: vel(:, :)
        real(wp), pointer :: force(:, :), mass(:), inv_mass(:)
        real(wp), pointer :: nm_vel(:), nm_vectors(:,:,:)
        integer, pointer :: iseed
        character                 :: num*9
        logical                   :: savefile,file_exists
        character(MaxFilename) :: nm_file

        ! use pointers
        natom = molecule%num_atoms
        nmodes = dynamics%nm_number
        coord => dynvars%coord
        coord_ref => dynvars%coord_ref
        force => dynvars%force
        vel => dynvars%velocity
        mass => molecule%mass
        iseed => dynamics%iseed
        simtim = dynamics%initial_time
        nm_vel                  =>dynvars%nm_vel
        nm_vectors              =>dynvars%nm_vectors
        dynvars%time = simtim
        dynvars%step = 0
        nm_file = dynamics%nm_file

        ! Read NMA vectors
        call read_nma(nm_vectors, nm_file, natom, nmodes)


        ! Init parameters to zero
        coord_ref(1:3,1:natom) = coord(1:3,1:natom)

        ! Initialize NMMD params
        dynvars%nm_amp(1:nmodes) = 0.0_wp
        dynvars%md_coord(1:3, 1:natom) = 0.0_wp

        nm_init_size = split_num(dynamics%nm_init)
        if (nm_init_size == nmodes) then
            call split(nmodes, nmodes, dynamics%nm_init, dynvars%nm_amp)
            write(MsgOut,'(A)') 'Initialize_NMA_Amp> Initialize NMA amplitudes'
            write(MsgOut, '(A)') '  nm_init    = ' //dynamics%nm_init
            write(MsgOut,'(A)') ''
        end if


        ! generate initial velocities
        !
        call generate_velocity(ensemble%temperature, molecule%num_atoms, &
                               molecule%mass, iseed, dynvars%velocity)
        

        ! init NMMD velocities
        do i = 1, natom
            do j = 1, nmodes
                nm_vel(j) = nm_vel(j) + nm_vectors(1, i, j)*vel(1, i)
                nm_vel(j) = nm_vel(j) + nm_vectors(2, i, j)*vel(2, i)
                nm_vel(j) = nm_vel(j) + nm_vectors(3, i, j)*vel(3, i)
            end do
        end do

        if (constraints%num_fixatm > 0) &
            call clear_fixatm_component(constraints, natom, dynvars%velocity)

        call stop_trans_rotation(molecule%num_atoms, molecule%mass, &
                                 dynamics%stop_com_translation, &
                                 dynamics%stop_com_rotation, &
                                 dynvars%coord, &
                                 constraints%fixatm, &
                                 dynvars%velocity)


        ! calculate energy(0) and forces(0)
        !
        call compute_energy_nmmd(molecule, dynamics, enefunc, &
                boundary, dynvars, pairlist, constraints, 0)
        
        
        ! output dynvars(0)
        !
        call compute_dynvars(molecule, enefunc, dynamics, boundary, ensemble, &
                             dynvars)
        call output_dynvars(output, enefunc, dynvars, ensemble)

        dynamics%restart = .true.

        ! at this point
        !   coord, velocity, and force are at 0

        return

    end subroutine initial_nmmd

    !======1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    compute_energy_nmmd
    !> @brief        compute energy of NMMD parameters
    !! @authors      RV
    !
    !======1=========2=========3=========4=========5=========6=========7=========8


    subroutine compute_energy_nmmd(molecule, dynamics, enefunc, &
        boundary, dynvars, pairlist, constraints, iteration)

        ! formal arguments
        type(s_dynamics), intent(in)    :: dynamics
        type(s_molecule), target, intent(inout)    :: molecule
        type(s_enefunc), target, intent(inout) :: enefunc
        type(s_pairlist), target, intent(inout) :: pairlist
        type(s_boundary), target, intent(in)    :: boundary
        type(s_dynvars), target, intent(inout) :: dynvars
        type(s_constraints), target, intent(inout) :: constraints
        integer, intent(in)     ::  iteration

        integer :: i,j ,nmodes, natom
        nmodes = dynamics%nm_number
        natom = molecule%num_atoms

        call compute_energy(molecule, enefunc, pairlist, boundary, &
                    mod(iteration, dynamics%eneout_period) == 0, &
                    enefunc%nonb_limiter, &
                    dynvars%coord, &
                    dynvars%trans, &
                    dynvars%coord_pbc, &
                    dynvars%energy, &
                    dynvars%temporary, &
                    dynvars%force, &
                    dynvars%force_omp, &
                    dynvars%virial, &
                    dynvars%virial_extern, &
                    constraints)


        ! update NMMD forces
        dynvars%nm_force(1:nmodes) = 0.0_wp
        dynvars%nm_random_force(1:nmodes) = 0.0_wp
        do i = 1, natom
            do j = 1, nmodes
                dynvars%nm_force(j) = dynvars%nm_force(j) + dynvars%nm_vectors(1, i, j)*dynvars%force(1, i)
                dynvars%nm_force(j) = dynvars%nm_force(j) + dynvars%nm_vectors(2, i, j)*dynvars%force(2, i)
                dynvars%nm_force(j) = dynvars%nm_force(j) + dynvars%nm_vectors(3, i, j)*dynvars%force(3, i)
                dynvars%nm_random_force(j) = dynvars%nm_random_force(j) + dynvars%nm_vectors(1, i, j)*dynvars%random_force(1, i)
                dynvars%nm_random_force(j) = dynvars%nm_random_force(j) + dynvars%nm_vectors(2, i, j)*dynvars%random_force(2, i)
                dynvars%nm_random_force(j) = dynvars%nm_random_force(j) + dynvars%nm_vectors(3, i, j)*dynvars%random_force(3, i)
            end do
        end do

        return

    end subroutine compute_energy_nmmd


    !======1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    langevin_thermostat_vv1
    !> @brief        control temperature using Langevin thermostat and barostat
    !! @authors      RV
    !! @param[in]    molecule    : molecule information
    !! @param[in]    dynamics    : dynamics information
    !! @param[in]    ensemble    : ensemble information
    !! @param[in]    boundary    : ensemble information
    !! @param[inout] constraints : bond constraints information
    !! @param[inout] dynvars     : dynamic variables information
    !
    !======1=========2=========3=========4=========5=========6=========7=========8


    subroutine langevin_thermostat_vv1(molecule, dynamics, ensemble, &
                                       boundary, constraints, dynvars)

        ! formal arguments
        type(s_dynamics), intent(in)    :: dynamics
        type(s_molecule), target, intent(in)    :: molecule
        type(s_ensemble), target, intent(inout) :: ensemble
        type(s_boundary),          intent(in)    :: boundary
        type(s_constraints), target, intent(inout) :: constraints
        type(s_dynvars), target, intent(inout) :: dynvars

        ! local variables
        real(wp)                  :: dt, h_dt, nm_inv_mass, nm_dt
        real(wp)                  :: temperature, gamma_t, scale_v
        real(wp)                  :: factor, sigma, temp0
        real(wp)                  :: v1, v2, rsq, grandom(1:3)
        real(wp)                  :: kBT, TWO_PI
        real(wp)                  :: tvel(3)
        integer                   :: i, j, natom, nmodes

        real(wp), pointer :: vel(:, :)
        real(wp), pointer :: coord(:, :), coord_ref(:, :), force(:, :)
        real(wp), pointer :: mass(:), inv_mass(:), temporary(:, :)
        real(wp), pointer :: viri(:, :), viri_const(:, :)
        real(wp), pointer :: random_force(:, :)
        logical, pointer :: fixatm(:)
        real(wp), pointer :: nm_amp(:), nm_vel(:), nm_force(:), nm_random_force(:)
        real(wp), pointer :: nm_coord(:,:), md_coord(:,:), nm_vectors(:,:,:)

        ! use pointers
        !
        dt = dynamics%timestep/AKMA_PS
        nm_dt = dynamics%nm_dt/AKMA_PS
        temp0 = ensemble%temperature
        gamma_t = ensemble%gamma_t*AKMA_PS
        natom = molecule%num_atoms
        nmodes = dynamics%nm_number
        nm_inv_mass= 1.0_wp/dynamics%nm_mass
        mass => molecule%mass
        inv_mass => molecule%inv_mass
        coord => dynvars%coord
        coord_ref => dynvars%coord_ref
        vel => dynvars%velocity
        force => dynvars%force
        random_force => dynvars%random_force
        temporary => dynvars%temporary
        viri => dynvars%virial
        viri_const => dynvars%virial_const
        fixatm => constraints%fixatm
        nm_amp                  =>dynvars%nm_amp
        nm_vel                  =>dynvars%nm_vel
        nm_force                =>dynvars%nm_force
        nm_random_force         =>dynvars%nm_random_force
        nm_coord                =>dynvars%nm_coord
        md_coord                =>dynvars%md_coord
        nm_vectors              =>dynvars%nm_vectors

        ! setup variables
        !
        kBT = KBOLTZ*temp0
        TWO_PI = 2.0_wp*PI
        h_dt = dt/2.0_wp

        ! random force
        !

        if (dynvars%step == 1) then

            factor = 2.0_wp*gamma_t*KBOLTZ*temp0/h_dt
            do j = 1, natom
                if (fixatm(j)) cycle

                sigma = sqrt(mass(j)*factor)

                rsq = 2.0_wp
                do while (rsq >= 1.0_wp)
                    v1 = 2.0_wp*random_get() - 1.0_wp
                    v2 = 2.0_wp*random_get() - 1.0_wp
                    rsq = v1*v1 + v2*v2
                end do
                rsq = sqrt(-2.0_wp*log(rsq)/rsq)
                grandom(1) = rsq*v1

                rsq = 2.0_wp
                do while (rsq >= 1.0_wp)
                    v1 = 2.0_wp*random_get() - 1.0_wp
                    v2 = 2.0_wp*random_get() - 1.0_wp
                    rsq = v1*v1 + v2*v2
                end do
                rsq = sqrt(-2.0_wp*log(rsq)/rsq)
                grandom(2) = rsq*v1

                rsq = 2.0_wp
                do while (rsq >= 1.0_wp)
                    v1 = 2.0_wp*random_get() - 1.0_wp
                    v2 = 2.0_wp*random_get() - 1.0_wp
                    rsq = v1*v1 + v2*v2
                end do
                rsq = sqrt(-2.0_wp*log(rsq)/rsq)
                grandom(3) = rsq*v1

                random_force(1:3, j) = sigma*grandom(1:3)

            end do

        end if

        scale_v = exp(-gamma_t*0.5_wp*nm_dt)

        ! thermostat nm
        nm_vel(1:nmodes) = nm_vel(1:nmodes)*scale_v 

        ! VV1 nm
        nm_vel(1:nmodes) = nm_vel(1:nmodes)+ 0.5_wp*nm_dt* (nm_force(1:nmodes)+nm_random_force(1:nmodes))* nm_inv_mass
        nm_amp(1:nmodes) = nm_amp(1:nmodes) + nm_dt*nm_vel(1:nmodes)



        scale_v = exp(-gamma_t*0.5_wp*dt)
        do i = 1, natom

            ! thermostat md
            vel(1:3,i) = vel(1:3,i)*scale_v 

            ! VV2 md            
            vel(1:3,i) = vel(1:3,i) + 0.5_wp*dt*(force(1:3,i)+random_force(1:3,i))*inv_mass(i) 
            md_coord(1:3, i) = md_coord(1:3, i) + vel(1:3, i)*dt
            
            ! nm_amp -> nm_coords
            nm_coord(1:3, i) = 0.0_wp
            do j = 1, nmodes
                nm_coord(1:3, i) = nm_coord(1:3, i) + (nm_amp(j)*nm_vectors(1:3, i, j))
            end do

            ! coord update
            coord(1:3, i) = coord_ref(1:3, i) + md_coord(1:3, i) + nm_coord(1:3, i)
        end do

        return

    end subroutine langevin_thermostat_vv1
    

    !======1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    langevin_thermostat_vv2
    !> @brief        control temperature using Langevin thermostat and barostat
  !! @authors      RV
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[in]    boundary    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] dynvars     : dynamic variables information
    !
    !======1=========2=========3=========4=========5=========6=========7=========8

    subroutine langevin_thermostat_vv2(molecule, dynamics, ensemble, &
                                       boundary, constraints, dynvars)

        ! formal arguments
        type(s_molecule), target, intent(in)    :: molecule
        type(s_dynamics), intent(in)    :: dynamics
        type(s_ensemble), target, intent(inout) :: ensemble
        type(s_boundary),  intent(in)    :: boundary
        type(s_constraints), target, intent(inout) :: constraints
        type(s_dynvars), target, intent(inout) :: dynvars

        ! local variables
        real(wp)                  :: dt, inv_dt, temp0, nm_inv_mass, nm_dt
        real(wp)                  :: factor, sigma
        real(wp)                  :: gamma_t
        real(wp)                  :: v1, v2, rsq, grandom(1:3), fcm(0:4)
        real(wp)                  :: h_dt, vel_scale
        integer                   :: j, natom, nmodes

        real(wp), pointer :: mass(:), inv_mass(:)
        real(wp), pointer :: temporary(:, :)
        real(wp), pointer :: vel(:, :),force(:, :)
        real(wp), pointer :: virial(:, :), viri_const(:, :)
        real(wp), pointer :: random_f(:, :), nm_vel(:), nm_force(:), nm_random_force(:)
        logical, pointer :: fixatm(:)

        ! use pointers
        !
        nmodes = dynamics%nm_number
        nm_inv_mass= 1.0_wp/dynamics%nm_mass
        dt = dynamics%timestep/AKMA_PS
        nm_dt = dynamics%nm_dt/AKMA_PS
        inv_dt = 1.0_wp/dt
        temp0 = ensemble%temperature
        gamma_t = ensemble%gamma_t*AKMA_PS
        mass => molecule%mass
        inv_mass => molecule%inv_mass
        natom = molecule%num_atoms
        vel => dynvars%velocity
        force => dynvars%force
        random_f => dynvars%random_force
        temporary => dynvars%temporary
        virial => dynvars%virial
        viri_const => dynvars%virial_const
        fixatm => constraints%fixatm

        nm_vel => dynvars%nm_vel
        nm_force => dynvars%nm_force
        nm_random_force => dynvars%nm_random_force

        ! time step
        !
        h_dt = dt/2.0_wp

        ! random force
        !
        fcm(0:3) = 0.0_wp
        factor = 2.0_wp*gamma_t*KBOLTZ*temp0/dt

        do j = 1, natom
            if (fixatm(j)) cycle

            sigma = sqrt(mass(j)*factor)

            rsq = 2.0_wp
            do while (rsq >= 1.0_wp)
                v1 = 2.0_wp*random_get() - 1.0_wp
                v2 = 2.0_wp*random_get() - 1.0_wp
                rsq = v1*v1 + v2*v2
            end do
            rsq = sqrt(-2.0_wp*log(rsq)/rsq)
            grandom(1) = rsq*v1

            rsq = 2.0_wp
            do while (rsq >= 1.0_wp)
                v1 = 2.0_wp*random_get() - 1.0_wp
                v2 = 2.0_wp*random_get() - 1.0_wp
                rsq = v1*v1 + v2*v2
            end do
            rsq = sqrt(-2.0_wp*log(rsq)/rsq)
            grandom(2) = rsq*v1

            rsq = 2.0_wp
            do while (rsq >= 1.0_wp)
                v1 = 2.0_wp*random_get() - 1.0_wp
                v2 = 2.0_wp*random_get() - 1.0_wp
                rsq = v1*v1 + v2*v2
            end do
            rsq = sqrt(-2.0_wp*log(rsq)/rsq)
            grandom(3) = rsq*v1

            random_f(1:3, j) = sigma*grandom(1:3)

        end do

        ! VV2
        !
        do j = 1, natom
            vel(1:3,j) = vel(1:3,j) + h_dt*(force(1:3,j)+random_f(1:3,j))*inv_mass(j)
        end do
        nm_vel(1:nmodes) = nm_vel(1:nmodes) + 0.5_wp*nm_dt*(nm_force(1:nmodes)+nm_random_force(1:nmodes))*nm_inv_mass

        ! termostat
        !
        vel_scale = exp(-gamma_t*h_dt)
        do j = 1, natom
            vel(1:3,j) = vel(1:3,j)*vel_scale
        end do
        vel_scale = exp(-gamma_t* 0.5_wp*nm_dt)
        nm_vel(1:nmodes) = nm_vel(1:nmodes)*vel_scale 

        return

    end subroutine langevin_thermostat_vv2

    !======1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    output_nmmd
    !> @brief        output NMMD
    !! @authors      RV
    !
    !======1=========2=========3=========4=========5=========6=========7=========8

    ! subroutine output_nmmd(dynvars, dynamics, natom, nmodes)
    !     type(s_dynvars), intent(inout) :: dynvars
    !     type(s_dynamics), intent(inout) :: dynamics
    !     integer, intent(inout) :: natom
    !     integer, intent(inout) :: nmodes

    !     integer :: unit_no
    !     character(MaxFilename)::restart_file, nm_prefix
    !     logical :: file_exists

    !     nm_prefix = dynamics%nm_prefix
    !     call include_id_to_nm_prefix(nm_prefix)

    !     if (dynamics%crdout_period > 0) then
    !         if (mod(dynvars%step,dynamics%crdout_period) == 0) then
                
    !             ! Write nm amplitudes
    !             inquire(file=trim(nm_prefix)//".nmmd", exist=file_exists)
    !             if (file_exists) then
    !                 call open_file (unit_no=unit_no, filename=trim(nm_prefix)//".nmmd", in_out=IOFileOutputReplace)
    !             else
    !                 call open_file (unit_no=unit_no, filename=trim(nm_prefix)//".nmmd", in_out=IOFileOutputNew)
    !             endif
    !             write(unit_no, *) dynvars%nm_amp
    !             call close_file(unit_no)

    !         endif
    !     endif

    !     return
        
    ! end subroutine output_nmmd

    !======1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    compute_nma
    !> @brief        compute normal mode analysis using Elnemo
    !! @authors      RV
    !
    !======1=========2=========3=========4=========5=========6=========7=========8

    ! subroutine compute_nma(molecule, dynvars, dynamics)
    !     type(s_molecule), intent(in) :: molecule
    !     type(s_dynvars), intent(in) :: dynvars
    !     type(s_dynamics),  intent(in) :: dynamics

    !     integer :: exitstatus, unit_no, nmodes
    !     logical :: file_exists
    !     character(MaxFilename)::nm_prefix

    !     nmodes = dynamics%nm_number
    !     nm_prefix = dynamics%nm_prefix
    !     call include_id_to_nm_prefix(nm_prefix)

    !     write(MsgOut,'(A)') 'Compute_NMA> Computing Normal Mode Analysis (ElNemo)'
    !     ! WRITE PDB
    !     open (unit=66, file=trim(nm_prefix)//"/run.pdb")
    !     call write_restart_pdb(66, molecule, dynvars%coord)
    !     close (66)

    !     ! ELNEMO PDBMAT
    !     inquire(file=trim(nm_prefix)//"/pdbmat.dat", exist=file_exists)
    !     if (file_exists) then
    !         call open_file (unit_no=unit_no, filename=trim(nm_prefix)//"/pdbmat.dat", in_out=IOFileOutputReplace)
    !     else
    !         call open_file (unit_no=unit_no, filename=trim(nm_prefix)//"/pdbmat.dat", in_out=IOFileOutputNew)
    !     endif
    !     write (unit_no, *) "Coordinate FILENAME        = run.pdb"
    !     write (unit_no, *) "MATRIx FILENAME            = pdbmat.sdijf"
    !     write (unit_no, '(A,F10.3)') "INTERACtion DISTance CUTOF = ", dynamics%elnemo_cutoff
    !     write (unit_no, *) "INTERACtion FORCE CONStant = 10.000000"
    !     write (unit_no, *) "Origin of MASS values      =    CON "
    !     write (unit_no, *) "Output PRINTing level      =      0 "
    !     call close_file (unit_no)

    !     ! ELNEMO DIAGRTB
    !     inquire(file=trim(nm_prefix)//"/diagrtb.dat", exist=file_exists)
    !     if (file_exists) then
    !         call open_file (unit_no=unit_no, filename=trim(nm_prefix)//"/diagrtb.dat", in_out=IOFileOutputReplace)
    !     else
    !         call open_file (unit_no=unit_no, filename=trim(nm_prefix)//"/diagrtb.dat", in_out=IOFileOutputNew)
    !     endif
    !     write (unit_no, *) "MATRix filename            = pdbmat.sdijf"
    !     write (unit_no, *) "COORdinates filename       = run.pdb"
    !     write (unit_no, *) "Eigenvector OUTPut filename= diagrtb.eigenfacs"
    !     write (unit_no, '(A, I10)') "Nb of VECTors required     = ", nmodes + 6
    !     write (unit_no, *) "EigeNVALues chosen         =       LOWE "
    !     write (unit_no, *) "Type of SUBStructuring     =       NONE "
    !     write (unit_no, '(A,I10)') "Nb of residues per BLOck   =        ", dynamics%elnemo_rtb_block
    !     write (unit_no, *) "Origin of MASS values      =       CONS "
    !     write (unit_no, *) "Temporary files cleaning   =       ALL  "
    !     write (unit_no, *) "MATRix FORMat              =       FREE "
    !     write (unit_no, *) "Output PRINting level      =          0 "
    !     call close_file (unit_no)

    !     ! ELNEMO RUN COMMAND
    !     call execute_command_line(                                          &
    !             "cd "// trim(nm_prefix) // " ; "                            &
    !             // trim(dynamics%elnemo_path)//"/nma_elnemo_pdbmat > pdbmat.log ; "  &
    !             // trim(dynamics%elnemo_path)//"/nma_diagrtb > diagrtb.log ;         &
    !             rm -f pdbmat.dat_run diagrtb.dat_run pdbmat.sdijf"          &
    !             , wait=.true., exitstat=exitstatus)
    !     if (exitstatus /= 0) then
    !         call error_msg('Mode PB')
    !     end if
    !     write(MsgOut,'(A)') ''

    !     return

    ! end subroutine compute_nma

    !======1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    read_nma
    !> @brief        read normal modes
    !! @authors      RV
    !
    !======1=========2=========3=========4=========5=========6=========7=========8

    subroutine read_nma(nm_vectors, nm_file, natom, nmodes)
        real(wp), dimension(:, :, :), intent(inout) :: nm_vectors
        character(MaxFilename), intent(in) :: nm_file
        integer, intent(in) :: nmodes
        integer, intent(in) :: natom

        integer :: i,j, unit_no


        call open_file(unit_no=unit_no, filename=nm_file, in_out=IOFileInput)

        do i = 1, nmodes
            read (unit_no, '(A)')
            read (unit_no, '(A)')
            do j = 1, natom
                    read (unit_no, *) nm_vectors(:,j,i)
            end do
        end do
        call close_file (unit_no)

        write(MsgOut,'(A)') 'Read_NMA> Read NMA'
        write(MsgOut, '(A)') '  nm_file    = '
        write(MsgOut, '(A)') '    '//nm_file
        write(MsgOut,'(A)') ''
        
        return

    end subroutine read_nma

end module at_nmmd_mod
