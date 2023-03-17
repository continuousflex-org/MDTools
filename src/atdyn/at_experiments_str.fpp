!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_experiments_str_mod
!> @brief   structure of experiments information
!! @authors Osamu Miyashita (OM), Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_experiments_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  type, public :: s_emfit
    integer               :: natom        ! number of atoms used to allocate work arrays
    integer               :: nx           ! number of voxel in x direction
    integer               :: ny           ! number of voxel in y direction
    integer               :: nz           ! number of voxel in z direction
    integer               :: cutoff       ! cutoff parameter for force calculation ...
    integer               :: n_grid_cut_x ! number of grids included for kernel calculation
    integer               :: n_grid_cut_y
    integer               :: n_grid_cut_z
    integer               :: emfit_period

    real(wp)              :: norm_exp     ! norm of target map for force calculation
    real(wp)              :: x0           ! x coordinate of the origin voxel
    real(wp)              :: y0           ! y coordinate of the origin voxel
    real(wp)              :: z0           ! z coordinate of the origin voxel
    real(wp)              :: dx           ! size of one voxel
    real(wp)              :: dy           ! size of one voxel
    real(wp)              :: dz           ! size of one voxel
    real(wp)              :: sigma        ! "resolution" parameter for generating simulated EM map
    real(wp)              :: tolerance    ! tolerance for error
    real(wp)              :: weight       ! weight parameter for fitting biasing potential
    real(wp)              :: corrcoeff    ! store correlation coefficient
    real(wp)              :: force_norm   ! store force from fitting energy term

    integer,  allocatable :: ig(:,:)
    real(wp), allocatable :: target_map(:,:,:)      ! memory for target
    real(wp), allocatable :: simulated_map(:,:,:)   ! memory for simulated map
    real(wp), allocatable :: bound_x(:), bound_y(:), bound_z(:)
  end type s_emfit

  type, public :: s_emfit_img
    ! Emfit Image
    real(wp), allocatable :: target_img(:,:)      ! memory for target
    real(wp), allocatable :: simulated_img(:,:)   ! memory for simulated map

    real(wp), allocatable :: rot_coord(:,:)
    integer , allocatable :: pixels(:,:)
    real(wp), allocatable :: gaussians_saved(:,:,:)
    real(wp), allocatable :: emfit_img_force(:,:)

    integer               :: image_size
    real(wp)              :: pixel_size 
    real(wp)              :: roll_angle 
    real(wp)              :: tilt_angle 
    real(wp)              :: yaw_angle  
    real(wp)              :: shift_x    
    real(wp)              :: shift_y   

    integer               :: period
    real(wp)              :: sigma     
    integer               :: cutoff   

    real(wp)              :: rot_matrix(3,3)
    real(wp)              :: inv_rot_matrix(3,3)

  end type s_emfit_img

  ! structures
  type, public :: s_experiments
    logical                       :: do_emfit
    integer                       :: emfit_type
    type(s_emfit)                 :: emfit
    type(s_emfit_img)             :: emfit_img
  end type s_experiments

  ! parameters for allocatable variables
  integer,      public, parameter :: ExperimentsEmfit = 1
  integer,      public, parameter :: ExperimentsEmfitImg = 2

  character(*), public, parameter :: ExperimentsTypes(2)  = (/&
                                                          'VOLUME',&
                                                          'IMAGE '/)

  ! subroutines
  public :: alloc_experiments
  public :: dealloc_experiments
  public :: dealloc_experiments_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_experiments
  !> @brief        allocate experiments information
  !! @authors      TM
  !! @param[inout] experiments      : information of experiments
  !! @param[in]    variable  : allocatable variable
  !! @param[in]    var_size1 : size of variables
  !! @param[in]    var_size2 : size of variables (optional)
  !! @param[in]    var_size3 : size of variables (optional)
  !! @param[in]    var_size4 : size of variables (optional)
  !! @param[in]    var_size5 : size of variables (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_experiments(experiments, variable, var_size1, var_size2, &
                               var_size3, var_size4, var_size5)

    ! formal arguments
    type(s_experiments),     intent(inout) :: experiments
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1
    integer,      optional,  intent(in)    :: var_size2
    integer,      optional,  intent(in)    :: var_size3
    integer,      optional,  intent(in)    :: var_size4
    integer,      optional,  intent(in)    :: var_size5

    ! local variables
    integer                  :: alloc_stat, dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(ExperimentsEmfit)

      if (allocated(experiments%emfit%target_map)) then
        if (size(experiments%emfit%target_map(:,0,0)) == var_size1) return
        deallocate(experiments%emfit%target_map,     &
                   experiments%emfit%simulated_map,  &
                   experiments%emfit%bound_x,        &
                   experiments%emfit%bound_y,        &
                   experiments%emfit%bound_z,        &
                   experiments%emfit%ig,             &
                   stat = dealloc_stat)
      end if

      allocate(experiments%emfit%target_map   (0:var_size1-1, 0:var_size2-1, 0:var_size3-1), &
               experiments%emfit%simulated_map(0:var_size1-1, 0:var_size2-1, 0:var_size3-1), &
               experiments%emfit%bound_x(0:var_size1),                &
               experiments%emfit%bound_y(0:var_size2),                &
               experiments%emfit%bound_z(0:var_size3),                &
               experiments%emfit%ig(1:3,1:var_size4),                 &
               stat = alloc_stat)

      experiments%emfit%target_map   (0:var_size1-1, 0:var_size2-1, 0:var_size3-1) = 0.0_wp
      experiments%emfit%simulated_map(0:var_size1-1, 0:var_size2-1, 0:var_size3-1) = 0.0_wp
      experiments%emfit%bound_x(0:var_size1)                = 0.0_wp
      experiments%emfit%bound_y(0:var_size2)                = 0.0_wp
      experiments%emfit%bound_z(0:var_size3)                = 0.0_wp
      experiments%emfit%ig(1:3,1:var_size4)                 = 0

    case(ExperimentsEmfitImg)

      if (allocated(experiments%emfit_img%target_img)) then
        if (size(experiments%emfit_img%target_img(:,1)) == var_size1) return
        deallocate(experiments%emfit_img%target_img,      &
                   experiments%emfit_img%simulated_img,   &
                   experiments%emfit_img%rot_coord,       &
                   experiments%emfit_img%pixels,          &
                   experiments%emfit_img%gaussians_saved, &
                   experiments%emfit_img%emfit_img_force, &
                   stat = dealloc_stat)
      end if

      allocate(experiments%emfit_img%target_img   (var_size1, var_size1), &
               experiments%emfit_img%simulated_img(var_size1, var_size1), &
               experiments%emfit_img%rot_coord    (3, var_size2),         &
               experiments%emfit_img%pixels       (2, var_size2),         &
               experiments%emfit_img%gaussians_saved(var_size1, var_size1, var_size2), & 
               experiments%emfit_img%emfit_img_force(3, var_size2),         &     
               stat = alloc_stat)

      experiments%emfit_img%target_img   (1:var_size1, 1:var_size1) = 0.0_wp
      experiments%emfit_img%simulated_img(1:var_size1, 1:var_size1) = 0.0_wp
      experiments%emfit_img%rot_coord(1:3, 1:var_size2) = 0.0_wp
      experiments%emfit_img%pixels(1:2, 1:var_size2) = 0
      experiments%emfit_img%gaussians_saved(1:var_size1, 1:var_size1, 1:var_size2) = 0.0_wp
      experiments%emfit_img%emfit_img_force(1:3, 1:var_size2) = 0.0_wp

    case default

      call error_msg('Alloc_Experiments> bad variable')

    end select


    if (alloc_stat /= 0) call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_experiments

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_experiments
  !> @brief        deallocate experiments information
  !! @authors      TM
  !! @param[inout] experiments : experiments information
  !! @param[in]    variable    : allocatable variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_experiments(experiments, variable)

    ! formal arguments
    type(s_experiments),     intent(inout) :: experiments
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case (ExperimentsEmfit)

      if (allocated(experiments%emfit%target_map)) then
        deallocate(experiments%emfit%target_map,     &
                   experiments%emfit%simulated_map,  &
                   experiments%emfit%bound_x,        &
                   experiments%emfit%bound_y,        &
                   experiments%emfit%bound_z,        &
                   experiments%emfit%ig,             &
                   stat = dealloc_stat)
      end if

    case(ExperimentsEmfitImg)

      if (allocated(experiments%emfit_img%target_img)) then
        deallocate(experiments%emfit_img%target_img,      &
                   experiments%emfit_img%simulated_img,   &
                   experiments%emfit_img%rot_coord,       &
                   experiments%emfit_img%pixels,          &
                   experiments%emfit_img%gaussians_saved, &
                   experiments%emfit_img%emfit_img_force, &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Experiments> bad variable')

    end select 


    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_experiments

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_experiments_all
  !> @brief        deallocate all experiments information
  !! @authors      TM
  !! @param[inout] experiments : information of experiments
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_experiments_all(experiments)

    ! formal arguments
    type(s_experiments),     intent(inout) :: experiments


    call dealloc_experiments(experiments, ExperimentsEmfit)
    call dealloc_experiments(experiments, ExperimentsEmfitImg)

    return

  end subroutine dealloc_experiments_all

end module at_experiments_str_mod
