  module module_pbl3d

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                         !
    ! Purpose: 1) Calculates turbulent fluxes of momentum, heat and moisture  !
    !                                                                         !
    !          2) Calculates tendencies due to 3D turbulent mixing            !
    !                                                                         !
    ! Author: Pedro A. Jimenez & Timothy W. Juliano                           !
    !                                                                         !
    ! Methods: Uses the Mellor and Yamada turbulence closure model            !
    !                                                                         !
    ! References: M73: Mellor, G. L., 1973: Analytic prediction of the        !
    !                  properties of stratified planetary surface layers. J.  !
    !                  Atmos. Sci., 30, 1061-1069.                            !
    !                                                                         !
    !             MY74: Mellor, G. L., and T. Yamada, 1974: A hierarchy of    !
    !                   turbulence closure models for planetary boundary      !
    !                   layers. J. Atmos. Sci., 31, 1791-1806.                !
    !                                                                         !
    !             YM75: Yamada, T., and G. Mellor, 1975: A simulation of the  !
    !                   Wangara atmospheric boundary layer data. J. Atmos.    !
    !                   Sci., 32, 2309-2329.                                  !
    !                                                                         !
    !              Y75: Yamada, T., The critical richardson number and the    !
    !                   ratio of the eddy transport coefficients obtained     !
    !                   from a turbulence closure model. JAS, 32, 926 - 933   !
    !                                                                         !
    !             MY82: Mellor, G. L., and T. Yamada, 1982: Development of a  !
    !                   turbulent closure model for geophysical fluid         !
    !                   problems. Rev. Geophys.                               !
    !                   Space Phys., 20, 851-875.                             !
    !                                                                         !
    !              Y83: Yamada, T., 1983: Simulations of nocturnal drainage   !
    !                   flows by a q2/l turbulence closure model              !
    !                                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use, intrinsic :: iso_fortran_env

    use module_model_constants, only: KARMAN, G, P608, T0, CP, RCP
    use module_state_description, only: param_first_scalar
    use module_big_step_utilities_em, only: grid_config_rec_type

    implicit none

    private

    public :: Init_pbl3d, &
              Init_temp_arrays_for_substep, Update_wrf_tends_temp_state_and_zero_tends, Allocate_temp_tends, Deallocate_temp_tends, &
              Calc_turb_fluxes, Vertical_turb_mix, Horizontal_turb_mix

    real, parameter :: G_OVER_TREF = G / T0
    real, parameter :: DWIND_DZ_MIN = 1.0E-10
    real            :: Q_SQ_MIN, TURB_FLUX_MIN
    logical, parameter :: SOLVE_WITH_SOR = .true.
    logical, parameter :: USE_HIGH_ORDER = .true.

    integer, parameter :: SP = REAL32
    integer, parameter :: DP = REAL64
    integer, parameter :: QP = REAL128
    integer, parameter :: SOLVER_PREC = DP
    integer, parameter :: EXTRA_PREC = QP
    integer, parameter :: WORKING_PREC = SP

      ! Model constants
    real :: a1, a2, b1, b2, c1, c2, c3, alpha, lambda

      ! For the level 2 model
    real :: rif_c, rf1, rf2, cm, ch2
    real :: fac1, fac2, fac3, fac4

    integer, parameter :: I_TO_TEST = 4, J_TO_TEST = 4, K_TO_TEST = 3

  contains

    subroutine Init_temp_arrays_for_substep (u_2, v_2, w_2, t_2, moist, rho, q_sq_2, &
              u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp, &
              ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Prepare for sub-stepping                                     !
      !              1. Save copy of state variables                          !
      !              2. Allocate temporary tendency arrays                    !
      !              3. Zero out temporary tendency arrays                    !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u_2, v_2, w_2, t_2, moist, rho, q_sq_2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp
      real, dimension (:,:,:), allocatable, intent (inout)      :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                   t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                                     intent (in)  :: ide, jde, kde,                &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte

      call Save_copy_state (u_2, v_2, w_2, t_2, moist, rho, q_sq_2, u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, &
              rho_tmp, q_sq_tmp, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Allocate_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
              t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, its, ite, jts, jte, kts, kte)

      call Zero_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              ide, jde, kde, its, ite, jts, jte, kts, kte)

    end subroutine Init_temp_arrays_for_substep


    subroutine Save_copy_state (u_2, v_2, w_2, t_2, moist, rho, q_sq_2, u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, &
              rho_tmp, q_sq_tmp, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Save copy of state arrays for sub-stepping                   !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u_2, v_2, w_2, t_2, moist, rho, q_sq_2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp
      integer,                                     intent (in)  :: ide, jde, kde,                &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte

        ! Local
      integer :: i, j, k, i_start, i_end, j_start, j_end, k_start, k_end

        ! Calculate indexes
      i_start = its - 1
      i_end = ite + 1
      j_start = jts - 1
      j_end = jte + 1
      k_start = kts

      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            u_tmp(i,k,j) = u_2(i,k,j)
          ENDDO
        ENDDO
      ENDDO

      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            v_tmp(i,k,j) = v_2(i,k,j)
          ENDDO
        ENDDO
      ENDDO

      k_end = kte

      DO j = j_start - 1, j_end + 1
        DO k = k_start, k_end
          DO i = i_start - 1, i_end + 1
            w_tmp(i,k,j) = w_2(i,k,j)
          ENDDO
        ENDDO
      ENDDO

      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            th_tmp(i,k,j) = t_2(i,k,j)
            qv_tmp(i,k,j) = moist(i,k,j)
            rho_tmp(i,k,j) = rho(i,k,j)
            q_sq_tmp(i,k,j) = q_sq_2(i,k,j)
          ENDDO
        ENDDO
      ENDDO

    end subroutine Save_copy_state


    subroutine Allocate_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
              t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Allocate temporary tendencies for sub-stepping               !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (:,:,:), allocatable, intent (inout)     :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                  t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                              intent (in)        :: its, ite, jts, jte, kts, kte

      IF (.NOT. ALLOCATED(ru_tendf_tmp))   ALLOCATE ( ru_tendf_tmp(its:ite,kts:kte,jts:jte)   )
      IF (.NOT. ALLOCATED(rv_tendf_tmp))   ALLOCATE ( rv_tendf_tmp(its:ite,kts:kte,jts:jte)   )
      IF (.NOT. ALLOCATED(rw_tendf_tmp))   ALLOCATE ( rw_tendf_tmp(its:ite,kts:kte,jts:jte)   )
      IF (.NOT. ALLOCATED(t_tendf_tmp))    ALLOCATE ( t_tendf_tmp(its:ite,kts:kte,jts:jte)    )
      IF (.NOT. ALLOCATED(moist_tend_tmp)) ALLOCATE ( moist_tend_tmp(its:ite,kts:kte,jts:jte) )
      IF (.NOT. ALLOCATED(q_sq_tend_tmp))  ALLOCATE ( q_sq_tend_tmp(its:ite,kts:kte,jts:jte) )

    end subroutine Allocate_temp_tends


    subroutine Deallocate_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
              t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Deallocate temporary tendencies for sub-stepping             !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (:,:,:), allocatable, intent (inout)     :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                  t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                              intent (in)        :: its, ite, jts, jte, kts, kte

      IF (ALLOCATED(ru_tendf_tmp))    DEALLOCATE ( ru_tendf_tmp   )
      IF (ALLOCATED(rv_tendf_tmp))    DEALLOCATE ( rv_tendf_tmp   )
      IF (ALLOCATED(rw_tendf_tmp))    DEALLOCATE ( rw_tendf_tmp   )
      IF (ALLOCATED(t_tendf_tmp))     DEALLOCATE ( t_tendf_tmp    )
      IF (ALLOCATED(moist_tend_tmp))  DEALLOCATE ( moist_tend_tmp )
      IF (ALLOCATED(q_sq_tend_tmp))   DEALLOCATE ( q_sq_tend_tmp  )

    end subroutine Deallocate_temp_tends


    subroutine Zero_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
              t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, ide, jde, kde, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Zero out temporary tendencies for sub-stepping               !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (its:ite, kts:kte, jts:jte), intent (out)  :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                    t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                              intent (in)          :: ide, jde, kde, its, ite, jts, jte, kts, kte

        ! Local
      integer :: i, j, k, i_start, i_end, j_start, j_end, k_start, k_end

        ! Calculate indexes for u
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            ru_tendf_tmp(i,k,j) = 0.0
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for v
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = jte
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            rv_tendf_tmp(i,k,j) = 0.0
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for w
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = kte

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            rw_tendf_tmp(i,k,j) = 0.0
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for theta/qv
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            t_tendf_tmp(i,k,j) = 0.0
            moist_tend_tmp(i,k,j) = 0.0
            q_sq_tend_tmp(i,k,j) = 0.0
          ENDDO
        ENDDO
      ENDDO

    end subroutine Zero_temp_tends


    subroutine Update_wrf_tends_temp_state_and_zero_tends (ru_tendf, rv_tendf, rw_tendf, t_tendf, &
              moist_tend, q_sq_tend, ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp, dt, c1h, c2h, c1f, c2f, muu, muv, mut, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                          !
      ! Purpose: Update arrays for sub-stepping                                  !
      !              1. Update WRF tendencies with temporary tendencies          !
      !              2. Update temporary state arrays with temporary tendencies  !
      !              3. Zero out temporary tendency arrays                       !
      !                                                                          !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                            !
      !                                                                          !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout)  :: ru_tendf, rv_tendf, rw_tendf, t_tendf, moist_tend, q_sq_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout)  :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                      t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout)  :: u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp
      real,                                        intent (in)     :: dt
      real, dimension (kms:kme),                   intent (in)     :: c1h, c2h, c1f, c2f
      real, dimension (ims:ime, jms:jme),          intent (in)     :: muu, muv, mut
      integer,                                     intent (in)     :: ide, jde, kde,                &
                                                                      ims, ime, jms, jme, kms, kme, &
                                                                      its, ite, jts, jte, kts, kte

      call Update_wrf_tends (ru_tendf, rv_tendf, rw_tendf, t_tendf, moist_tend, q_sq_tend, &
          ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Update_temp_state (u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp, &
          ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
          dt, c1h, c2h, c1f, c2f, muu, muv, mut, ide, jde, kde, ims, ime, jms, jme, kms, kme,  &
          its, ite, jts, jte, kts, kte)

      call Zero_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, &
          q_sq_tend_tmp, ide, jde, kde, its, ite, jts, jte, kts, kte)

    end subroutine Update_wrf_tends_temp_state_and_zero_tends


    subroutine Update_wrf_tends (ru_tendf, rv_tendf, rw_tendf, t_tendf, moist_tend, q_sq_tend, &
              ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Update WRF tendency arrays with temporary tendencies         !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout)  :: ru_tendf, rv_tendf, rw_tendf, t_tendf, moist_tend, q_sq_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)     :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                      t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                                     intent (in)     :: ide, jde, kde,                &
                                                                      ims, ime, jms, jme, kms, kme, &
                                                                      its, ite, jts, jte, kts, kte

        ! Local
      integer :: i, j, k, i_start, i_end, j_start, j_end, k_start, k_end

        ! Calculate indexes for u
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            ru_tendf(i,k,j) = ru_tendf(i,k,j) + ru_tendf_tmp(i,k,j)
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for v
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = jte
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            rv_tendf(i,k,j) = rv_tendf(i,k,j) + rv_tendf_tmp(i,k,j)
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for w
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = kte

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            rw_tendf(i,k,j) = rw_tendf(i,k,j) + rw_tendf_tmp(i,k,j)
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for theta/qv
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            t_tendf(i,k,j) = t_tendf(i,k,j) + t_tendf_tmp(i,k,j)
            moist_tend(i,k,j) = moist_tend(i,k,j) + moist_tend_tmp(i,k,j)
            q_sq_tend(i,k,j) = q_sq_tend(i,k,j) + q_sq_tend_tmp(i,k,j)
          ENDDO
        ENDDO
      ENDDO

    end subroutine Update_wrf_tends


    subroutine Update_temp_state (u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp, &
              ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              dt, c1h, c2h, c1f, c2f, muu, muv, mut, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Updates temporary state arrays for sub-stepping              !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout)  :: u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)     :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                      t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      real,                                        intent (in)     :: dt
      real, dimension (kms:kme),                   intent (in)     :: c1h, c2h, c1f, c2f
      real, dimension (ims:ime, jms:jme),          intent (in)     :: muu, muv, mut 
      integer,                                     intent (in)     :: ide, jde, kde, &
                                                                      ims, ime, jms, jme, kms, kme, &
                                                                      its, ite, jts, jte, kts, kte

        ! Local
      integer :: i, j, k, i_start, i_end, j_start, j_end, k_start, k_end

        ! Calculate indexes for u
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            u_tmp(i,k,j) = u_tmp(i,k,j) + (ru_tendf_tmp(i,k,j) * (dt / (c1h(k) * muu(i,j) + c2h(k))))
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for v
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = jte
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            v_tmp(i,k,j) = v_tmp(i,k,j) + (rv_tendf_tmp(i,k,j) * (dt / (c1h(k) * muv(i,j) + c2h(k))))
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for w
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = kte

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            w_tmp(i,k,j) = w_tmp(i,k,j) + (rw_tendf_tmp(i,k,j) * (dt / (c1f(k) * mut(i,j) + c2f(k))))
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for theta/qv
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            th_tmp(i,k,j) = th_tmp(i,k,j) + (t_tendf_tmp(i,k,j) * (dt / (c1h(k) * mut(i,j) + c2h(k))))
            qv_tmp(i,k,j) = qv_tmp(i,k,j) + (moist_tend_tmp(i,k,j) * (dt / (c1h(k) * mut(i,j) + c2h(k))))
            rho_tmp(i,k,j) = rho_tmp(i,k,j)
            q_sq_tmp(i,k,j) = max ( Q_SQ_MIN , q_sq_tmp(i,k,j) + (q_sq_tend_tmp(i,k,j) * (dt / (c1h(k) * mut(i,j) + c2h(k)))) )
          ENDDO
        ENDDO
      ENDDO

    end subroutine Update_temp_state


    subroutine Calc_turb_fluxes (config_flags, dz, z_at_mass, u, v, w, th, th2, qv, thetav, tsk, t2, q2, psfc, rho, ust, hfx, qfx, rmol, &
              du_dx_faces, du_dy_faces, du_dz_faces, &
              dv_dx_faces, dv_dy_faces, dv_dz_faces, &
              dw_dx_faces, dw_dy_faces, dw_dz_faces, &
              dqv_dx_faces, dqv_dy_faces, dqv_dz_faces, &
              dthetav_dx_faces, dthetav_dy_faces, dthetav_dz_faces, &
              du_dx_mass, du_dy_mass, du_dz_mass, &
              dv_dx_mass, dv_dy_mass, dv_dz_mass, &
              dw_dx_mass, dw_dy_mass, dw_dz_mass, &
              dthetav_dx_mass, dthetav_dy_mass, dthetav_dz_mass, &
              dqv_dx_mass, dqv_dy_mass, dqv_dz_mass, &
              itimestep, ss, n_tracer, &
              q_sq, q_sq_pa, q_sq_rat, q_sq_prog, q_sq_tend, &
              q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, q_sq_buoyancy_tend, q_sq_dissip_tend, &
              khx, khy, khz, l_master, l_master_pa, l_master_rat, l_master_at_mass, turb_flux_u2, turb_flux_v2, turb_flux_w2,   &
              turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
              turb_flux_theta2_v, &
              turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, &
              turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, &
              turb_flux_u2_pa, turb_flux_v2_pa, turb_flux_w2_pa,   &
              turb_flux_uv_pa, turb_flux_uw_pa, turb_flux_vw_pa, turb_flux_utheta_v_pa, turb_flux_vtheta_v_pa, turb_flux_wtheta_v_pa, &
              turb_flux_theta2_v_pa, &
              turb_flux_uqv_pa, turb_flux_vqv_pa, turb_flux_wqv_pa, &
              turb_flux_utheta_pa, turb_flux_vtheta_pa, turb_flux_wtheta_pa, &
              turb_flux_u2_rat, turb_flux_v2_rat, turb_flux_w2_rat,   &
              turb_flux_uv_rat, turb_flux_uw_rat, turb_flux_vw_rat, turb_flux_utheta_v_rat, turb_flux_vtheta_v_rat, &
              turb_flux_wtheta_v_rat, turb_flux_theta2_v_rat, &
              turb_flux_uqv_rat, turb_flux_vqv_rat, turb_flux_wqv_rat, &
              turb_flux_utheta_rat, turb_flux_vtheta_rat, turb_flux_wtheta_rat, mat_cond_heat, mat_cond_moist, &
              turb_flux_u2_mass, turb_flux_v2_mass, turb_flux_w2_mass, turb_flux_uv_mass, turb_flux_uw_mass, turb_flux_vw_mass, &
              turb_flux_wtheta_mass, turb_flux_wtheta_v_mass, msfux, msfuy, msfvx, msfvy, msftx, msfty, &
              dx, dy, rdx, rdy, dn, dnw, rdz, dt, id, &
              rdzw, fnm, fnp, cf1, cf2, cf3, zx, zy, mut, c1h, c2h, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its,    &
              ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Calculates turbulent fluxes                                  !
      !                                                                       !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                         !
      !                                                                       !
      ! Method: Mellor and Yamada model                                       !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags
      integer, intent(in) :: itimestep, id, ss, n_tracer

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: dz, z_at_mass, u, v, &
          w, th, th2, qv, rho, zx, zy, rdz, rdzw
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq, q_sq_pa
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq_prog
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_vdiff_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_hdiff_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_shear_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_buoyancy_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_dissip_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: khx, khy, khz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: thetav
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: du_dx_faces, du_dy_faces, du_dz_faces, &
          dv_dx_faces, dv_dy_faces, dv_dz_faces, dw_dx_faces, dw_dy_faces, dw_dz_faces, dthetav_dx_faces,   &
          dthetav_dy_faces, dthetav_dz_faces, dqv_dx_faces, dqv_dy_faces, dqv_dz_faces
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: du_dx_mass, du_dy_mass, du_dz_mass, &
          dv_dx_mass, dv_dy_mass, dv_dz_mass, dw_dx_mass, dw_dy_mass, dw_dz_mass, &
          dthetav_dx_mass, dthetav_dy_mass, dthetav_dz_mass, dqv_dx_mass, dqv_dy_mass, dqv_dz_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: l_master, l_master_pa, l_master_at_mass
      real, dimension (ims:ime, jms:jme),          intent (in)    :: msfux, msfuy, msfvx, msfvy, msftx, msfty
      real, dimension (ims:ime, jms:jme),          intent (in)    :: tsk, t2, q2, psfc, ust, hfx, qfx, rmol
      real, dimension (ims:ime, jms:jme),          intent (in)    :: mut
      real, dimension (kms:kme),                   intent (in)    :: c1h, c2h
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                                                                     turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                                                                     turb_flux_utheta_v, turb_flux_vtheta_v,   &
                                                                     turb_flux_wtheta_v, turb_flux_theta2_v,   &
                                                                     turb_flux_uqv, turb_flux_vqv,             &
                                                                     turb_flux_wqv, turb_flux_utheta,          &
                                                                     turb_flux_vtheta, turb_flux_wtheta
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: turb_flux_u2_pa, turb_flux_v2_pa, turb_flux_w2_pa, &
                                                                     turb_flux_uv_pa, turb_flux_uw_pa, turb_flux_vw_pa, &
                                                                     turb_flux_utheta_v_pa, turb_flux_vtheta_v_pa,   &
                                                                     turb_flux_wtheta_v_pa, turb_flux_theta2_v_pa,   &
                                                                     turb_flux_uqv_pa, turb_flux_vqv_pa,             &
                                                                     turb_flux_wqv_pa, turb_flux_utheta_pa,          &
                                                                     turb_flux_vtheta_pa, turb_flux_wtheta_pa
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: turb_flux_u2_rat, turb_flux_v2_rat, turb_flux_w2_rat, &
                                                                     turb_flux_uv_rat, turb_flux_uw_rat, turb_flux_vw_rat, &
                                                                     turb_flux_utheta_v_rat, turb_flux_vtheta_v_rat,   &
                                                                     turb_flux_wtheta_v_rat, turb_flux_theta2_v_rat,   &
                                                                     turb_flux_uqv_rat, turb_flux_vqv_rat,             &
                                                                     turb_flux_wqv_rat, turb_flux_utheta_rat,          &
                                                                     turb_flux_vtheta_rat, turb_flux_wtheta_rat, q_sq_rat, l_master_rat
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: mat_cond_heat, mat_cond_moist
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: turb_flux_u2_mass, turb_flux_v2_mass, turb_flux_w2_mass, &
                                                                     turb_flux_uv_mass, turb_flux_uw_mass, turb_flux_vw_mass, &
                                                                     turb_flux_wtheta_mass, turb_flux_wtheta_v_mass
      real, dimension (kms:kme),                   intent (in)    :: fnm, fnp, dn, dnw
      real,                                        intent (in)    :: dx, dy, rdx, rdy, cf1, cf2, cf3, dt
      integer,                                     intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                                                     ims, ime, jms, jme, kms, kme, &
                                                                     its, ite, jts, jte, kts, kte

        ! Local vars
!      real, dimension (its:ite, kts:kte, jts:jte) :: sm_l2, sh_l2, sm, sh
      real, dimension (its:ite, kts:kte, jts:jte) :: du_dx, du_dy, du_dz, &
                                                     dv_dx, dv_dy, dv_dz, &
                                                     dw_dx, dw_dy, dw_dz, &
                                                     dthetav_dx, dthetav_dy, dthetav_dz, &
                                                     dqv_dx, dqv_dy, dqv_dz
!      real, dimension ((ite-its+1)*(kte-kts+1)*(jte-jts+1)) :: mat_heat_hold, mat_moist_hold
!      real, allocatable :: mat_heat_pos(:), mat_moist_pos(:)
      logical, parameter :: DEBUG = .false.
      logical, parameter :: DO_FLUX_RATIO = .true.
      real, parameter :: SMALL = 1.E-18
      real, parameter :: min_rat = 0.0
      real, parameter :: max_rat = 2.0
      integer :: pbl_approx_opt, i, j, k, i_start, i_end, j_start, j_end, ktf
 !     integer :: mat_heat_count, mat_moist_count
      real :: mean_heat, mean_moist, median_heat, median_moist

        ! Calculate derivatives at the faces to prepare for fluxes
      call Calc_wind_derivatives (config_flags, u, v, w, du_dx, dv_dy, dw_dz,  &
          du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz,                            &
          msftx, msfty, msfux, msfuy, msfvx, msfvy,                            &
          rdx, rdy, dn, dnw, dz, rdz,                                          &
          rdzw, fnm, fnp, cf1, cf2, cf3, zx, zy,                               &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its,     &
          ite, jts, jte, kts, kte)

      call Calc_thetav (th2, qv, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, thetav)

      call Calc_mass_var_derivatives (thetav, qv,                      &
          msftx, msfty, msfux, msfuy, msfvx, msfvy, rdx, rdy,          &
          fnm, fnp, cf1, cf2, cf3, zx, zy, rdz, rdzw, dnw, dn,         &
          dthetav_dx, dthetav_dy, dthetav_dz,                          &
          dqv_dx, dqv_dy, dqv_dz,                                      &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme,  &
          its, ite, jts, jte, kts, kte)

        ! Calculate indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (DEBUG .and. itimestep == 1) then
        print *, 'its = ', its, 'ite = ', ite
        print *, 'ims = ', ims, 'ime = ', ime
        print *, 'jts = ', jts, 'jte = ', jte
        print *, 'jms = ', jms, 'jme = ', jme
      end if

        ! Save derivatives
      do j = j_start, j_end
        do k = kts, ktf + 1
          do i = i_start, i_end
           du_dx_faces(i, k, j) = du_dx(i, k, j)
           du_dy_faces(i, k, j) = du_dy(i, k, j)
           du_dz_faces(i, k, j) = du_dz(i, k, j)
           dv_dx_faces(i, k, j) = dv_dx(i, k, j)
           dv_dy_faces(i, k, j) = dv_dy(i, k, j)
           dv_dz_faces(i, k, j) = dv_dz(i, k, j)
           dw_dx_faces(i, k, j) = dw_dx(i, k, j)
           dw_dy_faces(i, k, j) = dw_dy(i, k, j)
           dw_dz_faces(i, k, j) = dw_dz(i, k, j)
           dqv_dx_faces(i, k, j) = dqv_dx(i, k, j)
           dqv_dy_faces(i, k, j) = dqv_dy(i, k, j)
           dqv_dz_faces(i, k, j) = dqv_dz(i, k, j)
           dthetav_dx_faces(i, k, j) = dthetav_dx(i, k, j)
           dthetav_dy_faces(i, k, j) = dthetav_dy(i, k, j)
           dthetav_dz_faces(i, k, j) = dthetav_dz(i, k, j)
          end do
        end do
      end do

      if (config_flags%pbl3d_opt < 2) then  ! PBL approx analytical solution only
        call Calc_fluxes_pbl_approx (config_flags, l_master, du_dz, dv_dz, dthetav_dz, dqv_dz, &
            th2, qv, dn, dnw, fnm, fnp, dx, dy, dz, zx, zy, msftx, msfty, rdz, rdzw, cf1, cf2, cf3, dt, id, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts,    &
            jte, kts, kte, q_sq, q_sq_prog, turb_flux_u2, turb_flux_v2, turb_flux_w2,     &
            turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                 &
            turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_uqv,    &
            turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, &
            u, v, w, thetav, rho, ust, hfx, qfx, rmol)
      else if (config_flags%pbl3d_opt == 2) then  ! PBL approx analytical + PBL approx numerical solution
        call Calc_fluxes_pbl_approx (config_flags, l_master_pa, du_dz, dv_dz, dthetav_dz, dqv_dz, &
            th2, qv, dn, dnw, fnm, fnp, dx, dy, dz, zx, zy, msftx, msfty, rdz, rdzw, cf1, cf2, cf3, dt, id, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts,    &
            jte, kts, kte, q_sq_pa, q_sq_prog, turb_flux_u2_pa, turb_flux_v2_pa, turb_flux_w2_pa,     &
            turb_flux_uv_pa, turb_flux_uw_pa, turb_flux_vw_pa, turb_flux_utheta_v_pa,                 &
            turb_flux_vtheta_v_pa, turb_flux_wtheta_v_pa, turb_flux_theta2_v_pa, turb_flux_uqv_pa,    &
            turb_flux_vqv_pa, turb_flux_wqv_pa, turb_flux_utheta_pa, turb_flux_vtheta_pa, turb_flux_wtheta_pa, &
            u, v, w, thetav, rho, ust, hfx, qfx, rmol)

        pbl_approx_opt = 1
        call Calc_fluxes (config_flags, pbl_approx_opt, id, config_flags%specified, q_sq, q_sq_prog, l_master, du_dx, du_dy, du_dz, dv_dx, dv_dy,   &
            dv_dz, dw_dx, dw_dy, dw_dz, dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy,  &
            dqv_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw,    &
            turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v,        &
            turb_flux_theta2_v, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,        &
            turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, mat_cond_heat, mat_cond_moist, th2, thetav, qv, rho, hfx, qfx, rmol, u, v, ust,  &
            dz, rdzw, fnm, fnp, cf1, cf2, cf3, dn, dnw, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme,  &
            its, ite, jts, jte, kts, kte)
      else if (config_flags%pbl3d_opt == 3) then  ! PBL approx analytical + full 3D numerical solution
        if (config_flags%replace_with_pa > 0) then
          call Calc_fluxes_pbl_approx (config_flags, l_master_pa, du_dz, dv_dz, dthetav_dz, dqv_dz, &
              th2, qv, dn, dnw, fnm, fnp, dx, dy, dz, zx, zy, msftx, msfty, rdz, rdzw, cf1, cf2, cf3, dt, id, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts,    &
              jte, kts, kte, q_sq_pa, q_sq_prog, turb_flux_u2_pa, turb_flux_v2_pa, turb_flux_w2_pa,     &
              turb_flux_uv_pa, turb_flux_uw_pa, turb_flux_vw_pa, turb_flux_utheta_v_pa,                 &
              turb_flux_vtheta_v_pa, turb_flux_wtheta_v_pa, turb_flux_theta2_v_pa, turb_flux_uqv_pa,    &
              turb_flux_vqv_pa, turb_flux_wqv_pa, turb_flux_utheta_pa, turb_flux_vtheta_pa, turb_flux_wtheta_pa, &
              u, v, w, thetav, rho, ust, hfx, qfx, rmol)
        end if

        pbl_approx_opt = 0
        call Calc_fluxes (config_flags, pbl_approx_opt, id, config_flags%specified, q_sq, q_sq_prog, l_master, du_dx, du_dy, du_dz, dv_dx, dv_dy,   &
            dv_dz, dw_dx, dw_dy, dw_dz, dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy,  &
            dqv_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw,    &
            turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v,        &
            turb_flux_theta2_v, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,        &
            turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, mat_cond_heat, mat_cond_moist, th2, thetav, qv, rho, hfx, qfx, rmol, u, v, ust,  &
            dz, rdzw, fnm, fnp, cf1, cf2, cf3, dn, dnw, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme,  &
            its, ite, jts, jte, kts, kte)
      end if

      if (config_flags%pbl3d_opt == 3 .and. config_flags%replace_with_pa > 0) then
        do j = j_start, j_end
          do k = kts + 1, ktf
            do i = i_start, i_end
              if (Sqrt(q_sq(i, k, j)) / l_master(i, k, j) < config_flags%replace_with_pa_rat) then
                turb_flux_u2_pa(i, k, j) = turb_flux_u2(i, k, j)
                turb_flux_v2_pa(i, k, j) = turb_flux_v2(i, k, j)
                turb_flux_w2_pa(i, k, j) = turb_flux_w2(i, k, j)
                turb_flux_uv_pa(i, k, j) = turb_flux_uv(i, k, j)
                turb_flux_uw_pa(i, k, j) = turb_flux_uw(i, k, j)
                turb_flux_vw_pa(i, k, j) = turb_flux_vw(i, k, j)
                turb_flux_utheta_pa(i, k, j) = turb_flux_utheta(i, k, j)
                turb_flux_vtheta_pa(i, k, j) = turb_flux_vtheta(i, k, j)
                turb_flux_wtheta_pa(i, k, j) = turb_flux_wtheta(i, k, j)
                turb_flux_utheta_v_pa(i, k, j) = turb_flux_utheta_v(i, k, j)
                turb_flux_vtheta_v_pa(i, k, j) = turb_flux_vtheta_v(i, k, j)
                turb_flux_wtheta_v_pa(i, k, j) = turb_flux_wtheta_v(i, k, j)
                turb_flux_uqv_pa(i, k, j) = turb_flux_uqv(i, k, j)
                turb_flux_vqv_pa(i, k, j) = turb_flux_vqv(i, k, j)
                turb_flux_wqv_pa(i, k, j) = turb_flux_wqv(i, k, j)
              end if
            end do
          end do
        end do
      end if

      if (DO_FLUX_RATIO) then
      ! Calculate flux ratio
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            q_sq_rat(i, k, j) = q_sq(i, k, j) / SIGN(MAX(ABS(q_sq_pa(i, k, j)),SMALL),q_sq_pa(i, k, j))
            l_master_rat(i, k, j) = l_master(i, k, j) / SIGN(MAX(ABS(l_master_pa(i, k, j)),SMALL),l_master_pa(i, k, j))
            turb_flux_u2_rat(i, k, j) = turb_flux_u2(i, k, j) / SIGN(MAX(ABS(turb_flux_u2_pa(i, k, j)),SMALL),turb_flux_u2_pa(i, k, j))
            turb_flux_v2_rat(i, k, j) = turb_flux_v2(i, k, j) / SIGN(MAX(ABS(turb_flux_v2_pa(i, k, j)),SMALL),turb_flux_v2_pa(i, k, j))
            turb_flux_w2_rat(i, k, j) = turb_flux_w2(i, k, j) / SIGN(MAX(ABS(turb_flux_w2_pa(i, k, j)),SMALL),turb_flux_w2_pa(i, k, j))
            turb_flux_uv_rat(i, k, j) = turb_flux_uv(i, k, j) / SIGN(MAX(ABS(turb_flux_uv_pa(i, k, j)),SMALL),turb_flux_uv_pa(i, k, j))
            turb_flux_uw_rat(i, k, j) = turb_flux_uw(i, k, j) / SIGN(MAX(ABS(turb_flux_uw_pa(i, k, j)),SMALL),turb_flux_uw_pa(i, k, j))
            turb_flux_vw_rat(i, k, j) = turb_flux_vw(i, k, j) / SIGN(MAX(ABS(turb_flux_vw_pa(i, k, j)),SMALL),turb_flux_vw_pa(i, k, j))
            turb_flux_utheta_v_rat(i, k, j) = turb_flux_utheta_v(i, k, j) / SIGN(MAX(ABS(turb_flux_utheta_v_pa(i, k, j)),SMALL),turb_flux_utheta_v_pa(i, k, j))
            turb_flux_vtheta_v_rat(i, k, j) = turb_flux_vtheta_v(i, k, j) / SIGN(MAX(ABS(turb_flux_vtheta_v_pa(i, k, j)),SMALL),turb_flux_vtheta_v_pa(i, k, j))
            turb_flux_wtheta_v_rat(i, k, j) = turb_flux_wtheta_v(i, k, j) / SIGN(MAX(ABS(turb_flux_wtheta_v_pa(i, k, j)),SMALL),turb_flux_wtheta_v_pa(i, k, j))
            turb_flux_theta2_v_rat(i, k, j) = turb_flux_theta2_v(i, k, j) / SIGN(MAX(ABS(turb_flux_theta2_v_pa(i, k, j)),SMALL),turb_flux_theta2_v_pa(i, k, j))
            turb_flux_utheta_rat(i, k, j) = turb_flux_utheta(i, k, j) / SIGN(MAX(ABS(turb_flux_utheta_pa(i, k, j)),SMALL),turb_flux_utheta_pa(i, k, j))
            turb_flux_vtheta_rat(i, k, j) = turb_flux_vtheta(i, k, j) / SIGN(MAX(ABS(turb_flux_vtheta_pa(i, k, j)),SMALL),turb_flux_vtheta_pa(i, k, j))
            turb_flux_wtheta_rat(i, k, j) = turb_flux_wtheta(i, k, j) / SIGN(MAX(ABS(turb_flux_wtheta_pa(i, k, j)),SMALL),turb_flux_wtheta_pa(i, k, j))
            turb_flux_uqv_rat(i, k, j) = turb_flux_uqv(i, k, j) / SIGN(MAX(ABS(turb_flux_uqv_pa(i, k, j)),SMALL),turb_flux_uqv_pa(i, k, j))
            turb_flux_vqv_rat(i, k, j) = turb_flux_vqv(i, k, j) / SIGN(MAX(ABS(turb_flux_vqv_pa(i, k, j)),SMALL),turb_flux_vqv_pa(i, k, j))
            turb_flux_wqv_rat(i, k, j) = turb_flux_wqv(i, k, j) / SIGN(MAX(ABS(turb_flux_wqv_pa(i, k, j)),SMALL),turb_flux_wqv_pa(i, k, j))
          end do
        end do
      end do
      end if

        ! Calculate eddy diffusivities for tracer mixing
      if ( ss .eq. config_flags%pbl3d_nsteps .and. n_tracer .ge. param_first_scalar ) then
        call Calc_kh (du_dz, dv_dz, dthetav_dz, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, &
            turb_flux_vtheta_v, turb_flux_wtheta_v, l_master, q_sq, khx, khy, khz, ids, ide, &
            jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

        ! Predict q_sq if level 2.5 model
      if ( config_flags%pbl3d_prog .gt. 0 ) then   ! Prognostic TKE
        call Calc_q_sq_l2p5 (config_flags, u, v, w, zx, zy, th2, qv, rho, thetav, tsk, t2, q2, psfc, ust, hfx, q_sq, q_sq_prog, q_sq_tend, &
            q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, q_sq_buoyancy_tend, q_sq_dissip_tend, l_master, l_master_at_mass, &
            turb_flux_u2, turb_flux_v2, turb_flux_w2, &
            turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta, turb_flux_wtheta_v, &
            turb_flux_u2_mass, turb_flux_v2_mass, turb_flux_w2_mass, turb_flux_uv_mass, turb_flux_uw_mass, turb_flux_vw_mass, turb_flux_wtheta_mass, &
            turb_flux_wtheta_v_mass, dn, dnw, fnm, fnp, msftx, msfty, msfux, msfuy, msfvx, msfvy, mut, c1h, c2h, &
            z_at_mass, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, dt, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, &
            du_dx_mass, du_dy_mass, du_dz_mass, &
            dv_dx_mass, dv_dy_mass, dv_dz_mass, &
            dw_dx_mass, dw_dy_mass, dw_dz_mass, &
            dthetav_dx_mass, dthetav_dy_mass, dthetav_dz_mass, &
            dqv_dx_mass, dqv_dy_mass, dqv_dz_mass)
      end if

      if (DEBUG) print *, 'Done Calc_turb_fluxes...'

    end subroutine Calc_turb_fluxes


    subroutine Calc_wind_derivatives (config_flags, u, v, w, du_dx, dv_dy, dw_dz, du_dy, dv_dx, &
        dw_dx, du_dz, dw_dy, dv_dz, &
        msftx, msfty, msfux, msfuy, msfvx, msfvy, &
        rdx, rdy, dn, dnw, dz, rdz,   &
        rdzw, fnm, fnp, cf1, cf2, cf3, zx, zy, &
        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its,    &
        ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                                         !
        ! Purpose: Calculate spatial derivatives of the wind variables at the center of the faces !
        !                                                                                         !
        ! Authors: Pedro A. Jimenez & Timothy W. Juliano                                          !
        !                                                                                         !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, v, w, zx, zy, dz, rdz, rdzw
      real, dimension (kms:kme), intent (in)                    :: fnm, fnp, dn, dnw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfty, msfux, msfuy, msfvx, msfvy
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dx, dv_dy, dw_dz, du_dy, dv_dx, &
                                                                   dw_dx, du_dz, dw_dy, dv_dz

      logical, parameter :: DEBUG = .false.

      call Calc_du_dz (u, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      call Calc_dv_dz (v, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      call Calc_dw_dz (w, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dw_dz)

      call Calc_du_dx (u, zx, rdz, dn, dnw, fnm, fnp, msftx, msfux, rdx, cf1, cf2, cf3,     &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dx  )

      call Calc_dv_dx (v, zx, rdz, dn, dnw, fnm, fnp, msftx, msfvx, rdx, cf1, cf2, cf3, ide, jde,  &
          kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dx     )

      call Calc_dw_dx (w, zx, rdz, msftx, rdx, ims, ime, jms, jme, kms, kme,  &
          ide, jde, kde, its, ite, jts, jte, kts, kte, config_flags, dw_dx    )

      call Calc_du_dy (u, zy, rdz, dn, dnw, fnm, fnp, msfty, msfuy, rdy, cf1, cf2, cf3,     &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dy  )

      call Calc_dv_dy (v, zy, rdz, dn, dnw, fnm, fnp, msfty, msfvy, rdy, cf1, cf2, cf3, ide, jds,    &
          jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dy  )

      call Calc_dw_dy (w, zy, rdz, msfty, rdy, ims, ime, jms, jme, kms, kme,  &
          ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dy                  )

      if (DEBUG) print *, 'Done Calc_wind_derivatives...'

    end subroutine Calc_wind_derivatives


    subroutine Calc_du_dz (u, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dz                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, rdz
      real, dimension (kms:kme),                   intent (in)  :: dnw
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dz

        ! Local vars
      integer :: i, j, k
      real, dimension (its:ite, kts:kte, jts:jte) :: u_at_mass
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Unstagger U
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            u_at_mass(i, k, j) = 0.5 * (u(i, k, j) + u(i + 1, k, j)) 
          end do
        end do
      end do

        ! Calculate derivative at the center of the faces
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            du_dz(i, k, j) = (u_at_mass(i, k, j) - u_at_mass(i, k - 1, j)) * rdz(i, k, j) 
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          du_dz(i, kts, j) = 0.0
            ! Top
          du_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dz


    subroutine Calc_dv_dz (v, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dz                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, rdz
      real, dimension (kms:kme),                   intent (in)  :: dnw
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dz

        ! Local vars
      integer :: i, j, k
      real, dimension (its:ite, kts:kte, jts:jte) :: v_at_mass
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Unstagger V
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            v_at_mass(i, k, j) = 0.5 * (v(i, k, j) + v(i, k, j + 1))
          end do
        end do
      end do

        ! Calculate derivative at the center of the faces
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dv_dz(i, k, j) = (v_at_mass(i, k, j) - v_at_mass(i, k - 1, j)) * rdz(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          dv_dz(i, kts, j) = 0.0
            ! Top face
          dv_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dv_dz


    subroutine Calc_dw_dz (w, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dw_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dz                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, rdz
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dz

        ! Local vars
      integer :: i, j, k
      real, dimension (its:ite, kts:kte, jts:jte) :: w_at_mass
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Calculate w at mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            w_at_mass(i, k, j) = 0.5 * (w(i, k, j) + w(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate partial w / partial z
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dw_dz(i, k, j) = (w_at_mass(i, k, j) - w_at_mass(i, k - 1, j)) * rdz(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          dw_dz(i, kts, j) =  0.0
            ! Top face
          dw_dz(i, kte, j) =  0.0
        end do
      end do

    end subroutine Calc_dw_dz


    subroutine Calc_du_dx (u, zx, rdz, dn, dnw, fnm, fnp, msftx, msfux, rdx, cf1, cf2, cf3,         &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dx      )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dx                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, zx, rdz
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfux
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dx

        ! Local vars
      integer :: i, j, k
      real    :: dz_dx, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_tx_ux
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: du_dz_dz_dx, u_mux, u_mux_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at u points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_tx_ux(i, j) = msftx(i, j) * msfux(i, j)
        end do
      end do

        ! Scale u component of vector to account for map projection
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            u_mux(i, k, j) = u(i, k, j) / msfux(i, j)
          end do
        end do
      end do

        ! Interpolate u_mux from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            u_mux_at_full(i, k, j) = fnm(k) * u_mux(i, k, j) + fnp(k) * u_mux(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          u_mux_at_full(i, 1, j) =  cf1 * u_mux(i, 1, j) + cf2 * u_mux(i, 2, j) + cf3 * u_mux(i, 3, j) 
            ! Top face
          u_mux_at_full(i, kte, j) =  cft1 *  u_mux(i, ktes1, j) +  cft2 * u_mux(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (du_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at center of faces
            dz_dx = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
              ! du_dz_dz_dx = (partial u_mux / partial z)|sigma (du_dz) * 
              ! (partial z / partial x)|sigma (dz_dx)
            du_dz_dz_dx(i, k, j) = (0.5 * (u_mux(i, k, j) + u_mux(i + 1, k, j)) - &
                0.5 * (u_mux(i, k - 1, j) + u_mux(i + 1, k - 1, j))) * rdz(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial u / partial x)|z = msf_tx_ux * (partial u_mux_at_full / partial x)|sigma - du_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            du_dx(i, k, j) = msf_tx_ux(i, j) * ((u_mux_at_full(i + 1, k, j) - u_mux_at_full(i, k, j)) * rdx -  &
                du_dz_dz_dx(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          du_dx(i, kts, j) = 0.0
            ! Top face
          du_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dx


    subroutine Calc_dv_dx (v, zx, rdz, dn, dnw, fnm, fnp, msftx, msfvx, rdx, cf1, cf2, cf3, ide, jde,         &
              kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dx            )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dx                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, zx, rdz
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfvx
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte
      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dx

        ! Local vars
      integer :: i, j, k
      real    :: dz_dx, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_tx_vx
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: dv_dz_dz_dx, v_mvx, v_mvx_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_tx_vx(i, j) = msftx(i, j) * msfvx(i, j)
        end do
      end do

        ! Scale v component of vector to account for map projection (at u points)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            v_mvx(i, k, j) = 0.25 * (v(i, k, j) + v(i, k, j + 1) + v(i - 1, k, j) + v(i - 1, k, j + 1)) / msfvx(i, j)
          end do
        end do
      end do

        ! Interpolate v_mvx from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            v_mvx_at_full(i, k, j) = fnm(k) * v_mvx(i, k, j) + fnp(k) * v_mvx(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          v_mvx_at_full(i, 1, j) = cf1 * v_mvx(i, 1, j) + cf2 * v_mvx(i, 2, j) + cf3 * v_mvx(i, 3, j)
            ! Top face
          v_mvx_at_full(i, kte, j) = cft1 * v_mvx(i, ktes1, j) + cft2 * v_mvx(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dv_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at center of faces
            dz_dx = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
              ! dv_dz_dz_dx = (partial v_mvx / partial z)|sigma (dv_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dv_dz_dz_dx(i, k, j) = (0.5 * (v_mvx(i, k, j) + v_mvx(i + 1, k, j)) - &
                0.5 * (v_mvx(i, k - 1, j) + v_mvx(i + 1, k - 1, j))) * rdz(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial v / partial x)|z = msf_tx_vx * (partial v_mvx_at_full / partial x)|sigma - dv_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dv_dx(i, k, j) = msf_tx_vx(i, j) * ((v_mvx_at_full(i + 1, k, j) - v_mvx_at_full(i, k, j)) * rdx -  &
                dv_dz_dz_dx(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dv_dx(i, 1, j) =  0.0
          dv_dx(i, kte, j) =  0.0
        end do
      end do

    end subroutine Calc_dv_dx


    subroutine Calc_dw_dx (w, zx, rdz, msftx, rdx, ims, ime, jms, jme, kms, kme,  &
              ide, jde, kde, its, ite, jts, jte, kts, kte, config_flags, dw_dx           )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dx                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, zx, rdz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx
      real,                                        intent (in)  :: rdx
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte
      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dx

        ! Local vars
      integer :: i, j, k, i_offset
      real    :: dz_dx
      real, dimension (its:ite, jts:jte) :: msf_tx_tx
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: dw_dz_dz_dx, w_mtx, w_mtx_at_half
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points squared
      do j = jts, jte
        do i = its, ite
          msf_tx_tx(i, j) = msftx(i, j) * msftx(i, j)
        end do
      end do

        ! Scale w component of vector to account for map projection
      if (USE_HIGH_ORDER) then
        i_offset = 2
      else
        i_offset = 1
      end if
      do j = j_start, j_end
        do k = kts, kte
          do i = i_start - i_offset, i_end + i_offset
            w_mtx(i, k, j) = w(i, k, j) / msftx(i, j)
          end do
        end do
      end do

        ! Average w_mtx to the mass points (half levels)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start - 1, i_end + 1
            w_mtx_at_half(i, k, j) = 0.5 * (w_mtx(i, k, j) + w_mtx(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dw_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at center of faces
            dz_dx = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
              ! dw_dz_dz_dx = (partial w_mtx / partial z)|sigma (dw_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dw_dz_dz_dx(i, k, j) = (w_mtx_at_half(i, k, j) - w_mtx_at_half(i, k - 1, j)) * rdz(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial w / partial x)|z = msf_tx_tx * (partial w_mtx / partial x)|sigma - dw_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            if (USE_HIGH_ORDER) then
              dw_dx(i, k, j) = msf_tx_tx(i, j) * (rdx * (8.0 * (w_mtx(i + 1, k, j) - w_mtx(i - 1, k, j)) - (w_mtx(i + 2, k, j) - w_mtx(i - 2, k, j)) ) / 12.0 - dw_dz_dz_dx(i, k, j))
            else
              dw_dx(i, k, j) = msf_tx_tx(i, j) * (rdx * 0.5 * (w_mtx(i + 1, k, j) - w_mtx(i - 1, k, j)) - dw_dz_dz_dx(i, k, j))
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dw_dx(i, 1, j) = 0.0
          dw_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dx


    subroutine Calc_du_dy (u, zy, rdz, dn, dnw, fnm, fnp, msfty, msfuy, rdy, cf1, cf2, cf3,       &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dy    )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dy                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, zy, rdz
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty, msfuy
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dy

        ! Local vars
      integer :: i, j, k
      real    :: dz_dy, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_ty_uy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: du_dz_dz_dy, u_muy, u_muy_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_ty_uy(i, j) = msfty(i, j) * msfuy(i, j)
        end do
      end do

        ! Scale u component of vector to account for map projection (at v points)
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            u_muy(i, k, j) = 0.25 * (u(i, k, j) + u(i + 1, k, j) + u(i, k, j - 1) + u(i + 1, k, j - 1)) / msfuy(i, j) 
          end do
        end do
      end do

        ! Interpolate u_muy from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            u_muy_at_full(i, k, j) = fnm(k) * u_muy(i, k, j) + fnp(k) * u_muy(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          u_muy_at_full(i, 1, j) =  cf1 * u_muy(i, 1, j) + cf2 * u_muy(i, 2, j) + cf3 * u_muy(i, 3, j)
            ! Top face
          u_muy_at_full(i, kte, j) = cft1 * u_muy(i, ktes1, j) +  cft2 * u_muy(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (du_dz_dz_dy)
      do j = j_start, j_end 
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at center of faces 
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
              ! du_dz_dz_dy = (partial u_muy / partial z)|sigma (du_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            du_dz_dz_dy(i, k, j) = (0.5 * (u_muy(i, k, j) + u_muy(i, k, j + 1)) - &
                0.5 * (u_muy(i, k - 1, j) + u_muy(i, k - 1, j + 1))) * rdz(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial u / partial y)|z = msf_ty_uy * (partial u_muy_at_full / partial y)|sigma - du_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            du_dy(i, k, j) = msf_ty_uy(i, j) * ((u_muy_at_full(i, k, j + 1) - u_muy_at_full(i, k, j)) * rdy -  &
                du_dz_dz_dy(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          du_dy(i, kts, j) = 0.0
            ! Top face
          du_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dy


    subroutine Calc_dv_dy (v, zy, rdz, dn, dnw, fnm, fnp, msfty, msfvy, rdy, cf1, cf2, cf3, ide, jds,     &
              jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dy   )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dy                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, zy, rdz
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty, msfvy
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jds, jde, kde,            &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte
      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dy

        ! Local vars
      integer :: i, j, k
      real    :: dz_dy, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_ty_vy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: dv_dz_dz_dy, v_mvy, v_mvy_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_ty_vy(i, j) = msfty(i, j) * msfvy(i, j)
        end do
      end do

        ! Scale v component of vector to account for map projection
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
              ! Because msfvy at the poles will be undefined (1./0.),
              ! we will have trouble. But we are OK since v at the
              ! poles is 0., and that takes precedence in this case
            if ((config_flags%polar) .and. ((j == jds) .or. (j == jde))) then
              v_mvy(i, k, j) = 0.0
            else
              v_mvy(i, k, j) = v(i, k, j) / msfvy(i, j)
            end if
          end do
        end do
      end do

        ! Interpolate v_mvy from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            v_mvy_at_full(i, k, j) = fnm(k) * v_mvy(i, k, j) + fnp(k) * v_mvy(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          v_mvy_at_full(i, 1, j) =  cf1 * v_mvy(i, 1, j) + cf2 * v_mvy(i, 2, j) + cf3 * v_mvy(i, 3, j)
            ! Top face
          v_mvy_at_full(i, kte, j) = cft1 * v_mvy(i, ktes1, j) + cft2 * v_mvy(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dv_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at center of faces
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
              ! dv_dz_dz_dy = (partial v_mvy / partial z)|sigma (dv_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dv_dz_dz_dy(i, k, j) = (0.5 * (v_mvy(i, k, j) + v_mvy(i, k, j + 1)) - &
                0.5 * (v_mvy(i, k - 1, j) + v_mvy(i, k - 1, j + 1))) * rdz(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial v / partial y)|z = msf_ty_vy * (partial v_mvy_at_full / partial y)|sigma - dv_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dv_dy(i, k, j) = msf_ty_vy(i, j) * ((v_mvy_at_full(i, k, j + 1) - v_mvy_at_full(i, k, j)) * rdy -  &
                dv_dz_dz_dy(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          dv_dy(i, 1, j) = 0.0
            ! Top face
          dv_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dv_dy


    subroutine Calc_dw_dy (w, zy, rdz, msfty, rdy, ims, ime, jms, jme, kms, kme,  &
              ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dy                  )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dy                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, zy, rdz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty
      real,                                        intent (in)  :: rdy
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dy

        ! Local vars
      integer :: i, j, k, j_offset
      real    :: dz_dy
      real, dimension (its:ite, jts:jte) :: msf_ty_ty
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: dw_dz_dz_dy, w_mty, w_mty_at_half
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1) 

        ! Map factor at mass points squared
      do j = jts, jte
        do i = its, ite
          msf_ty_ty(i, j) = msfty(i, j) * msfty(i, j)
        end do
      end do

        ! Scale w component of vector to account for map projection
      if (USE_HIGH_ORDER) then
        j_offset = 2
      else 
        j_offset = 1
      end if
      do j = j_start - j_offset, j_end + j_offset
        do k = kts, kte
          do i = i_start, i_end
            w_mty(i, k, j) = w(i, k, j) / msfty(i, j)
          end do
        end do
      end do

        ! Average w_mty to the mass points (half levels)
      do j = j_start - 1, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            w_mty_at_half(i,k,j) = 0.5 * (w_mty(i, k, j) + w_mty(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dw_dz_dz_dy)
      do j = j_start, j_end 
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at center of faces
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
              ! dw_dz_dz_dy = (partial w_mty / partial z)|sigma (dw_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dw_dz_dz_dy(i, k, j) = (w_mty_at_half(i, k, j) - w_mty_at_half(i, k - 1, j)) * rdz(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial w / partial y)|z = msf_ty_ty * (partial w_mty / partial y)|sigma - dw_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            if (USE_HIGH_ORDER) then
              dw_dy(i, k, j) = msf_ty_ty(i, j) * (rdy * (8.0 * (w_mty(i, k, j + 1) - w_mty(i, k, j - 1)) - (w_mty(i, k, j + 2) - w_mty(i, k, j - 2)) ) / 12.0 - dw_dz_dz_dy(i, k, j))
            else
              dw_dy(i, k, j) = msf_ty_ty(i, j) * (rdy * 0.5 * (w_mty(i, k, j + 1) - w_mty(i, k, j - 1)) - dw_dz_dz_dy(i, k, j))
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dw_dy(i, kts, j) = 0.0
          dw_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dy


    subroutine Calc_thetav (th, qv, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, thetav)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                   !
      ! Purpose: Calculate virtual potential temperature at mass points   !
      !                                                                   !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano                    !
      !                                                                   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: th, qv
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: thetav

        ! Local vars
      integer :: i, j, k
      real :: sqv
      integer :: ktf, i_start, i_end, j_start, j_end

      ktf = Min (kte, kde - 1)

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Calculate virtual potential temperature
      do j = j_start - 1, j_end + 1
        do k = kts, ktf
          do i = i_start - 1, i_end + 1
            sqv = qv(i, k, j) / (1.0 + qv(i, k, j))
            thetav(i, k, j) = th(i, k, j) * (1.0 + P608 * sqv)
          end do
        end do
      end do

    end subroutine Calc_thetav


    subroutine Calc_mass_var_derivatives (thetav, qv, msftx, msfty, msfux, msfuy, msfvx, msfvy,  &
              rdx, rdy, fnm, fnp, cf1, cf2, cf3, zx, zy, rdz, rdzw, dnw, dn, dthetav_dx,         &
              dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz, ids, ide, jds, jde, kds, kde,      &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte                         )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                     !
      ! Purpose: Calculate spatial derivatives of mass variables at the center of the faces !
      !                                                                                     !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano                                      !
      !                                                                                     !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (kms:kme),                   intent (in)  :: fnm, fnp, dn, dnw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfux, msfuy, msfvx, msfvy, msftx, msfty
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: thetav, qv, zx, zy
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: rdz, rdzw
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dthetav_dx, dthetav_dy, dthetav_dz, &
                                                                   dqv_dx,     dqv_dy,     dqv_dz

      logical, parameter :: DEBUG = .false.

        ! Calculate thetav derivatives
      call Calc_dmassvar_dx (thetav, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dx)

      call Calc_dmassvar_dy (thetav, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dy)

      call Calc_dmassvar_dz (thetav, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dz)

        ! Calculate qv derivatives
      call Calc_dmassvar_dx (qv, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dx)

      call Calc_dmassvar_dy (qv, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dy)

      call Calc_dmassvar_dz (qv, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dz)

      if (DEBUG) print *, 'Done Calc_mass_var_derivatives...'

    end subroutine Calc_mass_var_derivatives


    subroutine Calc_dmassvar_dx (var, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dx)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial x        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: var, zx, rdz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx
      real, dimension (kms:kme),                   intent (in)  :: fnm, fnp, dn, dnw
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dvar_dx

        ! Local vars
      integer :: i, j, k
      real :: dz_dx
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: var_at_full, dvar_dz_dz_dx
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktf = Min (kte, kde - 1)

      ktes1 = kte - 1
      ktes2 = kte - 2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Interpolate mass var from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            var_at_full(i, k, j) = 0.5 * (fnm(k) * (var(i - 1, k, j) + var(i, k, j)) + &
                fnp(k) * (var(i - 1, k - 1, j) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
          var_at_full(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
              cf3 * var(i, 3, j) + cf1 * var(i - 1, 1, j) +  &
              cf2 * var(i - 1, 2, j) + cf3 * var(i - 1, 3, j))
          var_at_full(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
              var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              var(i - 1, ktes1, j) + (var(i - 1, ktes1, j) - &
              var(i - 1, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dvar_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at center of faces
            dz_dx = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
              ! dvar_dz_dz_dx = (partial var / partial z)|sigma (dvar_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dvar_dz_dz_dx(i, k, j) =  (var(i, k, j) - var(i, k - 1, j)) * rdz(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial var / partial x)|z = msftx * (partial var_at_full / partial x)|sigma - dvar_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dvar_dx(i, k, j) = msftx(i, j) * (var_at_full(i + 1, k, j) -  var_at_full(i, k, j)) * rdx - dvar_dz_dz_dx(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dvar_dx(i, 1, j) = 0.0
          dvar_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dx


    subroutine Calc_dmassvar_dy (var, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dy)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial y        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: var, zy, rdz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty
      real, dimension (kms:kme),                   intent (in)  :: fnm, fnp, dn, dnw
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dvar_dy

        ! Local vars
      integer :: i, j, k
      real :: dz_dy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: var_at_full, dvar_dz_dz_dy
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktf = Min (kte, kde - 1)

      ktes1 = kte - 1
      ktes2 = kte - 2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Interpolate mass var from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1,   ktf
          do i = i_start, i_end
            var_at_full(i, k, j) = 0.5 * (fnm(k) * (var(i, k, j - 1) + var(i, k, j)) + &
                fnp(k) * (var(i, k - 1, j - 1) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
          var_at_full(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
              cf3 * var(i, 3, j) + cf1 * var(i, 1, j - 1) +  &
              cf2 * var(i, 2, j - 1) + cf3 * var(i, 3, j - 1))
          var_at_full(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
              var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              var(i, ktes1, j - 1) + (var(i, ktes1, j - 1) - &
              var(i, ktes2, j - 1)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dvar_dz_dz_dy)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at center of faces
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
              ! dvar_dz_dz_dy = (partial var / partial z)|sigma (dvar_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dvar_dz_dz_dy(i, k, j) = (var(i, k, j) - var(i, k - 1, j)) * rdz(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial var / partial y)|z = msfty * (partial var_at_full / partial y)|sigma - dvar_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dvar_dy(i, k, j) = msfty(i, j) * (var_at_full(i, k, j + 1) - var_at_full(i, k, j)) * rdy - dvar_dz_dz_dy(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dvar_dy(i, 1, j) = 0.0
          dvar_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dy


    subroutine Calc_dmassvar_dz (var, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial z        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: var, rdz
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dvar_dz

        ! Local vars
      integer :: i, j, k
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktf = Min (kte, kde - 1)

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Interpolate mass var from half levels to full levels (center of faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dvar_dz(i, k, j) = (var(i, k, j) - var(i, k - 1, j)) * rdz(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dvar_dz(i, kts, j) = 0.0
          dvar_dz(i, ktf + 1, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dz


    subroutine Fill_q_sq_with_q_sq_prog (q_sq_prog, q_sq, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq_prog
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      real,                                        intent (in)  :: cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq

        ! Local vars
      integer :: i, j, k
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end
      real    :: cft1, cft2

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            q_sq(i, k, j) = Max (Q_SQ_MIN, fnm(k) * q_sq_prog(i, k, j) + fnp(k) * q_sq_prog(i, k - 1, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          q_sq(i, kts, j) = Max (Q_SQ_MIN, q_sq(i, kts, j))
            ! Top face
          q_sq(i, ktf + 1, j) = Q_SQ_MIN
        end do
      end do

    end subroutine Fill_q_sq_with_q_sq_prog


    subroutine Fill_q_sq_prog_with_q_sq (q_sq_prog, q_sq, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_prog

        ! Local vars
      integer :: i, j, k
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end
      real    :: cft1, cft2

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_prog(i, k, j) = 0.5 * ( q_sq(i, k + 1, j) + q_sq(i, k, j) )
          end do
        end do
      end do

    end subroutine Fill_q_sq_prog_with_q_sq


    subroutine Calc_l_my_algebra_l2_or_l2p5 (my_level, pbl3d_l_opt, q_sq, z, dz, rdzw, rmol, wthv, dthetav_dz, l_master_l2, l_master_l2_in, l_master_l2p5, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                             !
        ! Purpose: Calculates the master length scale following MY74  !
        !                                                             !
        ! Author: Pedro A. Jimenez & Timothy W. Juliano               !
        !                                                             !
        ! Comments: Use an algebraic equation to calculate L          !
        !                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      character (len = *), intent(in) :: my_level
      integer,                                     intent (in), optional  :: pbl3d_l_opt
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in), optional  :: z, dz, rdzw
      real, dimension (ims:ime, jms:jme),          intent (in), optional  :: rmol
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in), optional  :: wthv
      real, dimension (its:ite, kts:kte, jts:jte), intent (in), optional  :: dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out), optional :: l_master_l2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in), optional :: l_master_l2_in
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out), optional :: l_master_l2p5
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte

        ! Local vars
      real, dimension (kts:kte) :: z_at_walls_1d
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real, parameter :: SMALL_VAL = 0.00001
      real, parameter :: alpha_1 = 0.23
        ! Original MYNN vals
      real, parameter :: alpha_2 = 1.0
      real, parameter :: alpha_3 = 5.0
      real, parameter :: alpha_4 = 100.0
      real, parameter :: cns = 2.7
        ! Operational MYNN vals
!      real, parameter :: alpha_2 = 0.65
!      real, parameter :: alpha_3 = 3.0
!      real, parameter :: alpha_4 = 20.0
!      real, parameter :: cns = 2.3
        ! Experimental MYNN vals
!      real, parameter :: alpha_2 = 0.3
!      real, parameter :: alpha_3 = 2.0
!      real, parameter :: alpha_4 = 10.0
!      real, parameter :: cns = 3.5
      real :: l0, l0_num, l0_den, qdz
      real :: q, q_dz, dz_wall
      real :: l_s, l_b, l_f, zet, kz, N, q_c

      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      level: if (my_level == 'l2') then   ! level 2 model

           ! Calculates l0 (Eq. 72 MY74)
        do j = j_start, j_end
          do i = i_start, i_end
            l0_num = SMALL_VAL
            l0_den = SMALL_VAL
            z_at_walls_1d(kts) = 0.0
            do k = kts + 1, ktf
              z_at_walls_1d(k) = z_at_walls_1d(k - 1) + 1.0 / rdzw(i, k - 1, j)
              dz_wall = 0.5 * (dz(i, k, j) + dz(i, k - 1, j))
              q = Sqrt (q_sq(i, k, j))
              q_dz = q * dz_wall
              l0_num = l0_num + q_dz * z_at_walls_1d(k)
              l0_den = l0_den + q_dz
            end do
            z_at_walls_1d(ktf + 1) = z_at_walls_1d(ktf) + 1.0 / rdzw(i, ktf, j)

            if (pbl3d_l_opt == 1) then  ! MY alpha
              l0 = alpha * l0_num / l0_den
            else if (pbl3d_l_opt == 2) then  ! MYNN alpha
              l0 = alpha_1 * l0_num / l0_den
              q_c = (G_OVER_TREF * Max (wthv(i, kts, j), 0.0) * l0) ** ( 1.0 / 3.0 )
            end if

            do k = kts, ktf
              if (pbl3d_l_opt == 1) then
                 ! Calculates master length scale (Eq. 71 MY74)
                l_master_l2(i, k, j) = l0 * KARMAN * z_at_walls_1d(k) / (KARMAN * z_at_walls_1d(k) + l0)
              else if (pbl3d_l_opt == 2) then
                 ! Calculates master length scale (Nakanishi 2001)
                 ! l surface
                kz  = KARMAN * z_at_walls_1d(k)
                zet = rmol(i, j) * z_at_walls_1d(k)
                if (zet > 0.0) then
                  l_s = kz / (1.0 + cns * MIN(1.0, zet))
                else
                  l_s = kz * (1.0 - alpha_4 * zet) ** 0.2
                end if
                 ! l buoyancy
                if (dthetav_dz(i, k, j) > 0.0) then
                  N = Sqrt (G_OVER_TREF * dthetav_dz(i, k, j))
                   ! note q_c is a velocity scale defined on the ground hence kts index
                  l_b = alpha_2 * Sqrt (q_sq(i, k, j)) / N * ( 1.0 + alpha_3 / alpha_2 * Sqrt (q_c / (l0 * N)) )
                  l_f = alpha_2 * Sqrt (q_sq(i, k, j)) / N
                else
                  l_b = 1.0E10
                  l_f = l_b
                end if

                l_master_l2(i, k, j) = MIN( (l_s * l0 * l_b) / (l_s * l0 + l_s * l_b + l0 * l_b), l_f )

              end if
            end do
            l_master_l2(i, kte, j) = l_master_l2(i, ktf, j)

            if (DEBUG) call Debug_l_my_algebra

          end do
        end do

      else if (my_level == 'l2p5') then   ! level 2.5 model

        do j = j_start, j_end
          do k = kts, ktf
            do i = i_start, i_end
              l_master_l2p5(i, k, j) = 0.5 * ( l_master_l2_in(i, k + 1, j) + l_master_l2_in(i, k, j) )
              if (l_master_l2p5(i, k, j) == 0.0) then
                print *, '--------------------------------------------------'
                print *, 'Bad l_master value at i, k, j = ', i, k, j, l_master_l2p5(i, k, j), l_master_l2_in(i, k + 1, j), l_master_l2_in(i, k, j)
              end if
            end do
          end do
        end do


           ! Calculates l0 (Eq. 72 MY74)
!        do j = j_start, j_end
!          do i = i_start, i_end
!            l0_num = SMALL_VAL
!            l0_den = SMALL_VAL
!            z_at_walls_1d(kts) = z(i, kts, j)
!            dz_wall = z(i, kts, j)
!            q = Sqrt (q_sq(i, kts, j))
!            q_dz = q * dz_wall
!            l0_num = l0_num + q_dz * z_at_walls_1d(kts)
!            l0_den = l0_den + q_dz
!            if (i==2 .and. j==1) then
!              print *, 'l2p5'
!              print *, q, q_dz, l0_num/l0_den
!            end if
!            do k = kts + 1, ktf
!              z_at_walls_1d(k) = z(i, k, j)
!              dz_wall = z(i, k, j) - z(i, k - 1, j)
!              q = Sqrt (q_sq(i, k, j))
!              q_dz = q * dz_wall
!              l0_num = l0_num + q_dz * z_at_walls_1d(k)
!              l0_den = l0_den + q_dz
!              if (i==2 .and. j==1) then
!                print *, 'l2p5'
!                print *, q, q_dz, l0_num/l0_den
!              end if
!            end do

!            l0 = alpha * l0_num / l0_den

             ! Calculates master length scale (Eq. 71 MY74)
!            do k = kts, ktf
!              l_master_l2p5(i, k, j) = l0 * KARMAN * z_at_walls_1d(k) / (KARMAN * z_at_walls_1d(k) + l0)
!            end do

!            if (DEBUG) then
!              print *, 'l_master_l2p5 = '
!              print *, l_master_l2p5(I_TO_TEST, :, J_TO_TEST)
!              print *, 'l0 =', l0
!              print *, 'l0_num =', l0_num
!              print *, 'l0_den =', l0_den
!              print *, 'z_at_walls = '
!              print *, z_at_walls_1d
!              print *, 'dz = '
!              print *, 0.5 * (z(i, 2:ktf, j) - z(i, 1:ktf-1, j))
!            end if

!          end do
!        end do

      end if level

      contains

        subroutine Debug_l_my_algebra

          if (i == I_TO_TEST .and. j == J_TO_TEST) then
              print *, 'l_master = '
              print *, l_master_l2(I_TO_TEST, :, J_TO_TEST)
              print *, 'l0 =', l0
              print *, 'l0_num =', l0_num
              print *, 'l0_den =', l0_den
              print *, 'z_at_walls = '
              print *, z_at_walls_1d
              print *, 'dz = '
              print *, dz(i, :, j)
              print *, 'rdzw = '
              print *, rdzw(i, :, j)
              print *, 'dzw = '
              print *, 1.0 / rdzw(i, :, j)
            end if

        end subroutine Debug_l_my_algebra

    end subroutine Calc_l_my_algebra_l2_or_l2p5


    subroutine Calc_fluxes_pbl_approx (config_flags, l_master, du_dz, dv_dz, dthetav_dz, dqv_dz, &
              t, qv, dn, dnw, fnm, fnp, dx, dy, dz, zx, zy, msftx, msfty, rdz, rdzw, cf1, cf2, cf3, dt, id, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts,    &
              jte, kts, kte, q_sq, q_sq_prog, turb_flux_u2, turb_flux_v2, turb_flux_w2,     &
              turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                 &
              turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_uqv,    &
              turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, &
              u, v, w, thetav, rho, ust, hfx, qfx, rmol)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                   !
      ! Purpose: Diagnose q_sq and turb fluxes                            !
      !                                                                   !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                     !
      !                                                                   !
      ! Comments: q_sq and fluxes are calculated at the faces             !
      !                                                                   !
      ! References: Y75                                                   !
      !                                                                   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)    :: config_flags

      real, dimension (kms:kme),                   intent (in)    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: dz, zx, zy, rdz, rdzw
      real,                                        intent (in)    :: cf1, cf2, cf3, dt, dx, dy
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz, dqv_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: l_master
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: turb_flux_u2, &
          turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw,        &
          turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v,                  &
          turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv,             &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, v, w, t, thetav, rho, qv
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfty, ust, hfx, qfx, rmol
      integer, intent (in)    :: id
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      real, dimension (its:ite, kts:kte, jts:jte) :: rif, sm_l2, sh_l2, sm, sh, q_ratio, q_sq_hl88
      logical, parameter :: use_hl88 = .true.
      character (4) :: my_level, mass_or_face

      mass_or_face = 'face'

        ! Calculate q_sq (level 2) or Fill q_sq with prognostic q_sq (level 2.5)
      if ( config_flags%pbl3d_prog .eq. 0 ) then   ! level 2 model
        my_level = 'l2'
          ! Calc stability functions for level 2 model
        call Calc_sm_sh_l2_or_l2p5 (mass_or_face=mass_or_face, my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
           ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
           its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm, sh=sh, rif=rif)
        call Calc_q_sq_l2_pbl_approx (mass_or_face=mass_or_face, l_master=l_master, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
            jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
            jts=jts, jte=jte, kts=kts, kte=kte, q_sq_face=q_sq, sm=sm, rif=rif)
          ! Relax q_sq at the boundaries (only for parent domain and specified B.C. runs)
        if ( id .eq. 1 .and. config_flags%specified) then
          call Taper_q_sq (q_sq, ids, jds, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        end if
      else if ( config_flags%pbl3d_prog .gt. 0 ) then   ! level 2.5 model
        if (use_hl88) then   ! calc sm and sh then q_sq for level 2 model to prepare for HL88 fix
          my_level = 'l2'
          call Calc_sm_sh_l2_or_l2p5 (mass_or_face=mass_or_face, my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm_l2, sh=sh_l2, rif=rif)
          call Calc_q_sq_l2_pbl_approx (mass_or_face=mass_or_face, l_master=l_master, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
              jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
              jts=jts, jte=jte, kts=kts, kte=kte, q_sq_hl88=q_sq_hl88, sm=sm_l2, rif=rif)
        end if
        call Fill_q_sq_with_q_sq_prog (q_sq_prog, q_sq, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        my_level = 'l2p5'
          ! Calc stability functions for level 2.5 model
        call Calc_sm_sh_l2_or_l2p5 (mass_or_face=mass_or_face, my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
          jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
          jts=jts, jte=jte, kts=kts, kte=kte, sm_l2=sm_l2, sh_l2=sh_l2, sm=sm, sh=sh, q_sq=q_sq, q_sq_hl88=q_sq_hl88, q_ratio=q_ratio, l_master_at_face=l_master)
      end if

        ! Diagnose master length scale at faces
      my_level = 'l2'
      call Calc_l_my_algebra_l2_or_l2p5 (my_level=my_level, pbl3d_l_opt=config_flags%pbl3d_l_opt, q_sq=q_sq, dz=dz, rdzw=rdzw, rmol=rmol, &
          wthv=turb_flux_wtheta_v, l_master_l2=l_master, dthetav_dz=dthetav_dz, &
          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

        ! Diagnose momentum, heat, moisture fluxes for pbl approx
      if ( config_flags%pbl3d_prog .eq. 0 ) then  ! For level 2 model, no need for realizability
        call Diagnose_fluxes_pbl_approx (fnm=fnm, fnp=fnp, u=u, v=v, t=t, thetav=thetav, rho=rho, qv=qv, ust=ust, hfx=hfx, qfx=qfx, &
            sm_tilde=sm, sh_tilde=sh, rif=rif, q_sq=q_sq, l_master=l_master, &
            du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, dqv_dz=dqv_dz, u2=turb_flux_u2, v2=turb_flux_v2, w2=turb_flux_w2, uv=turb_flux_uv, &
            uw=turb_flux_uw, vw=turb_flux_vw, uthv=turb_flux_utheta_v, vthv=turb_flux_vtheta_v, wthv=turb_flux_wtheta_v, th2v=turb_flux_theta2_v, &
            uqv=turb_flux_uqv, vqv=turb_flux_vqv, wqv=turb_flux_wqv, uth=turb_flux_utheta, vth=turb_flux_vtheta, wth=turb_flux_wtheta, ids=ids, ide=ide, &
            jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
      else if ( config_flags%pbl3d_prog .gt. 0 ) then  ! For level 2.5 model, include information for HL88 mod
        call Diagnose_fluxes_pbl_approx (fnm=fnm, fnp=fnp, u=u, v=v, t=t, thetav=thetav, rho=rho, qv=qv, ust=ust, hfx=hfx, qfx=qfx, &
            sm_tilde=sm, sh_tilde=sh, sm_l2=sm_l2, sh_l2=sh_l2, rif=rif, q_sq=q_sq, q_sq_hl88=q_sq_hl88, q_ratio=q_ratio, l_master=l_master, &
            du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, dqv_dz=dqv_dz, u2=turb_flux_u2, v2=turb_flux_v2, w2=turb_flux_w2, uv=turb_flux_uv, &
            uw=turb_flux_uw, vw=turb_flux_vw, uthv=turb_flux_utheta_v, vthv=turb_flux_vtheta_v, wthv=turb_flux_wtheta_v, th2v=turb_flux_theta2_v, &
            uqv=turb_flux_uqv, vqv=turb_flux_vqv, wqv=turb_flux_wqv, uth=turb_flux_utheta, vth=turb_flux_vtheta, wth=turb_flux_wtheta, ids=ids, ide=ide, &
            jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
      end if

      call Diagnose_fluxes_surface (config_flags%pbl3d_sfc_opt, u, v, t, thetav, rho, dz, qv, ust, hfx, qfx, &
          rmol, q_sq, l_master, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
          turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, &
          turb_flux_wtheta, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

    end subroutine Calc_fluxes_pbl_approx


    subroutine Calc_fluxes (config_flags, pbl_approx_opt, id, specified, q_sq, q_sq_prog, l_master, du_dx, du_dy, du_dz, dv_dx, dv_dy,  &
        dv_dz, dw_dx, dw_dy, dw_dz, dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, &
        dqv_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw,   &
        turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v,       &
        turb_flux_theta2_v, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,       &
        turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, mat_cond_heat, mat_cond_moist, t, thetav, qv, rho, hfx, qfx, rmol, u, v, ust, &
        dz, rdzw, fnm, fnp, cf1, cf2, cf3, dn, dnw, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
        its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                !
      ! Purpose: Diagnoses turb fluxes                                 !
      !                                                                !
      ! Author: Pedro A. Jimenez                                       !
      !                                                                !
      ! Comments: The fluxes are calculated at the faces.              !
      !                                                                !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)    :: config_flags

      logical,                                     intent (in)    :: specified
      integer,                                     intent (in)    :: pbl_approx_opt, id

      integer,                                     intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                                                     ims, ime, jms, jme, kms, kme, &
                                                                     its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dx, dv_dy, dw_dz, du_dy, dv_dx,      &
                                                                     dw_dx, du_dz, dw_dy, dv_dz, dthetav_dx, &
                                                                     dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, &
                                                                     dqv_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                                                                     turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                                                                     turb_flux_utheta_v, turb_flux_vtheta_v,   &
                                                                     turb_flux_wtheta_v, turb_flux_theta2_v,   &
                                                                     turb_flux_utheta, turb_flux_vtheta,       &
                                                                     turb_flux_wtheta, turb_flux_uqv,          &
                                                                     turb_flux_vqv, turb_flux_wqv,             &
                                                                     q_sq, l_master
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: mat_cond_heat, mat_cond_moist
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: u, v, rho, t, thetav, qv
      real, dimension (ims:ime, jms:jme),          intent (in)    :: ust, hfx, qfx, rmol
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: dz, rdzw
      real, dimension (kms:kme),                   intent (in)    :: fnm, fnp, dn, dnw
      real,                                        intent (in)    :: cf1, cf2, cf3

        ! Local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf, iter
      real :: qsq
      real, dimension (its:ite, kts:kte, jts:jte) :: sm, sh, sm_l2, sh_l2, rif, q_sq_hl88, q_ratio

      real (kind = SOLVER_PREC) :: tf_utv, tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv

      logical, parameter :: use_hl88 = .true.
      character (4) :: my_level, mass_or_face

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      mass_or_face = 'face'

        ! Calculate q_sq (level 2) or Fill q_sq with prognostic q_sq (level 2.5)
      if ( config_flags%pbl3d_prog .eq. 0 ) then   ! level 2 model
        my_level = 'l2'
        call Calc_sm_sh_l2_or_l2p5 (mass_or_face=mass_or_face, my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
           ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
           its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm, sh=sh, rif=rif)
        call Calc_q_sq_l2_pbl_approx (mass_or_face=mass_or_face, l_master=l_master, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
            jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
            jts=jts, jte=jte, kts=kts, kte=kte, q_sq_face=q_sq, sm=sm, rif=rif)
          ! Relax q_sq at the boundaries (only for parent domain and specified B.C. runs)
        if ( id .eq. 1 .and. config_flags%specified) then
          call Taper_q_sq (q_sq, ids, jds, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        end if
      else if ( config_flags%pbl3d_prog .gt. 0 ) then   ! level 2.5 model
        if (use_hl88) then   ! calc sm and sh then q_sq for level 2 model to prepare for HL88 fix
          my_level = 'l2'
          call Calc_sm_sh_l2_or_l2p5 (mass_or_face=mass_or_face, my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm_l2, sh=sh_l2, rif=rif)
          call Calc_q_sq_l2_pbl_approx (mass_or_face=mass_or_face, l_master=l_master, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
              jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
              jts=jts, jte=jte, kts=kts, kte=kte, q_sq_hl88=q_sq_hl88, sm=sm_l2, rif=rif)
        end if
        call Fill_q_sq_with_q_sq_prog (q_sq_prog, q_sq, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        my_level = 'l2p5'
          ! Calc stability functions for level 2.5 model
        call Calc_sm_sh_l2_or_l2p5 (mass_or_face=mass_or_face, my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
          jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
          jts=jts, jte=jte, kts=kts, kte=kte, sm_l2=sm_l2, sh_l2=sh_l2, sm=sm, sh=sh, q_sq=q_sq, q_sq_hl88=q_sq_hl88, q_ratio=q_ratio, l_master_at_face=l_master)
      end if

        ! Diagnose master length scale at faces
      my_level = 'l2'
      call Calc_l_my_algebra_l2_or_l2p5 (my_level=my_level, pbl3d_l_opt=config_flags%pbl3d_l_opt, q_sq=q_sq, dz=dz, rdzw=rdzw, rmol=rmol, &
          wthv=turb_flux_wtheta_v, l_master_l2=l_master, dthetav_dz=dthetav_dz, &
          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

      mat_cond_heat(its:ite,kts:kte,jts:jte) = -9999.
      mat_cond_moist(its:ite,kts:kte,jts:jte) = -9999.

      loop_j: do j = j_start, j_end
        loop_k: do k = kts + 1, ktf
          loop_i: do i = i_start, i_end

!            print *, 'i, k, j = ', i, k, j

            if_no_turb: if (q_sq(i, k, j) == Q_SQ_MIN) then

              call Set_fluxes_to_min (q_sq(i, k, j), turb_flux_u2(i, k, j),   &
                 turb_flux_v2(i, k, j), turb_flux_w2(i, k, j), turb_flux_uv(i, k, j), &
                 turb_flux_uw(i, k, j), turb_flux_vw(i, k, j), turb_flux_utheta_v(i, k, j), &
                 turb_flux_vtheta_v(i, k, j), turb_flux_wtheta_v(i, k, j), turb_flux_theta2_v(i, k, j), &
                 turb_flux_utheta(i, k, j), turb_flux_vtheta(i, k, j), turb_flux_wtheta(i, k, j),  &
                 turb_flux_uqv(i, k, j), turb_flux_vqv(i, k, j), turb_flux_wqv(i, k, j) )

            else

                ! Save temporary scalars
              if ( config_flags%pbl3d_prog .gt. 0 ) then  ! Level 2.5 model
                if ( q_ratio(i, k, j) <= 1.0 ) then  ! HL88 mod
                  qsq = q_sq_hl88(i, k, j)
                else
                  qsq = q_sq(i, k, j)
                end if
              else  ! Level 2 model
                qsq = q_sq(i, k, j)
              end if

              call Diagnose_fluxes (pbl_approx_opt, qsq, l_master(i, k, j), du_dx(i, k, j), du_dy(i, k, j),                          &
                  du_dz(i, k, j), dv_dx(i, k, j), dv_dy(i, k, j), dv_dz(i, k, j), dw_dx(i, k, j), dw_dy(i, k, j), dw_dz(i, k, j),    &
                  dthetav_dx(i, k, j), dthetav_dy(i, k, j), dthetav_dz(i, k, j), dqv_dx(i, k, j), dqv_dy(i, k, j), dqv_dz(i, k, j),  &
                  turb_flux_u2(i, k, j), turb_flux_v2(i, k, j), turb_flux_w2(i, k, j), turb_flux_uv(i, k, j), turb_flux_uw(i, k, j), &
                  turb_flux_vw(i, k, j), turb_flux_utheta_v(i, k, j), turb_flux_vtheta_v(i, k, j), turb_flux_wtheta_v(i, k, j),      &
                  turb_flux_theta2_v(i, k, j), turb_flux_uqv(i, k, j), turb_flux_vqv(i, k, j), turb_flux_wqv(i, k, j), tf_utv,       &
                  tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv, mat_cond_heat(i, k, j), mat_cond_moist(i, k, j)                            )

                ! Calc heat fluxes
              call Calc_heat_flux (t, qv, fnm, fnp, i, j, k, tf_utv, tf_vtv, tf_wtv,  tf_uqv, tf_vqv, tf_wqv, &
                  turb_flux_utheta(i, k, j), turb_flux_vtheta(i, k, j), turb_flux_wtheta(i, k, j), ims, ime, jms, jme, kms, kme)

                ! Apply HL88 mod for level 2.5 model
              if ( config_flags%pbl3d_prog .gt. 0 ) then
                if ( q_ratio(i, k, j) <= 1.0 ) then
                  turb_flux_u2(i, k, j) = turb_flux_u2(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_v2(i, k, j) = turb_flux_v2(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_w2(i, k, j) = turb_flux_w2(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_uv(i, k, j) = turb_flux_uv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_uw(i, k, j) = turb_flux_uw(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_vw(i, k, j) = turb_flux_vw(i, k, j) * q_ratio(i, k, j) ** 2.0

                  turb_flux_utheta_v(i, k, j) = turb_flux_utheta_v(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_vtheta_v(i, k, j) = turb_flux_vtheta_v(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_wtheta_v(i, k, j) = turb_flux_wtheta_v(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_utheta(i, k, j) = turb_flux_utheta(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_vtheta(i, k, j) = turb_flux_vtheta(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_wtheta(i, k, j) = turb_flux_wtheta(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_theta2_v(i, k, j) = turb_flux_theta2_v(i, k, j) * q_ratio(i, k, j) ** 2.0

                  turb_flux_uqv(i, k, j) = turb_flux_uqv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_vqv(i, k, j) = turb_flux_vqv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  turb_flux_wqv(i, k, j) = turb_flux_wqv(i, k, j) * q_ratio(i, k, j) ** 2.0
                end if
              end if

            end if if_no_turb

          end do loop_i
        end do loop_k
      end do loop_j

      call Diagnose_fluxes_surface (config_flags%pbl3d_sfc_opt, u, v, t, thetav, rho, dz, qv, ust, hfx, qfx, &
          rmol, q_sq, l_master, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
          turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, &
          turb_flux_wtheta, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

!      call Calc_sfc_fluxes (qv, hfx, qfx, rho, u, v, ust, turb_flux_uw, turb_flux_vw, turb_flux_wtheta, turb_flux_wqv, &
!          ims, ime, jms, jme, kms, kme, kts, i_start, i_end, j_start, j_end)
!      end do

    end subroutine Calc_fluxes


    subroutine Set_fluxes_to_min (q_sq, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, &
        tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v, tf_ut, tf_vt, tf_wt, tf_uqv, &
        tf_vqv, tf_wqv)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Set the turbulent fluxes to the minimum val  !
        !                                                       !
        ! Author: Pedro A. Jimenez                              !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, intent (inout) :: q_sq, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw,  &
          tf_utv, tf_vtv, tf_wtv, tf_t2v, tf_ut, tf_vt, tf_wt, tf_uqv, tf_vqv, &
          tf_wqv

      q_sq = Q_SQ_MIN
      tf_u2 = TURB_FLUX_MIN
      tf_v2 = TURB_FLUX_MIN
      tf_w2 = TURB_FLUX_MIN
      tf_uv = TURB_FLUX_MIN
      tf_uw = TURB_FLUX_MIN
      tf_vw = TURB_FLUX_MIN
      tf_utv = TURB_FLUX_MIN
      tf_vtv = TURB_FLUX_MIN
      tf_wtv = TURB_FLUX_MIN
      tf_t2v = TURB_FLUX_MIN
      tf_ut = TURB_FLUX_MIN
      tf_vt = TURB_FLUX_MIN
      tf_wt = TURB_FLUX_MIN
      tf_uqv = TURB_FLUX_MIN
      tf_vqv = TURB_FLUX_MIN
      tf_wqv = TURB_FLUX_MIN

    end subroutine Set_fluxes_to_min


    subroutine Diagnose_fluxes (pbl_approx_opt, q_sq_at_walls_1dpbl, l_master, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz,      &
        dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,        &
        turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_uqv, &
        turb_flux_vqv, turb_flux_wqv, tf_utv, tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv, mat_cond_heat, mat_cond_moist                )

      implicit none

      integer, intent(in)     :: pbl_approx_opt
      real, intent (in)    :: du_dx, dv_dy, dw_dz, du_dy, dv_dx,      &
                              dw_dx, du_dz, dw_dy, dv_dz,             &
                              dthetav_dx, dthetav_dy, dthetav_dz,     &
                              dqv_dx, dqv_dy, dqv_dz
      real, intent (in)    :: q_sq_at_walls_1dpbl, l_master
      real, intent (out)   :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                              turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                              turb_flux_utheta_v, turb_flux_vtheta_v,   &
                              turb_flux_wtheta_v, turb_flux_theta2_v,   &
                              turb_flux_uqv, turb_flux_vqv, turb_flux_wqv
      real, intent (out)   :: mat_cond_heat, mat_cond_moist
      real (kind = SOLVER_PREC), intent (out) :: tf_utv, tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv

        ! Local vars
      real (kind = SOLVER_PREC) :: q, q_sq, l, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_t2v
      real (kind = SOLVER_PREC) :: du_dx_tmp, dv_dy_tmp, dw_dz_tmp, du_dy_tmp, dv_dx_tmp,      &
                                   dw_dx_tmp, du_dz_tmp, dw_dy_tmp, dv_dz_tmp, dthetav_dx_tmp, &
                                   dthetav_dy_tmp, dthetav_dz_tmp, dqv_dx_tmp, dqv_dy_tmp, dqv_dz_tmp

      logical, parameter :: DEBUG = .false.

        ! Convert vars from single prec to double prec
      q_sq = Working_to_solver_prec(q_sq_at_walls_1dpbl)
      q = Sqrt ( q_sq )
      l = Working_to_solver_prec(l_master)
!      tf_u2 = Working_to_solver_prec(turb_flux_u2)
!      tf_v2 = Working_to_solver_prec(turb_flux_v2)
!      tf_w2 = Working_to_solver_prec(turb_flux_w2)
!      tf_uv = Working_to_solver_prec(turb_flux_uv)
!      tf_uw = Working_to_solver_prec(turb_flux_uw)
!      tf_vw = Working_to_solver_prec(turb_flux_vw)
!      tf_utv = Working_to_solver_prec(turb_flux_utheta_v)
!      tf_vtv = Working_to_solver_prec(turb_flux_vtheta_v)
!      tf_wtv = Working_to_solver_prec(turb_flux_wtheta_v)
!      tf_t2v = Working_to_solver_prec(turb_flux_theta2_v)

      !!! CHECK TO MAKE SURE WE CAN REPRODUCE PBL APPROXIMATION SOLUTION
      if (pbl_approx_opt == 1) then
        du_dx_tmp = 0.0
        du_dy_tmp = 0.0
        du_dz_tmp = Working_to_solver_prec(du_dz)
        dv_dx_tmp = 0.0
        dv_dy_tmp = 0.0
        dv_dz_tmp = Working_to_solver_prec(dv_dz)
        dw_dx_tmp = 0.0
        dw_dy_tmp = 0.0
        dw_dz_tmp = 0.0
        dthetav_dx_tmp = 0.0
        dthetav_dy_tmp = 0.0
        dthetav_dz_tmp = Working_to_solver_prec(dthetav_dz)
        dqv_dx_tmp = 0.0
        dqv_dy_tmp = 0.0
        dqv_dz_tmp = Working_to_solver_prec(dqv_dz)
      else
        du_dx_tmp = Working_to_solver_prec(du_dx)
        du_dy_tmp = Working_to_solver_prec(du_dy)
        du_dz_tmp = Working_to_solver_prec(du_dz)
        dv_dx_tmp = Working_to_solver_prec(dv_dx)
        dv_dy_tmp = Working_to_solver_prec(dv_dy)
        dv_dz_tmp = Working_to_solver_prec(dv_dz)
        dw_dx_tmp = Working_to_solver_prec(dw_dx)
        dw_dy_tmp = Working_to_solver_prec(dw_dy)
        dw_dz_tmp = Working_to_solver_prec(dw_dz)
        dthetav_dx_tmp = Working_to_solver_prec(dthetav_dx)
        dthetav_dy_tmp = Working_to_solver_prec(dthetav_dy)
        dthetav_dz_tmp = Working_to_solver_prec(dthetav_dz)
        dqv_dx_tmp = Working_to_solver_prec(dqv_dx)
        dqv_dy_tmp = Working_to_solver_prec(dqv_dy)
        dqv_dz_tmp = Working_to_solver_prec(dqv_dz)
      end if

        ! Solve system of equations for momentum fluxes
      call Solve_turb_system (q, q_sq, l, du_dx_tmp, du_dy_tmp, du_dz_tmp, dv_dx_tmp, dv_dy_tmp, dv_dz_tmp, dw_dx_tmp, dw_dy_tmp, dw_dz_tmp, &
          dthetav_dx_tmp, dthetav_dy_tmp, dthetav_dz_tmp, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v, mat_cond_heat)

        ! Solve system of equations for moisture fluxes
      call Solve_turb_system_moist (q, l, du_dx_tmp, du_dy_tmp, du_dz_tmp, dv_dx_tmp, dv_dy_tmp, dv_dz_tmp, dw_dx_tmp, dw_dy_tmp, dw_dz_tmp, dthetav_dx_tmp, &
          dthetav_dy_tmp, dthetav_dz_tmp, dqv_dx_tmp, dqv_dy_tmp, dqv_dz_tmp, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, &
          tf_uqv, tf_vqv, tf_wqv, mat_cond_moist)

        ! Convert vars from double prec to single prec
      turb_flux_u2 = Solver_to_working_prec(tf_u2)
      turb_flux_v2 = Solver_to_working_prec(tf_v2)
      turb_flux_w2 = Solver_to_working_prec(tf_w2)
      turb_flux_uv = Solver_to_working_prec(tf_uv)
      turb_flux_uw = Solver_to_working_prec(tf_uw)
      turb_flux_vw = Solver_to_working_prec(tf_vw)
      turb_flux_utheta_v = Solver_to_working_prec(tf_utv)
      turb_flux_vtheta_v = Solver_to_working_prec(tf_vtv)
      turb_flux_wtheta_v = Solver_to_working_prec(tf_wtv)
      turb_flux_theta2_v = Solver_to_working_prec(tf_t2v)
      turb_flux_uqv = Solver_to_working_prec(tf_uqv)
      turb_flux_vqv = Solver_to_working_prec(tf_vqv)
      turb_flux_wqv = Solver_to_working_prec(tf_wqv)

    end subroutine Diagnose_fluxes


    pure function Working_to_solver_prec (var)

      implicit none

      real (kind = WORKING_PREC), intent(in) :: var
      real (kind = SOLVER_PREC)              :: Working_to_solver_prec

      working_to_solver_prec = Real (var, kind = SOLVER_PREC)

    end function Working_to_solver_prec


    pure function Solver_to_working_prec (var)

      implicit none

      real (kind = SOLVER_PREC), intent(in) :: var
      real (kind = WORKING_PREC)            :: Solver_to_working_prec

      solver_to_working_prec = Real (var, kind = WORKING_PREC)

    end function Solver_to_working_prec


    subroutine Solve_turb_system (q, q_sq, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
              dthetav_dx, dthetav_dy, dthetav_dz, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v, mat_cond)

      implicit none

      real (kind = SOLVER_PREC), intent (in)    :: q, q_sq, l
      real (kind = SOLVER_PREC), intent (in)    :: du_dx, dv_dy, dw_dz, du_dy, dv_dx,      &
                                                   dw_dx, du_dz, dw_dy, dv_dz,             &
                                                   dthetav_dx, dthetav_dy, dthetav_dz
      real,                      intent (out)   :: mat_cond
      real (kind = SOLVER_PREC), intent (out)   :: tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v

        ! Local vars
      character (len = 1) :: equed
      integer, parameter   :: N_VARS = 10
      integer :: ipiv(N_VARS)
      integer :: info, info2, n_iter_sor
      real :: deter
      integer :: aa, niter, diter
      logical :: use_sor
      real (kind = SOLVER_PREC) :: rcond, cond, anorm, ferr, berr
      real (kind = SOLVER_PREC), dimension(4*N_VARS) :: work
      real (kind = SOLVER_PREC) :: swork(N_VARS*(N_VARS+N_VARS))
      real (kind = SOLVER_PREC), dimension(N_VARS) :: iwork, rsf, csf
      real (kind = SOLVER_PREC), allocatable, dimension(:) :: b(:), x(:)
      real (kind = SOLVER_PREC), allocatable, dimension(:,:) :: a(:, :), af(:, :)

      external dgesv
      external dsgesv
      external dgesvx
      real, external :: dlange
      external dgecon
      external dgetrf

      allocate (a(N_VARS, N_VARS), af(N_VARS, N_VARS), b(N_VARS), x(N_VARS))

        ! Fills in the B matrix
      call Fill_in_b_vector (b, q, q_sq, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)

        ! Fills in the A matrix
      call Fill_in_a_matrix (a, q, q_sq, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
          dthetav_dx, dthetav_dy, dthetav_dz)

!      if (SOLVE_WITH_SOR) then
!        use_sor =  Is_system_diag_dom (a, n_vars)
!        if (use_sor) then
!          call Solve_lin_sys_eqs_sor_gauss (a, b, n_vars, x, niter)
!        else
!          call Solve_lin_sys_eqs_gauss (a, b, n_vars, x)
!        end if
!      else
!        call Solve_lin_sys_eqs_gauss (a, b, n_vars, x)
!      end if

      call dgesvx ('N', 'N', N_VARS, 1, a, N_VARS, af, N_VARS, ipiv, equed, rsf, csf, b, N_VARS, x, N_VARS, rcond, ferr, berr, work, iwork, info)
      mat_cond = 1. / rcond

      tf_u2 = x(1)
      tf_v2 = x(2)
      tf_w2 = x(3)
      tf_uv = x(4)
      tf_uw = x(5)
      tf_vw = x(6)
      tf_utv = x(7)
      tf_vtv = x(8)
      tf_wtv = x(9)
      tf_t2v = x(10)

      deallocate (a, af, b, x)

    end subroutine Solve_turb_system


    function Is_system_diag_dom (coefs, n_vars) result(is_diag_dom)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                        !
      ! Purpose: Determines if a matrix is diagonally dominant !
      !                                                        !
      ! Author: Pedro A. Jimenez                               !
      !                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: n_vars
      real (kind = SOLVER_PREC), dimension(n_vars, n_vars), intent(in) :: coefs
      logical :: is_diag_dom

        ! Local vars
      integer :: i, j
      real (kind = SOLVER_PREC) :: sum_non_diag

      is_diag_dom = .true.

      do i = 1, n_vars
        sum_non_diag = 0.0
        do j = 1, n_vars
          if(i /= j) sum_non_diag = sum_non_diag + Abs (coefs(i, j))
        end do
        if(Abs (coefs(i, i)) < sum_non_diag) then
          is_diag_dom = .false.
          exit
        end if
      end do

    end function Is_system_diag_dom


    subroutine Solve_lin_sys_eqs_sor_gauss (a, b, n_vars, x, n_iter_sor)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                   !
      ! Purpose: Solves a linear system of eqs A * X = B  !
      !                                                   !
      ! Author: Pedro A. Jimenez                          !
      !                                                   !
      ! Method: Duccessive-over-relaxationa (SOR).        !
      !         and Gauss if it reaches the max iter.     !
      !         Values in the diagonal must be non zero.  !
      !                                                   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: n_vars
      real (kind = SOLVER_PREC), dimension(:,:), intent(inout) :: a(n_vars, n_vars)
      real (kind = SOLVER_PREC), dimension(:),   intent(inout) :: b(n_vars)
      real (kind = SOLVER_PREC), dimension(:),   intent(out)   :: x(n_vars)
      integer, intent(out)   :: n_iter_sor

      real (kind = SOLVER_PREC), parameter :: OMEGA = 1.2, EPS = 0.00001
      integer, parameter:: N_MAX_ITER = 100
      logical, parameter :: PRINT_INFO_ITER = .false.

      real (kind = SOLVER_PREC) :: conv, delta
      integer :: i, j, k


      x = 0.0
      do k = 1, N_MAX_ITER
        conv = 0.0
        do i = 1, n_vars
          delta = b(i)
          do j = 1, n_vars
            delta = delta - a(i, j) * x(j)
          end do
          x(i) = x(i) + OMEGA * delta / a(i, i)
          if (Abs (delta) > conv) conv = Abs (delta)
        end do
        if(conv < eps) exit
        if (PRINT_INFO_ITER) print *, 'Iteration.... ', k, x(1)
      end do

      n_iter_sor = k

      if (k > N_MAX_ITER) call Solve_lin_sys_eqs_gauss (a, b, n_vars, x)

    end subroutine Solve_lin_sys_eqs_sor_gauss


    subroutine Solve_lin_sys_eqs_gauss (a, b, N_VARS, x)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                   !
      ! Purpose: Solves a linear system of eqs A * X = B  !
      !                                                   !
      ! Author: Pedro A. Jimenez                          !
      !                                                   !
      ! Method: Gaussian elimination.                     !
      !         Values in the diagonal must be non zero.  !
      !                                                   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer,                          intent(in)    :: N_VARS
      real (kind = SOLVER_PREC), dimension(:,:), intent(inout) :: a(N_VARS, N_VARS)
      real (kind = SOLVER_PREC), dimension(:),   intent(inout) :: b(N_VARS)
      real (kind = SOLVER_PREC), dimension(:),   intent(out)   :: x(N_VARS)

        ! Local vars
      integer :: i, j, k
      real (kind = SOLVER_PREC) :: factor


      do k = 1, N_VARS - 1
        do i = k + 1, N_VARS
          factor = a(i, k) / a(k, k)
          a(i, k) = 0.0
          b(i) = b(i) - factor * b(k)
          do j = k + 1, N_VARS
            a(i, j) = a(i, j) - factor * a(k, j)
          end do
        end do
      end do

      x(N_VARS) = b(N_VARS) / a(N_VARS, N_VARS)
      do i = N_VARS - 1, 1, -1
        factor = 0.0
        do j = i + 1, N_VARS
          factor = factor + a(i, j) * x(j)
        end do
        x(i) = (b(i) - factor) / a(i, i)
      end do

    end subroutine Solve_lin_sys_eqs_gauss


    subroutine Solve_turb_system_moist (q, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dthetav_dx, &
              dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, &
              tf_uqv, tf_vqv, tf_wqv, mat_cond)

      implicit none

      real (kind = SOLVER_PREC), intent (in)    :: q, l
      real (kind = SOLVER_PREC), intent (in)    :: du_dx, dv_dy, dw_dz, du_dy, dv_dx,      &
                                                   dw_dx, du_dz, dw_dy, dv_dz,             &
                                                   dthetav_dx, dthetav_dy, dthetav_dz,     &
                                                   dqv_dx, dqv_dy, dqv_dz
      real (kind = SOLVER_PREC), intent (in)    :: tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv
      real (kind = SOLVER_PREC), intent (out)   :: tf_uqv, tf_vqv, tf_wqv
      real,                      intent (out)   :: mat_cond

        ! Local vars
      character (len = 1) :: equed
      integer, parameter   :: N_VARS = 4
      integer :: ipiv(N_VARS)
      integer :: info, info2, n_iter_sor
      real :: deter
      integer :: aa, niter, diter
      logical :: use_sor
      real (kind = SOLVER_PREC) :: rcond, cond, anorm, ferr, berr
      real (kind = SOLVER_PREC), dimension(4*N_VARS) :: work
      real (kind = SOLVER_PREC) :: swork(N_VARS*(N_VARS+N_VARS))
      real (kind = SOLVER_PREC), dimension(N_VARS) :: iwork, rsf, csf
      real (kind = SOLVER_PREC), allocatable, dimension(:) :: b(:), x(:)
      real (kind = SOLVER_PREC), allocatable, dimension(:,:) :: a(:, :), af(:, :)

      external dgesv
      external dsgesv
      external dgesvx
      external dgesvxx
      real, external :: dlange
      external dgecon
      external dgetrf

      allocate (a(N_VARS, N_VARS), af(N_VARS, N_VARS), b(N_VARS), x(N_VARS))

        ! Fills in the B matrix
      call Fill_in_b_vector_moist (b, dqv_dx, dqv_dy, dqv_dz, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv)

        ! Fills in the A matrix
      call Fill_in_a_matrix_moist (a, q, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
          dthetav_dx, dthetav_dy, dthetav_dz)

      call dgesvx ('N', 'N', N_VARS, 1, a, N_VARS, af, N_VARS, ipiv, equed, rsf, csf, b, N_VARS, x, N_VARS, rcond, ferr, berr, work, iwork, info)
      mat_cond = 1. / rcond

      tf_uqv = x(1)
      tf_vqv = x(2)
      tf_wqv = x(3)

      deallocate (a, af, b, x)

    end subroutine Solve_turb_system_moist


    subroutine Fill_in_b_vector (b, q, q_sq, l, u_x, u_y, u_z, v_x, &
        v_y, v_z, w_x, w_y, w_z)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Fills in the vector b with the indenpendent terms         !
      !                                                                    !
      ! Author: Pedro A. Jimenez                                           !
      !                                                                    !
      ! Comments: A * X = B wherein X are the tubulent fluxes of momentum, !
      !           theta_v, and theta2_v                                    !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      implicit none

      real (kind = SOLVER_PREC), intent(in) :: q, q_sq, l
      real (kind = SOLVER_PREC), intent(in) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
      real (kind = SOLVER_PREC), dimension(:), intent(inout) :: b(:)

      b = 0.0
      b(1) = Real (q ** 3.0 / (6.0 * a1 * l) + 3.0 * c1 * q_sq * u_x, kind = SOLVER_PREC)
      b(2) = Real (q ** 3.0 / (6.0 * a1 * l) + 3.0 * c1 * q_sq * v_y, kind = SOLVER_PREC)
      b(3) = Real (q ** 3.0 / (6.0 * a1 * l) + 3.0 * c1 * q_sq * w_z, kind = SOLVER_PREC)
      b(4) = Real (c1 * q_sq * (v_x + u_y), kind = SOLVER_PREC)
      b(5) = Real (c1 * q_sq * (w_x + u_z), kind = SOLVER_PREC)
      b(6) = Real (c1 * q_sq * (w_y + v_z), kind = SOLVER_PREC)

    end subroutine Fill_in_b_vector


    subroutine Fill_in_b_vector_moist (b, qv_x, qv_y, qv_z, tf_u2, tf_v2, tf_w2, tf_uv, &
        tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv)

        implicit none

      real (kind = SOLVER_PREC), intent(in) :: qv_x, qv_y, qv_z
      real (kind = SOLVER_PREC), intent(in) :: tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv
      real (kind = SOLVER_PREC), dimension(:), intent(inout) :: b(:)


      b(1) = Real (- (qv_x * tf_u2 + qv_y * tf_uv + qv_z * tf_uw), kind = SOLVER_PREC)
      b(2) = Real (- (qv_y * tf_v2 + qv_x * tf_uv + qv_z * tf_vw), kind = SOLVER_PREC)
      b(3) = Real (- (qv_z * tf_w2 + qv_x * tf_uw + qv_y * tf_vw), kind = SOLVER_PREC)
      b(4) = Real (- (qv_x * tf_utv + qv_y * tf_vtv + qv_z * tf_wtv), kind = SOLVER_PREC)

    end subroutine Fill_in_b_vector_moist


    subroutine Fill_in_a_matrix (a, q, q_sq, l, u_x, u_y, u_z, v_x, v_y, &
        v_z, w_x, w_y, w_z, tv_x, tv_y, tv_z)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Fills in the coefficients of the A Matrix                 !
      !                                                                    !
      ! Author: Pedro A. Jimenez                                           !
      !                                                                    !
      ! Comments: A * X = B wherein X are the tubulent fluxes of momentum, !
      !           theta_v, and theta2_v                                    !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!~!!!!!!

      implicit none

      real (kind = SOLVER_PREC), dimension(:, :), intent(inout) :: a(:, :)
      real (kind = SOLVER_PREC), intent(in) :: q, q_sq, l
      real (kind = SOLVER_PREC), intent(in) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, &
          tv_x, tv_y, tv_z

      a = 0.0

         ! u2 equation
      a(1, 1) = Real (q / (2.0 * a1 * l) + 2.0 * u_x, kind = SOLVER_PREC)
      a(1, 2) = Real (- v_y, kind = SOLVER_PREC)
      a(1, 3) = Real (- w_z, kind = SOLVER_PREC)
      a(1, 4) = Real (2.0 * u_y - v_x, kind = SOLVER_PREC)
      a(1, 5) = Real (2.0 * u_z - w_x, kind = SOLVER_PREC)
      a(1, 6) = Real (- v_z - w_y, kind = SOLVER_PREC)
      a(1, 9) = Real (G_OVER_TREF, kind = SOLVER_PREC)

         ! v2 equation
      a(2, 1) = Real (- u_x, kind = SOLVER_PREC)
      a(2, 2) = Real (q / (2.0 * a1 * l) + 2.0 * v_y, kind = SOLVER_PREC)
      a(2, 3) = Real (- w_z, kind = SOLVER_PREC)
      a(2, 4) = Real (2.0 * v_x - u_y, kind = SOLVER_PREC)
      a(2, 5) = Real (- u_z - w_x, kind = SOLVER_PREC)
      a(2, 6) = Real (2.0 * v_z - w_y, kind = SOLVER_PREC)
      a(2, 9) = Real (G_OVER_TREF, kind = SOLVER_PREC)

         ! w2 equation
      a(3, 1) = Real (- u_x, kind = SOLVER_PREC)
      a(3, 2) = Real (- v_y, kind = SOLVER_PREC)
      a(3, 3) = Real (q / (2.0 * a1 * l) + 2.0 * w_z, kind = SOLVER_PREC)
      a(3, 4) = Real (- u_y - v_x, kind = SOLVER_PREC)
      a(3, 5) = Real (2.0 * w_x - u_z, kind = SOLVER_PREC)
      a(3, 6) = Real (2.0 * w_y - v_z, kind = SOLVER_PREC)
      a(3, 9) = Real (- 2.0 * G_OVER_TREF, kind = SOLVER_PREC)

         ! uv equation
      a(4, 1) = Real (v_x, kind = SOLVER_PREC)
      a(4, 2) = Real (u_y, kind = SOLVER_PREC)
      a(4, 4) = Real (q / (3.0 * a1 * l) + u_x + v_y, kind = SOLVER_PREC)
      a(4, 5) = Real (v_z, kind = SOLVER_PREC)
      a(4, 6) = Real (u_z, kind = SOLVER_PREC)

       ! uw equation
      a(5, 1) = Real (w_x, kind = SOLVER_PREC)
      a(5, 3) = Real (u_z, kind = SOLVER_PREC)
      a(5, 4) = Real (w_y, kind = SOLVER_PREC)
      a(5, 5) = Real (q / (3.0 * a1 * l) + w_z + u_x, kind = SOLVER_PREC)
      a(5, 6) = Real (u_y, kind = SOLVER_PREC)
      a(5, 7) = Real (- G_OVER_TREF, kind = SOLVER_PREC)

         ! vw equation
      a(6, 2) = Real (w_y, kind = SOLVER_PREC)
      a(6, 3) = Real (v_z, kind = SOLVER_PREC)
      a(6, 4) = Real (w_x, kind = SOLVER_PREC)
      a(6, 5) = Real (v_x, kind = SOLVER_PREC)
      a(6, 6) = Real (q / (3.0 * a1 * l) + w_z + v_y, kind = SOLVER_PREC)
      a(6, 8) = Real (- G_OVER_TREF, kind = SOLVER_PREC)

         ! ut_v equation
      a(7, 1) = Real (tv_x, kind = SOLVER_PREC)
      a(7, 4) = Real (tv_y, kind = SOLVER_PREC)
      a(7, 5) = Real (tv_z, kind = SOLVER_PREC)
      a(7, 7) = Real (q / (3.0 * a2 * l) + u_x, kind = SOLVER_PREC)
      a(7, 8) = Real (u_y, kind = SOLVER_PREC)
      a(7, 9) = Real (u_z, kind = SOLVER_PREC)

         ! vt_v equation
      a(8, 2) = Real (tv_y, kind = SOLVER_PREC)
      a(8, 4) = Real (tv_x, kind = SOLVER_PREC)
      a(8, 6) = Real (tv_z, kind = SOLVER_PREC)
      a(8, 7) = Real (v_x, kind = SOLVER_PREC)
      a(8, 8) = Real (q / (3.0 * a2 * l) + v_y, kind = SOLVER_PREC)
      a(8, 9) = Real (v_z, kind = SOLVER_PREC)

         ! wt_v equation
      a(9, 3) = Real (tv_z, kind = SOLVER_PREC)
      a(9, 5) = Real (tv_x, kind = SOLVER_PREC)
      a(9, 6) = Real (tv_y, kind = SOLVER_PREC)
      a(9, 7) = Real (w_x, kind = SOLVER_PREC)
      a(9, 8) = Real (w_y, kind = SOLVER_PREC)
      a(9, 9) = Real (q / (3.0 * a2 * l) + w_z, kind = SOLVER_PREC)
      a(9, 10) = Real (- G_OVER_TREF, kind = SOLVER_PREC)

         ! tt_v equation
      a(10, 7) = Real (tv_x, kind = SOLVER_PREC)
      a(10, 8) = Real (tv_y, kind = SOLVER_PREC)
      a(10, 9) = Real (tv_z, kind = SOLVER_PREC)
      a(10, 10) = Real (q / (b2 * l), kind = SOLVER_PREC)

    end subroutine Fill_in_a_matrix


    subroutine Fill_in_a_matrix_moist (a, q, l, u_x, u_y, u_z, v_x, v_y, v_z, &
        w_x, w_y, w_z, tv_x, tv_y, tv_z)

      implicit none

      real (kind = SOLVER_PREC), dimension(:, :), intent(inout) :: a(:, :)
      real (kind = SOLVER_PREC), intent (in) :: q, l
      real (kind = SOLVER_PREC), intent (in) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, &
          w_y, w_z, tv_x, tv_y, tv_z

      a(1, 1) = Real (q / (3.0 * a2 * l) + u_x, kind = SOLVER_PREC)
      a(1, 2) = Real (u_y, kind = SOLVER_PREC)
      a(1, 3) = Real (u_z, kind = SOLVER_PREC)
      a(1, 4) = 0.0

        ! vr equation
      a(2, 1) = Real (v_x, kind = SOLVER_PREC)
      a(2, 2) = Real (q / (3.0 * a2 * l) + v_y, kind = SOLVER_PREC)
      a(2, 3) = Real (v_z, kind = SOLVER_PREC)
      a(2, 4) = 0.0

        ! wr equation
      a(3, 1) = Real (w_x, kind = SOLVER_PREC)
      a(3, 2) = Real (w_y, kind = SOLVER_PREC)
      a(3, 3) = Real (q / (3.0 * a2 * l) + w_z, kind = SOLVER_PREC)
      a(3, 4) = Real (-G_OVER_TREF, kind = SOLVER_PREC)

        ! wr equation
      a(4, 1) = Real (tv_x, kind = SOLVER_PREC)
      a(4, 2) = Real (tv_y, kind = SOLVER_PREC)
      a(4, 3) = Real (tv_z, kind = SOLVER_PREC)
      a(4, 4) = Real (2.0 * q / (b2 * l), kind = SOLVER_PREC)

    end subroutine Fill_in_a_matrix_moist


    subroutine Calc_heat_flux (th2, qv, fnm, fnp, i, j, k, tf_utv, &
        tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv, tf_ut, tf_vt, tf_wt, &
        ims, ime, jms, jme, kms, kme)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                             !
      ! Purpose: Calc heat fluxes based on the moist fluses and     !
      !          virtual potential temperature fluxes               !
      !                                                             !
      ! Author: Pedro A. Jimenez                                    !
      !                                                             !
      ! Reference: YM75 (Eq. A13)                                   !
      !                                                             !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: i, k, j, ims, ime, jms, jme, kms, kme
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: th2, qv
      real, dimension (kms:kme), intent (in) :: fnm, fnp
      real (kind = SOLVER_PREC), intent(in) :: tf_utv, tf_vtv, tf_wtv, &
          tf_uqv, tf_vqv, tf_wqv
      real, intent(out) :: tf_ut, tf_vt, tf_wt

        ! Local vars
      real :: th_wall, qv_wall

      th_wall = T0 + fnm(k) * th2(i, k, j) + fnp(k) * th2(i, k - 1, j)
      qv_wall = fnm(k) * qv(i, k, j) + fnp(k) * qv(i, k - 1, j)

      tf_ut = Real ((tf_utv - P608 * th_wall * tf_uqv) / (1.0 + P608 * qv_wall), &
          kind = WORKING_PREC)

      tf_vt = Real ((tf_vtv - P608 * th_wall * tf_vqv) / (1.0 + P608 * qv_wall), &
          kind = WORKING_PREC)

      tf_wt = Real ((tf_wtv - P608 * th_wall * tf_wqv) / (1.0 + P608 * qv_wall), &
          kind = WORKING_PREC)

    end subroutine Calc_heat_flux


    subroutine Diagnose_fluxes_pbl_approx (fnm, fnp, u, v, t, thetav, rho, qv, ust, hfx, qfx, &
              sm_tilde, sh_tilde, sm_l2, sh_l2, rif, q_sq, q_sq_hl88, q_ratio, l_master, du_dz, dv_dz, dthetav_dz, dqv_dz, &
              u2, v2, w2, uv, uw, vw, uthv, vthv, wthv, th2v, uqv, vqv, wqv, uth, vth, wth, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (kms:kme),                   intent (in)    :: fnm, fnp
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz, dqv_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: t, qv
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: u2, v2, w2, &
          uv, uw, vw, uthv, vthv, wthv, th2v, uqv, vqv, wqv, uth, vth, wth
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: sm_tilde, sh_tilde, rif
      real, dimension (its:ite, kts:kte, jts:jte), optional, intent (in)  :: q_ratio, q_sq_hl88, sm_l2, sh_l2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, v, thetav, rho
      real, dimension (ims:ime, jms:jme),          intent (in)  :: ust, hfx, qfx
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real :: qsq, l, sm, sh, ri, dudz, dvdz, dthvdz, dqvdz
      real :: t_wall, qv_wall
      real :: dxm, dym, tmpzx, tmpzy, sf_alpha
      real :: c_p, u1, v1, wind_speed, ust_sq
      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      j_loop: do j = j_start, j_end
        k_loop: do k = kts + 1, ktf
          i_loop: do i = i_start, i_end

            if_low_q_sq: if (q_sq(i, k, j) > Q_SQ_MIN) then

                ! Save temporary scalars
              if ( PRESENT(q_ratio) ) then  ! Level 2.5 model
                if ( q_ratio(i, k, j) <= 1.0 ) then  ! HL88 mod
                  qsq = q_sq_hl88(i, k, j)
                  sm = sm_l2(i, k, j)
                  sh = sh_l2(i, k, j)
                else
                  qsq = q_sq(i, k, j)
                  sm = sm_tilde(i, k, j)
                  sh = sh_tilde(i, k, j)
                end if
              else  ! Level 2 model
                qsq = q_sq(i, k, j)
                sm = sm_tilde(i, k, j)
                sh = sh_tilde(i, k, j)
              end if
              l = l_master(i, k, j)
              ri = rif(i, k, j)
              dudz = du_dz(i, k, j)
              dvdz = dv_dz(i, k, j)
              dthvdz = dthetav_dz(i, k, j)
              dqvdz = dqv_dz(i, k, j)

                ! Eq. 27 MY82
              uw(i, k, j) = Diag_uw_pbl_approx (qsq, l, sm, dudz)
              vw(i, k, j) = Diag_vw_pbl_approx (qsq, l, sm, dvdz)
              uv(i, k, j) = Diag_uv_pbl_approx (qsq, l, dudz, dvdz, uw(i, k, j), vw(i, k, j))

                ! Eq. 28 MY82
              uthv(i, k, j) = Diag_uthv_pbl_approx (l, sh, sm, dthvdz, dudz)
              vthv(i, k, j) = Diag_vthv_pbl_approx (l, sh, sm, dthvdz, dvdz)
              wthv(i, k, j) = Diag_wthv_pbl_approx (qsq, l, sh, dthvdz)
              wqv(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqvdz)

                ! Eq. 29 MY82
              th2v(i, k, j) = Diag_th2v_pbl_approx (l, sh, dthvdz)

                ! Eq. 26 MY82
              u2(i, k, j) = Diag_u2_pbl_approx (qsq, l, ri, dudz, dvdz, uw(i, k, j), vw(i, k, j), wthv(i, k, j))
              v2(i, k, j) = Diag_v2_pbl_approx (qsq, l, ri, dudz, dvdz, uw(i, k, j), vw(i, k, j), wthv(i, k, j))
              w2(i, k, j) = Diag_w2_pbl_approx (qsq, l, ri, dudz, dvdz, uw(i, k, j), vw(i, k, j), wthv(i, k, j))

              if (DEBUG) call Debug_negative_flux

              if (DEBUG) call Debug_big_qsq

                ! Diagnose the rest of qv fluxes
              uqv(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqvdz, uw(i, k, j), wqv(i, k, j))
              vqv(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqvdz, vw(i, k, j), wqv(i, k, j))

                ! Diagnose heat fluxes
              t_wall = T0 + fnm(k) * t(i, k, j) + fnp(k) * t(i, k - 1, j)
              qv_wall = fnm(k) * qv(i, k, j) + fnp(k) * qv(i, k - 1, j)

              uth(i, k, j) = Diag_xth_pbl_approx (t_wall, qv_wall, uthv(i, k, j), uqv(i, k, j))
              vth(i, k, j) = Diag_xth_pbl_approx (t_wall, qv_wall, vthv(i, k, j), vqv(i, k, j))
              wth(i, k, j) = Diag_xth_pbl_approx (t_wall, qv_wall, wthv(i, k, j), wqv(i, k, j))

                ! Apply HL88 mod for level 2.5 model
              if ( PRESENT(q_ratio) ) then
                if ( q_ratio(i, k, j) <= 1.0 ) then
                  u2(i, k, j) = u2(i, k, j) * q_ratio(i, k, j) ** 2.0
                  v2(i, k, j) = v2(i, k, j) * q_ratio(i, k, j) ** 2.0
                  w2(i, k, j) = w2(i, k, j) * q_ratio(i, k, j) ** 2.0
                  uv(i, k, j) = uv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  uw(i, k, j) = uw(i, k, j) * q_ratio(i, k, j) ** 2.0
                  vw(i, k, j) = vw(i, k, j) * q_ratio(i, k, j) ** 2.0

                  uthv(i, k, j) = uthv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  vthv(i, k, j) = vthv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  wthv(i, k, j) = wthv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  uth(i, k, j) = uth(i, k, j) * q_ratio(i, k, j) ** 2.0
                  vth(i, k, j) = vth(i, k, j) * q_ratio(i, k, j) ** 2.0
                  wth(i, k, j) = wth(i, k, j) * q_ratio(i, k, j) ** 2.0
                  th2v(i, k, j) = th2v(i, k, j) * q_ratio(i, k, j) ** 2.0

                  uqv(i, k, j) = uqv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  vqv(i, k, j) = vqv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  wqv(i, k, j) = wqv(i, k, j) * q_ratio(i, k, j) ** 2.0
                end if
              end if

            else

              call Set_turb_flux_min

            end if if_low_q_sq

          end do i_loop
        end do k_loop
      end do j_loop

        ! Zero fluxes at the suface
      call Set_turb_flux_min_surface

        ! Zero fluxes at the top
      call Set_turb_flux_min_top

      contains

        subroutine Debug_negative_flux

          if (u2(i, k, j) < 0 .or. v2(i, k, j) < 0 .or. w2(i, k, j) < 0) then
            print *, '------------------------------------------------------'
            print *, 'NEGATIVE FLUX'
            print *, 'i, k, j = ', i, k, j
            if (u2(i, k, j) < 0) then
              print *, 'u2, qsq, u2/qsq = ', u2(i, k, j), qsq, ABS( u2(i, k, j) / qsq )
            else if (v2(i, k, j) < 0) then
              print *, 'v2, qsq, v2/qsq = ', v2(i, k, j), qsq, ABS( v2(i, k, j) / qsq )
            else if (w2(i, k, j) < 0) then
              print *, 'w2, qsq, w2/qsq = ', w2(i, k, j), qsq, ABS( w2(i, k, j) / qsq )
            end if
          end if

        end subroutine Debug_negative_flux


        subroutine Debug_big_qsq

          if (u2(i, k, j) + v2(i, k, j) + w2(i, k, j) > (1.0 + epsilon(qsq)) * qsq) then
            print *, '------------------------------------------------------'
            print *, 'BIG QSQ'
            print *, 'i, k, j = ', i, k, j
            print *, 'u2, v2, w2, sum, = ', u2(i, k, j), v2(i, k, j), w2(i, k, j), (u2(i, k, j) + v2(i, k, j) + w2(i, k, j))
            print *, 'qsq, sum/qsq = ', qsq, (u2(i, k, j) + v2(i, k, j) + w2(i, k, j)) / qsq
            u2(i, k, j) = qsq / 3.0
            v2(i, k, j) = qsq / 3.0
            w2(i, k, j) = qsq / 3.0
          end if

        end subroutine Debug_big_qsq


        subroutine Set_turb_flux_min

          u2(i, k, j) = TURB_FLUX_MIN
          v2(i, k, j) = TURB_FLUX_MIN
          w2(i, k, j) = TURB_FLUX_MIN
          uv(i, k, j) = TURB_FLUX_MIN
          uw(i, k, j) = TURB_FLUX_MIN
          vw(i, k, j) = TURB_FLUX_MIN
          uqv(i, k, j) = TURB_FLUX_MIN
          vqv(i, k, j) = TURB_FLUX_MIN
          wqv(i, k, j) = TURB_FLUX_MIN
          uth(i, k, j) = TURB_FLUX_MIN
          uthv(i, k, j) = TURB_FLUX_MIN
          vth(i, k, j) = TURB_FLUX_MIN
          vthv(i, k, j) = TURB_FLUX_MIN
          wth(i, k, j) = TURB_FLUX_MIN
          wthv(i, k, j) = TURB_FLUX_MIN
          th2v(i, k, j) = TURB_FLUX_MIN

        end subroutine Set_turb_flux_min


        subroutine Set_turb_flux_min_surface

          do j = j_start, j_end
            do i = i_start, i_end
              u2(i, kts, j) = TURB_FLUX_MIN
              v2(i, kts, j) = TURB_FLUX_MIN
              w2(i, kts, j) = TURB_FLUX_MIN
              uv(i, kts, j) = TURB_FLUX_MIN
              uw(i, kts, j) = TURB_FLUX_MIN
              vw(i, kts, j) = TURB_FLUX_MIN
              uqv(i, kts, j) = TURB_FLUX_MIN
              vqv(i, kts, j) = TURB_FLUX_MIN
              wqv(i, kts, j) = TURB_FLUX_MIN
              uth(i, kts, j) = TURB_FLUX_MIN
              uthv(i, kts, j) = TURB_FLUX_MIN
              vth(i, kts, j) = TURB_FLUX_MIN
              vthv(i, kts, j) = TURB_FLUX_MIN
              wth(i, kts, j) = TURB_FLUX_MIN
              wthv(i, kts, j) = TURB_FLUX_MIN
              th2v(i, kts, j) = TURB_FLUX_MIN
            end do
          end do

        end subroutine Set_turb_flux_min_surface

        subroutine Set_turb_flux_min_top

          do j = j_start, j_end
            do i = i_start, i_end
              u2(i, ktf + 1, j) = TURB_FLUX_MIN
              v2(i, ktf + 1, j) = TURB_FLUX_MIN
              w2(i, ktf + 1, j) = TURB_FLUX_MIN
              uv(i, ktf + 1, j) = TURB_FLUX_MIN
              uw(i, ktf + 1, j) = TURB_FLUX_MIN
              vw(i, ktf + 1, j) = TURB_FLUX_MIN
              uqv(i, ktf + 1, j) = TURB_FLUX_MIN
              vqv(i, ktf + 1, j) = TURB_FLUX_MIN
              wqv(i, ktf + 1, j) = TURB_FLUX_MIN
              uth(i, ktf + 1, j) = TURB_FLUX_MIN
              uthv(i, ktf + 1, j) = TURB_FLUX_MIN
              vth(i, ktf + 1, j) = TURB_FLUX_MIN
              vthv(i, ktf + 1, j) = TURB_FLUX_MIN
              wth(i, ktf + 1, j) = TURB_FLUX_MIN
              wthv(i, ktf + 1, j) = TURB_FLUX_MIN
              th2v(i, ktf + 1, j) = TURB_FLUX_MIN
            end do
          end do

        end subroutine Set_turb_flux_min_top

    end subroutine Diagnose_fluxes_pbl_approx


    subroutine Diagnose_fluxes_surface (pbl3d_sfc_opt, u, v, t, thetav, rho, dz, qv, ust, hfx, qfx, rmol, &
              q_sq, l_master, u2, v2, w2, uv, uw, vw, uthv, vthv, wthv, th2v, uqv, vqv, wqv, uth, vth, wth, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                             !
        ! Purpose: Calculates the surface fluxes                      !
        !                                                             !
        ! Author: M. Eghdami & Pedro A. Jimenez & Timothy W. Juliano  !
        !                                                             !
        ! Comments: Nakanishi 2001 for 3D                             !
        !                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent (in)                                        :: pbl3d_sfc_opt
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: t, qv
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: u2, v2, w2, &
          uv, uw, vw, uthv, vthv, wthv, th2v, uqv, vqv, wqv, uth, vth, wth
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, v, thetav, rho, dz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: ust, hfx, qfx, rmol
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real :: zet, pmz1, pmz2, phz
      real :: qss, prt, gamma1, kz, elv
      real :: c_p, u1, v1, wind_speed
      real, dimension (its:ite, jts:jte) :: cosa, sina, ust_sq
      ! The parameters below depend on stability functions of module_sf_mynn
      REAL, PARAMETER :: cphm_st=5.0, cphm_unst=16.0, &
                         cphh_st=5.0, cphh_unst=16.0

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      do j = j_start, j_end
        do i = i_start, i_end
          c_p = CP * (1.0 + 0.84 * qv(i, kts, j))
          wth(i, kts, j) = hfx(i, j) / (rho(i, kts, j) * c_p)
          wthv(i, kts, j) = hfx(i, j) / (rho(i, kts, j) * CP)
          wqv(i, kts, j) = qfx(i, j) / rho(i, kts, j)
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          u1 = 0.5 * (u(i, 1, j) + u(i + 1, 1, j))
          v1 = 0.5 * (v(i, 1, j) + v(i, 1, j + 1))
          wind_speed = Max (Sqrt ((u1 ** 2.0) + (v1 ** 2.0)), 0.01)
          ust_sq(i, j) = ust(i, j) * ust(i, j)
          cosa(i, j) = u1 / wind_speed
          sina(i, j) = v1 / wind_speed
          uw(i, kts, j) = - ust_sq(i, j) * cosa(i, j)
          vw(i, kts, j) = - ust_sq(i, j) * sina(i, j)
        end do
      end do

        ! Calculate additional surface turbulent fluxes
      if (pbl3d_sfc_opt > 0) then
         ! Prandtl number
        gamma1 = 1.0 / 3.0 - 2.0 * (a1 / b1)
        prt = a1 * (gamma1 - c1) /  (a2 * gamma1)
        do j = j_start, j_end
          do i = i_start, i_end
            if (pbl3d_sfc_opt == 1) then
              zet  = 0.0
              pmz1 = 1.0
              pmz2 = 1.0
              phz  = prt
              elv  = 1.0 ! L/kz
            else if (pbl3d_sfc_opt == 2) then
              kz  = 0.5 * dz(i, kts, j) * KARMAN
              zet = 0.5 * dz(i, kts, j) * rmol(i, j)
              elv = 0.5 * l_master(i, kts + 1, j) / kz

               ! MYNN similarity functions (from Businger)
              if ( zet >= 0.0 ) then
                pmz1 = 1.0 + (cphm_st - 1.0) * zet ! This will always appear in the denomenator
                pmz2 = 1.0 + (cphm_st) * zet
                phz = prt + cphh_st * zet
              else
                pmz1 = 1.0 / (1.0 - cphm_unst * zet) ** 0.25 - zet
                pmz2 = 1.0 / (1.0 - cphm_unst * zet) ** 0.25
                phz = prt / sqrt(1.0 - cphh_unst * zet) ! cphh_unst does not match businger
              end if

!              pmz1 = psim(i, j) - zet ! This will always appear in the denomenator
!              pmz2 = psim(i, j)
!              phz = psih(i, j)

!              if (pmz1 .le. 0.0) then
!                write(90, *) '-----------------------------'
!                write(90, *) 'i =', i, 'j = ', j, 'k = ', k
!                write(90, *) 'kz = ', kz
!                write(90, *) 'zet = ', zet
!                write(90, *) 'elv = ', elv
!                write(90, *) 'pmz1 = ', pmz1
!                write(90, *) 'pmz1 - zet = ', pmz1 - zet
!                write(90, *) 'pmz2 = ', pmz2
!                write(90, *) 'phz = ', phz
!              end if

            end if

             ! Following equations derived from MY74
             ! Eq. 46
            qss = Max(Q_SQ_MIN, (b1 * elv * pmz1) ** (2.0 / 3.0) * ust_sq(i, j))
            q_sq(i, kts, j) = qss

             ! Other momentum Eq. 47a, b, c, d
            u2(i, kts, j) = qss * (1.0 / 3.0 + (a1 / b1) * &
                (4.0 * cosa(i, j) * cosa(i, j) * pmz2 - 2.0 * sina(i, j) * sina(i, j) * pmz2 + 2.0 * zet) / pmz1)
            v2(i, kts, j) = qss * (1.0 / 3.0 + (a1 / b1) * &
                (-2.0 * cosa(i, j) * cosa(i, j) * pmz2 + 4.0 * sina(i, j) * sina(i, j) * pmz2 + 2.0 * zet) / pmz1)
            w2(i, kts, j) = qss * (1.0 / 3.0 + (a1 / b1) * &
                (-2.0 * pmz2 - 4.0 * zet) / pmz1)
            uv(i, kts, j) = qss * (a1 / b1) * (6.0 * sina(i, j) * cosa(i, j) * pmz2) / pmz1

             ! Heat flux Eq. 48a, b
            uth(i, kts, j) = - wth(i, kts, j) * (3.0 * a2 / b1 ** (1.0 / 3.0)) * &
                (elv ** (2.0 / 3.0)) * (phz + pmz2) / pmz1 ** (1.0 / 3.0) * cosa(i, j)
            vth(i, kts, j) = - wth(i, kts, j) * (3.0 * a2 / b1 ** (1.0 / 3.0)) * &
                (elv ** (2.0 / 3.0)) * (phz + pmz2) / pmz1 ** (1.0 / 3.0) * sina(i, j)

             ! Moisture Eq. 48a, b
            uqv(i, kts, j) = - wqv(i, kts, j) * (3.0 * a2 / b1 ** (1.0 / 3.0)) * &
                (elv ** (2.0 / 3.0)) * (phz + pmz2) / pmz1 ** (1.0 / 3.0) * cosa(i, j)
            vqv(i, kts, j) = - wqv(i, kts, j) * (3.0 * a2 / b1 ** (1.0 / 3.0)) * &
                (elv ** (2.0 / 3.0)) * (phz + pmz2) / pmz1 ** (1.0 / 3.0) * sina(i, j)

             ! Heat variance Eq. 49
            th2v(i, kts, j) = wth(i, kts, j) ** 2.0 / MAX(ust_sq(i, j), 0.001) * (elv ** (2.0 / 3.0)) * &
                (b2 / b1 ** (1.0 / 3.0)) * phz / pmz1 ** (1.0 / 3.0)
          end do
        end do
      end if

    end subroutine Diagnose_fluxes_surface


    pure function Diag_uw_pbl_approx (q_sq, l_master, sm, du_dz)

      implicit none

      real, intent(in) :: q_sq, l_master, sm, du_dz
      real             :: diag_uw_pbl_approx

      diag_uw_pbl_approx = - Sqrt (q_sq) * l_master * sm * du_dz

    end function Diag_uw_pbl_approx


    pure function Diag_vw_pbl_approx (q_sq, l_master, sm, dv_dz)

      implicit none

      real, intent(in) :: q_sq, l_master, sm, dv_dz
      real             :: diag_vw_pbl_approx

      diag_vw_pbl_approx = - Sqrt (q_sq) * l_master * sm * dv_dz

    end function Diag_vw_pbl_approx


    pure function Diag_uv_pbl_approx (q_sq, l_master, du_dz, dv_dz, tf_uw, tf_vw)

      implicit none

      real, intent(in) :: q_sq, l_master, du_dz, dv_dz, tf_uw, tf_vw
      real             :: diag_uv_pbl_approx

      diag_uv_pbl_approx = 3.0 * a1 * l_master / Sqrt (q_sq) * (- tf_uw * dv_dz - tf_vw * du_dz)

    end function Diag_uv_pbl_approx


    pure function Diag_uthv_pbl_approx (l_master, sh, sm, dthetav_dz, du_dz)

      implicit none

      real, intent(in) :: l_master, sh, sm, dthetav_dz, du_dz
      real             :: diag_uthv_pbl_approx

      diag_uthv_pbl_approx = 3.0 * a2 * l_master ** 2.0 * (sh + sm) * dthetav_dz * du_dz

    end function Diag_uthv_pbl_approx


    pure function Diag_vthv_pbl_approx (l_master, sh, sm, dthetav_dz, dv_dz)

      implicit none

      real, intent(in) :: l_master, sh, sm, dthetav_dz, dv_dz
      real             :: diag_vthv_pbl_approx

      diag_vthv_pbl_approx = 3.0 * a2 * l_master ** 2.0 * (sh + sm) * dthetav_dz * dv_dz

    end function Diag_vthv_pbl_approx


    pure function Diag_wthv_pbl_approx (q_sq, l_master, sh, dthetav_dz)

      implicit none

      real, intent(in) :: q_sq, l_master, sh, dthetav_dz
      real             :: diag_wthv_pbl_approx

      diag_wthv_pbl_approx = - Sqrt (q_sq) * l_master * sh * dthetav_dz

    end function Diag_wthv_pbl_approx


    pure function Diag_wqv_pbl_approx (q_sq, l_master, sh, dqv_dz)

      implicit none

      real, intent(in) :: q_sq, l_master, sh, dqv_dz
      real             :: diag_wqv_pbl_approx

      diag_wqv_pbl_approx = - Sqrt (q_sq) * l_master * sh * dqv_dz

    end function Diag_wqv_pbl_approx


    pure function Diag_th2v_pbl_approx (l_master, sh, dthetav_dz)

      implicit none

      real, intent(in) :: l_master, sh, dthetav_dz
      real             :: diag_th2v_pbl_approx

      diag_th2v_pbl_approx = b2 * l_master ** 2.0 * sh * dthetav_dz ** 2.0

    end function Diag_th2v_pbl_approx


    pure function Diag_u2_pbl_approx (q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v)

      implicit none

      real, intent(in) :: q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v
      real             :: diag_u2_pbl_approx

      diag_u2_pbl_approx = q_sq / 3.0 + 2.0 * a1 * l_master / Sqrt (q_sq) *  &
          (- 2.0 * tf_uw * du_dz + tf_vw * dv_dz - G_OVER_TREF * tf_wtheta_v)

!      diag_u2_pbl_approx = q_sq / 3.0 + 2.0 * a1 * l_master / Sqrt (q_sq) *  &
!          (- (2.0 + rif) * tf_uw * du_dz + (1.0 - rif) * tf_vw * dv_dz)

    end function Diag_u2_pbl_approx


    pure function Diag_v2_pbl_approx (q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v)

      implicit none

      real, intent(in) :: q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v
      real             :: diag_v2_pbl_approx

      diag_v2_pbl_approx = q_sq / 3.0 + 2.0 * a1 * l_master / Sqrt (q_sq) *  &
          (tf_uw * du_dz - 2.0 * tf_vw * dv_dz - G_OVER_TREF * tf_wtheta_v)

!      diag_v2_pbl_approx = q_sq / 3.0 + 2.0 * a1 * l_master / Sqrt (q_sq) *  &
!          ((1 - rif) * tf_uw * du_dz - (2.0 + rif ) * tf_vw * dv_dz)

    end function Diag_v2_pbl_approx


    pure function Diag_w2_pbl_approx (q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v)

      implicit none

      real, intent(in) :: q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v
      real             :: diag_w2_pbl_approx

      diag_w2_pbl_approx = q_sq / 3.0 + 2.0 * a1 * l_master / Sqrt (q_sq) *  &
          (tf_uw * du_dz + tf_vw * dv_dz + 2.0 * G_OVER_TREF * tf_wtheta_v)

!      diag_w2_pbl_approx = q_sq * (1. / 3. - 2. * a1 / b1 - 6. * a1 / b1 * &
!          (1. - 2. / 3. * c2) * rif / (1. - rif))

!      diag_w2_pbl_approx = q_sq / 3.0 + 2.0 * a1 * l_master / Sqrt (q_sq) *  &
!          ((1.0 - 2.0 * rif) * tf_uw * du_dz + (1.0 - 2.0 * rif) * tf_vw * dv_dz)

    end function Diag_w2_pbl_approx


    pure function Diag_uqv_pbl_approx (q_sq, l_master, du_dz, dqv_dz, tf_uw, tf_wqv)

      implicit none

      real, intent(in) :: q_sq, l_master, du_dz, dqv_dz, tf_uw, tf_wqv
      real             :: diag_uqv_pbl_approx

      diag_uqv_pbl_approx = - (3.0 * a2 * l_master / Sqrt (q_sq)) * (tf_uw * dqv_dz + tf_wqv * du_dz)

    end function Diag_uqv_pbl_approx


    pure function Diag_vqv_pbl_approx (q_sq, l_master, dv_dz, dqv_dz, tf_vw, tf_wqv)

      implicit none

      real, intent(in) :: q_sq, l_master, dv_dz, dqv_dz, tf_vw, tf_wqv
      real             :: diag_vqv_pbl_approx

      diag_vqv_pbl_approx = - (3.0 * a2 * l_master / Sqrt (q_sq)) * (tf_vw * dqv_dz + tf_wqv * dv_dz)

    end function Diag_vqv_pbl_approx


    pure function Diag_xth_pbl_approx (t_wall, qv_wall, tf_xtheta_v, tf_xqv)

      implicit none

      real, intent(in) :: t_wall, qv_wall, tf_xtheta_v, tf_xqv
      real             :: diag_xth_pbl_approx

      diag_xth_pbl_approx = (tf_xtheta_v - P608 * t_wall * tf_xqv) / (1.0 + P608 * qv_wall)

    end function Diag_xth_pbl_approx


    subroutine Calc_kh (du_dz, dv_dz, dthetav_dz, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, &
              turb_flux_vtheta_v, turb_flux_wtheta_v, l_master, q_sq, khx, khy, khz, ids, ide, &
              jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                   !
      ! Purpose: Diagnose 3D eddy diffusivities for tracer mixing         !
      !                                                                   !
      ! Author: Timothy W. Juliano                                        !
      !                                                                   !
      ! Comments: kh* are calculated at the faces                         !
      !                                                                   !
      ! References: MY82                                                  !
      !                                                                   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: du_dz, dv_dz, dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: turb_flux_uw, turb_flux_vw, &
                                                                   turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master, q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: khx, khy, khz
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      do j = j_start, j_end
        do i = i_start, i_end
          khx(i, kts, j) = 0.0
          khy(i, kts, j) = 0.0
          khz(i, kts, j) = 0.0
        end do
      end do

        ! MY82, Eqs. 28 and 31
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! Khx
            if (dthetav_dz(i, k, j) == 0.0 .or. du_dz(i, k, j) == 0.0 .or. q_sq(i, k, j) == Q_SQ_MIN) then
              khx(i, k, j) = 0.0
            else
              khx(i, k, j) = Max( ( (turb_flux_utheta_v(i, k, j) / (3.0 * a2 * l_master(i, k, j) / sqrt(q_sq(i, k, j)))) &
                            + (turb_flux_uw(i, k, j) * dthetav_dz(i, k, j)) ) / (dthetav_dz(i, k, j) * du_dz(i, k, j)), 0.0 )
            end if
              ! Khy
            if (dthetav_dz(i, k, j) == 0.0 .or. dv_dz(i, k, j) == 0.0 .or. q_sq(i, k, j) == Q_SQ_MIN) then
              khy(i, k, j) = 0.0
            else
              khy(i, k, j) = Max( ( (turb_flux_vtheta_v(i, k, j) / (3.0 * a2 * l_master(i, k, j) / sqrt(q_sq(i, k, j)))) &
                            + (turb_flux_vw(i, k, j) * dthetav_dz(i, k, j)) ) / (dthetav_dz(i, k, j) * dv_dz(i, k, j)), 0.0 )
            end if
              ! Khz
            if (dthetav_dz(i, k, j) == 0.0 .or. q_sq(i, k, j) == Q_SQ_MIN) then
              khz(i, k, j) = 0.0
            else
              khz(i, k, j) = Max( -turb_flux_wtheta_v(i, k, j) / dthetav_dz(i, k, j), 0.0 )
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          khx(i, ktf + 1, j) = 0.0
          khy(i, ktf + 1, j) = 0.0
          khz(i, ktf + 1, j) = 0.0
        end do
      end do

    end subroutine Calc_kh


    subroutine Vertical_turb_mix (dth, du, dv, dw, dqv, khz, tracer, tracer_tend,   &
              pbl3d_opt, ss, config_flags, turb_flux_wtheta, turb_flux_uw, turb_flux_vw, &
              turb_flux_w2, turb_flux_wqv, rdz, dnw, dn, rho, fnm, fnp, n_tracer,    &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme,            &
              its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                               !
      ! Purpose: Performs vertical turbulent mixing   !
      !                                               !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano !
      !                                               !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: pbl3d_opt
      type (grid_config_rec_type), intent(in) :: config_flags
      integer,                     intent(in) :: ss, n_tracer
      integer,                     intent(in) :: ids, ide, jds, jde, kds, kde, &
                                                 ims, ime, jms, jme, kms, kme, &
                                                 its, ite, jts, jte, kts, kte

      real, dimension (kms:kme),                   intent(in) :: dnw, dn, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdz, khz, rho, &
          turb_flux_wtheta, turb_flux_uw, turb_flux_vw, turb_flux_w2,  &
          turb_flux_wqv
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) :: dth, &
          du, dv, dw, dqv
      real, dimension (ims:ime, kms:kme, jms:jme, n_tracer), intent(in) ::  tracer
      real, dimension (ims:ime, kms:kme, jms:jme, n_tracer), intent(inout) ::  tracer_tend
      
        ! Local vars
      integer :: im
      logical, parameter :: ADD_TEND_WIND=.true., ADD_TEND_T=.true., &
          ADD_TEND_Q=.true.
      
      if (ADD_TEND_WIND) then
        call Vertical_mixing_u (du, config_flags, turb_flux_uw, &
            dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms,   &
            jme, kms, kme, its, ite, jts, jte, kts, kte)

        call Vertical_mixing_v (dv, config_flags, turb_flux_vw, &
            dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms,   &
            jme, kms, kme, its, ite, jts, jte, kts, kte)

        if (pbl3d_opt > 0) then
          call Vertical_mixing_w (dw, config_flags, turb_flux_w2,   &
              dn, rho, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
        end if

      end if

      if (ADD_TEND_T) then
        call Vertical_mixing_s_with_flux (dth, config_flags, turb_flux_wtheta, &
            dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms,   &
            jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

      if (ADD_TEND_Q) then
        call Vertical_mixing_s_with_flux (dqv, config_flags, turb_flux_wqv,  &
            dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
            jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

      if (ss .eq. config_flags%pbl3d_nsteps .and. n_tracer .ge. param_first_scalar) then
        do im = param_first_scalar, n_tracer
          call Vertical_mixing_s (tracer_tend(ims, kms, jms, im), config_flags, tracer(ims, kms, jms, im), &
              khz, dnw, rdz, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        end do
      end if

    end subroutine Vertical_turb_mix


    subroutine Vertical_mixing_u (tendency, config_flags, turb_flux, dnw, rho, fnm, fnp,  & 
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                       !
      ! Purpose: Calculates tendencies of resolved variables  !
      !          at u points due to vertical turbulent        !
      !          mixing using turbulent fluxes.               !
      !                                                       !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano        !
      !                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer, intent(in ) :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: turb_flux, rho
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) :: tendency
      real, dimension (kms:kme), intent(in) ::    dnw, fnm, fnp

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real :: rho_avg
      real, dimension(its:ite, kts:kte, jts:jte) :: turb_flux_u
      real (kind = kind(turb_flux_u)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Define some indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
          config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
          config_flags%nested) i_end = Min (ide - 1, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
          config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
          config_flags%nested) j_end = Min (jde - 1, jte)

        ! Calc flux at u point in the faces
      do j = j_start, j_end
        do i = i_start, i_end
          rho_avg = 0.5 * (rho(i - 1, kts, j) + rho(i, kts, j)) 
          turb_flux_u(i, kts, j) = rho_avg * 0.5 * (turb_flux(i - 1, kts, j) + turb_flux(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf 
          do i = i_start, i_end
            rho_avg = 0.5 * (fnm(k) * (rho(i, k, j) + rho(i - 1, k, j)) + &
                fnp(k) * (rho(i, k - 1, j) + rho(i - 1, k - 1, j)))
            turb_flux_u(i, k, j) = rho_avg * 0.5 * (turb_flux(i - 1, k, j) + turb_flux(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          rho_avg = 0.5 * (rho(i - 1, ktf, j) + rho(i, ktf, j))
          turb_flux_u(i, ktf + 1, j) =  rho_avg * 0.5 * (turb_flux(i - 1, ktf + 1, j) + turb_flux(i, ktf + 1, j)) 
        end do
      end do

        ! Calc tendency
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tendency(i, k, j) = ( G / dnw(k) * (turb_flux_u(i, k + 1, j) - turb_flux_u(i, k, j)) ) / N_PBL3D
          end do
        end do
      end do

    end subroutine Vertical_mixing_u


    subroutine Vertical_mixing_v (tendency, config_flags, turb_flux, dnw, rho, fnm, fnp,  &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                       !
      ! Purpose: Calculates tendencies of resolved variables  !
      !          at v points due to vertical turbulent        !
      !          mixing using turbulent fluxes.               !
      !                                                       !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano        !
      !                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags
      integer, intent (in) :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte
      real, dimension (ims:ime, kms:kme, jms:jme ), intent(in) :: turb_flux, rho
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) ::tendency
      real, dimension (kms:kme), intent(in) :: dnw, fnm, fnp

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real :: rho_avg
      real, dimension(its:ite, kts:kte, jts:jte) :: turb_flux_v
      real (kind = kind(turb_flux_v)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Define indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      IF (config_flags%open_xs .or. config_flags%specified .or. &
          config_flags%nested) i_start = Max (ids + 1, its)
      IF (config_flags%open_xe .or. config_flags%specified .or. &
          config_flags%nested) i_end   = Min (ide - 1, ite)
      IF (config_flags%open_ys .or. config_flags%specified .or. &
          config_flags%nested) j_start = Max (jds + 1, jts)
      IF (config_flags%open_ye .or. config_flags%specified .or. &
          config_flags%nested) j_end   = Min (jde - 1, jte)

        ! Calc flux at v points in the faces
      do j = j_start, j_end
        do i = i_start, i_end
          rho_avg = 0.5 * (rho(i, kts, j - 1) + rho(i, kts, j))
          turb_flux_v(i, kts, j) = rho_avg * 0.5 * (turb_flux(i, kts, j - 1) + turb_flux(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            rho_avg = 0.5 * (fnm(k) * (rho(i, k, j - 1) + rho(i, k, j)) + &
                fnp(k) * (rho(i, k - 1, j - 1) + rho(i, k - 1, j)))
            turb_flux_v(i, k, j) = rho_avg * 0.5 * (turb_flux(i, k, j - 1) + turb_flux(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          rho_avg = 0.5 * (rho(i, ktf, j - 1) + rho(i, ktf, j))
          turb_flux_v(i, ktf + 1, j) =  rho_avg * 0.5 * (turb_flux(i, ktf + 1, j - 1) + turb_flux(i, ktf + 1, j))
        end do
      end do

        ! Calc tendency
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tendency(i, k, j) = ( G / dnw(k) * (turb_flux_v(i, k + 1, j) - turb_flux_v(i, k, j)) ) / N_PBL3D
          end do
        end do
      end do

    end subroutine Vertical_mixing_v


    subroutine Vertical_mixing_w (tendency, config_flags, turb_flux, dn, rho,  &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Calculates vertical velocity tendencies due  !
        !          to vertical turbulent mixing                 !
        !                                                       !
        ! Author: Pedro A. Jimenez                              !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer, intent(in) :: ids, ide, jds, jde, kds, kde, ims, &
          ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (kms:kme), intent(in) :: dn
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: turb_flux, rho
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) :: tendency

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real, dimension(its:ite, kts:kte, jts:jte) :: turb_flux_mass
      real (kind = kind(turb_flux_mass)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Define some dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end = Min (jde - 2, jte)

        ! Calc flux at mass points
      do j = j_start, j_end
        do k = kts + 1, ktf - 1
          do i = i_start, i_end
            turb_flux_mass(i, k, j) = 0.5 * (turb_flux(i, k + 1, j) + turb_flux(i, k, j)) * rho(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
          do i = i_start, i_end
            turb_flux_mass(i, kts, j) = 0.5 * (turb_flux(i, kts, j) + turb_flux(i, kts + 1, j)) * rho(i, kts, j)
            turb_flux_mass(i, ktf, j) = 0.5 * (turb_flux(i, ktf, j) + turb_flux(i, ktf + 1, j)) * rho(i, ktf, j)
          end do
      end do

        ! Calc tendency
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            tendency(i, k, j) = ( G * (turb_flux_mass(i, k, j) - turb_flux_mass(i, k - 1, j)) / dn(k) ) / N_PBL3D
          end do
        end do
      end do

    end subroutine Vertical_mixing_w


    subroutine Vertical_mixing_s (tendency, config_flags, var, khz, dnw, rdz, rho, fnm, fnp,  &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Calculates tendencies of resolved variables  !
        !          at mass points due to vertical turbulent     !
        !          mixing using turbulent fluxes.               !
        !          Code adapted from module_diffusion_em.F      !
        !                                                       !
        ! Authors: Timothy W. Juliano & Pedro A. Jimenez        !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer,                     intent(in) :: ids, ide, jds, jde, kds, kde, &
                                                 ims, ime, jms, jme, kms, kme, &
                                                 its, ite, jts, jte, kts, kte

      real, dimension (kms:kme),                   intent(in) :: dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: var, khz, rdz, rho
      real, dimension (ims:ime, kms:kme, jms:jme), intent(inout) :: tendency

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real, dimension(its:ite, kts:kte, jts:jte) :: H3, khz_avg

        ! Define some dimensions
      i_start = its
      i_end   = Min (ite, ide - 1)
      j_start = jts
      j_end   = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = Min (jde - 2, jte)

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            khz_avg(i, k, j) = khz(i, k, j) * (fnm(k) * rho(i, k, j) + fnp(k) * rho(i, k - 1, j))
            H3(i, k, j) = -khz_avg(i, k, j) * (var(i, k, j) - var(i, k - 1, j)) * rdz(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          H3(i, kts, j) = 0.
          H3(i, ktf + 1, j) = 0.
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tendency(i, k, j) = tendency(i, k, j) + ( G * (H3(i, k + 1, j) - H3(i, k, j)) / dnw(k) )
          end do
        end do
      end do

    end subroutine Vertical_mixing_s


    subroutine Vertical_mixing_s_with_flux (tendency, config_flags, turb_flux, dnw, rho, fnm, fnp,  &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Calculates tendencies of resolved variables  !
        !          at mass points due to vertical turbulent     !
        !          mixing using turbulent fluxes.               !
        !                                                       !
        ! Authors: Pedro A. Jimenez & Timothy W. Juliano        !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer,                     intent(in) :: ids, ide, jds, jde, kds, kde, &
                                                 ims, ime, jms, jme, kms, kme, &
                                                 its, ite, jts, jte, kts, kte

      real, dimension (kms:kme),                   intent(in) :: dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: turb_flux, rho
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) :: tendency

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real :: rho_kp1, rho_k
      real (kind = kind(turb_flux)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Define some dimensions
      i_start = its
      i_end   = Min (ite, ide - 1)
      j_start = jts
      j_end   = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = Min (jde - 2, jte)

        ! First layer
      do j = j_start, j_end
        do i = i_start, i_end
          rho_k = rho(i, kts, j)
          rho_kp1 = fnm(kts + 1) * rho(i, kts + 1, j) + fnp(kts + 1) * rho(i, kts, j)
          tendency(i, kts, j) = ( G * (rho_kp1 * turb_flux(i, kts + 1, j) - &
              rho_k * turb_flux(i, kts, j)) / dnw(kts) ) / N_PBL3D
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf - 1
          do i = i_start, i_end
            rho_k = fnm(k) * rho(i, k, j) + fnp(k) * rho(i, k - 1, j)
            rho_kp1 = fnm(k + 1) * rho(i, k + 1, j) + fnp(k + 1) * rho(i, k, j)
            tendency(i, k, j) = ( G * (rho_kp1 * turb_flux(i, k + 1, j) - &
                rho_k * turb_flux(i, k, j)) / dnw(k) ) / N_PBL3D
          end do
        end do
      end do

        ! Last layer
      do j = j_start, j_end
        do i = i_start, i_end
          rho_k = fnm(ktf) * rho(i, ktf, j) + fnp(ktf) * rho(i, ktf - 1, j)
          rho_kp1 = rho(i, ktf, j)
          tendency(i, ktf, j) = ( G * (rho_kp1 * turb_flux(i, ktf + 1, j) - &
              rho_k * turb_flux(i, ktf, j)) / dnw(ktf) ) / N_PBL3D
        end do
      end do

    end subroutine Vertical_mixing_s_with_flux

    subroutine Horizontal_turb_mix (rt_tendf, ru_tendf, rv_tendf, rw_tendf,    &
                                    moist_tendf, tracer_tend,                  &
                                    ss, config_flags,                          &
                                    turb_flux_u2, turb_flux_v2,                &
                                    turb_flux_uv, turb_flux_uw,                &
                                    turb_flux_vw, turb_flux_w2,                &
                                    turb_flux_utheta, turb_flux_vtheta,        &
                                    turb_flux_uqv, turb_flux_vqv,              &
                                    khx, khy, n_tracer, tracer,                &
                                    msfux, msfuy, msfvx, msfvy,                &
                                    msftx, msfty, dx, dy,                      &
                                    rdx, rdy, rdz, rdzw, fnm, fnp,             &
                                    cf1, cf2, cf3, zx, zy, dn, dnw, rho,       &
                                    ids, ide, jds, jde, kds, kde,              &
                                    ims, ime, jms, jme, kms, kme,              &
                                    its, ite, jts, jte, kts, kte               )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                               !
      ! Purpose: Performs horizontal turbulent mixing !
      !                                               !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano !
      !                                               !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

      TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
 
      INTEGER ,        INTENT(IN) :: ss 
      INTEGER ,        INTENT(IN) :: ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte

      REAL ,           INTENT(IN) :: cf1, cf2, cf3, dx, dy

      REAL , DIMENSION( kms:kme ) ,                 INTENT(IN) :: fnm
      REAL , DIMENSION( kms:kme ) ,                 INTENT(IN) :: fnp
      REAL , DIMENSION( kms:kme ) ,                 INTENT(IN) :: dnw
      REAL , DIMENSION( kms:kme ) ,                 INTENT(IN) :: dn

      REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN) ::  msfux, &
                                                                   msfuy, &
                                                                   msfvx, &
                                                                   msfvy, &
                                                                   msftx, &
                                                                   msfty

      REAL , DIMENSION( its:ite, kts:kte, jts:jte), INTENT(INOUT) :: rt_tendf, &
                                                                     ru_tendf, &
                                                                     rv_tendf, &
                                                                     rw_tendf, &
                                                                     moist_tendf

      REAL , DIMENSION( ims:ime, kms:kme, jms:jme, n_tracer), INTENT(INOUT) :: tracer_tend
      REAL , DIMENSION( ims:ime, kms:kme, jms:jme, n_tracer), INTENT(IN) :: tracer
      INTEGER                                               , INTENT(IN) :: n_tracer

      REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN) :: zx, zy, rdz, rdzw, &
                                                                  rho, khx, khy
                                                                  

      REAL, DIMENSION (ims:ime, kms:kme, jms:jme), intent(IN) :: turb_flux_u2,  &
          turb_flux_v2, turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_w2, &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_uqv, turb_flux_vqv

      REAL, INTENT(IN) :: rdx, rdy

        ! Local
      INTEGER :: im
      REAL, DIMENSION (its-1:ite, kts:kte, jts-1:jte) :: sf_alpha

      call Calc_slope_factor (config_flags, msftx, msfty, dx, dy, zx, zy, rdzw, sf_alpha, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_u (ru_tendf, config_flags, turb_flux_u2,  &
          turb_flux_uv, msfux, msfuy, rdx, rdy, fnm, fnp, dnw, zx, zy, &
          rdzw, rho, sf_alpha, ids, ide, jds, jde, kds, kde, ims, ime,    &
          jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_v (rv_tendf, config_flags, turb_flux_v2,  &
          turb_flux_uv, msfvx, msfvy, rdx, rdy, fnm, fnp, dnw, zx, zy, &
          rdzw, rho, sf_alpha, ids, ide, jds, jde, kds, kde, ims, ime,    &
          jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_w (rw_tendf, config_flags, turb_flux_uw, &
          turb_flux_vw, msftx, msfty, rdx, rdy, fnm, fnp, dn, zx, zy, &
          rdz, rho, sf_alpha, ids, ide, jds, jde, kds, kde, ims, ime,    &
          jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_s_with_flux (rt_tendf, config_flags,       &
          msftx, msfty, rdx, rdy, fnm, fnp, zx, zy, rdz, rdzw, dnw,     &
          rho, sf_alpha, turb_flux_utheta, turb_flux_vtheta,               &
          ids, ide, jds, jde, kds, kde, ims, ime,                       &
          jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_s_with_flux (moist_tendf, config_flags,            &
          msftx, msfty, rdx, rdy, fnm, fnp, zx, zy, rdz, rdzw, dnw, rho,        &
          sf_alpha, turb_flux_uqv, turb_flux_vqv, ids, ide, jds, jde, kds, kde,    &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      if (ss .eq. config_flags%pbl3d_nsteps .and. n_tracer .ge. param_first_scalar) then
        do im = param_first_scalar, n_tracer
          call Horizontal_mixing_s (tracer_tend(ims, kms, jms, im), config_flags, tracer(ims, kms, jms, im), &
              khx, khy, dn, dnw, rdzw, rho, fnm, fnp, cf1, cf2, cf3, rdx, rdy, zx, zy, msfux, msfvy, msftx, msfty, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        end do
      end if

    end subroutine Horizontal_turb_mix


    subroutine Calc_slope_factor (config_flags, msftx, msfty, dx, dy, zx, zy, rdzw, sf_alpha, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Calculate the slope factor to account for steep terrain   !
      !                                                                    !
      ! Author: Timothy W. Juliano                                         !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           Calculation is done at mass points                       !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent(in)  :: config_flags

      REAL, DIMENSION (ims:ime, jms:jme),          INTENT(IN)  :: msftx, msfty
      REAL,                                        INTENT(IN)  :: dx, dy
      REAL, DIMENSION (ims:ime, kms:kme, jms:jme), INTENT(IN)  :: zx, zy, rdzw
      REAL, DIMENSION (its-1:ite, kts:kte, jts-1:jte), INTENT(OUT) :: sf_alpha
      INTEGER,                                     INTENT(IN)  :: ide, jde, kde,                &
                                                                  ims, ime, jms, jme, kms, kme, &
                                                                  its, ite, jts, jte, kts, kte

        ! Local
      INTEGER :: i, j, k, i_start, i_end, j_start, j_end, ktf
      REAL    :: dxm, dym, tmpzx, tmpzy
      logical, parameter :: SLOPE_FACTOR = .true.

        ! Define dimensions
      i_start = its - 1
      i_end = ite
      j_start = jts - 1
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (SLOPE_FACTOR .and. config_flags%pbl3d_opt > 0) then
        do j = j_start, j_end
          do k = kts, ktf
            do i = i_start, i_end
              dxm = dx / msftx(i, j)
              dym = dy / msfty(i, j)
              tmpzx = (0.25 * ( abs(zx(i, k, j))+ abs(zx(i + 1, k, j)) + abs(zx(i, k + 1, j)) + abs(zx(i + 1, k + 1, j))) * rdzw(i, k, j) * dxm)
              tmpzy = (0.25 * ( abs(zy(i, k, j))+ abs(zy(i, k, j + 1)) + abs(zy(i, k + 1, j)) + abs(zy(i, k + 1, j + 1))) * rdzw(i, k, j) * dym)
              sf_alpha(i, k, j) = max(sqrt(tmpzx * tmpzx + tmpzy * tmpzy), 1.0)
            end do
          end do
        end do
      else
        do j = j_start, j_end
          do k = kts, ktf
            do i = i_start, i_end
              sf_alpha(i, k, j) = 1.0
            end do
          end do
        end do
      end if

    end subroutine Calc_slope_factor


    subroutine Horizontal_mixing_u (tendency, config_flags, turb_flux_u2,  &
        turb_flux_uv, msfux, msfuy, rdx, rdy, fnm, fnp, dnw, zx, zy, rdzw, &
        rho, sf_alpha, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme,    &
        kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Add the u tendency due to horizontal turb mix             !
      !                                                                    !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                      !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           The tendency is:                                         !
      !             du / dt = - [d (u'u') / dx + d (u'v') / dy]            !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, intent (in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: fnm, fnp, dnw
      real, dimension (ims:ime, jms:jme), intent(in) ::  msfux, msfuy 
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdzw, rho, &
          zx, zy, turb_flux_u2, turb_flux_uv
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha

        ! Local data
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      integer :: ie_ext, je_ext

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: &
          tf_u2_at_mass, tf_uv_at_corners

        ! Calc some indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
          config_flags%nested) i_start = Max(ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
          config_flags%nested) i_end = Min(ide - 1, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
          config_flags%nested) j_start = Max(jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
          config_flags%nested) j_end = Min(jde - 2, jte)

        ! tf_u2_at_mass = rho * turb_flux_u2 (move flux from center of faces to mass points)
      call Prep_turb_flux_u2 (config_flags, tf_u2_at_mass, turb_flux_u2, rho, ids, ide, &
          jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! tf_uv_at_corners = rho * turb_flux_uv (move flux from center of faces to half levels at corners)
      ie_ext = 0
      je_ext = 1
      call Prep_turb_flux_uv (config_flags, tf_uv_at_corners, turb_flux_uv, rho, ie_ext, je_ext, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! Calculate the u tendency due to horizontal turbulent mixing
      call Calc_htend_du (config_flags, tendency, turb_flux_u2, turb_flux_uv, tf_u2_at_mass, tf_uv_at_corners, &
          msfux, msfuy, rdx, rdy, dnw, zx, zy, rdzw, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

    end subroutine Horizontal_mixing_u


    subroutine Horizontal_mixing_v (tendency, config_flags, turb_flux_v2,  &
        turb_flux_uv, msfvx, msfvy, rdx, rdy, fnm, fnp, dnw, zx, zy, rdzw, &
        rho, sf_alpha, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme,    &
        kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Add the v tendency due to horizontal turb mix             !
      !                                                                    !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                      !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           The tendency is:                                         !
      !             dv / dt = - [d (u'v') / dx + d (v'v') / dy]            !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, intent (in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: fnm, fnp, dnw
      real, dimension (ims:ime, jms:jme), intent (in) :: msfvx, msfvy
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdzw, rho, &
          zx, zy, turb_flux_v2, turb_flux_uv
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      integer :: ie_ext, je_ext

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) ::  &
          tf_uv_at_corners, tf_v2_at_mass

        ! Indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end = Min (jde - 1, jte)

        ! tf_uv_at_corners = rho * turb_flux_uv (move flux from center of faces to half levels at corners)
      ie_ext = 1
      je_ext = 0
      call Prep_turb_flux_uv (config_flags, tf_uv_at_corners, turb_flux_uv, rho, ie_ext, je_ext, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! tf_v2_at_mass = rho * turb_flux_v2 (move flux from center of faces to mass points)
      call Prep_turb_flux_v2 (config_flags, tf_v2_at_mass, turb_flux_v2, rho, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! Calculate the v tendency due to horizontal turbulent mixing
      call Calc_htend_dv (config_flags, tendency, turb_flux_v2, turb_flux_uv, tf_v2_at_mass, tf_uv_at_corners, &
          msfvx, msfvy, rdx, rdy, dnw, zx, zy, rdzw, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

    end subroutine Horizontal_mixing_v


    subroutine Horizontal_mixing_w (tendency, config_flags, turb_flux_uw, &
        turb_flux_vw, msftx, msfty, rdx, rdy, fnm, fnp, dn, zx, zy, rdz,  &
        rho, sf_alpha, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme,     &
        kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Add the w tendency due to horizontal turb mix             !
      !                                                                    !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                      !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           The tendency is:                                         !
      !             dw / dt = - [d (u'w') / dx + d (v'w') / dy]            !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, intent(in) :: rdx, rdy

      real, dimension (kms:kme), intent(in) :: fnm, fnp, dn
      real, dimension (ims:ime, jms:jme), intent(in) ::  msftx, msfty

      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdz, rho, &
          zx, zy, turb_flux_uw, turb_flux_vw
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha

        ! Local data
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end, n_offset

      real, dimension(its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: &
          tf_uw_at_full, tf_vw_at_full

        ! Calc indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
          config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
          config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
          config_flags%nested) j_start = Max (jds + 1,jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
          config_flags%nested) j_end = Min (jde - 2,jte)

      if (.not. USE_HIGH_ORDER) then
          ! tf_uw_at_full = rho * turb_flux_uw (move flux from center of faces to u points in faces)
        call Prep_turb_flux_uw (config_flags, tf_uw_at_full, turb_flux_uw, fnm, fnp, rho, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

          ! tf_vw_at_full = rho * turb_flux_vw (move flux from center of faces to v points in faces)
        call Prep_turb_flux_vw (config_flags, tf_vw_at_full, turb_flux_vw, fnm, fnp, rho, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      else
        n_offset = 2
        do j = j_start - n_offset, j_end + n_offset
          do k = kts + 1, ktf
            do i = i_start - n_offset, i_end + n_offset
                ! tf_uw_at_full = rho * turb_flux_uw at the center of the faces
              tf_uw_at_full(i, k, j) = turb_flux_uw(i, k, j) * (fnm(k) * rho(i, k, j) + &
                  fnp(k) * rho(i, k - 1, j))
                ! tf_vw_at_full = rho * turb_flux_vw at the center of the faces
              tf_vw_at_full(i, k, j) = turb_flux_vw(i, k, j) * (fnm(k) * rho(i, k, j) + &
                  fnp(k) * rho(i, k - 1, j))
            end do
          end do
        end do
      end if

        ! Calculate the w tendency due to horizontal turbulent mixing
      call Calc_htend_dw (config_flags, tendency, turb_flux_uw, turb_flux_vw, tf_uw_at_full, tf_vw_at_full, &
          msftx, msfty, rdx, rdy, dn, zx, zy, rdz, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

    end subroutine Horizontal_mixing_w


    subroutine Horizontal_mixing_s_with_flux (tendency, config_flags, msftx, msfty, rdx, rdy, fnm, fnp,   &
              zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_x, turb_flux_y, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Add the scalar tendency due to horizontal turb mix        !
      !                                                                    !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                      !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           The tendency is:                                         !
      !             ds / dt = - [d (u's') / dx + d (v's') / dy]            !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, intent (in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: fnm, fnp, dnw
      real, dimension (ims:ime, jms:jme), intent (in) ::  &
          msftx, msfty
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: rdz, &
          rdzw, zx, zy, turb_flux_x, turb_flux_y, rho
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: tendency

        ! Local vars
      integer :: ktf, i_start, i_end, j_start, j_end
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) ::    &
          tf_x_at_u, tf_y_at_v

        ! Define indexes
      i_start = its
      i_end   = Min (ite, ide - 1)
      j_start = jts
      j_end   = Min (jte, jde - 1)
      ktf = Min (kte, kde-1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
           config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
           config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
           config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
           config_flags%nested) j_end = Min (jde - 2, jte)

        ! Scalar turbulent flux in x direction (at u points)
      call Prep_turb_flux_sx (turb_flux_x, tf_x_at_u, rho, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

        ! Scalar turbulent flux in y direction (at v points)
      call Prep_turb_flux_sy (turb_flux_y, tf_y_at_v, rho, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

        ! Calculate scalar tendency due to horizontal turbulent mixing (at mass points)
      call Calc_htend_ds (config_flags, tendency, turb_flux_x, turb_flux_y, tf_x_at_u, tf_y_at_v, rho, sf_alpha, &
          fnm, fnp, zx, zy, dnw, rdzw, msfty, msftx, rdx, rdy, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, &
          kts, kte, i_start, i_end, j_start, j_end, ktf)

    end subroutine Horizontal_mixing_s_with_flux


    subroutine Horizontal_mixing_s (tendency, config_flags, var, khx, khy, dn, dnw, rdzw, rho, fnm, fnp,  &
              cf1, cf2, cf3, rdx, rdy, zx, zy, msfux, msfvy, msftx, msfty, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Calculates tendencies of resolved variables  !
        !          at mass points due to vertical turbulent     !
        !          mixing using turbulent fluxes.               !
        !          Code adapted from module_diffusion_em.F      !
        !                                                       !
        ! Authors: Timothy W. Juliano & Pedro A. Jimenez        !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer,                     intent(in) :: ids, ide, jds, jde, kds, kde, &
                                                 ims, ime, jms, jme, kms, kme, &
                                                 its, ite, jts, jte, kts, kte

      real,                                        intent(in) :: cf1, cf2, cf3, rdx, rdy
      real, dimension (kms:kme),                   intent(in) :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent(in) :: msfux, msfvy, msftx, msfty
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: var, khx, khy, zx, zy, rdzw, rho
      real, dimension (ims:ime, kms:kme, jms:jme), intent(inout) :: tendency

        ! Local vars
      integer :: i, j, k, ktf, ktes1, ktes2
      integer :: i_start, i_end, j_start, j_end
      real :: mrdx, mrdy, tmpzx, tmpzy, rdzu, rdzv
      real, dimension(its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: H1, H2, H1avg, H2avg, &
                                                                    khx_avg, khy_avg, zx_at_m, zy_at_m

        ! Define some dimensions
      i_start = its
      i_end   = Min (ite, ide - 1)
      j_start = jts
      j_end   = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)
      ktes1 = kte - 1
      ktes2 = kte - 2

      if (config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = Min (jde - 2, jte)

! H1 = partial var over partial x

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            khx_avg(i, k, j) = 0.25 * (khx(i - 1, k, j) + khx(i, k, j) + khx(i - 1, k + 1, j) + khx(i, k + 1, j)) * &
                                            0.5 * (rho(i - 1, k, j) + rho(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            H1avg(i, k, j) = 0.5 * (fnm(k) * (var(i - 1, k, j) + var(i, k, j)) + fnp(k) * (var(i - 1, k - 1, j) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
          H1avg(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
                                    cf3 * var(i, 3, j) + cf1 * var(i - 1, 1, j) + &
                                    cf2 * var(i - 1, 2, j) + cf3 * var(i - 1, 3, j))
          H1avg(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
                                        var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
                                        var(i - 1, ktes1, j) + (var(i - 1, ktes1, j) - &
                                        var(i - 1, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            tmpzx = 0.5 * (zx(i, k, j) + zx(i, k + 1, j))
            rdzu = 2. / (1. / rdzw(i, k, j) + 1. / rdzw(i - 1, k, j))
            H1(i, k, j) = -msfux(i, j) * khx_avg(i, k, j) * ( rdx * (var(i, k, j) - var(i - 1, k, j)) - tmpzx * &
                           (H1avg(i, k + 1, j) - H1avg(i, k, j)) * rdzu )
          end do
        end do
      end do

! H2 = partial var over partial y

      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            khy_avg(i, k, j) = 0.25 * (khy(i, k, j - 1) + khy(i, k, j) + khy(i, k + 1, j - 1) + khy(i, k + 1, j)) * &
                                            0.5 * (rho(i, k, j - 1) + rho(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            H2avg(i, k, j) = 0.5 * (fnm(k) * (var(i, k, j - 1) + var(i, k, j)) + fnp(k) * (var(i, k - 1, j - 1) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
          H2avg(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
                                    cf3 * var(i, 3, j) + cf1 * var(i, 1, j - 1) + &
                                    cf2 * var(i, 2, j - 1) + cf3 * var(i, 3, j - 1))
          H2avg(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
                                        var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
                                        var(i, ktes1, j - 1) + (var(i, ktes1, j - 1) - &
                                        var(i, ktes2, j - 1)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            tmpzy = 0.5 * (zy(i, k, j) + zy(i, k + 1, j))
            rdzv = 2. / (1. / rdzw(i, k, j) + 1. / rdzw(i, k, j - 1))
            H2(i, k, j) = -msfvy(i, j) * khy_avg(i, k, j) * ( rdy * (var(i, k, j) - var(i, k, j - 1)) - tmpzy * &
                           (H2avg(i, k + 1, j) - H2avg(i, k, j)) * rdzv )
          end do
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            H1avg(i, k, j) = 0.5 * (fnm(k) * (H1(i + 1, k, j) + H1(i, k, j)) + &
                                    fnp(k) * (H1(i + 1, k - 1, j) + H1(i, k - 1, j)))
            H2avg(i, k, j) = 0.5 * (fnm(k) * (H2(i, k, j + 1) + H2(i, k, j)) + &
                                    fnp(k) * (H2(i, k - 1, j + 1) + H2(i, k - 1, j)))
            zx_at_m(i, k, j) = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
            zy_at_m(i, k, j) = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          H1avg(i, kts, j) = 0.
          H1avg(i, ktf + 1, j) = 0.
          H2avg(i, kts, j) = 0.
          H2avg(i, ktf + 1, j) = 0.
          zx_at_m(i, kts, j) = 0.25 * (zx(i, kts, j) + zx(i + 1, kts, j) + zx(i, kts + 1, j) + zx(i + 1, kts + 1, j))
          zy_at_m(i, kts, j) = 0.25 * (zy(i, kts, j) + zy(i, kts, j + 1) + zy(i, kts + 1, j) + zy(i, kts + 1, j + 1))
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msftx(i, j) * rdx
            mrdy = msfty(i, j) * rdy

            tendency(i, k, j) = tendency(i, k, j) + g / (dnw(k) * rdzw(i, k, j)) * &
                 ( mrdx * (H1(i + 1, k, j) - H1(i, k, j)) + mrdy * (H2(i, k, j + 1) - H2(i, k, j)) - &
                   msftx(i, j) * zx_at_m(i, k, j) * (H1avg(i, k + 1, j) - H1avg(i, k, j)) * rdzw(i, k, j) - &
                   msfty(i, j) * zy_at_m(i, k, j) * (H2avg(i, k + 1, j) - H2avg(i, k, j)) * rdzw(i, k, j) )
          end do
        end do
      end do

    end subroutine Horizontal_mixing_s


    subroutine Prep_turb_flux_u2 (config_flags, tf_at_mass, turb_flux_u2, rho, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent (in) :: ids, ide, jds, jde, kds, kde, ims, &
         ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (inout) :: tf_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: &
          turb_flux_u2, rho

        ! Local vars
      integer :: i, j, k, ktf, i_start, i_end, j_start, j_end

        ! Define indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end = Min (ide - 1, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end = Min (jde - 1, jte)

      i_start = i_start - 1

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tf_at_mass(i, k, j) = rho(i, k, j) * 0.5 * (turb_flux_u2(i, k + 1, j) + turb_flux_u2(i, k, j))
          end do
        end do
      end do

    end subroutine Prep_turb_flux_u2


    subroutine Prep_turb_flux_v2 (config_flags, tf_at_mass, turb_flux_v2, rho, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent (in) :: ids, ide, jds, jde, kds, kde, ims, &
         ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (inout) :: tf_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: &
          turb_flux_v2, rho

        ! Local vars
      integer :: i, j, k, ktf, i_start, i_end, j_start, j_end

        ! Define indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end = Min (ide - 1, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end = Min (jde - 1, jte)

      j_start = j_start - 1

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tf_at_mass(i, k, j) = rho(i, k, j) * 0.5 * (turb_flux_v2(i, k + 1, j) + turb_flux_v2(i, k, j))
          end do
        end do
      end do

    end subroutine Prep_turb_flux_v2


    subroutine Prep_turb_flux_uv (config_flags, tf_at_corners, turb_flux_uv, rho, ie_ext, je_ext, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent (in) :: ids, ide, jds, jde, kds, kde,  &
                              ims, ime, jms, jme, kms, kme,  &
                              its, ite, jts, jte, kts, kte

      integer, intent (in) :: ie_ext, je_ext

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1), &
          intent (inout) :: tf_at_corners
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in)  &
          ::  turb_flux_uv, rho

        ! Local vars
      integer :: i, j, k, ktf, i_start, i_end, j_start, j_end
      real  :: turb_flux_avg

        ! Calc indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end = Min (ide - 1 - ie_ext, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end = Min (jde - 1 - je_ext, jte)

      i_end = i_end + ie_ext
      j_end = j_end + je_ext

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            turb_flux_avg = 0.125 * (turb_flux_uv(i - 1, k, j) + turb_flux_uv(i, k, j) +  &
                turb_flux_uv(i - 1, k, j - 1) + turb_flux_uv(i, k, j - 1) + turb_flux_uv(i - 1, k + 1, j) &
                + turb_flux_uv(i, k + 1, j) + turb_flux_uv(i - 1, k + 1, j - 1) + turb_flux_uv(i, k + 1, j - 1))
            tf_at_corners(i, k, j) = turb_flux_avg * 0.25 * (rho(i - 1, k, j) + rho(i, k, j) + &
                rho(i - 1, k, j - 1) + rho(i, k, j - 1))
          end do
        end do
      end do

    end subroutine Prep_turb_flux_uv


    subroutine Prep_turb_flux_uw (config_flags, tf_at_full, turb_flux_uw,  &
        fnm, fnp, rho, ids, ide, jds, jde, kds, kde, ims, ime, jms,   &
        jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent (in) :: ids, ide, jds, jde, kds, kde, ims, ime, &
          jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (kms:kme), intent (in) :: fnm, fnp
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2), & 
          intent (inout) :: tf_at_full
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in):: rho, &
          turb_flux_uw

        ! Local vars
      integer  :: i, j, k, ktf, i_start, i_end, j_start, j_end
      real :: turb_flux_avg

        ! Define indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1 )
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .OR. config_flags%specified .OR. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .OR. config_flags%specified .OR. &
         config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .OR. config_flags%specified .OR. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .OR. config_flags%specified .OR. &
         config_flags%nested) j_end = Min (jde - 2, jte)

      i_end = i_end + 1

       ! Bottom face
      do j = j_start, j_end
        do i = i_start, i_end
          turb_flux_avg = 0.5 * (turb_flux_uw(i, kts, j) +  &
              turb_flux_uw(i - 1, kts, j))
          tf_at_full(i, kts, j) = turb_flux_avg * 0.5 * &
              (rho(i - 1, kts, j) + rho(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            turb_flux_avg = 0.5 * (turb_flux_uw(i, k, j) +  &
                turb_flux_uw(i - 1, k, j)) 
            tf_at_full(i, k, j) = turb_flux_avg * 0.5 * (fnm(k) *  &
                (rho(i - 1, k, j) + rho(i, k, j) ) + fnp(k) * &
                (rho(i - 1, k - 1, j) + rho(i, k - 1, j)))
          end do
        end do
      end do

       ! Top face
      do j = j_start, j_end
        do i = i_start, i_end
          tf_at_full(i, ktf + 1, j) = 0.0
        end do
      end do

    end subroutine Prep_turb_flux_uw


    subroutine Prep_turb_flux_vw (config_flags, tf_at_full, turb_flux_vw, &
        fnm, fnp, rho, ids, ide, jds, jde, kds, kde, ims, ime, jms,  &
        jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent(in)  :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, ims, ime, &
         jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (kms:kme), intent(in) :: fnm, fnp
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2), &
          intent (inout) :: tf_at_full
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rho, &
          turb_flux_vw

        ! Local vars
      integer :: i, j, k, ktf, i_start, i_end, j_start, j_end
      real :: turb_flux_avg

        ! Define indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      IF (config_flags%open_xs .OR. config_flags%specified .OR. &
         config_flags%nested) i_start = Max (ids + 1, its)
      IF (config_flags%open_xe .OR. config_flags%specified .OR. &
         config_flags%nested) i_end = Min (ide - 2, ite)
      IF (config_flags%open_ys .OR. config_flags%specified .OR. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      IF (config_flags%open_ye .OR. config_flags%specified .OR. &
         config_flags%nested) j_end = Min (jde - 2, jte)

      j_end = j_end + 1

       ! Bottom face
      do j = j_start, j_end
        do i = i_start, i_end
          turb_flux_avg = 0.5 * (turb_flux_vw(i, kts, j) +  &
              turb_flux_vw(i, kts, j - 1))
          tf_at_full(i, kts, j) = turb_flux_avg * 0.5 * &
              (rho(i, kts, j) + rho(i, kts, j - 1))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            turb_flux_avg = 0.5 * (turb_flux_vw(i, k, j) + turb_flux_vw(i, k, j - 1)) 
            tf_at_full(i, k, j) = turb_flux_avg * 0.5 * (fnm(k) * (rho(i, k, j) + &
                rho(i, k, j - 1)) + fnp(k) * (rho(i, k - 1, j) + rho(i, k - 1, j - 1)))
          end do
        end do
      end do

        ! Top face
      do j = j_start, j_end
        do i = i_start, i_end
          tf_at_full(i, ktf + 1, j) = 0.0
        end do
      end do

    end subroutine Prep_turb_flux_vw


    subroutine Prep_turb_flux_sx (turb_flux_x, tf_x_at_u, rho, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: turb_flux_x, rho
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (inout) :: tf_x_at_u

      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      real :: turb_flux_average

        ! tf_x_at_u = turb flux in x direction (move flux from center of faces to u points)
      do j = j_start, j_end
        do k = kts, ktf - 1
          do i = i_start, i_end + 1
            turb_flux_average = 0.25 * (turb_flux_x(i - 1, k, j) + turb_flux_x(i, k, j) +  &
               turb_flux_x(i - 1, k + 1, j) + turb_flux_x(i, k + 1, j))
            tf_x_at_u(i, k, j) = turb_flux_average  * 0.5 * (rho(i - 1, k, j) + rho(i, k, j))
          end do
        end do
      end do

       ! Top face
      do j = j_start, j_end
          do i = i_start, i_end + 1
            tf_x_at_u(i, ktf, j) = 0.25 * (turb_flux_x(i - 1, ktf, j) + turb_flux_x(i, ktf, j)) * 0.5 * (rho(i - 1, ktf, j) + rho(i, ktf, j))
          end do
      end do

    end subroutine Prep_turb_flux_sx


    subroutine Prep_turb_flux_sy (turb_flux_y, tf_y_at_v, rho, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: turb_flux_y, rho
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (inout) :: tf_y_at_v

      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      real :: turb_flux_average

        ! tf_y_at_v = turb flux in y direction (move flux from center of faces to v points)
      do j = j_start, j_end + 1
        do k = kts, ktf - 1
          do i = i_start, i_end
            turb_flux_average = 0.25 * (turb_flux_y(i, k, j - 1) + turb_flux_y(i, k, j) + &
                turb_flux_y(i, k + 1, j - 1) + turb_flux_y(i, k + 1, j))
            tf_y_at_v(i, k, j) = turb_flux_average * 0.5 * (rho(i, k, j - 1) + rho(i, k, j))
          end do
        end do
      end do

       ! Top face
      do j = j_start, j_end + 1
        do i = i_start, i_end
          tf_y_at_v(i, ktf, j) = 0.25 * (turb_flux_y(i, ktf, j - 1) + turb_flux_y(i, ktf, j)) * 0.5 * (rho(i, ktf, j - 1) + rho(i, ktf, j))
        end do
      end do

    end subroutine Prep_turb_flux_sy


    subroutine Calc_htend_du (config_flags, tendency, turb_flux_u2, turb_flux_uv, tf_u2_at_mass, tf_uv_at_corners, &
              msfux, msfuy, rdx, rdy, dnw, zx, zy, rdzw, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf
      real,    intent(in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: dnw
      real, dimension (ims:ime, jms:jme), intent(in) ::  msfux, msfuy
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdzw, rho, &
          zx, zy, turb_flux_u2, turb_flux_uv
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1), intent(in) :: &
          tf_u2_at_mass, tf_uv_at_corners

        ! Local data
      integer :: i, j, k
      real    :: mrdx, mrdy, tmpdz, tmpalpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: &
          tf_u2_at_u, tf_uv_at_u, zy_at_u, zx_at_u
      real (kind = kind(tf_u2_at_u)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Calculate component related to slope in sigma surfaces
        ! tf_u2_at_u = rho * turb_flux_u2 at u points in the faces
        ! tf_uv_at_u = rho * turb_flux_uv at u points in the faces
      do j = j_start, j_end
        do i = i_start, i_end
          tf_u2_at_u(i, kts, j) = 0.5 * (turb_flux_u2(i - 1, kts, j) + &
              turb_flux_u2(i, kts, j)) * 0.5 * (rho(i - 1, kts, j) +   &
              rho(i, kts, j))
          tf_uv_at_u(i, kts, j) = 0.5 * (turb_flux_uv(i - 1, kts, j) + &
              turb_flux_uv(i, kts, j)) * 0.5 * (rho(i - 1, kts, j) +   &
              rho(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            tf_u2_at_u(i, k, j) = 0.5 * (turb_flux_u2(i - 1, k, j) + &
                turb_flux_u2(i, k, j)) * 0.25 * (rho(i - 1, k, j) +   &
                rho(i, k, j) + rho(i - 1, k - 1, j) + rho(i, k - 1, j))
            tf_uv_at_u(i, k, j) = 0.5 * (turb_flux_uv(i - 1, k, j) + &
                turb_flux_uv(i, k, j)) * 0.25 * (rho(i - 1, k, j) +   &
                rho(i, k, j) + rho(i - 1, k - 1, j) + rho(i, k - 1, j))
              ! Horizontal derivatives of z at u points
            zy_at_u(i, k, j) = 0.125 * (zy(i - 1, k, j) + zy(i, k, j) + &
                zy(i - 1, k, j + 1) + zy(i, k, j + 1) + zy(i - 1, k + 1, j) + &
                zy(i, k + 1, j) + zy(i - 1, k + 1, j + 1) + zy(i, k + 1, j + 1))
            zx_at_u(i, k, j) = 0.5 * (zx(i, k, j) + zx(i, k + 1 , j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          tf_u2_at_u(i, ktf + 1, j) = 0.0
          tf_uv_at_u(i, ktf + 1, j) = 0.0
          zx_at_u(i, kts, j) = 0.5 * (zx(i, kts, j) + zx(i, kts + 1 , j))
          zy_at_u(i, kts, j) = 0.125 * (zy(i - 1, kts, j) + zy(i, kts, j) + &
              zy(i - 1, kts, j + 1) + zy(i, kts, j + 1) + zy(i - 1, kts + 1, j) + &
              zy(i, kts + 1, j) + zy(i - 1, kts + 1, j + 1) + zy(i, kts + 1, j + 1))
        end do
      end do

        ! Calc tendency using coordinate transform (sigma to cartesian)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msfux(i, j) * rdx
            mrdy = msfuy(i, j) * rdy
            tmpdz = (1.0 / rdzw(i, k, j) + 1.0 / rdzw(i - 1, k, j)) / 2.0
            tmpalpha = 0.5 * (sf_alpha(i - 1, k, j) + sf_alpha(i, k, j))
            tendency(i, k, j) = tendency(i, k, j) + ( G * tmpdz / dnw(k) *             &
               (mrdx * (tf_u2_at_mass(i, k, j) - tf_u2_at_mass(i - 1, k, j)) +                      &
                mrdy * (tf_uv_at_corners(i, k, j + 1) - tf_uv_at_corners(i, k, j)) -                      &
                msfux(i, j) * zx_at_u(i, k, j) * (tf_u2_at_u(i, k + 1, j) - tf_u2_at_u(i, k, j)) / tmpdz - &
                msfuy(i, j) * zy_at_u(i, k, j) * (tf_uv_at_u(i, k + 1, j) - tf_uv_at_u(i, k, j)) / tmpdz) ) / (N_PBL3D * tmpalpha)
          end do
        end do
      end do

    end subroutine Calc_htend_du


    subroutine Calc_htend_dv (config_flags, tendency, turb_flux_v2, turb_flux_uv, tf_v2_at_mass, tf_uv_at_corners, &
              msfvx, msfvy, rdx, rdy, dnw, zx, zy, rdzw, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf
      real,    intent(in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: dnw
      real, dimension (ims:ime, jms:jme), intent(in) ::  msfvx, msfvy
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdzw, rho, &
          zx, zy, turb_flux_v2, turb_flux_uv
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1), intent(in) :: &
          tf_v2_at_mass, tf_uv_at_corners

        ! Local data
      integer :: i, j, k
      real    :: mrdx, mrdy, tmpdz, tmpalpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: &
          tf_v2_at_v, tf_uv_at_v, zy_at_v, zx_at_v
      real (kind = kind(tf_v2_at_v)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Calculate component related to slope in sigma surfaces
        ! tf_uv_at_v = rho * turb_flux_uv at v points in the faces
        ! tf_v2_at_v = rho * turb_flux_v2 at v points in the faces
      do j = j_start, j_end
        do i = i_start, i_end
          tf_uv_at_v(i, kts, j) = 0.5 * (turb_flux_uv(i, kts, j - 1) + &
              turb_flux_uv(i, kts, j)) * 0.5 * (rho(i, kts, j - 1) +   &
              rho(i, kts, j))
          tf_v2_at_v(i, kts, j) = 0.5 * (turb_flux_v2(i, kts, j - 1) + &
              turb_flux_v2(i, kts, j)) * 0.5 * (rho(i, kts, j - 1) +   &
              rho(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            tf_uv_at_v(i, k, j) = 0.5 * (turb_flux_uv(i, k, j - 1) + &
                turb_flux_uv(i, k, j)) * 0.25 * (rho(i, k, j - 1) +   &
                rho(i, k, j) + rho(i, k - 1, j - 1) + rho(i, k - 1, j))
            tf_v2_at_v(i, k, j) = 0.5 * (turb_flux_v2(i, k, j - 1) + &
                turb_flux_v2(i, k, j)) * 0.25 * (rho(i, k, j - 1) +   &
                rho(i, k, j) + rho(i, k - 1, j - 1) + rho(i, k - 1, j))
              ! Horiontal derivatives of z at v points
            zx_at_v(i, k, j) = 0.125 * (zx(i, k, j) + zx(i + 1, k, j) +     &
                zx(i, k, j - 1) + zx(i + 1, k, j - 1) + zx(i, k + 1, j) +   &
                zx(i + 1, k + 1, j) + zx(i, k + 1, j - 1) + zx(i + 1, k + 1, j - 1))
            zy_at_v(i, k, j) = 0.5 * (zy(i, k, j) + zy(i, k + 1 , j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          tf_uv_at_v(i, ktf + 1,j)= 0.0
          tf_v2_at_v(i, ktf + 1, j) = 0.0
          zx_at_v(i, kts, j) = 0.125 * (zx(i, kts, j) + zx(i + 1, kts, j) +     &
              zx(i, kts, j - 1) + zx(i + 1, kts, j - 1) + zx(i, kts + 1, j) + &
              zx(i + 1, kts + 1, j) + zx(i, kts + 1, j - 1) + zx(i + 1, kts + 1, j - 1))
          zy_at_v(i, kts, j) = 0.5 * (zy(i, kts, j) + zy(i, kts + 1 , j))
        end do
      end do

        ! Calc tendency using coordinate transform (sigma to cartesian)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msfvx(i, j) * rdx
            mrdy = msfvy(i, j) * rdy
            tmpdz = (1.0 / rdzw(i, k, j) + 1.0 / rdzw(i, k, j - 1)) / 2.0
            tmpalpha = 0.5 * (sf_alpha(i, k, j - 1) + sf_alpha(i, k, j))
            tendency(i, k, j) = tendency(i, k, j) + ( G * tmpdz / dnw(k) * &
                (mrdy * (tf_v2_at_mass(i, k, j) - tf_v2_at_mass(i, k, j - 1)) +       &
                 mrdx * (tf_uv_at_corners(i + 1, k, j) - tf_uv_at_corners(i, k, j)) -       &
                 msfvx(i, j) * zx_at_v(i, k, j) * (tf_uv_at_v(i, k + 1, j) - tf_uv_at_v(i, k, j)) / tmpdz - &
                 msfvy(i, j) * zy_at_v(i, k, j) * (tf_v2_at_v(i, k + 1, j) - tf_v2_at_v(i, k, j)) / tmpdz) ) / (N_PBL3D * tmpalpha)
          end do
        end do
      end do

    end subroutine Calc_htend_dv


    subroutine Calc_htend_dw (config_flags, tendency, turb_flux_uw, turb_flux_vw, tf_uw_at_full, tf_vw_at_full, &
              msftx, msfty, rdx, rdy, dn, zx, zy, rdz, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf
      real,    intent(in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: dn
      real, dimension (ims:ime, jms:jme), intent(in) ::  msftx, msfty
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdz, rho, &
          zx, zy, turb_flux_uw, turb_flux_vw
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2), intent(in) :: &
          tf_uw_at_full, tf_vw_at_full

        ! Local data
      integer :: i, j, k
      real    :: mrdx, mrdy, tmpalpha
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: &
          tf_uw_at_mass, tf_vw_at_mass, zy_at_w, zx_at_w
      real (kind = kind(tf_uw_at_mass)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Calculate component related to slope in sigma surfaces
        ! tf_uw_at_mass = turb_flux_uw * rho at mass points
        ! tf_vw_at_mass = turb_flux_vw * rho at mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tf_uw_at_mass(i, k, j) = 0.5 * (turb_flux_uw(i, k, j) + turb_flux_uw(i, k + 1, j)) *  rho(i, k, j)
            tf_vw_at_mass(i, k, j) = 0.5 * (turb_flux_vw(i, k, j) + turb_flux_vw(i, k + 1, j)) *  rho(i, k, j)
              ! Horiontal derivatives of z at center of faces
            zx_at_w(i, k, j) = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
            zy_at_w(i, k, j) = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
          end do
        end do
      end do

        ! Lid
      do j = j_start, j_end
        do i = i_start, i_end
          tf_uw_at_mass(i, ktf + 1, j) = 0.0
          tf_vw_at_mass(i, ktf + 1, j) = 0.0
        end do
      end do

        ! Calc tendency using coordinate transform (sigma to cartesian)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            mrdx = msftx(i, j) * rdx
            mrdy = msfty(i, j) * rdy
            tmpalpha = 0.5 * (sf_alpha(i, k - 1, j) + sf_alpha(i, k, j))
            if (USE_HIGH_ORDER) then
              tendency(i, k, j) = tendency(i, k, j) + ( G / (dn(k) * rdz(i, k, j)) * (   &
                  mrdx * ((8.0 * (tf_uw_at_full(i + 1, k, j) - tf_uw_at_full(i - 1, k, j))) -  &
                  (tf_uw_at_full(i + 2, k, j) - tf_uw_at_full(i - 2, k, j))) / 12.0 +          &
                  mrdy * (8.0 * (tf_vw_at_full(i, k, j + 1) - tf_vw_at_full(i, k, j - 1)) -    &
                  (tf_vw_at_full(i, k, j + 2) - tf_vw_at_full(i, k, j - 2))) / 12.0 -          &
                  rdz(i, k, j) * (msftx(i, j) * zx_at_w(i, k, j) * (tf_uw_at_mass(i, k, j) - tf_uw_at_mass(i, k - 1, j)) +  &
                                  msfty(i, j) * zy_at_w(i, k, j) * (tf_vw_at_mass(i, k, j) - tf_vw_at_mass(i, k - 1, j))) ) ) / (N_PBL3D * tmpalpha)
            else
              tendency(i, k, j) = tendency(i, k, j) + ( G / (dn(k) * rdz(i, k, j)) * (   &
                  mrdx * (tf_uw_at_full(i + 1, k, j) - tf_uw_at_full(i, k, j)) +       &
                  mrdy * (tf_vw_at_full(i, k, j + 1) - tf_vw_at_full(i, k, j)) -       &
                  rdz(i, k, j) * (msftx(i, j) * zx_at_w(i, k, j) * (tf_uw_at_mass(i, k, j) - tf_uw_at_mass(i, k - 1, j)) +  &
                                  msfty(i, j) * zy_at_w(i, k, j) * (tf_vw_at_mass(i, k, j) - tf_vw_at_mass(i, k - 1, j))) ) ) / (N_PBL3D * tmpalpha)
            end if
          end do
        end do
      end do

    end subroutine Calc_htend_dw


    subroutine Calc_htend_ds (config_flags, tendency, turb_flux_x, turb_flux_y, tf_x_at_u, tf_y_at_v, rho, sf_alpha, &
              fnm, fnp, zx, zy, dnw, rdzw, msfty, msftx, rdx, rdy, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, &
              kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, intent (in) :: rdx, rdy
      real, dimension (kms:kme), intent (in) :: fnm, fnp, dnw
      real, dimension (ims:ime, jms:jme), intent (in) :: msftx, msfty
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: turb_flux_x, &
          turb_flux_y, rho, rdzw, zx, zy
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (in) :: tf_x_at_u, tf_y_at_v
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: tendency

      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      real :: mrdx, mrdy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ) ::  &
          tf_x_at_full, tf_y_at_full, zx_at_m, zy_at_m
      real (kind = kind(tf_x_at_full)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Calculate component related to slope in sigma surfaces
        ! tf_x_at_full = turb_flux_x * rho at the center of the faces
        ! tf_y_at_full = turb_flux_y * rho at the center of the faces
      do j = j_start, j_end
        do i = i_start, i_end
          tf_x_at_full(i, kts, j) = turb_flux_x(i, kts, j) * rho(i, kts, j)
          tf_y_at_full(i, kts, j) = turb_flux_y(i, kts, j) * rho(i, kts, j)
          end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            tf_x_at_full(i, k, j) = turb_flux_x(i, k, j) * (fnm(k) * rho(i, k, j) + &
                fnp(k) * rho(i, k - 1, j))
            tf_y_at_full(i, k, j) = turb_flux_y(i, k, j) * (fnm(k) * rho(i, k, j) + &
                fnp(k) * rho(i, k - 1, j))
            zx_at_m(i, k, j) = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + &
                zx(i, k + 1, j) + zx(i + 1, k + 1, j))
            zy_at_m(i, k, j) = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + &
                zy(i, k + 1, j) + zy(i, k + 1, j + 1))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          tf_x_at_full(i, ktf + 1, j) = 0.0
          tf_y_at_full(i, ktf + 1, j) = 0.0
          zx_at_m(i, kts, j) = 0.25 * (zx(i, kts, j) + zx(i + 1, kts, j) + &
              zx(i, kts + 1, j) + zx(i + 1, kts + 1, j))
          zy_at_m(i, kts, j) = 0.25 * (zy(i, kts, j) + zy(i, kts, j + 1) + &
              zy(i, kts + 1, j) + zy(i, kts + 1, j + 1))
        end do
      end do

        ! Calc tendency using coordinate transform (sigma to cartesian)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msftx(i, j) * rdx
            mrdy = msfty(i, j) * rdy
            tendency(i, k, j) = tendency(i, k, j) + ( G / (dnw(k) * rdzw(i, k, j)) * (  &
                mrdx * (tf_x_at_u(i + 1, k, j) - tf_x_at_u(i, k, j)) +                &
                mrdy * (tf_y_at_v(i, k, j + 1) - tf_y_at_v(i, k, j)) -                &
                msftx(i, j) * zx_at_m(i, k, j) * (tf_x_at_full(i, k + 1, j) - tf_x_at_full(i, k, j)) * rdzw(i, k, j) - &
                msfty(i, j) * zy_at_m(i, k, j) * (tf_y_at_full(i, k + 1, j) - tf_y_at_full(i, k, j)) * rdzw(i, k, j)) ) / (N_PBL3D * sf_alpha(i, k, j))
          end do
        end do
      end do

    end subroutine Calc_htend_ds


    subroutine Calc_sm_sh_l2_or_l2p5 (mass_or_face, my_level, du_dz, dv_dz, dthetav_dz, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, &
              kms, kme, its, ite, jts, jte, kts, kte, sm_l2, sh_l2, sm, sh, q_sq, q_sq_hl88, q_ratio, l_master_at_mass, l_master_at_face, rif)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                        !
      ! Purpose: Calculates Sm and Sh for level 2 model following Yamada (1975)                !
      !                                                                                        !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                                          !
      !                                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      character (len = *), intent(in)                             :: mass_or_face
      character (len = *), intent(in)                             :: my_level
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz
      real, dimension (its:ite, kts:kte, jts:jte), intent (out)   :: sm, sh
      real, dimension (its:ite, kts:kte, jts:jte), intent (out), optional :: rif, q_ratio
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in),  optional :: q_sq
      real, dimension (its:ite, kts:kte, jts:jte), intent (in),  optional :: sm_l2, sh_l2, q_sq_hl88
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in),  optional :: l_master_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in),  optional :: l_master_at_face
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf, k_start, n_iter_sor
      real :: ri, turb_prtl_num, gm, gh, phi1, phi2, phi3, phi4, phi5, N, D
      real, dimension (its:ite, kts:kte, jts:jte) :: l_master

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      level: if (my_level == 'l2') then   ! level 2 model

        if (mass_or_face == 'face') then
          k_start = kts + 1
        else
          k_start = kts
        end if

        j_loop_l2: do j = j_start, j_end
          k_loop_l2: do k = k_start, ktf
            i_loop_l2: do i = i_start, i_end

              if (abs(du_dz(i, k, j)) < DWIND_DZ_MIN .and. abs(dv_dz(i, k, j)) < DWIND_DZ_MIN) then
                ri = G_OVER_TREF * dthetav_dz(i, k, j) / (2.0 * DWIND_DZ_MIN * DWIND_DZ_MIN)
              else
                ri = G_OVER_TREF * dthetav_dz(i, k, j) / (du_dz(i, k, j) ** 2.0 + dv_dz(i, k, j) ** 2.0)
              end if

                ! Eq. on page 929 Y75
              rif(i, k, j) = fac1 * (ri + fac2 - sqrt(ri * ri + fac3 * ri + fac4))
              if (rif(i, k, j) /= rif(i, k, j)) then
                write(90, *) 'i =', i, 'j = ', j, 'k = ', k
                write(90, *) '    dt_dz = ', dthetav_dz(i, k, j)
                write(90, *) '    du_dz = ', du_dz(i, k, j)
                write(90, *) '    dv_dz = ', dv_dz(i, k, j)
              end if

              if_turb: if (rif(i, k, j) < rif_c) then

                  ! Eq. 9a Y75
                sm(i, k, j) = cm * (rif_c - rif(i, k, j)) * (rf1 - rif(i, k, j)) / ((1.0 - rif(i, k, j)) * (rf2 - rif(i, k, j)))
                  ! Eq. 11 Y75
                turb_prtl_num = (ch2 / cm) * (rf2 - rif(i, k, j)) / (rf1 - rif(i, k, j))
                sh(i, k, j) = turb_prtl_num * sm(i, k, j)

              else

                sm(i, k, j) = TURB_FLUX_MIN
                sh(i, k, j) = TURB_FLUX_MIN

              end if if_turb

            end do i_loop_l2
          end do k_loop_l2
        end do j_loop_l2

        do j = j_start, j_end
          do i = i_start, i_end
            if (mass_or_face == 'face') then
              sm(i, kts, j) = TURB_FLUX_MIN
              sh(i, kts, j) = TURB_FLUX_MIN
            end if
            sm(i, ktf + 1, j) = TURB_FLUX_MIN
            sh(i, ktf + 1, j) = TURB_FLUX_MIN
          end do
        end do

      else if (my_level == 'l2p5') then   ! level 2.5 model

        if (PRESENT(l_master_at_mass)) then

          k_start = kts

          do j = j_start, j_end
            do k = k_start, ktf
              do i = i_start, i_end
                l_master(i, k, j) = l_master_at_mass(i, k, j)
              end do
            end do
          end do

        else if (PRESENT(l_master_at_face)) then

          k_start = kts + 1

          do j = j_start, j_end
            do k = k_start, ktf
              do i = i_start, i_end
                l_master(i, k, j) = l_master_at_face(i, k, j)
              end do
            end do
          end do

        end if

        j_loop_l2p5: do j = j_start, j_end
          k_loop_l2p5: do k = k_start, ktf
            i_loop_l2p5: do i = i_start, i_end

                  ! Realizability criteria from HL88
                if ( q_sq(i, k, j) < q_sq_hl88(i, k, j) ) then
                  q_ratio(i, k, j) = sqrt( q_sq(i, k, j) / q_sq_hl88(i, k, j) )
                  sm(i, k, j) = sm_l2(i, k, j) * q_ratio(i, k, j)
                  sh(i, k, j) = sh_l2(i, k, j) * q_ratio(i, k, j)
                else
                  q_ratio(i, k, j) = sqrt( q_sq(i, k, j) / q_sq_hl88(i, k, j) )
                    ! Eqs. 33a and 33b on page 855 MY82
                  gh = - ( (l_master(i, k, j) ** 2.0) / Max (q_sq(i, k, j), Q_SQ_MIN) ) * G_OVER_TREF * dthetav_dz(i, k, j)
                  gm = ( (l_master(i, k, j) ** 2.0) / Max (q_sq(i, k, j), Q_SQ_MIN) ) * ( Max( du_dz(i, k, j) ** 2.0 + dv_dz(i, k, j) ** 2.0, 1.0E-10 ) )

                    ! Eqs. 2.31, 2.33, 2.35 - 2.41 from HL88
                  phi1 = 1.0 - 3.0 * a2 * b2 * gh
                  phi2 = 1.0 - 9.0  * a1 * a2 * gh
                  phi3 = phi1 + 9.0 * a2 * a2 * gh
                  phi4 = phi1 - 12.0 * a1 * a2 * gh
                  phi5 = phi1 - 18.0 * a1 * a2 * gh
                  N = phi2 + 18.0 * c1 * a1 * a1 * gm
                  D = phi2 * phi4 + 6.0 * a1 * a1 * phi3 * gm

                  sm(i, k, j) = a1 * (phi3 - 3.0 * c1 * phi4) / D
                  sh(i, k, j) = a2 * N / D

                    ! Eqs. 3.17, 3.19, 3.20 from HL88
                  sm(i, k, j) = max ( min ( sm(i, k, j), sqrt ( 0.44 ** 2.0 / gm ) ), 0.0 )
                  sh(i, k, j) = max ( min ( sh(i, k, j), 0.76 * b2 ), 0.0 )
                end if

            end do i_loop_l2p5
          end do k_loop_l2p5
        end do j_loop_l2p5

        do j = j_start, j_end
          do i = i_start, i_end
            if (PRESENT(l_master_at_face)) then
              sm(i, kts, j) = 0.0 !TURB_FLUX_MIN
              sh(i, kts, j) = 0.0 !TURB_FLUX_MIN
            end if
            sm(i, ktf + 1, j) = 0.0 !TURB_FLUX_MIN
            sh(i, ktf + 1, j) = 0.0 !TURB_FLUX_MIN
          end do
        end do

      end if level

    end subroutine Calc_sm_sh_l2_or_l2p5


    subroutine Calc_q_sq_l2_pbl_approx (mass_or_face, l_master, du_dz, dv_dz, dthetav_dz, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, q_sq_face, q_sq_hl88, sm, rif)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                        !
      ! Purpose: Calculates diagnostic TKE equation (level 2)                                  !
      !                                                                                        !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                                          !
      !                                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      character (len = *), intent(in)                             :: mass_or_face
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out), optional :: q_sq_face
      real, dimension (its:ite, kts:kte, jts:jte), intent (out), optional :: q_sq_hl88
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: sm, rif
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, kts1, ktf
      real :: ri, turb_prtl_num
      real, dimension (its:ite, kts:kte, jts:jte) :: q_sq
      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      if (mass_or_face == 'mass') then
        kts1 = kts
      else
        kts1 = kts + 1
      end if

      ktf = Min (kte, kde - 1)

      j_loop: do j = j_start, j_end
        k_loop: do k = kts1, ktf
          i_loop: do i = i_start, i_end

            if_turb: if (rif(i, k, j) < rif_c) then

              q_sq(i, k, j) =  b1 * l_master(i, k, j) ** 2.0 *        &
                  (du_dz(i, k, j) ** 2.0 + dv_dz(i, k, j) ** 2.0) *  &
                  (1.0 - rif(i, k, j)) * sm(i, k, j)

              q_sq(i, k, j) = Max (q_sq(i, k, j), Q_SQ_MIN)

              if (DEBUG) call Debug_big_qsq

            else

              q_sq(i, k, j) = Q_SQ_MIN

            end if if_turb

          end do i_loop
        end do k_loop
      end do j_loop

        ! Zero q_sq at the suface and top
      do j = j_start, j_end
        do i = i_start, i_end
          if (mass_or_face == 'face') then
            q_sq(i, kts, j) = Q_SQ_MIN
          end if
          q_sq(i, ktf + 1, j) = Q_SQ_MIN
        end do
      end do

      if (PRESENT(q_sq_face)) then

        do j = j_start, j_end
          do k = kts, ktf + 1
            do i = i_start, i_end
              q_sq_face(i, k, j) = q_sq(i, k, j)
            end do
          end do
        end do

      end if

      if (PRESENT(q_sq_hl88)) then

        do j = j_start, j_end
          do k = kts, ktf + 1
            do i = i_start, i_end
              q_sq_hl88(i, k, j) = q_sq(i, k, j)
            end do
          end do
        end do

      end if

      contains

        subroutine Debug_big_qsq

          if (q_sq(i, k, j) > 20.0) then
            print *, '-----------------------------'
            print *, 'i, k, j = ', i, k, j
            print *, 'q_sq = ', q_sq(i, k, j)
            print *, 'du_dz = ', du_dz(i, k, j)
            print *, 'dv_dz = ', dv_dz(i, k, j)
            print *, 'dthetav_dz = ', dthetav_dz(i, k, j)
            print *, 'l master = ', l_master(i, k, j)
            print *, 'Ri = ', ri
            print *, 'Ri flux = ', rif(i, k, j)
            print *, 'sm = ', sm(i, k, j)
          end if

        end subroutine Debug_big_qsq

    end subroutine Calc_q_sq_l2_pbl_approx


    subroutine Calc_q_sq_l2p5 (config_flags, u, v, w, zx, zy, t, qv, rho, thetav, tsk, t2, q2, psfc, ust, hfx, q_sq, q_sq_prog, q_sq_tend, &
              q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, q_sq_buoyancy_tend, q_sq_dissip_tend, l_master_at_face, &
              l_master_at_mass, turb_flux_u2_face, turb_flux_v2_face, turb_flux_w2_face, &
              turb_flux_uv_face, turb_flux_uw_face, turb_flux_vw_face, turb_flux_wtheta_face, turb_flux_wtheta_v_face, &
              turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta, &
              turb_flux_wtheta_v, dn, dnw, fnm, fnp, msftx, msfty, msfux, msfuy, msfvx, msfvy, mut, c1, c2, &
              z_at_mass, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, dt, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, &
              du_dx_mass, du_dy_mass, du_dz_mass, dv_dx_mass, dv_dy_mass, dv_dz_mass, dw_dx_mass, dw_dy_mass, dw_dz_mass, &
              dthetav_dx_mass, dthetav_dy_mass, dthetav_dz_mass, dqv_dx_mass, dqv_dy_mass, dqv_dz_mass)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                        !
      ! Purpose: Calculates derivatives and fluxes at mass points                              !
      !          to prepare for prognostic TKE equation (level 2.5)                            !
      !                                                                                        !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                                          !
      !                                                                                        !
      ! Method: 1) Follows module_diffusion_em to calc shear, buoyancy, and dissipation terms  !
      !                                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: u, v, w, zx, zy, t, qv, rho, thetav
      real, dimension (ims:ime, jms:jme),          intent (in)    :: tsk, t2, q2, psfc, ust, hfx
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_vdiff_tend,    &
                                                                     q_sq_hdiff_tend,    &
                                                                     q_sq_shear_tend,    &
                                                                     q_sq_buoyancy_tend, &
                                                                     q_sq_dissip_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master_at_face
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: l_master_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                                                                     turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                                                                     turb_flux_wtheta, turb_flux_wtheta_v
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: turb_flux_u2_face, turb_flux_v2_face, turb_flux_w2_face, &
                                                                     turb_flux_uv_face, turb_flux_uw_face, turb_flux_vw_face, &
                                                                     turb_flux_wtheta_face, turb_flux_wtheta_v_face
      real, dimension (kms:kme),                   intent (in)    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)    :: msftx, msfty, msfux, msfuy, msfvx, msfvy
      real, dimension (ims:ime, jms:jme),          intent (in)    :: mut
      real, dimension (kms:kme),                   intent (in)    :: c1, c2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: z_at_mass, rdz, rdzw
      real,                                        intent (in)    :: rdx, rdy, cf1, cf2, cf3, dt
      integer, intent (in) :: ids, ide, jds, jde, kds, kde,  &
                              ims, ime, jms, jme, kms, kme,  &
                              its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: du_dx_mass, du_dy_mass, du_dz_mass, &
          dv_dx_mass, dv_dy_mass, dv_dz_mass, dw_dx_mass, dw_dy_mass, dw_dz_mass, &
          dthetav_dx_mass, dthetav_dy_mass, dthetav_dz_mass, dqv_dx_mass, dqv_dy_mass, dqv_dz_mass

        ! Local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real, dimension (its:ite, kts:kte, jts:jte) :: sm_l2, sh_l2, sm, sh, rif, q_sq_hl88, q_ratio
      real, dimension (its:ite, kts:kte, jts:jte) :: du_dz, dv_dz, dw_dz, du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, &
           dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz
      real, dimension (its:ite, kts:kte, jts:jte) :: turb_flux_wqv
      real, dimension (its:ite, jts:jte)          :: thetav_sk
      logical, parameter :: use_hl88 = .true.
      character (4) :: my_level, mass_or_face

      ktf     = MIN( kte, kde-1 )
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      call Calc_wind_derivatives_at_mass (config_flags, u, v, w, zx, zy, du_dz, dv_dz, dw_dz, &
          du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, rdzw, dn, dnw, fnm, fnp, msftx, msfty, &
          msfux, msfuy, msfvx, msfvy, rdx, rdy, cf1, cf2, cf3, ids, ide, jds, jde, kds, kde, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Calc_surface_thetav (tsk, t2, q2, psfc, thetav_sk, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Calc_mass_var_derivatives_at_mass (thetav, thetav_sk, q2, qv, dthetav_dx, dthetav_dy, dthetav_dz, &
          dqv_dx, dqv_dy, dqv_dz, rdzw, zx, zy, msftx, msfty, dn, dnw, fnm, fnp, rdx, rdy, cf1, cf2, cf3, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! Save derivatives
      do j = j_start, j_end
        do k = kts, ktf + 1
          do i = i_start, i_end
           du_dx_mass(i, k, j) = du_dx(i, k, j)
           du_dy_mass(i, k, j) = du_dy(i, k, j)
           du_dz_mass(i, k, j) = du_dz(i, k, j)
           dv_dx_mass(i, k, j) = dv_dx(i, k, j)
           dv_dy_mass(i, k, j) = dv_dy(i, k, j)
           dv_dz_mass(i, k, j) = dv_dz(i, k, j)
           dw_dx_mass(i, k, j) = dw_dx(i, k, j)
           dw_dy_mass(i, k, j) = dw_dy(i, k, j)
           dw_dz_mass(i, k, j) = dw_dz(i, k, j)
           dthetav_dx_mass(i, k, j) = dthetav_dx(i, k, j)
           dthetav_dy_mass(i, k, j) = dthetav_dy(i, k, j)
           dthetav_dz_mass(i, k, j) = dthetav_dz(i, k, j)
           dqv_dx_mass(i, k, j) = dqv_dx(i, k, j)
           dqv_dy_mass(i, k, j) = dqv_dy(i, k, j)
           dqv_dz_mass(i, k, j) = dqv_dz(i, k, j)
          end do
        end do
      end do

      if (use_hl88) then
        mass_or_face = 'mass'
        my_level = 'l2'
        call Calc_sm_sh_l2_or_l2p5 (mass_or_face=mass_or_face, my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
            ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm_l2, sh=sh_l2, rif=rif)
        call Calc_q_sq_l2_pbl_approx (mass_or_face=mass_or_face, l_master=l_master_at_mass, du_dz=du_dz, dv_dz=dv_dz, &
            dthetav_dz=dthetav_dz, ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, &
            kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, q_sq_hl88=q_sq_hl88, sm=sm_l2, rif=rif)
      end if

        ! Calculate Sm and Sh for level 2.5 according to Mellor and Yamada (1984)
      my_level = 'l2p5'
      call Calc_sm_sh_l2_or_l2p5 (mass_or_face=mass_or_face, my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm_l2=sm_l2, sh_l2=sh_l2, sm=sm, sh=sh, &
          q_sq=q_sq_prog, q_sq_hl88=q_sq_hl88, q_ratio=q_ratio, l_master_at_mass=l_master_at_mass)

        ! Diagnose master length scale at mass points for prognostic TKE
      my_level = 'l2p5'
      call Calc_l_my_algebra_l2_or_l2p5 (my_level=my_level, q_sq=q_sq_prog, z=z_at_mass, l_master_l2_in=l_master_at_face, &
          l_master_l2p5=l_master_at_mass, ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
          ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

      if ( config_flags%pbl3d_opt .lt. 2 ) then
        call Diagnose_fluxes_l2p5_pbl_approx (u, v, t, qv, rho, ust, hfx, du_dz, dv_dz, dthetav_dz, dqv_dz, &
            q_sq_prog, q_sq_hl88, q_ratio, sm, sh, sm_l2, sh_l2, rif, l_master_at_mass, &
            turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, &
            turb_flux_wtheta, turb_flux_wtheta_face, turb_flux_wtheta_v, turb_flux_wtheta_v_face, turb_flux_wqv, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      else
        call Diagnose_fluxes_l2p5 (du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
            dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, &
            turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta, turb_flux_wtheta_v, &
            turb_flux_wtheta_v_face, turb_flux_wqv, q_sq_prog, q_sq_hl88, q_ratio, l_master_at_mass, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

        ! Calc TKE tendency
      call Calc_q_sq_rhs (config_flags, sm, sh, zx, zy, du_dz, dv_dz, dw_dz, du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, &
          dthetav_dz, dqv_dz, q_sq, q_sq_prog, q_sq_tend, q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, &
          q_sq_buoyancy_tend, q_sq_dissip_tend, l_master_at_face, l_master_at_mass, turb_flux_u2, turb_flux_v2, turb_flux_w2, &
          turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta, turb_flux_wtheta_v, dn, dnw, fnm, fnp, &
          msftx, msfty, msfux, msfuy, msfvx, msfvy, mut, c1, c2, z_at_mass, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, dt, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

    end subroutine Calc_q_sq_l2p5


    subroutine Calc_wind_derivatives_at_mass (config_flags, u, v, w, zx, zy, du_dz, dv_dz, dw_dz, &
              du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, rdzw, dn, dnw, fnm, fnp, msftx, msfty, &
              msfux, msfuy, msfvx, msfvy, rdx, rdy, cf1, cf2, cf3, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                                         !
        ! Purpose: Calculate spatial derivatives of the wind variables at mass points             !
        !                                                                                         !
        ! Authors: Timothy W. Juliano & Pedro A. Jimenez                                          !
        !                                                                                         !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, v, w, zx, zy, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfty, msfux, msfuy, msfvx, msfvy
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dz, dv_dz, dw_dz, &
                                                                   du_dx, dv_dx, dw_dx, &
                                                                   du_dy, dv_dy, dw_dy
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3

      logical, parameter :: DEBUG = .false.

      call Calc_du_dz_at_mass (u, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      call Calc_dv_dz_at_mass (v, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      call Calc_dw_dz_at_mass (w, rdzw, dn, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dw_dz)

      call Calc_du_dx_at_mass (u, zx, rdzw, dn, dnw, fnm, fnp, msftx, msfux, rdx, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dx)

      call Calc_dv_dx_at_mass (v, zx, rdzw, dn, dnw, fnm, fnp, msftx, msfvx, rdx, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dx)

      call Calc_dw_dx_at_mass (w, zx, rdzw, msftx, rdx, ims, ime, jms, jme, kms, kme,  &
          ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dx)

      call Calc_du_dy_at_mass (u, zy, rdzw, dn, dnw, fnm, fnp, msfty, msfuy, rdy, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dy)

      call Calc_dv_dy_at_mass (v, zy, rdzw, dn, dnw, fnm, fnp, msfty, msfvy, rdy, cf1, cf2, cf3, ide, jds, &
              jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dy)

      call Calc_dw_dy_at_mass (w, zy, rdzw, msfty, rdy, ims, ime, jms, jme, kms, kme, &
              ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dy)

      if (DEBUG) print *, 'Done Calc_wind_derivatives_at_mass...'

    end subroutine Calc_wind_derivatives_at_mass


    subroutine Calc_du_dz_at_mass (u, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dz                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dz
      real,                                        intent (in)  :: cf1, cf2, cf3

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) :: u_at_full
      real, dimension (its:ite, kts:kte, jts:jte)   :: u_at_face
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            u_at_full(i, k, j) = fnm(k) * u(i, k, j) + fnp(k) * u(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          u_at_full(i, 1, j) =  0.0
            ! Top face
          u_at_full(i, kte, j) =  cft1 * u(i, ktes1, j) + cft2 * u(i, ktes2, j)
        end do
      end do

        ! Calculate u at the center of the faces
      do j = j_start, j_end
        do k = kts, kte
          do i = i_start, i_end
            u_at_face(i, k, j) = 0.5 * (u_at_full(i, k, j) + u_at_full(i + 1, k, j))
          end do
        end do
      end do

        ! Calculate derivative at the mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            du_dz(i, k, j) = (u_at_face(i, k + 1, j) - u_at_face(i, k, j)) * rdzw(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          du_dz(i, kts, j) = 0.0
            ! Top
          du_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dz_at_mass


    subroutine Calc_dv_dz_at_mass (v, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dz                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dz
      real,                                        intent (in)  :: cf1, cf2, cf3

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) :: v_at_full
      real, dimension (its:ite, kts:kte, jts:jte)   :: v_at_face
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            v_at_full(i, k, j) = fnm(k) * v(i, k, j) + fnp(k) * v(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          v_at_full(i, 1, j) =  0.0
            ! Top face
          v_at_full(i, kte, j) =  cft1 * v(i, ktes1, j) + cft2 * v(i, ktes2, j)
        end do
      end do

        ! Calculate v at the center of the faces
      do j = j_start, j_end
        do k = kts, kte
          do i = i_start, i_end
            v_at_face(i, k, j) = 0.5 * (v_at_full(i, k, j) + v_at_full(i, k, j + 1))
          end do
        end do
      end do

        ! Calculate derivative at the mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dv_dz(i, k, j) = (v_at_face(i, k + 1, j) - v_at_face(i, k, j)) * rdzw(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          dv_dz(i, kts, j) = 0.0
            ! Top
          dv_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dv_dz_at_mass


    subroutine Calc_dw_dz_at_mass (w, rdzw, dn, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dw_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dz                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dz

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Calculate derivative at the mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dw_dz(i, k, j) = (w(i, k + 1, j) - w(i, k, j)) * rdzw(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          dw_dz(i, kts, j) = 0.0
            ! Top
          dw_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dz_at_mass


    subroutine Calc_du_dx_at_mass (u, zx, rdzw, dn, dnw, fnm, fnp, msftx, msfux, rdx, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dx)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dx                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, zx, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfux
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dx

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2, dz_dx
      real, dimension (its:ite, jts:jte) :: msf_tx_ux
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) :: du_dz_dz_dx, u_mux, u_mux_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

      logical, parameter :: DEBUG = .false.

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at u points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_tx_ux(i, j) = msftx(i, j) * msfux(i, j)
        end do
      end do

        ! Scale u component of vector to account for map projection
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            u_mux(i, k, j) = u(i, k, j) / msfux(i, j)
          end do
        end do
      end do

        ! Interpolate u_mux from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            u_mux_at_full(i, k, j) = fnm(k) * u_mux(i, k, j) + fnp(k) * u_mux(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          u_mux_at_full(i, 1, j) =  cf1 * u_mux(i, 1, j) + cf2 * u_mux(i, 2, j) + cf3 * u_mux(i, 3, j)
            ! Top face
          u_mux_at_full(i, kte, j) =  cft1 *  u_mux(i, ktes1, j) +  cft2 * u_mux(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (du_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at mass points
            dz_dx = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
              ! du_dz_dz_dx = (partial u_mux / partial z)|sigma (du_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            du_dz_dz_dx(i, k, j) = (0.5 * (u_mux_at_full(i, k + 1, j) + u_mux_at_full(i + 1, k + 1, j)) - &
                0.5 * (u_mux_at_full(i, k, j) + u_mux_at_full(i + 1, k, j))) * rdzw(i, k, j) * dz_dx

            if (DEBUG) then
              if ((i == 49 .or. i == 50 .or. i == 51) .and. (k == 3 .or. k == 4) .and. j == 54) then
                print *, '-----------------------------'
                print *, 'i = ', i, 'j = ', j, 'k = ', k
                print *, 'dz_dx: '
                print *, dz_dx
                print *, 'du_dz_dz_dx: '
                print *, du_dz_dz_dx(i, k, j)
              end if
            end if
          end do
        end do
      end do

        ! (partial u / partial x)|z = msf_tx_ux * (partial u_mux_at_full / partial x)|sigma - du_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            du_dx(i, k, j) = msf_tx_ux(i, j) * ((u_mux(i + 1, k, j) - u_mux(i, k, j)) * rdx - du_dz_dz_dx(i, k, j))

            if (DEBUG) then
              if ((i == 49 .or. i == 50 .or. i == 51) .and. (k == 3 .or. k == 4) .and. j == 54) then
                print *, '-----------------------------'
                print *, 'i = ', i, 'j = ', j, 'k = ', k
                print *, 'u_mux_at_full: '
                print *, u_mux_at_full(i + 1, k, j), u_mux_at_full(i, k, j)
                print *, 'rdx: '
                print *, rdx
                print *, 'du_dx: '
                print *, du_dx(i, k, j)
              end if
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          du_dx(i, kts, j) = 0.0
            ! Top
          du_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dx_at_mass


    subroutine Calc_dv_dx_at_mass (v, zx, rdzw, dn, dnw, fnm, fnp, msftx, msfvx, rdx, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dx)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dx                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, zx, rdzw
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfvx
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dx

        ! Local vars
      integer :: i, j, k
      real    :: dz_dx, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_tx_vx
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: dv_dz_dz_dx, v_mvx, v_mvx_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_tx_vx(i, j) = msftx(i, j) * msfvx(i, j)
        end do
      end do

        ! Scale v component of vector to account for map projection (at u points)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            v_mvx(i, k, j) = 0.25 * (v(i, k, j) + v(i, k, j + 1) + v(i - 1, k, j) + v(i - 1, k, j + 1)) / msfvx(i, j)
          end do
        end do
      end do

        ! Interpolate v_mvx from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            v_mvx_at_full(i, k, j) = fnm(k) * v_mvx(i, k, j) + fnp(k) * v_mvx(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          v_mvx_at_full(i, 1, j) = cf1 * v_mvx(i, 1, j) + cf2 * v_mvx(i, 2, j) + cf3 * v_mvx(i, 3, j)
            ! Top face
          v_mvx_at_full(i, kte, j) = cft1 * v_mvx(i, ktes1, j) + cft2 * v_mvx(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dv_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at mass points
            dz_dx = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
              ! dv_dz_dz_dx = (partial v_mvx / partial z)|sigma (dv_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dv_dz_dz_dx(i, k, j) = (0.5 * (v_mvx_at_full(i, k + 1, j) + v_mvx_at_full(i + 1, k + 1, j)) - &
                0.5 * (v_mvx_at_full(i, k, j) + v_mvx_at_full(i + 1, k, j))) * rdzw(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial v / partial x)|z = msf_tx_vx * (partial v_mvx_at_full / partial x)|sigma - dv_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dv_dx(i, k, j) = msf_tx_vx(i, j) * ((v_mvx(i + 1, k, j) - v_mvx(i, k, j)) * rdx - dv_dz_dz_dx(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!          dv_dx(i, 1, j) =  0.0
          dv_dx(i, kte, j) =  0.0
        end do
      end do

    end subroutine Calc_dv_dx_at_mass


    subroutine Calc_dw_dx_at_mass (w, zx, rdzw, msftx, rdx, ims, ime, jms, jme, kms, kme,  &
              ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dx)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dx                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, zx, rdzw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx
      real,                                        intent (in)  :: rdx
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dx

        ! Local vars
      integer :: i, j, k, i_offset
      real    :: dz_dx
      real, dimension (its:ite, jts:jte) :: msf_tx_tx
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: dw_dz_dz_dx, w_mtx, w_mtx_at_half
      integer :: ktf, i_start, i_end, j_start, j_end

      logical, parameter :: DEBUG = .false.

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points squared
      do j = jts, jte
        do i = its, ite
          msf_tx_tx(i, j) = msftx(i, j) * msftx(i, j)
        end do
      end do

        ! Scale w component of vector to account for map projection
      if (USE_HIGH_ORDER) then
        i_offset = 2
      else
        i_offset = 1
      end if
      do j = j_start, j_end
        do k = kts, kte
          do i = i_start - i_offset, i_end + i_offset
            w_mtx(i, k, j) = w(i, k, j) / msftx(i, j)
          end do
        end do
      end do

        ! Average w_mtx to the mass points (half levels)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start - i_offset, i_end + i_offset
            w_mtx_at_half(i, k, j) = 0.5 * (w_mtx(i, k, j) + w_mtx(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dw_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at mass points
            dz_dx = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
              ! dw_dz_dz_dx = (partial w_mtx / partial z)|sigma (dw_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dw_dz_dz_dx(i, k, j) = (w_mtx(i, k + 1, j) - w_mtx(i, k, j)) * rdzw(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial w / partial x)|z = msf_tx_tx * (partial w_mtx / partial x)|sigma - dw_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end

            if (DEBUG) then
              if (i == 3 .and. k == 2 .and. j == 3) then
                print *, '-----------------------------'
                print *, 'ims = ', ims, 'ime = ', ime
                print *, 'i_start = ', i_start, 'i_end = ', i_end
                print *, 'jms = ', jms, 'jme = ', jme
                print *, 'j_start = ', j_start, 'j_end = ', j_end
                print *, 'Row: '
                print *, w(:, k, j)
                print *, 'Column: '
                print *, w(i, k, :)
              end if
            end if

            if (USE_HIGH_ORDER) then
              dw_dx(i, k, j) = msf_tx_tx(i, j) * (rdx * (8.0 * (w_mtx_at_half(i + 1, k, j) - w_mtx_at_half(i - 1, k, j)) - (w_mtx_at_half(i + 2, k, j) - w_mtx_at_half(i - 2, k, j)) ) / 12.0 - dw_dz_dz_dx(i, k, j))
            else
              dw_dx(i, k, j) = msf_tx_tx(i, j) * (rdx * 0.5 * (w_mtx_at_half(i + 1, k, j) - w_mtx_at_half(i - 1, k, j)) - dw_dz_dz_dx(i, k, j))
            end if

          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!          dw_dx(i, 1, j) = 0.0
          dw_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dx_at_mass


    subroutine Calc_du_dy_at_mass (u, zy, rdzw, dn, dnw, fnm, fnp, msfty, msfuy, rdy, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dy)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dy                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, zy, rdzw
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty, msfuy
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dy

        ! Local vars
      integer :: i, j, k
      real    :: dz_dy, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_ty_uy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: du_dz_dz_dy, u_muy, u_muy_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

      logical, parameter :: DEBUG = .false.

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_ty_uy(i, j) = msfty(i, j) * msfuy(i, j)
        end do
      end do

        ! Scale u component of vector to account for map projection (at v points)
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            u_muy(i, k, j) = 0.25 * (u(i, k, j) + u(i + 1, k, j) + u(i, k, j - 1) + u(i + 1, k, j - 1)) / msfuy(i, j)
          end do
        end do
      end do

        ! Interpolate u_muy from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            u_muy_at_full(i, k, j) = fnm(k) * u_muy(i, k, j) + fnp(k) * u_muy(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          u_muy_at_full(i, 1, j) =  cf1 * u_muy(i, 1, j) + cf2 * u_muy(i, 2, j) + cf3 * u_muy(i, 3, j)
            ! Top face
          u_muy_at_full(i, kte, j) = cft1 * u_muy(i, ktes1, j) +  cft2 * u_muy(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (du_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at mass points
            dz_dy = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
              ! du_dz_dz_dy = (partial u_muy / partial z)|sigma (du_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            du_dz_dz_dy(i, k, j) = (0.5 * (u_muy_at_full(i, k + 1, j) + u_muy_at_full(i, k + 1, j + 1)) - &
                0.5 * (u_muy_at_full(i, k, j) + u_muy_at_full(i, k, j + 1))) * rdzw(i, k, j) * dz_dy

            if (DEBUG) then
              if (i == 54 .and. (k == 3 .or. k == 4) .and. (j == 59 .or. j == 60 .or. j == 61 .or. j == 62)) then
                print *, '-----------------------------'
                print *, 'i = ', i, 'j = ', j, 'k = ', k
                print *, 'dz_dy: '
                print *, dz_dy
                print *, 'du_dz_dz_dy: '
                print *, du_dz_dz_dy(i, k, j)
              end if
            end if

          end do
        end do
      end do

        ! (partial u / partial y)|z = msf_ty_uy * (partial u_muy_at_full / partial y)|sigma - du_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            du_dy(i, k, j) = msf_ty_uy(i, j) * ((u_muy(i, k, j + 1) - u_muy(i, k, j)) * rdy - du_dz_dz_dy(i, k, j))

            if (DEBUG) then
              if (i == 54 .and. (k == 3 .or. k == 4) .and. (j == 59 .or. j == 60 .or. j == 61 .or. j == 62)) then
                print *, '-----------------------------'
                print *, 'i = ', i, 'j = ', j, 'k = ', k
                print *, 'u_muy_at_full: '
                print *, u_muy_at_full(i, k, j + 1), u_muy_at_full(i, k, j)
                print *, 'rdy: '
                print *, rdy
                print *, 'du_dy: '
                print *, du_dy(i, k, j)
              end if
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
!          du_dy(i, kts, j) = 0.0
            ! Top face
          du_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dy_at_mass


    subroutine Calc_dv_dy_at_mass (v, zy, rdzw, dn, dnw, fnm, fnp, msfty, msfvy, rdy, cf1, cf2, cf3, ide, jds, &
              jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dy)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dy                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, zy, rdzw
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty, msfvy
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jds, jde, kde,            &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte
      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dy

        ! Local vars
      integer :: i, j, k
      real    :: dz_dy, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_ty_vy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: dv_dz_dz_dy, v_mvy, v_mvy_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_ty_vy(i, j) = msfty(i, j) * msfvy(i, j)
        end do
      end do

        ! Scale v component of vector to account for map projection
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
              ! Because msfvy at the poles will be undefined (1./0.),
              ! we will have trouble. But we are OK since v at the
              ! poles is 0., and that takes precedence in this case
            if ((config_flags%polar) .and. ((j == jds) .or. (j == jde))) then
              v_mvy(i, k, j) = 0.0
            else
              v_mvy(i, k, j) = v(i, k, j) / msfvy(i, j)
            end if
          end do
        end do
      end do

        ! Interpolate v_mvy from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            v_mvy_at_full(i, k, j) = fnm(k) * v_mvy(i, k, j) + fnp(k) * v_mvy(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          v_mvy_at_full(i, 1, j) =  cf1 * v_mvy(i, 1, j) + cf2 * v_mvy(i, 2, j) + cf3 * v_mvy(i, 3, j)
            ! Top face
          v_mvy_at_full(i, kte, j) = cft1 * v_mvy(i, ktes1, j) + cft2 * v_mvy(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dv_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at mass points
            dz_dy = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
              ! dv_dz_dz_dy = (partial v_mvy / partial z)|sigma (dv_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dv_dz_dz_dy(i, k, j) = (0.5 * (v_mvy_at_full(i, k + 1, j) + v_mvy_at_full(i, k + 1, j + 1)) - &
                0.5 * (v_mvy_at_full(i, k, j) + v_mvy_at_full(i, k, j + 1))) * rdzw(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial v / partial y)|z = msf_ty_vy * (partial v_mvy_at_full / partial y)|sigma - dv_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dv_dy(i, k, j) = msf_ty_vy(i, j) * ((v_mvy(i, k, j + 1) - v_mvy(i, k, j)) * rdy - dv_dz_dz_dy(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
!          dv_dy(i, 1, j) = 0.0
            ! Top face
          dv_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dv_dy_at_mass


    subroutine Calc_dw_dy_at_mass (w, zy, rdzw, msfty, rdy, ims, ime, jms, jme, kms, kme, &
              ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dy)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dy                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, zy, rdzw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty
      real,                                        intent (in)  :: rdy
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dy

        ! Local vars
      integer :: i, j, k, j_offset
      real    :: dz_dy
      real, dimension (its:ite, jts:jte) :: msf_ty_ty
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: dw_dz_dz_dy, w_mty, w_mty_at_half
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points squared
      do j = jts, jte
        do i = its, ite
          msf_ty_ty(i, j) = msfty(i, j) * msfty(i, j)
        end do
      end do

        ! Scale w component of vector to account for map projection
      if (USE_HIGH_ORDER) then
        j_offset = 2
      else
        j_offset = 1
      end if
      do j = j_start - j_offset, j_end + j_offset
        do k = kts, kte
          do i = i_start, i_end
            w_mty(i, k, j) = w(i, k, j) / msfty(i, j)
          end do
        end do
      end do

        ! Average w_mty to the mass points (half levels)
      do j = j_start - j_offset, j_end + j_offset
        do k = kts, ktf
          do i = i_start, i_end
            w_mty_at_half(i,k,j) = 0.5 * (w_mty(i, k, j) + w_mty(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dw_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at mass points
            dz_dy = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
              ! dw_dz_dz_dy = (partial w_mty / partial z)|sigma (dw_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dw_dz_dz_dy(i, k, j) = (w_mty(i, k + 1, j) - w_mty(i, k, j)) * rdzw(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial w / partial y)|z = msf_ty_ty * (partial w_mty / partial y)|sigma - dw_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            if (USE_HIGH_ORDER) then
              dw_dy(i, k, j) = msf_ty_ty(i, j) * (rdy * (8.0 * (w_mty_at_half(i, k, j + 1) - w_mty_at_half(i, k, j - 1)) - (w_mty_at_half(i, k, j + 2) - w_mty_at_half(i, k, j - 2)) ) / 12.0 - dw_dz_dz_dy(i, k, j))
            else
              dw_dy(i, k, j) = msf_ty_ty(i, j) * (rdy * 0.5 * (w_mty_at_half(i, k, j + 1) - w_mty_at_half(i, k, j - 1)) - dw_dz_dz_dy(i, k, j))
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!          dw_dy(i, kts, j) = 0.0
          dw_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dy_at_mass


    subroutine Calc_surface_thetav (tsk, t2, q2, psfc, thetav_sk, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, jms:jme),          intent (in)  :: tsk, t2, q2, psfc
      real, dimension (its:ite, jts:jte),          intent (out) :: thetav_sk
      integer,                                     intent (in)  :: ide, jde, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
        ! Local vars
      integer :: i, j
      integer :: i_start, i_end, j_start, j_end
      real    :: sqv, theta_sk

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      do j = j_start, j_end
        do i = i_start, i_end
          sqv = q2(i, j) / (1.0 + q2(i, j))
          theta_sk = t2(i, j) * (100000. / psfc(i, j)) ** rcp
          thetav_sk(i, j) = theta_sk * (1.0 + P608 * sqv)
        end do
      end do

    end subroutine Calc_surface_thetav


    subroutine Calc_mass_var_derivatives_at_mass (thetav, thetav_sk, q2, qv, dthetav_dx, dthetav_dy, dthetav_dz, &
              dqv_dx, dqv_dy, dqv_dz, rdzw, zx, zy, msftx, msfty, dn, dnw, fnm, fnp, rdx, rdy, cf1, cf2, cf3, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                     !
      ! Purpose: Calculate spatial derivatives of mass variables at mass points             !
      !                                                                                     !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano                                      !
      !                                                                                     !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: thetav, qv
      real, dimension (its:ite, jts:jte),          intent (in)  :: thetav_sk
      real, dimension (ims:ime, jms:jme),          intent (in)  :: q2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: rdzw, zx, zy
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfty
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dthetav_dx, dthetav_dy, dthetav_dz, &
                                                                   dqv_dx, dqv_dy, dqv_dz
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3

      logical, parameter :: DEBUG = .false.

        ! Calculate thetav derivatives
      call Calc_dmassvar_dx_at_mass (thetav, zx, rdzw, msftx, dn, dnw, fnm, fnp, rdx, cf1, cf2, cf3, ide, jde, kde, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dx)

      call Calc_dmassvar_dy_at_mass (thetav, zy, rdzw, msfty, dn, dnw, fnm, fnp, rdy, cf1, cf2, cf3, ide, jde, kde, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dy)

      call Calc_dmassvar_dz_at_mass (var=thetav, thetav_sk=thetav_sk, rdzw=rdzw, dn=dn, dnw=dnw, fnm=fnm, fnp=fnp, &
          cf1=cf1, cf2=cf2, cf3=cf3, ide=ide, jde=jde, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, dvar_dz=dthetav_dz)

        ! Calculate qv derivatives
      call Calc_dmassvar_dx_at_mass (qv, zx, rdzw, msftx, dn, dnw, fnm, fnp, rdx, cf1, cf2, cf3, ide, jde, kde, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dx)

      call Calc_dmassvar_dy_at_mass (qv, zy, rdzw, msfty, dn, dnw, fnm, fnp, rdy, cf1, cf2, cf3, ide, jde, kde, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dy)

      call Calc_dmassvar_dz_at_mass (var=qv, q2=q2, rdzw=rdzw, dn=dn, dnw=dnw, fnm=fnm, fnp=fnp, &
          cf1=cf1, cf2=cf2, cf3=cf3, ide=ide, jde=jde, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, dvar_dz=dqv_dz)

      if (DEBUG) print *, 'Done Calc_mass_var_derivatives_at_mass...'

    end subroutine Calc_mass_var_derivatives_at_mass


    subroutine Calc_dmassvar_dx_at_mass (var, zx, rdzw, msftx, dn, dnw, fnm, fnp, rdx, cf1, cf2, cf3,  &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dx)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial x        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: var, zx, rdzw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx
      real, dimension (kms:kme),                   intent (in)  :: fnm, fnp, dn, dnw
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dvar_dx

        ! Local vars
      integer :: i, j, k
      real :: dz_dx
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: var_at_u, var_at_full, dvar_dz_dz_dx
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktf = Min (kte, kde - 1)

      ktes1 = kte - 1
      ktes2 = kte - 2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ! Interpolate mass var to u points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            var_at_u(i, k, j) = 0.5 * (var(i - 1, k, j) + var(i, k, j))
          end do
        end do
      end do

        ! Interpolate mass var from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            var_at_full(i, k, j) = 0.5 * (fnm(k) * (var(i - 1, k, j) + var(i, k, j)) + &
                fnp(k) * (var(i - 1, k - 1, j) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
          var_at_full(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
              cf3 * var(i, 3, j) + cf1 * var(i - 1, 1, j) +  &
              cf2 * var(i - 1, 2, j) + cf3 * var(i - 1, 3, j))
          var_at_full(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
              var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              var(i - 1, ktes1, j) + (var(i - 1, ktes1, j) - &
              var(i - 1, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dvar_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at mass points
            dz_dx = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
              ! dvar_dz_dz_dx = (partial var / partial z)|sigma (dvar_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dvar_dz_dz_dx(i, k, j) =  (0.5 * (var_at_full(i, k + 1, j) + var_at_full(i + 1, k + 1, j)) - &
                0.5 * (var_at_full(i, k, j) + var_at_full(i + 1, k, j))) * rdzw(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial var / partial x)|z = msftx * (partial var_at_u / partial x)|sigma - dvar_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dvar_dx(i, k, j) = msftx(i, j) * (var_at_u(i + 1, k, j) -  var_at_u(i, k, j)) * rdx - dvar_dz_dz_dx(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!          dvar_dx(i, 1, j) = 0.0
          dvar_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dx_at_mass


    subroutine Calc_dmassvar_dy_at_mass (var, zy, rdzw, msfty, dn, dnw, fnm, fnp, rdy, cf1, cf2, cf3,  &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dy)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial y        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: var, zy, rdzw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty
      real, dimension (kms:kme),                   intent (in)  :: fnm, fnp, dn, dnw
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dvar_dy

        ! Local vars
      integer :: i, j, k
      real :: dz_dy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: var_at_v, var_at_full, dvar_dz_dz_dy
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktf = Min (kte, kde - 1)

      ktes1 = kte - 1
      ktes2 = kte - 2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Interpolate mass var from to v points
      do j = j_start, j_end + 1
        do k = kts,   ktf
          do i = i_start, i_end
            var_at_v(i, k, j) = 0.5 * (var(i, k, j - 1) + var(i, k, j))
          end do
        end do
      end do

        ! Interpolate mass var from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1,   ktf
          do i = i_start, i_end
            var_at_full(i, k, j) = 0.5 * (fnm(k) * (var(i, k, j - 1) + var(i, k, j)) + &
                fnp(k) * (var(i, k - 1, j - 1) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
          var_at_full(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
              cf3 * var(i, 3, j) + cf1 * var(i, 1, j - 1) +  &
              cf2 * var(i, 2, j - 1) + cf3 * var(i, 3, j - 1))
          var_at_full(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
              var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              var(i, ktes1, j - 1) + (var(i, ktes1, j - 1) - &
              var(i, ktes2, j - 1)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dvar_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at mass points
            dz_dy = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
              ! dvar_dz_dz_dy = (partial var / partial z)|sigma (dvar_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dvar_dz_dz_dy(i, k, j) = (0.5 * (var_at_full(i, k + 1, j) + var_at_full(i, k + 1, j + 1)) - &
                0.5 * (var_at_full(i, k, j) + var_at_full(i, k, j + 1))) * rdzw(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial var / partial y)|z = msfty * (partial var_at_full / partial y)|sigma - dvar_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dvar_dy(i, k, j) = msfty(i, j) * (var_at_v(i, k, j + 1) - var_at_v(i, k, j)) * rdy - dvar_dz_dz_dy(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!          dvar_dy(i, 1, j) = 0.0
          dvar_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dy_at_mass


    subroutine Calc_dmassvar_dz_at_mass (var, thetav_sk, q2, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial z        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme),  intent (in)  :: var, rdzw
      real, dimension (its:ite, jts:jte), optional, intent (in)  :: thetav_sk
      real, dimension (ims:ime, jms:jme), optional, intent (in)  :: q2
      real, dimension (kms:kme),                    intent (in)  :: dn, dnw, fnm, fnp
      integer,                                      intent (in)  :: ide, jde, kde,                  &
                                                                    ims, ime, jms, jme, kms, kme,   &
                                                                    its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte),  intent (out) :: dvar_dz
      real,                                         intent (in)  :: cf1, cf2, cf3

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end
      real, dimension (its:ite, kts:kte, jts:jte) :: var_at_full

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            var_at_full(i, k, j) = fnm(k) * var(i, k, j) + fnp(k) * var(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          if (PRESENT(thetav_sk)) then
            var_at_full(i, 1, j) =  thetav_sk(i, j) - 300.
          end if
          if (PRESENT(q2)) then
            var_at_full(i, 1, j) =  q2(i, j)
          end if
            ! Top face
          var_at_full(i, kte, j) =  cft1 * var(i, ktes1, j) + cft2 * var(i, ktes2, j)
        end do
      end do

        ! Calculate derivative at mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dvar_dz(i, k, j) = (var_at_full(i, k + 1, j) - var_at_full(i, k, j)) * rdzw(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          dvar_dz(i, kts, j) = 0.0
            ! Top
          dvar_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dz_at_mass


    subroutine Diagnose_fluxes_l2p5_pbl_approx (u, v, t, qv, rho, ust, hfx, du_dz, dv_dz, dthetav_dz, dqv_dz, &
              q_sq, q_sq_hl88, q_ratio, sm_tilde, sh_tilde, sm_l2, sh_l2, rif, l_master, &
              u2, v2, w2, uv, uw, vw, wth, wth_face, wthv, wthv_face, wqv, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)   :: u, v, t, qv, rho
      real, dimension (ims:ime, jms:jme),          intent (in)   :: ust, hfx
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)   :: du_dz, dv_dz, dthetav_dz, dqv_dz

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)   :: q_sq
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)   :: sm_tilde, sh_tilde, sm_l2, sh_l2, rif, q_sq_hl88, q_ratio
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)   :: l_master, wth_face, wthv_face
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)  :: u2, v2, w2, &
                                                                    uv, uw, vw, wth, wthv
      real, dimension (its:ite, kts:kte, jts:jte), intent (out)  :: wqv
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real :: qsq, l, sm, sh, ri, dudz, dvdz, dthvdz, dqvdz
      real :: c_p, u1, v1, wind_speed, ust_sq
      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      j_loop: do j = j_start, j_end
        k_loop: do k = kts, ktf
          i_loop: do i = i_start, i_end

                  ! Save temporary scalars
                l = l_master(i, k, j)
                if ( q_ratio(i, k, j) <= 1.0 ) then  ! HL88 realizability
                  qsq = q_sq_hl88(i, k, j)
                  sm = sm_l2(i, k, j)
                  sh = sh_l2(i, k, j)
                else
                  qsq = q_sq(i, k, j)
                  sm = sm_tilde(i, k, j)
                  sh = sh_tilde(i, k, j)
                end if
                qsq = Max (qsq, Q_SQ_MIN)
                l = l_master(i, k, j)
                ri = rif(i, k, j)
                dudz = du_dz(i, k, j)
                dvdz = dv_dz(i, k, j)
                dthvdz = dthetav_dz(i, k, j)
                dqvdz = dqv_dz(i, k, j)

                  ! Eq. 27 MY82
                uw(i, k, j) = Diag_uw_pbl_approx (qsq, l, sm, dudz)
                vw(i, k, j) = Diag_vw_pbl_approx (qsq, l, sm, dvdz)
                uv(i, k, j) = Diag_uv_pbl_approx (qsq, l, dudz, dvdz, uw(i, k, j), vw(i, k, j))

                  ! Eq. 28 MY82
                wqv(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqvdz)

                if (k == 1) then
                  wthv(i, k, j) = ( wthv_face(i, k, j) + wthv_face(i, k + 1, j) ) / 2.
                  wth(i, k, j) = ( wth_face(i, k, j) + wth_face(i, k + 1, j) ) / 2.
                else
                  wthv(i, k, j) = Diag_wthv_pbl_approx (qsq, l, sh, dthvdz)
                  wth(i, k, j) = Diag_xth_pbl_approx (t(i, k, j), qv(i, k, j), wthv(i, k, j), wqv(i, k, j))
                end if

                  ! Eq. 26 MY82
                u2(i, k, j) = Diag_u2_pbl_approx (qsq, l, ri, dudz, dvdz, uw(i, k, j), vw(i, k, j), wthv(i, k, j))
                v2(i, k, j) = Diag_v2_pbl_approx (qsq, l, ri, dudz, dvdz, uw(i, k, j), vw(i, k, j), wthv(i, k, j))
                w2(i, k, j) = Diag_w2_pbl_approx (qsq, l, ri, dudz, dvdz, uw(i, k, j), vw(i, k, j), wthv(i, k, j))

                  ! Apply HL88 mod
                if ( q_ratio(i, k, j) <= 1.0) then
                  u2(i, k, j) = u2(i, k, j) * q_ratio(i, k, j) ** 2.0
                  v2(i, k, j) = v2(i, k, j) * q_ratio(i, k, j) ** 2.0
                  w2(i, k, j) = w2(i, k, j) * q_ratio(i, k, j) ** 2.0
                  uv(i, k, j) = uv(i, k, j) * q_ratio(i, k, j) ** 2.0
                  uw(i, k, j) = uw(i, k, j) * q_ratio(i, k, j) ** 2.0
                  vw(i, k, j) = vw(i, k, j) * q_ratio(i, k, j) ** 2.0

                  if (k > 1) then
                    wthv(i, k, j) = wthv(i, k, j) * q_ratio(i, k, j) ** 2.0
                    wth(i, k, j) = wth(i, k, j) * q_ratio(i, k, j) ** 2.0
                  end if

                  wqv(i, k, j) = wqv(i, k, j) * q_ratio(i, k, j) ** 2.0
                end if


                if (DEBUG) call Debug_turb_flux

          end do i_loop
        end do k_loop
      end do j_loop

      contains

        subroutine Debug_turb_flux

          if (i == I_TO_TEST .and. j == J_TO_TEST) then
            print *, '------------------------------------------'
            print *, 'i =', i, 'j = ', j, 'k = ', k
            print *, '     q_sq = ', qsq
            print *, '    dt_dz = ', dthvdz
            print *, '    du_dz = ', dudz
            print *, '    dv_dz = ', dvdz
            print *, ' l_master = ', l
            print *, ' turb_flux_uw = ', uw(i, k, j)
            print *, ' turb_flux_vw = ', vw(i, k, j)
            print *, ' turb_flux_wtheta_v = ', wthv(i, k, j)
          end if

        end subroutine Debug_turb_flux


        subroutine Set_turb_flux_min

          u2(i, k, j) = TURB_FLUX_MIN
          v2(i, k, j) = TURB_FLUX_MIN
          w2(i, k, j) = TURB_FLUX_MIN
          uv(i, k, j) = TURB_FLUX_MIN
          uw(i, k, j) = TURB_FLUX_MIN
          vw(i, k, j) = TURB_FLUX_MIN
          wthv(i, k, j) = TURB_FLUX_MIN

        end subroutine Set_turb_flux_min

    end subroutine Diagnose_fluxes_l2p5_pbl_approx


    subroutine Diagnose_fluxes_l2p5 (du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dthetav_dx, dthetav_dy, dthetav_dz, &
        dqv_dx, dqv_dy, dqv_dz, u2, v2, w2, uv, uw, vw, wth, wthv, wthv_face, tf_wqv, q_sq, q_sq_hl88, q_ratio, l_master, &
        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (its:ite, kts:kte, jts:jte), intent (in)   :: du_dx, du_dy, du_dz, &
                                                                    dv_dx, dv_dy, dv_dz, &
                                                                    dw_dx, dw_dy, dw_dz
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)   :: dthetav_dx, dthetav_dy, dthetav_dz
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)   :: dqv_dx, dqv_dy, dqv_dz

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)   :: q_sq, l_master, wthv_face
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)   :: q_sq_hl88, q_ratio
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)  :: u2, v2, w2, &
                                                                    uv, uw, vw, wth, wthv
      real, dimension (its:ite, kts:kte, jts:jte), intent (out)  :: tf_wqv
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real :: qsq
      real (kind = SOLVER_PREC) :: utv, vtv, wtv, uqv, vqv, wqv
      real, dimension (its:ite, kts:kte, jts:jte) :: uthv, vthv, th2v, tf_uqv, tf_vqv
      real, dimension (its:ite, kts:kte, jts:jte) :: mat_cond_heat, mat_cond_moist
      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      j_loop: do j = j_start, j_end
        k_loop: do k = kts, ktf
          i_loop: do i = i_start, i_end

            if_no_turb: if (q_sq(i, k, j) == Q_SQ_MIN) then

              call Set_turb_flux_min

              if (k == 1) then
                wthv(i, k, j) = ( wthv_face(i, k, j) + wthv_face(i, k + 1, j) ) / 2.
              end if

            else

                ! Save temporary scalars
              if ( q_ratio(i, k, j) <= 1.0 ) then  ! HL88 mod
                qsq = q_sq_hl88(i, k, j)
              else
                qsq = q_sq(i, k, j)
              end if

              call Diagnose_fluxes (0, qsq, l_master(i, k, j), du_dx(i, k, j), du_dy(i, k, j), du_dz(i, k, j),                       &
                  dv_dx(i, k, j), dv_dy(i, k, j), dv_dz(i, k, j), dw_dx(i, k, j), dw_dy(i, k, j), dw_dz(i, k, j),                    &
                  dthetav_dx(i, k, j), dthetav_dy(i, k, j), dthetav_dz(i, k, j), dqv_dx(i, k, j), dqv_dy(i, k, j), dqv_dz(i, k, j),  &
                  u2(i, k, j), v2(i, k, j), w2(i, k, j), uv(i, k, j), uw(i, k, j), vw(i, k, j), uthv(i, k, j), vthv(i, k, j),        &
                  wthv(i, k, j), th2v(i, k, j), tf_uqv(i, k, j), tf_vqv(i, k, j), tf_wqv(i, k, j), utv, vtv, wtv, uqv, vqv, wqv,     &
                  mat_cond_heat(i, k, j), mat_cond_moist(i, k, j)                                                                    )

                ! Apply HL88 mod for level 2.5 model
              if ( q_ratio(i, k, j) <= 1.0 ) then
                u2(i, k, j) = u2(i, k, j) * q_ratio(i, k, j) ** 2.0
                v2(i, k, j) = v2(i, k, j) * q_ratio(i, k, j) ** 2.0
                w2(i, k, j) = w2(i, k, j) * q_ratio(i, k, j) ** 2.0
                uv(i, k, j) = uv(i, k, j) * q_ratio(i, k, j) ** 2.0
                uw(i, k, j) = uw(i, k, j) * q_ratio(i, k, j) ** 2.0
                vw(i, k, j) = vw(i, k, j) * q_ratio(i, k, j) ** 2.0
                wthv(i, k, j) = wthv(i, k, j) * q_ratio(i, k, j) ** 2.0
                tf_wqv(i, k, j) = tf_wqv(i, k, j) * q_ratio(i, k, j) ** 2.0
              end if

              if (k == 1) then
                wthv(i, k, j) = ( wthv_face(i, k, j) + wthv_face(i, k + 1, j) ) / 2.
              end if

            end if if_no_turb

          end do i_loop
        end do k_loop
      end do j_loop

      contains

        subroutine Set_turb_flux_min

          u2(i, k, j) = TURB_FLUX_MIN
          v2(i, k, j) = TURB_FLUX_MIN
          w2(i, k, j) = TURB_FLUX_MIN
          uv(i, k, j) = TURB_FLUX_MIN
          uw(i, k, j) = TURB_FLUX_MIN
          vw(i, k, j) = TURB_FLUX_MIN
          wthv(i, k, j) = TURB_FLUX_MIN
          tf_wqv(i, k, j) = TURB_FLUX_MIN

        end subroutine Set_turb_flux_min

    end subroutine Diagnose_fluxes_l2p5


    subroutine Calc_q_sq_rhs (config_flags, sm, sh, zx, zy, du_dz, dv_dz, dw_dz, du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, &
              dthetav_dz, dqv_dz, q_sq, q_sq_prog, q_sq_tend, q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, &
              q_sq_buoyancy_tend, q_sq_dissip_tend, l_master_at_face, l_master_at_mass, turb_flux_u2, turb_flux_v2, turb_flux_w2, &
              turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta, turb_flux_wtheta_v, dn, dnw, fnm, fnp, &
              msftx, msfty, msfux, msfuy, msfvx, msfvy, mut, c1, c2, z_at_mass, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, dt, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                        !
      ! Purpose: Calculates RHS of prognostic TKE equation (level 2.5)                         !
      !                                                                                        !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                                          !
      !                                                                                        !
      ! Method: 1) Follows module_diffusion_em to calc shear, buoyancy, and dissipation terms  !
      !                                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: zx, zy
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: sm, sh
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dw_dz, &
                            du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, dthetav_dz, dqv_dz
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_vdiff_tend,    &
                                                                     q_sq_hdiff_tend,    &
                                                                     q_sq_shear_tend,    &
                                                                     q_sq_buoyancy_tend, &
                                                                     q_sq_dissip_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master_at_face
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                                                                     turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                                                                     turb_flux_wtheta, turb_flux_wtheta_v
      real, dimension (kms:kme),                   intent (in)    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)    :: msftx, msfty, msfux, msfuy, msfvx, msfvy
      real, dimension (ims:ime, jms:jme),          intent (in)    :: mut
      real, dimension (kms:kme),                   intent (in)    :: c1, c2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: z_at_mass, rdz, rdzw
      real,                                        intent (in)    :: rdx, rdy, cf1, cf2, cf3, dt
      integer, intent (in) :: ids, ide, jds, jde, kds, kde,  &
                              ims, ime, jms, jme, kms, kme,  &
                              its, ite, jts, jte, kts, kte

        ! Local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      logical, parameter :: DEBUG = .false.
      real (kind = kind(q_sq_tend)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

      ktf     = MIN( kte, kde-1 )
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      IF ( config_flags%open_xs .or. config_flags%specified .or. &
           config_flags%nested) i_start = MAX(ids+1,its)
      IF ( config_flags%open_xe .or. config_flags%specified .or. &
           config_flags%nested) i_end   = MIN(ide-2,ite)
      IF ( config_flags%open_ys .or. config_flags%specified .or. &
           config_flags%nested) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. config_flags%specified .or. &
           config_flags%nested) j_end   = MIN(jde-2,jte)

        !!! Diffusion !!!
      call Calc_q_sq_vertical_diffusion (N_PBL3D, sm, q_sq_vdiff_tend, q_sq_tend, l_master_at_face, q_sq_prog, q_sq, &
          fnm, fnp, dnw, rdz, rdzw, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      if ( config_flags%pbl3d_opt .gt. 1 ) then
        call Calc_q_sq_horizontal_diffusion (N_PBL3D, q_sq_hdiff_tend, q_sq_tend, l_master_at_face, q_sq_prog, q_sq, &
            msftx, msfty, msfux, msfvy, fnm, fnp, dn, dnw, zx, zy, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, &
            ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)
      end if

        !!! Shear !!!
      if ( config_flags%pbl3d_opt .lt. 2 ) then
        call Calc_q_sq_shear_pbl_approx (N_PBL3D, sm, l_master_at_mass, q_sq_prog, q_sq_shear_tend, q_sq_tend, du_dz, dv_dz, turb_flux_uw, turb_flux_vw, &
            mut, c1, c2, fnm, fnp, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)
      else
        call Calc_q_sq_shear (N_PBL3D, q_sq_shear_tend, q_sq_tend, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
            turb_flux_u2, turb_flux_uw, turb_flux_v2, turb_flux_vw, turb_flux_w2, turb_flux_uv, mut, c1, c2, fnm, fnp, &
            ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)
      end if

        !!! Buoyancy !!!
      call Calc_q_sq_buoyancy (N_PBL3D, q_sq_buoyancy_tend, q_sq_tend, turb_flux_wtheta_v, dthetav_dz, q_sq_prog, l_master_at_mass, sh, mut, c1, c2, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

        !!! Dissipation !!!
      call Calc_q_sq_dissip (N_PBL3D, q_sq_dissip_tend, q_sq_tend, l_master_at_mass, q_sq_prog, mut, c1, c2, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

        !!! Bound tendency !!!
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_tend(i, k, j) = max ( q_sq_tend(i, k, j), -(c1(k) * mut(i,j) + c2(k)) * max ( 0.0 , N_PBL3D * q_sq_prog(i, k, j) ) / dt )

            if (DEBUG) call Debug_q_sq_tend
          end do
        end do
      end do

      contains

      subroutine Debug_q_sq_tend

        if (i == I_TO_TEST .and. j == J_TO_TEST ) then
          print *, '------------------------------------------'
          print *, '  q_sq_tend final = ', q_sq_tend(i, k, j)
        end if

      end subroutine Debug_q_sq_tend

    end subroutine Calc_q_sq_rhs


    subroutine Calc_q_sq_vertical_diffusion (N_PBL3D, sm, q_sq_vdiff_tend, q_sq_tend, l_master, &
              q_sq_mass, q_sq_face, fnm, fnp, dnw, rdz, rdzw, ims, ime, jms, jme, kms, kme,     &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf             )

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                !
        ! Purpose: Calculate vertical diffusion term for prognostic TKE  !
        !          Follows MY82 Eq. 24                                   !
        !                                                                !
        ! Authors: Timothy W. Juliano & Pedro A. Jimenez                 !
        !                                                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: sm
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_vdiff_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq_mass, q_sq_face
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master
      real, dimension (kms:kme),                   intent (in)  :: dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: rdz, rdzw
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      real, parameter :: Sq = 3.0
      real :: q_sq_vdiff_tend_tmp
      real, dimension (its:ite, kts:kte, jts:jte) :: q_sq_vdiff_hold
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            q_sq_vdiff_hold(i, k, j) = - ( Sq * sm(i, k, j) * l_master(i, k, j) * Sqrt(q_sq_face(i, k, j)) * &
                                         ( q_sq_mass(i, k, j) - q_sq_mass(i, k - 1, j) ) * rdz(i, k, j) )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  l_master = ', l_master(i, k, j)
              print *, '  q_sq_face = ', q_sq_face(i, k, j)
              print *, '  q_sq_mass = ', q_sq_mass(i, k, j), q_sq_mass(i, k - 1, j)
              print *, '  q_sq_vdiff_hold = ', q_sq_vdiff_hold(i, k, j)
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          q_sq_vdiff_hold(i, kts, j) = 0.0
          q_sq_vdiff_hold(i, ktf + 1, j) = 0.0
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_vdiff_tend(i, k, j) = rdzw(i, k, j) * ( q_sq_vdiff_hold(i, k + 1, j) - q_sq_vdiff_hold(i, k, j) )
            q_sq_vdiff_tend_tmp = ( G / dnw(k) ) *  &
                                  ( q_sq_vdiff_hold(i, k + 1, j) - q_sq_vdiff_hold(i, k, j) )
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( q_sq_vdiff_tend_tmp / N_PBL3D )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '  q_sq_vdiff_tend = ', q_sq_vdiff_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if 
          end do
        end do
      end do

    end subroutine Calc_q_sq_vertical_diffusion


    subroutine Calc_q_sq_horizontal_diffusion (N_PBL3D, q_sq_hdiff_tend, q_sq_tend, l_master, q_sq_mass, q_sq_face, &
              msftx, msfty, msfux, msfvy, fnm, fnp, dn, dnw, zx, zy, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_hdiff_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq_mass, q_sq_face
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfty, msfux, msfvy
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: zx, zy, rdz, rdzw
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k, ktes1, ktes2
      real :: cft1, cft2, tmpzx, rdzu, rdzv, dz_dx, dz_dy, &
              dq_sq_dz_dz_dx, dq_sq_dz_dz_dy, mrdx, mrdy
      real :: q_sq_hdiff_tend_tmp
      real, parameter :: Sq = 0.20
      real, dimension (its:ite, kts:kte, jts:jte) :: q_sq_hdiff_hold
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) :: xkxavg, q_sq_at_full_u, q_sq_at_full_v, &
           dq_sq_dx_at_u, dq_sq_dy_at_v, dq_sq_dx_at_w, dq_sq_dy_at_w, zx_at_mass, zy_at_mass

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

        !!! Do dq_sq/dx first !!!

         ! Interpolate q_sq from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            q_sq_at_full_u(i, k, j) = 0.5 * (fnm(k) * (q_sq_mass(i - 1, k, j) + q_sq_mass(i, k, j)) + &
                fnp(k) * (q_sq_mass(i - 1, k - 1, j) + q_sq_mass(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
          q_sq_at_full_u(i, kts, j) = 0.5 * (cf1 * q_sq_mass(i, 1, j) + cf2 * q_sq_mass(i, 2, j) + &
              cf3 * q_sq_mass(i, 3, j) + cf1 * q_sq_mass(i - 1, 1, j) +  &
              cf2 * q_sq_mass(i - 1, 2, j) + cf3 * q_sq_mass(i - 1, 3, j))
          q_sq_at_full_u(i, ktf + 1, j) = 0.5 * (q_sq_mass(i, ktes1, j) + (q_sq_mass(i, ktes1, j) - &
              q_sq_mass(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              q_sq_mass(i - 1, ktes1, j) + (q_sq_mass(i - 1, ktes1, j) - &
              q_sq_mass(i - 1, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate average eddy diffusivity at u points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            xkxavg(i, k, j) = Sq * 0.25 * ( l_master(i - 1, k, j) * Sqrt(q_sq_face(i - 1, k, j)) +  &
                                            l_master(i, k, j) * Sqrt(q_sq_face(i, k, j))         +  &
                                            l_master(i - 1, k + 1, j) * Sqrt(q_sq_face(i - 1, k + 1, j)) +  &
                                            l_master(i, k + 1, j) * Sqrt(q_sq_face(i, k + 1, j))         )
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
          xkxavg(i, ktf + 1, j) = xkxavg(i, ktf, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dq_sq_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
              ! dz_dx = (partial z / partial x)|sigma at u points
            dz_dx = 0.5 * (zx(i, k, j) + zx(i, k + 1, j))
            rdzu = 2. / (1. / rdzw(i, k, j) + 1. / rdzw(i - 1, k, j))
              ! dq_sq_dz_dz_dx = (partial q_sq / partial z)|sigma (dq_sq_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dq_sq_dz_dz_dx = (q_sq_at_full_u(i, k + 1, j) - q_sq_at_full_u(i, k, j)) * rdzu * dz_dx
              ! (partial q_sq / partial x)|z = msfux * (partial q_sq_at_mass / partial x)|sigma - dq_sq_dz_dz_dx)
              ! at u points
            dq_sq_dx_at_u(i, k, j) = - msfux(i, j) * xkxavg(i, k, j) *  &
                                       0.5 * ( q_sq_at_full_u(i, k + 1, j) + q_sq_at_full_u(i, k, j) ) * &
                                       (q_sq_mass(i, k, j) - q_sq_mass(i - 1, k, j)) * rdx - dq_sq_dz_dz_dx
          end do
        end do
      end do

        !!! Do dq_sq/dy second !!!

         ! Interpolate q_sq from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            q_sq_at_full_v(i, k, j) = 0.5 * (fnm(k) * (q_sq_mass(i, k, j - 1) + q_sq_mass(i, k, j)) + &
                fnp(k) * (q_sq_mass(i, k - 1, j - 1) + q_sq_mass(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
          q_sq_at_full_v(i, kts, j) = 0.5 * (cf1 * q_sq_mass(i, 1, j) + cf2 * q_sq_mass(i, 2, j) + &
              cf3 * q_sq_mass(i, 3, j) + cf1 * q_sq_mass(i, 1, j - 1) +  &
              cf2 * q_sq_mass(i, 2, j - 1) + cf3 * q_sq_mass(i, 3, j - 1))
          q_sq_at_full_v(i, ktf + 1, j) = 0.5 * (q_sq_mass(i, ktes1, j) + (q_sq_mass(i, ktes1, j) - &
              q_sq_mass(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              q_sq_mass(i, ktes1, j - 1) + (q_sq_mass(i, ktes1, j - 1) - &
              q_sq_mass(i, ktes2, j - 1)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate average eddy diffusivity at v points
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            xkxavg(i, k, j) = Sq * 0.25 * ( l_master(i, k, j - 1) * Sqrt(q_sq_face(i, k, j - 1)) +  &
                                            l_master(i, k, j) * Sqrt(q_sq_face(i, k, j))         +  &
                                            l_master(i, k + 1, j - 1) * Sqrt(q_sq_face(i, k + 1, j - 1)) +  &
                                            l_master(i, k + 1, j) * Sqrt(q_sq_face(i, k + 1, j))         )
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
          xkxavg(i, ktf + 1, j) = xkxavg(i, ktf, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dq_sq_dz_dz_dy)
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial y)|sigma at y points
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k + 1, j))
            rdzv = 2. / (1. / rdzw(i, k, j) + 1. / rdzw(i, k, j - 1))
              ! dq_sq_dz_dz_dy = (partial q_sq / partial z)|sigma (dq_sq_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dq_sq_dz_dz_dy = (q_sq_at_full_v(i, k + 1, j) - q_sq_at_full_v(i, k, j)) * rdzv * dz_dy
              ! (partial q_sq / partial y)|z = msfvy * (partial q_sq_at_mass / partial y)|sigma - dq_sq_dz_dz_dy)
              ! at v points
            dq_sq_dy_at_v(i, k, j) = - msfvy(i, j) * xkxavg(i, k, j) * &
                                       0.5 * ( q_sq_at_full_v(i, k + 1, j) + q_sq_at_full_v(i, k, j) ) * &
                                       (q_sq_mass(i, k, j) - q_sq_mass(i, k, j - 1)) * rdy - dq_sq_dz_dz_dy
          end do
        end do
      end do

        !!! Do full horizontal diffusion !!!

        ! Interpolate dq_sq_dx and dq_sq_dy to w points
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dq_sq_dx_at_w(i, k, j) = 0.5 * (fnm(k) * (dq_sq_dx_at_u(i + 1, k, j) + dq_sq_dx_at_u(i, k, j)) +  &
                                            fnp(k) * (dq_sq_dx_at_u(i + 1, k - 1, j) + dq_sq_dx_at_u(i, k - 1, j)))
            dq_sq_dy_at_w(i, k, j) = 0.5 * (fnm(k) * (dq_sq_dy_at_v(i, k, j + 1) + dq_sq_dy_at_v(i, k, j)) +  &
                                            fnp(k) * (dq_sq_dy_at_v(i, k - 1, j + 1) + dq_sq_dy_at_v(i, k - 1, j)))

            zx_at_mass(i, k, j) = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
            zy_at_mass(i, k, j) = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dq_sq_dx_at_w(i, kts, j) = 0.
          dq_sq_dx_at_w(i, ktf + 1, j) = 0.
          dq_sq_dy_at_w(i, kts, j) = 0.
          dq_sq_dy_at_w(i, ktf + 1,j) = 0.

          zx_at_mass(i, kts, j) = 0.25 * (zx(i, kts, j) + zx(i + 1, kts ,j) + zx(i, kts + 1, j) + zx(i + 1, kts + 1, j))
          zy_at_mass(i, kts, j) = 0.25 * (zy(i, kts, j) + zy(i, kts, j + 1) + zy(i, kts + 1, j) + zy(i, kts + 1, j + 1))
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msftx(i, j) * rdx
            mrdy = msfty(i, j) * rdy

            q_sq_hdiff_tend(i, k, j) = 2. * ( mrdx * (dq_sq_dx_at_u(i + 1, k, j) - dq_sq_dx_at_u(i, k, j)) + &
                                              mrdy * (dq_sq_dy_at_v(i, k, j + 1) - dq_sq_dy_at_v(i, k, j)) - &
                                              msftx(i, j) * zx_at_mass(i, k, j) * (dq_sq_dx_at_w(i, k + 1, j) - dq_sq_dx_at_w(i, k, j)) * rdzw(i, k, j) - &
                                              msfty(i, j) * zy_at_mass(i, k, j) * (dq_sq_dy_at_w(i, k + 1, j) - dq_sq_dy_at_w(i, k, j)) * rdzw(i, k, j) )
            q_sq_hdiff_tend_tmp = G / (dnw(k) * rdzw(i, k, j)) * &
                                  2. * ( mrdx * (dq_sq_dx_at_u(i + 1, k, j) - dq_sq_dx_at_u(i, k, j)) + &
                                         mrdy * (dq_sq_dy_at_v(i, k, j + 1) - dq_sq_dy_at_v(i, k, j)) - &
                                         msftx(i, j) * zx_at_mass(i, k, j) * (dq_sq_dx_at_w(i, k + 1, j) - dq_sq_dx_at_w(i, k, j)) * rdzw(i, k, j) - &
                                         msfty(i, j) * zy_at_mass(i, k, j) * (dq_sq_dy_at_w(i, k + 1, j) - dq_sq_dy_at_w(i, k, j)) * rdzw(i, k, j) )
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( q_sq_hdiff_tend_tmp / N_PBL3D )
          end do
        end do
      end do

    end subroutine Calc_q_sq_horizontal_diffusion


    subroutine Calc_q_sq_shear_pbl_approx (N_PBL3D, sm, l_master, q_sq, q_sq_shear_tend, q_sq_tend, du_dz, dv_dz, turb_flux_uw, turb_flux_vw, &
              mut, c1, c2, fnm, fnp, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master, q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_shear_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: sm, du_dz, dv_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: turb_flux_uw, &
                                                                   turb_flux_vw

      real, dimension (ims:ime, jms:jme), intent (in) :: mut
      real, dimension (kms:kme),          intent (in) :: c1, c2, fnm, fnp
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_shear_tend(i, k, j) = -2. * ( turb_flux_uw(i, k, j) * du_dz(i, k, j) +  &
                                               turb_flux_vw(i, k, j) * dv_dz(i, k, j) )
!            if (sm(i, k, j) == 0.0) then
!              q_sq_shear_tend(i, k, j) = 0.0
!            else
!              q_sq_shear_tend(i, k, j) = ( c1(k) * mut(i, j) + c2(k) ) *             &
!                                         2. * ( turb_flux_uw(i, k, j) * turb_flux_uw(i, k, j) / (l_master(i, k, j) * q_sq(i, k, j) * sm(i, k, j)) +  &
!                                                turb_flux_vw(i, k, j) * turb_flux_vw(i, k, j) / (l_master(i, k, j) * q_sq(i, k, j) * sm(i, k, j)) )
!            end if
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( (c1(k) * mut(i, j) + c2(k)) * q_sq_shear_tend(i, k, j) / N_PBL3D )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq_shear_tend = ', q_sq_shear_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if
          end do
        end do
      end do

    end subroutine Calc_q_sq_shear_pbl_approx


    subroutine Calc_q_sq_shear (N_PBL3D, q_sq_shear_tend, q_sq_tend, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
              turb_flux_u2, turb_flux_uw, turb_flux_v2, turb_flux_vw, turb_flux_w2, turb_flux_uv, mut, c1, c2, fnm, fnp, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_shear_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: du_dx, du_dy, du_dz,  &
                                                                   dv_dx, dv_dy, dv_dz,  &
                                                                   dw_dx, dw_dy, dw_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: turb_flux_u2,         &
          turb_flux_uw, turb_flux_v2, turb_flux_vw, turb_flux_w2, turb_flux_uv

      real, dimension (ims:ime, jms:jme), intent (in) :: mut
      real, dimension (kms:kme),          intent (in) :: c1, c2, fnm, fnp
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_shear_tend(i, k, j) = -2. * ( turb_flux_u2(i, k, j) * du_dx(i, k, j) +  &
                                               turb_flux_uv(i, k, j) * du_dy(i, k, j) +  &
                                               turb_flux_uw(i, k, j) * du_dz(i, k, j) +  &
                                               turb_flux_uv(i, k, j) * dv_dx(i, k, j) +  &
                                               turb_flux_v2(i, k, j) * dv_dy(i, k, j) +  &
                                               turb_flux_vw(i, k, j) * dv_dz(i, k, j) +  &
                                               turb_flux_uw(i, k, j) * dw_dx(i, k, j) +  &
                                               turb_flux_vw(i, k, j) * dw_dy(i, k, j) +  &
                                               turb_flux_w2(i, k, j) * dw_dz(i, k, j) )
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( (c1(k) * mut(i, j) + c2(k)) * q_sq_shear_tend(i, k, j) / N_PBL3D )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq_shear_tend = ', q_sq_shear_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if
          end do
        end do
      end do

    end subroutine Calc_q_sq_shear


    subroutine Calc_q_sq_buoyancy (N_PBL3D, q_sq_buoyancy_tend, q_sq_tend, turb_flux_wtheta_v, dthv_dz, q_sq, l_master, sh, mut, c1, c2, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_buoyancy_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq, l_master, turb_flux_wtheta_v
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: dthv_dz, sh

      real, dimension (ims:ime, jms:jme), intent (in) :: mut
      real, dimension (kms:kme),          intent (in) :: c1, c2
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_buoyancy_tend(i, k, j) = 2. * G_OVER_TREF * turb_flux_wtheta_v(i, k, j)
!            q_sq_buoyancy_tend(i, k, j) = ( c1(k) * mut(i, j) + c2(k) ) *  &
!                                          -2. * G_OVER_TREF * l_master(i, k, j) * Sqrt(q_sq(i, k, j)) * sh(i, k, j) * dthv_dz(i, k, j)
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( (c1(k) * mut(i, j) + c2(k)) * q_sq_buoyancy_tend(i, k, j) / N_PBL3D )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq_buoyancy_tend = ', q_sq_buoyancy_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if
          end do
        end do
      end do

!      k = kts
!      do j = j_start, j_end
!        do i = i_start, i_end
!          q_sq_buoyancy_tend(i, k, j) = ( c1(k) * mut(i, j) + c2(k) ) *  &
!                                        -2. * G_OVER_TREF * ((l_master(i, k, j) * Sqrt(q_sq(i, k, j)) * sh(i, k, j) * dthv_dz(i, k, j)) - turb_flux_wtheta(i, k, j)) / 2.
!          q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( q_sq_buoyancy_tend(i, k, j) / N_PBL3D )
!
!          if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
!            print *, '------------------------------------------'
!            print *, '  i =', i, 'j = ', j, 'k = ', k
!            print *, '  q_sq_buoyancy_tend = ', q_sq_buoyancy_tend(i, k, j)
!            print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
!          end if
!        end do
!      end do

    end subroutine Calc_q_sq_buoyancy


    subroutine Calc_q_sq_dissip (N_PBL3D, q_sq_dissip_tend, q_sq_tend, l_master, q_sq, mut, c1, c2, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_dissip_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master

      real, dimension (ims:ime, jms:jme), intent (in) :: mut
      real, dimension (kms:kme),          intent (in) :: c1, c2
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            if ( DEBUG ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq = ', q_sq(i, k, j)
              print *, '  l_master = ', l_master(i, k, j)
            end if
            q_sq_dissip_tend(i, k, j) = 2. * ( q_sq(i, k, j) ** 1.5 ) /  &
                                             ( b1 * l_master(i, k, j) )
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) - ( (c1(k) * mut(i, j) + c2(k)) * q_sq_dissip_tend(i, k, j) / N_PBL3D )

            if ( DEBUG ) then
              print *, '------------------------------------------'
!              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq_dissip_tend = ', q_sq_dissip_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if
          end do
        end do
      end do

    end subroutine Calc_q_sq_dissip


    subroutine Taper_q_sq (q_sq, ids, jds, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      integer,                                     intent (in)    :: ids, jds, ide, jde, kde,      &
                                                                     ims, ime, jms, jme, kms, kme, &
                                                                     its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      integer :: bdy_dist, taper_dist
      real    :: taper

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            bdy_dist = min( min((i-ids),(ide-i)) , min((j-jds),(jde-j)) )
            bdy_dist = max( 1, min( bdy_dist, 5 ) )
            taper_dist = 5
            if (bdy_dist .LE. taper_dist) then
              taper = max(0., min( 1., real(bdy_dist) / real(taper_dist) ) )
              q_sq(i, k, j) = q_sq(i, k, j) * taper
            end if
          end do
        end do
      end do

    end subroutine Taper_q_sq


    function Calc_median (list, num)

      implicit none

      real, dimension(num), intent(in) :: list
      integer,              intent(in) :: num
      real                             :: Calc_median

      if (MOD(num,2) .NE. 0) then
        calc_median = REAL(list(num/2 + 1))
      else
        calc_median = REAL(list(num/2) + list(num/2 + 1)) / 2.0
      end if

    end function Calc_median


    subroutine Init_pbl3d (restart, pbl3d_opt, pbl3d_sfc_opt, pbl3d_l_opt, pbl3d_prog, pbl3d_constants, &
              pbl3d_nsteps, fnm, fnp, dnw, ph, phb, t_2, qv, thetav, alt, u, v, ust, hfx, qfx, rmol, &
              q_sq, q_sq_prog, q_sq_vdiff, q_sq_hdiff, q_sq_shear, q_sq_buoyancy, q_sq_dissip,                       &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme,               &
              its, ite, jts, jte, kts, kte, l_master, l_master_at_mass, turb_flux_u2,                     &
              turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw,                   &
              turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v,                     &
              turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_utheta,                 &
              turb_flux_vtheta, turb_flux_wtheta, turb_flux_uqv,                        &
              turb_flux_vqv, turb_flux_wqv                                              )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Initializes the PBL3D                       !
      !                                                      !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano        !
      !                                                      !
      ! Method: 1) User selects set of closure constants     !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      logical, intent(in) :: restart
      integer, intent(in) :: pbl3d_opt, pbl3d_sfc_opt, pbl3d_l_opt, pbl3d_prog, pbl3d_nsteps
      character (len = *), intent(in) :: pbl3d_constants

        ! Grid vars
      real, dimension(kms:kme),                   intent(in)  :: fnm, fnp, dnw
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)  :: ph, phb

        ! Thermo/kinematic state vars
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)  :: t_2, qv, u, v, alt
      real, dimension(ims:ime, jms:jme),          intent(in)  :: ust, hfx, qfx, rmol
      real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: thetav

        ! Turb state vars
      real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: q_sq_prog, q_sq,     &
          q_sq_vdiff, q_sq_hdiff, q_sq_shear, q_sq_buoyancy, q_sq_dissip
      real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: l_master, l_master_at_mass, &
          turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,          &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                  &
          turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v,      &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,            &
          turb_flux_uqv, turb_flux_vqv, turb_flux_wqv

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      logical, parameter :: DEBUG = .true.

      if (DEBUG) then
        write (OUTPUT_UNIT, *) 'Inside subroutine Init_turb_mixing_3d:'
        write (OUTPUT_UNIT, *) 'PBL3D option = ', pbl3d_opt
        write (OUTPUT_UNIT, *) 'Model constants = ', trim(pbl3d_constants)
        write (OUTPUT_UNIT, *) 'PBL3D number of substeps = ', pbl3d_nsteps
      end if

      call Set_model_constants (pbl3d_constants)

      call Set_my_l2_global_vars

      call Set_my_qsq_min (pbl3d_prog)

      call Set_init_turb_state (restart,                                 &
          pbl3d_opt, pbl3d_prog, pbl3d_sfc_opt, pbl3d_l_opt,             &
          ids, ide, jds, jde, kds, kde,                                  &
          ims, ime, jms, jme, kms, kme,                                  &
          its, ite, jts, jte, kts, kte,                                  &
          fnm, fnp, dnw, ph, phb,                                        &
          t_2, qv, thetav, alt, u, v, ust, hfx, qfx, rmol,               &
          q_sq, q_sq_prog, q_sq_vdiff, q_sq_hdiff, q_sq_shear,           &
          q_sq_buoyancy, q_sq_dissip,                                    &
          l_master, l_master_at_mass,                                    &
          turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,        &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                &
          turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v,    &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,          &
          turb_flux_uqv, turb_flux_vqv, turb_flux_wqv                    )

    end subroutine Init_pbl3d


    subroutine Set_model_constants (pbl3d_constants)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Assigns values to the model constants                        !
      !                                                                       !
      ! Author: Pedro A. Jimenez                                              !
      !                                                                       !
      ! Method: 1) M73 constants                                              !
      !                                                                       !
      !         2) MY74 constants                                             !
      !                                                                       !
      !         3) MY82 constants                                             !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      character (len = *), intent(in) :: pbl3d_constants

      logical, parameter :: DEBUG = .false.

      real, parameter :: A1_M73 = 0.78, &
                         A2_M73 = 0.79, &
                         B1_M73 = 15.0, &
                         B2_M73 = 8.0,  &
                         C1_M73 = 0.056,&
                         C2_M73 = 0.0  ,&
                         C3_M73 = 0.0  ,&
                         ALPHA_M73 = 0.1, &
                         LAMBDA_M73 = 0.23

      real, parameter :: A1_MY74 = 0.78, &
                         A2_MY74 = 0.78, &
                         B1_MY74 = 15.0, &
                         B2_MY74 = 8.0,  &
                         C1_MY74 = 0.056,&
                         C2_MY74 = 0.0  ,&
                         C3_MY74 = 0.0  ,&
                         ALPHA_MY74 = 0.1, &
                         LAMBDA_MY74 = 0.23

      real, parameter :: A1_MY82 = 0.92, &
                         A2_MY82 = 0.74, &
                         B1_MY82 = 16.6, &
                         B2_MY82 = 10.1, &
                         C1_MY82 = 0.08 ,&
                         C2_MY82 = 0.0  ,&
                         C3_MY82 = 0.0  ,&
                         ALPHA_MY82 = 0.1, &
                         LAMBDA_MY82 = 0.23

      real, parameter :: A1_MYJ = 0.659888514560862645,    &
                         A2_MYJ = 0.6574209922667784586,   &
                         B1_MYJ = 11.87799326209552761,    &
                         B2_MYJ = 7.226971804046074028,    &
                         C1_MYJ = 0.000830955950095854396, &
                         C2_MYJ = 0.0,                     &
                         C3_MYJ = 0.0,                     &
                         ALPHA_MYJ = 0.3,                  &
                         LAMBDA_MYJ = 0.23

      real, parameter :: A1_Boulac = 0.3,     &
                         A2_Boulac = 0.33,    &
                         B1_Boulac = 8.4,     &
                         B2_Boulac = 6.4,     &
                         C1_Boulac = 0.08,    &
                         C2_Boulac = 0.0,     &
                         C3_Boulac = 0.0,     &
                         ALPHA_Boulac = 0.1,  &
                         LAMBDA_Boulac = 0.23

        ! Set model contants
      select case (trim(pbl3d_constants))

        case('M73')

          a1 = A1_M73
          a2 = A2_M73
          b1 = B1_M73
          b2 = B2_M73
          c1 = C1_M73
          c2 = C2_M73
          c3 = C3_M73
          alpha = ALPHA_M73
          lambda = LAMBDA_M73

        case('MY74')

          a1 = A1_MY74
          a2 = A2_MY74
          b1 = B1_MY74
          b2 = B2_MY74
          c1 = C1_MY74
          c2 = C2_MY74
          c3 = C3_MY74
          alpha = ALPHA_MY74
          lambda = LAMBDA_MY74

        case('MY82')

          a1 = A1_MY82
          a2 = A2_MY82
          b1 = B1_MY82
          b2 = B2_MY82
          c1 = C1_MY82
          c2 = C2_MY82
          c3 = C3_MY82
          alpha = ALPHA_MY82
          lambda = LAMBDA_MY82

        case('MYJ')

          a1 = A1_MYJ
          a2 = A2_MYJ
          b1 = B1_MYJ
          b2 = B2_MYJ
          c1 = C1_MYJ
          c2 = C2_MYJ
          c3 = C3_MYJ
          alpha = ALPHA_MYJ
          lambda = LAMBDA_MYJ

        case('Boulac')

          a1 = A1_Boulac
          a2 = A2_Boulac
          b1 = B1_Boulac
          b2 = B2_Boulac
          c1 = C1_Boulac
          c2 = C2_Boulac
          c3 = C3_Boulac
          alpha = ALPHA_Boulac
          lambda = LAMBDA_Boulac

        case default

          write (OUTPUT_UNIT, *) 'ERROR: Unknown definition of model constants'
          write (OUTPUT_UNIT, *) 'Valid values of turbmix3d_constants are Y73, MY74 and MY82'
          stop

      end select

      if (DEBUG) then
        write (OUTPUT_UNIT, *) 'Inside subroutine Set_model_constants:'
        write (OUTPUT_UNIT, *) 'Model constants:'
        write (OUTPUT_UNIT, *) 'a1 = ', a1
        write (OUTPUT_UNIT, *) 'a2 = ', a2
        write (OUTPUT_UNIT, *) 'b1 = ', b1
        write (OUTPUT_UNIT, *) 'b2 = ', b2
        write (OUTPUT_UNIT, *) 'c1 = ', c1
        write (OUTPUT_UNIT, *) 'c2 = ', c2
        write (OUTPUT_UNIT, *) 'c3 = ', c3
        write (OUTPUT_UNIT, *) 'Alpha = ', alpha
        write (OUTPUT_UNIT, *) 'Lambda = ', lambda
      end if

      if (DEBUG) print *, 'Done Set_model_constants...'

    end subroutine Set_model_constants


    subroutine Set_my_l2_global_vars

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Set global variables for Level 2                             !
      !                                                                       !
      ! Author: Pedro A. Jimenez (adapted from MYNN)                          !
      !                                                                       !
      ! Method: 1) M75 global vars                                            !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      logical, parameter :: DEBUG = .false.

      real              :: e1, e2, e3, e4, e5

        ! For the level 2 model (Y75)
      e1 = b1 - 6.0 * a1
      e2 = b1 + 12.0 * a1 * (1.0 - c2) + 3.0 * b2 * (1.0 - c3)
      e3 = b1 * (1.0 - 3.0 * c1) - 6.0 * a1
      e4 = b1 * (1.0 - 3.0 * c1) + 12.0 * a1 * (1.0 - c2) + 9.0 * a2 * (1.0 - c2)
      e5 = b1 + 3.0 * a1 * (1.0 - c2) + 3.0 * b2 * (1.0 - c3)

      rif_c = e1 / e2
      rf1   = e3 / e4
      rf2   = e1 / e5
      cm    = a1 * e2 * e4 / (b1 * e5)
      ch2   = a2 * e2 / b1

      fac1 = a2 * e5 / (2.0 * a1 * e4)
      fac2 = a1 * e3 / (a2 * e5)
      fac3 = 2.0 * (a1 / a2) * (e3 * e5 - 2.0 * e1 * e4 ) / (e5 * e5)
      fac4 = (a1 * e3 / (a2 * e5) ) ** 2.0

      if (DEBUG) print *, 'Done Set_my_l2_global_vars...'

    end subroutine Set_my_l2_global_vars


    subroutine Set_my_qsq_min (pbl3d_prog)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Set minimum values for Q_SQ and TURB_FLUX                    !
      !                                                                       !
      ! Author: Timothy W. Juliano                                            !
      !                                                                       !
      ! Method: Depending on level 2 or 2.5 model                             !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: pbl3d_prog

      logical, parameter :: DEBUG = .false.

      if (pbl3d_prog == 0) then
        Q_SQ_MIN = 1.0E-5
        TURB_FLUX_MIN = 1.0E-12
      else
        Q_SQ_MIN = 1.0E-5
        TURB_FLUX_MIN = 0.0
      end if

      if (DEBUG) print *, 'Done Set_my_qsq_min...'

    end subroutine Set_my_qsq_min


    subroutine Set_init_turb_state (restart,                                   &
              pbl3d_opt, pbl3d_prog, pbl3d_sfc_opt, pbl3d_l_opt,               &
              ids, ide, jds, jde, kds, kde,                                    &
              ims, ime, jms, jme, kms, kme,                                    &
              its, ite, jts, jte, kts, kte,                                    &
              fnm, fnp, dnw, ph, phb,                                          &
              t_2, qv, thetav, alt, u, v, ust, hfx, qfx, rmol,                 &
              q_sq, q_sq_prog, q_sq_vdiff, q_sq_hdiff, q_sq_shear,             &
              q_sq_buoyancy, q_sq_dissip,                                      &
              l_master, l_master_at_mass,                                      &
              turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,          &
              turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                  &
              turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v,      &
              turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,            &
              turb_flux_uqv, turb_flux_vqv, turb_flux_wqv                      )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Initialize turbulence state                                  !
      !                                                                       !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano (adapted from MYNN)     !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      logical, parameter :: DEBUG = .true.

      logical, intent(in) :: restart
      integer, intent(in) :: pbl3d_opt, pbl3d_prog, pbl3d_sfc_opt, pbl3d_l_opt
      integer :: i, j, k, itf, jtf, ktf
      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, dimension (kms:kme),                  intent(in)  :: fnm, fnp, dnw
      real, dimension(ims:ime, jms:jme),          intent(in)  :: ust, hfx, qfx, rmol
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)  :: t_2, qv, u, v, &
                                                                 alt, ph, phb
      real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: thetav
      real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: q_sq_prog, q_sq,     &
          q_sq_vdiff, q_sq_hdiff, q_sq_shear, q_sq_buoyancy, q_sq_dissip
      real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: l_master, l_master_at_mass, &
          turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,          &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                  &
          turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v,      &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,            &
          turb_flux_uqv, turb_flux_vqv, turb_flux_wqv

        ! Local vars
      real, dimension (kts:kte) :: sqv
      real, dimension (ims:ime, kms:kme, jms:jme) :: rdz, rdzw, z_at_mass, z_at_w, dz
      real, dimension (its:ite, kts:kte, jts:jte) :: sm, sh, rif
      real, dimension (its:ite, kts:kte, jts:jte) :: du_dz, dv_dz, &
                                                     dthetav_dz, dqv_dz
      real, dimension (ims:ime, kms:kme, jms:jme) :: rho
      integer :: i_start, i_end, j_start, j_end, init_iter

      if (restart) return

        ! Defines mass point at faces indexes
      jtf = Min (jte, jde - 1)
      itf = Min (ite, ide - 1)
      ktf = Min (kte, kde - 1)

      if (DEBUG) write (OUTPUT_UNIT, *) "ktf = ", ktf, "itf = ", itf, "jtf = ", jtf , "kte = ", kte, "kde = ", kde, "kme = ", kme

        ! Initializes thetav, TKE, L, and the turbulent fluxes
      j_loop_init: do j = jts, jtf
        i_loop_init: do i = its, itf
          sqv(kts:ktf) = qv(i, kts:ktf, j) / (1.0 + qv(i, kts:ktf, j))
          thetav(i, kts:ktf, j) = t_2(i, kts:ktf, j) * (1.0 + P608 * sqv(kts:ktf))
          q_sq(i, kts:ktf, j) = Q_SQ_MIN
          l_master(i, kts:ktf, j) = Q_SQ_MIN
          if (pbl3d_prog == 1) then
            l_master_at_mass(i, kts:ktf, j) = Q_SQ_MIN
          end if
          turb_flux_u2(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_v2(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_w2(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_uv(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_uw(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_vw(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_utheta_v(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_vtheta_v(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_wtheta_v(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_theta2_v(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_utheta(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_vtheta(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_wtheta(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_uqv(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_vqv(i, kts:ktf, j) = TURB_FLUX_MIN
          turb_flux_wqv(i, kts:ktf, j) = TURB_FLUX_MIN
          rho(i, kts:ktf, j) = 1. / alt(i, kts:ktf, j)
        end do i_loop_init
      end do j_loop_init

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      do j = j_start, j_end

        do k = 1, kte
          do i = i_start, i_end
            z_at_w(i, k, j) = ( ph(i, k, j) + phb(i, k, j) ) / g
          end do
        end do

        do k = 1, ktf
          do i = i_start, i_end
            z_at_mass(i, k, j) = ( z_at_w(i, k + 1, j) + z_at_w(i, k, j) ) / 2.
          end do
        end do

        do k = 1, ktf
          do i = i_start, i_end
            rdzw(i, k, j) = 1.0 / ( z_at_w(i, k + 1, j) - z_at_w(i, k, j) )
            dz(i, k, j) = z_at_w(i, k + 1, j) - z_at_w(i, k, j)
          end do
        end do

        do k = 2, ktf
          do i = i_start, i_end
            rdz(i, k, j) = 2.0 / ( z_at_w(i, k + 1, j) - z_at_w(i, k - 1, j) )
          end do
        end do

        do i = i_start, i_end
          rdz(i, 1, j) = 2. / (z_at_w(i, 2, j) - z_at_w(i, 1, j))
        end do

      end do

        ! Calc vertical derivatives
      call Calc_du_dz (u, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      call Calc_dv_dz (v, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      call Calc_dmassvar_dz (thetav, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dz)

      call Calc_dmassvar_dz (qv, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dz)

        ! Calc stability functions
      call Calc_sm_sh_l2_or_l2p5 (mass_or_face='face', my_level='l2', du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm, sh=sh, rif=rif)

      k = 1
      do j = j_start, j_end
        do i = i_start, i_end
          q_sq(i, k, j) = Max(Q_SQ_MIN, b1 ** (2.0 / 3.0) * ust(i, j) * ust(i, j))
        end do
      end do

!      do init_iter = 1, 4
          ! Calc q_sq
        call Calc_q_sq_l2_pbl_approx (mass_or_face='face', l_master=l_master, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
            jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
            jts=jts, jte=jte, kts=kts, kte=kte, q_sq_face=q_sq, sm=sm, rif=rif)

          ! Diagnose master length scale
        call Calc_l_my_algebra_l2_or_l2p5 (my_level='l2', pbl3d_l_opt=pbl3d_l_opt, q_sq=q_sq, dz=dz, rdzw=rdzw, rmol=rmol, &
            wthv=turb_flux_wtheta_v, l_master_l2=l_master, dthetav_dz=dthetav_dz, &
            ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

          ! Diagnose momentum, heat, moisture fluxes for pbl approx
      call Diagnose_fluxes_pbl_approx (fnm=fnm, fnp=fnp, u=u, v=v, t=t_2, thetav=thetav, rho=rho, qv=qv, ust=ust, hfx=hfx, qfx=qfx, &
          sm_tilde=sm, sh_tilde=sh, rif=rif, q_sq=q_sq, l_master=l_master, &
          du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, dqv_dz=dqv_dz, u2=turb_flux_u2, v2=turb_flux_v2, w2=turb_flux_w2, uv=turb_flux_uv, &
          uw=turb_flux_uw, vw=turb_flux_vw, uthv=turb_flux_utheta_v, vthv=turb_flux_vtheta_v, wthv=turb_flux_wtheta_v, th2v=turb_flux_theta2_v, &
          uqv=turb_flux_uqv, vqv=turb_flux_vqv, wqv=turb_flux_wqv, uth=turb_flux_utheta, vth=turb_flux_vtheta, wth=turb_flux_wtheta, ids=ids, ide=ide, &
          jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

!        call Diagnose_fluxes_surface (pbl3d_sfc_opt, u, v, t_2, thetav, rho, dz, qv, ust, hfx, qfx, &
!          rmol, q_sq, l_master, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, &
!          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
!          turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, &
!          turb_flux_wtheta, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
!      end do

      if (pbl3d_prog == 1) then
        call Fill_q_sq_prog_with_q_sq (q_sq_prog, q_sq, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          ! Diagnose master length scale at mass points for prognostic TKE
        call Calc_l_my_algebra_l2_or_l2p5 (my_level='l2p5', q_sq=q_sq_prog, z=z_at_mass, l_master_l2_in=l_master, &
            l_master_l2p5=l_master_at_mass, ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
            ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
      end if

      if (DEBUG) print *, 'Done Set_init_turb_state...'

    end subroutine Set_init_turb_state

  end module module_pbl3d
