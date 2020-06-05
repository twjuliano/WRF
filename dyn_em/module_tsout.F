!WRF
!
MODULE module_tsout

USE module_model_constants

CONTAINS

SUBROUTINE compute_vertical_slice( u_slice, v_slice,                          &
                                   w_slice, t_slice,                          &
                                   u, v, w, t, z, z_at_w,                     &
                                   hgt,                                       &
                                   ph,phb,                                    &
                                   slice_height,                              &
                                   ids, ide, jds, jde, kds, kde,              &
                                   ims, ime, jms, jme, kms, kme,              &
                                   its, ite, jts, jte, kts, kte               )

!=======================================================================
!
!  This subroutine calculates interpolated horizontal slices at specified heights.
!
!======================================================================= 

   IMPLICIT NONE

   REAL, DIMENSION(ims:ime,jms:jme), INTENT( OUT ) :: u_slice, v_slice, &
                                                      w_slice, t_slice
   
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(  IN ) :: u, v, w, t, z, z_at_w
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(  IN ) :: ph,phb
   REAL, DIMENSION(ims:ime,jms:jme), INTENT(  IN ) :: hgt
   
   REAL, INTENT(  IN ) :: slice_height
   
   INTEGER , INTENT( IN  ) :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte

   ! local variables
   

   INTEGER :: i, j, k, k_max

   REAL :: z_1, z_2, d_1, d_2, d_z, fac_1, fac_2
   REAL :: slice_height_abs
   REAL :: u_m_1, u_m_2, v_m_1, v_m_2
   REAL :: z_w_1,z_w_2

   
!-----------------------------------------------------------------------
! executable code starts here
!
   
   
!   DO i = its, ite
!      DO j = jts, MIN(jte,jde-1)
!
!         DO k = kts+1, MIN(kte,kde-1)  
!
!            !z_2 = 0.5*( z(i-1,k,j) + z(i,k,j) ) ! project z to u points
!            z_2 =  z(i,k,j) ! slices are interpolated to cell centers
!           
!            slice_height_abs = slice_height + hgt(i,j)
!            IF ( z_2 .GE. slice_height_abs ) THEN
!
!               z_1 = z(i,k-1,j) 
!
!               d_z = z_2 - z_1                  ! linear interpolation in height
!               d_2 = z_2 - slice_height_abs
!!               u_slice(i,j) = d_2 
!               d_1 = slice_height_abs - z_1
!               fac_2 = d_2/d_z
!               fac_1 = d_1/d_z
!               u_m_2 = 0.5*(u(i-1,k,j) + u(i,k,j) ) ! u at cell center level k
!               u_m_1 = 0.5*(u(i-1,k-1,j) + u(i,k-1,j) ) ! u at cell center level k-1
!               
!               u_slice(i,j) = fac_2*u_m_1 + fac_1*u_m_2
!                
!               GOTO 10
!               
!            ENDIF
!            
!         END DO
!
!10       CONTINUE
!         
!      END DO
!   END DO
!
!   
!   DO i = its, MIN(ite,ide-1)
!      DO j = jts, jte
!
!         DO k = kts+1, MIN(kte,kde-1)
!
!            z_2 =  z(i,k,j) ! slices are interpolated to cell centers
!
!            slice_height_abs = slice_height + hgt(i,j)
!            IF ( z_2 .GE. slice_height_abs ) THEN
!               z_1 = z(i,k-1,j)
!
!               d_z = z_2 - z_1
!               d_2 = z_2 - slice_height_abs
!               d_1 = slice_height_abs - z_1
!               fac_2 = d_2/d_z
!               fac_1 = d_1/d_z
!               
!               v_m_2 = 0.5*(v(i,k,j-1) + v(i,k,j) ) ! v at cell center level k
!               v_m_1 = 0.5*(v(i,k-1,j-1) + v(i,k-1,j) ) ! v at cell center level k-1
!               v_slice(i,j) = fac_2*v_m_1 + fac_1*v_m_2
!                
!               GOTO 20
!               
!            ENDIF
!
!         END DO
!
!20       CONTINUE
!         
!      END DO
!   END DO
!
!   
   DO i = its, MIN(ite,ide-1)
      DO j = jts, MIN(jte,jde-1)

         DO k = kts+1, MIN(kte,kde-1)

            z_2 =  z(i,k,j) ! slices are interpolated to cell centers
            z_w_2 = 0.5*( z(i,k,j) + z(i,k+1,j) ) 

            slice_height_abs = slice_height + hgt(i,j)


            IF ( z_w_2 .GE. slice_height_abs ) THEN
               z_w_1 = 0.5*(z(i,k-1,j) + z(i,k,j)) 

               d_z = z_w_2 - z_w_1
               d_2 = z_w_2 - slice_height_abs
               d_1 = slice_height_abs - z_w_1
               fac_2 = d_2/d_z
               fac_1 = d_1/d_z
               
               w_slice(i,j) = fac_2*w(i,k-1,j) + fac_1*w(i,k,j)

               GOTO 20
               
            ENDIF

20       CONTINUE

            IF ( z_2 .GE. slice_height_abs ) THEN
               z_1 = z(i,k-1,j)

               d_z = z_2 - z_1
               d_2 = z_2 - slice_height_abs
               d_1 = slice_height_abs - z_1
               fac_2 = d_2/d_z
               fac_1 = d_1/d_z
               
               u_m_2 = 0.5*(u(i+1,k,j) + u(i,k,j) ) ! u at cell center level k
               u_m_1 = 0.5*(u(i+1,k-1,j) + u(i,k-1,j) ) ! u at cell center level k-1
               v_m_2 = 0.5*(v(i,k,j+1) + v(i,k,j) ) ! v at cell center level k
               v_m_1 = 0.5*(v(i,k-1,j+1) + v(i,k-1,j) ) ! v at cell center level k-1

               u_slice(i,j) = fac_2*u_m_1 + fac_1*u_m_2
               v_slice(i,j) = fac_2*v_m_1 + fac_1*v_m_2
               t_slice(i,j) = fac_2*t(i,k-1,j) + fac_1*t(i,k,j)
                
               GOTO 30
               
            ENDIF

         END DO

30       CONTINUE
         
      END DO
   END DO




   DO i = its, MIN(ite,ide-1)
      DO j = jts, MIN(jte,jde-1)

         DO k = kts+1, MIN(kte,kde-1)

            z_2 =  z(i,k,j) ! slices are interpolated to cell centers
            z_w_2 = 0.5*( z(i,k,j) + z(i,k+1,j) ) 

            slice_height_abs = slice_height + hgt(i,j)


            IF ( z_w_2 .GE. slice_height_abs ) THEN
               z_w_1 = 0.5*(z(i,k-1,j) + z(i,k,j)) 

               d_z = z_w_2 - z_w_1
               d_2 = z_w_2 - slice_height_abs
               d_1 = slice_height_abs - z_w_1
               fac_2 = d_2/d_z
               fac_1 = d_1/d_z
               
               w_slice(i,j) = fac_2*w(i,k-1,j) + fac_1*w(i,k,j)

               GOTO 40
               
            ENDIF

         END DO

40       CONTINUE
         
      END DO
   END DO
     
   
   RETURN

END SUBROUTINE compute_vertical_slice


  
SUBROUTINE cell_centered_calc( uts, vts, wts, tts, pts, tkets, zts,       &
                               t11ts, t22ts, t33ts, t12ts, t13ts, t23ts,  &
                               h3ts, usts, z0ts,                          &
                               u, v, w, t, p, pb, tke, ph, phb,           &
                               t11, t22, t33, t12, t13, t23,              &
                               h3, ust, z0,                               &
                               ids, ide, jds, jde, kds, kde,              &
                               ims, ime, jms, jme, kms, kme,              &
                               its, ite, jts, jte, kts, kte                )

!=======================================================================
!
!  This subroutine calculates cell centered parameters for reduced output
!
!======================================================================= 

   IMPLICIT NONE

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT( OUT ) :: uts, vts, wts, &
                                                       tts, pts, tkets, zts, &
                                                       t11ts, t22ts, t33ts,  &
                                                       t12ts, t13ts, t23ts,  &
                                                       h3ts 

   REAL, DIMENSION(ims:ime,jms:jme), INTENT( OUT ) :: usts, z0ts

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT( IN ) :: u, v, w,           &
                                                             t, p, pb,          &
                                                             tke, ph, phb, t11, &
                                                             t22, t33, t12,     & 
                                                             t13, t23, h3 
                                                         
   REAL, DIMENSION(ims:ime,jms:jme), INTENT( IN ) :: ust, z0 
   
   INTEGER , INTENT( IN  ) :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte

! local variables
   INTEGER :: i, j, k

!-----------------------------------------------------------------------
! executable starts here
!

   DO i = its, ite
      DO j = jts, jte
         
         DO k=kts,MIN(kte,kde-1)
            uts(i,k,j)  = 0.5 * ( u(i,k,j) + u(i+1,k,j) )
            vts(i,k,j)  = 0.5 * ( v(i,k,j) + v(i,k,j+1) )
            wts(i,k,j)  = 0.5 * ( w(i,k,j) + w(i,k+1,j) )
            tts(i,k,j)  = t(i,k,j)
            pts(i,k,j)  = p(i,k,j) + pb(i,k,j)
            tkets(i,k,j) = tke(i,k,j)
            t11ts(i,k,j) = t11(i,k,j)
            t22ts(i,k,j) = t22(i,k,j)
            t33ts(i,k,j) = t33(i,k,j)
            t12ts(i,k,j) = 0.25 * ( t12(i,k,j)   + t12(i+1,k,j) + t12(i,k,j+1) + t12(i+1,k,j+1) )
            t13ts(i,k,j) = 0.5 * ( t13(i,k,j) + t13(i+1,k,j) ) 
            t23ts(i,k,j) = 0.5 * ( t23(i,k,j) + t23(i,k,j+1) ) 
            h3ts(i,k,j)  = h3(i,k,j)
            zts(i,k,j)   = ( phb(i,k,j) + ph(i,k,j) )/g
         ENDDO
         
         t13ts(i,kde,j) = 0.5 * ( t13(i,kde,j) + t13(i+1,kde,j) ) 
         t23ts(i,kde,j) = 0.5 * ( t23(i,kde,j) + t23(i,kde,j+1) ) 
         h3ts(i,kde,j)  = h3(i,kde,j) 
         zts(i,kde,j)   = ( phb(i,kde,j) + ph(i,kde,j) )/g
                   
         usts(i,j) = ust(i,j)  
         z0ts(i,j) = z0(i,j)  
         
      END DO
   END DO
   
   RETURN

END SUBROUTINE cell_centered_calc

END MODULE module_tsout
