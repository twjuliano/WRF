!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_sf_spec_ideal_tsk

#ifdef DM_PARALLEL
    USE module_dm
#endif

CONTAINS

!-------------------------------------------------------------------
  
SUBROUTINE SPEC_IDEAL_TSK_HFX_LOC( TSK,                       &
                                   FLHC,T3D,                  &
                                   P3D,PSFC,                  &
                                   P1000mb,ROVCP,             &
                                   spec_hfx,                  &
                                   spec_sf_heatrate,          &
                                   dt,                        &
                                   spec_dumvar_2,             &
                                   ids,ide, jds,jde, kds,kde, &
                                   ims,ime, jms,jme, kms,kme, &
                                   its,ite, jts,jte, kts,kte  )
          
!-------------------------------------------------------------------
      IMPLICIT NONE

REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::  TSK

REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::   T3D
REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::   P3D

REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) ::  PSFC
REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) ::  FLHC

REAL,     INTENT(IN   ) :: P1000mb
REAL,     INTENT(IN   ) :: ROVCP

REAL,     INTENT(IN   ) :: spec_hfx        
REAL,     INTENT(IN   ) :: spec_sf_heatrate   
REAL,     INTENT(IN   ) :: dt                      
REAL,     INTENT(INOUT) :: spec_dumvar_2

INTEGER,  INTENT(IN   ) :: ids,ide, jds,jde, kds,kde, &
                           ims,ime, jms,jme, kms,kme, &
                           its,ite, jts,jte, kts,kte

! Local variables -------------------------------------------------------------

REAL, DIMENSION( its:ite, jts:jte ) ::   THX

INTEGER ::  i, j, i_end, j_end

REAL :: tconst, var1

!------------------------------------------------------------------------------

i_end = MIN(ite,ide-1)
j_end = MIN(jte,jde-1)

IF (spec_sf_heatrate .NE. 0.0 ) THEN

   IF (spec_hfx .NE. 0.0) THEN

      PRINT*,'spec_hfx must be zero when specifying a cooling rate'
      STOP

   ENDIF 
   
   spec_dumvar_2 = spec_dumvar_2 + dt

   IF (spec_dumvar_2 .GE. 0.0) THEN

      tconst = 60.0
      
      spec_dumvar_2 = spec_dumvar_2 - tconst

      DO J=jts,j_end
         DO I=its,i_end

            TSK(i,j) = TSK(i,j) + (spec_sf_heatrate/tconst)

         ENDDO
      ENDDO

   ENDIF

ELSE

   DO j=jts,j_end 
      DO i=its,i_end

         THX(i,j) = T3D(i,1,j)*( P1000mb/P3D(i,1,j) )**rovcp
         
      ENDDO
   ENDDO

   IF ( spec_hfx .EQ. 0.0 ) THEN
      
     DO j=jts,j_end 
        DO i=its,i_end
 
          TSK(i,j) = THX(i,j)/( P1000mb/PSFC(i,j) )**rovcp  
 
        ENDDO
     ENDDO

   ELSE
      
     DO j=jts,j_end 
        DO i=its,i_end

           var1 = ( spec_hfx/( MAX(FLHC(i,j),0.001) ) + THX(i,j) )/( P1000mb/PSFC(i,j) )**rovcp

           TSK(I,J) = MIN( MAX( var1,TSK(I,J)-1.0 ),TSK(I,J)+1.0 ) ! limit temperature change to 1K per time step

        ENDDO
    ENDDO

  ENDIF
  
ENDIF   


END SUBROUTINE SPEC_IDEAL_TSK_HFX_LOC
  
!-------------------------------------------------------------------

SUBROUTINE SPEC_IDEAL_TSK_HFX_MEAN( TSK,HFX,                   &
                                    FLHC,T3D,                  &
                                    P3D,PSFC,                  &
                                    P1000mb,ROVCP,             &
                                    spec_hfx,                  &
                                    spec_sf_heatrate,          &
                                    dt,                        &
                                    spec_dumvar_1,             &
                                    spec_dumvar_2,             &
                                    ids,ide, jds,jde, kds,kde, &
                                    ims,ime, jms,jme, kms,kme, &
                                    its,ite, jts,jte, kts,kte  )
          
!-------------------------------------------------------------------
      IMPLICIT NONE

REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::  TSK

REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::   T3D
REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::   P3D

REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) ::  PSFC
REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) ::  FLHC
REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) ::  HFX

REAL,     INTENT(IN   ) :: P1000mb
REAL,     INTENT(IN   ) :: ROVCP

REAL,     INTENT(IN   ) :: spec_hfx            
REAL,     INTENT(IN   ) :: spec_sf_heatrate   
REAL,     INTENT(IN   ) :: dt                      
REAL,     INTENT(INOUT) :: spec_dumvar_1
REAL,     INTENT(INOUT) :: spec_dumvar_2

INTEGER,  INTENT(IN   ) :: ids,ide, jds,jde, kds,kde, &
                           ims,ime, jms,jme, kms,kme, &
                           its,ite, jts,jte, kts,kte

! Local variables -------------------------------------------------------------

REAL, DIMENSION( its:ite, jts:jte ) ::   THX

INTEGER ::  i, j, i_end, j_end

REAL(8) :: flhc_patch_sum
REAL(8) :: flhc_dom_avg

REAL(8) :: ps_patch_sum
REAL(8) :: ps_dom_avg

REAL(8) :: hfx_patch_sum
REAL(8) :: hfx_dom_avg

REAL(8) :: thx_patch_sum
REAL(8) :: thx_dom_avg

REAL(8) :: tsk_patch_sum
REAL(8) :: tsk_dom_avg

REAL(8) :: fac

REAL :: tconst

REAL :: dumvar_1, dumvar_2

!------------------------------------------------------------------------------

i_end = MIN(ite,ide-1)
j_end = MIN(jte,jde-1)

fac = DBLE((ide-ids)*(jde-jds))

IF (spec_sf_heatrate .NE. 0.0 ) THEN

   IF (spec_hfx .NE. 0.0) THEN

      PRINT*,'spec_hfx must be zero when specifying a cooling rate'
      STOP

   ENDIF 
   
   spec_dumvar_2 = spec_dumvar_2 + dt
   
   IF (spec_dumvar_2 .GE. 0.0) THEN

      tconst = 60.0
      
      spec_dumvar_2 = spec_dumvar_2 - tconst

      tsk_patch_sum = 0.0d0
   
      DO j=jts,j_end 
         DO i=its,i_end

            tsk_patch_sum = tsk_patch_sum + TSK(i,j)
         
         ENDDO
      ENDDO

      tsk_dom_avg = tsk_patch_sum/fac

#ifdef DM_PARALLEL
   
      tsk_dom_avg  = wrf_dm_sum_real8(tsk_patch_sum/fac)

#endif   

      DO J=jts,j_end
         DO I=its,i_end

            TSK(i,j) =  tsk_dom_avg + (spec_sf_heatrate/tconst)

         ENDDO
      ENDDO

   ENDIF

ELSE

   DO j=jts,j_end 
      DO i=its,i_end

         THX(i,j) = T3D(i,1,j)*(P1000mb/P3D(i,1,j))**rovcp
         
      ENDDO
   ENDDO
   
   flhc_patch_sum = 0.0d0
   ps_patch_sum = 0.0d0
   hfx_patch_sum = 0.0d0
   thx_patch_sum = 0.0d0

   DO j=jts,j_end 
      DO i=its,i_end

         flhc_patch_sum = flhc_patch_sum + FLHC(i,j) 

         ps_patch_sum = ps_patch_sum + PSFC(i,j) 

         hfx_patch_sum = hfx_patch_sum + HFX(i,j) 
         
         thx_patch_sum = thx_patch_sum + THX(i,j) 

      ENDDO
   ENDDO
      
   flhc_dom_avg = flhc_patch_sum/fac
   ps_dom_avg = ps_patch_sum /fac
   hfx_dom_avg = hfx_patch_sum/fac
   thx_dom_avg = thx_patch_sum/fac

#ifdef DM_PARALLEL
   
      flhc_dom_avg  = wrf_dm_sum_real8(flhc_patch_sum/fac)

      ps_dom_avg  = wrf_dm_sum_real8(ps_patch_sum/fac)

      hfx_dom_avg  = wrf_dm_sum_real8(hfx_patch_sum/fac)

      thx_dom_avg = wrf_dm_sum_real8(thx_patch_sum/fac)
      
#endif

      
   IF (flhc_dom_avg .EQ. 0.0) THEN
         
      ! do nothing: flhc is zero at startup

   ELSE IF ( spec_hfx .EQ. 0.0 ) THEN

      dumvar_2 = thx_dom_avg/( P1000mb/ps_dom_avg )**rovcp 
      
      DO j=jts,j_end 
         DO i=its,i_end
            
            TSK(i,j) = dumvar_2 
            
         ENDDO
      ENDDO
         
   ELSE
         
      dumvar_1 = spec_hfx - hfx_dom_avg
         
      IF ( spec_dumvar_1 .EQ. 0.0) spec_dumvar_1 = spec_hfx
      
      spec_dumvar_1 = spec_dumvar_1 + dumvar_1   
      
      dumvar_2 = (spec_dumvar_1/flhc_dom_avg + thx_dom_avg)/(P1000mb/(ps_dom_avg))**rovcp 
      
      DO j=jts,j_end 
         DO i=its,i_end
            
            TSK(i,j) = dumvar_2 
            
         ENDDO
      ENDDO
      
   ENDIF

ENDIF   

   
END SUBROUTINE SPEC_IDEAL_TSK_HFX_MEAN

!-------------------------------------------------------------------

SUBROUTINE SET_MOL_UST_Z0 (ids,ide, jds,jde, kds,kde,                    &
                           ims,ime, jms,jme, kms,kme,                    &
                           its,ite, jts,jte, kts,kte,                    &
                           znt,                                          &
                           spec_z0,                                      &
                           ust,mol,                                      &
                           ustt,molt                                     )
   

!-------------------------------------------------------------------

IMPLICIT NONE


INTEGER,  INTENT(IN   ) :: ids,ide, jds,jde, kds,kde, &
                           ims,ime, jms,jme, kms,kme, &
                           its,ite, jts,jte, kts,kte

REAL,     INTENT(IN   ) :: spec_z0
REAL,     DIMENSION( ims:ime, jms:jme ),        INTENT(INOUT) :: znt
REAL,     DIMENSION( ims:ime, jms:jme ),        INTENT(INOUT) :: ust
REAL,     DIMENSION( ims:ime, jms:jme ),        INTENT(INOUT) :: mol
REAL,     DIMENSION( ims:ime, jms:jme ),        INTENT(INOUT) :: ustt
REAL,     DIMENSION( ims:ime, jms:jme ),        INTENT(INOUT) :: molt

INTEGER ::  i, j, i_start, i_end, j_start, j_end

DO j = jms, jme
   DO i = ims, ime
      
      znt(i,j) = spec_z0
      ust(i,j) = ustt(i,j)
      mol(i,j) = molt(i,j)
      
   ENDDO
ENDDO


END SUBROUTINE SET_MOL_UST_Z0


SUBROUTINE SET_MOLT_USTT (ids,ide, jds,jde, kds,kde,                    &
                          ims,ime, jms,jme, kms,kme,                    &
                          its,ite, jts,jte, kts,kte,                    &
                          ust,mol,                                      &
                          ustt,molt                                     )
   

!-------------------------------------------------------------------

IMPLICIT NONE


INTEGER,  INTENT(IN   ) :: ids,ide, jds,jde, kds,kde, &
                           ims,ime, jms,jme, kms,kme, &
                           its,ite, jts,jte, kts,kte

REAL,     DIMENSION( ims:ime, jms:jme ),        INTENT(INOUT) :: ust
REAL,     DIMENSION( ims:ime, jms:jme ),        INTENT(INOUT) :: mol
REAL,     DIMENSION( ims:ime, jms:jme ),        INTENT(INOUT) :: ustt
REAL,     DIMENSION( ims:ime, jms:jme ),        INTENT(INOUT) :: molt

INTEGER ::  i, j

DO j = jms, jme
   DO i = ims, ime
      
      ustt(i,j) = ust(i,j)
      molt(i,j) = mol(i,j)
      
   ENDDO
ENDDO


END SUBROUTINE SET_MOLT_USTT
  
  
END MODULE module_sf_spec_ideal_tsk


    
