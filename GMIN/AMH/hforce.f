
C     --------------------- HFORCE.F ----------------------

      SUBROUTINE HFORCE(H_CORD,NITCORD,IDX1,IDX2,R1,R2,POT,FACTOR,
     *  I_CLASS,F_CORD)
 
C     --------------------------------------------------

C     HDRGN FINDS THE  POTENTIAL DUE TO HYDROGEN BONDS BETWEEN N AND O     

C     ---------------------------------------------------

      USE AMHGLOBALS,  ONLY:MAXSIZ,MAXCRD, PRCORD,HO_ZERO,NO_ZERO,SIGMA_NO,SIGMA_H,HBSCL

      IMPLICIT NONE


C     ARGUMENT DECLARATIONS:

          
              DOUBLE PRECISION  H_CORD(MAXSIZ,3),F_CORD(MAXSIZ,3,MAXCRD)

C     INTERNAL VARIABLES:

         INTEGER IDX1,IDX2 
C        --- DO LOOP INDICES ---

         INTEGER I_AXIS,I_CLASS 


         DOUBLE PRECISION  R1,FACTOR,
     *         POT,
     *         NITCORD(MAXSIZ,3),
     *         R2,
     *         DV_DRNO,DV_DRHO,DRNO_DO(3),
     *         DRNO_DN(3),DRHO_DO(3),DRHO_DH(3)

C     --------------------- BEGIN -----------------------


C     FIND FORCE DUE TO HBONDS

C DO *NOT* ZERO F_CORD HERE BECAUSE THE CALLER (HDRGN) IS WORKING OUT RUNNING TOTAL

        DV_DRNO = -HBSCL(I_CLASS)*POT*((R1 - NO_ZERO)/
     *                                    (SIGMA_NO**2))*FACTOR
        DV_DRHO = -HBSCL(I_CLASS)*POT*((R2 - HO_ZERO)/
     *                                    (SIGMA_H**2))*FACTOR

        DO I_AXIS = 1,3
        
        DRNO_DO(I_AXIS) = 
     *  (PRCORD(IDX1,I_AXIS,1,3)-NITCORD(IDX2,I_AXIS)) /R1 

        DRNO_DN(I_AXIS) = 
     *  -(PRCORD(IDX1,I_AXIS,1,3)-NITCORD(IDX2,I_AXIS)) /R1 

        DRHO_DO(I_AXIS) =
     *  (PRCORD(IDX1,I_AXIS,1,3)-H_CORD(IDX2,I_AXIS))/R2

        DRHO_DH(I_AXIS) = 
     *  -(PRCORD(IDX1,I_AXIS,1,3)-H_CORD(IDX2,I_AXIS))/R2

C        VEC1 IS THE FORCE VECTOR ACTING ON N

C        VEC1(I_AXIS) = -1.0*DRNO_DN(I_AXIS)*DV_DRNO
C        VEC2(I_AXIS) = -1.0*DV_DRHO*DRHO_DH(I_AXIS)
C        VEC3(I_AXIS)=
C    *  PRCORD(IDX1,I_AXIS,1,3)-NITCORD(IDX2,I_AXIS)

        F_CORD(IDX2,I_AXIS,1) = F_CORD(IDX2,I_AXIS,1) -
     *  DV_DRNO*DRNO_DN(I_AXIS)*0.7032820        
        F_CORD(IDX2,I_AXIS,1) = F_CORD(IDX2,I_AXIS,1) -
     *  DV_DRHO*DRHO_DH(I_AXIS)*0.8929599

        F_CORD(IDX2-1,I_AXIS,1) = F_CORD(IDX2-1,I_AXIS,1) -     
     *  DV_DRNO*DRNO_DN(I_AXIS)*0.4831806   
        F_CORD(IDX2-1,I_AXIS,1) = F_CORD(IDX2-1,I_AXIS,1) -
     *  DV_DRHO*DRHO_DH(I_AXIS)*0.8409657

        F_CORD(IDX1,I_AXIS,3) = F_CORD(IDX1,I_AXIS,3) -     
     *  DV_DRNO*DRNO_DO(I_AXIS)   
        F_CORD(IDX1,I_AXIS,3) = F_CORD(IDX1,I_AXIS,3) -
     *  DV_DRHO*DRHO_DO(I_AXIS)

        F_CORD(IDX2-1,I_AXIS,3) = F_CORD(IDX2-1,I_AXIS,3) + 
     *  DV_DRNO*DRNO_DN(I_AXIS)*0.1864626   
        F_CORD(IDX2-1,I_AXIS,3) = F_CORD(IDX2-1,I_AXIS,3) +
     *  DV_DRHO*DRHO_DH(I_AXIS)*0.7338894

        ENDDO  ! I_AXIS        


C     ----------------------- DONE -----------------------

      RETURN
      END
