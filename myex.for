C******************** START FILE EXPERT.FOR ; GROUP EXPERT ******************
C*****************************************************************
C
C
      SUBROUTINE EXPERT(KCLASS,KSUB,KPOINT)
C
C	0.04	MODIFY STANDARD OPERATION
C
C
      use trcom
      use Controller
C
      real(kind=8), dimension(6) :: y
      real(kind=8) :: u
      LOGICAL IEXPER
C
C-----------------------------------------------------------------
C
      IEXPER(INCALL,IPCALL,INTES,IPTES) =
     >     ((INCALL.EQ.INTES).AND.(IPCALL.EQ.IPTES))
C
C
C-----------------------------------------------------------------
C
C
C	SET TRACER VARIABLES
C
      NCLASS=KCLASS
      NSUB=KSUB
      NPOINT=KPOINT
C
C	ARE DIAGNOSTICS REQUIRED?
C
      IF (NLREPT)CALL REPORT(KCLASS,KSUB,KPOINT)
C
C
C	DO SPECIAL FUNCTIONS
C
      INUM=1000*KCLASS+KSUB
C
C
C  EXAMPLE:  TEST TO DO SOMETHING ON EXPERT CALL FROM STEPON, POINT 2.
C
CSAMPLE      IF(IEXPER(INUM,KPOINT,L3STEPON,2)) THEN
CSAMPLE          ....
CSAMPLE      ENDIF

      if (inum == L3POSTEP .and. kpoint == 11) then

         y = computeOutput(OMEGA)
         
		 u = computeControl(y)
         
         TQANOM = applyControl(OMEGA, u)

      end if
C
C
C
C
      RETURN
C
C
C
C
      END
C******************** END FILE EXPERT.FOR ; GROUP EXPERT ******************
