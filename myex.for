      module Controller
      	use trcom
C
      	implicit none
C
      	real(kind=8), dimension(10,10) :: A0
      	real(kind=8), dimension(6,60) :: C
      	real(kind=8), dimension(10) :: K
      	real(kind=8), dimension(10,6) :: L
      	real(kind=8), dimension(60) :: gauss
C
C
      contains
C
      	function computeOutput(omega) result (y)
      		real(kind=8), dimension(60), intent(in) :: omega
      		real(kind=8), dimension(6) :: y
C
      		y = matmul(C, omega)
C
      	end function computeOutput
C
      	function computeControl(y) result (u)
      		real(kind=8), dimension(6), intent(in) :: y
      		real(kind=8) :: u
      		real(kind=8), dimension(10), save :: xhat = 0
C
      		u = -dot_product(K, xhat)
      		xhat = xhat + DT0 * (matmul(A0, xhat) + matmul(L, y))
C
      	end function computeControl
C
      	function applyControl(omega, I2) result (torque)
      		real(kind=8), dimension(60), intent(in) :: omega
      		real(kind=8), intent(in) :: I2
      		real(kind=8), dimension(60) :: torque
      		
      		torque = gauss * I2 * omega
C
      	end function applyControl
      	
      	subroutine loadData()
      		integer :: i
C
      		open(unit=101, file='/u/igoumiri/transp_kaye_method/control_data/A.dat', status='old')
      		do i=1,10
      			read(101,*) A(i,:)
      		end do
      		close(101)
C
      		open(unit=102, file='/u/igoumiri/transp_kaye_method/control_data/C.dat', status='old')
      		do i=1,6
      			read(102,*) C(i,:)
      		end do
      		close(102)
C
      		open(unit=103, file='/u/igoumiri/transp_kaye_method/control_data/K.dat', status='old')
      		do i=1,10
      			read(103,*) K(i,:)
      		end do
      		close(103)
C
      		open(unit=104, file='/u/igoumiri/transp_kaye_method/control_data/L.dat', status='old')
      		do i=1,10
      			read(104,*) L(i,:)
      		end do
      		close(104)
C
      		open(unit=105, file='/u/igoumiri/transp_kaye_method/control_data/gauss.dat', status='old')
      		do i=1,60
      			read(105,*) gauss(i,:)
      		end do
      		close(105)
C
      	end subroutine
C
      end module Controller
C
C
C
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
C
      if (inum == L3AUX .and. kpoint == 6) then
C
      	call loadData()
C
      end if
C
      if (inum == L3POSTEP .and. kpoint == 11) then
C
      	y = computeOutput(OMEGA)
C
      	u = computeControl(y)
C
      	TQANOM = applyControl(OMEGA, u)
C
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
