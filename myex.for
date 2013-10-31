      module Controller
      	use trcom
C
      	implicit none
C
      	real(kind=8), dimension(10,10) :: A0
      	real(kind=8), dimension(10) :: B
      	real(kind=8), dimension(6,60) :: C
      	real(kind=8), dimension(10) :: K
      	real(kind=8), dimension(10,6) :: L
      	real(kind=8), dimension(6) :: F
      	real(kind=8), dimension(60) :: gauss
      	real(kind=8) :: u0
      	real(kind=8), dimension(6) :: r0
      	real(kind=8), dimension(6) :: rd
C
C
      contains
C
      	function computeOutput(omega) result (y)
      		real(kind=8), dimension(60), intent(in) :: omega
      		real(kind=8), dimension(6) :: y
C
      		y = matmul(C, omega) - r0
C
      	end function computeOutput
C
      	function computeControl(y) result (u)
      		real(kind=8), dimension(6), intent(in) :: y
      		real(kind=8) :: u, Frd
      		real(kind=8), dimension(10), save :: xhat = 0
C
      		Frd = dot_product(F, rd)
      		u = Frd - dot_product(K, xhat)
      		xhat = xhat + DT0 * (matmul(A0, xhat) + matmul(L, y) + B * Frd)
C
      	end function computeControl
C
      	function applyControl(omega, u) result (torque)
      		real(kind=8), dimension(60), intent(in) :: omega
      		real(kind=8), intent(in) :: u
      		real(kind=8), dimension(60) :: torque
C
      		torque = gauss * (u + u0) * omega
C
      	end function applyControl
C
      	subroutine loadData()
      		integer :: status, i
      		character(len=8) :: name
C
      		open(unit=100, status='old', form='formatted',
     &			file='/u/igoumiri/transp_kaye_method/parameters.dat')
      		status = 0
      		do while (status == 0)
      			read(100, *, iostat=status) name
      			select case (name)
      				case ('#', '!', '//', ';')
      					continue
      				case ('A')
      					do i=1,10
      						read(100,*) A0(i,:)
      					end do
      				case ('B')
      					do i=1,10
      						read(100,*) B(i)
      					end do
      				case ('C')
      					do i=1,6
      						read(100,*) C(i,:)
      					end do
      				case ('K')
      					do i=1,10
      						read(100,*) K(i)
      					end do
      				case ('L')
      					do i=1,10
      						read(100,*) L(i,:)
      					end do
      				case ('F')
      					do i=1,6
      						read(100,*) F(i)
      					end do
      				case ('gauss')
      					do i=1,60
      						read(100,*) gauss(i)
      					end do
      				case ('u0')
      					read(100,*) u0
      				case ('r0')
      					do i=1,6
      						read(100,*) r0(i)
      					end do
      				case ('rd')
      					do i=1,6
      						read(100,*) rd(i)
      					end do
      				case default
      					write(0,*) "Imene's controller: error in parameter file"
      			end select
      		end do
      		close(100)
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
