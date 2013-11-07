      module Controller
      use trcom, only: TA, DT0, NZONES

      implicit none

      save     ! save module variables

      integer, parameter :: r8 = selected_real_kind(12,100)
      integer, parameter :: kunit = 100 ! probably safe unit number
      integer, parameter :: max_path_len = 512 ! probably long enough
      character(len=45) :: parameters =
     &   '/u/igoumiri/transp_kaye_method/parameters.nml'
      character(len=max_path_len) :: datafile
      character(len=max_path_len) :: statefile, controlfile, timefile
      namelist /input_files/ datafile
      namelist /output_files/ statefile, controlfile, timefile

      ! --- state changes during run ---
      real(kind=r8), dimension(:), allocatable :: xhat

      ! --- fixed during run ---
      real(kind=r8), dimension(:,:), allocatable :: A0
      real(kind=r8), dimension(:,:), allocatable :: B
      real(kind=r8), dimension(:,:), allocatable :: C
      real(kind=r8), dimension(:,:), allocatable :: K
      real(kind=r8), dimension(:,:), allocatable :: L
      real(kind=r8), dimension(:,:), allocatable :: F
      real(kind=r8), dimension(:), allocatable :: u0
      real(kind=r8), dimension(:), allocatable :: r0
      real(kind=r8), dimension(:), allocatable :: rd
      real(kind=r8), dimension(:,:), allocatable :: gauss

      contains

      ! ----------------------------------------------------------------
      function computeOutput(omega) result (y)
      real(kind=r8), dimension(NZONES), intent(in) :: omega
      real(kind=r8), dimension(size(C, 1)) :: y

      y = matmul(C, omega) - r0

      end function computeOutput

      ! ----------------------------------------------------------------
      function computeControl(y) result (u)
      real(kind=r8), dimension(size(C, 1)), intent(in) :: y
      real(kind=r8), dimension(size(B, 2)) :: u, Frd

      Frd = matmul(F, rd)
      u = Frd - matmul(K, xhat)
      xhat = xhat
     &   + DT0 * (matmul(A0, xhat) + matmul(L, y) + matmul(B, Frd))

      end function computeControl

      ! ----------------------------------------------------------------
      function applyControl(omega, u) result (torque)
      real(kind=r8), dimension(NZONES), intent(in) :: omega
      real(kind=r8), dimension(size(B, 2)), intent(in) :: u
      real(kind=r8), dimension(NZONES) :: torque

      torque = matmul(gauss, u + u0) * omega

      end function applyControl

      ! ----------------------------------------------------------------
      subroutine initController(u, y, restart)
      real(kind=r8), dimension(:), allocatable, save :: u, y
      logical, intent(in), optional :: restart

      ! Read paths from namelists
      open(unit=kunit, file=parameters, status='old')
      read(kunit, nml=input_files)
      rewind(kunit)
      read(kunit, nml=output_files)
      close(kunit)

      call loadData

      ! Right now we load directly A-BK-LC from the datafile but later
      ! it might be more natural to load A and do the operation here
      !A0 = A0 - matmul(B, K) - matmul(L, C)

      allocate(xhat(size(A0, 1)))
      allocate(u(size(B, 2)))
      allocate(y(size(C, 1)))

      ! If restarting, read state
      ! Otherwise, init xhat and create output files
      if (present(restart) .and. restart) then
         ! %->inform, !->warn, ?->fatal error
         print *, "%controller: "
     &      // "reading state from '" // trim(statefile)  // "'"
         call readState
      else
         print *, '%controller: initializing state to zero'
         xhat = 0._r8

         open(unit=kunit, file=statefile, status='replace')
         write(kunit,*) '! state - size: ', size(xhat)
         close(kunit)

         open(unit=kunit, file=controlfile, status='replace')
         write(kunit,*) '! control - size: ', size(u)
         close(kunit)

         open(unit=kunit, file=timefile, status='replace')
         write(kunit,*) '! time - size: 1'
         close(kunit)

         call saveAll(u)
      end if

      end subroutine initController

      ! ----------------------------------------------------------------
      subroutine saveAll(u)
      real(kind=r8), dimension(size(B, 2)), intent(in) :: u
      integer :: i

      ! save time
      open(unit=kunit, file=timefile,
     &   status='old', position='append', action='write')
      write(kunit,*) TA
      close(kunit)

      ! save xhat
      open(unit=kunit, file=statefile,
     &   status='old', position='append', action='write')
      write(kunit,*) xhat
      close(kunit)

      ! save control
      open(unit=kunit, file=controlfile,
     &   status='old', position='append', action='write')
      write(kunit,*) u
      close(kunit)

      end subroutine saveAll

      ! ----------------------------------------------------------------
      subroutine readState

      open(unit=kunit, file=statefile,
     &   status='old', position='append', action='read')
      backspace(kunit)
      read(kunit,*) xhat
      close(kunit)

      end subroutine readState

      ! ----------------------------------------------------------------
      subroutine loadData
      integer :: i, rows, cols
      character(len=5) :: name

      open(unit=kunit, file=datafile, status='old')
      do
         read(kunit, *, end=50) name

         select case (name)
         case ('A')
            backspace(kunit)
            read(kunit, *) name, rows, cols
            allocate(A0(rows, cols))
            do i=1,rows
               read(kunit,*) A0(i,:)
            end do
         case ('B')
            backspace(kunit)
            read(kunit,*) name, rows, cols
            allocate(B(rows, cols))
            do i=1,rows
               read(kunit,*) B(i,:)
            end do
         case ('C')
            backspace(kunit)
            read(kunit,*) name, rows, cols
            allocate(C(rows, cols))
            do i=1,rows
               read(kunit,*) C(i,:)
            end do
         case ('K')
            backspace(kunit)
            read(kunit,*) name, rows, cols
            allocate(K(rows, cols))
            do i=1,rows
               read(kunit,*) K(i,:)
            end do
         case ('L')
            backspace(kunit)
            read(kunit,*) name, rows, cols
            allocate(L(rows, cols))
            do i=1,rows
               read(kunit,*) L(i,:)
            end do
         case ('F')
            backspace(kunit)
            read(kunit,*) name, rows, cols
            allocate(F(rows, cols))
            do i=1,rows
               read(kunit,*) F(i,:)
            end do
         case ('gauss')
            backspace(kunit)
            read(kunit,*) name, rows, cols
            allocate(gauss(rows, cols))
            do i=1,rows
               read(kunit,*) gauss(i,:)
            end do
         case ('u0')
            backspace(kunit)
            read(kunit,*) name, rows
            allocate(u0(rows))
            do i=1,rows
               read(kunit,*) u0(i)
            end do
         case ('r0')
            backspace(kunit)
            read(kunit,*) name, rows
            allocate(r0(rows))
            do i=1,rows
               read(kunit,*) r0(i)
            end do
         case ('rd')
            backspace(kunit)
            read(kunit,*) name, rows
            allocate(rd(rows))
            do i=1,rows
               read(kunit,*) rd(i)
            end do
         case default
            if (scan(name, '#!;') == 1) then ! a comment
               cycle
            else
               write(0,*) "?controller: "
     &            // "error in data file '" // trim(datafile) // "'"
               call bad_exit
            end if
         end select
      end do
 50   close(kunit)

      ! check that:
      ! - nothing is missing
      ! - the size of all parameters are consistent
      call checkConsistency

      end subroutine loadData

      ! ----------------------------------------------------------------
      subroutine checkConsistency

         ! check that all variables were provided
         call assert(allocated(A0), 'A must be provided')
         call assert(allocated(B), 'B must be provided')
         call assert(allocated(C), 'C must be provided')
         call assert(allocated(K), 'K must be provided')
         call assert(allocated(L), 'L must be provided')
         call assert(allocated(F), 'F must be provided')
         call assert(allocated(u0), 'u0 must be provided')
         call assert(allocated(r0), 'r0 must be provided')
         call assert(allocated(rd), 'rd must be provided')
         call assert(allocated(gauss), '"gauss" must be provided')

         ! check dimensions
         call assert(size(A0, 1) == size(A0, 2), 'A must be square')
         call assert(size(B, 1) == size(A0, 1),
     &      'B must have as many rows as A')
         call assert(size(C, 2) == NZONES,
     &      'the number of columns of C must be equal to NZONES')
         call assert(size(K, 1) == size(B, 2),
     &      'K must have as many rows as the number of columns of B')
         call assert(size(K, 2) == size(A0, 2),
     &      'K must have as many columns as A')
         call assert(size(L, 1) == size(A0, 1),
     &      'L must have as many rows as A')
         call assert(size(L, 2) == size(C, 1),
     &      'L must have as many columns as the number of rows of C')
         call assert(size(F, 1) == size(B, 2),
     &      'L must have as many rows as the number of columns of B')
         call assert(size(F, 2) == size(C, 1),
     &      'L must have as many columns as the number of rows of C')
         call assert(size(u0) == size(B, 2),
     &   'the length of u0 must be equal to the number of columns of B')
         call assert(size(r0) == size(C, 1),
     &      'the length of r0 must be equal to the number of rows of C')
         call assert(size(rd) == size(C, 1),
     &      'the length of rd must be equal to the number of rows of C')
         call assert(size(gauss, 1) == size(B, 2),
     &   '"gauss" must have as many rows as the number of columns of B')
         call assert(size(gauss, 2) == NZONES,
     &      'the number of columns of "gauss" must be equal to NZONES')

      end subroutine checkConsistency

      ! ----------------------------------------------------------------
      subroutine assert(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not.condition) then
         write(0,*) '?controller: ' // message
         call bad_exit
      end if

      end subroutine assert

      end module Controller
C
C
C
C******************** START FILE EXPERT.FOR ; GROUP EXPERT *************
C***********************************************************************
C
C
      subroutine expert(KCLASS,KSUB,KPOINT)
C
C     0.04	MODIFY STANDARD OPERATION
C
C
      use trcom, only: NCLASS,NSUB,NPOINT, L3POSTEP, OMEGA, TQANOM
      use trcom, only: NLREPT, LCENTR, LEDGE, LCM1, LEP1, DVOL, NSTEP
      use trcom, only: L3AUXVAL, NMODVPH, NZONES

      use Controller
C
      implicit none

      real(kind=r8), dimension(:), allocatable :: u
      real(kind=r8), dimension(:), allocatable :: y

      integer :: INUM, INCALL, IPCALL, INTES, IPTES
      integer :: KCLASS, KSUB, KPOINT
      logical :: IEXPER

      logical :: idone = .false.  ! flag if initialization has been performed
C
C-----------------------------------------------------------------------
C
      IEXPER(INCALL,IPCALL,INTES,IPTES) =
     >     ((INCALL.EQ.INTES).AND.(IPCALL.EQ.IPTES))
C
C
C-----------------------------------------------------------------------
C
C
C     SET TRACER VARIABLES
C
      NCLASS=KCLASS
      NSUB=KSUB
      NPOINT=KPOINT
C
C     ARE DIAGNOSTICS REQUIRED?
C
      if (NLREPT) call REPORT(KCLASS,KSUB,KPOINT)
C
C
C     DO SPECIAL FUNCTIONS
C
      INUM=1000*KCLASS+KSUB
C
C
C     EXAMPLE:  TEST TO DO SOMETHING ON EXPERT CALL FROM STEPON, POINT 2.
C
C     SAMPLE      IF(IEXPER(INUM,KPOINT,L3STEPON,2)) THEN
C     SAMPLE          ....
C     SAMPLE      ENDIF
C
      if (inum == L3AUXVAL .and. kpoint == 6) then
         ! This part is run once at the very first step

         call initController(u, y)

         idone = .true.   ! on first call only
      end if

      if (inum == L3POSTEP .and. kpoint == 11 .and. NSTEP>1) then
         if (.not. idone) then
            ! This part is run once after a restart

            call initController(u, y, restart=.true.)

            idone = .true.   ! on first call only
         end if

         if (NMODVPH/=0) then  ! only call when predicting rotation

            y = computeOutput(OMEGA(LCENTR:LEDGE,2))

            u = computeControl(y)

            ! TQANOM(I) is the torque applied to a zone
            ! multiply (Nt*m/cm^3) by zone volume
            TQANOM(LCENTR:LEDGE) = NZONES * DVOL(LCENTR:LEDGE,2)
     &         * applyControl(OMEGA(LCENTR:LEDGE,2), u)
            TQANOM(LCM1) = TQANOM(LCENTR)
            TQANOM(LEP1) = TQANOM(LEDGE)
         end if

         call saveAll(u) ! saving each time step is not good
                         ! TODO: need to change this
      end if
C
C
C
C
      return
C
C
C
C
      end
C********************END FILE EXPERT.FOR ; GROUP EXPERT ****************
