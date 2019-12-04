MODULE parameters_test

IMPLICIT NONE

  REAL(8), PARAMETER :: pi = 4*ATAN(1.d0)
 
  INTEGER(8), PARAMETER :: lightspeed = 3e10 !cm/s

  REAL(8), PARAMETER :: dt = 0.04 ! Time Step
  INTEGER, PARAMETER :: maxn = 1000 ! Default # of particles
  INTEGER, PARAMETER :: maxt = 10 ! Max time

CONTAINS

! ---------------------- BEGIN weighf -----------------------------------
!
! The weighting function/kernel normalized to unity which converges to a 
! delta function as h -> 0.
!
! -----------------------

FUNCTION weighf(pos1, pos2, h) ! Monaghan (1992)

! ---------------------- Variables  
!
! pos1, pos2: {real - input}
!    Position 1 and position 2 
!
! h: {real - input}
!    Smoothing length.
!
! ----------------------

IMPLICIT NONE

  REAL(8), intent(in) :: pos1, pos2, h
!  CHARACTER(10), intent(in) :: method
  REAL(8) :: pos
  REAL(8), PARAMETER :: pi = 4*ATAN(1.d0)
  REAL(8) :: weighf

  weighf = 0.d0

  pos = pos1-pos2

!  if (method == "gauss") then

  weighf = exp(-0.5*pos*pos/(h*h))/(h*h*pi) ! 2D Gaussian

!  else  ! Alternative Kernel also from Monaghan (1992)
!  q = abs(pos1-pos2)/h

!  weighf = 10/(h*h*7*pi)

!  if (0 <= q .and. q <= 1) then

!    weighf = weighf*(1-1.5*q*q+0.75*q*q*q*q)

!  elseif (0 <= q .and. q <= 1) then

!    weighf = weighf*0.25*(2-q)*(2-q)*(2-q)

!  elseif (2 <= q) then

!    weighf = 0

!  else

!    print*, "Something went wrong with weighting kernel! Aborting."

!    STOP

!  endif 

!  endif

END FUNCTION  

!------------------------ END weighf ---------------------------------

END MODULE
