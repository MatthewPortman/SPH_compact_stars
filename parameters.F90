MODULE parameters

IMPLICIT NONE

  REAL(8), PARAMETER :: pi = 4*ATAN(1.d0)
 
  INTEGER(8), PARAMETER :: G = 6.67e-8 !cm^2 g^-1 s^-2
  INTEGER(8), PARAMETER :: lightspeed = 3e10 !cm/s
  REAL(8), PARAMETER :: gas_const = 50.1e-15 ! g cm^2 s^-2 K^-1 

  REAL(8), PARAMETER :: star_mass = 4e33 !2 mstar in grams  
!  REAL(8), PARAMETER :: part_mass = 1 !n = 1.675e-24 g
!  REAL(8), PARAMETER :: part_dens = 1 !
  REAL(8), PARAMETER :: dt = 0.04 ! Time Step
  INTEGER, PARAMETER :: maxn = 10000 ! Default # of particles
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
!    Position 1 and position 2 (distance not squared... consider keeping squared to avoid another complication if nothing else uses regular r).
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

!  else
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
!------------------------ BEGIN gweighf ---------------------------------
!
! The gradient of the weighting function/kernel.
!
! Insert finite difference here...
! -----------------------

attributes(global) SUBROUTINE gweighf(pos, h, gradx, grady) ! Monaghan (1992)

use cudafor
! ---------------------- Variables  
!
! x/ypos1, x/ypos2: {real - input}
!    Coord positions of particles 1 and 2 
!
! h: {real - input}
!    Smoothing length.
!
! ----------------------

IMPLICIT NONE

  REAL(8), intent(inout) :: pos(:,:)
  REAL(8), intent(in), value :: h
  REAL(8), intent(out) :: gradx(:,:), grady(:,:)

  REAL(8) :: xpos(maxn,maxn), ypos(maxn,maxn), rx(maxn,maxn), pos1(maxn), pos2(maxn), ry(maxn,maxn)
  REAL(8), PARAMETER :: pi = 4*ATAN(1.d0)
  INTEGER :: aa, bb, i, maxn

  maxn = 10000

  bb = blockDim%x * (blockIdx%x - 1) + threadIDx%x
  aa = blockDim%y * (blockIdx%y - 1) + threadIDx%x

  if ((aa .ge. maxn) .or. (bb .ge. maxn)) then !Stop running on hitting # of total threads

    return

  endif

  pos1(bb) = pos(1,bb)
  pos2(aa) = pos(2,aa)

  !xpos(:,bb) = pos1(:) - pos1(bb)
  !ypos(:,aa) = pos2(:) - pos2(aa)

  CALL syncthreads()

  do i = 1, maxn

    xpos(i,bb) = pos1(i) - pos1(bb)
    ypos(i,aa) = pos2(i) - pos2(aa)

    rx(bb,i) = xpos(bb,i)*xpos(bb,i) + ypos(bb,i)*ypos(bb,i) 
    ry(aa,i) = xpos(aa,i)*xpos(aa,i) + ypos(aa,i)*ypos(aa,i) 

    gradx(bb,i) = -exp(-0.5*rx(bb,i)*rx(bb,i)/(h*h))/(h*h*h*h*pi)*xpos(bb,i)
    grady(aa,i) = -exp(-0.5*ry(aa,i)*ry(aa,i)/(h*h))/(h*h*h*h*pi)*ypos(aa,i)
!  CALL syncthreads()
  enddo

!  gradx(aa,aa) = 1
!  grady(bb,bb) = 1
  !ry(:,:) = 1

END SUBROUTINE

!------------------------ END gweighf ---------------------------------

END MODULE
