! Authorship Info

! export PGI=/opt/pgi; export PATH=/opt/pgi/linux86-64/18.4/bin:$PATH; export LM_LICENSE_FILE=$LM_LICENSE_FILE:/opt/pgi/license.dat-COMMUNITY-18.4;

! To do:
!
! ---------- Take derivative of smoothing kernel (might switch to Gaussian)
!
! ---------- Figure out coordinates. 
!
! ---------- Fill star with initial particle soup
!
! ---------- Mess with units
!
! ---------- BOUNDARIES
!
! Figure out how to talk about PDEs in here (mixing from rotation? Momentum according to Jared)
!
! ---------- Implement large scale rigid rotation
!
! Implement differential rotation
!
! Parallelize some of these loops
!
! ---------- Implement gravity
!
! ---------- Artifical viscosity? Goes away with rigid rotation
!
!---------- Switch all the indices -_-
!
! Import things as vectors rather than individual coordinates

program sph_test

use cudafor
use parameters

IMPLICIT NONE

  REAL(8) :: c ! speed of sound in cm/s

  REAL(8) :: pos(1:2,1:maxn) = 0 ! Particle position r,theta
  REAL(8) :: vel(1:2,1:maxn) = 0 ! Particle velocity v_r, theta
  REAL(8) :: accel(1:2,1:maxn) = 0 ! Particle Acceleration r,theta (gravitational)
  REAL(8) :: m(1:maxn) ! Particle mass
  REAL(8) :: dist(1:maxn) ! Particle distance
  REAL(8) :: rho(1:maxn) ! Particle Density (local - scalar) g/cm^2
  REAL(8) :: P(1:maxn) ! Pressure at particle position (local - scalar) g/(cm*s^2)
  REAL(8) :: ut(1:maxn) = 0 ! Particle Energy (local - scalar)
  REAL(8) :: hsm ! Smoothing length 
  REAL(8) :: diffrot ! Differential Rotation Value

!  REAL(8) :: weighf, gweighf ! Function declarations

  REAL(8) :: dr, dtheta ! Initial particle spacing in polar

  INTEGER :: i, j, t, it_t
  REAL(8) :: lambda, k, nu, A, ut_old, start, finish, timing, rot
  REAL(8) :: ratioi, ratioj, vol, temp, dist2, vel0, avgP, avgrho
  REAL(8) :: grad(2,maxn,maxn)

  REAL(8), ALLOCATABLE :: tempo(:,:)

  REAL(8), ALLOCATABLE, device :: pos_d(:,:) !pos_d(2,maxn)
  REAL(8), ALLOCATABLE, device :: gradx_d(:,:), grady_d(:,:), ry_d(:,:) !gradx_d(maxn,maxn), grady_d(maxn,maxn), ry_d(maxn,maxn)
  type(dim3) :: grid, tblock

  CHARACTER(50) :: filename1, filename2, filename3

  k = 0.1
  nu = 10
  rot = 5
  A = 1!40

  dr = 1.d0/maxn 
  dtheta = 2.d0*pi/maxn ! unitless 

! Assuming intial quantities are evenly distributed among particles
  m(:) = 2.d0/maxn ! Normalized to 1 particle = 1 unit

  hsm = dt*sqrt(1000.d0/maxn) ! Mocz 2011
  lambda = 16.d0*k/(0.75*0.75*pi)

  CALL loadvalues(2,maxn,pos(:,:))

  do i = 1, maxn

    dist(i) = sqrt(pos(1,i)*pos(1,i)+pos(2,i)*pos(2,i))
    pos(1,i) = dist(i)*cos(i*dtheta) ! Unitless via parameters
    pos(2,i) = dist(i)*sin(i*dtheta) 

  enddo

  it_t = maxt/dt

  if (maxn == 1000) then

    filename1 = "./1000particles/pos_1000.dat"
    filename2 = "./1000particles/ut_1000.dat"
    filename3 = "./1000particles/dense_1000.dat"

  elseif (maxn == 10000) then

    filename1 = "./10000particles/pos_10000.dat"
    filename2 = "./10000particles/ut_10000.dat"
    filename3 = "./10000particles/dense_10000.dat"

  else

    filename1 = "./pos.dat"
    filename2 = "./ut.dat"
    filename3 = "./dense.dat"

  endif

  OPEN (unit=14, file=filename1, action="write",status="replace")
  OPEN (unit=15, file=filename2, action="write",status="replace")
  OPEN (unit=16, file=filename3, action="write",status="replace")

  CALL CPU_TIME(start)

  ALLOCATE(tempo(maxn,maxn))

do t = 1, it_t !leapfrog
  
  do i = 1, maxn

! Thanks https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec05.pdf
! Assuming no angular acceleration for everyone's sake.

! Differential rotation -- Eriguchi & Muller (1985)
    diffrot = rot*A*A/(A*A+dist(i)*dist(i))  
!    diffrot = rot + diffrot

! Updating positions
      pos(:,i) = pos(:,i) + vel(:,i)*dt + 0.5*accel(:,i)*dt*dt

! Velocities via forward differencing
      vel(:,i) = vel(:,i) + accel(:,i)*dt

    if (t > 1) then

! Updating positions with rotation consideration
     pos(1,i) = pos(1,i) + dist(i)*cos(diffrot*dt)
     pos(2,i) = pos(2,i) + dist(i)*sin(diffrot*dt)

! Velocity due to rigid rotation
      vel(1,i) = vel(1,i) - dist(i)*diffrot*sin(diffrot*dt) !pos(2,i)/dist(i)
      vel(2,i) = vel(2,i) + dist(i)*diffrot*cos(diffrot*dt) !rot*pos(1,i)/dist(i)

    endif

    !else

! Updating positions
     ! pos(:,i) = pos(:,i)+vel(:,i)*dt+0.5*accel(:,i)*dt*dt 

!Averaging the previous step and current to get the half time step (Mocz 2011)/wikipedia
      !vel(:,i) = 0.5*(vel(:,i)+dt*accel(:,i))

  dist(i) = sqrt(pos(1,i)*pos(1,i)+pos(2,i)*pos(2,i))

  enddo 

  do i = 1, maxn

! Differential rotation -- Eriguchi & Muller (1985)
    diffrot = rot*A*A/(A*A+dist(i)*dist(i))  

! Initializing Density

  rho(i) = m(i)*weighf(dist(i), 0.d0, hsm)
  if (rho(i) < 0.001) rho(i) = 2.d0/maxn

    do j = i+1, maxn

      temp = m(j)*weighf(dist(i),dist(j),hsm)

! Density at particle position 
      rho(i) = rho(i) + temp
      rho(j) = rho(j) + temp

    enddo

! Polytropic Pressure/EOS (Monaghan & Price 2004)
    P(i) = k*rho(i)*rho(i) 

! Lambda: Acceleration of particle due to gravity -- Monaghan & Price (2004)
    accel(:,i) = -nu*vel(:,i) - lambda*pos(:,i) 

! Can finite difference this to approximate derivative instead of just dividing by dt

    accel(1,i) = accel(1,i) - dist(i)*diffrot*diffrot*cos(diffrot*dt)
    accel(2,i) = accel(2,i) - dist(i)*diffrot*diffrot*sin(diffrot*dt) !+ diffrot/dt*asin(pos(2,i)/dist(i))

enddo

! PARALLELIZING -------------------------------- 

  ALLOCATE(pos_d(2,maxn))
  ALLOCATE(gradx_d(maxn,maxn))
  ALLOCATE(grady_d(maxn,maxn))

  tblock = dim3(1024,1,1)
  grid = dim3(ceiling(real(maxn)/real(tblock%x)),ceiling(real(maxn)/real(tblock%x)),1)

  pos_d(:,:) = pos(:,:)

  CALL gweighf<<<grid, tblock>>>(pos_d, hsm, gradx_d, grady_d)

  grad(1,:,:) = gradx_d(:,:)
  grad(2,:,:) = grady_d(:,:)

  DEALLOCATE(pos_d)
  DEALLOCATE(gradx_d)
  DEALLOCATE(grady_d)

  do i = 1,maxn

    ratioi = P(i)/(rho(i)*rho(i))
    
    do j = i+1, maxn

      ratioj = P(j)/(rho(j)*rho(j))

! Monaghan & Price (2004) 
! Acceleration due to pressure

      temp = m(j)*(ratioi + ratioj)

      accel(:,i) = accel(:,i) - temp*grad(:,i,j)
      accel(:,j) = accel(:,j) + temp*grad(:,i,j)

! Thermal Energy per unit mass/dt. 

      temp = temp*abs(dot_product(vel(:,i)-vel(:,j),grad(:,i,j)))

      ut(i) = ut(i) + temp
      ut(j) = ut(j) + temp
! Assume no thermal conduction between adjacent particles. Just another term that would have to be added.
    enddo

  enddo

  CALL avgdist(P(:), m(:), rho(:), dist(:), hsm, avgP)
  CALL avgdist(rho(:), m(:), rho(:), dist(:), hsm, avgrho)

!print*, "avgrho", avgrho
!print*, "avgP", avgP
  CALL causality(avgP, avgrho)

  do i = 1,maxn

    if (t > 1) pos(:,i) = pos(:,i)/maxval(dist(:))
    WRITE(14,100, advance='no') pos(1,i)
    WRITE(14,100) pos(2,i)

  enddo

  if (mod(t,10) == 0) print*, "Time Iteration", t, "Max Distance", maxval(dist(:))

  ut_old = sum(ut(:))

  dist(:) = dist(:)/maxval(dist(:)) !Re-normalizing to get properly bounded results
! Theoretically we could avoid this if I could find the 'force' that's causing this issue. I have
! been so far unable.
!  rho(:) = rho(:)/maxval(rho(:))

  WRITE(15,100) ut(:)/maxval(ut(:))
  WRITE(16,100) rho(:)/maxval(rho(:)) !ut(:)/maxval(ut(:))


enddo ! Time loop

  CALL CPU_TIME(finish)

  timing = finish-start

  print*
  print*, "Time to run: ", timing, "s"

  100  FORMAT (F16.8)

  CLOSE(14)
  CLOSE(15)
  CLOSE(16)

  DEALLOCATE(tempo)

END

! ---------------------- END MAIN PROGRAM -------------------------------
!
! ---------------------- SUBFUNCTIONS/SUBROUTINES -----------------------
! ---------------------- BEGIN avgdist ---------------------------------
!
! This subroutine calculates the average of a property over a distribution
! of particles.
!
! -----------------------

SUBROUTINE avgdist(param, m_i, rho_i, p1, hh, out_param) ! Monaghan (1992)

! ---------------------- Variables  
!
! param: {real - input}
!    Property of the system to be averaged over.
!
! m_i: {real - input}
!    Particle mass.
!
! rho_i: {real - input}
!    Particle mass.
!
! p1: {real - input}
!    Position 1 (see weighf).
!
! hh: {real - input}
!    Smoothing length.
!
! ----------------------

use parameters

IMPLICIT NONE

  REAL(8), intent(in) :: param(1:maxn), p1(1:maxn), hh, m_i(1:maxn), rho_i(1:maxn)
  INTEGER :: i
  REAL(8), intent(out) :: out_param

  out_param = 0.d0

  do i = 1, maxn

    out_param = out_param + weighf(p1(i), 0.d0, hh)*m_i(i)*param(i)/rho_i(i)

  enddo

END SUBROUTINE  

!---------------------------- END weighf ---------------------------------
!
! ---------------------- BEGIN causality ---------------------------------
!
! The subroutine causality checks if causality has been violated.
!
! -----------------------

SUBROUTINE causality(PP, Rho)

! ---------------------- Variables  
!
! n: {integer - input}
!    Problem (matrix and vector) size. 
!
! ----------------------

use parameters

IMPLICIT NONE

  REAL(8), INTENT(IN) :: PP, Rho
  REAL(8) :: c

! Causality check AKA is the speed of sound less than the speed of light?

  c = PP/Rho

  if (c > lightspeed) then

    print*, "CAUSALITY HAS BEEN BREACHED."
    print*, "Density: ", Rho, "Pressure: ", PP
    print*, "The model does not converge."

    STOP

  endif

end  

!------------------------ END Causality ---------------------------------
!
! ---------------------- BEGIN loadvalues -----------------------------
!
! The function loadvalues loads vector xx with values. 
!
! -----------------------

SUBROUTINE loadvalues(m, n, xx)

! ---------------------- Variables  
!
! m, n: {integer - input}
!    Size of matrix to be filled.
!
! xx(m,n): {float - output}
!    Matrix to be output of size (m,n) to be filled.    
!
! ----------------------

IMPLICIT NONE

  INTEGER, intent(in) :: m,n
  REAL(8), intent(out) :: xx(m,n)

! To get a random number every time, we have to restart the seed. Per
! the suggestion above, we have restarted the seed according to the 
! millisecond count of the system clock which should give a random
! number sequence that varies.

  INTEGER :: values(1:8), k
  INTEGER, ALLOCATABLE :: seed(:) 
 
  CALL DATE_AND_TIME(values = values)

  CALL RANDOM_SEED(size = k)
  ALLOCATE(seed(1:k))
  seed(:) = values(8)
  CALL random_seed(put = seed)

! AA and BB are generated by the intrinsic RANDOM_NUMBER routine.

  CALL RANDOM_NUMBER(xx)

end
! ---------------------- END loadvalues -------------------------------
!
! ---------------------- BEGIN AllocCheck -------------------------------
!
! AllocCheck checks whether or not memory allocation was successful.
!
! ----------------------

SUBROUTINE AllocCheck(Estate)

! ---------------------- Variable
!
! Estate: {integer}
!    Integer indicating allocation success or fail.
!
! ---------------------- 

IMPLICIT NONE

INTEGER, intent(in) :: Estate

! If the allocation is unsuccessful, then Estate will be greater than 0. 
! If this is so, we write to STDOUT to let the user know before stopping
! the program.

  if (Estate > 0) then

    WRITE(*,*) "Dynamic memory allocation failed, quitting."

    stop

  end if

END

! ---------------------- END AllocCheck ---------------------------------
!
! ---------------------- END SUBROUTINES/SUBFUNCTIONS -----------------  



