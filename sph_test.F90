! Matthew Portman
! 5/10/18
! sph_test.F90
! --------------------------------------------------------------------------------
!
! This code simulates a toy star by iterating forward in time and with Smoothed 
! Particle Hydrodynamics (SPH) to model particle interactions.
!
! ---------------------------------------------------------------------------------

program sph_test

use parameters_test

IMPLICIT NONE

  REAL(8) :: c ! speed of sound in cm/s
  REAL(8) :: ang_vel = 40 ! revs/s, downscaled to appropro time step

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

!  REAL(8) :: weighf, gweighff ! Function declarations

  REAL(8) :: dr, dtheta ! Initial particle spacing in polar

  INTEGER :: i, j, t, it_t
  REAL(8) :: lambda, k, nu, A, ut_old, start, finish, timing, rot
  REAL(8) :: ratioi, ratioj, vol, temp, dist2, vel0, avgP, avgrho
  REAL(8) :: grad(1:2,1)

  CHARACTER(30) :: filename1, filename2, filename3

  k = 0.1
  nu = 10
  rot = 5
  A = 1!40

!  star_moment = 0.5*star_mass*star_rad*star_rad !g*cm^2
!  star_E = 0.5*I*ang_vel*ang_vel ! Rotational KE g*cm^2*s^-2

  dr = 1.d0/maxn ! cm -> unitless (units of star_radius)... 'units' taken from mocz
  dtheta = 2.d0*pi/maxn ! unitless 

! Assuming intial quantities are evenly distributed among particles
  m(:) = 2.d0/maxn ! Normalized to 1 particle = 1 unit

  hsm = dt*sqrt(1000.d0/maxn) ! Mocz 2011
  lambda = 16.d0*k/(0.75*0.75*pi)

  ! Loading in random positions for all particles
  CALL loadvalues(2,maxn,pos(:,:))

  do i = 1, maxn

  ! Caclulating starting distances and places the particles at angles.
    dist(i) = sqrt(pos(1,i)*pos(1,i)+pos(2,i)*pos(2,i))
    pos(1,i) = dist(i)*cos(i*dtheta) 
    pos(2,i) = dist(i)*sin(i*dtheta) 

  enddo

  it_t = maxt/dt

  ! Print naming
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

do t = 1, it_t 
  
  do i = 1, maxn

! Thanks https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec05.pdf
! Assuming no angular acceleration for everyone's sake.

! Differential rotation -- Eriguchi & Muller (1985)
    diffrot = rot*A*A/(A*A+dist(i)*dist(i))  

! Updating positions
      pos(:,i) = pos(:,i) + vel(:,i)*dt + 0.5*accel(:,i)*dt*dt

! Velocities via forward differencing
      vel(:,i) = vel(:,i) + accel(:,i)*dt

    if (t > 1) then

! Updating positions with rotation consideration
     pos(1,i) = pos(1,i) + dist(i)*cos(diffrot*dt)
     pos(2,i) = pos(2,i) + dist(i)*sin(diffrot*dt)

! Velocity due to rigid rotation
      vel(1,i) = vel(1,i) - dist(i)*diffrot*sin(diffrot*dt) 
      vel(2,i) = vel(2,i) + dist(i)*diffrot*cos(diffrot*dt) 

    endif

    !else ! Leapfrog... Not sure if my implementation was right so I'll keep it commented out.

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

! Including acceleration due to rotation
    accel(1,i) = accel(1,i) - dist(i)*diffrot*diffrot*cos(diffrot*dt)
    accel(2,i) = accel(2,i) - dist(i)*diffrot*diffrot*sin(diffrot*dt)

enddo

  do i = 1,maxn

    ratioi = P(i)/(rho(i)*rho(i))
    
    do j = i+1, maxn

      ratioj = P(j)/(rho(j)*rho(j))
      temp = -m(j)*(ratioi + ratioj)

      CALL gweighff(pos(1,i), pos(2,i), pos(1,j), pos(2,j), dist(i), dist(j), hsm, grad(:,1))
	  
! Monaghan & Price (2004) 
      accel(:,i) = accel(:,i) + temp*grad(:,1)
      accel(:,j) = accel(:,j) - temp*grad(:,1)

! Assume no thermal conduction between adjacent particles. Just another term that would have to be added.
! Thermal Energy per unit mass/dt. 
      temp = -temp*abs(dot_product(vel(:,i)-vel(:,j),grad(:,1)))
      ut(i) = ut(i) + temp 
      ut(j) = ut(j) + temp

    enddo

  enddo

! Checking causality
  CALL avgdist(P(:), m(:), rho(:), dist(:), hsm, avgP)
  CALL avgdist(rho(:), m(:), rho(:), dist(:), hsm, avgrho)

  CALL causality(avgP, avgrho)

! Printing. Here's effectively where I use my Boundary Conditions to 'keep' particles inside the star.
! It's a little bit crude to re-normalize but until I can fine tune what is driving them out, it'll do.

  do i = 1,maxn

    if (t > 1) pos(:,i) = pos(:,i)/maxval(dist(:))
    WRITE(14,100, advance='no') pos(1,i)
    WRITE(14,100) pos(2,i)

  enddo

  if (mod(t,10) == 0) print*, "Time Iteration", t, "Max dist", maxval(dist(:))
  ut_old = sum(ut(:))

  dist(:) = dist(:)/maxval(dist(:)) 

  WRITE(15,100) ut(:)/maxval(ut(:))
  WRITE(16,100) rho(:)/maxval(rho(:))


enddo ! Time loop

  100  FORMAT (F16.8)

  CLOSE(14)
  CLOSE(15)
  CLOSE(16)

  CALL CPU_TIME(finish)

  timing = finish - start

  print*
  print*, "Timing: ", timing, "s"

END

! ---------------------- END MAIN PROGRAM -------------------------------
!
! ---------------------- SUBFUNCTIONS/SUBROUTINES -----------------------
!
!------------------------ BEGIN gweighff ---------------------------------
!
! The gradient of the weighting function/kernel.
!
! -----------------------


SUBROUTINE gweighff(xpos1, ypos1, xpos2, ypos2, r1, r2, h, grad) ! Monaghan (1992)

! ---------------------- Variables  
!
! x/ypos1, x/ypos2: {real - input}
!    Coord positions of particles 1 and 2 
!
! h: {real - input}
!    Smoothing length.
!
! ----------------------

use parameters_test

IMPLICIT NONE

  REAL(8), intent(in) :: xpos1, ypos1, xpos2, ypos2, r1, r2, h
  REAL(8), intent(out) :: grad(1:2,1)

  REAL(8) :: temp, xpos, ypos 

  ! Difference in positions
  xpos = xpos1-xpos2
  ypos = ypos1-ypos2

  ! Analytical solution to gaussian kernel
  temp = -weighf(r1,r2,h)/(h*h)

  grad(1,1) = temp*xpos
  grad(2,1) = temp*ypos

END SUBROUTINE

!------------------------ END gweighff ---------------------------------
!
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

use parameters_test

IMPLICIT NONE

  REAL(8), intent(in) :: param(1:maxn), p1(1:maxn), hh, m_i(1:maxn), rho_i(1:maxn)
  INTEGER :: i
  REAL(8), intent(out) :: out_param

  out_param = 0.d0

  ! Averaging over the system property
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

use parameters_test

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
! Didn't end up allocating any arrays but I always keep this around just 
! in case. 
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



