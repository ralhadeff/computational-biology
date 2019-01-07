! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ LD/MC code for implementation @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! this is naken code that runs 3D (in A coordinates) MC or LD
! the code is empty and is designed to be implemented by users

      program main

!	written by Raphael Alhadeff Jan 2018

      implicit none

!	conversion constant
      real, parameter :: bz_cu = 0.83146, kcal2cu = 418.3965
      double precision, parameter :: piValue = 3.14159265359
!	write frequency, i for local loops, type of run (0-MC,1-Overdamp LD,2-Underdamp LD)
!	p - particle index
      integer :: nParticles, logWriteFreq, method, p, i
!	total number of step and current step index
      integer(8) :: nsteps, istep
!	object for particle
!	position, velocity, force xyz 
!	sigma, mass, total energy 
!	note that for dimensionless systems you will have to split sigma, mass and friction to the different coordinates used (xyz in space can use the same values)
      type Particle
         double precision :: x, y, z, xV, yV, zV, xF, yF, zF, eTotal, sigma
         real :: mass
      end type
!	particles arrays
      type (Particle), dimension(:), allocatable :: particles
!       total energy, previous total energy, temperature(K), dt and friction if LD
!	current roll for MC
      double precision :: totalEnergy, previousEnergy, roll
      real :: temperature, dt, friction
!	input file name
      character(len=30) :: fileName
!	random number array (pregenerated for speed)
!	counter and size for rnd array
!	current seed for separate calls (without using array)
!	seed numbers to make different runs every execution
      double precision, dimension(:), allocatable :: rnd
      integer :: rndCounter, rndSize, currentSeed
      integer, dimension(1:12) :: seed
!       MC maximal step size (A) and energy barrier peak plateau width
      integer :: stepSize
!	rerun with specific seeds or not (using seeds file)
!	whether to print output files
!	whether to stop running (for early stopping)
      logical :: rerun, printOutput, kill
!	format for printing pdb's
      character(len=128) :: pdbFormat="('ATOM',2x,i5,2x,a4,a3,i6,4x,3f8.3)"
!##################################
! add all FF variables necessary ##
!##################################

! @@@@@@@@@@@@@@@@
! @@ read input @@
! @@@@@@@@@@@@@@@@
!       read input file name from the command line
      call get_command_argument(1,fileName)
      open(31,file=fileName)
      print *, "Reading in: """, trim(fileName), """"
      read(31,*) nsteps, logWriteFreq
!	store method in temp for translation
      read(31,*) temperature, temp
      if (temp.eq."MC") then
         method = 0
      elseif (temp.eq."OD") then
         method = 1
      else if (temp.eq."UD") then
         method = 2
      else
         stop "Method invalid, please specify MC, OD or UD"
      end if
      if (method.ne.0) then
         read(31,*) dt, friction
      end if
      if (method.eq.0) then
         text="Monte Carlo"
      elseif (method.eq.1) then
         text = "overdamped Langevin Dynamics"
      elseif (method.eq.2) then 
         text = "underdamped Langevin Dynamics"
      else 
         stop " Your method is invalid. Please choose 0, 1 or 2"
      end if
      print *, "Running ", trim(text)
      read(31,*) stepSize, rerun, printOutput
      read(31,*) nParticles
!	QA
      if (nParticles.lt.0) then 
         stop " Cannot use negative number of particles"
      end if
      allocate(particles(nParticles))
!######################
! add all FF readins ##
!######################
!	read in coordinates for the particles
!######################################################
! modify based on arguments needed for each particle ##
!######################################################
      do i = 0,nParticles
         if (method.eq.0) then
!	read MC input
            read(31,*) myosin(p)%x,myosin(p)%y,myosin(p)%z    
         else 
!	read LD input (with mass)
            read(31,*) myosin(p)%x,myosin(p)%y,myosin(p)%z, myosin(p)%mass
         end if
      end do
      close(31)

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ random number generator @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!	allocate rnd array
!	set maximum
      if (nsteps<= 10000000) then
         rndSize = nsteps;
      else
!       size of maximal array for random number
         rndSize = 10000000
      end if
!	5 because I need 1 number to choose the particle, 3 for xyz, and then the 5th for the roll
      allocate(rnd(rndSize*5))
      rndCounter = 0;
      if (rerun) then
         open(555,file='seeds')
         read(555,*) seed
         close(555)
      else 
!       open the urandom file (got this off the internet)
         open(89, file='/dev/urandom',access='stream' , form='UNFORMATTED')
!       get 12 seed integers (for the random_number subroutine seeds)
!       it takes 12 for some reason
         read(89) seed
         close(89)
!	write out seeds in case I want to rerun the same run again
         open(444,file='seeds',status='unknown')
         write(444,*) seed
         close(444)
      end if
!       feed the seeds
      call random_seed(put=seed)
!       generate random numbers in the array
      call random_number(rnd)

! @@@@@@@@@@@@@@@@@@@@@
! @@ print out files @@
! @@@@@@@@@@@@@@@@@@@@@
!	movie.pdb is the output movie
!	start is the first frame of the movie for visualization purposes (you can insert CONECT here)
      open(200,file='movie.pdb',status='unknown')
      open(201,file='start.pdb',status='unknown')
!	only for UD, print temperatures
      if (method.eq.2) then
         open(203,file='temperature.dat',status='unknown')
      end if
      do p=1,nParticles
!####################
! modify as needed ##
!####################
         write(201,pdbFormat) p,'A   ','RES',p,particles(p)%x,particles(p)%y,particles(p)%z
      end do
!#############################
! add connections, optional ##
!#############################
!         write(201,'(''CONECT'',i5,i5)') i,j
      close(201)

! @@@@@@@@@@@@@@@@@@@@
! @@ Initialization @@
! @@@@@@@@@@@@@@@@@@@@
!	set all velocities to 0
      particles%xV = 0
      particles%yV = 0
      particles%zV = 0
      if (method.eq.0) then
!	calculate and save previous energy
         call applyFF(particles)
         previousEnergy = totalEnergy
      end if
      if (printOutput) then
!#############################################
! initialize your output arrays or whatever ##
!#############################################
      end if
      if (method.ne.0) then
!	precalculate sigma's and get first seed
         call calculateSigma(particles)
!	seed has to be negative
         if (seed(1) .gt. 0) then
            currentSeed = -seed(1)
         else
            currentSeed = seed(1)
         end if
      end if

! @@@@@@@@@@@@@@@@@@@@@
! @@ simulation loop @@
! @@@@@@@@@@@@@@@@@@@@@
      do istep = 1,nsteps
!	apply method
         if (method.eq.0) then
            call doMC(particles)
         else
            call doLD(particles)
         end if
!       print frame to pdb using the user provided frequency
         if (mod(istep,logWriteFreq).eq.0) then
            write(200,'(''MODEL'',i12)') istep
            do p=1,nParticles
!####################
! modify as needed ##
!####################
               write(200,pdbFormat) p,'A   ','RES',p,particles(p)%x,particles(p)%y,particles(p)%z
            enddo
            write(200,'(''TER'')')
            write(200,'(''ENDMDL'',i12)') istep
            if (method.eq.2) then
               write(203,*) getTemperature(particles)
            end if
         endif
!	do these checks only every 1000 steps, to save on computation time
         if (mod(iStep,1000).eq.0) then
!	check whether there is a kill 'command'
!	the kill command is simply the existence of a 'kill' file in the folder
            inquire (file="kill",exist=kill)
            if(kill) then
!######################################################
! implement whatever you want done upon kill request ##
!######################################################
               print *, "A kill has been requested"
!	delete kill file to prevent next submittion from stopping immediately
               open(unit=54322,file="kill",status='old')
               close(unit=54322,status='delete')
!	exit loop and print outputs
               exit
            end if
!	check whether there is a midpoint histogram request 'command'
!	reuse kill
            inquire (file="histograms",exist=kill)
            if(kill) then
!       delete kill file to prevent next iteration from printing again and again
               open(unit=54322,file="histograms",status='old')
               close(unit=54322,status='delete')
!       print outputs
!##################################################################
! implement whatever you want done upon output-interrupt request ##
!##################################################################
               kill = .false.
            end if
         end if
!	end of simulation loop
      end do
!       close pdb files
      close(200)
      if (method.eq.2) then
         close(203)
      end if

      if (printOutput) then
!#####################################################
! implement whatever you want for output at the end ##
!#####################################################
!         call printOutPut
      end if

!	end of main routine
      contains
! @@@@@@@@@@@@@@@@@@@@@@@@@
! @@ End of main program @@
! @@@@@@@@@@@@@@@@@@@@@@@@@

! @@@@@@@@@@@@@@@@@@@@@@@
! @@ perform a LD step @@
! @@@@@@@@@@@@@@@@@@@@@@@
      subroutine doLD(particles)

      implicit none

      type(Particle), dimension(:), intent(inout) :: particles

!	reset forces and apply potential
      call applyFF(particles)
!	generate random forces
!	myosin
      do p = 1,nParticles
         if (method.eq.1) then
            particles(p)%x = particles(p)%x + particles(p)%sigma*gasdev(currentSeed)
            particles(p)%y = particles(p)%y + particles(p)%sigma*gasdev(currentSeed)
            particles(p)%z = particles(p)%z + particles(p)%sigma*gasdev(currentSeed)
         else 
            particles(p)%xF = particles(p)%xF + particles(p)%sigma*gasdev(currentSeed)
            particles(p)%yF = particles(p)%yF + particles(p)%sigma*gasdev(currentSeed)
            particles(p)%zF = particles(p)%zF + particles(p)%sigma*gasdev(currentSeed)
         end if
      end do 
!       do mechanics
      do p = 1,nParticles
         call doMechanics(particles(p))
      end do

      end subroutine

! @@@@@@@@@@@@@@@@@@
! @@ do mechanics @@
! @@@@@@@@@@@@@@@@@@
      subroutine doMechanics(body)
      
      implicit none
  
      type (Particle) :: body

!	note that the random forces/translations have been applied already
      if (method.eq.1) then
!	overdamp
!	translate
         body%x = body%x + body%xF*dt/(body%mass*friction)
         body%y = body%y + body%yF*dt/(body%mass*friction)
         body%z = body%z + body%zF*dt/(body%mass*friction)
      else
!	underdamp
!	add drag force
         body%xF = body%xF - friction*body%mass*body%xV
         body%yF = body%yF - friction*body%mass*body%yV
         body%zF = body%zF - friction*body%mass*body%zV
!	translate
         body%x  = body%x  + body%xV*dt+0.5d0*(body%xF/body%mass)*dt**2
         body%y  = body%y  + body%yV*dt+0.5d0*(body%yF/body%mass)*dt**2
         body%z  = body%z  + body%zV*dt+0.5d0*(body%zF/body%mass)*dt**2
!	update velocity
         body%xV = body%xV + (body%xF/body%mass)*dt
         body%yV = body%yV + (body%yF/body%mass)*dt
         body%zV = body%zV + (body%zF/body%mass)*dt
      end if

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@
! @@ perform a MC step @@
! @@@@@@@@@@@@@@@@@@@@@@@
      subroutine doMC(particles)


      type(Particle), dimension(:), intent(inout) :: particles

!       x,y,z_temp is the vector coordinates of the MC step
!       chance is the boltzmann probability for the metropolis criteria
!	roll is a global variable because it is sometimes used later in checkATP
      double precision :: xTemp, yTemp, zTemp, chance

!       get a random particle
!       rnd has 6*rndSize elements
!       for each step it needs a particle number, x,y,z and the metropolis number (in that order)
      i = floor(rnd((rndCounter)*5)*(nParticles))+1
!	from -stepSize/2 to +stepSize/2
!	note that this results in steps that are 0 to 0.87*stepSize in 3D
      xTemp = (rnd((rndCounter)*6+1)-0.5d0)*stepSize
      yTemp = (rnd((rndCounter)*6+2)-0.5d0)*stepSize
      zTemp = (rnd((rndCounter)*6+3)-0.5d0)*stepSize
      roll   = rnd((rndCounter)*6+4)
!	move the particle
      call moveMC(particles,i,xTemp,yTemp,zTemp)
!       get energy after change
      call applyFF(particles)
!       calculate chance to pass (metropolis)
      chance = exp((totalEnergy-previousEnergy)/(-0.002d0*temperature))
!       if accepted:
!       save the previous energy, and update boolean
      if (chance.gt.roll) then
         previousEnergy = totalEnergy
      else
!       if rejected revert the change
         call moveMC(particles,i,-xTemp,-yTemp,-zTemp)
      end if
!       generate new rnd values if needed (only relevent for very long runs)
      rndCounter=rndCounter+1
      if (rndCounter>=rndSize) then
         call random_number(rnd)
         rndCounter=0
      end if

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ subroutine to preform the MC moving @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine moveMC(particles,i,x,y,z)
      
      implicit none

      type(Particle), dimension(:), intent(inout) :: particles
!	i is the index of the particle.
      integer :: i, a, b
      double precision :: x,y,z

!	move particle 
      particles(i)%x=particles(i)%x+x
      particles(i)%y=particles(i)%y+y
      particles(i)%z=particles(i)%z+z

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ subroutine to calculate the total energy of the system @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine applyFF(particles)

      implicit none

      type(Particle), dimension(:), intent(inout) :: particles

!	have to use a private i because global i might be in use
      integer :: i, j 

      if (method.eq.0) then
!	Set all energies to 0
         totalEnergy = 0.d0
         particle%eTotal = 0.d0
      else
!	set forces to 0
         particles%xF = 0.d0
         particles%yF = 0.d0
         particles%zF = 0.d0
      end if

! @@@@@@@@@@@@@@
! @@ bonds    @@
! @@ angles   @@
! @@ torsions @@
! @@@@@@@@@@@@@@
!#################################
! optional - pre written methods##
!#################################
!         call getBondEnergy(kb,r0,particleI,particleJ)
!         call getAngleEnergy(ka,aMin,aMax,pI,pJ,pK)
!         call getTorsion(kt,a0,I,J,K,L)
!         call getRejection(kb,rMin,I,J) ; harmonic rejection below rMin

!	Accumulate the total energy
      if (method.eq.0) then
         do p = 1,nParticles
            totalEnergy = totalEnergy + particles(p)%eTotal
         end do
      end if

!######################################
! implement all of your FF here      ##
! this is the core of the simulation ##
!######################################

      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Calculate the energy of marcus parbolas @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double precision function markus(x, startX, endX, startE, barrierE, endE, gap,flatten)
!	modified markus with different skewness, and no h12

!	startX and endX are the points in real space where the potential should begin and end
!	startE and endE are the energies at the beginning and ending
!	barrier is the energy at the peak
!	gap is the length of the gap (the plateau at the peak)

      implicit none

      double precision :: f, g, t, barrier, a1, a2, output
      double precision, intent(in) :: x, barrierE
      real, intent(in) :: startX, startE, endX, endE, gap
!	whether to flatten beyond the range or not?
      logical, intent(in) :: flatten

!       NOTE - the barrier HAS to be higher than both startE and endE
!	otherwise - barrier will be set to the higher of the two
      if (endE.gt.startE .AND. barrierE.lt.endE) then 
         barrier = endE
      else if (startE.gt.endE .AND. barrierE.lt.startE) then
         barrier = startE
      else
         barrier = barrierE
      end if

      a1 = 4*(barrier - startE)/(endX-startX-gap)**2
      a2 = 4*(barrier - endE)  /(startX-endX+gap)**2

      t = x
      if (method.eq.0) then
!	energy
!	if flatten is true, then set to edges' energy at the beginning and end
!	otherwise, leave as is	
         if (t.lt.startX .AND. flatten) then
            output = startE
         else if (t.gt.(startX+endX-gap)/2 .AND. t.lt.(startX+endX+gap)/2) then
            output = barrier
         else if (t.ge.endX .AND. flatten) then
            output = endE
         else
!	calculate the energy in the scaled coordinates
!	energy is NOT to be scaled, return as is
            if (t.lt.(startX+endX)/2) then
               output = a1*(t-startX)**2+startE
            else
               output = a2*(t-endX)**2+endE
            end if
         end if
      else 
!	force (using a scale function rather than the markus
!	treat the truncations
         if (t.lt.startX .AND. flatten) then
            output = 0
         else if (t.gt.(startX+endX-gap)/2 .AND. t.lt.(startX+endX+gap)/2) then
            output = 0
         else if (t.ge.endX .AND. flatten) then
            output = 0
         else if (t.lt.(startX+endX)/2) then
            output = -2*a1*(t-startX)
         else
            output = -2*a2*(t-endX)
         end if
         output = output*kcal2cu
      end if

      markus = output

      end function markus     

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Calculate the energy of marcus parbolas @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double precision function footAngleEnergy(theta, footDown, footUp, k, fDEnergy, barrierE, fUEnergy, gap)
!	modified of the modified markus for the foot angle
!	the reason I need a different function is that the edges should be more wall-like
!	IMPORTANT NOTE:	
!	for reasons of ease, the lever down angle is negative and the lever up is positive
!	lever down is going forward (so tailing foot) 

      implicit none

!	k for harmonic constraint at the edges
      real, intent(in) :: footDown, footUp, fDEnergy, fUEnergy, barrierE, gap, k
      double precision, intent(in) :: theta
      double precision :: barrier, a1, a2, output

!       NOTE - the barrier HAS to be higher than both energies
!       otherwise - barrier will be set to the higher of the two
      if (fDEnergy.gt.fUEnergy .AND. barrierE.lt.fDEnergy) then
         barrier = fDEnergy
      else if (fUEnergy.gt.fDEnergy .AND. barrierE.lt.fUEnergy) then
         barrier = fUEnergy
      else
         barrier = barrierE
      end if

      a1 = 4*(barrier - fDEnergy)/(footUp-footDown-gap)**2
      a2 = 4*(barrier - fUEnergy)/(footDown-footUp+gap)**2

      if (method.eq.0) then
!       energy
         if (theta.lt.footDown) then
            output = fDEnergy+k*(theta-footDown)**2
         else if (theta.gt.footUp) then
            output = fUEnergy+k*(theta-footUp)**2
         else if (theta.gt.(footDown+footUp-gap)/2 .AND. theta.lt.(footDown+footUp+gap)/2) then
            output = barrier
         else if (theta.lt.(footDown+footUp)/2) then 
            output = a1*(theta-footDown)**2+fDEnergy
         else
            output = a2*(theta-footUp)**2+fUEnergy
         end if
      else
         if (theta.lt.footDown) then
            output = 2*k*(theta-footDown)
         else if (theta.gt.footUp) then
            output = 2*k*(theta-footUp)
         else if (theta.gt.(footDown+footUp-gap)/2 .AND. theta.lt.(footDown+footUp+gap)/2) then
            output = 0
         else if (theta.lt.(footDown+footUp)/2) then
            output = 2*a1*(theta-footDown)
         else
            output = 2*a2*(theta-footUp)
         end if
         output = output*kcal2cu
      end if

      footAngleEnergy = output

      end function footAngleEnergy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Get rejection energy (half harmonic) @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getRejection(kb,rMin,pi,pj)

      implicit none
!       bond constant and equilibrium length
      real, intent(in) ::  kb, rMin
!       particles i and j
      type (Particle), intent(inout) :: pi, pj
!       local helper variables
!       vector i to j, length and base force (to avoid recalculating)
      double precision, dimension(:) :: vec(3)
      double precision :: dr, f

      call vector(pi,pj,vec)
      dr = length(vec)
      if (dr.lt.rMin) then
         if (method.ne.0) then
            f = -2.d0 * kb * (rMin-dr)/dr
            f = f * kcal2cu
!       add force to pi
            call addForce(pi,f*vec)
!       add force (negative because of dxyz direction) to pj
            call addForce(pj,-f*vec)
         else
!       add energies            
            pi%eTotal = pi%eTotal + kb * (dr-rMin)**2
         end if
      end if

      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Get bond energy and force (harmonic) @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getBondEnergy(kb,r0,pi,pj)

      implicit none
!	bond constant and equilibrium length
      real, intent(in) ::  kb, r0
!	particles i and j
      type (Particle), intent(inout) :: pi, pj
!	local helper variables
      double precision :: dr, f
      double precision, dimension(:) :: vec(3)

      call vector(pi,pj,vec)
      dr = length(vec)
      if (method.ne.0) then
         f = -2.d0 * kb * (r0-dr)/dr 
         f = f * kcal2cu
!	add force to pi
         call addForce(pi,f*vec)
!	add force (negative because of dxyz direction) to pj
         call addForce(pj,-f*vec)
      else 
!	add energies            
         pi%eTotal = pi%eTotal + kb * (dr-r0)**2
      end if

      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Get angle energy (harmonic) @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getAngleEnergy(ka,thetaMin,thetaMax,pi,pj,pk)

      implicit none
!	angle constant and equilibrium theta (radians)
      real, intent(in) :: ka, thetaMin, thetaMax
!	particles i, j and k
!	j is the 'hinge' particle
      type (Particle), intent(inout) :: pi, pj, pk
!	local variables
!	different delta x,y,z for the vector calculations
!	vector sizes
!	cosine of Theta, theta and the base force
!	st is sqrt(1-cos**2) - derivative of acos
!	aa and bb are just precalculated stuff
      double precision :: rij,rkj
      double precision, dimension(3) :: vecIJ, vecKJ, forceI, forceK
      double precision :: cosTheta, theta, f, sinTheta, theta0
!       tempF are the final forces
      double precision :: tempF(6)

      call vector(pi,pj,vecIJ)
      rij = length(vecIJ)
      call vector(pk,pj,vecKJ)
      rkj = length(vecKJ)
      cosTheta = dot(vecIJ,vecKJ)/(rij*rkj)
!	in case of machine precision overload
      if (cosTheta.gt.1.0) then
         theta = 0.d0
      else if (cosTheta.lt.-1.0) then 
         theta = piValue
      else 
         theta = acos(cosTheta)
      end if
!	only proceed if theta is within requested range
      if (theta.gt.thetaMax.OR.theta.le.thetaMin) then
!	update reference based on which side of the potential theta is
         if (theta.le.thetaMin) then
            theta0 = thetaMin
         else if (theta.gt.thetaMax) then
            theta0 = thetaMax
         else
            theta0 = theta
         end if
         if (method.ne.0) then
!	do the derivations
            f = -2*ka*(theta-theta0)*kcal2cu
            sinTheta = sqrt(1-cosTheta**2)
!	if forces will overflow skip (force = 0)
            if (sinTheta.gt.1.0d-8) then
!	force on pi
               forceI = (-f/sinTheta)*(-vecKJ+(rkj/rij)*cosTheta*vecIJ)/(rij*rkj)
               call addForce(pi, forceI)
!	force on pk
               forceK = (-f/sinTheta)*(-vecIJ+(rij/rkj)*cosTheta*vecKJ)/(rij*rkj)
               call addForce(pk, forceK)
!	force on pj
               call addForce(pj, -forceI-forceK)
            end if
         else
!	convert to energy and pass it to particle j
            pj%eTotal=pj%eTotal+(ka*(theta - theta0)**2)
         end if
      end if

      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Get distance between two particles @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double precision function distance(pi,pj)

      implicit none
!	particles i and j
      type (Particle), intent(in) :: pi,pj
      double precision, dimension(:) :: vec(3)
      call vector(pj,pi,vec)
      distance = length(vec)
      end function distance

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ align normal of a plane (3 particles) to a given vector @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine alignPlane(k,theta0,pi,pj,pk,ref,angle)

      implicit none

      type (Particle), intent(inout) :: pi,pj,pk
      double precision, dimension(3), intent(in) :: ref
      real, intent(in) :: theta0
      double precision, intent(in) :: k
      double precision, intent(out) :: angle
!       local variables
!	theta is the angle between n (normal of ijk plane) and m
      double precision :: rm, rn, cosTheta, theta, sinTheta, f
      double precision, dimension(:) :: VecJI(3), vecJK(3), n(3), vecKI(3), fI(3),fJ(3),fK(3), helper1(3),helper2(3),force(3)

      call vector(pj,pi,vecJI)
      call vector(pj,pk,vecJK)
      call vector(pk,pi,vecKI)
      call cross(vecJI,vecJK,n)
      rm = length(ref)
      rn = length(n)
      cosTheta = dot(n,ref)/(rm*rn)
      theta = acos(cosTheta)
      angle = theta
      if (method.eq.0) then 
         pj%eTotal = pj%eTotal+k*(theta-theta0)**2
      else 
         sinTheta = sqrt(1-cosTheta**2)
         if (sinTheta.gt.1.0d-8) then
            f = 2*k*((theta-theta0)/(sinTheta*rm*rn))*kcal2cu
!	force on i
            call cross(vecJK,ref,helper1)
            call cross(vecJK,n,helper2)
            force = -f*(helper1-(rm/rn)*cosTheta*helper2)
            call addForce(pi, force)
!	force on j
            call cross(ref,vecKI,helper1)
            call cross(n,vecKI,helper2)
            force = f*(helper1-(rm/rn)*cosTheta*helper2)
            call addForce(pj, force)
!	force on k
            call cross(vecJI,ref,helper1)
            call cross(vecJI,n,helper2)
            force = f*(helper1-(rm/rn)*cosTheta*helper2)
            call addForce(pk, force)
         end if
      end if

      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Apply forces on four particles around a torsion angle @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getTorsion(kt, phi0, pi, pj, pk, pl, iType, iParameter)

      implicit none

      real, intent(in) :: kt, phi0
      type (Particle), intent(inout) :: pi, pj, pk, pl
!	iType is the type of call (get angle, input force etc ) and iParameter is a case specific argument
      integer, optional, intent(in) :: iType
!	types:
!	1 - return angle and exit
!	2 - input parameter into f and apply forces
!	3 - apply a periodic force (1+cos(nPhi+Phi0)
!	    in this case, iParameter is equal to 10*n+phi0 (is explicit since n is integer and phi0 <= 2*pi <10)
      double precision, optional, intent(inout) :: iParameter
!       vectors ij, jk etc. cross products (n1 and n2)
      double precision, dimension(3) :: vecJI, vecJK, vecLK, jixjk, kjxkl, n1xn2
!	lengths, angles, direction is a helper dot product to determine the angle's characteristics
      double precision :: n1, n2, cosPhi, phi, direction, dPhi
!	derivative and helper stuff
      double precision :: sinPhi
      double precision,dimension(3) :: temp1, temp2, force1, force2, force3, force4
!	force/energy
      double precision :: f, phase
      integer :: i, n

      call vector(pj,pi, vecJI)
      call vector(pj,pk, vecJK)
      call vector(pl,pk, vecLK)
!	ji cross jk
      call cross(vecJI,vecJK,jixjk)
!	lk cross jk is also kj cross kl
      call cross(vecLK,vecJK, kjxkl)
      n1 = length(jixjk)
      n2 = length(kjxkl)
      if (n1.gt.1.d-8 .AND. n2.gt.1.d-8) then
         cosPhi = dot(jixjk,kjxkl)/(n1*n2)
!	I need this condition because I was having overlow of cosPhi being too close to +-1
         if (cosPhi.gt.1.0) then
            phi = 0
         else if (cosPhi.lt.-1.0) then
            phi = piValue
         else
            phi = acos(cosPhi)
         end if
         call cross(jixjk,kjxkl,n1xn2)
!	convension stuff to get the correct angle range and sign
         direction = dot(vecJK,n1xn2)
         if (direction.gt.0) then
            phi = -phi
         end if
!	phi was the angle between the normals, we're interested in the angle between the planes
         phi = piValue -phi
         if (phi.gt.piValue) then
            phi = phi - 2*piValue
         end if
!	output angle and quit if requested
         if (present(iType) .AND. iType.eq.1) then
            iParameter = phi
         else
!	get delta Phi
            dPhi = phi0 - phi
!	normalize dPhi to select shortest path to phi0
            if (abs(dPhi).gt.piValue) then
               dPhi = sign(2*piValue-abs(dPhi),dPhi)
            end if
!	get energy/forces
            if (method.eq.0) then
               if (present(iType)) then
                  if (iType.eq.2) then
                     f = iParameter
                  end if
               else
                  f = kt*dPhi**2
               end if
!	split energy between the two particles of the bond itself
               pj%eTotal=pj%eTotal+f/2
               pk%eTotal=pk%eTotal+f/2
            else
               if (present(iType)) then
                  if (iType.eq.2) then
                     f = iParameter
                  end if
               else
                  f = -2*kt*dPhi*kcal2cu
               end if
!	apply force based on partial derivatives
!	derivative of aCos
               sinPhi = sqrt(1-cosPhi**2)
               if (direction.gt.0) then
                  sinPhi = -sinPhi
               end if
!	if dAcosdPhi is too small, the force will overflow (and should be zero)
               if (abs(sinPhi).gt.1.d-8) then
                  temp1 = -1*((n1/n2)*kjxkl-cosPhi*jixjk)/(n1**2*sinPhi)
                  temp2 = -1*((n2/n1)*jixjk-cosPhi*kjxkl)/(n2**2*sinPhi)
                  call cross(vecJK,temp1,force1)
                  call cross(vecJI,temp1,force2)
                  call cross(temp2,vecLK,force3)
                  call cross(temp2,vecJK,force4)
!	actualy forces
                  call addForce(pi, f*force1)
                  call addForce(pj, f*(-force1+force2-force3))
                  call addForce(pk, f*(-force2-force4+force3))
                  call addForce(pl, f*force4)
               end if
            end if
         end if
      end if
      
      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Add vector to force values of particle @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine addForce(p, forceV)

      implicit none

      type (Particle), intent(inout) :: p
      double precision, dimension(:), intent(in) :: forceV(3)

      p%xF = p%xF+forceV(1)
      p%yF = p%yF+forceV(2)
      p%zF = p%zF+forceV(3)

      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Get vector between two particles @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine vector(pi, pj, output)

      implicit none

      type (Particle), intent(in) :: pi, pj
      double precision, dimension(:), intent(out) :: output(3)

      output(1) = pj%x-pi%x
      output(2) = pj%y-pi%y
      output(3) = pj%z-pi%z

      end subroutine

!@@@@@@@@@@@@@@@@@
!@@ Vector size @@
!@@@@@@@@@@@@@@@@@
      double precision function length(vec)

      implicit none

      double precision, dimension(:), intent(in) :: vec

      length = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)

      end function

!@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Vector dot product @@
!@@@@@@@@@@@@@@@@@@@@@@@@
      double precision function dot(vec1, vec2)

      implicit none

      double precision, dimension(:), intent(in) :: vec1, vec2

      dot = vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)

      end function

!@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Vector cross product @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine cross(vec1, vec2, output)

      implicit none

      double precision, dimension(:), intent(in) :: vec1, vec2
      double precision, dimension(:), intent(out) :: output(3)

      output(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
      output(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
      output(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)

      end subroutine

!@@@@@@@@@@@@@@@@@
!@@ Vector size @@
!@@@@@@@@@@@@@@@@@
      double precision function getTemperature(myosin,adp,pi)

      implicit none

      type(Particle), dimension(:), intent(in) :: myosin, adp, pi
      integer :: i
      double precision :: totalK

      totalK = 0.d0

!	sum the kinetic energy of all particles
      do i = 1, pLast
         totalK = totalK + 0.5d0*myosin(i)%mass*(myosin(i)%xV**2+myosin(i)%yV**2+myosin(i)%zV**2)
      end do
      do i = 1, 4
         totalK = totalK + 0.5d0*adp(i)%mass*(adp(i)%xV**2+adp(i)%yV**2+adp(i)%zV**2)
         totalK = totalK + 0.5d0*pi(i)%mass*(pi(i)%xV**2+pi(i)%yV**2+pi(i)%zV**2)
      end do

!	divide by Bolzmann and n of particles (pLast + 4 adp + 4 pi)
      getTemperature = totalK/(bz_cu*pLast+8)

      end function

!@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Precalculate sigma's @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine calculateSigma(particles)

      implicit none

      integer :: i

      type(Particle), dimension(:), intent(inout) :: particles

!       particles
      if (method.eq.1) then
         particles%sigma = sqrt(2*bz_cu*temperature*dt/(particles%mass*friction))
      else
         particles%sigma = sqrt(2*bz_cu*temperature*particles%mass*friction/dt)
      end if

      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Generate random float from 0 to 1 @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double precision function randomFloat(seed)

      implicit none

      integer seed, IA, IM, IQ, IR, NTAB, NDIV
      double precision EPS, RNMX

      parameter (IA=16807, IM=2147483647, IQ=127773, IR=2836, NTAB=32, NDIV=1+(IM-1)/NTAB, EPS=1.2e-7, RNMX=1.-EPS)

      integer j,k,iv(NTAB),iy
      save iv,iy
      DATA iv /NTAB*0/, iy /0/

      if (seed.le.0.or.iy.eq.0)  then
         seed=max(-seed,1)
         do j=NTAB+8, 1, -1
            k=seed/IQ
            seed=IA*(seed-k*IQ)-IR*k
            if (seed.lt.0) seed=seed+IM
            if (j.le.NTAB) iv(j)=seed
         end do
         iy=iv(1)
      end if

      k=seed/IQ
      seed=IA*(seed-K*IQ) - IR*k

      if (seed.lt.0) seed=seed+IM

      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=seed

      randomFloat=min((1.d0/IM)*iy, RNMX)

      end function randomFloat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Generate random number with normal distribution @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double precision function gasdev(seed)

      implicit none

      integer seed
      integer iset
      double precision fac,gset,rsq,v1,v2, output
      logical loop
      save iset, gset
      data iset/0/

      loop = .true.
      if (seed.lt.0) iset = 0
      if (iset.eq.0) then
         do while (loop)
            v1 = 2.d0*randomFloat(seed)-1
            v2 = 2.d0*randomFloat(seed)-1
            rsq = v1**2 + v2**2
            if (rsq.ge.1..or.rsq.eq.0.) then
               loop = .true.
            else
               loop = .false.
            end if
         end do
         fac = sqrt(-2.d0*log(rsq)/rsq)
         gset = v1*fac
         output = v2*fac
         iset = 1
      else
         output = gset
         iset = 0
      endif

      gasdev = output
      end function

      end program
