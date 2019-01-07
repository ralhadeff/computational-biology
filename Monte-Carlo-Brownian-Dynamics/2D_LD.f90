! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ LD/MC code for adenylate kinase renormalization @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      program main

!	written by Raphael Alhadeff March 2018

      implicit none

!	conversion constant
      real, parameter :: bz_cu = 0.83146, kcal2cu = 418.3965
!	write frequency, i for local loops, type of run (0-MC,1-Overdamp LD,2-Underdamp LD)
!	p - particle index
      integer :: nParticles, logWriteFreq, method, p, i, j
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
      real :: temperature, dt
      real, dimension(2) :: friction
!	input file name
      character(len=30) :: fileName, runType, text
!	random number array (pregenerated for speed)
!	counter and size for rnd array
!	current seed for separate calls (without using array)
!	seed numbers to make different runs every execution
      double precision, dimension(:), allocatable :: rnd
      integer :: rndCounter, rndSize, currentSeed
      integer, dimension(1:12) :: seed
!       MC maximal step size (A) and energy barrier peak plateau width
      real :: stepSize
!	rerun with specific seeds or not (using seeds file)
!	whether to print output files
!	whether to stop running (for early stopping)
      logical :: rerun, printOutput, kill
!	format for printing pdb's
      character(len=128) :: pdbFormat="('ATOM',2x,i5,2x,a4,a3,i6,4x,3f8.3)"
!	FF variables
      real :: xStart, xEnd, xExtra, yStart, yEnd, yExtra
      real :: boundaryK
!	modified markus energy at start, end and barrier
      real :: x0, x1, xb
      real :: w0, w1, wb
      real :: y0, y1, yb
!	histogram array
      integer, dimension(:,:), allocatable :: histogram
!	free energy, potential energy and the gap
      real :: eneF, eneP, eneG
!	max gap for missing points
      real :: eneMax
      integer :: eneCount
!	state variables
      integer :: state, iState, destination, reason
      real :: threshold
! @@@@@@@@@@@@@@@@
! @@ read input @@
! @@@@@@@@@@@@@@@@
!       read input file name from the command line
      call get_command_argument(1,fileName)
      print *, "You are running ADK-LD/MC v.6 (May 2018)"
      open(31,file=fileName)
      print *, "Reading in: """, trim(fileName), """"
      read(31,*) nsteps, logWriteFreq
!	store method in temp for translation
      read(31,*) temperature, runType
      if (runType.eq."MC") then
         method = 0
      elseif (runType.eq."OD") then
         method = 1
      else if (runType.eq."UD") then
         method = 2
      else
         stop "Method invalid, please specify MC, OD or UD"
      end if
      if (method.ne.0) then
         read(31,*) dt, friction(1), friction(2)
      else 
         read(31,*) dt
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
      read(31,*) stepSize, threshold, rerun, printOutput
!	QA
      if (threshold.gt.1 .OR. threshold.lt.0.01) then
         stop " Please change threshold to be larger than 0 and smaller than 1.0. Recommended: 0.2-0.5"
      end if
      threshold = 1-threshold
      read(31,*) nParticles
!	QA
      if (nParticles.lt.0) then 
         stop " Cannot use negative number of particles"
      end if
      allocate(particles(nParticles))
!	read in boundaries
      read(31,*) xStart, xEnd, xExtra
      read(31,*) yStart, yEnd, yExtra
      allocate(histogram(int(xStart-2*xExtra):int(xEnd+2*xExtra),int(yStart-2*yExtra):int(yEnd+2*yExtra)))
      read(31,*) boundaryK
      read(31,*) x0, xb, x1
      read(31,*) w0, wb, w1
      read(31,*) y0, yb, y1
!	read in coordinates for the particles
      do i = 1,nParticles
         p = i
         if (method.eq.0) then
!	read MC input
            read(31,*) particles(p)%x,particles(p)%y,particles(p)%z    
         else 
!	read LD input (with mass)
            read(31,*) particles(p)%x,particles(p)%y,particles(p)%z, particles(p)%mass
         end if
      end do
!	check for a destination
      read(31,*,IOSTAT=reason) destination
      if (reason < 0) then
         destination = -1
      else
         print *, "This simulation will stop immediately upon reaching state ", destination
      end if
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
      allocate(rnd(rndSize*6))
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
!	only for UD, print temperatures
      if (method.eq.2) then
         open(203,file='temperature.dat',status='unknown')
      end if
      open(207,file='state.dat',status='unknown')

! @@@@@@@@@@@@@@@@@@@@
! @@ Initialization @@
! @@@@@@@@@@@@@@@@@@@@
!	set all velocities to 0
      particles%xV = 0
      particles%yV = 0
      if (method.eq.0) then
!	calculate and save previous energy
         call applyFF(particles)
         previousEnergy = totalEnergy
      end if
      if (printOutput) then
         histogram = 0
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
!	start state as illegal -1
      state = -1
      state = getState(particles(1)%x,particles(1)%y)
      write(207,*) 0, state

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
!	populate histogram file
         histogram(int(particles(1)%x),int(particles(1)%y)) = histogram(int(particles(1)%x),int(particles(1)%y)) + 1
!	update state if needed
        iState = getState(particles(1)%x,particles(1)%y)
        if (iState.ne.state) then
           state=iState
           write(207,*) istep, state
           if (state .eq. destination) then
              print *, "Destination state reached. Exiting on step", istep
              exit
           end if
        end if
!       print frame to pdb using the user provided frequency
         if (mod(istep,logWriteFreq).eq.0) then
            write(200,'(''MODEL'',i12)') istep
            do p=1,nParticles
!	use z as the energy of the system
               particles(p)%z = getFunctionAt(particles(p)%x,particles(p)%y)
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
      close(207)

      if (printOutput) then
         open(unit=5438,file="energy.pdb",status='unknown')
         open(unit=5439,file="potential.pdb",status='unknown')
         open(unit=5440,file="dd(G-E).pdb",status='unknown')
!	first pass only get the delta energy
         eneG = 0
         do i = xStart-xExtra, xEnd+xExtra
            do j = yStart-xExtra, yEnd+xExtra
               if (histogram(i,j).gt.0) then
                  eneF = -log(real(histogram(i,j)))/0.6
                  eneP = getFunctionAt(dble(i),dble(j))
                  eneG = eneG + (eneF-eneP)
                  if (abs(eneF-eneP).gt.eneMax) then
                     eneMax = abs(eneF-eneP)
                  end if
                  eneCount = eneCount + 1
               end if
            end do
         end do
         eneG = eneG/eneCount
!	second pass print out free energy
         do i = xStart-2*xExtra, xEnd+2*xExtra
            do j = yStart-2*yExtra, yEnd+2*yExtra
               if (histogram(i,j).gt.0) then
                  write(5438,pdbFormat) 1,'A   ','RES',1,real(i), real(j), -log(real(histogram(i,j)))/0.6-eneG
               end if
            end do
         end do
!	print potential only inside boundary
         do i = xStart-xExtra, xEnd+xExtra
            do j = yStart-yExtra, yEnd+yExtra
               write(5439,pdbFormat) 1,'A   ','RES',1,real(i), real(j), getFunctionAt(dble(i),dble(j))
               if (histogram(i,j).gt.0) then

                  write(5440,pdbFormat) 1,'A   ','RES',1,real(i), real(j), -log(real(histogram(i,j)))/0.6-eneG - getFunctionAt(dble(i),dble(j))
               else
                  write(5440,pdbFormat) 1,'A   ','RES',1,real(i), real(j), eneMax+2
               end if
            end do
         end do
         close(5438)
         close(5439)
         close(5440)
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
         else 
            particles(p)%xF = particles(p)%xF + particles(p)%sigma*gasdev(currentSeed)
            particles(p)%yF = particles(p)%yF + particles(p)%sigma*gasdev(currentSeed)
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
         body%x = body%x + body%xF*dt/(body%mass*friction(1))
         body%y = body%y + body%yF*dt/(body%mass*friction(2))
      else
!	underdamp
!	add drag force
         body%xF = body%xF - friction(1)*body%mass*body%xV
         body%yF = body%yF - friction(2)*body%mass*body%yV
!	translate
         body%x  = body%x  + body%xV*dt+0.5d0*(body%xF/body%mass)*dt**2
         body%y  = body%y  + body%yV*dt+0.5d0*(body%yF/body%mass)*dt**2
!	update velocity
         body%xV = body%xV + (body%xF/body%mass)*dt
         body%yV = body%yV + (body%yF/body%mass)*dt
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
      double precision :: xTemp, yTemp, chance

!       get a random particle
!       rnd has 6*rndSize elements
!       for each step it needs a particle number, x,y,z and the metropolis number (in that order)
      i = floor(rnd((rndCounter)*5)*(nParticles))+1
!	from -stepSize/2 to +stepSize/2
!	note that this results in steps that are 0 to 0.87*stepSize in 3D
      xTemp = (rnd((rndCounter)*6+1)-0.5d0)*stepSize
      yTemp = (rnd((rndCounter)*6+2)-0.5d0)*stepSize
      roll   = rnd((rndCounter)*6+4)
!	move the particle
      call moveMC(particles,i,xTemp,yTemp)
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
         call moveMC(particles,i,-xTemp,-yTemp)
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
      subroutine moveMC(particles,i,x,y)
      
      implicit none

      type(Particle), dimension(:), intent(inout) :: particles
      integer :: i 
      double precision :: x,y

!	move particle 
      particles(i)%x=particles(i)%x+x
      particles(i)%y=particles(i)%y+y

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ subroutine to calculate the total energy of the system @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine applyFF(particles)

      implicit none

      type(Particle), dimension(:), intent(inout) :: particles
      double precision :: dx, dy, dummy

      if (method.eq.0) then
!	Set all energies to 0
         totalEnergy = 0.d0
         particles%eTotal = 0.d0
      else
!	set forces to 0
         particles%xF = 0.d0
         particles%yF = 0.d0
      end if

!	add the boundaries in xy
      do p = 1,nParticles
         call boundary(particles(p), .true., xStart-xExtra, .true.)
         call boundary(particles(p), .true., xEnd+xExtra, .false.)
         call boundary(particles(p), .false., yStart-yExtra, .true.)
         call boundary(particles(p), .false., yEnd+yExtra, .false.)
         if (method.eq.0) then
            particles(p)%eTotal = particles(p)%eTotal + getFunctionAt(particles(p)%x,particles(p)%y) 
         else
            dummy = getFunctionAt(particles(p)%x,particles(p)%y,dx,dy)
            particles(p)%xF = particles(p)%xF + dx
            particles(p)%yF = particles(p)%yF + dy
         end if
      end do

!	Accumulate the total energy
      if (method.eq.0) then
         do p = 1,nParticles
            totalEnergy = totalEnergy + particles(p)%eTotal
         end do
      end if

      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Apply boundary force @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine boundary(p, isX, ref, lt)

      implicit none
!	particle
      type(Particle), intent(inout) :: p
!	reference point
      real :: ref
!	x or y?
!	less than or greater than?
      logical :: isX, lt
!	helper - value or x or y
      double precision :: var, f

      f = 0.d0

      if (isX) then
         var = p%x
      else
         var = p%y
      end if

      if ( ( lt .AND. (var .lt. ref)) .OR. (.NOT.lt .AND. (var .gt. ref)) ) then
         if (method.eq.0) then
            f = boundaryK * (ref-var)**2
            p%eTotal = p%eTotal + f
         else
            f = 2.d0 * boundaryK * (ref-var)
            f = f * kcal2cu
            if (isX) then
               p%xF = p%xF + f
            else 
               p%yF = p%yF + f
            end if
         end if
      end if
   
      end subroutine boundary

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Calculate the energy of the system at x y @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double precision function getFunctionAt(x,y,dx,dy)

      implicit none

      double precision, intent(in) :: x,y
      double precision, optional, intent(out) :: dx, dy
      double precision :: f
!	to save the original method
      double precision :: me
!	interpolated x curves (lam interpolant, steep is the steepness of dxdy)
      real :: a0, a1, ab, lam, steep

      f = 0.d0
      
      if (y.le.yStart) then
         a0 = x0
         ab = xb
         a1 = x1
      elseif (y.ge.yEnd) then
         a0 = w0
         ab = wb
         a1 = w1
      else
!	there an issue with the tilted plateau - currenty I don't have a good solution
         lam = (y-yStart)/(yEnd-yStart)
         a0 = x0*(1-lam)+w0*lam
         ab = xb*(1-lam)+wb*lam
         a1 = x1*(1-lam)+w1*lam
      end if

      if (present(dx) .AND. present(dy)) then
         dx = markus(x,xStart,xEnd,a0,ab,a1,stepSize*2,.false.)
         dy = markus(y,yStart,yEnd,y0,yb,y1,stepSize*2,.false.)
!	add force from gradient between the two potentials used
         if (y.gt.yStart .AND. y.lt.yEnd) then
            if (x.lt.(xStart+xEnd)/2-stepSize) then
               steep = w0-x0
            elseif (x.gt.(xStart+xEnd)/2+stepSize) then
               steep = w1-x1
            else
               steep = wb-xb
            end if
            dy = dy - kcal2cu*steep/(yEnd-yStart)
         end if
      else 
         me = method
         method = 0
         f = 0.d0
!	markus for x
         f = markus(x,xStart,xEnd,a0,ab,a1,stepSize*2,.false.)
!	markus for y
         f = f + markus(y,yStart,yEnd,y0,yb,y1,stepSize*2,.false.)
         method = me
      end if

      getFunctionAt = f
      
      end function getFunctionAt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Get the state of the system @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      integer function getState(x,y)

      implicit none
 
      double precision, intent(in) :: x, y 
      integer :: st
      real :: xS, yS

!	states will be - 11 for +x +y
!                        01 for -x +y
!                        10 for +x -y
!                        00 for -x -y
!	for simplicity 0.d0 will be treated as +
! a state is only treated as changed if the difference from the median is at least threshold*size of box
      st = 0
      xS = threshold*((xEnd-xStart)/2)
      yS = threshold*((yEnd-yStart)/2)

!       starting from positive x
      if (state.eq.10 .OR. state.eq.11) then
         if (x.ge.(xStart+xEnd)/2-xS) then
            st = st+10
         end if
      end if
!       starting from negative x
      if (state.eq.1 .OR. state.eq.0) then
         if (x.ge.(xStart+xEnd)/2+xS) then
            st = st+10
         end if
      end if
!       starting from positive y
      if (state.eq.1 .OR. state.eq.11) then
         if (y.ge.(yStart+yEnd)/2-yS) then
            st = st+1
         end if
      end if
!       starting from negative y
      if (state.eq.0 .OR. state.eq.10) then
         if (y.ge.(yStart+yEnd)/2+yS) then
            st = st+1
         end if
      end if
      
      getState = st

      end function getState

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Calculate the energy of marcus parbolas @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double precision function markus(x, startX, endX, startE, barrierE, endE, gap, flatten)
!	modified markus with different skewness, and no h12

!	startX and endX are the points in real space where the potential should begin and end
!	startE and endE are the energies at the beginning and ending
!	barrier is the energy at the peak
!	gap is the length of the gap (the plateau at the peak)

      implicit none

      double precision :: t, barrier, a1, a2, output
      double precision, intent(in) :: x
      real, intent(in) :: startX, startE, endX, endE, barrierE, gap
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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Temperature (kinetic) @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double precision function getTemperature(myosin)

      implicit none

      type(Particle), dimension(:), intent(in) :: myosin
      integer :: i
      double precision :: totalK

      totalK = 0.d0

!	sum the kinetic energy of all particles
      do i = 1, nParticles
         totalK = totalK + 0.5d0*myosin(i)%mass*(myosin(i)%xV**2+myosin(i)%yV**2+myosin(i)%zV**2)
      end do

!	divide by Bolzmann and n of particles (pLast)
      getTemperature = totalK/(bz_cu*nParticles)

      end function

!@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Precalculate sigma's @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine calculateSigma(particles)

      implicit none

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
