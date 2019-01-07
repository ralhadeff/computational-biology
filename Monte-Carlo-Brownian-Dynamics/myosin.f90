! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ LD/MC code for the myosin-V system @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! important comment: MC code might have bugs, don't use it without proper testing!!!

      program main

!	written by Raphael Alhadeff April 2017
!     version updated Jan 2018

      implicit none

!	conversion constant
      real, parameter :: bz_cu = 0.83146, kcal2cu = 418.3965
      double precision, parameter :: piValue = 3.14159265359
!	write frequency, i for local loops, type of run (0-MC,1-Overdamp LD,2-Underdamp LD)
!	iLast to make the code more readable
!	I will use p for index of particles and d for indice of domains
      integer :: nMyosin, nActin, logWriteFreq, method, p, d, i, pLast, dLast
!	total number of step and current step index
      integer(8) :: nsteps, istep
!	object for particle
!	position, velocity, force xyz 
!	index of bound entity (actin for myosin, myosin for adp and adp for pi)
!	sigma, mass, total energy 
      type Particle
         double precision :: x, y, z, xV, yV, zV, xF, yF, zF, eTotal, sigma
         real :: mass
         integer :: bound
      end type
!	particles arrays: pi is phosphate. ATP will be adp+p adjacent
      type (Particle), dimension(:), allocatable :: myosin, actin
      type (Particle), dimension(4) :: adp, pi
!	domain type is an octahedron
!	center and plus XYZ or minus XYZ
!	all indices of myosin particle from myosin array
      type DomainT
         integer :: center, pX, pY, pZ
      end type
      type (DomainT), dimension(:), allocatable :: domain
!	foot type 
!	index of domain and indices of adp and atp
      type Foot
         integer :: domain, nuc1, nuc2, knee
      end type
      type (Foot) :: feet(2)
!       total energy, previous total energy, temperature(K), dt and friction if LD
!	current roll for MC
      double precision :: totalEnergy, previousEnergy, roll
      real :: temperature, dt, friction
!	xyz coordinate for first actin particle, dx,dy,dz for each consequent particle in the filament
      real :: xActin, dxActin, yActin, dyActin, zActin, dzActin
!	input file name
      character(len=30) :: fileName
!	FF: in main program so it can be read from input file
!	two sides truncated harmonic angle ka and max/min values
!	harmonic bond, and minimum distance for non bonded rejection
      real :: ka, thetaMax, thetaMin, kb, r0, rMin
!	torsion angle parameters
      real :: kt, phi0
!	foot angle thresholds, dG and barrier
!	footSpacer is the size of the transition zone for lambda (see footAngle subroutine)
      real :: thetaFootUp, thetaFootDown, footDG, footBarrier, footSpacer
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
!	whether to pring energies and histogram files
!	whether to calculate and print states
!	whether to stop running (for early stopping)
      logical :: rerun, printOutput, printStates, kill
!	actin myosin interactions
      real :: actMyoSpacer, actMyoLength, couplingRatio, actMyoBarrierHigh, actMyoBarrierLow, actMyoBind, couplingMultiplyer
!	myosin nucleotide interactions
      integer :: nucRegime
      real :: myoNucSpacer, nucLen, myoNucHighBarrier, myoNucLowBarrier, myoNucBind
      real :: nucleotideBoundary, nucleotideBoundaryK, phosphateBoundary
!	ADP P interactions
      real :: adpPiSpacer, adpPiBarrierW, adpPiBarrierP, adpPi2ndBarrierW, adpPi2ndBarrierP, adpPiBind, adpPi2ndBind
!	force on feet, force constant is linear in the positive x direction (starting from x=0)
      real :: pullForce
!	state of the system stuff
!	will use integer codes as follows:
!	myoAngle = -1 down, +1 up
!	myoBind = 1 bound, 0 unbound
!	myoNuc = 0 nothing, 2 adp, 3 atp, 4 adp+p
      integer, dimension(2) :: myoAngleS, myoBindS, myoNucS
!	format for printing pdb's
      character(len=128) :: pdbFormat="('ATOM',2x,i5,2x,a4,a3,i6,4x,3f8.3)"
!	histogram counters
!	first index is the particle index, second in then the distance bin
!	angle is 0.00 to 3.50 in 0.01 increments, rest is 1 A increments
      integer(8) :: PiDistance(4,0:100), ADPDistance(4,0:100), MyoDistance(2,0:900), MyoAngle(2,-350:350)
!	text string for printouts
      character(len=128) :: text, temp, vFormat

! @@@@@@@@@@@@@@@@
! @@ read input @@
! @@@@@@@@@@@@@@@@
!       read input file from the command line
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
      read(31,*) stepSize, rerun, printOutput, printStates
      read(31,*) nMyosin, nActin
!	QA
      if (nMyosin.lt.0) then 
         stop " Cannot use negative number of joints"
      else
         write(temp,*), nMyosin
         print *, "Using ", trim(adjustl(temp)), " lever particles (each lever), excluding the joint and heads"
      end if
      allocate(myosin(4*(2*nMyosin+3)))
      allocate(domain(2*nMyosin+3))
      if (nActin.lt.0) then
         stop " Cannot use negative number of actin particles"
      end if
      allocate(actin(nActin))
!	keep reading
      read(31,*) kb, r0, rMin
      read(31,*) ka, thetaMin, thetaMax
!	store the symbol for the regime in temp
      read(31,*) thetaFootDown, thetaFootUp, footSpacer, temp
      if (temp.eq."D") then
         nucRegime = -1
         print *, "Nucleotide relase is allowed below the down angle only"
      elseif (temp.eq."U") then
         nucRegime = 1
         print *, "Nucleotide relase is allowed above the up angle only"
      elseif (temp.eq."UC") then
         nucRegime = 0
         print *, "Nucleotide relase is allowed only when both myosins are bound (uncoupled)"
      else
         stop "Nucleotide release/bind regime invalide, please specify U, D or UC"
      end if
      read(31,*) kt, phi0
      read(31,*) footDG, footBarrier
      read(31,*) actMyoSpacer, actMyoLength, couplingRatio, couplingMultiplyer
      read(31,*) actMyoBarrierLow, actMyoBarrierHigh
      read(31,*) actMyoBind
      read(31,*) myoNucSpacer, nucLen 
      read(31,*) myoNucLowBarrier, myoNucHighBarrier, myoNucBind
      read(31,*) adpPiSpacer, adpPiBarrierW, adpPiBarrierP
      read(31,*) adpPi2ndBarrierW, adpPi2ndBarrierP
      read(31,*) adpPiBind, adpPi2ndBind
      read(31,*) nucleotideBoundary, phosphateBoundary, nucleotideBoundaryK
      read(31,*) pullForce
      if (pullForce.ne.0) then
         if (pullForce.gt.0) then
            text = "A pulling force of "
            write(temp,'(f8.3)'), pullForce
         else
            text = "A pushing force of "
            i = floor(log10(abs(pullForce)))+6
            write(vFormat,*), i
            write(vFormat,*), "(f", trim(adjustl(vFormat)), ".3)"
            write(temp,vFormat), abs(pullForce)
         end if 
         print *, trim(text), " ",trim(temp), " is applied to the joint particle"
      end if
      read(31,*) xActin, dxActin, yActin, dyActin, zActin, dzActin
!	to make things clear:
!	particles 1-4 - one foot, then 5...4*(n+1) are the first lever (could be 0) then 4*(n+1)+1...4*(n+2) is the joint
!	then 4*(n+2)+1...4*(2n+2) is the second lever and lastly 4*(2n+2)+1...4*(2n+3) is the second foot
!	reading in ONLY the domains' centers, the levers and joint (linearly)
      do i = 0,2*nMyosin+2
         p = 4*i+1
         if (method.eq.0) then
            read(31,*) myosin(p)%x,myosin(p)%y,myosin(p)%z    
         else 
            read(31,*) myosin(p)%x,myosin(p)%y,myosin(p)%z, myosin(p)%mass
         end if
      end do
!	nucleotides
      do i = 1,4
         if (method.eq.0) then
            read(31,*) adp(i)%x, adp(i)%y, adp(i)%z
            read(31,*) pi(i)%x,   pi(i)%y,   pi(i)%z
         else
            read(31,*) adp(i)%x, adp(i)%y, adp(i)%z, adp(i)%mass
            read(31,*) pi(i)%x,   pi(i)%y,   pi(i)%z, pi(i)%mass
         end if
      end do
      close(31)
!       setup lasts
      pLast = 4*(2*nMyosin+3)
      dLast = 2*nMyosin+3
!       setup all domains and feet indices
      do d = 1,dLast
         p = 4*(d-1)+1
         domain(d)%center = p
         domain(d)%pX = p+1
         domain(d)%pY = p+2
         domain(d)%pZ = p+3
      end do
      feet(1)%domain = 1
      feet(1)%nuc1 = 1
      feet(1)%nuc2 = 2
      feet(1)%knee = 2
      feet(2)%domain = dLast
      feet(2)%nuc1 = 3
      feet(2)%nuc2 = 4
      feet(2)%knee = dLast-1
!	update coordinates for the sub particles
!	TODO start aligned to angles
      do d = 1,dLast
!	+x
         myosin(domain(d)%pX)%x = myosin(domain(d)%center)%x+r0/2
         myosin(domain(d)%pX)%y = myosin(domain(d)%center)%y
         myosin(domain(d)%pX)%z = myosin(domain(d)%center)%z
         myosin(domain(d)%pX)%mass = myosin(domain(d)%center)%mass
!       +y
         myosin(domain(d)%pY)%x = myosin(domain(d)%center)%x
         myosin(domain(d)%pY)%y = myosin(domain(d)%center)%y+r0/2
         myosin(domain(d)%pY)%z = myosin(domain(d)%center)%z
         myosin(domain(d)%pY)%mass = myosin(domain(d)%center)%mass
!       +z
         myosin(domain(d)%pZ)%x = myosin(domain(d)%center)%x
         myosin(domain(d)%pZ)%y = myosin(domain(d)%center)%y
         myosin(domain(d)%pZ)%z = myosin(domain(d)%center)%z+r0/2
         myosin(domain(d)%pZ)%mass = myosin(domain(d)%center)%mass
      end do
!       Create X, Y and Z coordinates of nActin particles
      do i = 1,nActin
         actin(i)%x = xActin + (i-1) * dxActin
         actin(i)%y = yActin + (i-1) * dyActin
         actin(i)%z = zActin + (i-1) * dzActin
      end do

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
!	6 because I need 1 number to choose the particle, 3 for xyz, and then the 5th for the roll
!	6th for the occasional roll for ATP/ADP stuff
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
!	start is the first frame of the movie for visualization purposes
!	actin.pdb is an actin string for visualization purposes
      open(200,file='movie.pdb',status='unknown')
      open(201,file='start.pdb',status='unknown')
      open(199,file='actin.pdb',status='unknown')
!	only for UD, print temperatures
      if (method.eq.2) then
         open(203,file='temperature.dat',status='unknown')
      end if
      if (printStates) then
         open(202,file='states.dat',status='unknown')
         open(333,file='cases.pdb',status='unknown')
      end if
      do d=1,dLast
         p = domain(d)%center
         write(201,pdbFormat) p,'H   ','MYO',p,myosin(p)%x,myosin(p)%y,myosin(p)%z
         p = domain(d)%pX
         write(201,pdbFormat) p,'OH  ','MYO',p,myosin(p)%x,myosin(p)%y,myosin(p)%z
         p = domain(d)%pY
         write(201,pdbFormat) p,'S   ','MYO',p,myosin(p)%x,myosin(p)%y,myosin(p)%z
         p = domain(d)%pZ
         write(201,pdbFormat) p,'NA  ','MYO',p,myosin(p)%x,myosin(p)%y,myosin(p)%z
      end do
      do i = 1,4
         write(201,pdbFormat) pLast+2*i-1,'CA  ','ADP',i,adp(i)%x,adp(i)%y,adp(i)%z
         write(201,pdbFormat) pLast+2*i,'P   ','PHO',i,pi(i)%x,pi(i)%y,pi(i)%z
      end do
!	connections in myosin file
      do d=1,dLast-1
!       between domains
         write(201,'(''CONECT'',i5,i5)') domain(d)%center,domain(d+1)%center
      end do
      do d=1,dLast
!       inside each domain
         write(201,'(''CONECT'',i5,i5)') domain(d)%center,domain(d)%pX
         write(201,'(''CONECT'',i5,i5)') domain(d)%center,domain(d)%pY
         write(201,'(''CONECT'',i5,i5)') domain(d)%center,domain(d)%pZ
      end do
!	add bonds to ADP+P
      do i=1,4
         write(201,'(''CONECT'',i5,i5)') pLast+2*i-1,pLast+2*i
      end do
!	write actin particles to file
      do i=1,nActin
         write(199,pdbFormat) i,'OH  ','ACT', i, actin(i)%x, actin(i)%y, actin(i)%z
      enddo
!	write connect lines
      do i=1,nActin-1
         write(199,'(''CONECT'',i5,i5)') i,i+1
      enddo
      close(199)
      close(201)

! @@@@@@@@@@@@@@@@@@@@
! @@ Initialization @@
! @@@@@@@@@@@@@@@@@@@@
!	set all velocities to 0
      myosin%xV = 0
      myosin%yV = 0
      myosin%zV = 0
      adp%xV = 0
      adp%yV = 0
      adp%zV = 0
      pi%xV = 0
      pi%yV = 0
      pi%zV = 0
      if (method.eq.0) then
!	calculate and save previous energy
         call applyFF(myosin,adp,pi)
         previousEnergy = totalEnergy
      end if
      if (printOutput) then
!	set histogram arrays to zero
         PiDistance  = 0
         ADPDistance = 0
         MyoDistance = 0
         MyoAngle    = 0
      end if
      if (method.ne.0) then
!	precalculate sigma's and get first seed
         call calculateSigma(myosin,adp,pi)
         if (seed(1) .gt. 0) then
            currentSeed = -seed(1)
         else
            currentSeed = seed(1)
         end if
      end if
!	TODO crap writing
!	initialize bound
      call markBound(myosin)
!	TODO add state initialization
      do i = 1,2
         myoBindS(i) = -1
         myoAngleS(i) = 0
         myoNucS(i)  = -1
      end do

! @@@@@@@@@@@@@@@@@@@@@
! @@ simulation loop @@
! @@@@@@@@@@@@@@@@@@@@@
      do istep = 1,nsteps
!	apply method
         if (method.eq.0) then
            call doMC(myosin,adp,pi)
         else
            call doLD(myosin,adp,pi)
         end if
!       print frame to pdb using the user provided frequency
         if (mod(istep,logWriteFreq).eq.0) then
            write(200,'(''MODEL'',i12)') istep
            do p=1,pLast
               write(200,pdbFormat) p,'OH  ','MYO',p,myosin(p)%x,myosin(p)%y,myosin(p)%z
            enddo
            do i=1,4
               write(200,pdbFormat) i*2-1+pLast,'NA  ','ADP',i,adp(i)%x,adp(i)%y,adp(i)%z
               write(200,pdbFormat) i*2+pLast,'P   ','PHO',i,pi(i)%x,pi(i)%y,pi(i)%z
            end do
            write(200,'(''TER'')')
            write(200,'(''ENDMDL'',i12)') istep
            if (method.eq.2) then
               write(203,*) getTemperature(myosin,adp,pi)
            end if
         endif
!	state stuff
         if (printStates) then
            call checkState(myoAngleS, myoBindS, myoNucS)
         end if
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
               call printHistograms 
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
      if (printStates) then
         close(202)
         close(333)
      end if

      if (printOutput) then
!	output subroutines use energy and not forces
!	for markus subroutine to work properly set method to MC
         method = 0
         call printOutLandscapes
         call printHistograms
      end if

!	end of main routine
      contains
! @@@@@@@@@@@@@@@@@@@@@@@@@
! @@ End of main program @@
! @@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine markBound(myosin)

      implicit none
 
      type (Particle), dimension(:), intent(inout) :: myosin
      real :: actMyoDist, mDist
      integer :: i, iBound, f

      do f = 1,2
         actMyoDist=100000.0
         do i=1,nActin
            mDist = distance(myosin(domain(feet(f)%domain)%center),actin(i))
            if(mDist.lt.actMyoDist)then
               actMyoDist=mDist
               iBound=i
            endif
         end do
         myosin(domain(feet(f)%domain)%center)%bound = iBound
      end do

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@
! @@ perform a LD step @@
! @@@@@@@@@@@@@@@@@@@@@@@
      subroutine doLD(myosin,adp,pi)

      implicit none

      type(Particle), dimension(:), intent(inout) :: myosin,adp,pi

!	reset forces and apply potential
      call applyFF(myosin,adp,pi)
!	generate random forces
!	myosin
      do p = 1,pLast
         if (method.eq.1) then
            myosin(p)%x = myosin(p)%x + myosin(p)%sigma*gasdev(currentSeed)
            myosin(p)%y = myosin(p)%y + myosin(p)%sigma*gasdev(currentSeed)
            myosin(p)%z = myosin(p)%z + myosin(p)%sigma*gasdev(currentSeed)
         else 
            myosin(p)%xF = myosin(p)%xF + myosin(p)%sigma*gasdev(currentSeed)
            myosin(p)%yF = myosin(p)%yF + myosin(p)%sigma*gasdev(currentSeed)
            myosin(p)%zF = myosin(p)%zF + myosin(p)%sigma*gasdev(currentSeed)
         end if
      end do 
!	adp and pi
      do i = 1,4
         if (method.eq.1) then
            adp(i)%x = adp(i)%x + adp(i)%sigma*gasdev(currentSeed)
            adp(i)%y = adp(i)%y + adp(i)%sigma*gasdev(currentSeed)
            adp(i)%z = adp(i)%z + adp(i)%sigma*gasdev(currentSeed)
            pi(i)%x = pi(i)%x + pi(i)%sigma*gasdev(currentSeed)
            pi(i)%y = pi(i)%y + pi(i)%sigma*gasdev(currentSeed)
            pi(i)%z = pi(i)%z + pi(i)%sigma*gasdev(currentSeed)
         else
            adp(i)%xF = adp(i)%xF + adp(i)%sigma*gasdev(currentSeed)
            adp(i)%yF = adp(i)%yF + adp(i)%sigma*gasdev(currentSeed)
            adp(i)%zF = adp(i)%zF + adp(i)%sigma*gasdev(currentSeed)
            pi(i)%xF = pi(i)%xF + pi(i)%sigma*gasdev(currentSeed)
            pi(i)%yF = pi(i)%yF + pi(i)%sigma*gasdev(currentSeed)
            pi(i)%zF = pi(i)%zF + pi(i)%sigma*gasdev(currentSeed)
         end if
      end do
!       do mechanics
      do p = 1,pLast
         call doMechanics(myosin(p))
      end do
!       adp and pi
      do i = 1,4
         call doMechanics(adp(i))
         call doMechanics(pi(i))
      end do

!	reconstruct ATP or ADP if needed
      call checkATPLD(myosin,adp,pi)

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
      subroutine doMC(myosin,adp,pi)


      type(Particle), dimension(:), intent(inout) :: myosin,adp,pi

!       x,y,z_temp is the vector coordinates of the MC step
!       chance is the boltzmann probability for the metropolis criteria
!	roll is a global variable because it is sometimes used later in checkATP
      double precision :: xTemp, yTemp, zTemp, chance

!       get a random particle
!       rnd has 6*rndSize elements
!       for each step it needs a particle number, x,y,z and the metropolis number (in that order)
!       number of myosin particles (2 levers + 2 heads + 1 joint) +4 ADP +4 P +1 to start from 1
      i = floor(rnd((rndCounter)*6)*(pLast+8))+1
!	from -stepSize/2 to +stepSize/2
!	note that this results in steps that are 0 to 0.87*stepSize in 3D
      xTemp = (rnd((rndCounter)*6+1)-0.5d0)*stepSize
      yTemp = (rnd((rndCounter)*6+2)-0.5d0)*stepSize
      zTemp = (rnd((rndCounter)*6+3)-0.5d0)*stepSize
      roll   = rnd((rndCounter)*6+4)
!	move the particle
      call moveMC(myosin,adp,pi,i,xTemp,yTemp,zTemp)
!       get energy after change
      call applyFF(myosin,adp,pi)
!       calculate chance to pass (metropolis)
      chance = exp((totalEnergy-previousEnergy)/(-0.002d0*temperature))
!       if accepted:
!       save the previous energy, and update boolean
      if (chance.gt.roll) then
         previousEnergy = totalEnergy
!       if the accepted move was an ADP particle
!       check if an ATP particle needs to be restored
         if (i.gt.pLast .AND. i.lt.pLast+5) then
!	determine which foot this is and pass the foot
            if (i.eq.feet(1)%nuc1+pLast .OR. i.eq.feet(1)%nuc2+pLast) then
               call checkATP(feet(1),myosin,adp,pi)
            else
               call checkATP(feet(2),myosin,adp,pi)
            end if
         end if
      else
!       if rejected revert the change
         call moveMC(myosin,adp,pi,i,-xTemp,-yTemp,-zTemp)
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
      subroutine moveMC(myosin,adp,pi,i,x,y,z)
      
      implicit none

      type(Particle), dimension(:), intent(inout) :: myosin,adp,pi
!	i is the index of the particle.
!	1-nMyosin is myosin
!	next 4 are ADP and next 4 are Pi
!	ADP and Pi will move with myosin feet for ATP binding smoothness
!	Pi will also move with ADP
!	ADP and Pi will move at 1/10th the stepsize to scale down their coordinates (for visual purposes)
      integer :: i, a, b
      double precision :: x,y,z

!	move particle 
      if (i.eq.domain(feet(1)%domain)%center) then
!	first foot, move ADP and Pi with it
         myosin(i)%x=myosin(i)%x+x
         myosin(i)%y=myosin(i)%y+y
         myosin(i)%z=myosin(i)%z+z
         a = feet(1)%nuc1
         b = feet(1)%nuc2
         adp(a)%x=adp(a)%x+x
         adp(a)%y=adp(a)%y+y
         adp(a)%z=adp(a)%z+z
         pi(a)%x=pi(a)%x+x
         pi(a)%y=pi(a)%y+y
         pi(a)%z=pi(a)%z+z
         adp(b)%x=adp(b)%x+x
         adp(b)%y=adp(b)%y+y
         adp(b)%z=adp(b)%z+z
         pi(b)%x=pi(b)%x+x
         pi(b)%y=pi(b)%y+y
         pi(b)%z=pi(b)%z+z
      else if (i.eq.domain(feet(2)%domain)%center) then
!       last foot, move ADP and Pi with it
         myosin(i)%x=myosin(i)%x+x
         myosin(i)%y=myosin(i)%y+y
         myosin(i)%z=myosin(i)%z+z
         a = feet(2)%nuc1
         b = feet(2)%nuc2
         adp(a)%x=adp(a)%x+x
         adp(a)%y=adp(a)%y+y
         adp(a)%z=adp(a)%z+z
         pi(a)%x=pi(a)%x+x
         pi(a)%y=pi(a)%y+y
         pi(a)%z=pi(a)%z+z
         adp(b)%x=adp(b)%x+x
         adp(b)%y=adp(b)%y+y
         adp(b)%z=adp(b)%z+z
         pi(b)%x=pi(b)%x+x
         pi(b)%y=pi(b)%y+y
         pi(b)%z=pi(b)%z+z
      else if (i.le.pLast) then
!       only myosin particle
         myosin(i)%x=myosin(i)%x+x
         myosin(i)%y=myosin(i)%y+y
         myosin(i)%z=myosin(i)%z+z
      else if (i.lt.pLast+5)then
!       adp, move Pi with it; 1/10 step
         adp(i-pLast)%x=adp(i-pLast)%x+x/10
         adp(i-pLast)%y=adp(i-pLast)%y+y/10
         adp(i-pLast)%z=adp(i-pLast)%z+z/10
         pi(i-pLast)%x=pi(i-pLast)%x+x/10
         pi(i-pLast)%y=pi(i-pLast)%y+y/10
         pi(i-pLast)%z=pi(i-pLast)%z+z/10
      else
!       pi 1/10 step
         pi(i-pLast-4)%x=pi(i-pLast-4)%x+x/10
         pi(i-pLast-4)%y=pi(i-pLast-4)%y+y/10
         pi(i-pLast-4)%z=pi(i-pLast-4)%z+z/10
      end if

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ subroutine to calculate the total energy of the system @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine applyFF(myosin,adp,pi)

      implicit none

      type(Particle), dimension(:), intent(inout) :: myosin,adp,pi

!	have to use a private i because global i might be in use
      integer :: i, j, pMiddle, dMiddle 
!	indices to make code more readable
      integer, dimension(4) :: myo

      if (method.eq.0) then
!	Set all energies to 0
         totalEnergy = 0.d0
         myosin%eTotal = 0.d0
         adp%eTotal = 0.d0
         pi%eTotal = 0.d0
      else
!	set forces to 0
         myosin%xF = 0.d0
         myosin%yF = 0.d0
         myosin%zF = 0.d0
         adp%xF = 0.d0
         adp%yF = 0.d0
         adp%zF = 0.d0
         pi%xF = 0.d0
         pi%yF = 0.d0
         pi%zF = 0.d0
      end if

!       define indices of the middle
      dMiddle = nMyosin+2
      pMiddle = domain(nMyosin+2)%center

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ actin myosin interactions @@
! @@ also feet angles          @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      do i = 1, 2
         call footEnergy(myosin, adp, pi, feet(i))
      end do
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ Biasing force to accelerate diffusion @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      call biasDiffusion(myosin, adp, pi)

! @@@@@@@@@@@@@@@@@@@@@
! @@ Tweezers energy @@
! @@@@@@@@@@@@@@@@@@@@@
      call pull(myosin(pMiddle),pullForce)

! @@@@@@@@@@@@@@
! @@ bonds    @@
! @@ angles   @@
! @@ torsions @@
! @@@@@@@@@@@@@@
!	inter domain bonds
      do d = 1, dLast-1
         call getBondEnergy(kb,r0,myosin(domain(d)%center),myosin(domain(d+1)%center))
      end do
!       internal particles
      do d = 1, dLast
         call getBondEnergy(kb,r0/2,myosin(domain(d)%center),myosin(domain(d)%pX))
         call getBondEnergy(kb,r0/2,myosin(domain(d)%center),myosin(domain(d)%pY))
         call getBondEnergy(kb,r0/2,myosin(domain(d)%center),myosin(domain(d)%pZ))
         call getAngleEnergy(ka,real(piValue/2),real(piValue/2),myosin(domain(d)%pX),myosin(domain(d)%center),myosin(domain(d)%pY))
         call getAngleEnergy(ka,real(piValue/2),real(piValue/2),myosin(domain(d)%pX),myosin(domain(d)%center),myosin(domain(d)%pZ))
         call getAngleEnergy(ka,real(piValue/2),real(piValue/2),myosin(domain(d)%pY),myosin(domain(d)%center),myosin(domain(d)%pZ))
!       for chirality
         call getTorsion(kt,real(piValue/2),myosin(domain(d)%pX),myosin(domain(d)%center),myosin(domain(d)%pZ),myosin(domain(d)%pY))
      end do
!	inter domains angles
!	split because the two feet are inverted
!	first foot
      d = 1
!       to prevent collapse of the plane (affecting torsions below)
!	I had to increase the constant because it was still causing problems
         call getAngleEnergy(5*ka,real(0.2*piValue),real(0.8*piValue),myosin(domain(d+1)%center),myosin(domain(d)%center),myosin(domain(d)%pY))
!	skip first angle because it is treated in the foot energy (lever-up or lever-down)
      call getAngleEnergy(ka,thetaMin,thetaMax,myosin(domain(d)%center),myosin(domain(d+1)%center),myosin(domain(d+1)%pZ))
      call getTorsion(kt,0.0,myosin(domain(d)%pX),myosin(domain(d)%center),myosin(domain(d+1)%center),myosin(domain(d+1)%pX))
      myo(1) = domain(d)%pX
      myo(2) = domain(d)%center
      myo(3) = domain(d+1)%center
      myo(4) = domain(d+1)%pY
      call getTorsion(kt,real(piValue/2),myosin(myo(1)),myosin(myo(2)),myosin(myo(3)),myosin(myo(4)))
!	first lever
      do d = 2, dMiddle-2
         call getAngleEnergy(ka,real(piValue),real(piValue),myosin(domain(d)%center),myosin(domain(d)%pZ),myosin(domain(d+1)%center))
         call getAngleEnergy(ka,thetaMin,thetaMax,myosin(domain(d)%center),myosin(domain(d+1)%center),myosin(domain(d+1)%pZ))
         call getTorsion(kt,0.0,myosin(domain(d)%pX),myosin(domain(d)%center),myosin(domain(d+1)%center),myosin(domain(d+1)%pX))
         myo(1) = domain(d)%pX
         myo(2) = domain(d)%center
         myo(3) = domain(d+1)%center
         myo(4) = domain(d+1)%pY
         call getTorsion(kt,real(piValue/2),myosin(myo(1)),myosin(myo(2)),myosin(myo(3)),myosin(myo(4)))         
      end do
!	angle to joint
      d = dMiddle-1
      call getAngleEnergy(ka,thetaMin,thetaMax,myosin(domain(d)%center),myosin(domain(d)%pZ),myosin(domain(d+1)%center))
!	to couple between the two lever's XYZ's systems
      call getTorsion(kt,0.0,myosin(domain(d)%pY),myosin(domain(d)%center),myosin(domain(d+1)%center),myosin(domain(d+1)%pY))
!	last foot
      d = dLast
!       to prevent collapse of the plane (affecting torsions below)
         call getAngleEnergy(5*ka,real(0.2*piValue),real(0.8*piValue),myosin(domain(d-1)%center),myosin(domain(d)%center),myosin(domain(d)%pY))
!       skip first angle because it is treated in the foot energy (lever-up or lever-down)
      call getAngleEnergy(ka,thetaMin,thetaMax,myosin(domain(d)%center),myosin(domain(d-1)%center),myosin(domain(d-1)%pZ))
      call getTorsion(kt,0.0,myosin(domain(d)%pX),myosin(domain(d)%center),myosin(domain(d-1)%center),myosin(domain(d-1)%pX))
      myo(1) = domain(d)%pX
      myo(2) = domain(d)%center
      myo(3) = domain(d-1)%center
      myo(4) = domain(d-1)%pY
      call getTorsion(kt,real(piValue/2),myosin(myo(1)),myosin(myo(2)),myosin(myo(3)),myosin(myo(4)))
!	second lever
      do d = dLast-1, dMiddle+2, -1
         call getAngleEnergy(ka,real(piValue),real(piValue),myosin(domain(d)%center),myosin(domain(d)%pZ),myosin(domain(d-1)%center))
         call getAngleEnergy(ka,thetaMin,thetaMax,myosin(domain(d)%center),myosin(domain(d-1)%center),myosin(domain(d-1)%pZ))
        call getTorsion(kt,0.0,myosin(domain(d)%pX),myosin(domain(d)%center),myosin(domain(d-1)%center),myosin(domain(d-1)%pX))
         myo(1) = domain(d)%pX
         myo(2) = domain(d)%center
         myo(3) = domain(d-1)%center
         myo(4) = domain(d-1)%pY
         call getTorsion(kt,real(piValue/2),myosin(myo(1)),myosin(myo(2)),myosin(myo(3)),myosin(myo(4)))
      end do
!	angle to joint
      d = dMiddle+1
      call getAngleEnergy(ka,thetaMin,thetaMax,myosin(domain(d)%center),myosin(domain(d)%pZ),myosin(domain(d-1)%center))
!       to couple between the two lever's XYZ's systems
      call getTorsion(kt,0.0,myosin(domain(d)%pY),myosin(domain(d)%center),myosin(domain(d-1)%center),myosin(domain(d-1)%pY))

! @@@@@@@@@@@@@@@
! @@ Rejection @@
! @@@@@@@@@@@@@@@
!	rejection between the lever pairs and feet domains (w.r.t. joint)
      do i = 1,nMyosin+1
         call getRejection(kb,rMin,myosin(domain(i)%center),myosin(domain(dLast-(i-1))%center))
      end do

!	Accumulate the total energy
      if (method.eq.0) then
         do p = 1,pLast
            totalEnergy = totalEnergy + myosin(p)%eTotal
         end do
      end if

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@
! @@ Foot pull energy @@ 
! @@@@@@@@@@@@@@@@@@@@@@
      subroutine pull(myosin,force)

      implicit none
 
      type (Particle), intent(inout) :: myosin
      real :: force

      if (method.ne.0) then
!	force
         myosin%xF = myosin%xF - force*kcal2cu
      else 
!	energy
         myosin%eTotal = myosin%eTotal + myosin%x * force    
      end if 

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ Foot energy, including:   @@ 
! @@ myosin-actin, myosin-ADP  @@
! @@ ADP-Pi, ADP/Pi boundaries @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine footEnergy(myosin, adp, pi, indices)

      implicit none

      type(Particle), dimension(:), intent(inout) :: myosin,adp,pi
      type(Foot), intent(in) :: indices

!	distances to be used
      double precision :: mDist, pDist, myoAdpDist, myoAdpDist2, adpPDist, adpPDist2, actMyoDist
!	for nucleotide regime (optional)
      double precision :: actMyo1, actMyo2
!	closest's, and distance between ADPs
      double precision :: adpClosest, piClosest, interADP
!	angle to be used, energy to be used, e to be recycled, lambda's to be recycled
      double precision :: footTheta, energy, e, lambda, lambda2
!	energy for the coupling term
      double precision :: fd, alpha, Ua
!	index of actin and loops
      integer iBound, i, iP
!	helper indices to make code more readable
      integer, dimension(4) :: myo
!	histogram stuff
      integer :: floored
!	foot angle force stuff
      double precision :: dx, dy, dz, r, cosTheta, st
      double precision, dimension(:) :: df(3)

!	startup
      energy = 0
!	ankle binds the nucleotides
      myoAdpDist = distance(myosin(domain(indices%domain)%center),adp(indices%nuc1))
      myoAdpDist2 = distance(myosin(domain(indices%domain)%center),adp(indices%nuc2))
      adpPDist = distance(adp(indices%nuc1), pi(indices%nuc1))
      adpPDist2 = distance(adp(indices%nuc2), pi(indices%nuc2))
!	get the foot angle
      myo(1) = domain(indices%knee)%center
      myo(2) = domain(indices%domain)%center
      myo(3) = domain(indices%domain)%pY
      myo(4) = domain(indices%domain)%pZ
      call getTorsion(0.0,0.0,myosin(myo(1)),myosin(myo(2)),myosin(myo(3)),myosin(myo(4)),1,footTheta)
!	distance between ADPs
      interADP = distance(adp(indices%nuc1), adp(indices%nuc2))
!       find the closest actin particle to particle of Myosin
!	for simplicity distance will be between the tail particle (heel) and the back actin particle from the doublets
      actMyoDist=100000.d0
      do i=1,nActin
         mDist = distance(myosin(domain(indices%domain)%center),actin(i))
         if(mDist.lt.actMyoDist)then
            actMyoDist=mDist
            iBound=i
         endif
      end do
      myosin(domain(indices%domain)%center)%bound = iBound

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ Myosin actin interactions @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!       TODO rubbish writing
      actMyo1 = distance(myosin(domain(feet(1)%domain)%center),actin(myosin(domain(feet(1)%domain)%center)%bound))
      actMyo2 = distance(myosin(domain(feet(2)%domain)%center),actin(myosin(domain(feet(2)%domain)%center)%bound))
!	energy only affected by the closest ADP
      if (myoAdpDist.lt.myoAdpDist2) then
         adpClosest = myoAdpDist
         piClosest  = adpPDist
      else
         adpClosest = myoAdpDist2
         piClosest  = adpPDist2
      end if
!	myosin - actin binding energies
      if ((method.eq.0 .AND. &
              adpClosest.lt.myoNucSpacer+nucLen/2 .AND. &
              piClosest.lt.adpPiSpacer+nucLen+nucLen/2) .OR. &
          (method.ne.0 .AND. &
              adpClosest.le.myoNucSpacer .AND. & 
              piClosest.le.adpPiSpacer+nucLen)) then
!       if foot binds ATP or ADP+Pi
         lambda = actMyoBarrierLow
      else if (method.eq.0 .OR. (method.ne.0 .AND. &
                                    (adpClosest.ge.myoNucSpacer+nucLen .OR. &
                                    piClosest.ge.adpPiSpacer+nucLen*2))) then
         lambda = actMyoBarrierHigh
      else
!	LD combination
         if (adpClosest.lt.myoNucSpacer) then
            lambda = (piClosest-(adpPiSpacer+nucLen))/nucLen
         else if (piClosest.lt.adpPiSpacer+nucLen) then
            lambda = (adpClosest-myoNucSpacer)/nucLen
         else
            lambda = 1-   (1-(piClosest-(adpPiSpacer+nucLen))/nucLen) &
                         * (1-(adpClosest-myoNucSpacer)/nucLen)
         end if
!	use lambda to get the combined barrier
!	notice that the combination generates a barrier that is 1/2 at the threshold point
         lambda = actMyoBarrierLow*(1-lambda) + actMyoBarrierHigh * lambda
      end if
!       TODO rubbish again
      if (nucRegime.eq.0) then
!       fix barrier to be high in case this foot is trying to detatch while other foot is already detached
         if (actMyoDist.eq.actMyo1) then
!       if myo1 is bound an myo2 is unbound, force high barrier
            if (actMyo1.lt.actMyoSpacer+actMyoLength/2 .AND. actMyo2.gt.actMyoSpacer+actMyoLength/2) then
               lambda = actMyoBarrierHigh
            end if
         else if (actMyoDist.eq.actMyo2) then
            if (actMyo2.lt.actMyoSpacer+actMyoLength/2 .AND. actMyo1.gt.actMyoSpacer+actMyoLength/2) then
               lambda = actMyoBarrierHigh
            end if
         end if
      end if
!	NOTE
!	the energy function I am using is going to be the markus term, m(d), for the energy profile (barrier and attraction to the actin
!	then the energy to keep the foot aligned to the xyz coordinate system of the actin
!	and the coupling term f(d)
!	the alignment energy is simply Ua = k*(alpha**2+beta**2+gamma**2) where alpha, beta and gamma are:
!	   the angles between the x'y', x'z', and y'z' planes (of the foot particle) and the X Y Z axis of the system (absolute coordinates)
!	then the coupling term is f(d) = (spacer-d)/(length/2)+1 (and =1 below spacer and =0 above spacer+length/2)
!	so the energy function is m(d)+f(d)*Ua
!	the derivative of m(d) is already calculated by markus subroutine, and the dAngle/dxyz is calculated by alignPlane subroutine
!	so I only need to explicitly code the interaction terms
!	m(d)
      e = markus(actMyoDist, actMyoSpacer, actMyoSpacer+actMyoLength, actMyoBind, lambda, 0.0, 1.75*stepSize,.true.)
!	f(d)
      if (actMyoDist.lt.actMyoSpacer) then 
         fd = 1
      else if (actMyoDist.lt.actMyoSpacer+couplingRatio*actMyoLength) then
         fd = (actMyoSpacer-actMyoDist)/(couplingRatio*actMyoLength) + 1
      else 
         fd = 0
      end if
!	Ua only if f(d) is not equal to 0
      if (actMyoDist.lt.actMyoSpacer+couplingRatio*actMyoLength) then
!	k has to be multiplied by f(d) to get the correct energy/force
         i = domain(indices%domain)%center
         df(1) = 0
         df(2) = 0
         df(3) = 1
         call alignPlane(fd*ka*couplingMultiplyer,0.0,myosin(i+2),myosin(i),myosin(i+1),df,alpha)
         Ua = couplingMultiplyer*ka*alpha**2         
         df(1) = 0
         df(2) = 1
         df(3) = 0
!	reuse alpha for beta
         call alignPlane(fd*ka*couplingMultiplyer,0.0,myosin(i+1),myosin(i),myosin(i+3),df,alpha)
         Ua = Ua + couplingMultiplyer*ka*alpha**2
         df(1) = 1
         df(2) = 0
         df(3) = 0
!	resuse alpha for gamma
         call alignPlane(fd*ka*couplingMultiplyer,0.0,myosin(i+3),myosin(i),myosin(i+2),df,alpha)
         Ua = Ua + couplingMultiplyer*ka*alpha**2
      end if
!	note that m(d) and Ua get the energy or the force based on method
      if (method.eq.0) then
!	add energy of m(d)
         energy = energy + e
!	energy of Ua is already incorporated in the subroutine
!	TODO double check that statement above
      else
!	apply force to myosin
!	force of dm(d)
         call vector(myosin(domain(indices%domain)%center),actin(myosin(domain(indices%domain)%center)%bound),df)
         df = -e*df/actMyoDist
         call addForce(myosin(domain(indices%domain)%center),df)
!	force of f(d)*dUa is incorporated in subroutine
!	force of df(d)*Ua
!	use fd for the derivative
         if (actMyoDist.lt.actMyoSpacer) then
            fd = 0
         else if (actMyoDist.lt.actMyoSpacer+couplingRatio*actMyoLength) then
            fd = -1/(couplingRatio*actMyoLength)
         else
            fd = 0
         end if             
         call vector(myosin(domain(indices%domain)%center),actin(myosin(domain(indices%domain)%center)%bound),df) 
!	time dd/dr time Ua
         df = -fd*df*Ua/actMyoDist
         call addForce(myosin(domain(indices%domain)%center),df)
      end if

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ ADP myosin interactions @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!	between ADP and myosin
      do i=1,2  
!	assign pointer and value
         if (i.eq.1) then
            iP = indices%nuc1
            mDist = myoAdpDist
         else
            iP = indices%nuc2
            mDist = myoAdpDist2
         end if
         if (actMyo1.gt.actMyoSpacer+actMyoLength .AND. actMyo2.gt.actMyoSpacer+actMyoLength) then
!       if both are myosins are unbound, kill program. myosin has detached
            open(unit=54322,file="kill",status='unknown')
            close(54322)
         end if
!	based on regime, establish whether a nucleotide should be able to bind/release
         if (nucRegime.eq.0) then
!       establish whether both myosins are bound to actin
            if (actMyo1.lt.actMyoSpacer+actMyoLength/2 &
                .AND. actMyo2.lt.actMyoSpacer+actMyoLength/2) then
!       both are bound, allow nucleotides to bind/unbind
               lambda = 0
            else
!       one is bound, continue diffusion
               lambda = 1
            end if
         else if (((footTheta.lt.thetaFootDown .AND. nucRegime.eq.-1) &
                 .OR. (footTheta.gt.thetaFootUp .AND. nucRegime.eq.1) ) &
                .AND. ( &
             (method.eq.0 .AND. actMyoDist.lt.actMyoSpacer+actMyoLength/2) .OR. &
             (method.ne.0 .AND. actMyoDist.le.actMyoSpacer))) then
!	foot is bound at the down angle
            lambda = 0
         else if (method.eq.0 .OR. &
                 (method.ne.0 .AND. &
                    (actMyoDist.ge.actMyoSpacer+actMyoLength .OR. & 
                    ((footTheta.ge.thetaFootDown+footSpacer .AND. nucRegime.eq.-1) .OR. &
                     (footTheta.le.thetaFootUp-footSpacer .AND. nucRegime.eq.1))))) then
            lambda = 1
         else
!       LD combination
            if (actMyoDist.lt.actMyoSpacer) then
               if (nucRegime.eq.-1) then
                  lambda = (footTheta-thetaFootDown)/footSpacer
               else
                  lambda = (thetaFootUp-footTheta)/footSpacer
               end if
            else if ((footTheta.lt.thetaFootDown .AND. nucRegime.eq.-1) .OR. &
                     (footTheta.gt.thetaFootUp .AND. nucRegime.eq.1)) then
               lambda = (actMyoDist-actMyoSpacer)/actMyoLength
            else
               if (nucRegime.eq.-1) then
                  lambda = 1-(1-(footTheta-thetaFootDown)/footSpacer)*(1-(actMyoDist-actMyoSpacer)/actMyoLength)
               else
                  lambda = 1-(1-(thetaFootUp-footTheta)/footSpacer)*(1-(actMyoDist-actMyoSpacer)/actMyoLength)
               end if
            end if
         end if
         lambda = myoNucLowBarrier*(1-lambda) + myoNucHighBarrier * lambda
         e = markus(mDist, myoNucSpacer, myoNucSpacer+nucLen, myoNucBind, lambda, 0.0, 0.175*stepSize,.true.)
         if (method.eq.0) then
!       add energy
            energy = energy + e
         else
!       apply force to adp
            adp(iP)%xF = adp(iP)%xF - e * (myosin(domain(indices%domain)%center)%x-adp(iP)%x)/mDist
            adp(iP)%yF = adp(iP)%yF - e * (myosin(domain(indices%domain)%center)%y-adp(iP)%y)/mDist
            adp(iP)%zF = adp(iP)%zF - e * (myosin(domain(indices%domain)%center)%z-adp(iP)%z)/mDist
         end if
!	add boundary for ADP and myosin
         if (mDist.gt.nucleotideBoundary) then
            if (method.eq.0) then
               energy = energy + nucleotideBoundaryK*(mDist-nucleotideBoundary)**2
            else
!	recycle e
               e = 2*nucleotideBoundaryK*(mDist-nucleotideBoundary)*kcal2cu
               adp(iP)%xF = adp(iP)%xF + e * (myosin(domain(indices%domain)%center)%x-adp(iP)%x)/mDist
               adp(iP)%yF = adp(iP)%yF + e * (myosin(domain(indices%domain)%center)%y-adp(iP)%y)/mDist
               adp(iP)%zF = adp(iP)%zF + e * (myosin(domain(indices%domain)%center)%z-adp(iP)%z)/mDist         
            end if
         end if
      end do

! @@@@@@@@@@@@@@@@@@@@@@@@@
! @@ ADP Pi interactions @@
! @@@@@@@@@@@@@@@@@@@@@@@@@
!	ADP Pi have 2 barriers. One for the 'bond' and one for the Pi release
      do i = 1,2
!       assign pointer and value
         if (i.eq.1) then
               iP = indices%nuc1
               mDist = myoAdpDist
               pDist = adpPDist
         else
               iP = indices%nuc2
               mDist = myoAdpDist2
               pDist = adpPDist2
         end if
         if ((method.eq.0 .AND. & 
                 mDist.lt.myoNucSpacer+nucLen/2 .AND. actMyoDist.ge.actMyoSpacer+actMyoLength/2) .OR. &
             (method.ne.0 .AND. &
                 mDist.le.myoNucSpacer .AND. actMyoDist.ge.actMyoSpacer+actMyoLength) ) then
!       Nuc is bound to foot and foot is not bound to actin - hydrolysis but no release:
            lambda = adpPiBarrierP
            lambda2 = adpPi2ndBarrierW
         else if ((method.eq.0 .AND. &
                      mDist.lt.myoNucSpacer+nucLen/2 .AND. actMyoDist.lt.actMyoSpacer+actMyoLength/2) .OR. &
                  (method.ne.0 .AND. &
                      mDist.le.myoNucSpacer .AND. actMyoDist.le.actMyoSpacer) ) then
!       ADP+P is bound to foot and foot is bound to actin - release but not hydrolysis
            lambda = adpPiBarrierW
            lambda2 = adpPi2ndBarrierP
         else if (method.eq.0 .OR. &
                 (method.ne.0 .AND. mDist.ge.myoNucSpacer+nucLen) ) then
!       Free nucleotide - no release and no hydrolysis
            lambda = adpPiBarrierW
            lambda2 = adpPi2ndBarrierW
!       combinations
         else 
            if (mDist.le.myoNucSpacer) then
               lambda  = 1 - (actMyoDist - actMyoSpacer)/actMyoLength
               lambda2 = (actMyoDist - actMyoSpacer)/actMyoLength
            else if (actMyoDist.le.actMyoSpacer) then
               lambda  = 1
               lambda2 = (mDist-myoNucSpacer)/nucLen
            else if (actMyoDist.ge.actMyoSpacer+actMyoLength) then
               lambda  = (mDist-myoNucSpacer)/nucLen
               lambda2 = 1
            else
               lambda  = 1-((actMyoDist - actMyoSpacer)/actMyoLength)*(1-(mDist-myoNucSpacer)/nucLen)
               lambda2 = 1-(1-(actMyoDist - actMyoSpacer)/actMyoLength)*(1-(mDist-myoNucSpacer)/nucLen)
            end if
!	calculated weighted barriers then energy
            lambda  = adpPiBarrierP *(1-lambda)  + adpPiBarrierW  * lambda
            lambda2 = adpPi2ndBarrierP*(1-lambda2) + adpPi2ndBarrierW * lambda2
         end if
!       bond energy
            e = markus(pDist, adpPiSpacer, adpPiSpacer+nucLen, adpPiBind, lambda, 0.0, 0.175*stepSize,.true.) &
!       release energy
              + markus(pDist, adpPiSpacer+nucLen, adpPiSpacer+nucLen*2, adpPi2ndBind, lambda2, 0.0, 0.175*stepSize,.true.)
         if (method.eq.0) then
!       add energy
            energy = energy + e
         else
!       apply force to pi
            call vector (pi(iP),adp(iP),df)
            df = e*df/pDist
            call addForce(pi(iP),-df)
         end if
!       add boundary for ADP and Pi
         if (pDist.gt.phosphateBoundary) then
            if (method.eq.0) then
               energy = energy + nucleotideBoundaryK*(pDist-phosphateBoundary)**2
            else
!       recycle e
               e = 2*nucleotideBoundaryK*(pDist-phosphateBoundary)*kcal2cu
               call vector(pi(iP),adp(iP),df)
               df = e*df/pDist
               call addForce(pi(iP),df)
            end if
         end if
      end do

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ Rejection between ADPs @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@
!	threshold is the max radius for relevent interactions between the nuc and the myo (so its the diameter)
      if (interADP.lt.(myoNucSpacer+nucLen)) then
         if (method.eq.0) then
            energy = energy + (nucleotideBoundaryK/2)*(interADP-(myoNucSpacer+nucLen))**2
         else
!	use half K to make the force less violent
            e = nucleotideBoundaryK*(interADP-(myoNucSpacer+nucLen))*kcal2cu/interADP
            call vector(adp(indices%nuc1),adp(indices%nuc2),df)
            df = e*df
!	apply force only to the unbound nucleotide
            if (myoAdpDist.gt.myoAdpDist2) then
               call addForce(adp(indices%nuc1),df)
            else 
               call addForce(adp(indices%nuc2),-df)   
            end if
         end if
      end if

! @@@@@@@@@@@@@@@@@@@@@@@
! @@ foot angle energy @@
! @@@@@@@@@@@@@@@@@@@@@@@
!	Foot angle energy
      e = footAngleEnergy(footTheta, thetaFootDown, thetaFootUp,kt, footDG, footBarrier, 0.0, 0.02*stepSize)
      if (method.eq.0) then
         energy = energy + e
      else
         myo(1) = domain(indices%knee)%center
         myo(2) = domain(indices%domain)%center
         myo(3) = domain(indices%domain)%pY
         myo(4) = domain(indices%domain)%pZ
         call getTorsion(kt,0.0,myosin(myo(1)),myosin(myo(2)),myosin(myo(3)),myosin(myo(4)),2,e)
      end if 
!	ensure the foot stays in the xz plane
      myo(1) = domain(indices%knee)%center
      myo(2) = domain(indices%domain)%center
      myo(3) = domain(indices%domain)%pX
      myo(4) = domain(indices%domain)%pZ
      call getTorsion(kt,0.0,myosin(myo(1)),myosin(myo(2)),myosin(myo(3)),myosin(myo(4)))
!	energy will be inserted in the myosin
      myosin(domain(indices%domain)%center)%eTotal = myosin(domain(indices%domain)%center)%eTotal + energy

! @@@@@@@@@@@@@@@@@@@@@
! @@ Histogram stuff @@
! @@@@@@@@@@@@@@@@@@@@@
!	determine which particle this is
      if (printOutput) then
         if (indices%domain.eq.feet(1)%domain) then
            floored = floor(adpPDist)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.100) then
               floored = 100
            end if
            PiDistance(1,floored) = PiDistance(1,floored)+1
            floored = floor(adpPDist2)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.100) then
               floored = 100
            end if
            PiDistance(2,floored) = PiDistance(2,floored)+1
            floored = floor(myoAdpDist)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.100) then
               floored = 100
            end if
            ADPDistance(1,floored) = ADPDistance(1,floored)+1
            floored = floor(myoAdpDist2)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.100) then
               floored = 100
            end if
            ADPDistance(2,floored) = ADPDistance(2,floored)+1
            floored = floor(actMyoDist)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.900) then
               floored = 900
            end if
            MyoDistance(1,floored) = MyoDistance(1,floored)+1
            floored = floor(footTheta*100)
            if (floored.le.-350) then
               floored = -350
            else if (floored.ge.350) then
               floored = 350
            end if
            MyoAngle(1,floored) = MyoAngle(1,floored)+1
         else
            floored = floor(adpPDist)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.100) then
               floored = 100
            end if
            PiDistance(3,floored) = PiDistance(3,floored)+1
            floored = floor(adpPDist2)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.100) then
               floored = 100
            end if
            PiDistance(4,floored) = PiDistance(4,floored)+1
            floored = floor(myoAdpDist)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.100) then
               floored = 100
            end if
            ADPDistance(3,floored) = ADPDistance(3,floored)+1
            floored = floor(myoAdpDist2)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.100) then
               floored = 100
            end if
            ADPDistance(4,floored) = ADPDistance(4,floored)+1
            floored = floor(actMyoDist)
            if (floored.le.0) then
               floored = 0
            else if (floored.ge.900) then
               floored = 900
            end if
            MyoDistance(2,floored) = MyoDistance(2,floored)+1
            floored = floor(footTheta*100)
            if (floored.le.-350) then
               floored = -350
            else if (floored.ge.350) then
               floored = 350
            end if
            MyoAngle(2,floored) = MyoAngle(2,floored)+1
         end if
      end if

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ Apply a linear force to help the myosin find the actin @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine biasDiffusion(myosin, adp, pi)

      implicit none

      type(Particle), dimension(:), intent(inout) :: myosin,adp,pi
      double precision, dimension(2) :: distances
!	indices will be presaved to make code more redable
!	index i and index j is the opposite index
      integer :: i, j, myo(2), act(2)
      double precision, dimension(3) :: vec
      double precision :: f

!	get indices and calculate distances
      do i = 1,2
         myo(i) = domain(feet(i)%domain)%center
         act(i) = myosin(myo(i))%bound
         distances(i) = distance(myosin(myo(i)),actin(act(i)))
      end do

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ Check if an ATP particle needs to be reconstructed to ADP @@ 
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine checkATP(indices,myosin,adp,pi)
!	MC version

      implicit none

      type(Particle), dimension(:), intent(inout) :: myosin,adp,pi
      type(Foot) :: indices
      integer ::  a, b, c;

      a = indices%nuc1
      b = indices%nuc2
      c = domain(indices%domain)%center

!	determine whether a reconstruction is needed
      if ((distance(adp(a),myosin(c)).gt.myoNucSpacer+nucLen)*   &
         ((distance(adp(b),myosin(c)).gt.myoNucSpacer+nucLen))) then
!	recycle c, get a random number from 0 to 9
         c = (mod(floor(roll*1000),10))         
         if ((distance(adp(a),pi(a)).lt.adpPiSpacer+nucLen/2)*(distance(adp(b),pi(b)).lt.adpPiSpacer+nucLen/2)) then
!       both are ATP:
!	transform ATP to ADP
            if (c.lt.5) then
               pi(a)%x = adp(a)%x+50
               pi(a)%y = adp(a)%y
               pi(a)%z = adp(a)%z
            else
               pi(b)%x = adp(b)%x+50
               pi(b)%y = adp(b)%y
               pi(b)%z = adp(b)%z
            end if
         else if ((distance(adp(a),pi(a)).gt.adpPiSpacer+nucLen)*(distance(adp(b),pi(b)).gt.adpPiSpacer+nucLen)) then
!	both are ADP: transform to ATP
            if (c.lt.5) then
               pi(a)%x = adp(a)%x
               pi(a)%y = adp(a)%y
               pi(a)%z = adp(a)%z
            else
               pi(b)%x = adp(b)%x
               pi(b)%y = adp(b)%y
               pi(b)%z = adp(b)%z
            end if
         end if
!	recalculate energy for the next step
         call applyFF(myosin,adp,pi)
         previousEnergy = totalEnergy
      end if

      end subroutine

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ Check if an ATP particle needs to be reconstructed to ADP @@ 
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine checkATPLD(myosin,adp,pi)
!	LD version

      implicit none

      type(Particle), dimension(:), intent(inout) :: myosin,adp,pi

!	index of the ADP moved
!	a and b are the adp's. c is the myosin head
      integer ::  a, b, c
!	index
      integer :: i
!	saved distances
      double precision :: distA, distB

      do i = 1,2
         a = feet(i)%nuc1
         b = feet(i)%nuc2
         c = domain(feet(i)%domain)%center
!	determine whether both nucleotides are unbound
         distA = distance(adp(a),myosin(c))
         distB = distance(adp(b),myosin(c))
         if ((distA.gt.myoNucSpacer+nucLen)*   &
            ((distB.gt.myoNucSpacer+nucLen))) then
            if ((distance(adp(a),pi(a)).lt.adpPiSpacer+nucLen/2)*(distance(adp(b),pi(b)).lt.adpPiSpacer+nucLen/2)) then
!       both are ATP:
!	transform ATP to ADP, the farthest
!	+1 to avoid overflows with distance 0 or something
               if (distA.gt.distB) then
                  pi(a)%x = adp(a)%x+phosphateBoundary-nucLen
                  pi(a)%y = adp(a)%y+1
                  pi(a)%z = adp(a)%z+1
                  pi(a)%xV = 0
                  pi(a)%yV = 0
                  pi(a)%zV = 0
               else
                  pi(b)%x = adp(b)%x+phosphateBoundary-nucLen
                  pi(b)%y = adp(b)%y+1
                  pi(b)%z = adp(b)%z+1
                  pi(b)%xV = 0
                  pi(b)%yV = 0
                  pi(b)%zV = 0
               end if
            else if ((distance(adp(a),pi(a)).gt.adpPiSpacer+nucLen)*(distance(adp(b),pi(b)).gt.adpPiSpacer+nucLen)) then
!	both are ADP: transform to ATP
               if (distA.gt.distB) then
                  pi(a)%x = adp(a)%x+1
                  pi(a)%y = adp(a)%y+1
                  pi(a)%z = adp(a)%z+1
                  pi(a)%xV = 0
                  pi(a)%yV = 0
                  pi(a)%zV = 0
               else
                  pi(b)%x = adp(b)%x+1
                  pi(b)%y = adp(b)%y+1
                  pi(b)%z = adp(b)%z+1
                  pi(b)%xV = 0
                  pi(b)%yV = 0
                  pi(b)%zV = 0
               end if
            end if
         end if
      end do

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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Print out potential curves @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine printOutLandscapes

      implicit none

      integer :: i,j
      double precision :: energy, x, y, e

      open(301,file='myoActATPorADPPi.ene',status='unknown')
      open(302,file='myoActElse.ene',status='unknown')
      open(304,file='myoNucATPorDownBound.ene',status='unknown')
      open(305,file='myoNucElse.ene',status='unknown')
      open(306,file='nucPLooseFoot.ene',status='unknown')
      open(307,file='nucPBoundFoot.ene',status='unknown')
      open(308,file='nucPFree.ene',status='unknown')
      open(309,file='footAngle.ene',status='unknown')

      do i = 0,2500
         x = i/10.d0
         e = markus(x, actMyoSpacer, actMyoSpacer+actMyoLength, actMyoBind, 1.d0*actMyoBarrierLow, 0.0, 1.75*stepSize,.true.)
         write(301,'(f8.3,f8.3)') x, e
         e = markus(x, actMyoSpacer, actMyoSpacer+actMyoLength, actMyoBind, 1.d0*actMyoBarrierHigh, 0.0, 1.75*stepSize,.true.)
         write(302,'(f8.3,f8.3)') x, e
         x = i/49.d0
         energy = markus(x, myoNucSpacer, myoNucSpacer+nucLen, myoNucBind, 1.d0*myoNucLowBarrier, 0.0, 0.175*stepSize,.true.)
         if (x.gt.nucleotideBoundary) then
            energy = energy +nucleotideBoundaryK*(x-nucleotideBoundary)**2
         end if
         write(304,'(f8.3,f8.3)') x, energy
         energy = markus(x, myoNucSpacer, myoNucSpacer+nucLen, myoNucBind, 1.d0*myoNucHighBarrier, 0.0, 0.175*stepSize,.true.)
         if (x.gt.nucleotideBoundary) then
            energy = energy +nucleotideBoundaryK*(x-nucleotideBoundary)**2
         end if
         write(305,'(f8.3,f8.3)') x, energy 
         energy = markus(x, adpPiSpacer, adpPiSpacer+nucLen, adpPiBind, 1.d0*adpPiBarrierP, 0.0, 0.175*stepSize,.true.) &
                + markus(x, adpPiSpacer+nucLen, adpPiSpacer+nucLen*2, adpPi2ndBind, 1.d0*adpPi2ndBarrierW, 0.0, 0.175*stepSize,.true.)
         if (x.gt.phosphateBoundary) then
            energy = energy +nucleotideBoundaryK*(x-phosphateBoundary)**2
         end if
         write(306,'(f8.3,f8.3)') x, energy
         energy = markus(x, adpPiSpacer, adpPiSpacer+nucLen, adpPiBind, 1.d0*adpPiBarrierW, 0.0, 0.175*stepSize,.true.) &
                + markus(x, adpPiSpacer+nucLen, adpPiSpacer+nucLen*2, adpPi2ndBind, 1.d0*adpPi2ndBarrierP, 0.0, 0.175*stepSize,.true.)
         if (x.gt.phosphateBoundary) then
            energy = energy +nucleotideBoundaryK*(x-phosphateBoundary)**2
         end if
         write(307,'(f8.3,f8.3)') x, energy
         energy = markus(x, adpPiSpacer, adpPiSpacer+nucLen, adpPiBind, 1.d0*adpPiBarrierW, 0.0, 0.175*stepSize,.true.) &
                + markus(x, adpPiSpacer+nucLen, adpPiSpacer+nucLen*2, adpPi2ndBind, 1.d0*adpPi2ndBarrierW, 0.0, 0.175*stepSize,.true.)
         if (x.gt.phosphateBoundary) then
            energy = energy +nucleotideBoundaryK*(x-phosphateBoundary)**2
         end if
         write(308,'(f8.3,f8.3)') x, energy
         x = i/600.d0
      end do
      do i = -4000,4000
         x = i/1000.0
         write(309,'(f10.3,f10.3)') x, footAngleEnergy(x, thetaFootDown, thetaFootUp, kt, footDG, footBarrier, 0.0, 0.02*stepSize)
      end do
      
      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Print out histograms @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine printHistograms

      implicit none

      open(301,file='pi1.histo',status='unknown')
      open(302,file='pi2.histo',status='unknown')
      open(303,file='pi3.histo',status='unknown')
      open(304,file='pi4.histo',status='unknown')
      open(401,file='adp1.histo',status='unknown')
      open(402,file='adp2.histo',status='unknown')
      open(403,file='adp3.histo',status='unknown')
      open(404,file='adp4.histo',status='unknown')
      open(501,file='myo1.histo',status='unknown')
      open(502,file='myo2.histo',status='unknown')
      open(601,file='angle1.histo',status='unknown')
      open(602,file='angle2.histo',status='unknown')

      do i = 0,100
         write(301,*) i, PiDistance(1,i)
         write(302,*) i, PiDistance(2,i)
         write(303,*) i, PiDistance(3,i)
         write(304,*) i, PiDistance(4,i)
         write(401,*) i, ADPDistance(1,i)
         write(402,*) i, ADPDistance(2,i)
         write(403,*) i, ADPDistance(3,i)
         write(404,*) i, ADPDistance(4,i)
      end do

      do i = 0,900
         write(501,*) i, MyoDistance(1,i)
         write(502,*) i, MyoDistance(2,i)
      end do

      do i = -350,350
         write(601,*) i/100.d0, MyoAngle(1,i)
         write(602,*) i/100.d0, MyoAngle(2,i)
      end do

!	close files
      close(301)
      close(302)
      close(303)
      close(304)
      close(401)
      close(402)
      close(403)
      close(404)
      close(501)
      close(502)
      close(601)
      close(602)

      end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@ Precalculate sigma's @@
!@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine calculateSigma(myosin,adp,p)

      implicit none

      integer :: i

      type(Particle), dimension(:), intent(inout) :: myosin,adp,p

!       myosin particles
      if (method.eq.1) then
         myosin%sigma = sqrt(2*bz_cu*temperature*dt/(myosin%mass*friction))
      else
         myosin%sigma = sqrt(2*bz_cu*temperature*myosin%mass*friction/dt)
      end if
!       adp and pi
      if (method.eq.1) then
         adp%sigma = sqrt(2*bz_cu*temperature*dt/(adp%mass*friction))
         p%sigma   = sqrt(2*bz_cu*temperature*dt/(p%mass*friction))
      else
         adp%sigma = sqrt(2*bz_cu*temperature*adp%mass*friction/dt)
         p%sigma   = sqrt(2*bz_cu*temperature*p%mass*friction/dt)
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

!@@@@@@@@@@@@@@@@@
!@@ State stuff @@
!@@@@@@@@@@@@@@@@@
      subroutine checkState(myoAngleS, myoBindS, myoNucS)

      implicit none

      integer, dimension(2), intent(INOUT) :: myoAngleS, myoBindS, myoNucS
      integer :: f, i, j, n
      integer, dimension(4) :: indices
!	temp value 
      double precision :: dist, distP
!	1 for adp1, 2 for adp2, 3 adp-pi1, 4 adp-pi2
      double precision, dimension(4) :: innerDist
!	flag for extra printout
      logical :: pdb

      pdb = .false.

!	myosin binding to actin
      do f = 1,2
         i = domain(feet(f)%domain)%center
         dist = distance(myosin(i),actin(myosin(i)%bound))
         if (myoBindS(f).eq.0 .AND. dist.lt.actMyoSpacer) then
            write(202,*), iStep, "ACT:Myo",f , "now bound to ", myosin(i)%bound
            pdb = .true.
            myoBindS(f) = 1
         else if (myoBindS(f).eq.1 .AND. dist.gt.actMyoSpacer+actMyoLength) then
            write(202,*), iStep, "ACT:Myo",f ," is now unbound"
            pdb = .true.
            myoBindS(f) = 0
!	flagged for initiation
         else if (myoBindS(f).eq.-1) then
            if (dist.lt.actMyoSpacer+actMyoLength/2) then
               myoBindS(f) = 1
            else
               myoBindS(f) = 0
            end if
         end if
      end do

!       myosin foot angle
      do f = 1,2
         indices(1) = domain(feet(f)%knee)%center
         indices(2) = domain(feet(f)%domain)%center
         indices(3) = domain(feet(f)%domain)%pY
         indices(4) = domain(feet(f)%domain)%pZ
         call getTorsion(0.0,0.0,myosin(indices(1)),myosin(indices(2)),myosin(indices(3)),myosin(indices(4)),1,dist)
         if (myoAngleS(f).eq.-1 .AND. dist.gt.thetaFootUp - footSpacer) then
            write(202,*), iStep, "ANG:Myo",f , " is now angle Up"
            pdb = .true.
            myoAngleS(f) = 1
         else if (myoAngleS(f).eq.1 .AND. dist.lt.thetaFootDown + footSpacer) then
            write(202,*), iStep, "ANG:Myo",f ," is now angle Down"
            pdb = .true.
            myoAngleS(f) = -1
         else if (myoAngleS(f).eq.0) then
            if (dist.lt.(thetaFootUp+thetaFootDown)/2) then
               myoAngleS(f) = -1
            else
               myoAngleS(f) = 1
            end if
         end if
      end do

!       myosin binding to nucleotide
      do f = 1,2
         i = domain(feet(f)%domain)%center
         do n = 1,2
            if (n.eq.1) then
               j = feet(f)%nuc1
               innerDist(1) = distance(myosin(i),adp(j))
               innerDist(3) = distance(adp(j),pi(j))
            else
               j = feet(f)%nuc2
               innerDist(2) = distance(myosin(i),adp(j))
               innerDist(4) = distance(adp(j),pi(j))
            end if
         end do
!	mark closest
         if (innerDist(1).lt.innerDist(2)) then
            dist = innerDist(1)
            distP = innerDist(3)
         else
            dist = innerDist(2)
            distP = innerDist(4)
         end if
!       myoNuc = 0 nothing, 2 adp, 3 atp, 4 adp+p
         if (myoNucS(f).eq.0) then
            if (dist.lt.myoNucSpacer .AND. distP.lt.adpPiSpacer+nucLen/2) then
               write(202,*), iStep, "NUC:Myo",f , "now bound to ATP"
               pdb = .true.
               myoNucS(f) = 3
            else if (dist.lt.myoNucSpacer .AND. distP.gt.adpPiSpacer+(3*nucLen)/2) then
               write(202,*), iStep, "NUC:Myo",f , "now bound to ADP"
               pdb = .true.
               myoNucS(f) = 2
            else if (dist.lt.myoNucSpacer) then
               write(202,*), iStep, "NUC:Myo",f , "now bound to ADP+P"
               pdb = .true.
               myoNucS(f) = 4       
            end if
         else if (myoNucS(f).gt.0 .AND. dist.gt.myoNucSpacer+nucLen) then
            write(202,*), iStep, "NUC:Myo",f , "released the nucleotide"
            pdb = .true.
            myoNucS(f) = 0
         else if (myoNucS(f).eq.2 .AND. distP.lt.myoNucSpacer+nucLen) then
            write(202,*), iStep, "NUC:Myo",f , "has bound a phosphate"
            pdb = .true.
            myoNucS(f) = 4
         else if (myoNucS(f).eq.4 .AND. distP.lt.myoNucSpacer) then
            write(202,*), iStep, "NUC:Myo",f , "has created ATP"
            pdb = .true.
            myoNucS(f) = 3            
         else if (myoNucS(f).eq.3 .AND. distP.gt.myoNucSpacer+nucLen) then
            write(202,*), iStep, "NUC:Myo",f , "has hydrolyzed ATP"
            pdb = .true.
            myoNucS(f) = 4
         else if (myoNucS(f).eq.4 .AND. distP.gt.myoNucSpacer+2*nucLen) then
            write(202,*), iStep, "NUC:Myo",f , "has released a phosphate"
            pdb = .true.
            myoNucS(f) = 2
!	initialize
         else if (myoNucS(f).eq.-1) then
            if (dist.gt.myoNucSpacer+nucLen/2) then
               myoNucS(f) = 0
            else if (distP.lt.myoNucSpacer+nucLen/2) then
               myoNucS(f) = 3
            else if (distP.gt.myoNucSpacer+(3*nucLen)/2) then
               myoNucS(f) = 2
            else 
               myoNucS(f) = 4              
            end if
         end if
      end do

!	print special pdb
      if (pdb) then
         write(333,'(''MODEL'',i12)') iStep
         do n=1,pLast
            write(333,pdbFormat) n,'OH  ','MYO',n,myosin(n)%x,myosin(n)%y,myosin(n)%z
         enddo
         do i=1,4
            write(333,pdbFormat) i*2-1+pLast,'NA  ','ADP',i,adp(i)%x,adp(i)%y,adp(i)%z
            write(333,pdbFormat) i*2+pLast,'P   ','PHO',i,pi(i)%x,pi(i)%y,pi(i)%z
         end do
         write(333,'(''TER'')')
         write(333,'(''ENDMDL'',i12)') istep
      end if

      end subroutine

      end program
