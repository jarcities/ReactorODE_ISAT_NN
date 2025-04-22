! filepath: /home/jay/Code/reactorode_isat_nn/driver_samples.f90
!-----------------------------------------------------------------------
! Module definition for external Fortran interfaces
module extFGH

    interface
       !--------------------------------------------------------------------
       ! Declaration of the C-bound subroutine "myfgh":
       ! - Evaluates the residual function (and optionally Jacobian) for the reactor model.
       ! - Arguments: need (flag array), nx (size of state vector),
       !   x (state vector), nf, nh (dimensions), iusr and rusr (user controlled arrays),
       !   and outputs f, g, h.
       subroutine myfgh ( need, nx, x, nf, nh, iusr, rusr, f, g, h ) bind ( c )
          use iso_c_binding
          integer (c_int) :: need(*)  ! flag array (e.g., whether to compute Jacobian)
          integer (c_int) :: nx       ! number of state variables
          integer (c_int) :: nf       ! size of output f (residual)
          integer (c_int) :: nh       ! size of h
          integer (c_int) :: iusr(*)  ! integer user array
          real (c_double) :: x(*)     ! input state vector (transformed)
          real (c_double) :: rusr(*)  ! double user parameters (includes tolerances, etc.)
          real (c_double) :: f(*)     ! output residual vector
          real (c_double) :: g(*)     ! output Jacobian (if requested)
          real (c_double) :: h(*)     ! additional output
       end subroutine myfgh
 
       !--------------------------------------------------------------------
       ! Declaration of the C-bound subroutine "mymix":
       ! - Mixes two state vectors (x1 and x2) using a specified mixing parameter.
       subroutine mymix ( nx, x1, x2, alpha, iusr, rusr ) bind ( c )
          use iso_c_binding
          integer (c_int) :: nx      ! number of state variables
          integer (c_int) :: iusr(*) ! user integer array
          real (c_double) :: x1(*)   ! first state vector
          real (c_double) :: x2(*)   ! second state vector
          real (c_double) :: rusr(*) ! user double array
          real (c_double) :: alpha(*)! mixing parameter (e.g., dt/mixTime)
       end subroutine mymix
 
       !--------------------------------------------------------------------
       ! Declaration of the C-bound subroutine "toxhat":
       ! - Transforms a physical state vector (ptcl) into the internal (x) representation.
       subroutine toxhat ( ptcl, x, nx, rusr ) bind ( c )
          use iso_c_binding
          integer (c_int) :: nx      ! number of state variables
          real (c_double) :: ptcl(*) ! physical state vector (input)
          real (c_double) :: x(*)    ! transformed state vector (output)
          real (c_double) :: rusr(*) ! user parameters for transformation
       end subroutine toxhat
 
       !--------------------------------------------------------------------
       ! Declaration of the C-bound subroutine "myfnn":
       ! - Evaluates the neural network function on the state vector x.
       subroutine myfnn ( nx, x, fnn ) bind ( c )
          use iso_c_binding
          integer (c_int) :: nx      ! number of state variables
          real (c_double) :: x(*)    ! input state vector
          real (c_double) :: fnn(*)  ! output from the neural network
       end subroutine myfnn
 
       !--------------------------------------------------------------------
       ! Declaration of the C-bound subroutine "fromxhat":
       ! - Converts the internal state (x) back to the physical state (ptcl).
       subroutine fromxhat ( x, ptcl, nx, rusr ) bind ( c )
          use iso_c_binding
          integer (c_int) :: nx       ! number of state variables
          real (c_double) :: x(*)     ! input state in x-space
          real (c_double) :: ptcl(*)  ! output physical state vector
          real (c_double) :: rusr(*)  ! user parameters for transformation
       end subroutine fromxhat
 
    end interface
 
 end module extFGH
 !-----------------------------------------------------------------------
 ! End of module "extFGH"
 !-----------------------------------------------------------------------
 
 program main
 
    ! Bring in the interfaces from the module extFGH.
    use extFGH
 
    implicit none
 
    ! Declaration of pointer arrays and scalar variables.
    integer, pointer :: need(:), need2(:), iusr(:)
    integer :: ii, jj, kk
    integer :: nx = 11, nf = 11, nh = 1, ns = 10, info(100), nPtcl = 250, idd, nSteps = 4000, counter, i1, i2
    double precision, pointer :: x(:), x1(:), x2(:), rusr(:), f(:), fisat(:), fisatsum(:), fnn(:)
    double precision, pointer :: g(:,:), h(:), ptcls(:,:), alpha(:)
    double precision :: rinfo(70), stats(100), dt, start, finish
    double precision :: flowThroughTime, mixTime
 
    ! Index identifiers for species (e.g., hydrogen, oxygen, nitrogen)
    integer :: nH2 = 2, nO2 = 5, nN2 = 11
 
    ! Parameters for particle initialization – fuel fraction and temperatures.
    double precision :: mFrac = 0.1, Tfuel = 1050.0, Tox = 1050.0
 
    double precision :: rr  ! variable for random number
 
    ! Allocate arrays with the specified sizes.
    allocate(need(3), need2(3), iusr(3), x(nx), x1(nx), x2(nx), rusr(2*nx+5), f(nx), fisat(nx), fisatsum(nx), fnn(nx), g(nx,nx), h(1), ptcls(nx+1, nPtcl), alpha(1))
 
    ! Set simulation parameters.
    dt = 5e-7                  ! time step for integration
    flowThroughTime = 1e-4       ! characteristic time for particle replacement (flow-through)
    mixTime = 1e-4             ! characteristic mixing time
 
    ! Initialize flag arrays:
    need = (/ 1, 0, 0 /)       ! "need" array; 1 means at least function evaluation is required
    need2 = (/ 1, 1, 0 /)      ! "need2" array; here the second element is 1 meaning request for Jacobian computation
 
    ! Initialize the rusr (user double array) parameters:
    rusr(1:2*nx) = 0.0         ! first 2*nx entries set to 0 (no normalization initially)
    rusr(1:nx+1) = (/ 1050.0, 1e-5, 7.5967e-08, 1.4358e-07, 2.1008e-06, 1.5337e-07, 2.1554e-06, 1.4550e-09, 9.1008e-11, 1.0000e-35, 7.8992e-06, 358.9534 /)
       ! The first nx+1 values include initial conditions (e.g., Tfuel) and scaling factors.
    rusr(2*nx+1) = 1e-8       ! absolute tolerance (atol) used by solvers
    rusr(2*nx+2) = 1e-8       ! relative tolerance (rtol)
    rusr(2*nx+3) = dt         ! pass the time step (dt) to routines
    rusr(2*nx+4) = 1e-6       ! finite difference perturbation (dx)
    rusr(2*nx+5) = 101325.0   ! pressure (p)
 
    ! Initialize ISAT control flags.
    info  = 0                 ! set the info array to 0 initially
    rinfo = 0.0               ! set the rinfo array to 0 initially (for real-valued control parameters)
    info(2)  = 0              ! if_g flag, not used (set to 0)
    info(12) = 2              ! option for ISAT operation (e.g., table operation mode)
    info(28) = 0              ! idites, possibly related to iteration control (set to 0)
    rinfo(1) = 1e-3           ! error tolerance for ISAT (first tolerance)
    rinfo(2) = 1e-3           ! error tolerance for ISAT (second tolerance)
    rinfo(3) = 1e2            ! another parameter (possibly maximum number of entries)
    !-- special settings: Table 2
    !info(21,2) = 1   ! ifull – commented-out special setting
    rinfo(8) = 1000.          ! another ISAT parameter (e.g., maximum number of lookups or similar)
    !-- special settings: Table 3
 
    ! Set the seed for random number generation for reproducibility.
    call random_seed(put=(/10,11/))
 
    ! Record CPU start time.
    call cpu_time(start)
 
    ! Initialize a particle counter.
    counter = 0
 
    !-------------------------------
    ! Initialize the particle states.
    do ii = 1, nPtcl ! Loop over each of the 250 particles
 
       call random_number(rr)  ! Generate a random number between 0 and 1
 
       if ( rr < mFrac ) then
          ! If the random number is below the fuel fraction, set this particle as fuel.
          ptcls(1,ii) = Tfuel            ! Set temperature to Tfuel.
          ptcls(2:nx,ii) = 0.0           ! Set the remainder of the state vector to 0.
          ptcls(nH2,ii) = 1.0            ! Set the hydrogen mass fraction to 1 (fuel inflow).
       else
          ! Otherwise, set the particle as oxidizer.
          ptcls(1,ii) = Tox              ! Set temperature to Tox.
          ptcls(2:nx,ii) = 0.0           ! Set the rest of the state vector to 0.
          ptcls(nO2,ii) = 1.0/4.76       ! Set oxygen fraction.
          ptcls(nN2,ii) = 3.76/4.76      ! Set nitrogen fraction (oxidizer inflow).
       end if
 
       counter = counter + 1            ! Increment the particle counter.
       ptcls(nx+1,ii) = counter         ! Store the particle identifier (ID).
 
    end do ! End of particle initialization.
    !-------------------------------
 
    !-------------------------------
    ! Main time-stepping loop (simulate reactor dynamics for nSteps iterations).
    do jj = 1, nSteps ! Loop over time steps
 
       if (modulo(jj,100).eq.0) then
          ! Every 100 steps, print current particle states.
          print *, ptcls
       end if
 
       !-------------------------------
       ! Reaction step for each particle.
       do ii = 1, nPtcl
          ! Transform the particle's physical state to the internal state representation.
          call toxhat( ptcls(1:nx,ii), x, nx, rusr )
 
          ! Optionally, one could directly call myfgh for a reaction step.
          ! The following call is commented out:
          ! call myfgh( need, nx, x, nf, nh, iusr, rusr, f, g, h )
 
          ! Instead, call the ISAT routine which uses myfgh internally if needed:
          call isatab( idd, 0, nx, x, nf, nh, nh, myfgh, iusr, rusr, info, rinfo, f, g, h, stats)
             ! idd: ISAT table identifier.
             ! 0: possibly a flag indicating a new/ongoing lookup.
             ! nx, x: state vector size and its current value.
             ! nf, nh, nh: sizes for function outputs and additional arrays.
             ! myfgh: pointer to the subroutine for function evaluation.
             ! iusr, rusr: user-supplied parameters.
             ! info, rinfo: ISAT control and tolerance arrays.
             ! f, g, h: outputs (residual, Jacobian if requested, extra data).
             ! stats: performance statistics array.
          !----------------------------------------------------------------
 
          ! Evaluate the neural network correction for the current state.
          call myfnn( nx, x, fnn )
 
          ! Update the state: add the residual and NN correction to the current state.
          x = x + f + fnn
 
          ! Transform the updated state back to the physical variables for the particle.
          call fromxhat( x, ptcls(1:nx,ii), nx, rusr )
       end do ! End of reaction step loop.
 
       if (modulo(jj,100).eq.0) then
          ! Print particle states again every 100 time steps.
          print *, ptcls
       end if
 
       !-------------------------------
       ! Flow-through (replacement) step for particles:
       do ii = 1, nPtcl
          call random_number(rr)          ! Generate a random number.
          if ( rr < dt/flowThroughTime ) then
             ! With probability dt/flowThroughTime, replace this particle.
             call random_number(rr)       ! Get another random decision.
             if ( rr < mFrac ) then
                ! Initialize as fuel.
                ptcls(1,ii) = Tfuel
                ptcls(2:nx,ii) = 0.0
                ptcls(nH2,ii) = 1.0       ! Fuel inflow.
             else
                ! Initialize as oxidizer.
                ptcls(1,ii) = Tox
                ptcls(2:nx,ii) = 0.0
                ptcls(nO2,ii) = 1.0/4.76
                ptcls(nN2,ii) = 3.76/4.76  ! Oxidizer inflow.
             end if
             counter = counter + 1         ! Update counter.
             ptcls(nx+1,ii) = counter       ! Assign new ID.
          end if
       end do
       !-------------------------------
 
       !-------------------------------
       ! Mixing step: blend states of randomly selected particle pairs.
       alpha(1) = dt/mixTime     ! Determine the mixing parameter.
       do ii = 1, nPtcl
          call random_number(rr)
          i1 = int(ceiling(rr * real(nPtcl)))   ! Randomly select first particle index.
          if (i1 .eq. 0) i1 = 1                   ! Ensure index is valid.
          call random_number(rr)
          i2 = int(ceiling(rr * real(nPtcl)))   ! Randomly select second particle index.
          if (i2 .eq. 0) i2 = 1                   ! Ensure valid index.
 
          x1 = ptcls(1:nx, i1)    ! Extract state of particle i1.
          x2 = ptcls(1:nx, i2)    ! Extract state of particle i2.
 
          call mymix( nx, x1, x2, dt/mixTime, iusr, rusr )
             ! Mix the two states using the given mixing parameter.
          ptcls(1:nx, i1) = x1    ! Write back the updated state for particle i1.
          ptcls(1:nx, i2) = x2    ! Write back the updated state for particle i2.
       end do
       !-------------------------------
    end do ! End of main time-stepping loop.
    !-----------------------------------------------------------------------
 
    ! Get the CPU time at the end of the simulation.
    call cpu_time(finish)
 
    ! A final call to myfgh (perhaps for diagnostics or final evaluation).
    call myfgh( need, nx, x, nf, nh, iusr, rusr, f, g, h )
 
    !-------------------------------
    ! Performance testing: repeated ISAT evaluations.
    do ii = 1, 100
       call cpu_time(start)  ! Record start time.
       call isatab( idd, 0, nx, x, nf, nh, nh, myfgh, iusr, rusr, info, rinfo, f, g, h, stats)
          ! Evaluate ISAT for the state x.
       call cpu_time(finish) ! Record finish time.
       print '("Time = ",e10.3," seconds.")', finish - start
          ! Print the elapsed time for the ISAT call.
       print *, f(1) * rusr(7+1)
          ! Print a scaled version of the first component of the residual f.
       print *, ii       ! Print the current iteration number.
    end do
    !-------------------------------
 
 end program main