module extFGH

   interface

      subroutine myfgh ( need, nx, x, nf, nh, iusr, rusr, f, g, h ) bind ( c )
         use iso_c_binding
         integer (c_int) :: need(*),nx,nf,nh,iusr(*)
         real (c_double) :: x(*),rusr(*),f(*),g(*),h(*)
      end subroutine myfgh

      subroutine mymix ( nx, x1, x2, alpha, iusr, rusr ) bind ( c )
         use iso_c_binding
         integer (c_int) :: nx,iusr(*)
         real (c_double) :: x1(*),x2(*),rusr(*),alpha(*)
      end subroutine mymix

      subroutine toxhat ( ptcl, x, nx, rusr ) bind ( c )
         use iso_c_binding
         integer (c_int) :: nx
         real (c_double) :: x(*),rusr(*),ptcl(*)
      end subroutine toxhat

      subroutine myfnn ( nx, x, fnn ) bind ( c )
         use iso_c_binding
         integer (c_int) :: nx
         real (c_double) :: x(*),fnn(*)
      end subroutine myfnn

      subroutine fromxhat ( x, ptcl, nx, rusr ) bind ( c )
         use iso_c_binding
         integer (c_int) :: nx
         real (c_double) :: x(*),rusr(*),ptcl(*)
      end subroutine fromxhat

   end interface

end module extFGH

program main

   !use isat_rnu
   use extFGH


   implicit none

   integer, pointer :: need(:),need2(:),iusr(:)
   integer :: ii,jj,kk
   integer :: nx=11,nf=11,nh=1, ns=10, info(100), nPtcl = 250, idd, nSteps = 4000, counter, i1, i2
   double precision, pointer    :: x(:), x1(:), x2(:), rusr(:), f(:), fisat(:), fisatsum(:), fnn(:), g(:,:), h(:), ptcls(:,:), alpha(:)
   double precision :: rinfo(70), stats(100), dt, start, finish

   double precision :: flowThroughTime, mixTime

   integer :: nH2 = 2, nO2 = 5, nN2 = 11

   double precision :: mFrac = 0.1, Tfuel = 1050.0, Tox = 1050.0

   double precision :: rr

   allocate(need(3),need2(3),iusr(3),x(nx),x1(nx),x2(nx),rusr(2*nx+5),f(nx),fisat(nx),fisatsum(nx),fnn(nx),g(nx,nx),h(1),ptcls(nx+1,nPtcl),alpha(1))

   dt = 5e-7
   flowThroughTime = 1e-4
   mixTime = 1e-4

   need = (/ 1, 0, 0 /)
   need2 = (/ 1, 1, 0 /)
   rusr(1:2*nx) = 0.0 ! no normalization initially
   rusr(1:nx+1) = (/ 1050.0, 1e-5, 7.5967e-08, 1.4358e-07, 2.1008e-06, 1.5337e-07, 2.1554e-06, 1.4550e-09, 9.1008e-11, 1.0000e-35, 7.8992e-06, 358.9534 /)
   rusr(2*nx+1) = 1e-8 ! atol
   rusr(2*nx+2) = 1e-8 ! rtol
   rusr(2*nx+3) = dt ! dt
   rusr(2*nx+4) = 1e-6 ! dx
   rusr(2*nx+5) = 101325.0 ! p

   info  = 0
   rinfo = 0
   info(2)  = 0	! if_g
   info(12) = 2	! isat op
   info(28) = 0   ! idites
   rinfo(1) = 1e-3 !errtol
   rinfo(2) = 1e-3 !errtol
   rinfo(3) = 1e2
!-- special settings: Table 2
   !info(21,2) = 1   ! ifull
   rinfo(8) = 1000.  !  stomby
!-- special settings: Table 3



   call random_seed(put=(/10,11/))

   call cpu_time(start)

   counter = 0

   do ii = 1,nPtcl ! begin initialization

      call random_number(rr)

      if ( rr < mFrac ) then
         ptcls(1,ii) = Tfuel
         ptcls(2:nx,ii) = 0.0
         ptcls(nH2,ii) = 1.0		! fuel inflow
      else
         ptcls(1,ii) = Tox
         ptcls(2:nx,ii) = 0.0
         ptcls(nO2,ii) = 1.0/4.76
         ptcls(nN2,ii) = 3.76/4.76	! oxidizer inflow
      end if

      counter = counter + 1
      ptcls(nx+1,ii) = counter

   end do ! end initialization



   do jj = 1,nSteps ! begin time step

      if (modulo(jj,100).eq.0) then
         print *, ptcls
      end if


      do ii = 1,nPtcl
         call toxhat( ptcls(1:nx,ii), x, nx, rusr )

         !print *, x

         !call myfgh( need, nx, x, nf, nh, iusr, rusr, f, g, h )

         !print *, f

         call isatab( idd,0,nx,x,nf,nh,nh, myfgh, iusr,rusr, info, rinfo, f,g,h,stats)

         !fisatsum = fisatsum + fisat

         call myfnn( nx, x, fnn )

         x = x + f + fnn

         call fromxhat( x, ptcls(1:nx,ii), nx, rusr )
      end do ! reaction

      if (modulo(jj,100).eq.0) then
         print *, ptcls
      end if

      do ii = 1,nPtcl
         call random_number(rr)
         if ( rr < dt/flowThroughTime ) then

            call random_number(rr)

            if ( rr < mFrac ) then
               ptcls(1,ii) = Tfuel
               ptcls(2:nx,ii) = 0.0
               ptcls(nH2,ii) = 1.0		! fuel inflow
            else
               ptcls(1,ii) = Tox
               ptcls(2:nx,ii) = 0.0
               ptcls(nO2,ii) = 1.0/4.76
               ptcls(nN2,ii) = 3.76/4.76	! oxidizer inflow
            end if

            counter = counter + 1
            ptcls(nx+1,ii) = counter

         end if
      end do

      alpha(1) = dt/mixTime

      do ii = 1,nPtcl
         call random_number(rr)
         i1 = int(ceiling(rr*real(nPtcl)))
         if (i1 .eq. 0) i1 = 1
         call random_number(rr)
         i2 = int(ceiling(rr*real(nPtcl)))
         if (i2 .eq. 0) i2 = 1

         x1 = ptcls(1:nx,i1)
         x2 = ptcls(1:nx,i2)

         call mymix( nx, x1, x2, dt/mixTime, iusr, rusr )

         ptcls(1:nx,i1) = x1
         ptcls(1:nx,i2) = x2

      end do

   end do ! end time step



   call cpu_time(finish)

   !call myfgh( need, nx, x, nf, nh, iusr, rusr, f, g, h )

   !do ii = 1,100

   !  call cpu_time(start)
   !  call isatab( idd,0,nx,x,nf,nh,nh, myfgh, iusr,rusr, info, rinfo, f,g,h,stats)
   !  call cpu_time(finish)
   !  print '("Time = ",e10.3," seconds.")',finish-start

   !  print *, f(1)*rusr(7+1)
   !  print *, ii

   !end do



end program main
