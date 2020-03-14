module edmf

! MKW 2020/03/12: I pulled the edmf plume routine and dependent subroutines out
! of vertical_diffusion.F90 to make it easier to implement elsewhere.

contains

! =============================================================================== !
!  Eddy-diffusivity mass-flux routine                                                                               !
! =============================================================================== !

subroutine edmf( kts,   kte,  dt, zm, pw, iexner, & ! phis,  &
               nup,    u,    v,  th, thl,  thv, qt,                 &
               ust, wthl, wqt, pblh, qc,                             &
         ! outputs updraft properties for diagnostics
               dry_a,   moist_a,                                        &
               dry_w,   moist_w,                                        &
               dry_qt,  moist_qt,                                       &
               dry_thl, moist_thl,                                      &
               dry_u,   moist_u,                                        &
               dry_v,   moist_v,                                        &
                         moist_qc,                                        &
          ! outputs - variables needed for solver
               ae, aw, awthl, awqv, awql, awqi, awu, awv,           &
               thlflx, qtflx, sflx )

! Original author: Marcin Kurowski, JPL
! Modified heavily by Mikael Witte, UCLA/JPL for implementation in CESM2/E3SM
! Last modified 11 Feb. 2020

! MKW TODO: CLUBB operates on thl instead of s or sl; so maybe get rid of sflx and aws
! MKW TODO: make sure all s, sl variables are removed; can also remove phis as input then

! Variables needed for solver
! ae = sum_i (1-a_i)
! aw3 = sum (a_i w_i)
! aws3 = sum (a_i w_i*s_i); s=thl*cp
! aws3,awqv3,awql3,awqi3,awu3,awv3 similar as above except for different variables
!

! MKW: which z grid is the right one? seems like zm since fluxes are calculated on "edges"/interfaces
!      should I therefore input all thermo variables interpolated to momentum grid?

! - mass flux variables are computed on edges:
!  upa,upw,upqt,... kts:kte+1
!  dry_a,moist_a,dry_w,moist_w, ... kts:kte+1


     integer, intent(in) :: kts,kte,nup
     real(r8), dimension(kts:kte), intent(in) :: u,v,th,thl,qt,qc,thv,iexn
     real(r8), dimension(kts-1:kte), intent(in) :: zm,p
     real(r8), dimension(kts:kte), intent(in) :: iexner

     !real(r8), intent(in) :: phis ! surface geopotential
     real(r8), intent(in) :: ust, wthl, wqt,pblh ! MKW TODO: figure out where ust and pblh come - surface_varnce and ???
     real(r8),value :: dt ! only needed for random number generator

! outputs - updraft properties
     real(r8),dimension(kts:kte+1), intent(out) :: dry_a, moist_a, dry_w, moist_w, &
                 dry_qt, moist_qt, dry_thl, moist_thl,          &
                 dry_u,  moist_u,  dry_v,   moist_v,    moist_qc
     real(r8),dimension(kts:kte+1), intent(out) :: ae,aw,awthl,awqv,awql,awqi,awu,awv
     real(r8),dimension(kts:kte+1), intent(out) :: thlflx, qtflx

! INTERNVAL VARIABLES
! sums over all plumes
     real(r8), dimension(kts:kte+1) :: moist_th, dry_th, awqt, awth
! updraft properties
     real(r8), dimension(kts:kte+1,1:nup) ::                &
                 upw, upthl, upqt, upqc, upth, upqv, upql,  &
                 upqi, upa, upu, upv, upthv, ups
! entrainment variables
     real(r8), dimension(kts+1:kte,1:nup) :: ent,entf
     integer,  dimension(kts+1:kte,1:nup) :: enti
! internal variables
     integer :: k,i,ih
     real(r8) :: wthv, wstar, qstar, thstar, sigmaw, sigmaqt, sigmath, z0, &
                 wmin, wmax, wlv, wtv, wp
     real(r8) :: b, qtn, thln, thvn, thn, qcn, qln, qin, un, vn, wn2, &
                 entexp, entexpu, entw

! internal surface cont

     real(r8) :: iexh

! w parameters
     real(r8),parameter :: &
       wa = 1., &
       wb = 1.5

! entrainment parameters
     real(r8),parameter :: &
!      L0   = 150.,&
!       ENT0 = .5
!!       L0   = 150.,&
!!       ENT0 = .8
!        L0   = 100.,&
!        ENT0 = .42
      L0   = 50.,&
      ENT0 = .22
!!      L0   = 50.,&
!!       ENT0 = .18
!       L0   = 25.,&
!       ENT0 = .11

!! parameters defining initial conditions for updrafts
     real(r8),parameter :: &
!       pwmin = 1.55,&
     pwmin = 1.5,&
     pwmax = 3.

! min values to avoid singularities
     real(r8),parameter :: &
     wstarmin = 1.e-3, &
     pblhmin  = 100.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! INITIALIZE OUTPUT VARIABLES
! set updraft properties to zero
     dry_a     = 0.
     moist_a   = 0.
     dry_w     = 0.
     moist_w   = 0.
     dry_qt    = 0.
     moist_qt  = 0.
     dry_thl   = 0.
     moist_thl = 0.
     dry_u     = 0.
     moist_u   = 0.
     dry_v     = 0.
     moist_v   = 0.
     moist_qc  = 0.
! outputs - variables needed for solver
     aw        = 0.
     aws       = 0.
     awqv      = 0.
     awql      = 0.
     awqi      = 0.
     awu       = 0.
     awv       = 0.

! this is the environmental area - by default 1.
     ae = 1.

! START MAIN COMPUTATION
     upw   = 0.
     upthl = 0.
     upthv = 0.
     upqt  = 0.
     upa   = 0.
     upu   = 0.
     upv   = 0.
     upqc  = 0.
     ent   = 0.
     upth  = 0.
     upql  = 0.
     upqi  = 0.
     upqv  = 0.

     pblh = max(pblh,pblhmin)
     wthv = wthl+zvir*thv(kts)*wqt

     print *,'zm=',zm
     print *,'pw=',p
     print *,'wthl',wthl
     print *,'wqt',wqt
     print *,'wthv',wthv

     print*,'iexn 1D=',iexn
     ! if surface buoyancy is positive then do mass-flux, otherwise not
     if (.false.) then !wthv > 0.0 ) then

       ! compute entrainment coefficient
       ! get dz/L0
       do i=1,nup
         do k=kts+1,kte
           entf(k,i) = (zm(k)-zm(k-1)) / L0 ! MKW TODO: replace all instances of zm(k)-zm(k-1) with dzm variable?
!               print*,'k,zmu,zmd=',k,zm(k),zm(k-1)
         enddo
       enddo

       ! get Poisson P(dz/L0)
       call Poisson( 1, nup, kts+1, kte, entf, enti)

       ! entrainent: Ent=Ent0/dz*P(dz/L0)
       do i=1,nup
         do k=kts+1,kte
           ent(k,i) = real( enti(k,i))*ent0/(zm(k)-zm(k-1) )
         enddo
       enddo

       ! surface conditions
!          wstar  = max( wstarmin, (gravit/thv(1)*wthv*pblh)**(1./3.) )
       wstar  = max( wstarmin, (gravit/thv(kts)*wthv*pblh)**(1./3.) )
       qstar  = wqt / wstar
       thstar = wthl / wstar

       print*,'wstar=',wstar
       print*,'qstar=',qstar
       print*,'thstar=',thstar

       sigmaw  = 0.572 * wstar     / 1.
       sigmaqt = 2.89 * abs(qstar) / 1.
       sigmath = 2.89 * abs(thstar)/ 1.

       wmin = sigmaw * pwmin
       wmax = sigmaw * pwmax

       do i=1,nup

         wlv = wmin + (wmax-wmin) / (real(nup)) * (real(i)-1.)
         wtv = wmin + (wmax-wmin) / (real(nup)) * real(i)

         upw(kts,i) = 0.5 * (wlv+wtv)
         upa(kts,i) = 0.5 * erf( wtv/(sqrt(2.)*sigmaw) ) &
                    - 0.5 * erf( wlv/(sqrt(2.)*sigmaw) )

         upu( kts, i) = u(kts)
         upv( kts, i) = v(kts)

         upqc(kts,i)  = 0.
         upqt(kts,i)  = qt(kts)  + 0.32 * upw(kts,i) * sigmaqt/sigmaw
         upthv(kts,i) = thv(kts) + 0.58 * upw(kts,i) * sigmath/sigmaw
         upthl(kts,i) = upthv(kts,i) / (1.+zvir*upqt(kts,i))
         upth(kts,i)  = upthl(kts,i)
         upqv(kts,i)  = upqt(kts,i)

       enddo


       ! integrate updrafts
       do i=1,nup
         do k=kts+1,kte+1

           entexp  = exp(-ent(k,i)*(zm(k)-zm(k-1)))
           entexpu = exp(-ent(k,i)*(zm(k)-zm(k-1))/3.)  !orig
!              entexpu = exp(-ent(k,i)*(zm(k)-zm(k-1))/1.)

           qtn  = qt(k-1) *(1.-entexp ) + upqt (k-1,i)*entexp
           thln = thl(k-1)*(1.-entexp ) + upthl(k-1,i)*entexp
           un   = u(k-1)  *(1.-entexpu) + upu  (k-1,i)*entexpu
           vn   = v(k-1)  *(1.-entexpu) + upv  (k-1,i)*entexpu

           iexh = (1.e5/p(k))**(rair/cpair) ! MKW NOTE: why not just use CLUBB exner??
           call condensation_edmf(qtn, thln, p(k), iexh, &
                                 thvn, qcn, thn, qln, qin)

!Condensation within updrafts, input/output at full levels:
!               call condensation_edmf(qtn, thln, p(k), (iexn(k-1)+iexn(k))/2., &
!                                     thvn, qcn, thn, qln, qin)

           b=gravit*(0.5*(thvn+upthv(k-1,i))/thv(k-1)-1.)
           !b=mapl_grav*(thvn/thv(k-1)-1.)
           ! Wn^2
           ! to avoid singularities w equation has to be computed diferently if wp==0

           wp = wb*ent(k,i)
           if (wp==0.) then
             wn2 = upw(k-1,i)**2+2.*wa*b*(zm(k)-zm(k-1))
           else
             entw = exp(-2.*wp*(zm(k)-zm(k-1)))
             wn2 = entw*upw(k-1,i)**2+wa*b/(wb*ent(k,i))*(1.-entw)
           end if

           if (wn2>0.) then
             upw(k,i)   = sqrt(wn2)
             upthv(k,i) = thvn
             upthl(k,i) = thln
             upqt(k,i)  = qtn
             upqc(k,i)  = qcn
             upu(k,i)   = un
             upv(k,i)   = vn
             upa(k,i)   = upa(k-1,i)
             upth(k,i)  = thn
             upql(k,i)  = qln
             upqi(k,i)  = qin
             upqv(k,i)  = qtn - qcn
           else
             exit
           end if
         enddo
       enddo

       ! writing updraft properties for output
       ! all variables, except Areas are now multipled by the area
       ! to confirm with WRF grid setup we do not save the first and the last row
       do k=kts,kte+1

         ! first sum over all i-updrafts
         do i=1,nup
           if (upqc(k,i)>0.) then
             moist_a(k)   = moist_a(k)   + upa(k,i)
             moist_w(k)   = moist_w(k)   + upa(k,i)*upw(k,i)
             moist_qt(k)  = moist_qt(k)  + upa(k,i)*upqt(k,i)
             moist_thl(k) = moist_thl(k) + upa(k,i)*upthl(k,i)
             moist_th(k)  = moist_th(k)  + upa(k,i)*upth(k,i)
             moist_u(k)   = moist_u(k)   + upa(k,i)*upu(k,i)
             moist_v(k)   = moist_v(k)   + upa(k,i)*upv(k,i)
             moist_qc(k)  = moist_qc(k)  + upa(k,i)*upqc(k,i)
           else
             dry_a(k)     = dry_a(k)     + upa(k,i)
             dry_w(k)     = dry_w(k)     + upa(k,i)*upw(k,i)
             dry_qt(k)    = dry_qt(k)    + upa(k,i)*upqt(k,i)
             dry_thl(k)   = dry_thl(k)   + upa(k,i)*upthl(k,i)
             dry_th(k)    = dry_th(k)    + upa(k,i)*upth(k,i)
             dry_u(k)     = dry_u(k)     + upa(k,i)*upu(k,i)
             dry_v(k)     = dry_v(k)     + upa(k,i)*upv(k,i)
           endif
         enddo

         if ( dry_a(k) > 0. ) then
           dry_w(k)   = dry_w(k)   / dry_a(k)
           dry_qt(k)  = dry_qt(k)  / dry_a(k)
           dry_thl(k) = dry_thl(k) / dry_a(k)
           dry_th(k)  = dry_th(k)  / dry_a(k)
           dry_u(k)   = dry_u(k)   / dry_a(k)
           dry_v(k)   = dry_v(k)   / dry_a(k)
         else
           dry_w(k)   = 0.
           dry_qt(k)  = 0.
           dry_thl(k) = 0.
           dry_th(k)  = 0.
           dry_u(k)   = 0.
           dry_v(k)   = 0.
         endif

         if ( moist_a(k) > 0. ) then
           moist_w(k)   = moist_w(k)   / moist_a(k)
           moist_qt(k)  = moist_qt(k)  / moist_a(k)
           moist_thl(k) = moist_thl(k) / moist_a(k)
           moist_th(k)  = moist_th(k)  / moist_a(k)
           moist_u(k)   = moist_u(k)   / moist_a(k)
           moist_v(k)   = moist_v(k)   / moist_a(k)
           moist_qc(k)  = moist_qc(k)  / moist_a(k)
         else
           moist_w(k)   = 0.
           moist_qt(k)  = 0.
           moist_thl(k) = 0.
           moist_th(k)  = 0.
           moist_u(k)   = 0.
           moist_v(k)   = 0.
           moist_qc(k)  = 0.
         endif

       enddo

       ! print*,'PHIS(i)=',phis
! Dry static energy
! phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
!                          + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)


       do k=kts,kte+1
!             iexh = (p(1)/p(k))**(rair/cpair) !as in dp_coupling.F90
!if you modify dp_coupling.F90 to use 1.e5 as reference pressure:
         iexh = (1.e5/p(k))**(rair/cpair)
         print*,'IEXH=',iexh

!CHECK SURF GEOPOTENTIAL
!         do i=1,nup
!           ups(k,i)=cpair*upth(k,i)/iexh  + gravit*zm(k) + phis(ih) !surface geopotential
!         enddo
       enddo

       do k=kts,kte+1
         do i=1,nup
           aw  (k) = aw  (k) + upa(k,i)*upw(k,i)
           awu (k) = awu (k) + upa(k,i)*upw(k,i)*upu(k,i)
           awv (k) = awv (k) + upa(k,i)*upw(k,i)*upv(k,i)
           !aws (k) = aws (k) + upa(k,i)*upw(k,i)*upth(k,i)*cpair
           !aws (k) = aws (k) + upa(k,i)*upw(k,i)*ups(k,i)
           awthl(k)= awthl(k)+ upa(k,i)*upw(k,i)*upthl(k,i) !*cpair/iexh
           awth(k) = awth(k) + upa(k,i)*upw(k,i)*upth(k,i) !*cpair/iexh
           awqt(k) = awqt(k) + upa(k,i)*upw(k,i)*upqt(k,i)
           awqv(k) = awqv(k) + upa(k,i)*upw(k,i)*upqv(k,i)
           awql(k) = awql(k) + upa(k,i)*upw(k,i)*upql(k,i)
           awqi(k) = awqi(k) + upa(k,i)*upw(k,i)*upqi(k,i)
         enddo
       enddo

       do k=kts+1,kte+1
         iexh = (1.e5/p(k))**(rair/cpair)
         thlflx(k)= (awthl(k) - aw(k)*0.5*(thl(k-1)+thl(k)) ) ! MKW NOTE: used to be slflx, but CLUBB works on thl
         !sflx( k)= (awth(k) - aw(k)*0.5*(th(k-1)+th(k)) )*cpair/iexh
         qtflx(k)= awqt(k)  - aw(k)*0.5*(qt(k-1)+qt(k))
       enddo
       !iexh = (1.e5/p(kts))**(rair/cpair)
       thlflx(kts) = 0.
       !sflx(kts)  = 0.
       qtflx(kts) = 0.

     end if  ! ( wthv > 0.0 )


end subroutine edmf


subroutine condensation_edmf( qt, thl, p, iex, thv, qc, th, ql, qi)
!
! zero or one condensation for edmf: calculates thv and qc
!
     use wv_saturation,      only : qsat

     real(r8),intent(in) :: qt,thl,p,iex
     real(r8),intent(out):: thv,qc,th,ql,qi

     !local variables
     integer :: niter,i
     real(r8) :: diff,t,qs,qcold,es,wf

     ! max number of iterations
     niter=50
     ! minimum difference
     diff=2.e-5

     qc=0.
     t=thl/iex

!by definition:
! T   = Th*Exner, Exner=(p/p0)^(R/cp)   (1)
! Thl = Th - L/cp*ql/Exner              (2)
!so:
! Th  = Thl + L/cp*ql/Exner             (3)
! T   = Th*Exner=(Thl+L/cp*ql/Exner)*Exner    (4)
!     = Thl*Exner + L/cp*ql
     do i=1,niter
       wf = get_watf(t)
       t = thl/iex+get_alhl(wf)/cpair*qc   !as in (4)

       ! qsat, p is in pascal (check!)
       call qsat(t,p,es,qs)
       qcold = qc
       qc = max(0.5*qc+0.5*(qt-qs),0.)
       if (abs(qc-qcold)<diff) exit
     enddo

     wf = get_watf(t)
     t = thl/iex+get_alhl(wf)/cpair*qc


     call qsat(t,p,es,qs)
     qc = max(qt-qs,0.)
     thv = (thl+get_alhl(wf)/cpair*iex*qc)*(1.+zvir*(qt-qc)-qc)
     th = t*iex
     qi = qc*(1.-wf)
     ql = qc*wf

     contains

     function get_watf(t)
       real(r8) :: t,get_watf,tc
       real(r8), parameter :: &
       tmax=-10., &
       tmin=-40.

       tc=t-273.16

       if (tc>tmax) then
         get_watf=1.
       else if (tc<tmin) then
         get_watf=0.
       else
         get_watf=(tc-tmin)/(tmax-tmin);
       end if

     end function get_watf


     function get_alhl(wf)
     !latent heat of the mixture based on water fraction
       use physconst,        only : latvap , latice

       real(r8) :: get_alhl,wf

       get_alhl = wf*latvap+(1.-wf)*(latvap+latice)

     end function get_alhl

end subroutine condensation_edmf

subroutine Poisson(istart,iend,jstart,jend,mu,poi)
       integer, intent(in) :: istart,iend,jstart,jend
       real(r8), dimension(istart:iend,jstart:jend),intent(in) :: mu
       integer, dimension(istart:iend,jstart:jend), intent(out) :: poi

       integer :: i,j

       ! do this only once
       ! call init_random_seed

       do i=istart,iend
         do j=jstart,jend
           call   random_poisson(mu(i,j),.true.,poi(i,j))
         enddo
       enddo

end subroutine Poisson


subroutine init_random_seed()
       ! use iso_fortran_env, only: int64
       ! use ifport, only: getpid
       use,intrinsic :: iso_fortran_env
       implicit none
       integer, allocatable :: seed(:)
       integer :: i, n, un, istat, dt(8), pid
       integer(int64) :: t

       call random_seed(size = n)
       allocate(seed(n))
       ! First try if the OS provides a random number generator
       open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
       if (istat == 0) then
         read(un) seed
         close(un)
       else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
         call system_clock(t)
         if (t == 0) then
           call date_and_time(values=dt)
           t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
         end if
         !pid = getpid()
         ! for distributed memory jobs we need to fix this
         pid = 1
         t = ieor(t, int(pid, kind(t)))
         do i=1,n
           seed(i) = lcg(t)
         end do
       end if
       call random_seed(put=seed)

       contains
       ! This simple PRNG might not be good enough for real work, but is
       ! sufficient for seeding a better PRNG.
       function lcg(s)
         integer :: lcg
         integer(int64) :: s
         if (s == 0) then
           s = 104729
         else
           s = mod(s, 4294967296_int64)
         end if
         s = mod(s * 279470273_int64, 4294967291_int64)
         lcg = int(mod(s, int(huge(0), int64)), kind(0))
       end function lcg

end subroutine init_random_seed


subroutine random_Poisson(mu,first,ival)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                           RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                    Compiled and Written by:
!
!                         Barry W. Brown
!                          James Lovato
!
!             Department of Biomathematics, Box 237
!             The University of Texas, M.D. Anderson Cancer Center
!             1515 Holcombe Boulevard
!             Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate POIsson random deviate
!                            Function
! Generates a single random deviate from a Poisson distribution with mean mu.
!                            Arguments
!     mu --> The mean of the Poisson distribution from which
!            a random deviate is to be generated.
!                              REAL mu
!                              Method
!     For details see:
!               Ahrens, J.H. and Dieter, U.
!               Computer Generation of Poisson Deviates
!               From Modified Normal Distributions.
!               ACM Trans. Math. Software, 8, 2
!               (June 1982),163-179
!     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
!     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
!     SEPARATION OF CASES A AND B

!     .. Scalar Arguments ..
	REAL(r8), INTENT(IN)    :: mu
	LOGICAL, INTENT(IN) :: first
INTEGER             :: ival
!     ..
!     .. Local Scalars ..
	REAL(r8)          :: b1, b2, c, c0, c1, c2, c3, del, difmuk, e, fk, fx, fy, g,  &
                    omega, px, py, t, u, v, x, xx
	REAL(r8), SAVE    :: s, d, p, q, p0
        INTEGER       :: j, k, kflag
	LOGICAL, SAVE :: full_init
        INTEGER, SAVE :: l, m
!     ..
!     .. Local Arrays ..
	REAL(r8), SAVE    :: pp(35)
!     ..
!     .. Data statements ..
	REAL(r8), PARAMETER :: a0 = -.5, a1 = .3333333, a2 = -.2500068, a3 = .2000118,  &
                a4 = -.1661269, a5 = .1421878, a6 = -.1384794,   &
                 a7 = .1250060

	REAL(r8), PARAMETER :: fact(10) = (/ 1., 1., 2., 6., 24., 120., 720., 5040.,  &
            40320., 362880. /)

!     ..
!     .. Executable Statements ..
   IF (mu > 10.0) THEN
!     C A S E  A. (RECALCULATION OF S, D, L IF MU HAS CHANGED)

  IF (first) THEN
s = SQRT(mu)
d = 6.0*mu*mu

!             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
!             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
!             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .

l = mu - 1.1484
full_init = .false.
  END IF


!     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE

	  g = mu + s*random_normal()
	  IF (g > 0.0) THEN
		ival = g

	!     STEP I. IMMEDIATE ACCEPTANCE IF ival IS LARGE ENOUGH

		IF (ival>=l) RETURN

	!     STEP S. SQUEEZE ACCEPTANCE - SAMPLE U

		fk = ival
		difmuk = mu - fk
		CALL RANDOM_NUMBER(u)
		IF (d*u >= difmuk*difmuk*difmuk) RETURN
	  END IF

	!     STEP P. PREPARATIONS FOR STEPS Q AND H.
	!             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
	!             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
	!             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
	!             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
	!             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.

	  IF (.NOT. full_init) THEN
		omega = .3989423/s
		b1 = .4166667E-1/mu
		b2 = .3*b1*b1
		c3 = .1428571*b1*b2
		c2 = b2 - 15.*c3
		c1 = b1 - 6.*b2 + 45.*c3
		c0 = 1. - b1 + 3.*b2 - 15.*c3
		c = .1069/mu
		full_init = .true.
	  END IF

	  IF (g < 0.0) GO TO 50

	!             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)

	  kflag = 0
	  GO TO 70

	!     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)

	  40 IF (fy-u*fy <= py*EXP(px-fx)) RETURN

	!     STEP E. EXPONENTIAL SAMPLE - random_exponential() FOR STANDARD EXPONENTIAL
	!             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
	!             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)

	  50 e = random_exponential()
	  CALL RANDOM_NUMBER(u)
	  u = u + u - 1
	  t = 1.8 + SIGN(e, u)
	  IF (t <= (-.6744)) GO TO 50
	  ival = mu + s*t
	  fk = ival
	  difmuk = mu - fk

	!             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)

	  kflag = 1
	  GO TO 70

	!     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)

	  60 IF (c*ABS(u) > py*EXP(px+e) - fy*EXP(fx+e)) GO TO 50
	  RETURN

	!     STEP F. 'SUBROUTINE' F. CALCULATION OF PX, PY, FX, FY.
	!             CASE ival < 10 USES FACTORIALS FROM TABLE FACT

	  70 IF (ival>=10) GO TO 80
	  px = -mu
	  py = mu**ival/fact(ival+1)
	  GO TO 110

	!             CASE ival >= 10 USES POLYNOMIAL APPROXIMATION
	!             A0-A7 FOR ACCURACY WHEN ADVISABLE
	!             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)

	  80 del = .8333333E-1/fk
	  del = del - 4.8*del*del*del
	  v = difmuk/fk
	  IF (ABS(v)>0.25) THEN
		px = fk*LOG(1. + v) - difmuk - del
	  ELSE
		px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
	  END IF
	  py = .3989423/SQRT(fk)
	  110 x = (0.5 - difmuk)/s
	  xx = x*x
	  fx = -0.5*xx
	  fy = omega* (((c3*xx + c2)*xx + c1)*xx + c0)
	  IF (kflag <= 0) GO TO 40
	  GO TO 60

	!---------------------------------------------------------------------------
	!     C A S E  B.    mu < 10
	!     START NEW TABLE AND CALCULATE P0 IF NECESSARY

	ELSE
	  IF (first) THEN
		m = MAX(1, INT(mu))
		l = 0
		p = EXP(-mu)
		q = p
		p0 = p
	  END IF

	!     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD

	  DO
		CALL RANDOM_NUMBER(u)
		ival = 0
		IF (u <= p0) RETURN

	!     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
	!             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
	!             (0.458=PP(9) FOR MU=10)

		IF (l == 0) GO TO 150
		j = 1
		IF (u > 0.458) j = MIN(l, m)
		DO k = j, l
		  IF (u <= pp(k)) GO TO 180
		END DO
		IF (l == 35) CYCLE

	!     STEP C. CREATION OF NEW POISSON PROBABILITIES P
	!             AND THEIR CUMULATIVES Q=PP(K)

		150 l = l + 1
		DO k = l, 35
		  p = p*mu / k
		  q = q + p
		  pp(k) = q
		  IF (u <= q) GO TO 170
		END DO
		l = 35
	  END DO

	  170 l = k
	  180 ival = k
	  RETURN
	END IF

	RETURN
	END subroutine random_Poisson



	FUNCTION random_normal() RESULT(fn_val)

	! Adapted from the following Fortran 77 code
	!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
	!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
	!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

	!  The function random_normal() returns a normally distributed pseudo-random
	!  number with zero mean and unit variance.

	!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
	!  and J.F. Monahan augmented with quadratic bounding curves.

	REAL(r8) :: fn_val

	!     Local variables
	REAL(r8)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
				r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

	!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

	DO
	  CALL RANDOM_NUMBER(u)
	  CALL RANDOM_NUMBER(v)
	  v = 1.7156 * (v - 0.5 )

	!     Evaluate the quadratic form
	  x = u - s
	  y = ABS(v) - t
	  q = x**2 + y*(a*y - b*x)

	!     Accept P if inside inner ellipse
	  IF (q < r1) EXIT
	!     Reject P if outside outer ellipse
	  IF (q > r2) CYCLE
	!     Reject P if outside acceptance region
	  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
	END DO

	!     Return ratio of P's coordinates as the normal deviate
	fn_val = v/u
	RETURN

	END FUNCTION random_normal





	FUNCTION random_exponential() RESULT(fn_val)

	! Adapted from Fortran 77 code from the book:
	!     Dagpunar, J. 'Principles of random variate generation'
	!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

	! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
	! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
	! TO EXP(-random_exponential), USING INVERSION.

	REAL(r8)  :: fn_val

	!     Local variable
	REAL(r8)  :: r

	DO
	  CALL RANDOM_NUMBER(r)
	  IF (r > 0.) EXIT
	END DO

	fn_val = -LOG(r)
	RETURN

	END FUNCTION random_exponential

end module edmf
