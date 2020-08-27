module edmf

  use physics_utils, only: rtype, rtype8, itype!, btype ! MKW: btype not currently used

  use physconst,     only: rgas => rair, cp => cpair, ggr => gravit, &
                           lcond => latvap, lice => latice, eps => zvir
  ! use shoc,          only: linear_interp

  implicit none

  public :: integrate_mf, init_random_seed, calc_mf_vertflux, compute_tmpi3

  private

  !=========================================================
  ! Physical constants used for mass flux plumes
  !=========================================================

  !! MKW: is it easier to inherit these from SHOC for now?
  !use physconst,     only: rair => rgas, cpair => cp, gravit => ggr, &
  !                         latvap => lcond, latice => lice, zvir => eps
  !! MKW: this is how these variables are initialized in module shoc - values are set in shoc_init which is called from shoc_init_e3sm
  ! ! These are set in initialization and should be set to
  ! !  to the values used in whatever host model SHOC is
  ! !  implemented in
  ! real(rtype) :: ggr   ! gravity [m/s^2]
  ! real(rtype) :: rgas  ! dry air gas constant [J/kg.K]
  ! real(rtype) :: cp    ! specific heat of dry air [J/kg.K]
  ! real(rtype) :: lcond ! latent heat of vaporization [J/kg]
  ! real(rtype) :: lice  ! latent heat of fusion [J/kg]
  ! real(rtype) :: eps   ! rh2o/rair - 1 [-]


contains

  ! TODO: Consider an "edmf_init" routine to initialize physical constants

  ! =============================================================================== !
  !  Eddy-diffusivity mass-flux routine                                                                               !
  ! =============================================================================== !

  subroutine integrate_mf(shcol, nz, nzi, dt,                      & ! input
                 zt_in, zi_in, dz_zt_in, p_in,                     & ! input - MKW 20200804 removed iex and dz_zi_in
                 nup,    u_in,   v_in,   thl_in,   thv_in, qt_in,  & ! input
                 ust,    wthl,   wqt,    pblh,     qc_in,          & ! input
                 dry_a_out,   moist_a_out,                         & ! output: updraft properties for diagnostics
                 dry_w_out,   moist_w_out,                         & ! output: updraft properties for diagnostics
                 dry_qt_out,  moist_qt_out,                        & ! output: updraft properties for diagnostics
                 dry_thl_out, moist_thl_out,                       & ! output: updraft properties for diagnostics
                 dry_u_out,   moist_u_out,                         & ! output: updraft properties for diagnostics
                 dry_v_out,   moist_v_out,                         & ! output: updraft properties for diagnostics
                              moist_qc_out,                        & ! output: updraft properties for diagnostics
                 ae_out, aw_out,                                   & ! output: variables needed for  diffusion solver
                 awthl_out, awqt_out,                              & ! output: variables needed for  diffusion solver
                 awql_out, awqi_out,                               & ! output: variables needed for  diffusion solver
                 awu_out, awv_out)                                   ! output: variables needed for  diffusion solver
                 !thlflx_out, qtflx_out )                             ! output: MF turbulent flux diagnostics

  ! Original author: Marcin Kurowski, JPL
  ! Modified heavily by Mikael Witte and Maria Chinita Candeias, UCLA/JPL for implementation in E3SM
  ! Last modified 4 Aug. 2020

  !
  ! Variables needed for solver:
  ! ae = sum_i (1-a_i)
  ! aw = sum (a_i w_i)
  ! awthl = sum(a_i w_i*thl_i)
  ! awqt  = sum(a_i w_i*qt_i)
  ! awql,awqi,awu,awv similar to above except for different variables - not currently coupled to SHOC diffusion solver
  !

  ! - mass flux variables are computed on edges (i.e. momentum grid):
  !  upa,upw,upqt,... 1:nzi
  !  dry_a,moist_a,dry_w,moist_w, ... 1:nzi
       ! Variable(s)
       integer, intent(in) :: shcol,nz,nzi,nup
       real(rtype), dimension(shcol,nz),  intent(in) :: zt_in,dz_zt_in,p_in !,iex_in
       real(rtype), dimension(shcol,nzi), intent(in) :: zi_in
       real(rtype), dimension(shcol,nz),  intent(in) :: u_in,v_in,thl_in,qt_in,qc_in,thv_in  ! all on thermodynamic/midpoint levels

       real(rtype), dimension(shcol), intent(in) :: ust, wthl, wqt
       real(rtype), dimension(shcol), intent(in) :: pblh
       real(rtype), value :: dt ! only needed for random number generator

  ! outputs - updraft properties
       real(rtype),dimension(shcol,nzi), intent(out) :: &
              dry_a_out, moist_a_out, dry_w_out, moist_w_out,                &
              dry_qt_out, moist_qt_out, dry_thl_out, moist_thl_out,          &
              dry_u_out,  moist_u_out,  dry_v_out,   moist_v_out,    moist_qc_out
  ! outputs - variables needed for diffusion solver
       real(rtype),dimension(shcol,nzi), intent(out) :: &
              ae_out,aw_out,awthl_out,awqt_out,awql_out,awqi_out,awu_out,awv_out
  ! outputs - flux diagnostics
       !real(rtype),dimension(shcol,nzi), intent(out) :: thlflx_out, qtflx_out

  ! INTERNAL VARIABLES
  ! flipped variables (i.e. index 1 is at surface)
       real(rtype), dimension(shcol,nz)  :: zt, dz_zt, iexner, p
       real(rtype), dimension(shcol,nzi) :: zi
       real(rtype), dimension(shcol,nz)  :: u,v,thl,qt,qc,thv
  ! flipped updraft properties (i.e. index 1 is at surface)
       real(rtype), dimension(shcol,nzi) :: dry_a, moist_a, dry_w, moist_w, &
                                   dry_qt, moist_qt, dry_thl, moist_thl, &
                                   dry_u, moist_u, dry_v, moist_v, moist_qc
       real(rtype), dimension(shcol,nzi) :: ae, aw, awthl, awqt, awql, awqi, awu, awv
       !real(rtype), dimension(shcol,nzi) :: thlflx, qtflx

  ! sums over all plumes
       real(rtype), dimension(shcol,nz) :: moist_th, dry_th, awqv, awth

  ! updraft properties
       real(rtype), dimension(nzi,nup) ::                &
                   upw, upthl, upqt, upqc, upth, upqv, upql,  &
                   upqi, upa, upu, upv, upthv, ups
  ! entrainment variables
       real(rtype), dimension(nz,nup) :: ent,entf
       integer,  dimension(nz,nup) :: enti
  ! internal variables
       integer :: k,j,i,ih
       real(rtype) :: wthv, wstar, qstar, thstar, sigmaw, sigmaqt, sigmath, z0, &
                   wmin, wmax, wlv, wtv, wp
       real(rtype) :: pbj, b, qtn, thln, thvn, thn, qcn, qln, qin, un, vn, wn2, &
                   entexp, entexpu, entw

  ! internal surface cont

       real(rtype) :: iexh
       real(rtype) :: dzt(nz)!, dzi(nzi)
       real(rtype) :: thl_zi(nzi),qt_zi(nzi)

  ! w parameters
       real(rtype),parameter :: &
         wa = 1._rtype, &
         wb = 1.5_rtype

  ! entrainment parameters
       real(rtype),parameter :: &
  !      L0   = 150._rtype,&
  !       ent0 = .5_rtype
  !!       L0   = 150._rtype,&
  !!       ent0 = .8_rtype
  !        L0   = 100._rtype,&
  !        ent0 = .42_rtype
           L0   = 50._rtype,&
           ent0 = .22_rtype
  !!      L0   = 50._rtype,&
  !!       ent0 = .18_rtype
  !       L0   = 25._rtype,&
  !       ent0 = .11_rtype

  !! parameters defining initial conditions for updrafts
       real(rtype),parameter :: &
  !       pwmin = 1.55_rtype,&
          pwmin = 1.5_rtype,&
          pwmax = 3._rtype

  ! min values to avoid singularities
       real(rtype),parameter :: &
          wstarmin = 1.e-3, &
          pblhmin  = 100.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  ! Flip vertical coordinates and all input variables
     do k=1,nzi
       ! thermodynamic grid variables
       if (k<nzi) then
         zt(:,k) = zt_in(:,nz-k+1)
         dz_zt(:,k) = dz_zt_in(:,nz-k+1)
         ! iexner(:,k) = iex_in(:,nz-k+1)
         p(:,k) = p_in(:,nz-k+1)

         u(:,k) = u_in(:,nz-k+1)
         v(:,k) = v_in(:,nz-k+1)
         thl(:,k) = thl_in(:,nz-k+1)
         qt(:,k) = qt_in(:,nz-k+1)
         qc(:,k) = qc_in(:,nz-k+1)
         thv(:,k) = thv_in(:,nz-k+1)
       endif

       ! momentum altitude grid
       zi(:,k) = zi_in(:,nzi-k+1)
     enddo


  ! INITIALIZE OUTPUT VARIABLES
  ! set updraft properties to zero
     dry_a     = 0._rtype
     moist_a   = 0._rtype
     dry_w     = 0._rtype
     moist_w   = 0._rtype
     dry_qt    = 0._rtype
     moist_qt  = 0._rtype
     dry_thl   = 0._rtype
     moist_thl = 0._rtype
     dry_u     = 0._rtype
     moist_u   = 0._rtype
     dry_v     = 0._rtype
     moist_v   = 0._rtype
     moist_qc  = 0._rtype
  ! outputs - variables needed for solver
     aw        = 0._rtype
     ! aws       = 0._rtype
     awthl     = 0._rtype
     awqt      = 0._rtype
     awqv      = 0._rtype
     awql      = 0._rtype
     awqi      = 0._rtype
     awu       = 0._rtype
     awv       = 0._rtype
  ! outputs - diagnostics
     !thlflx    = 0._rtype
     !qtflx     = 0._rtype

  ! this is the environmental area - by default 1.
     ae = 1._rtype

  ! START MAIN COMPUTATION
  ! NOTE: SHOC does not invert the vertical coordinate, which by default is ordered from lowest to highest pressure
  !     (i.e. top of atmosphere to bottom) so surface-based do loops are performed in reverse (i.e. from nz to 1)
     do j=1,shcol
       ! zero out plume properties
       upw   = 0._rtype
       upthl = 0._rtype
       upthv = 0._rtype
       upqt  = 0._rtype
       upa   = 0._rtype
       upu   = 0._rtype
       upv   = 0._rtype
       upqc  = 0._rtype
       ent   = 0._rtype
       upth  = 0._rtype
       upql  = 0._rtype
       upqi  = 0._rtype
       upqv  = 0._rtype

       pbj = max(pblh(j),pblhmin)
       wthv = wthl(j)+eps*thv(j,1)*wqt(j)

       ! if surface buoyancy is positive then do mass-flux, otherwise not
       if (wthv>0.0) then
         dzt = dz_zt(j,:)

         ! interpolate thl and qt to interface grid
         call linear_interp(zt,zi,thl(j,:),thl_zi,nz,nzi,shcol,0._rtype)
         call linear_interp(zt,zi,qt(j,:),qt_zi,nz,nzi,shcol,0._rtype)

         ! compute entrainment coefficient
         ! get dz/L0
         do i=1,nup
           do k=2,nz
             entf(k,i) = dzt(k) / L0
           enddo
         enddo

         ! get Poisson P(dz/L0)
         call Poisson( 1, nz, 1, nup, entf, enti)

         ! entrainment: Ent=Ent0/dz*P(dz/L0)
         do i=1,nup
           do k=2,nz
             ent(k,i) = real( enti(k,i))*ent0/dzt(k)
           enddo
         enddo

         ! surface conditions
         wstar  = max( wstarmin, (ggr/thv(j,1)*wthv*pbj)**(1._rtype/3._rtype) )
         qstar  = wqt(j) / wstar
         thstar = wthl(j) / wstar

  !       print*,'wstar=',wstar
  !       print*,'qstar=',qstar
  !       print*,'thstar=',thstar

         sigmaw  = 0.572_rtype * wstar     / 1._rtype
         sigmaqt = 2.89_rtype * abs(qstar) / 1._rtype
         sigmath = 2.89_rtype * abs(thstar)/ 1._rtype

         wmin = sigmaw * pwmin
         wmax = sigmaw * pwmax

         do i=1,nup

           wlv = wmin + (wmax-wmin) / (real(nup)) * (real(i)-1._rtype)
           wtv = wmin + (wmax-wmin) / (real(nup)) * real(i)

           upw(1,i) = 0.5_rtype * (wlv+wtv)
           upa(1,i) = 0.5_rtype * erf( wtv/(sqrt(2._rtype)*sigmaw) ) &
                      - 0.5_rtype * erf( wlv/(sqrt(2._rtype)*sigmaw) )

           upu(1, i) = u(j,1)
           upv(1, i) = v(j,1)

           upqc(1,i)  = 0._rtype
           upqt(1,i)  = qt(j,1)  + 0.32_rtype * upw(1,i) * sigmaqt/sigmaw
           upthv(1,i) = thv(j,1) + 0.58_rtype * upw(1,i) * sigmath/sigmaw
           upthl(1,i) = upthv(1,i) / (1._rtype+eps*upqt(1,i))
           upth(1,i)  = upthl(1,i)
           upqv(1,i)  = upqt(1,i)

         enddo

         ! integrate updrafts
         do i=1,nup
           do k=2,nzi

             entexp  = exp(-ent(k,i)*dzt(k))
             entexpu = exp(-ent(k,i)*dzt(k)/3._rtype)

             qtn  = qt(j,k-1) *(1._rtype-entexp ) + upqt (k-1,i)*entexp
             thln = thl(j,k-1)*(1._rtype-entexp ) + upthl(k-1,i)*entexp
             un   = u(j,k-1)  *(1._rtype-entexpu) + upu  (k-1,i)*entexpu
             vn   = v(j,k-1)  *(1._rtype-entexpu) + upv  (k-1,i)*entexpu
             iexh = (1.e5_rtype / p(j,k))**(rgas/cp) ! MKW NOTE: why not just use SHOC exner??

             !Condensation within updrafts, input/output at full levels:
             call condensation_mf(qtn, thln, p(j,k), iexh, &
                                   thvn, qcn, thn, qln, qin)

             ! To avoid singularities w equation has to be computed diferently if wp==0
             b=ggr*(0.5_rtype*(thvn+upthv(k-1,i))/thv(j,k-1)-1._rtype)
             wp = wb*ent(k,i)
             if (wp==0.) then
               wn2 = upw(k-1,i)**2._rtype+2._rtype*wa*b*dzt(k)
             else
               entw = exp(-2._rtype*wp*dzt(k))
               wn2 = entw*upw(k-1,i)**2._rtype+wa*b/(wb*ent(k,i))*(1._rtype-entw)
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
         ! all variables, except areas (moist_a and dry_a) are now multipled by the area
         do k=1,nzi

           ! first sum over all i-updrafts
           do i=1,nup
             if (upqc(k,i)>0.) then
               moist_a(j,k)   = moist_a(j,k)   + upa(k,i)
               moist_w(j,k)   = moist_w(j,k)   + upa(k,i)*upw(k,i)
               moist_qt(j,k)  = moist_qt(j,k)  + upa(k,i)*upqt(k,i)
               moist_thl(j,k) = moist_thl(j,k) + upa(k,i)*upthl(k,i)
               moist_th(j,k)  = moist_th(j,k)  + upa(k,i)*upth(k,i)
               moist_u(j,k)   = moist_u(j,k)   + upa(k,i)*upu(k,i)
               moist_v(j,k)   = moist_v(j,k)   + upa(k,i)*upv(k,i)
               moist_qc(j,k)  = moist_qc(j,k)  + upa(k,i)*upqc(k,i)
             else
               dry_a(j,k)     = dry_a(j,k)     + upa(k,i)
               dry_w(j,k)     = dry_w(j,k)     + upa(k,i)*upw(k,i)
               dry_qt(j,k)    = dry_qt(j,k)    + upa(k,i)*upqt(k,i)
               dry_thl(j,k)   = dry_thl(j,k)   + upa(k,i)*upthl(k,i)
               dry_th(j,k)    = dry_th(j,k)    + upa(k,i)*upth(k,i)
               dry_u(j,k)     = dry_u(j,k)     + upa(k,i)*upu(k,i)
               dry_v(j,k)     = dry_v(j,k)     + upa(k,i)*upv(k,i)
             endif
           enddo

           if ( dry_a(j,k) > 0. ) then
             dry_w(j,k)   = dry_w(j,k)   / dry_a(j,k)
             dry_qt(j,k)  = dry_qt(j,k)  / dry_a(j,k)
             dry_thl(j,k) = dry_thl(j,k) / dry_a(j,k)
             dry_th(j,k)  = dry_th(j,k)  / dry_a(j,k)
             dry_u(j,k)   = dry_u(j,k)   / dry_a(j,k)
             dry_v(j,k)   = dry_v(j,k)   / dry_a(j,k)
           else
             dry_w(j,k)   = 0._rtype
             dry_qt(j,k)  = 0._rtype
             dry_thl(j,k) = 0._rtype
             dry_th(j,k)  = 0._rtype
             dry_u(j,k)   = 0._rtype
             dry_v(j,k)   = 0._rtype
           endif

           if ( moist_a(j,k) > 0._rtype ) then
             moist_w(j,k)   = moist_w(j,k)   / moist_a(j,k)
             moist_qt(j,k)  = moist_qt(j,k)  / moist_a(j,k)
             moist_thl(j,k) = moist_thl(j,k) / moist_a(j,k)
             moist_th(j,k)  = moist_th(j,k)  / moist_a(j,k)
             moist_u(j,k)   = moist_u(j,k)   / moist_a(j,k)
             moist_v(j,k)   = moist_v(j,k)   / moist_a(j,k)
             moist_qc(j,k)  = moist_qc(j,k)  / moist_a(j,k)
           else
             moist_w(j,k)   = 0._rtype
             moist_qt(j,k)  = 0._rtype
             moist_thl(j,k) = 0._rtype
             moist_th(j,k)  = 0._rtype
             moist_u(j,k)   = 0._rtype
             moist_v(j,k)   = 0._rtype
             moist_qc(j,k)  = 0._rtype
           endif

         enddo

         do k=1,nzi
           do i=1,nup
             ae  (j,k) = ae  (j,k) - upa(k,i)
             aw  (j,k) = aw  (j,k) + upa(k,i)*upw(k,i)
             awu (j,k) = awu (j,k) + upa(k,i)*upw(k,i)*upu(k,i)
             awv (j,k) = awv (j,k) + upa(k,i)*upw(k,i)*upv(k,i)
             !aws (k) = aws (k) + upa(k,i)*upw(k,i)*upth(k,i)*cpair
             !aws (k) = aws (k) + upa(k,i)*upw(k,i)*ups(k,i)
             awthl(j,k)= awthl(j,k)+ upa(k,i)*upw(k,i)*upthl(k,i) !*cpair/iexh
             awth(j,k) = awth(j,k) + upa(k,i)*upw(k,i)*upth(k,i) !*cpair/iexh
             awqt(j,k) = awqt(j,k) + upa(k,i)*upw(k,i)*upqt(k,i)
             awqv(j,k) = awqv(j,k) + upa(k,i)*upw(k,i)*upqv(k,i)
             awql(j,k) = awql(j,k) + upa(k,i)*upw(k,i)*upql(k,i)
             awqi(j,k) = awqi(j,k) + upa(k,i)*upw(k,i)*upqi(k,i)
           enddo
         enddo

         ! MKW: not diagnosing fluxes here at present, comment all this out
         !do k=2,nzi
         !  thlflx(j,k)= awthl(j,k) - aw(j,k)*thl(j,k-1) ! MKW NOTE: used to be slflx, but CLUBB works on thl
         !  !sflx( k)= (awth(k) - aw(k)*0.5*(th(k-1)+th(k)) )*cpair/iexh ! not using this since all s/sl stuff is handled in clubb_cam_tend
         !  qtflx(j,k)= awqt(j,k)  - aw(j,k)*qt(j,k-1)
         !enddo
         !thlflx(j,1) = 0._rtype
         !!sflx(kts)  = 0.
         !qtflx(j,1) = 0._rtype

         print*,'max(1-ae)=',maxval(1._rtype-ae(j,:))

       end if  ! ( wthv > 0.0 )
     end do ! j=1,shcol

     ! flip output variables so index 1 = model top (i.e. lowest pressure)
     do k=1,nzi
       dry_a_out(:,nzi-k+1) = dry_a(:,k)
       dry_w_out(:,nzi-k+1) = dry_w(:,k)
       dry_qt_out(:,nzi-k+1) = dry_qt(:,k)
       dry_thl_out(:,nzi-k+1) = dry_thl(:,k)
       dry_u_out(:,nzi-k+1) = dry_u(:,k)
       dry_v_out(:,nzi-k+1) = dry_v(:,k)

       moist_a_out(:,nzi-k+1) = moist_a(:,k)
       moist_w_out(:,nzi-k+1) = moist_w(:,k)
       moist_qt_out(:,nzi-k+1) = moist_qt(:,k)
       moist_thl_out(:,nzi-k+1) = moist_thl(:,k)
       moist_u_out(:,nzi-k+1) = moist_u(:,k)
       moist_v_out(:,nzi-k+1) = moist_v(:,k)
       moist_qc_out(:,nzi-k+1) = moist_qc(:,k)

       ae_out(:,nzi-k+1) = ae(:,k)
       aw_out(:,nzi-k+1) = aw(:,k)
       awthl_out(:,nzi-k+1) = awthl(:,k)
       awqt_out(:,nzi-k+1) = awqt(:,k)
       awql_out(:,nzi-k+1) = awql(:,k)
       awqi_out(:,nzi-k+1) = awqi(:,k)
       awu_out(:,nzi-k+1) = awu(:,k)
       awv_out(:,nzi-k+1) = awv(:,k)

       !thlflx_out(:,nzi-k+1) = thlflx(:,k)
       !qtflx_out(:,nzi-k+1) = qtflx(:,k)
     end do


  end subroutine integrate_mf

  !==============================================================
  ! Linear interpolation to get values on various grids
  ! MKW: copying this from shoc.F90 since I can't get it to work with a "use" statement

subroutine linear_interp(x1,x2,y1,y2,km1,km2,ncol,minthresh)
    implicit none

    integer, intent(in) :: km1, km2
    integer, intent(in) :: ncol
    real(rtype), intent(in) :: x1(ncol,km1), y1(ncol,km1)
    real(rtype), intent(in) :: x2(ncol,km2)
    real(rtype), intent(in) :: minthresh
    real(rtype), intent(out) :: y2(ncol,km2)

    integer :: k1, k2, i

#if 1
    !i = check_grid(x1,x2,km1,km2,ncol)
    if (km1 .eq. km2+1) then
       do k2 = 1,km2
          k1 = k2+1
          do i = 1,ncol
             y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
          end do
       end do
    elseif (km2 .eq. km1+1) then
       k2 = 1
       do i = 1,ncol
          y2(i,k2) = y1(i,1) + (y1(i,2)-y1(i,1))*(x2(i,k2)-x1(i,1))/(x1(i,2)-x1(i,1))
       end do
       do k2 = 2, km2-1
          k1 = k2
          do i = 1,ncol
             y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
          end do
       end do
       k2 = km2
       do i = 1,ncol
          y2(i,k2) = y1(i,km1) + (y1(i,km1)-y1(i,km1-1))*(x2(i,k2)-x1(i,km1))/(x1(i,km1)-x1(i,km1-1))
       end do
    else
       print *,km1,km2
    end if
    do k2 = 1,km2
       do i = 1,ncol
          if (y2(i,k2) .lt. minthresh) then
             y2(i,k2) = minthresh
          endif
       end do
    end do
#else
    do i=1,ncol
       do k2=1,km2
          if( x2(i,k2) <= x1(i,1) ) then
             y2(i,k2) = y1(i,1) + (y1(i,2)-y1(i,1))*(x2(i,k2)-x1(i,1))/(x1(i,2)-x1(i,1))
          elseif( x2(i,k2) >= x1(i,km1) ) then
             y2(i,k2) = y1(i,km1) + (y1(i,km1)-y1(i,km1-1))*(x2(i,k2)-x1(i,km1))/(x1(i,km1)-x1(i,km1-1))
          else
             do k1 = 2,km1
                if( (x2(i,k2)>=x1(i,k1-1)).and.(x2(i,k2)<x1(i,k1)) ) then
                   y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
                endif
             enddo ! end k1 loop
          endif

          if (y2(i,k2) .lt. minthresh) then
             y2(i,k2) = minthresh
          endif

       enddo ! end k2 loop
    enddo ! i loop
#endif

    return

end subroutine linear_interp


  subroutine condensation_mf( qt, thl, p, iex, thv, qc, th, ql, qi)
  !
  ! zero or one condensation for edmf: calculates thv and qc
  !
       use wv_saturation,      only : qsat

       real(rtype),intent(in) :: qt,thl,p,iex
       real(rtype),intent(out):: thv,qc,th,ql,qi

       !local variables
       integer :: niter,i
       real(rtype) :: diff,t,qs,qcold,es,wf

       ! max number of iterations
       niter=50
       ! minimum difference
       diff=2.e-5_rtype

       qc=0._rtype
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
         t = thl/iex+get_alhl(wf)/cp*qc   !as in (4)

         ! qsat, p is in pascal (check!)
         call qsat(t,p,es,qs)
         qcold = qc
         qc = max(0.5_rtype*qc+0.5_rtype*(qt-qs),0._rtype)
         if (abs(qc-qcold)<diff) exit
       enddo

       wf = get_watf(t)
       t = thl/iex+get_alhl(wf)/cp*qc


       call qsat(t,p,es,qs)
       qc = max(qt-qs,0._rtype)
       thv = (thl+get_alhl(wf)/cp*iex*qc)*(1.+eps*(qt-qc)-qc)
       th = t*iex
       qi = qc*(1.-wf)
       ql = qc*wf

       contains

       function get_watf(t)
         real(rtype) :: t,get_watf,tc
         real(rtype), parameter :: &
         tmax=-10._rtype, &
         tmin=-40._rtype

         tc=t-273.16_rtype

         if (tc>tmax) then
           get_watf=1._rtype
         else if (tc<tmin) then
           get_watf=0._rtype
         else
           get_watf=(tc-tmin)/(tmax-tmin);
         end if

       end function get_watf


       function get_alhl(wf)
       !latent heat of the mixture based on water fraction
         real(rtype) :: get_alhl,wf

         get_alhl = wf*lcond+(1._rtype-wf)*(lcond+lice)

       end function get_alhl

  end subroutine condensation_mf

  subroutine calc_mf_vertflux(shcol,nlev,nlevi,aw,awvar,var,varflx)

    implicit none

  ! INPUT VARIABLES
    ! number of SHOC columns
    integer, intent(in) :: shcol
    ! number of midpoint levels
    integer, intent(in) :: nlev
    ! number of interface levels
    integer, intent(in) :: nlevi
    ! Sum plume (a_i*w_i) [m/s]
    real(rtype), intent(in) :: aw(shcol,nlevi)
    ! Sum plume vertical flux of generic variable var (a_i*w_i*var_i) [units vary]
    real(rtype), intent(in) :: awvar(shcol,nlevi)
    ! Input variable on thermo/full grid [units vary]
    real(rtype), intent(in) :: var(shcol,nlevi) ! NOTE: var is interpolated to zi, so has dim nzi

  ! OUTPUT VARIABLE
    real(rtype), intent(out) :: varflx(shcol,nlevi)

  ! INTERNAL VARIABLES
    integer :: i,k

    ! MKW TODO: SHOC has separate subroutines for lower (k=nlevi) and
    !   upper (k=1) boundary conditions. Make these off later if SCREAM
    !   folks want that. Should be very quick.

    ! diagnose MF fluxes
    varflx(:shcol,1) = 0._rtype;
    do k=2,nlev
      do i=1,shcol
        ! MKW NOTE: we may change this to
        varflx(i,k)= awvar(i,k) - aw(i,k)*var(i,k)
      end do
    end do
    varflx(:shcol,nlevi) = 0._rtype;
  end subroutine calc_mf_vertflux

  subroutine compute_tmpi3(nlevi, shcol, dtime, rho_zi, tmpi3)

    !intent-ins
    integer,     intent(in) :: nlevi, shcol
    !time step [s]
    real(rtype), intent(in) :: dtime
    !air density at interfaces [kg/m3]
    real(rtype), intent(in) :: rho_zi(shcol,nlevi)

    !intent-out
    real(rtype), intent(out) :: tmpi3(shcol,nlevi)

    !local vars
    integer :: i, k

    tmpi3(:,1) = 0._rtype
    ! eqn: tmpi3 = dt*g*rho
    do k = 2, nlevi
      do i = 1, shcol
         tmpi3(i,k) = dtime *  ggr*rho_zi(i,k)
      enddo
    enddo

  end subroutine compute_tmpi3

  subroutine Poisson(istart,iend,jstart,jend,mu,poi)
         ! Variable(s)
         integer, intent(in) :: istart,iend,jstart,jend
         real(rtype), dimension(istart:iend,jstart:jend),intent(in) :: mu
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

         implicit none
         integer, allocatable :: seed(:)
         integer :: i, n, un, istat, dt(8), pid
         integer(itype) :: t

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
             t = (dt(1) - 1970) * 365_itype * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_itype * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_itype * 60 * 60 * 1000 &
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
           integer(itype) :: s
           if (s == 0) then
             s = 104729
           else
             s = mod(s, 4294967296_itype)
           end if
           s = mod(s * 279470273_itype, 4294967291_itype)
           lcg = int(mod(s, int(huge(0), itype)), kind(0))
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
  	REAL(rtype), INTENT(IN)    :: mu
  	LOGICAL, INTENT(IN) :: first
    INTEGER             :: ival
  !     ..
  !     .. Local Scalars ..
  	REAL(rtype)          :: b1, b2, c, c0, c1, c2, c3, del, difmuk, e, fk, fx, fy, g,  &
                      omega, px, py, t, u, v, x, xx
  	REAL(rtype), SAVE    :: s, d, p, q, p0
          INTEGER       :: j, k, kflag
  	LOGICAL, SAVE :: full_init
          INTEGER, SAVE :: l, m
  !     ..
  !     .. Local Arrays ..
  	REAL(rtype), SAVE    :: pp(35)
  !     ..
  !     .. Data statements ..
  	REAL(rtype), PARAMETER :: a0 = -.5_rtype, a1 = .3333333_rtype, a2 = -.2500068_rtype, a3 = .2000118_rtype,  &
                  a4 = -.1661269_rtype, a5 = .1421878_rtype, a6 = -.1384794_rtype,   &
                   a7 = .1250060_rtype

  	REAL(rtype), PARAMETER :: fact(10) = (/ 1._rtype, 1._rtype, 2._rtype, 6._rtype, 24._rtype, 120._rtype, 720._rtype, 5040._rtype,  &
              40320._rtype, 362880._rtype /)

  !     ..
  !     .. Executable Statements ..
     IF (mu > 10.0_rtype) THEN
  !     C A S E  A. (RECALCULATION OF S, D, L IF MU HAS CHANGED)

    IF (first) THEN
  s = SQRT(mu)
  d = 6.0_rtype*mu*mu

  !             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
  !             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
  !             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .

  l = mu - 1.1484
  full_init = .false.
    END IF


  !     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE

  	  g = mu + s*random_normal()
  	  IF (g > 0.0_rtype) THEN
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
  		omega = .3989423_rtype/s
  		b1 = .4166667E-1_rtype/mu
  		b2 = .3_rtype*b1*b1
  		c3 = .1428571_rtype*b1*b2
  		c2 = b2 - 15._rtype*c3
  		c1 = b1 - 6._rtype*b2 + 45._rtype*c3
  		c0 = 1._rtype - b1 + 3._rtype*b2 - 15._rtype*c3
  		c = .1069_rtype/mu
  		full_init = .true.
  	  END IF

  	  IF (g < 0.0_rtype) GO TO 50

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
  	  t = 1.8_rtype + SIGN(e, u)
  	  IF (t <= (-.6744_rtype)) GO TO 50
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

  	  80 del = .8333333E-1_rtype/fk
  	  del = del - 4.8_rtype*del*del*del
  	  v = difmuk/fk
  	  IF (ABS(v)>0.25_rtype) THEN
  		px = fk*LOG(1._rtype + v) - difmuk - del
  	  ELSE
  		px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
  	  END IF
  	  py = .3989423_rtype/SQRT(fk)
  	  110 x = (0.5_rtype - difmuk)/s
  	  xx = x*x
  	  fx = -0.5_rtype*xx
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
          REAL(rtype) :: fn_val

  	!     Local variables
  	REAL(rtype)     :: s = 0.449871_rtype, t = -0.386595_rtype, a = 0.19600_rtype, b = 0.25472_rtype,           &
  				r1 = 0.27597_rtype, r2 = 0.27846_rtype, u, v, x, y, q

  	!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

  	DO
  	  CALL RANDOM_NUMBER(u)
  	  CALL RANDOM_NUMBER(v)
  	  v = 1.7156_rtype * (v - 0.5_rtype )

  	!     Evaluate the quadratic form
  	  x = u - s
  	  y = ABS(v) - t
  	  q = x**2._rtype + y*(a*y - b*x)

  	!     Accept P if inside inner ellipse
  	  IF (q < r1) EXIT
  	!     Reject P if outside outer ellipse
  	  IF (q > r2) CYCLE
  	!     Reject P if outside acceptance region
  	  IF (v**2._rtype < -4.0_rtype*LOG(u)*u**2) EXIT
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
          REAL(rtype)  :: fn_val

  	!     Local variable
  	REAL(rtype)  :: r

  	DO
  	  CALL RANDOM_NUMBER(r)
  	  IF (r > 0._rtype) EXIT
  	END DO

  	fn_val = -LOG(r)
  	RETURN

  	END FUNCTION random_exponential

end module edmf

! End of EDMF module. Thanks for visiting.
