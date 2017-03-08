c----------------------------------------------------------------------
c miscellaneous numerical routines for CFE4D
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c specific x first derivative routine used by excision routines 
c-----------------------------------------------------------------------
        subroutine df1_int_x(f,f_x,dx,i,chr,ex,Nx)

        implicit none
        integer Nx,i
        real*8 f(Nx),chr(Nx),ex,f_x,dx
        logical first
        save first
        data first/.true./

        if (i.eq.1.or.chr(i-1).eq.ex) then
           if (i.le.(Nx-1).and.chr(i+1).ne.ex) then
              f_x=(-f(i)+f(i+1))/dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_x: error in chr stencil (A)'
                 write(*,*) '    i,Nx,dx=',i,Nx,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        else if (i.eq.Nx.or.chr(i+1).eq.ex) then
           if (i.ge.2.and.chr(i-1).ne.ex) then
              f_x=(f(i)-f(i-1))/dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int: error in chr stencil (B)'
                 write(*,*) '    i,Nx,dx=',i,Nx,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        else
           if (chr(i+1).ne.ex.and.chr(i-1).ne.ex) then
              f_x=(f(i+1)-f(i-1))/2/dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int: error in chr stencil (C)'
                 write(*,*) '    i,Nx,dx=',i,Nx,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        end if

        return
        end

c----------------------------------------------------------------------
c the following computes all first derivatives of f
c at a point i at time level n
c
c stencil reduces to first order on excision surface
c----------------------------------------------------------------------
        subroutine df1_int(f_np1,f_n,f_nm1,f_t,f_x,
     &                     dx,dt,i,chr,ex,Nx,name)
        implicit none
        integer Nx,i
        real*8 f_np1(Nx),f_n(Nx),f_nm1(Nx)
        real*8 f_t,f_x,dx,dt,ex,chr(Nx)
        character*(*) name

        logical ltrace,first
        parameter (ltrace=.false.) 
        save first
        data first/.true./
        real*8 f_x_np1,f_x_nm1

        if (chr(i).eq.ex) then
         write(*,*) 'df1_int: error ... point excised'
         stop
        end if

        f_t=(f_np1(i)-f_nm1(i))/2/dt

        call df1_int_x(f_n,f_x,dx,i,chr,ex,Nx)

        if (ltrace) then
           write(*,*) 'df1_int for ',name
           write(*,*) ' f_t=',f_t
           write(*,*) ' f_x=',f_x
        end if

        return
        end

c----------------------------------------------------------------------
c same as df1_int above, but computes all second derivatives as well 
c----------------------------------------------------------------------
        subroutine df2_int(f_np1,f_n,f_nm1,
     &                     f_t,f_x,
     &                     f_tt,f_tx,f_xx,
     &                     dx,dt,i,chr,ex,Nx,name)
        implicit none
        integer Nx,i
        real*8 f_np1(Nx),f_n(Nx),f_nm1(Nx)
        real*8 f_t,f_x,f_tt,f_tx,f_xx
        real*8 dx,dt,ex,chr(Nx)
        character*(*) name
        logical first
        save first
        data first/.true./

        logical ltrace
        parameter (ltrace=.false.)
        real*8 f_x_np1,f_x_nm1

        call df1_int(f_np1,f_n,f_nm1,f_t,f_x,
     &               dx,dt,i,chr,ex,Nx,name)

        f_tt=(f_np1(i)-2*f_n(i)+f_nm1(i))/dt/dt 

        f_xx=0

        if (chr(i).eq.ex) then
          write(*,*) 'df2_int: error ... point excised'
          stop
        end if

        call df1_int_x(f_np1,f_x_np1,dx,i,chr,ex,Nx)
        call df1_int_x(f_nm1,f_x_nm1,dx,i,chr,ex,Nx)
        f_tx=(f_x_np1-f_x_nm1)/2/dt

        !i

        if (i.eq.1.or.(chr(i-1).eq.ex)) then
          if (i.ge.(Nx-1).or.
     &        chr(i+1).eq.ex.or.chr(i+2).eq.ex) then
            if (first) then
              first=.false.
              write(*,*) 'df2_int: error in chr (A)'
              write(*,*) '    i,Nx,dx=',i,Nx,dx
              write(*,*) '    (first error only)'
            end if
            return
          end if
          f_xx=(f_n(i+2)-2*f_n(i+1)+f_n(i))/dx/dx 

        else if (i.eq.Nx.or.(chr(i+1).eq.ex)) then
          if (i.le.2.or.
     &        chr(i-1).eq.ex.or.chr(i-2).eq.ex) then
            if (first) then
              first=.false.
              write(*,*) 'df2_int: error in chr (B)'
              write(*,*) '    i,Nx,dx=',i,Nx,dx
              write(*,*) '    (first error only)'
            end if
            return
          end if
          f_xx=(f_n(i)-2*f_n(i-1)+f_n(i-2))/dx/dx 

        else if (chr(i+1).ne.ex.and.chr(i-1).ne.ex) then
          f_xx=(f_n(i+1)-2*f_n(i)+f_n(i-1))/dx/dx 

        else
          if (first) then
            first=.false.
            write(*,*) 'df2_int: error in chr (C)'
            write(*,*) '    i,Nx,dx=',i,Nx,dx
            write(*,*) '    (first error only)'
          end if
          return
        end if

        if (ltrace) then
           write(*,*) 'df2_int for ',name
           write(*,*) ' f_tt=',f_tt
           write(*,*) ' f_tx=',f_tx
           write(*,*) ' f_xx=',f_xx
        end if

        return
        end


c----------------------------------------------------------------------
c in polar coordinates t,x, for x in [0,1]
c
c initializes f with a 1d gaussian-like profile ... multiplication
c by (1-x^2) is for correct asymptotics for AdS 
c
c f = amp*(1-x^2)*exp (- (r-r0)^2/delta^2) + A  , x > r0
c   = amp*(1-x^2) + A		                , x < r0
c
c where the profile is constructed to be phi1~(1-x^2)
c for correct AdS asymptotics (1-x^2)^2*phi1~(1-x)^3 and regularity at x=0 
c----------------------------------------------------------------------
        subroutine gauss1d(f,amp,r0,delta,xu0,ax,
     &                     amp2,r02,delta2,xu02,ax2,
     &                     L,x,Nx)
        implicit none
        integer i,Nx
        real*8 f(Nx),x(Nx),L
        real*8 amp,r0,delta,ax,xu0,r
        real*8 amp2,r02,delta2,ax2,xu02,r2
        real*8 x0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------

        do i=1,Nx
          f(i)=0
          x0=x(i)
          
          r=sqrt((x0-xu0)**2/ax**2)
          r2=sqrt((x0-xu02)**2/ax2**2)

          ! Gaussian phi=amp*exp(-r^2)/r^3 profile 
          ! remember that phi=phi1*(1-x^2)^2
          if (r.ge.1) then
             f(i)=0
          else if (r.gt.r0) then
             f(i)=amp*exp(-((r-r0)/delta)**2)*(1-x0**2)
          else
             f(i)=amp*(1-x0**2)
          end if

        end do

        return
        end

c-----------------------------------------------------------------------
c for variables with lin_zero_bnd ... zeros residual there
c-----------------------------------------------------------------------
        subroutine lin_zero_bnd_res(f,phys_bdy,all,Nx)
        implicit none
        integer Nx,all
        real*8 f(Nx)
        integer phys_bdy(2)

        integer i,is,ie

        ! initialize fixed-size variables
        data i,is,ie/0,0,0/

        !--------------------------------------------------------------

        if (phys_bdy(1).eq.1.or.all.eq.1) then
          f(2)=0
        end if

        if (phys_bdy(2).eq.1.or.all.eq.1) then
          f(Nx-1)=0
        end if

        return
        end

c-----------------------------------------------------------------------
c in polar coordinates t,x, for x in [0,1]
c
c The following initializes the metric for AdS-Schwarzschild initial data
c
c WARNING: non-horizon-penetrating coordinates
c-----------------------------------------------------------------------
        subroutine init_sch(gb_tt,gb_tx,gb_xx,psi,ief_bh_r0,
     &                      L,phys_bdy,x,chr,ex,Nx)
        implicit none
        integer Nx
        integer phys_bdy(2)
        real*8 gb_tt(Nx),gb_tx(Nx),gb_xx(Nx),psi(Nx)
        real*8 ief_bh_r0
        real*8 chr(Nx),ex,L
        real*8 x(Nx)

        real*8 dx
        
        integer i
        real*8 x0
        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------

        dx=(x(2)-x(1))

        ! initialize gbars 
        do i=2,Nx-1
 
          x0=x(i)

          if (chr(i).ne.ex) then 

            gb_tt(i)=ief_bh_r0*(1-x0**2)/x0
            gb_xx(i)=L**4*ief_bh_r0*(1-x0**2)*(1+x0**2)**2
     &               /(x0**2+L**2*(1-x0**2)**2)
     &               /(x0**3+L**2*(1-x0**2)**2*(x0-ief_bh_r0*(1-x0**2)))
            psi(i)=0

          endif

        end do

        ! x=0 regularization of gbars
        call axi_reg_g(gb_tt,gb_tx,gb_xx,psi,chr,ex,L,x,Nx)

        return
        end

c-----------------------------------------------------------------------
c in polar coordinates t,x, for x in [0,1]
c
c The following initializes the metric for pure AdS
c-----------------------------------------------------------------------
        subroutine init_ads(gb_tt,gb_tx,gb_xx,psi,
     &                      L,phys_bdy,x,chr,ex,Nx)
        implicit none
        integer Nx
        integer phys_bdy(2)
        real*8 gb_tt(Nx),gb_tx(Nx),gb_xx(Nx),psi(Nx)
        real*8 chr(Nx),ex,L
        real*8 x(Nx)

        real*8 dx
        
        integer i
        real*8 x0
        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------

        dx=(x(2)-x(1))

        ! initialize gbars 
        do i=2,Nx-1
 
          x0=x(i)

          if (chr(i).ne.ex) then 

            gb_tt(i)=0
            gb_xx(i)=0
            psi(i)=0

          endif

        end do

        ! x=0 regularization of gbars
        call axi_reg_g(gb_tt,gb_tx,gb_xx,psi,chr,ex,L,x,Nx)

        return
        end

c----------------------------------------------------------------------
c calculates all the tensorial objects in x coordinates at point i
c----------------------------------------------------------------------
        subroutine tensor_init(
     &                  db_txtx_np1,db_txtx_n,db_txtx_nm1,
     &                  db_txty_np1,db_txty_n,db_txty_nm1,
     &                  db_txtz_np1,db_txtz_n,db_txtz_nm1,
     &                  db_txxy_np1,db_txxy_n,db_txxy_nm1,
     &                  db_txxz_np1,db_txxz_n,db_txxz_nm1,
     &                  db_txyz_np1,db_txyz_n,db_txyz_nm1,
     &                  db_tyty_np1,db_tyty_n,db_tyty_nm1,
     &                  db_tytz_np1,db_tytz_n,db_tytz_nm1,
     &                  db_tyxy_np1,db_tyxy_n,db_tyxy_nm1,
     &                  db_tyxz_np1,db_tyxz_n,db_tyxz_nm1,
     &                  db_tyyz_np1,db_tyyz_n,db_tyyz_nm1,
     &                  db_tztz_np1,db_tztz_n,db_tztz_nm1,
     &                  db_tzxy_np1,db_tzxy_n,db_tzxy_nm1,
     &                  db_tzxz_np1,db_tzxz_n,db_tzxz_nm1,
     &                  db_tzyz_np1,db_tzyz_n,db_tzyz_nm1,
     &                  db_xyxy_np1,db_xyxy_n,db_xyxy_nm1,
     &                  db_xyxz_np1,db_xyxz_n,db_xyxz_nm1,
     &                  db_xyyz_np1,db_xyyz_n,db_xyyz_nm1,
     &                  db_xzxz_np1,db_xzxz_n,db_xzxz_nm1,
     &                  db_xzyz_np1,db_xzyz_n,db_xzyz_nm1,
     &                  db_yzyz_np1,db_yzyz_n,db_yzyz_nm1,
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                  gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                  h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                  A_l,A_l_x,Hads_l,
     &                  gamma_ull,gamma_ull_x,
     &                  riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                  einstein_ll,set_ll,
     &                  phi10_x,phi10_xx,
     &                  x,dt,chr,L,ex,Nx,i)
        implicit none

        integer Nx
        real*8 chr(Nx),ex
        real*8 x(Nx),dt,L

        real*8 db_txtx_np1(Nx),db_txtx_n(Nx),db_txtx_nm1(Nx)
        real*8 db_txty_np1(Nx),db_txty_n(Nx),db_txty_nm1(Nx)
        real*8 db_txtz_np1(Nx),db_txtz_n(Nx),db_txtz_nm1(Nx)
        real*8 db_txxy_np1(Nx),db_txxy_n(Nx),db_txxy_nm1(Nx)
        real*8 db_txxz_np1(Nx),db_txxz_n(Nx),db_txxz_nm1(Nx)
        real*8 db_txyz_np1(Nx),db_txyz_n(Nx),db_txyz_nm1(Nx)
        real*8 db_tyty_np1(Nx),db_tyty_n(Nx),db_tyty_nm1(Nx)
        real*8 db_tytz_np1(Nx),db_tytz_n(Nx),db_tytz_nm1(Nx)
        real*8 db_tyxy_np1(Nx),db_tyxy_n(Nx),db_tyxy_nm1(Nx)
        real*8 db_tyxz_np1(Nx),db_tyxz_n(Nx),db_tyxz_nm1(Nx)
        real*8 db_tyyz_np1(Nx),db_tyyz_n(Nx),db_tyyz_nm1(Nx)
        real*8 db_tztz_np1(Nx),db_tztz_n(Nx),db_tztz_nm1(Nx)
        real*8 db_tzxy_np1(Nx),db_tzxy_n(Nx),db_tzxy_nm1(Nx)
        real*8 db_tzxz_np1(Nx),db_tzxz_n(Nx),db_tzxz_nm1(Nx)
        real*8 db_tzyz_np1(Nx),db_tzyz_n(Nx),db_tzyz_nm1(Nx)
        real*8 db_xyxy_np1(Nx),db_xyxy_n(Nx),db_xyxy_nm1(Nx)
        real*8 db_xyxz_np1(Nx),db_xyxz_n(Nx),db_xyxz_nm1(Nx)
        real*8 db_xyyz_np1(Nx),db_xyyz_n(Nx),db_xyyz_nm1(Nx)
        real*8 db_xzxz_np1(Nx),db_xzxz_n(Nx),db_xzxz_nm1(Nx)
        real*8 db_xzyz_np1(Nx),db_xzyz_n(Nx),db_xzyz_nm1(Nx)
        real*8 db_yzyz_np1(Nx),db_yzyz_n(Nx),db_yzyz_nm1(Nx)
        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 Hb_t_np1(Nx),Hb_t_n(Nx),Hb_t_nm1(Nx)
        real*8 Hb_x_np1(Nx),Hb_x_n(Nx),Hb_x_nm1(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx)

        integer a,b,c,d,e,f,g,h

        integer i

        real*8 x0     
        real*8 dx

        real*8 grad_phi1_sq

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,chi,theta)
        !--------------------------------------------------------------
        real*8 g0_ll(4,4),g0_uu(4,4)
        real*8 g0_ll_x(4,4,4),g0_uu_x(4,4,4),g0_ll_xx(4,4,4,4)
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4),gads_ll_xx(4,4,4,4)
        real*8 h0_ll(4,4),h0_uu(4,4)
        real*8 h0_ll_x(4,4,4),h0_uu_x(4,4,4),h0_ll_xx(4,4,4,4)
        real*8 gamma_ull(4,4,4),gamma_ull_x(4,4,4,4)
        real*8 riemann_ulll(4,4,4,4)
        real*8 ricci_ll(4,4),ricci_lu(4,4),ricci
        real*8 einstein_ll(4,4),set_ll(4,4)
        real*8 Hads_l(4),A_l(4),A_l_x(4,4)
        real*8 phi10_x(4),phi10_xx(4,4)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t,gb_tt_x
        real*8 gb_tt_tt,gb_tt_tx,gb_tt_xx
        real*8 gb_tx_t,gb_tx_x
        real*8 gb_tx_tt,gb_tx_tx,gb_tx_xx
        real*8 gb_xx_t,gb_xx_x
        real*8 gb_xx_tt,gb_xx_tx,gb_xx_xx
        real*8 psi_t,psi_x
        real*8 psi_tt,psi_tx,psi_xx
        real*8 phi1_t,phi1_x
        real*8 phi1_tt,phi1_tx,phi1_xx
  
        real*8 gb_tt0,gb_tx0,gb_xx0,psi0,phi10

        real*8 g0_tt_ads_x,g0_tt_ads_xx
        real*8 g0_xx_ads_x,g0_xx_ads_xx
        real*8 g0_psi_ads_x,g0_psi_ads_xx

        real*8 g0_tt_ads0,g0_xx_ads0,g0_psi_ads0

        real*8 Hb_t_t,Hb_t_x
        real*8 Hb_x_t,Hb_x_x

        real*8 Hb_t0,Hb_x0

        real*8 lambda4

!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))

        x0=x(i)

        ! AdS4D cosmological constant
        !(lambda4=-(n-1)(n-2)/2/L^2) for n=4 dimensional AdS)
        lambda4=-3/L/L

        ! set gads values
        g0_tt_ads0 =-cosh(2*x0/(1-x0**2))**2
        g0_xx_ads0 =(-3d0/lambda4)*4*(1+x0**2)**2/(1-x0**2)**4
        g0_psi_ads0=(-3d0/lambda4)*sinh(2*x0/(1-x0**2))**2

        ! set gbar values
        gb_tt0=gb_tt_n(i)
        gb_tx0=gb_tx_n(i)
        gb_xx0=gb_xx_n(i)
        psi0  =psi_n(i)

        ! set hbar values
        Hb_t0=Hb_t_n(i)
        Hb_x0=Hb_x_n(i)

        ! set phi1 value
        phi10=phi1_n(i)

        ! set gads derivatives
        g0_tt_ads_x  =(-2*(1+x0**2)*sinh((4*x0)/(1-x0**2)))
     &               /(-1+x0**2)**2
        g0_tt_ads_xx =(-4*(2*(1+x0**2)**2*cosh((4*x0)/(-1+x0**2))
     &               +x0*(-3 + 2*x0**2 + x0**4)
     &               *sinh((4*x0)/(-1 + x0**2))))/(-1 + x0**2)**4
        g0_xx_ads_x  =(-16*L**2*x0*(1+x0**2)*(3+x0**2))
     &               /(-1+x0**2)**5
        g0_xx_ads_xx =(16*L**2*(3+39*x0**2+33*x0**4+5*x0**6))
     &               /(-1+x0**2)**6
        g0_psi_ads_x =(2*L**2*(1+x0**2)*sinh((4*x0)/(1-x0**2)))
     &               /(-1+x0**2)**2
        g0_psi_ads_xx=(4*L**2*(2*(1+x0**2)**2
     &               *cosh((4*x0)/(-1+x0**2))
     &               +x0*(-3 + 2*x0**2 + x0**4)
     &               *sinh((4*x0)/(-1 + x0**2))))/(-1 + x0**2)**4
 
        ! calculate gbar derivatives
        call df2_int(gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &       gb_tt_t,gb_tt_x,
     &       gb_tt_tt,gb_tt_tx,gb_tt_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')
        call df2_int(gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &       gb_tx_t,gb_tx_x,
     &       gb_tx_tt,gb_tx_tx,gb_tx_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')
        call df2_int(gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &       gb_xx_t,gb_xx_x,
     &       gb_xx_tt,gb_xx_tx,gb_xx_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')
        call df2_int(psi_np1,psi_n,psi_nm1,
     &       psi_t,psi_x,
     &       psi_tt,psi_tx,psi_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')

        ! calculate hbar derivatives
        call df1_int(Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &       Hb_t_t,Hb_t_x,
     &       dx,dt,i,chr,ex,Nx,'Hb_t')
        call df1_int(Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &       Hb_x_t,Hb_x_x,
     &       dx,dt,i,chr,ex,Nx,'Hb_x')

        ! calculate phi1 derivatives
        call df2_int(phi1_np1,phi1_n,phi1_nm1,
     &       phi1_t,phi1_x,
     &       phi1_tt,phi1_tx,phi1_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')

!NEED TO UPDATE BELOW GIVEN NEW G0_ADS ABOVE!
        ! give values to the metric, using sin(chi)=1,sin(theta)=1 w.l.o.g 
        !(considering chi,theta-independent case, so chi=pi/2,theta=pi/2 will do)
        g0_ll(1,1)=g0_tt_ads0+gb_tt0
        g0_ll(1,2)=           gb_tx0
        g0_ll(2,2)=g0_xx_ads0+gb_xx0
        g0_ll(3,3)=g0_psi_ads0+psi0*x0**2
        g0_ll(4,4)=g0_psi_ads0+psi0*x0**2

        g0_uu(1,1)=
     &          (g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4))
        g0_uu(1,2)=
     &         (-g0_ll(1,2)*g0_ll(3,3)*g0_ll(4,4))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4))
        g0_uu(2,2)=
     &          (g0_ll(1,1)*g0_ll(3,3)*g0_ll(4,4))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4))
        g0_uu(3,3)=
     &          (g0_ll(1,1)*g0_ll(2,2)*g0_ll(4,4)
     &                  -g0_ll(1,2)**2*g0_ll(4,4))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4))
        g0_uu(4,4)=
     &          (g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)
     &                  -g0_ll(1,2)**2*g0_ll(3,3))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4))

        g0_ll_x(1,1,1)   =0
     &                   +gb_tt_t
        g0_ll_x(1,1,2)   =g0_tt_ads_x
     &                   +gb_tt_x
        g0_ll_xx(1,1,1,1)=0
     &                   +gb_tt_tt
        g0_ll_xx(1,1,1,2)=0
     &                   +gb_tt_tx
        g0_ll_xx(1,1,2,2)=g0_tt_ads_xx
     &                   +gb_tt_xx

        g0_ll_x(1,2,1)   =0
     &                   +gb_tx_t
        g0_ll_x(1,2,2)   =0
     &                   +gb_tx_x
        g0_ll_xx(1,2,1,1)=0
     &                   +gb_tx_tt
        g0_ll_xx(1,2,1,2)=0
     &                   +gb_tx_tx
        g0_ll_xx(1,2,2,2)=0
     &                   +gb_tx_xx

        g0_ll_x(2,2,1)   =0
     &                   +gb_xx_t
        g0_ll_x(2,2,2)   =g0_xx_ads_x
     &                   +gb_xx_x
        g0_ll_xx(2,2,1,1)=0
     &                   +gb_xx_tt
        g0_ll_xx(2,2,1,2)=0
     &                   +gb_xx_tx
        g0_ll_xx(2,2,2,2)=g0_xx_ads_xx
     &                   +gb_xx_xx

        g0_ll_x(3,3,1)   =0
     &                   +psi_t*x0**2
        g0_ll_x(3,3,2)   =g0_psi_ads_x
     &                   +psi_x*x0**2
     &                   +psi0*(2*x0)
        g0_ll_xx(3,3,1,1)=0
     &                   +psi_tt*x0**2
        g0_ll_xx(3,3,1,2)=0
     &                   +psi_tx*x0**2
     &                   +psi_t*(2*x0)
        g0_ll_xx(3,3,2,2)=g0_psi_ads_xx
     &                   +psi_xx*x0**2
     &                   +psi_x*(2*x0)
     &                   +psi_x*(2*x0)
     &                   +psi0*(2)

        g0_ll_x(4,4,1)   =0
     &                   +psi_t*x0**2
        g0_ll_x(4,4,2)   =g0_psi_ads_x
     &                   +psi_x*x0**2
     &                   +psi0*(2*x0)
        g0_ll_xx(4,4,1,1)=0
     &                   +psi_tt*x0**2
        g0_ll_xx(4,4,1,2)=0
     &                   +psi_tx*x0**2
     &                   +psi_t*(2*x0)
        g0_ll_xx(4,4,2,2)=g0_psi_ads_xx
     &                   +psi_xx*x0**2
     &                   +psi_x*(2*x0)
     &                   +psi_x*(2*x0)
     &                   +psi0*(2)
        g0_ll_xx(4,4,3,3)=-2*g0_psi_ads0          ! WARNING: from sin^2chi factor in pure ads term

        do a=1,3
          do b=a+1,4
            g0_ll(b,a)=g0_ll(a,b)
            g0_uu(b,a)=g0_uu(a,b) 
            do c=1,4
              g0_ll_x(b,a,c)=g0_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              g0_uu_x(a,b,c)=0
              do d=1,4
                do e=1,4
                  g0_uu_x(a,b,c)=g0_uu_x(a,b,c)
     &                          -g0_ll_x(d,e,c)
     &                           *g0_uu(a,d)*g0_uu(b,e)
                end do
     &  
              end do
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              do d=1,4
                g0_ll_xx(a,b,c,d)=
     &             g0_ll_xx(min(a,b),max(a,b),min(c,d),max(c,d))
              end do
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              gamma_ull(a,b,c)=0
              do d=1,4
                gamma_ull(a,b,c)=gamma_ull(a,b,c)
     &                          +0.5d0*g0_uu(a,d)
     &                                *(g0_ll_x(c,d,b)
     &                                 -g0_ll_x(b,c,d)
     &                                 +g0_ll_x(d,b,c))
              end do
            end do
          end do
        end do

        ! give values to the ads metric, using sin(chi)=1,sin(theta)=1 w.l.o.g 
        !(considering chi,theta-independent case, so chi=pi/2,theta=pi/2 will do)
        gads_ll(1,1)=g0_tt_ads0
        gads_ll(2,2)=g0_xx_ads0
        gads_ll(3,3)=g0_psi_ads0
        gads_ll(4,4)=g0_psi_ads0

        gads_uu(1,1)=1/g0_tt_ads0
        gads_uu(2,2)=1/g0_xx_ads0
        gads_uu(3,3)=1/g0_psi_ads0
        gads_uu(4,4)=1/g0_psi_ads0

        gads_ll_x(1,1,2)   =g0_tt_ads_x  
        gads_ll_xx(1,1,2,2)=g0_tt_ads_xx 
        gads_ll_x(2,2,2)   =g0_xx_ads_x  
        gads_ll_xx(2,2,2,2)=g0_xx_ads_xx 
        gads_ll_x(3,3,2)   =g0_psi_ads_x 
        gads_ll_xx(3,3,2,2)=g0_psi_ads_xx
        gads_ll_x(4,4,2)   =g0_psi_ads_x 
        gads_ll_xx(4,4,2,2)=g0_psi_ads_xx
        gads_ll_xx(4,4,3,3)=-2*g0_psi_ads0          ! WARNING: from sin^2chi factor in pure ads term
                
        do a=1,3
          do b=a+1,4
            gads_ll(b,a)=gads_ll(a,b)
            gads_uu(b,a)=gads_uu(a,b)
            do c=1,4
              gads_ll_x(b,a,c)=gads_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              gads_uu_x(a,b,c)=
     &              -gads_ll_x(1,1,c)*gads_uu(a,1)*gads_uu(b,1)
     &              -gads_ll_x(1,2,c)*(gads_uu(a,1)*gads_uu(b,2)
     &                               +gads_uu(a,2)*gads_uu(b,1))
     &              -gads_ll_x(1,3,c)*(gads_uu(a,1)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,1))
     &              -gads_ll_x(1,4,c)*(gads_uu(a,1)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,1))
     &              -gads_ll_x(2,2,c)*gads_uu(a,2)*gads_uu(b,2)
     &              -gads_ll_x(2,3,c)*(gads_uu(a,2)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,2))
     &              -gads_ll_x(2,4,c)*(gads_uu(a,2)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,2))
     &              -gads_ll_x(3,3,c)*gads_uu(a,3)*gads_uu(b,3)
     &              -gads_ll_x(3,4,c)*(gads_uu(a,3)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,3))
     &              -gads_ll_x(4,4,c)*gads_uu(a,4)*gads_uu(b,4)
            end do
          end do
        end do

        ! give values to the metric deviation, using sin(chi)=1,sin(theta)=1 w.l.o.g 
        !(considering chi,theta-independent case, so chi=pi/2,theta=pi/2 will do)
        h0_ll(1,1)=gb_tt0
        h0_ll(1,2)=gb_tx0
        h0_ll(2,2)=gb_xx0
        h0_ll(3,3)=psi0*x0**2
        h0_ll(4,4)=psi0*x0**2
        
        h0_uu(1,1)=g0_uu(1,1)-gads_uu(1,1)
        h0_uu(1,2)=g0_uu(1,2)
        h0_uu(2,2)=g0_uu(2,2)-gads_uu(2,2)
        h0_uu(3,3)=g0_uu(3,3)-gads_uu(3,3)
        h0_uu(4,4)=g0_uu(4,4)-gads_uu(4,4)

        h0_ll_x(1,1,1)   =g0_ll_x(1,1,1)-gads_ll_x(1,1,1)
        h0_ll_x(1,1,2)   =g0_ll_x(1,1,2)-gads_ll_x(1,1,2)
        h0_ll_xx(1,1,1,1)=g0_ll_xx(1,1,1,1)-gads_ll_xx(1,1,1,1)
        h0_ll_xx(1,1,1,2)=g0_ll_xx(1,1,1,2)-gads_ll_xx(1,1,1,2)
        h0_ll_xx(1,1,2,2)=g0_ll_xx(1,1,2,2)-gads_ll_xx(1,1,2,2)

        h0_ll_x(1,2,1)   =g0_ll_x(1,2,1)-gads_ll_x(1,2,1)
        h0_ll_x(1,2,2)   =g0_ll_x(1,2,2)-gads_ll_x(1,2,2)
        h0_ll_xx(1,2,1,1)=g0_ll_xx(1,2,1,1)-gads_ll_xx(1,2,1,1)
        h0_ll_xx(1,2,1,2)=g0_ll_xx(1,2,1,2)-gads_ll_xx(1,2,1,2)
        h0_ll_xx(1,2,2,2)=g0_ll_xx(1,2,2,2)-gads_ll_xx(1,2,2,2)

        h0_ll_x(2,2,1)   =g0_ll_x(2,2,1)-gads_ll_x(2,2,1)
        h0_ll_x(2,2,2)   =g0_ll_x(2,2,2)-gads_ll_x(2,2,2)
        h0_ll_xx(2,2,1,1)=g0_ll_xx(2,2,1,1)-gads_ll_xx(2,2,1,1)
        h0_ll_xx(2,2,1,2)=g0_ll_xx(2,2,1,2)-gads_ll_xx(2,2,1,2)
        h0_ll_xx(2,2,2,2)=g0_ll_xx(2,2,2,2)-gads_ll_xx(2,2,2,2)

        h0_ll_x(3,3,1)   =g0_ll_x(3,3,1)-gads_ll_x(3,3,1)
        h0_ll_x(3,3,2)   =g0_ll_x(3,3,2)-gads_ll_x(3,3,2)
        h0_ll_xx(3,3,1,1)=g0_ll_xx(3,3,1,1)-gads_ll_xx(3,3,1,1)
        h0_ll_xx(3,3,1,2)=g0_ll_xx(3,3,1,2)-gads_ll_xx(3,3,1,2)
        h0_ll_xx(3,3,2,2)=g0_ll_xx(3,3,2,2)-gads_ll_xx(3,3,2,2)

        h0_ll_x(4,4,1)   =g0_ll_x(4,4,1)-gads_ll_x(4,4,1)
        h0_ll_x(4,4,2)   =g0_ll_x(4,4,2)-gads_ll_x(4,4,2)
        h0_ll_xx(4,4,1,1)=g0_ll_xx(4,4,1,1)-gads_ll_xx(4,4,1,1)
        h0_ll_xx(4,4,1,2)=g0_ll_xx(4,4,1,2)-gads_ll_xx(4,4,1,2)
        h0_ll_xx(4,4,2,2)=g0_ll_xx(4,4,2,2)-gads_ll_xx(4,4,2,2)
        h0_ll_xx(4,4,3,3)=g0_ll_xx(4,4,3,3)-gads_ll_xx(4,4,3,3)

        h0_uu_x(1,1,1)=g0_uu_x(1,1,1)-gads_uu_x(1,1,1)
        h0_uu_x(1,1,2)=g0_uu_x(1,1,2)-gads_uu_x(1,1,2)

        h0_uu_x(1,2,1)=g0_uu_x(1,2,1)-gads_uu_x(1,2,1)
        h0_uu_x(1,2,2)=g0_uu_x(1,2,2)-gads_uu_x(1,2,2)

        h0_uu_x(2,2,1)=g0_uu_x(2,2,1)-gads_uu_x(2,2,1)
        h0_uu_x(2,2,2)=g0_uu_x(2,2,2)-gads_uu_x(2,2,2)

        h0_uu_x(3,3,1)=g0_uu_x(3,3,1)-gads_uu_x(3,3,1)
        h0_uu_x(3,3,2)=g0_uu_x(3,3,2)-gads_uu_x(3,3,2)

        h0_uu_x(4,4,1)=g0_uu_x(4,4,1)-gads_uu_x(4,4,1)
        h0_uu_x(4,4,2)=g0_uu_x(4,4,2)-gads_uu_x(4,4,2)

        do a=1,3
          do b=a+1,4
            h0_ll(b,a)=h0_ll(a,b)
            h0_uu(b,a)=h0_uu(a,b)
            do c=1,4
              h0_ll_x(b,a,c)=h0_ll_x(a,b,c)
              h0_uu_x(b,a,c)=h0_uu_x(a,b,c)
            end do
          end do
        end do

        ! give values to the gh source functions
        A_l(1)=Hb_t0*(1-x0**2)
        A_l(2)=Hb_x0*(1-x0**2)

        A_l_x(1,1)=Hb_t_t*(1-x0**2)
        A_l_x(1,2)=Hb_t_x*(1-x0**2)
     &            -2*x0*Hb_t0

        A_l_x(2,1)=Hb_x_t*(1-x0**2)
        A_l_x(2,2)=Hb_x_x*(1-x0**2)
     &            -2*x0*Hb_x0

        ! give values to the ads gh source functions
        Hads_l(1)=0
        Hads_l(2)=(2-2*x0**2+8*x0**4)/(1-x0**2+x0**6-x0**8)/x0

        ! give values to the scalar field
        phi10_x(1)=phi1_t*(1-x0**2)**2
        phi10_x(2)=phi1_x*(1-x0**2)**2
     &            +phi10*(-4*x0)*(1-x0**2)
        phi10_x(3)=0
        phi10_x(4)=0

        phi10_xx(1,1)=phi1_tt*(1-x0**2)**2
        phi10_xx(1,2)=phi1_tx*(1-x0**2)**2
     &               +phi1_t*(-4*x0)*(1-x0**2)
        phi10_xx(1,3)=0
        phi10_xx(1,4)=0
        phi10_xx(2,2)=phi1_xx*(1-x0**2)**2
     &               +phi1_x*(2)*(-4*x0)*(1-x0**2)
     &               +phi10*(-4+12*x0**2)
        phi10_xx(2,3)=0
        phi10_xx(2,4)=0
        phi10_xx(3,3)=0                
        phi10_xx(3,4)=0
        phi10_xx(4,4)=0

        do a=1,3
          do b=a+1,4
            phi10_xx(b,a)=phi10_xx(a,b)
          end do
        end do

        ! calculate Christoffel symbol derivatives at point i
        !(gamma^a_bc,e = 1/2 g^ad_,e(g_bd,c  + g_cd,b  - g_bc,d)
        !              +   1/2 g^ad(g_bd,ce + g_cd,be - g_bc,de))
        !(WARNING: 
        ! this second derivative needs more info on chi,theta, so 
        ! cannot use the g0_ll_x, g0_uu_x, g0_ll_xx set by chi=pi/2,theta=pi/2, 
        ! only (from pure AdS piece):
        ! gamma_ull(4,4,3)=cotchi
        ! gamma_ull(2,4,4)=-1/L^2*x/(1-x)*(L^2*(1-x)^2+x^2)*sinchi^2
        ! gamma_ull(3,4,4)=-coschi*sinchi
        ! gamma_ull(4,3,4)=cotchi
        ! are affected in this AdS4D case, so hardcoded these at the end)
        do a=1,4
          do b=1,4
            do c=1,4
              do e=1,4
                gamma_ull_x(a,b,c,e)=0
                do d=1,4
                  gamma_ull_x(a,b,c,e)=gamma_ull_x(a,b,c,e)
     &              +0.5d0*g0_uu_x(a,d,e)*(g0_ll_x(b,d,c)+
     &                     g0_ll_x(c,d,b)-g0_ll_x(b,c,d))
     &              +0.5d0*g0_uu(a,d)*(g0_ll_xx(b,d,c,e)+
     &                     g0_ll_xx(c,d,b,e)-g0_ll_xx(b,c,d,e))
                end do
              end do
            end do
          end do
        end do
        gamma_ull_x(4,4,3,3)=-1 
        gamma_ull_x(2,4,4,3)=0
        gamma_ull_x(3,4,4,3)=1
        gamma_ull_x(4,3,4,3)=-1

        ! calculate Riemann tensor at point i
        !(R^a_bcd =gamma^a_bd,c - gamma^a_bc,d
        !          +gamma^a_ce gamma^e_bd - gamma^a_de gamma^e_bc)
        do a=1,4
          do b=1,4
            do c=1,4
              do d=1,4
                riemann_ulll(a,b,c,d)=
     &                gamma_ull_x(a,b,d,c)-gamma_ull_x(a,b,c,d)
                do e=1,4
                   riemann_ulll(a,b,c,d)=riemann_ulll(a,b,c,d)
     &               +gamma_ull(a,c,e)*gamma_ull(e,b,d)
     &               -gamma_ull(a,d,e)*gamma_ull(e,b,c)
                end do
              end do
            end do
          end do
        end do

        ! calculate Ricci tensor at point i
        !(R_bd = R^a_bad)
        do b=1,4
          do d=1,4
            ricci_ll(b,d)=0
            do a=1,4
              ricci_ll(b,d)=ricci_ll(b,d)+riemann_ulll(a,b,a,d)
            end do
          end do
        end do

        ! calculate raised Ricci tensor at point i
        !(R_a^b = R_ad g^db)
        do a=1,4
          do b=1,4
            ricci_lu(a,b)=0
            do d=1,4
              ricci_lu(a,b)=ricci_lu(a,b)+ricci_ll(a,d)*g0_uu(d,b)
            end do
          end do
        end do

        ! calculate Ricci scalar
        !(R = R_a^a)
        ricci=0
        do a=1,4
          ricci=ricci+ricci_lu(a,a)
        end do
  
        ! calculates Einstein tensor at point i
        !(G_ab = R_ab - 1/2 R g_ab)
        do a=1,4
          do b=1,4
            einstein_ll(a,b)=ricci_ll(a,b)-0.5d0*ricci*g0_ll(a,b)
          end do
        end do

        ! calculates stress-energy tensor at point i 
        !(T_ab = 2*phi1,a phi1,b - (phi1,c phi1,d) g^cd g_ab + ...)
        grad_phi1_sq=0
        do a=1,4
          do b=1,4
            grad_phi1_sq=grad_phi1_sq
     &                  +phi10_x(a)*phi10_x(b)*g0_uu(a,b)
          end do
        end do

        do a=1,4
          do b=1,4
            set_ll(a,b)=
     &            phi10_x(a)*phi10_x(b)
     &           -g0_ll(a,b)*(grad_phi1_sq/2)
          end do
        end do

        return
        end
