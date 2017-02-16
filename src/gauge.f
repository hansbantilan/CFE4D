c----------------------------------------------------------------------
c in polar coordinates t,x, for x in [0,1]
c
c routines associated with the source functions 
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c smooth (C2) transition function from 0 at x=rho1, to 1 at x=rho2
c i.e. trans=0 everywhere for rho1=1.0, rho2=1.0 
c      trans=1 everywhere for rho1=0.0, rho2=0.0
c----------------------------------------------------------------------
        real*8 function trans(x,rho1,rho2)
        implicit none
        real*8 x,rho1,rho2

        real*8 xb

        ! initialize fixed-size variables 
        data xb/0.0/

        !--------------------------------------------------------------
 
        xb=(rho2-x)/(rho2-rho1)

        if (x.ge.rho2) then
          trans=1
        else if (x.ge.rho1) then
          trans=1-xb**3*(6*xb**2-15*xb+10)
        else
          trans=0
        end if

        return
        end

c----------------------------------------------------------------------
c smooth (C-inf) transition function from 0 at x=rho1, to 1 at x=rho2
c i.e. trans2=0 everywhere for rho1=1.0, rho2=1.0 
c      trans2=1 everywhere for rho1=0.0, rho2=0.0
c----------------------------------------------------------------------
        real*8 function trans2(x,rho1,rho2)
        implicit none
        real*8 x,rho1,rho2

        real*8 xb

        ! initialize fixed-size variables 
        data xb/0.0/

        !--------------------------------------------------------------
 
        xb=(x-rho1)/(rho2-rho1)

        if (x.ge.rho2) then
          trans2=1
        else if (x.ge.rho1) then
          trans2=exp(-1/xb)/(exp(-1/xb)+exp(-1/(1-xb)))
        else
          trans2=0
        end if

        return
        end

c----------------------------------------------------------------------
c Evolution of Hb is split into two routines, hb_t_evo and hb_i_evo.
c
c parameters:
c
c gauge : integer controlling which scheme
c
c rho1,rho2: real constants, gauge dependent
c
c current schemes:
c
c Hb_t:
c
c gauge = 0 : fixed gauge
c gauge = 1 : AdS-asymptotics gauge
c
c Hb_i:
c
c gauge = 0 : fixed gauge
c gauge = 1 : AdS-asymptotics gauge
c
c----------------------------------------------------------------------
        subroutine hb_t_evo(res,
     &                    gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                    gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                    gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                    psi_np1,psi_n,psi_nm1,
     &                    Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                    Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                    phi1_np1,phi1_n,phi1_nm1,
     &                    L,x,dt,chr,ex,
     &                    phys_bdy,ghost_width,Nx,
     &                    Hb_t_0,Hb_x_0,
     &                    gauge,t_n,rho1,rho2,rho3,rho4,xi1,xi2,cbulk)
        implicit none
        integer Nx,gauge,phys_bdy(2),ghost_width(2)
        real*8 res(Nx),t_n,t_np1
        real*8 chr(Nx),ex,L
        real*8 x(Nx),dt,rho1,rho2,rho3,rho4,xi1,xi2,cbulk
        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 Hb_t_np1(Nx),Hb_t_n(Nx),Hb_t_nm1(Nx)
        real*8 Hb_x_np1(Nx),Hb_x_n(Nx),Hb_x_nm1(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx)

        real*8 Hb_t_0(Nx),Hb_x_0(Nx)

        integer i

        real*8 x0
        real*8 dx

        real*8 F_t_np1,F_x_np1,F_y_np1
        real*8 Hb_t0,Hb_x0
        real*8 trans
        real*8 f0,g0
        real*8 f1

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------

        dx=x(2)-x(1)

        t_np1=t_n+dt

        do i=1,Nx
          x0=x(i)
          res(i)=0
          if (chr(i).eq.ex) then
            Hb_t_np1(i)=0
          end if
        end do

        !--------------------------------------------------------------
        ! gauge 0
        !--------------------------------------------------------------

        if (gauge.eq.0) return

        !--------------------------------------------------------------
        ! gauge 1
        !--------------------------------------------------------------

        if (gauge.eq.1) then
          do i=2,Nx-1
            if (chr(i).ne.ex) then
                x0=x(i)
                Hb_t0=Hb_t_0(i)

                f0=trans(x0,rho1,rho2)
                f1=trans(x0,rho3,rho4)
                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4

                F_t_np1=gb_tx_np1(i)*1.5d0*f1
     &                 +gb_tx_np1(i)*cbulk*(1-f1)

                if (xi2.le.1e-16) then
                  Hb_t_np1(i)=F_t_np1
                else
                  Hb_t_np1(i)=F_t_np1+(Hb_t0-F_t_np1)*exp(-g0)
                end if
            end if
          end do
          return
        end if

        !--------------------------------------------------------------
        ! otherwise
        !--------------------------------------------------------------

        write(*,*) 'hb_t_evo : error, gauge,',gauge,' unknown'

        return
        end

c-----------------------------------------------------------------------
        subroutine hb_i_evo(res,
     &                    gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                    gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                    gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                    psi_np1,psi_n,psi_nm1,
     &                    Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                    Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                    phi1_np1,phi1_n,phi1_nm1,
     &                    L,x,dt,chr,ex,
     &                    phys_bdy,ghost_width,Nx,
     &                    Hb_t_0,Hb_x_0,
     &                    gauge,t_n,rho1,rho2,rho3,rho4,xi1,xi2,cbulk)

        implicit none
        integer Nx,gauge,phys_bdy(2),ghost_width(2)
        real*8 res(Nx),t_n,t_np1
        real*8 chr(Nx),ex,L
        real*8 x(Nx),dt,rho1,rho2,rho3,rho4,xi1,xi2,cbulk
        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 Hb_t_np1(Nx),Hb_t_n(Nx),Hb_t_nm1(Nx)
        real*8 Hb_x_np1(Nx),Hb_x_n(Nx),Hb_x_nm1(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx)

        real*8 Hb_t_0(Nx),Hb_x_0(Nx)

        integer i

        real*8 x0
        real*8 dx

        real*8 F_t_np1,F_x_np1,F_y_np1
        real*8 Hb_t0,Hb_x0
        real*8 trans
        real*8 f0,g0
        real*8 f1

        real*8 PI
        parameter (PI=3.141592653589793d0) 

        !--------------------------------------------------------------

        dx=x(2)-x(1)

        t_np1=t_n+dt

        do i=1,Nx
          x0=x(i)
          res(i)=0
          if (chr(i).eq.ex) then
            Hb_x_np1(i)=0
          end if
        end do

        !--------------------------------------------------------------
        ! gauge 0
        !--------------------------------------------------------------

        if (gauge.eq.0) return

        !--------------------------------------------------------------
        ! gauge 1
        !--------------------------------------------------------------

        if (gauge.eq.1) then 
          do i=2,Nx-1
            if (chr(i).ne.ex) then
                x0=x(i)
                Hb_x0=Hb_x_0(i)

                f0=trans(x0,rho1,rho2)
                f1=trans(x0,rho3,rho4)
                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4

                F_x_np1=gb_xx_np1(i)*1.5d0*f1
     &                 +gb_xx_np1(i)*cbulk*(1-f1)

                if (xi2.le.1e-16) then
                  Hb_x_np1(i)=F_x_np1
                else
                  Hb_x_np1(i)=F_x_np1+(Hb_x0-F_x_np1)*exp(-g0)
                end if
            end if
          end do
          return
        end if

        !--------------------------------------------------------------
        ! otherwise
        !--------------------------------------------------------------

        write(*,*) 'hb_i_evo : error, gauge,',gauge,' unknown'

        return
        end

c-----------------------------------------------------------------------
