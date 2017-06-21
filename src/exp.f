c----------------------------------------------------------------------
c in cartesian coordinates t,x,y,z for x in [0,1], y in [0,1], z in [-1,1]
c
c support routines for apph.c
c
c (fixed using correct g=gads+(1-x^2)gb factorization)
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c returns true if a given point is not excised, and can be updated if
c in the *ghost region* of the grid
c----------------------------------------------------------------------
        logical function can_calc_ex(chr,i,j,Nx,Ny,Nz,ex)
        implicit none
        integer i,j,Nx,Ny,Nz
        real*8 chr(Nx,Ny),ex

        integer i1,j1

        ! initialize fixed-size variables
        data i1,j1/0,0/

        !--------------------------------------------------------------

        can_calc_ex=.true.

        if (chr(i,j).ne.ex.and.(i.gt.2.and.i.lt.Nx-1.and.
     &                          j.gt.2.and.j.lt.Ny-1)
     &                         ) return

        do i1=max(1,min(Nx-2,i-1)),min(Nx,max(3,i+1))
           do j1=max(1,min(Ny-2,j-1)),min(Ny,max(3,j+1))
              if (chr(i1,j1).eq.ex) then
                 can_calc_ex=.false.
                 return
              end if
           end do
        end do

        return 
        end

c----------------------------------------------------------------------
c the following computes the outgoing null expansion (theta) 
c
c f is the level-surface function, and [is0..ie0,js0..je0] specifies
c the region overwhich theta should be computed
c
c is_ex is set to 1 if couldn't compute expansion an some point due
c to closeness of an excised region.
c----------------------------------------------------------------------
        subroutine calc_exp(theta,f,is0,ie0,js0,je0,is_ex,
     &                    gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                    gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                    gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                    gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                    gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                    gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                    psi_np1,psi_n,psi_nm1,
     &                    L,x,y,z,dt,chr,ex,do_ex,
     &                    Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,is0,ie0,js0,je0,do_ex
        integer is_ex
        real*8 theta(Nx,Ny),f(Nx,Ny),chr(Nx,Ny),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L

        real*8 gb_tt_np1(Nx,Ny),gb_tt_n(Nx,Ny),gb_tt_nm1(Nx,Ny)
        real*8 gb_tx_np1(Nx,Ny),gb_tx_n(Nx,Ny),gb_tx_nm1(Nx,Ny)
        real*8 gb_ty_np1(Nx,Ny),gb_ty_n(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xx_n(Nx,Ny),gb_xx_nm1(Nx,Ny)
        real*8 gb_xy_np1(Nx,Ny),gb_xy_n(Nx,Ny),gb_xy_nm1(Nx,Ny)
        real*8 gb_yy_np1(Nx,Ny),gb_yy_n(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 psi_np1(Nx,Ny),psi_n(Nx,Ny),psi_nm1(Nx,Ny)

        integer i,j,is,ie,js,je
        integer a,b,c,d,e

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 x0,y0
        real*8 dx,dy

        real*8 zeros(Nx,Ny)

        real*8 f_x,f_y,f_xx,f_xy,f_yy
        real*8 tmp1,tmp2,tmp3,tmp4,tmp5

        real*8 n_l(5),s_l(5)
        real*8 n_u(5),s_u(5)
        real*8 n_l_x(5,5),n_u_x(5,5),s_l_x(5,5)
        real*8 f0_x(5) 
        real*8 f0_xx(5,5)
        real*8 gam_uu(5,5),sig_uu(5,5)
        real*8 gam_uu_x(5,5,5)
        real*8 normsusq

        real*8 g0gamfx(5)
        real*8 nufx,nuxfx(5),gamxfxfx(5)

        real*8 theta_ads

        integer do_ex_n
        logical ltrace,any_ex,can_calc_ex
        parameter (ltrace=.false.,do_ex_n=0)

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,w,y,z)
        !--------------------------------------------------------------
        real*8 g0_ll(5,5),g0_uu(5,5)
        real*8 g0_ll_x(5,5,5),g0_uu_x(5,5,5),g0_ll_xx(5,5,5,5)
        real*8 gads_ll(5,5),gads_uu(5,5)
        real*8 gads_ll_x(5,5,5),gads_uu_x(5,5,5),gads_ll_xx(5,5,5,5)
        real*8 h0_ll(5,5),h0_uu(5,5)
        real*8 h0_ll_x(5,5,5),h0_uu_x(5,5,5),h0_ll_xx(5,5,5,5)
        real*8 gamma_ull(5,5,5),gamma_ull_x(5,5,5,5)
        real*8 riemann_ulll(5,5,5,5)
        real*8 ricci_ll(5,5),ricci_lu(5,5),ricci
        real*8 einstein_ll(5,5),set_ll(5,5)
        real*8 Hads_l(5),A_l(5),A_l_x(5,5)
        real*8 phi10_x(5),phi10_xx(5,5)

        !--------------------------------------------------------------
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data i,j/0,0/
        data is,ie,js,je/0,0,0,0/

        data x0,y0/0.0,0.0/
        data dx,dy/0.0,0.0/

        data f_x,f_y/0.0,0.0/
        data f_xx,f_xy/0.0,0.0/
        data f_yy/0.0/
        data tmp1,tmp2,tmp3,tmp4,tmp5/0.0,0.0,0.0,0.0,0.0/

        data n_l,s_l/5*0.0,5*0.0/
        data n_u,s_u/5*0.0,5*0.0/
        data n_l_x,n_u_x,s_l_x/25*0.0,25*0.0,25*0.0/
        data f0_x/5*0.0/
        data f0_xx/25*0.0/
        data gam_uu,sig_uu/25*0.0,25*0.0/
        data gam_uu_x/125*0.0/
        data normsusq/0.0/

        data g0gamfx/5*0.0/
        data nufx,nuxfx,gamxfxfx/0.0,5*0.0,5*0.0/

        data g0_ll,g0_uu/25*0.0,25*0.0/
        data gads_ll,gads_uu/25*0.0,25*0.0/
        data h0_ll,h0_uu/25*0.0,25*0.0/
        data gamma_ull/125*0.0/
        data gamma_ull_x/625*0.0/

        data g0_ll_x,g0_uu_x/125*0.0,125*0.0/
        data gads_ll_x,gads_uu_x/125*0.0,125*0.0/
        data h0_ll_x,h0_uu_x/125*0.0,125*0.0/

        data g0_ll_xx/625*0.0/
        data gads_ll_xx/625*0.0/
        data h0_ll_xx/625*0.0/

        data ricci/0.0/
        data ricci_ll,ricci_lu/25*0.0,25*0.0/
        data einstein_ll,set_ll/25*0.0,25*0.0/
        data riemann_ulll/625*0.0/

        data A_l,Hads_l/5*0.0,5*0.0/
        data A_l_x/25*0.0/

        data phi10_x/5*0.0/
        data phi10_xx/25*0.0/

!----------------------------------------------------------------------

        if (ltrace) write(*,*) 'calc_exp ... N=',Nx,Ny,Nz

        is_ex=0

        do i=1,Nx
          do j=1,Ny
            zeros(i,j)=0
          end do
        end do

        do i=is0,ie0
          do j=js0,je0
           theta(i,j)=1.0d0
          end do
        end do

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        is=max(2,is0)
        ie=min(Nx-1,ie0)
        js=max(2,js0)
        je=min(Ny-1,je0)
!        if (x(1) .lt.0+dx/2 .and.is.eq.2     ) is=3     !REMOVE
!        if (x(Nx).gt.1-dx/2 .and.ie.eq.(Nx-1)) ie=Nx-2 
!        if (y(1) .lt.0+dy/2 .and.js.eq.2     ) js=3 
!        if (y(Ny).gt.1-dy/2 .and.je.eq.(Ny-1)) je=Ny-2

        do i=is,ie
          do j=js,je
            if (ltrace) write(*,*) 'i,j:',i,j

            any_ex=.false.

            if (.not.(do_ex.eq.0.or.
     &          can_calc_ex(chr,i,j,Nx,Ny,Nz,ex))) then
                write(*,*) ' can_calc_ex: i,j has excised neighbors,'
                write(*,*) '              so cannot be updated '
                write(*,*) ' i,j=',i,j
                write(*,*) ' x(i),y(j)=',x(i),y(j)
                write(*,*) ' chr(i,j)=',chr(i,j)
                write(*,*) ' Nx,Ny=',Nx,Ny
                write(*,*) ' dx,dy=',dx,dy
               any_ex=.true.
               is_ex=1
            else
               x0=x(i)
               y0=y(j)              
             
!               ! computes tensors at point i,j 
!               call tensor_init(
!     &                 gb_tt_np1,gb_tt_n,gb_tt_nm1,
!     &                 gb_tx_np1,gb_tx_n,gb_tx_nm1,
!     &                 gb_ty_np1,gb_ty_n,gb_ty_nm1,
!     &                 gb_xx_np1,gb_xx_n,gb_xx_nm1,
!     &                 gb_xy_np1,gb_xy_n,gb_xy_nm1,
!     &                 gb_yy_np1,gb_yy_n,gb_yy_nm1,
!     &                 psi_np1,psi_n,psi_nm1,
!     &                 zeros,zeros,zeros,
!     &                 zeros,zeros,zeros,
!     &                 zeros,zeros,zeros,
!     &                 zeros,zeros,zeros,
!     &                 g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
!     &                 gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
!     &                 h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
!     &                 A_l,A_l_x,Hads_l,
!     &                 gamma_ull,gamma_ull_x,
!     &                 riemann_ulll,ricci_ll,ricci_lu,ricci,
!     &                 einstein_ll,set_ll,
!     &                 phi10_x,phi10_xx,
!     &                 x,y,dt,chr,L,ex,Nx,Ny,i,j)

            end if

            if (.not.any_ex) then

              ! computes needed derivatives of f=r-R(chi,phi)
              call df2_int(f,f,f,
     &             tmp1, f_x, f_y, 
     &             tmp2,tmp3,tmp4,
     &             f_xx,f_xy,f_yy,
     &             x,y,dt,i,j,chr,ex,Nx,Ny,'f')
 
              ! define unit time-like vector n, normal to t=const
              ! surfaces
              n_l(1)=-1/sqrt(-g0_uu(1,1))
              do a=1,5
                n_u(a)=n_l(1)*g0_uu(a,1)+
     &                 n_l(2)*g0_uu(a,2)+
     &                 n_l(3)*g0_uu(a,3)+
     &                 n_l(4)*g0_uu(a,4)+
     &                 n_l(5)*g0_uu(a,5)
              end do
              do b=1,5
                n_l_x(1,b)=-1/2.0d0/sqrt(-g0_uu(1,1))**3*g0_uu_x(1,1,b)
              end do
              do a=1,5
                do b=1,5
                  n_u_x(a,b)=n_l_x(1,b)*g0_uu(a,1)+
     &                       n_l_x(2,b)*g0_uu(a,2)+
     &                       n_l_x(3,b)*g0_uu(a,3)+
     &                       n_l_x(4,b)*g0_uu(a,4)+
     &                       n_l_x(5,b)*g0_uu(a,5)+
     &                       n_l(1)*g0_uu_x(a,1,b)+
     &                       n_l(2)*g0_uu_x(a,2,b)+
     &                       n_l(3)*g0_uu_x(a,3,b)+
     &                       n_l(4)*g0_uu_x(a,4,b)+
     &                       n_l(5)*g0_uu_x(a,5,b)
                end do
              end do

              ! define gradients of the flow field f=r-AH_R(chi,phi)
              f0_x(1)=0
              f0_x(2)=f_x
              f0_x(3)=f_y
              f0_x(4)=0
              f0_x(5)=0
              f0_xx(1,1)=0
              f0_xx(1,2)=0
              f0_xx(1,3)=0
              f0_xx(1,4)=0
              f0_xx(1,5)=0
              f0_xx(2,2)=f_xx
              f0_xx(2,3)=f_xy
              f0_xx(2,4)=0
              f0_xx(2,5)=0
              f0_xx(3,3)=f_yy
              f0_xx(3,4)=0
              f0_xx(3,5)=0
              f0_xx(4,4)=0
              f0_xx(4,5)=0
              f0_xx(5,5)=0

              do a=1,4
                do b=a+1,5
                  f0_xx(b,a)=f0_xx(a,b)
                end do
              end do

              ! define metric on codimension-1 surfaces 
              do a=1,5
                do b=1,5
                  gam_uu(a,b)=g0_uu(a,b)+n_u(a)*n_u(b)
                end do
              end do
              do a=1,5
                do b=1,5
                  do c=1,5
                    gam_uu_x(a,b,c)=g0_uu_x(a,b,c)
     &                             +n_u_x(a,c)*n_u(b)
     &                             +n_u(a)*n_u_x(b,c)
                  end do
                end do
              end do

              ! define unit space-like vector s, orthogonal to n and
              ! projected gradient of the flow field f
              do a=1,5
                s_u(a)=0.0d0
                do b=1,5
                  s_u(a)=s_u(a)+gam_uu(a,b)*f0_x(b)
                end do
              end do
              normsusq=0.0d0
              do a=1,5
                do b=1,5
                  normsusq=normsusq+gam_uu(a,b)*f0_x(a)*f0_x(b)
                end do
              end do
              do a=1,5
                s_u(a)=s_u(a)/sqrt(normsusq)
              end do
              do a=1,5
                s_l(a)=s_u(1)*g0_ll(a,1)+
     &                 s_u(2)*g0_ll(a,2)+
     &                 s_u(3)*g0_ll(a,3)+
     &                 s_u(4)*g0_ll(a,4)+
     &                 s_u(5)*g0_ll(a,5)
              end do

              nufx=0
              do a=1,5
                nufx=nufx
     &              +n_u(a)*f0_x(a)
                nuxfx(a)=0
                gamxfxfx(a)=0
                do c=1,5
                  nuxfx(a)=nuxfx(a)
     &                    +n_u_x(c,a)*f0_x(c)
                  do d=1,5
                    gamxfxfx(a)=gamxfxfx(a)
     &                        +gam_uu_x(c,d,a)*f0_x(c)*f0_x(d)
     &                        +gam_uu(c,d)*f0_xx(c,a)*f0_x(d)
     &                        +gam_uu(c,d)*f0_x(c)*f0_xx(d,a)
                  end do
                end do
              end do
              do a=1,5
                do b=1,5
                  s_l_x(a,b)=
     &                     (f0_xx(a,b)+n_l_x(a,b)*nufx+n_l(a)*nuxfx(b))
     &                     /sqrt(normsusq)
     &                    -(f0_x(a)+n_l(a)*nufx)*gamxfxfx(b)
     &                     /2.0d0/sqrt(normsusq)**3
                end do
              end do

              ! define metric on codimension-2 surfaces 
              do a=1,5
                do b=1,5
                  sig_uu(a,b)=g0_uu(a,b)+n_u(a)*n_u(b)-s_u(a)*s_u(b)
                end do
              end do

              ! for theta: outward null expansion 
              theta(i,j)=0.0d0
              do c=1,5
                do d=1,5
                  theta(i,j)=theta(i,j)
     &                   +sig_uu(c,d)*(n_l_x(c,d)+s_l_x(c,d))
                  do e=1,5
                    theta(i,j)=theta(i,j)
     &                   -sig_uu(c,d)*gamma_ull(e,c,d)*(n_l(e)+s_l(e))
                  end do
                end do
              end do

!              theta_ads=3.0d0/2.0d0/sqrt(x0**2+y0**2)*(1+x0**2+y0**2)
!              theta(i,j)=theta(i,j)/theta_ads

              ! if theta is nan, then set it to the following instead
              if (.not.(theta(i,j).le.abs(theta(i,j)))) then
                 theta(i,j)=0.0d0
              end if
 
            end if

          end do
        end do

        !--------------------------------------------------------------
        ! fill in the borders if desired via first order extrapolation
        !--------------------------------------------------------------
 
        do j=js,je
          if (is0.eq.1) theta(1,j)=theta(2,j)
          if (ie0.eq.Nx) theta(Nx,j)=theta(Nx-1,j)
        end do

        do i=is-1,ie+1
          if (js0.eq.1) theta(i,1)=theta(i,2)
          if (je0.eq.Ny) theta(i,Ny)=theta(i,Ny-1)
        end do

        return
        end

c-----------------------------------------------------------------------
c given AH_R,AH_xc, computes the corresponding coordinate center
c and principle axis radii for excision ex_r0,ex_xc0
c-----------------------------------------------------------------------
        subroutine fill_ex_params(AH_R,AH_xc,ex_r0,ex_xc0,
     &                            AH_Nchi,AH_Nphi,dx,dy,axisym)
        implicit none
        integer axisym
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(2),ex_r0(2),ex_xc0(2)
        
        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer i,j
        real*8 x,y,AH_chi,AH_phi,dahchi,dahphi
        real*8 dx,dy
        real*8 xmin,xmax,ymin,ymax

        ! initialize fixed-size variables
        data i,j/0,0/

        data x,y,AH_chi,AH_phi/0.0,0.0,0.0,0.0/
        data dahchi,dahphi/0.0,0.0/
        data xmin,xmax,ymin,ymax/0.0,0.0,0.0,0.0/

        !--------------------------------------------------------------

        xmin=1
        xmax=0
        ymin=1
        ymax=0

        if (AH_xc(2).lt.dy) then
          dahchi=PI/(AH_Nchi-1)
        else
          dahchi=2*PI/(AH_Nchi-1)
        end if

        do i=1,AH_Nchi
          do j=1,AH_Nphi

            AH_chi=(i-1)*dahchi
            AH_phi=(j-1)*dahphi

            ! AH_R,AH_chi: polar coordinates of point on AH, wrt center of AH
            ! AH_xc(1),AH_xc(2): cartesian coordinates of center of AH, wrt origin
            ! x,y cartesian coordinates of point on AH, wrt origin
            x=AH_R(i,j)*cos(AH_chi)+AH_xc(1)
            y=AH_R(i,j)*sin(AH_chi)+AH_xc(2)

            xmin=min(x,xmin)
            ymin=min(y,ymin)

            xmax=max(x,xmax)
            ymax=max(y,ymax)

          end do
        end do

        !sets excision coordinate center as AH coordinate center
        ex_xc0(1)=AH_xc(1)
        ex_xc0(2)=AH_xc(2)

        !sets excision principal axis radii as ellipse semi-axes
        ex_r0(1)=(xmax-xmin)/2
        if (AH_xc(2).lt.dy) then
          ex_r0(2)=ymax
        else
          ex_r0(2)=(ymax-ymin)/2
        end if

        return
        end

c-----------------------------------------------------------------------
c the follow checks whether the coordinate ellipse a, 
c given by [r_a,xc_a], is entirely contained in b, given by [r_b,xc_b]
c If so, a_in_b and a_int_b is set to 1 (else 0)
c If a intersects b (but does *not* entirely contain b), 
c then a_int_b is set to 1 (else 0)
c
c NOTE: only the 2*dim "corner" points along the principle axis
c are checked, and so there are situations where a_int_b might
c incorrectly be cleared. However, for the kinds of AH's we're
c dealing with, this should not happen
c-----------------------------------------------------------------------
        subroutine is_inside(a_in_b,a_int_b,r_a,xc_a,r_b,xc_b,dim)
        implicit none
        integer a_in_b,dim,a_int_b,side
        real*8 r_a(dim),xc_a(dim)
        real*8 r_b(dim),xc_b(dim)

        real*8 d2,d0
        logical ltrace

        ! initialize fixed-size variables
        data d2,d0/0.0,0.0/

        !--------------------------------------------------------------
        ! we only check points along the principle axis
        !--------------------------------------------------------------

        a_in_b=1
        a_int_b=0

        d0=(xc_a(2) - xc_b(2))**2/r_b(2)**2
        if (dim.gt.2) d0=d0+(xc_a(3) - xc_b(3))**2/r_b(3)**2
        d2=(xc_a(1)+r_a(1) - xc_b(1))**2/r_b(1)**2 + d0
        if (d2.gt.1) a_in_b=0
        if (d2.lt.1) a_int_b=1
        d2=(xc_a(1)-r_a(1) - xc_b(1))**2/r_b(1)**2 + d0
        if (d2.gt.1) a_in_b=0
        if (d2.lt.1) a_int_b=1

        d0=(xc_a(1) - xc_b(1))**2/r_b(1)**2
        if (dim.gt.2) d0=d0+(xc_a(3) - xc_b(3))**2/r_b(3)**2
        d2=(xc_a(2)+r_a(2) - xc_b(2))**2/r_b(2)**2 + d0
        if (d2.gt.1) a_in_b=0
        if (d2.lt.1) a_int_b=1
        d2=(xc_a(2)-r_a(2) - xc_b(2))**2/r_b(2)**2 + d0
        if (d2.gt.1) a_in_b=0
        if (d2.lt.1) a_int_b=1

        if (dim.gt.2) then
           d0=(xc_a(1) - xc_b(1))**2/r_b(1)**2
           d0=d0+(xc_a(2) - xc_b(2))**2/r_b(2)**2
           d2=(xc_a(3)+r_a(3) - xc_b(3))**2/r_b(3)**2 + d0
           if (d2.gt.1) a_in_b=0
           if (d2.lt.1) a_int_b=1
           d2=(xc_a(3)-r_a(3) - xc_b(3))**2/r_b(3)**2 + d0
           if (d2.gt.1) a_in_b=0
           if (d2.lt.1) a_int_b=1
        end if

        !--------------------------------------------------------------
        ! the above would have missed an intersection where a
        ! is large, and b is small and away from a's principle axis.
        ! to check for an intersection, all of b's points must
        ! be one side or another of a:
        !--------------------------------------------------------------

        side=0

        d0=(xc_b(2) - xc_a(2))**2/r_a(2)**2
        if (dim.gt.2) d0=d0+(xc_b(3) - xc_a(3))**2/r_a(3)**2
        d2=(xc_b(1)+r_b(1) - xc_a(1))**2/r_a(1)**2 + d0
        if (d2.gt.1) side=side+1
        if (d2.lt.1) side=side-1
        d2=(xc_b(1)-r_b(1) - xc_a(1))**2/r_a(1)**2 + d0
        if (d2.gt.1) side=side+1
        if (d2.lt.1) side=side-1

        d0=(xc_b(1) - xc_a(1))**2/r_a(1)**2
        if (dim.gt.2) d0=d0+(xc_b(3) - xc_a(3))**2/r_a(3)**2
        d2=(xc_b(2)+r_b(2) - xc_a(2))**2/r_a(2)**2 + d0
        if (d2.gt.1) side=side+1
        if (d2.lt.1) side=side-1
        d2=(xc_b(2)-r_b(2) - xc_a(2))**2/r_a(2)**2 + d0
        if (d2.gt.1) side=side+1
        if (d2.lt.1) side=side-1

        if (dim.gt.2) then
           d0=(xc_b(1) - xc_a(1))**2/r_a(1)**2
           d0=d0+(xc_b(2) - xc_a(2))**2/r_a(2)**2
           d2=(xc_b(3)+r_b(3) - xc_a(3))**2/r_a(3)**2 + d0
           if (d2.gt.1) side=side+1
           if (d2.lt.1) side=side-1
           d2=(xc_b(3)-r_b(3) - xc_a(3))**2/r_a(3)**2 + d0
           if (d2.gt.1) side=side+1
           if (d2.lt.1) side=side-1
        end if

        if (abs(side).ne.2*dim)  a_int_b=1

        return
        end

c-----------------------------------------------------------------------
c is a point on the hypersurface r=R(chi) (relative to the center
c AH_xc) interior to the cartesian bounding box?? 
c-----------------------------------------------------------------------
        subroutine ah_is_int(is_int,AH_R,AH_xc,i,j,bbox,dx,dy,
     &                       AH_Nchi,AH_Nphi,axisym)
        implicit none
        integer axisym
        integer AH_Nchi,AH_Nphi,i,j,is_int
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(2)
        real*8 dx,dy,bbox(4)
        
        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 x,y,AH_chi,AH_phi,dahchi,dahphi

!LOOKING!
        logical ltrace
        parameter (ltrace=.false.)
!        parameter (ltrace=.true.)

        ! initialize fixed-size variables 
        data x,y,AH_chi,AH_phi/0.0,0.0,0.0,0.0/
        data dahchi,dahphi/0.0,0.0/

        !--------------------------------------------------------------

        if (AH_xc(2).lt.dy) then
          dahchi=PI/(AH_Nchi-1)
        else
          dahchi=2*PI/(AH_Nchi-1)
        end if

        AH_chi=(i-1)*dahchi
        AH_phi=(j-1)*dahphi

        ! AH_R,AH_chi,AH_phi: polar coordinates of point on AH, wrt center of AH
        ! AH_xc(1),AH_xc(2): cartesian coordinates of center of AH, wrt origin
        ! x,y,z cartesian coordinates of point on AH, wrt origin
        x=AH_R(i,j)*cos(AH_chi)+AH_xc(1)
        y=AH_R(i,j)*sin(AH_chi)+AH_xc(2)

        if (AH_xc(2).lt.dy) then
          if ((x-bbox(1)).gt.(2.5*dx).and.(bbox(2)-x).gt.(2.5*dx).and.
     &        ((y-bbox(3)).gt.(2.5*dy).or.bbox(3).lt.(0+dy/2)).and. !allow y=0
     &        (bbox(4)-y).gt.(2.5*dy))
     &    then
             is_int=1
          else
             is_int=0
          end if
        else 
          if ((x-bbox(1)).gt.(2.5*dx).and.(bbox(2)-x).gt.(2.5*dx).and.
     &        (y-bbox(3)).gt.(2.5*dy).and.(bbox(4)-y).gt.(2.5*dy))
     &    then
             is_int=1
          else
             is_int=0
          end if
        end if

        if (ltrace) then
           write(*,*)
           write(*,*) 'is_int: ',is_int
           write(*,*) 'bbox(1),bbox(2)=',bbox(1),bbox(2)
           write(*,*) 'bbox(3),bbox(4)=',bbox(3),bbox(4)
           write(*,*) 'x,y=',x,y
           write(*,*) 'R,AH_chi',AH_R(i,j),AH_chi
           write(*,*) 'xc,yc',AH_xc(1),AH_xc(2)
        end if

        return
        end

c-----------------------------------------------------------------------
c marks unowned points (-1) as owned by node=rank if expansion can be
c calculated via interior stencils (2 *cells* away from boundaries,
c as calculated above) on node rank ... fills in AH_lev(..)=L as well,
c note that physical boundaries are exceptions; see ah_is_int() 
c-----------------------------------------------------------------------
        subroutine ah_fill_own(AH_R,AH_xc,AH_own,AH_lev,bbox,dx,dy,
     &                         rank,L,AH_Nchi,AH_Nphi,axisym)
        implicit none
        integer axisym
        integer rank,AH_Nchi,AH_Nphi,L
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(2)
        integer AH_own(AH_Nchi,AH_Nphi),AH_lev(AH_Nchi,AH_Nphi)
        real*8 dx,dy,bbox(4)
        
        integer i,j,is_int

        logical ltrace
        parameter (ltrace=.false.)

        ! initialize fixed-size variables
        data i,j,is_int/0,0,0/

        !--------------------------------------------------------------

        do i=1,AH_Nchi
           do j=1,AH_Nphi
              if (AH_own(i,j).eq.-1) then
                 call ah_is_int(is_int,AH_R,AH_xc,i,j,bbox,dx,dy,
     &                          AH_Nchi,AH_Nphi,axisym)
                 if (is_int.eq.1) then
                    AH_own(i,j)=rank
                    AH_lev(i,j)=L
                 end if
                 if (ltrace) then
                    write(*,*) 'ah_fill_own: bbox=',bbox
                    write(*,*) 'i,j=',i,j,' own=',AH_own(i,j)
                 end if
              end if
           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c evaluates the hypersurface function f = r - R(chi,phi) in
c the range [is..ie,js..je,ks..ke], with
c
c x=AH_R(chi,phi) + AH_xc(1)
c y=chi/PI + AH_xc(2) 
c
c AH_chi goes from 0 ("x" axis) to PI ("x" axis)
c AH_phi goes from 0 ("z" axis) to PI ("z" axis) 
c
c AH_chi and AH_phi are measured relative to center of AH, which in
c cartesian coordinates is given by AH_xc
c-----------------------------------------------------------------------
        subroutine ah_fill_f(AH_R,AH_xc,f,is,ie,js,je,x,y,z,
     &                       AH_Nchi,AH_Nphi,Nx,Ny,Nz,axisym)
        implicit none
        integer axisym
        integer is,ie,js,je,AH_Nchi,AH_Nphi,Nx,Ny,Nz
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(2)
        real*8 x(Nx),y(Ny),z(Nz),f(Nx,Ny)
        
        integer i,j,i0,j0,is0,ie0,js0,je0
        real*8 AH_chi,AH_phi,r,xb,yb,zb,dahchi,dahphi,ft,fp,rtp
        real*8 xb0,yb0,zb0

        real*8 dx,dy

        real*8 PI
        parameter (PI=3.141592653589793d0)

        logical ltrace
        parameter (ltrace=.false.)

        ! initialize fixed-size variables
        data i,j,i0,j0,is0,ie0/0,0,0,0,0,0/
        data js0,je0/0,0/

        data AH_chi,AH_phi,r,xb,yb,zb/0.0,0.0,0.0,0.0,0.0,0.0/
        data dahchi,dahphi,ft,fp,rtp/0.0,0.0,0.0,0.0,0.0/
        data xb0,yb0,zb0/0.0,0.0,0.0/

        data dx,dy/0.0,0.0/ 
        !--------------------------------------------------------------

        if (ltrace) then
           write(*,*) 'ah_fill_f:'
           write(*,*) 'is,ie,js,je:',is,ie,js,je
        end if

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (AH_xc(2).lt.dy) then
          dahchi=PI/(AH_Nchi-1)
        else
          dahchi=2*PI/(AH_Nchi-1)
        end if

        xb0=AH_xc(1)
        yb0=AH_xc(2)

        is0=is
        js0=js
        ie0=ie
        je0=je

!        if (is0.eq.1.and.abs(x(1)+1).lt.1.0d-10) is0=2      !REMOVE
!        if (js0.eq.1.and.abs(y(1)+1).lt.1.0d-10) js0=2
!        if (ie0.eq.Nx.and.abs(x(Nx)-1).lt.1.0d-10) ie0=Nx-1
!        if (je0.eq.Ny.and.abs(y(Ny)-1).lt.1.0d-10) je0=Ny-1

        do i=is0,ie0
           xb=x(i)
           do j=js0,je0
              yb=y(j)

                 !xb,yb: cartesian coordinates of some point "not on, but near" AH, wrt origin
                 !AH_xc(1),AH_xc(2): cartesian coordinates of center of AH, wrt origin
                 !r,AH_chi,AH_phi polar coordinates of the point "not on, but near" AH, wrt center of AH
                 !AH_R,AH_chi,AH_phi: polar coordinates of point "projected" on AH, wrt center of AH)
                 r=sqrt( (xb-AH_xc(1))**2
     &                  +(yb-AH_xc(2))**2 )
                 AH_chi=atan2( yb-AH_xc(2),xb-AH_xc(1) )
                 AH_phi=0
                 if (AH_chi.lt.0) AH_chi=AH_chi+2*PI

                 !-----------------------------------------------------
                 ! bi-linear interpolate in AH_chi and AH_phi
                 !-----------------------------------------------------
                 i0=min(int(AH_chi/dahchi+1),AH_Nchi-1)
                 if (axisym.eq.0) then
                    j0=min(int(AH_phi/dahphi+1),AH_Nphi-1)
                 else
                    j0=1
                 end if

                 if (r.eq.0) then
!                   write(*,*) 'WARNING!!! i0 being set to 1'
!                   write(*,*) 'center of AH is sampled by ah_fill_f()'
!                   write(*,*) '(AH points too close to center of AH)'
!                   write(*,*) '(... increase AH_r0 or resolution)'
!                    write(*,*) ' i0,j0=',i0,j0
!                    write(*,*) 'xb,yb=',xb,yb
!                    write(*,*) 'is,ie=',is,ie
!                    write(*,*) 'js,je=',js,je
!                    write(*,*) 'i,j=',i,j
!                    write(*,*) 'AH_R(i0,j0)',AH_R(i0,j0)
!                    write(*,*) 'AH_chi,AH_phi=',AH_chi,AH_phi
!                    write(*,*) 'r=',r
!                    write(*,*) 'rtp=',rtp
!                    write(*,*) '---------------------------'
                   i0=1
                 end if

                 if (ltrace) write(*,*) 'i0,j0=',i0,j0, ' Nchi,Nphi=',
     &              AH_Nchi,AH_Nphi

                 ft=((AH_chi-(i0-1)*dahchi))/dahchi
                 fp=((AH_phi-(j0-1)*dahphi))/dahphi

                 if (axisym.eq.0) then
                    rtp=AH_R(i0,j0)*(1-ft)*(1-fp)+
     &               AH_R(i0+1,j0)*(ft)*(1-fp)+
     &               AH_R(i0,j0+1)*(1-ft)*(fp)+
     &               AH_R(i0+1,j0+1)*(ft)*(fp)
                 else
                    rtp=AH_R(i0,j0)*(1-ft)+
     &                  AH_R(i0+1,j0)*(ft)
                 end if

                 f(i,j)=r-rtp

!                 !TEST in ah_fill_f()!
!                 if (rtp.lt.0.29d0.or.rtp.gt.0.31d0) then
!                   write(*,*) 'i0,j0=',i0,j0
!                   write(*,*) 'AH_R(i0,j0)=',AH_R(i0,j0)
!                   write(*,*) 'AH_R(i0+1,j0)=',AH_R(i0+1,j0)
!                   write(*,*) 'AH_R(i0,j0+1)=',AH_R(i0,j0+1)
!                   write(*,*) 'AH_R(i0+1,j0+1)=',AH_R(i0+1,j0+1)
!                   write(*,*) '-------------------'
!                 end if
!                 do i=1,AH_Nchi
!                   do j=1,AH_Nphi
!                     if (AH_R(i,j).lt.0.89d0.or.AH_R(i,j).gt.0.91d0) then
!                       write(*,*) 'AH_R,i,j=',AH_R(i,j),i,j
!                     end if
!                   end do
!                 end do

                 if (ltrace) then
                    write(*,*) 'xb,yb=',xb,yb
                    write(*,*) 'AH_R(i0,j0)',AH_R(i0,j0)
                    write(*,*) 'AH_chi,AH_phi=',AH_chi,AH_phi
                    write(*,*) 'r=',r
                    write(*,*) 'rtp=',rtp
                    write(*,*) '---------------------------'
                 end if
           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c the following calculates the expansion at the single point (i0,j0)
c on the hypersurface f=r-AH_R(chi,phi) (with offset AH_xc)
c
c also returns the area element da at (i0,j0);
c and if (i0,j0) corresponds to chi=Pi/2, the equatorial circumference
c element d_ceq, and if (i0,j0) corresponds to phi=0, the polar
c circumference element d_cp 
c
c NOTE: if Nchi is even, will miss equatorial plane by dahchi
c
c area element is set to -1000 if point was too close to excision boundary
c to compute.
c
c NOTE: this routine modifies f and chi, but only in the vicinity of (i0,j0)
c
c if use_AH_new!=0 (in AH_new.inc), then computes chi using new tensor routines
c-----------------------------------------------------------------------
        subroutine calc_exp0(AH_R,AH_xc,AH_theta,i0,j0,AH_Nchi,
     &                    AH_Nphi,theta,f,da,d_ceq,d_cp,
     &                    gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                    gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                    gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                    gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                    gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                    gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                    psi_np1,psi_n,psi_nm1,
     &                    L,x,y,z,dt,chr,ex,do_ex,
     &                    Nx,Ny,Nz,axisym)
        implicit none
        integer axisym
        integer Nx,Ny,Nz,i0,j0,AH_Nchi,AH_Nphi,do_ex
        real*8 theta(Nx,Ny),f(Nx,Ny),AH_xc(2),da,d_ceq,d_cp
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_theta(AH_Nchi,AH_Nphi)
        real*8 chr(Nx,Ny),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L

        real*8 gb_tt_np1(Nx,Ny),gb_tt_n(Nx,Ny),gb_tt_nm1(Nx,Ny)
        real*8 gb_tx_np1(Nx,Ny),gb_tx_n(Nx,Ny),gb_tx_nm1(Nx,Ny)
        real*8 gb_ty_np1(Nx,Ny),gb_ty_n(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xx_n(Nx,Ny),gb_xx_nm1(Nx,Ny)
        real*8 gb_xy_np1(Nx,Ny),gb_xy_n(Nx,Ny),gb_xy_nm1(Nx,Ny)
        real*8 gb_yy_np1(Nx,Ny),gb_yy_n(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 psi_np1(Nx,Ny),psi_n(Nx,Ny),psi_nm1(Nx,Ny)

        real*8 cosx(Nx),cosy(Ny),cosz(Nz)
        real*8 sinx(Nx),siny(Ny),sinz(Nz)

        integer i,j,i1,j1,is,ie,js,je,is_ex

        real*8 x0,y0
        real*8 PI
        parameter (PI=3.141592653589793d0)
        real*8 dx,dy,dahchi,dahphi,AH_chi,AH_phi,fx,fy,fz
        real*8 r_t,r_p,st,sp,ct,cp,r
        real*8 dxb_dt,dyb_dt,dzb_dt
        real*8 dxb_dp,dyb_dp,dzb_dp
        real*8 dtt,dpp,dpt
        real*8 gb_xx0,gb_xy0,gb_yy0,psi0
        real*8 f_x,f_y,f_xx,f_xy,f_yy

        real*8 drh_dch,dx_dch,dy_dch
        real*8 dx_dch_rh,dy_dch_rh
        real*8 g_xx0,g_xy0,g_yy0,g_thth0,g_phph0
        real*8 g_chch0,g_const0

        real*8 rho0,f0

        real*8 g0_xx_ads0,g0_xy_ads0,g0_yy_ads0,g0_psi_ads0

        integer is_bad

!LOOKING!
        logical ltrace
        parameter (ltrace=.false.)
!        parameter (ltrace=.true.)

        ! initialize fixed-size variables 
        data i,j,i1,j1/0.0,0.0,0.0,0.0/
        data is,ie,js,je/0.0,0.0,0.0,0.0/
        data is_ex/0.0/

        data x0,y0/0.0,0.0/
        data dx,dy/0.0,0.0/
        data dahchi,dahphi/0.0,0.0/
        data AH_chi,AH_phi,fx,fy,fz/0.0,0.0,0.0,0.0,0.0/
        data r_t,r_p,st,sp,ct,cp,r/0.0,0.0,0.0,0.0,0.0,0.0,0.0/
        data dxb_dt,dyb_dt,dzb_dt/0.0,0.0,0.0/
        data dxb_dp,dyb_dp,dzb_dp/0.0,0.0,0.0/
        data dtt,dpp,dpt/0.0,0.0,0.0/
        data gb_xx0,gb_xy0,gb_yy0/0.0,0.0,0.0/
        data psi0/0.0/
        data f_x,f_y/0.0,0.0/
        data f_xx,f_xy,f_yy/0.0,0.0,0.0/

        data drh_dch,dx_dch,dy_dch/0.0,0.0,0.0/
        data dx_dch_rh,dy_dch_rh/0.0,0.0/
        data g_xx0,g_xy0,g_yy0,g_thth0,g_phph0/0.0,0.0,0.0,0.0,0.0/
        data g_chch0,g_const0/0.0,0.0/

        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (AH_xc(2).lt.dy) then
          dahchi=PI/(AH_Nchi-1)
        else
          dahchi=2*PI/(AH_Nchi-1)
        end if
        dahphi=1.0d0

        AH_chi=(i0-1)*dahchi
        AH_phi=(j0-1)*dahphi

        !AH_R,AH_chi,AH_phi polar coordinates of point on AH, wrt center of AH
        !AH_xc(1),AH_xc(2) cartesian coordinates of center of AH, wrt origin
        !x0,y0 cartesian coordinates of point on AH, wrt origin
        x0=AH_R(i0,j0)*cos(AH_chi)+AH_xc(1)
        y0=AH_R(i0,j0)*sin(AH_chi)+AH_xc(2)

        da=0
        d_ceq=0
        d_cp=0

        ! extract (i,j) cartesian grid point closest to x0,y0
        i=(x0-x(1))/dx+1
        j=(y0-y(1))/dy+1

        ! for bilinear interpolation of theta from 
        ! (i,j),(i+1,j),(i,j+1),(i+1,j+1)
        !(NOTE: since, x0,y0 does not necessarily lie on a grid point these fx,fy are *not* identically zero)
        fx=(x0-((i-1)*dx+x(1)))/dx
        fy=(y0-((j-1)*dy+y(1)))/dy

        if (ltrace) then
           write(*,*) 'calc_exp0: i0,j0=',i0,j0
           write(*,*) ' AH_R,AH_chi,AH_phi=',AH_R(i0,j0),AH_chi,AH_phi
           write(*,*) ' i,j=',i,j
           write(*,*) ' x0,y0=',x0,y0
           write(*,*) ' x(1),y(1)=',x(1),y(1)
           write(*,*) ' dx,dy=',dx,dy
           write(*,*) '----------------------------'
        end if

        !------------------------------------------------------------
        ! we will use bilinear interpolation to calculate theta,
        ! hence we need f in a sufficient box around to caculate theta
        ! in a 2^2 box adjacent to (i,j)
        !------------------------------------------------------------

        is_bad=0
        if (i.lt.2.or.j.lt.2.or.i.ge.Nx.or.j.ge.Ny)
     &  then
          is_bad=1
          if (i0.eq.1.or.i0.eq.AH_Nchi) 
     &    then !these are allowed ... reg_ah_r() will fill it in 
            AH_theta(i0,j0)=0
            return
          end if              
        end if

        if (is_bad.eq.1) then
!           write(*,*) '--------------------------------------------'
!           write(*,*) 'WARNING!!! theta set to 100 at i0,j0=',i0,j0
!           write(*,*) 'and regularity will *not* take care of it!'
!           write(*,*) '(AH points too bunched up near i0 or j0 endpts)'
!           write(*,*) '(... reduce AH_Nchi or AH_Nphi)'
!           write(*,*) ' i0,j0=',i0,j0
!           write(*,*) ' AH_R,AH_chi,AH_phi=',AH_R(i0,j0),AH_chi,AH_phi
!           write(*,*) ' xc(1),xc(2)=',AH_xc(1),AH_xc(2)
!           write(*,*) ' i,j=',i,j
!           write(*,*) ' x0,y0=',x0,y0
!           write(*,*) ' x(1),y(1)=',x(1),y(1)
!           write(*,*) ' dx,dy=',dx,dy
!           write(*,*) ' dahchi,dahphi=',dahchi,dahphi
!           write(*,*) '--------------------------------------------'
           AH_theta(i0,j0)=100
        end if

        is=max(1,i-2)
        ie=min(Nx,i+3)
        js=max(1,j-2)
        je=min(Ny,j+3)

        call ah_fill_f(AH_R,AH_xc,f,is,ie,js,je,x,y,z,
     &              AH_Nchi,AH_Nphi,Nx,Ny,Nz,axisym)

        is=i
        ie=i+1
        js=j
        je=j+1

        if (do_ex.ne.0) then
           do i1=is,ie   
              do j1=js,je   
                 if (chr(i1,j1).eq.ex) then
                    write(*,*) ' calc_exp0: pt i1,j1 is excised'
                    write(*,*) ' AH_chi,AH_phi=',AH_chi,AH_phi
                    write(*,*) ' x0,y0=',x0,y0
                    write(*,*) ' dx,dy=',dx,dy
                    write(*,*) ' i,j=',i,j
                    write(*,*) ' i1,j1=',i1,j1
                    write(*,*) ' x(i1),y(j1)=',x(i1),y(j1)
                    write(*,*) ' chr(i1,j1)=',chr(i1,j1)
                    write(*,*) ' Nx,Ny=',Nx,Ny
                    da=-10000
                    AH_theta(i0,j0)=0
                    return
                 end if
              end do
           end do
        end if

        call calc_exp(theta,f,is,ie,js,je,is_ex,
     &             gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &             gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &             gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &             gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &             gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &             gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &             psi_np1,psi_n,psi_nm1,
     &             L,x,y,z,dt,chr,ex,do_ex,Nx,Ny,Nz)

        ! interpolate theta from 
        ! (i,j),(i+1,j),(i,j+1),(i+1,j+1)
        AH_theta(i0,j0) = (1-fx)*(1-fy)*theta(i,j)+
     &                 (  fx)*(1-fy)*theta(i+1,j)+
     &                 (1-fx)*(  fy)*theta(i,j+1)+
     &                 (  fx)*(  fy)*theta(i+1,j+1)

        !--------------------------------------------------------------
        ! proper area element
        ! 
        ! NOTE: some definitions below are copied from tensor_init() in
        ! misc.f, ideally should just call tensor_init()
        !--------------------------------------------------------------

        rho0=sqrt(x0**2+y0**2)
        f0=(1-rho0**2)**2+4*rho0**2/L**2

        drh_dch=(AH_R(i0+1,j0)-AH_R(i0-1,j0))/2/dahchi

        dx_dch=cos(AH_chi)*drh_dch-rho0*sin(AH_chi)  !when rho=rho(chi)
        dy_dch=sin(AH_chi)*drh_dch+rho0*cos(AH_chi)  

        dx_dch_rh=rho0*sin(AH_chi)  !when rho=const
        dy_dch_rh=rho0*cos(AH_chi)

        g0_xx_ads0 =(x0**2*(1+rho0**2)**2/f0+y0**2)
     &             /(1-rho0**2)**2
     &             /rho0**2
     &             *4
        g0_xy_ads0 =((1+rho0**2)**2/f0-1)
     &             /(1-rho0**2)**2
     &             /rho0**2
     &             *x0*y0
     &             *4
        g0_yy_ads0 =(y0**2*(1+rho0**2)**2/f0+x0**2)
     &             /(1-rho0**2)**2
     &             /rho0**2
     &             *4
        g0_psi_ads0=(y0**2)
     &             /(1-rho0**2)**2
     &             *4

        gb_xx0=((1-fx)*(1-fy)*gb_xx_n(i,j)+
     &       (  fx)*(1-fy)*gb_xx_n(i+1,j)+
     &       (1-fx)*(  fy)*gb_xx_n(i,j+1)+
     &       (  fx)*(  fy)*gb_xx_n(i+1,j+1))

        gb_xy0=((1-fx)*(1-fy)*gb_xy_n(i,j)+
     &       (  fx)*(1-fy)*gb_xy_n(i+1,j)+
     &       (1-fx)*(  fy)*gb_xy_n(i,j+1)+
     &       (  fx)*(  fy)*gb_xy_n(i+1,j+1))

        gb_yy0=((1-fx)*(1-fy)*gb_yy_n(i,j)+
     &       (  fx)*(1-fy)*gb_yy_n(i+1,j)+
     &       (1-fx)*(  fy)*gb_yy_n(i,j+1)+
     &       (  fx)*(  fy)*gb_yy_n(i+1,j+1))

        psi0=((1-fx)*(1-fy)*psi_n(i,j)+
     &       (  fx)*(1-fy)*psi_n(i+1,j)+
     &       (1-fx)*(  fy)*psi_n(i,j+1)+
     &       (  fx)*(  fy)*psi_n(i+1,j+1))

        g_xx0   =g0_xx_ads0+gb_xx0*(1-rho0**2)
        g_xy0   =g0_xy_ads0+gb_xy0*(1-rho0**2)
        g_yy0   =g0_yy_ads0+gb_yy0*(1-rho0**2)
        g_thth0 =g0_psi_ads0+psi0*(1-rho0**2)*y0**2
        g_phph0 =g0_psi_ads0+psi0*(1-rho0**2)*y0**2

        g_chch0 =g_xx0*dx_dch_rh**2+2*g_xy0*dx_dch_rh*dy_dch_rh
     &          +g_yy0*dy_dch_rh**2
        g_const0=g_xx0*dx_dch**2+2*g_xy0*dx_dch*dy_dch+g_yy0*dy_dch**2

        da=4*PI*sqrt(g_const0*g_thth0*g_phph0)*dahchi

        if (i0.eq.((AH_Nchi-1)/2)+1) then
           d_ceq=2*PI*sqrt(g_phph0)
        end if
        d_cp=2*sqrt(g_chch0)*dahchi

        if (is_ex.eq.1) then
           da=-10000
        end if

        return
        end

c-----------------------------------------------------------------------
c Kreiss-Oliger-like smoothing of R (or theta, or whatever)
c
c + removes *isolated* "spikes" that sometimes arise in late time
c distorted holes prior to merger with rem_spikes option
c-----------------------------------------------------------------------
        subroutine smooth_ah_r(AH_R,AH_w1,AH_eps,AH_Nchi,AH_Nphi)
        implicit none
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_eps
        real*8 AH_w1(AH_Nchi,AH_Nphi)
        
        integer i,j,jp1,jp2,jm1,jm2,j0
        real*8 r0,r1,rip1,rip2,rim1,rim2,rjp1,rjp2,rjm1,rjm2,r
        real*8 avg_var,v1,v2,v3,v4

        logical rem_spikes
        parameter (rem_spikes=.false.)
        real*8 max_var
        parameter (max_var=5)

        ! initialize fixed-size variables 
        data i,j,jp1,jp2,jm1,jm2,j0/0,0,0,0,0,0,0/

        data r0,r1,rip1,rip2/0.0,0.0,0.0,0.0/
        data rim1,rim2,rjp1,rjp2/0.0,0.0,0.0,0.0/
        data rjm1,rjm2,r/0.0,0.0,0.0/
        data avg_var,v1,v2,v3,v4/0.0,0.0,0.0,0.0,0.0/
 
        !--------------------------------------------------------------

        if (AH_Nphi.gt.1.and.rem_spikes) then
           avg_var=0
           do i=1,AH_Nchi
              do j=1,AH_Nphi
                 AH_w1(i,j)=0
              end do
           end do
           do i=2,AH_Nchi-1
              do j=1,AH_Nphi-1
                 jm1=j-1
                 if (jm1.lt.1) jm1=jm1+AH_Nphi-1
                 jp1=j+1
                 if (jp1.gt.AH_Nphi) jp1=jp1-AH_Nphi+1
                 r=AH_R(i,j)
                 v1=(r-AH_R(i+1,j))**2
                 v2=(r-AH_R(i-1,j))**2
                 v3=(r-AH_R(i,jp1))**2
                 v4=(r-AH_R(i,jm1))**2
                 !-----------------------------------------------------
                 ! by subtracting the max, we are only sensitive
                 ! to isolated spikes (i.e, the points adjacent to 
                 ! the spike will now have a small variance)
                 !-----------------------------------------------------
                 AH_w1(i,j)=sqrt(v1+v2+v3+v4-max(v1,v2,v3,v4))
                 avg_var=avg_var+AH_w1(i,j)
              end do
           end do
           avg_var=avg_var/(AH_Nphi-1)/(AH_Nchi-2)
           do i=2,AH_Nchi-1
              do j=1,AH_Nphi-1
                 jm1=j-1
                 if (jm1.lt.1) jm1=jm1+AH_Nphi-1
                 jp1=j+1
                 if (jp1.gt.AH_Nphi) jp1=jp1-AH_Nphi+1
                 if (AH_w1(i,j)/avg_var.gt.max_var.and.
     &               AH_w1(i+1,j)/avg_var.lt.max_var.and.
     &               AH_w1(i-1,j)/avg_var.lt.max_var.and.
     &               AH_w1(i,jp1)/avg_var.lt.max_var.and.
     &               AH_w1(i,jm1)/avg_var.lt.max_var) then
                    AH_R(i,j)=0.25d0*(AH_R(i+1,j)+AH_R(i-1,j)+
     &                                AH_R(i,jp1)+AH_R(i,jm1))
                 end if
              end do
           end do
           do i=1,AH_Nchi
              AH_R(i,AH_Nphi)=AH_R(i,1)
           end do
        end if

        r0=0
        r1=0
        do i=1,AH_Nphi
           r0=r0+AH_R(1,i)
           r1=r1+AH_R(AH_Nchi,i)
        end do
        r0=r0/AH_Nphi
        r1=r1/AH_Nphi
        do i=1,AH_Nphi
           AH_R(1,i)=r0
           AH_R(AH_Nchi,i)=r1
        end do

        do i=1,AH_Nchi
           do j=1,AH_Nphi
              AH_w1(i,j)=AH_R(i,j)
           end do
        end do

        do i=2,AH_Nchi-1
           do j=1,AH_Nphi
              r=AH_w1(i,j)
              rip1=AH_w1(i+1,j)
              rim1=AH_w1(i-1,j)
              if (i.eq.(AH_Nchi-1)) then
                 rip2=r
              else
                 rip2=AH_w1(i+2,j)
              end if
              if (i.eq.2) then
                 rim2=r
              else
                 rim2=AH_w1(i-2,j)
              end if
              if (AH_Nphi.lt.5) then
                 jp1=j
                 jm1=j
                 jp2=j
                 jm2=j
              else
                 jp1=j+1
                 jp2=j+2
                 jm1=j-1
                 jm2=j-2
                 if (jp1.gt.AH_Nphi) jp1=jp1-AH_Nphi+1
                 if (jp2.gt.AH_Nphi) jp2=jp2-AH_Nphi+1
                 if (jm1.lt.1) jm1=jm1+AH_Nphi-1
                 if (jm2.lt.1) jm2=jm2+AH_Nphi-1
              end if
              rjp1=AH_w1(i,jp1)
              rjp2=AH_w1(i,jp2)
              rjm1=AH_w1(i,jm1)
              rjm2=AH_w1(i,jm2)

              AH_R(i,j)=r-AH_eps/16*(
     &          rip2+rim2-4*(rip1+rim1)+6*r +
     &          rjp2+rjm2-4*(rjp1+rjm1)+6*r)
           end do
        end do
        
        do i=1,AH_Nchi
           AH_R(i,1)=AH_R(i,AH_Nphi)
        end do

        return 
        end

c-----------------------------------------------------------------------
c enforce on-axis regularity of R (or theta, or whatever)
c-----------------------------------------------------------------------
        subroutine reg_ah_r(AH_R,AH_Nchi,AH_Nphi)
        implicit none
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi)
        
        integer i,j
        real*8 r0,r1

        ! initialize fixed-size variables 
        data i,j/0,0/
        data r0,r1/0.0,0.0/

        !--------------------------------------------------------------

!        ! set by average value
!        r0=0
!        r1=0
!        do j=1,AH_Nphi
!           r0=r0+AH_R(1,j)
!           r1=r1+AH_R(AH_Nchi,j)
!        end do
!        r0=r0/AH_Nphi
!        r1=r1/AH_Nphi
!        do j=1,AH_Nphi
!           AH_R(1,j)=r0
!           AH_R(AH_Nchi,j)=r1
!        end do

!        ! testing pureads
!        do i=1,AH_Nchi
!           AH_R(i,1)=3.0d0
!           AH_R(i,AH_Nphi)=3.0d0
!        end do
!        do j=2,AH_Nphi-1 !i=1, i=AH_Nchi are the poles i.e. the same pt for all j
!           AH_R(1,j)=AH_R(1,1)
!           AH_R(AH_Nchi,j)=AH_R(AH_Nchi,1)
!        end do
!
!        ! set by zero-derivative extrapolation
!        do i=1,AH_Nchi
!           AH_R(i,1)=(AH_R(i,2)*4-AH_R(i,3))/3
!           AH_R(i,AH_Nphi)=(AH_R(i,AH_Nphi-1)*4-AH_R(i,AH_Nphi-1))/3
!        end do
!        do j=1,AH_Nphi
!           AH_R(1,j)=(AH_R(2,j)*4-AH_R(3,j))/3
!           AH_R(AH_Nchi,j)=(AH_R(AH_Nchi-1,j)*4-AH_R(AH_Nchi-2,j))/3
!        end do

        ! set by zero-derivative extrapolation
        AH_R(1,1)=(AH_R(2,1)*4-AH_R(3,1))/3
        AH_R(AH_Nchi,1)=(AH_R(AH_Nchi-1,1)*4-AH_R(AH_Nchi-2,1))/3

        return 
        end

c-----------------------------------------------------------------------
c the following changes AH_xc so that median(AH_R) is zero, and
c adjusts AH_R conversely. 
c
c NOTE: the manner in which AH_R is adjusted assumes that the change 
c in AH_xc is *small* !
c-----------------------------------------------------------------------
        subroutine adjust_ah_xc(AH_R,AH_xc,AH_Nchi,AH_Nphi,
     &                          dx,dy,dz,axisym)
        implicit none
        integer axisym
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(2)
        
        integer i,j
        real*8 dahchi,dahphi,AH_chi,AH_phi
        real*8 xmin,xmax,ymin,ymax,zmin,zmax,x0,y0,z0,dx0,dy0,dz0
        real*8 dx,dy,dz
        real*8 csx,csy,csz

        real*8 PI
        parameter (PI=3.141592653589793d0)

!LOOKING!
        logical ltrace
        parameter (ltrace=.false.)
!        parameter (ltrace=.true.)

        ! initialize fixed-size variables
        data i,j/0,0/

        data dahchi,dahphi,AH_chi,AH_phi/0.0,0.0,0.0,0.0/
        data xmin,xmax,ymin,ymax/0.0,0.0,0.0,0.0/
        data zmin,zmax,x0,y0,z0/0.0,0.0,0.0,0.0,0.0/
        data dx0,dy0,dz0/0.0,0.0,0.0/
        data csx,csy,csz/0.0,0.0,0.0/

        !--------------------------------------------------------------

        xmin=0
        ymin=0
        zmin=0
        xmax=0
        ymax=0
        zmax=0

        if (AH_xc(2).lt.dy) then
          dahchi=PI/(AH_Nchi-1)
        else
          dahchi=2*PI/(AH_Nchi-1)
        end if

        do i=1,AH_Nchi
           AH_chi=(i-1)*dahchi
           do j=1,AH_Nphi
              AH_phi=(j-1)*dahphi

              !x0,y0 cartesian coordinates of point on AH, wrt center of AH
              !AH_R,AH_chi,AH_phi polar coordinates of point on AH, wrt center of AH
              x0=AH_R(i,j)*cos(AH_chi)
              y0=AH_R(i,j)*sin(AH_chi)

              xmin=min(xmin,x0)
              xmax=max(xmax,x0)
              ymin=min(ymin,y0)
              ymax=max(ymax,y0)
           end do
        end do

        ! for now, don't move xc (yc) when xc=0 (yc=0)
        if (AH_xc(1).eq.0) then
          dx0=0
        else
          dx0=(xmax+xmin)/2
        end if
        if (AH_xc(2).eq.0) then
          dy0=0
        else
          dy0=(ymax+ymin)/2
        end if

        AH_xc(1)=AH_xc(1)+dx0
        AH_xc(2)=AH_xc(2)+dy0

        if (ltrace) then
          write(*,*) '-----------------------------------'
          write(*,*) 'dx0,AH_xc(1)=',dx0,AH_xc(1)
          write(*,*) 'dy0,AH_xc(2)=',dy0,AH_xc(2)
        end if

        do i=1,AH_Nchi
           AH_chi=(i-1)*dahchi
           do j=1,AH_Nphi
              AH_phi=(j-1)*dahphi

              !x0,y0 cartesian coordinates of point on AH, wrt center of AH
              !AH_R,AH_chi,AH_phi polar coordinates of point on AH, wrt center of AH
              ! FIX: WHY IS CSX,CSY USED HERE? WHY NOT JUST
              ! AH_R=SQRT((ABS(X)-DX)^2+(ABS(Y)-DY)^2)
              x0=AH_R(i,j)*cos(AH_chi)
              y0=AH_R(i,j)*sin(AH_chi)
              csx=x0/AH_R(i,j)
              csy=y0/AH_R(i,j)
              AH_R(i,j)=sqrt((abs(x0)-csx*dx0)**2+
     &                       (abs(y0)-csy*dy0)**2)

           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c More agressive averaging than smooth_ah_r() above
c
c eps not used yet
c 
c-----------------------------------------------------------------------
        subroutine smooth_ah_r_b(AH_R,AH_w1,AH_eps,AH_Nchi,AH_Nphi)
        implicit none
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_eps
        real*8 AH_w1(AH_Nchi,AH_Nphi)
        
        integer i,j,i1,j1,i0,j0
        integer avg_rad
        real*8 maxr,sum,r0,r1,num,r
        parameter (avg_rad=2)

        ! initialize fixed-size variables
        data i,j,i1,j1,i0,j0/0,0,0,0,0,0/

        data maxr,sum,r0,r1,num,r/0.0,0.0,0.0,0.0,0.0,0.0/

        !--------------------------------------------------------------

        do i=1,AH_Nchi
           do j=1,AH_Nphi
              AH_w1(i,j)=AH_R(i,j)
           end do
        end do

        do i=1,AH_Nchi
           do j=1,AH_Nphi
              maxr=0
              sum=0
              num=0
              do i1=i-avg_rad,i+avg_rad
                 do j1=j-avg_rad,j+avg_rad
                    r=((i1-i)**2+(j1-j)**2)**(0.5d0)
                    if (r.lt.(avg_rad+0.01)) then
                       i0=i1
                       j0=j1
                       if (i0.lt.1) i0=2-i0
                       if (i0.gt.AH_Nchi) i0=AH_Nchi-(i0-AH_Nchi)
                       if (j0.lt.1) j0=j0+AH_Nphi-1
                       if (j0.gt.AH_Nphi) j0=j0-AH_Nphi+1
                       sum=sum+(avg_rad+1-r)*AH_w1(i0,j0)
                       num=num+(avg_rad+1-r)
                       if (abs(AH_w1(i0,j0)).gt.abs(maxr)) 
     &                    maxr=AH_w1(i0,j0)
                    end if
                 end do
              end do
              sum=sum/num
              !sum=(sum-maxr)/(num-1)
              AH_R(i,j)=sum
           end do
        end do

        r0=0
        r1=0
        do i=1,AH_Nphi
           r0=r0+AH_R(1,i)
           r1=r1+AH_R(AH_Nchi,i)
        end do
        r0=r0/AH_Nphi
        r1=r1/AH_Nphi
        do i=1,AH_Nphi
           AH_R(1,i)=r0
           AH_R(AH_Nchi,i)=r1
        end do

        do i=1,AH_Nchi
           AH_R(i,1)=AH_R(i,AH_Nphi)
        end do

        return 
        end

