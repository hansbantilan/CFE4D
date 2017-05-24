c----------------------------------------------------------------------
c in cartesian coordinates t,x,y
c
c routines for computing y=0 axis regularity conditions for
c metric and scalar fied
c----------------------------------------------------------------------
        subroutine axi_reg_g(gb_tt,gb_tx,gb_ty,
     &                       gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
     &                       L,x,y,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        integer regtype
        real*8 gb_tt(Nx,Ny),gb_tx(Nx,Ny),gb_ty(Nx,Ny)
        real*8 gb_xx(Nx,Ny),gb_xy(Nx,Ny),gb_yy(Nx,Ny)
        real*8 psi(Nx,Ny),tfunction(Nx,Ny),chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)

        integer i
        real*8 dx,dy
        real*8 PI
        parameter (PI=3.141592653589793d0)

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (abs(y(1)).gt.dy/2) return
 
        do i=1,Nx
           if (chr(i,1).ne.ex) then

              ! 1-pt axis bcs, using 3-pt zero derivative to get axis pt; 
              ! set psi=gb_yy at y=0
              if (regtype.eq.1) then 
                gb_tt(i,1)=(4*gb_tt(i,2)-gb_tt(i,3))/3
                gb_tx(i,1)=(4*gb_tx(i,2)-gb_tx(i,3))/3
                gb_xx(i,1)=(4*gb_xx(i,2)-gb_xx(i,3))/3
                gb_yy(i,1)=(4*gb_yy(i,2)-gb_yy(i,3))/3
                psi(i,1)=gb_yy(i,1)
                gb_xy(i,1)=0
                gb_ty(i,1)=0
              ! 1-pt axis bcs, using 3-pt zero derivative to get axis pt; 
              ! set gb_yy=psi at y=0
              else if (regtype.eq.2) then 
                gb_tt(i,1)=(4*gb_tt(i,2)-gb_tt(i,3))/3
                gb_tx(i,1)=(4*gb_tx(i,2)-gb_tx(i,3))/3
                gb_xx(i,1)=(4*gb_xx(i,2)-gb_xx(i,3))/3
                psi(i,1)=(4*psi(i,2)-psi(i,3))/3
                gb_yy(i,1)=psi(i,1)
                gb_xy(i,1)=0
                gb_ty(i,1)=0
              !(experimental 2-pt axis bcs)!
              else if (regtype.eq.3) then  
                gb_tt(i,1)=(4*gb_tt(i,3)-gb_tt(i,5))/3
                gb_tx(i,1)=(4*gb_tx(i,3)-gb_tx(i,5))/3
                gb_xx(i,1)=(4*gb_xx(i,3)-gb_xx(i,5))/3
                gb_yy(i,1)=(4*gb_yy(i,3)-gb_yy(i,5))/3
                psi(i,1)=gb_yy(i,1)
                gb_xy(i,1)=0
                gb_ty(i,1)=0
                gb_tt(i,2)= gb_tt(i,1)/4 + 3*gb_tt(i,3)/2
     &                     -gb_tt(i,4)   +   gb_tt(i,5)/4
                gb_tx(i,2)= gb_tx(i,1)/4 + 3*gb_tx(i,3)/2
     &                     -gb_tx(i,4)   +   gb_tx(i,5)/4
                gb_xx(i,2)= gb_xx(i,1)/4 + 3*gb_xx(i,3)/2
     &                     -gb_xx(i,4)   +   gb_xx(i,5)/4
                gb_yy(i,2)= gb_yy(i,1)/4 + 3*gb_yy(i,3)/2
     &                     -gb_yy(i,4)   +   gb_yy(i,5)/4
                psi(i,2)  = psi(i,1)/4   + 3*psi(i,3)/2
     &                     -psi(i,4)     +   psi(i,5)/4
                gb_xy(i,2)=0.5*gb_xy(i,3)
                gb_ty(i,2)=0.5*gb_ty(i,3)
              ! 2-pt axis bcs, using 5-pt zero derivative to get non-axis pt; 
              ! set psi=gb_yy at y=0
              else if (regtype.eq.4) then  
                gb_tt(i,1)=(4*gb_tt(i,3)-gb_tt(i,5))/3
                gb_tx(i,1)=(4*gb_tx(i,3)-gb_tx(i,5))/3
                gb_xx(i,1)=(4*gb_xx(i,3)-gb_xx(i,5))/3
                gb_yy(i,1)=(4*gb_yy(i,3)-gb_yy(i,5))/3
                psi(i,1)=gb_yy(i,1)
                gb_xy(i,1)=0
                gb_ty(i,1)=0
                gb_tt(i,2)=(
     &            25*gb_tt(i,1)+36*gb_tt(i,3)-16*gb_tt(i,4)+3*gb_tt(i,5)
     &            )/48
                gb_tx(i,2)=(
     &            25*gb_tx(i,1)+36*gb_tx(i,3)-16*gb_tx(i,4)+3*gb_tx(i,5)
     &            )/48
                gb_xx(i,2)=(
     &            25*gb_xx(i,1)+36*gb_xx(i,3)-16*gb_xx(i,4)+3*gb_xx(i,5)
     &            )/48
                gb_yy(i,2)=(
     &            25*gb_yy(i,1)+36*gb_yy(i,3)-16*gb_yy(i,4)+3*gb_yy(i,5)
     &            )/48
                psi(i,2)=(
     &            25*psi(i,1)+36*psi(i,3)-16*psi(i,4)+3*psi(i,5)
     &            )/48
                gb_xy(i,2)=0.5*gb_xy(i,3)
                gb_ty(i,2)=0.5*gb_ty(i,3)
              ! 2-pt axis bcs, using 5-pt zero derivative to get non-axis pt; 
              ! set gb_yy=psi at y=0
              else if (regtype.eq.5) then  
                gb_tt(i,1)=(4*gb_tt(i,3)-gb_tt(i,5))/3
                gb_tx(i,1)=(4*gb_tx(i,3)-gb_tx(i,5))/3
                gb_xx(i,1)=(4*gb_xx(i,3)-gb_xx(i,5))/3
                gb_yy(i,1)=(4*gb_yy(i,3)-gb_yy(i,5))/3
                gb_yy(i,1)=psi(i,1)
                gb_xy(i,1)=0
                gb_ty(i,1)=0
                gb_tt(i,2)=(
     &            25*gb_tt(i,1)+36*gb_tt(i,3)-16*gb_tt(i,4)+3*gb_tt(i,5)
     &            )/48
                gb_tx(i,2)=(
     &            25*gb_tx(i,1)+36*gb_tx(i,3)-16*gb_tx(i,4)+3*gb_tx(i,5)
     &            )/48
                gb_xx(i,2)=(
     &            25*gb_xx(i,1)+36*gb_xx(i,3)-16*gb_xx(i,4)+3*gb_xx(i,5)
     &            )/48
                gb_yy(i,2)=(
     &            25*gb_yy(i,1)+36*gb_yy(i,3)-16*gb_yy(i,4)+3*gb_yy(i,5)
     &            )/48
                psi(i,2)=(
     &            25*psi(i,1)+36*psi(i,3)-16*psi(i,4)+3*psi(i,5)
     &            )/48
                gb_xy(i,2)=0.5*gb_xy(i,3)
                gb_ty(i,2)=0.5*gb_ty(i,3)
              ! 3-pt axis bcs, using 5-pt zero derivative to get non-axis pts; 
              ! set psi=gb_yy at y=0
              else if (regtype.eq.6) then  
                gb_tt(i,1)=(4*gb_tt(i,4)-gb_tt(i,7))/3
                gb_tx(i,1)=(4*gb_tx(i,4)-gb_tx(i,7))/3
                gb_xx(i,1)=(4*gb_xx(i,4)-gb_xx(i,7))/3
                gb_yy(i,1)=(4*gb_yy(i,4)-gb_yy(i,7))/3
                psi(i,1)=gb_yy(i,1)
                gb_xy(i,1)=0
                gb_ty(i,1)=0
                gb_tt(i,2)=( 107*gb_tt(i,1) + 100*gb_tt(i,4)
     &                       -75*gb_tt(i,5) + 18*gb_tt(i,6) )/150
                gb_tx(i,2)=( 107*gb_tx(i,1) + 100*gb_tx(i,4)
     &                       -75*gb_tx(i,5) + 18*gb_tx(i,6) )/150
                gb_xx(i,2)=( 107*gb_xx(i,1) + 100*gb_xx(i,4)
     &                       -75*gb_xx(i,5) + 18*gb_xx(i,6) )/150
                gb_yy(i,2)=( 107*gb_yy(i,1) + 100*gb_yy(i,4)
     &                       -75*gb_yy(i,5) + 18*gb_yy(i,6) )/150
                psi(i,2)  =( 107*psi(i,1)   + 100*psi(i,4)
     &                       -75*psi(i,5)   + 18*psi(i,6) )/150
                gb_xy(i,2)=gb_xy(i,4)/3
                gb_ty(i,2)=gb_ty(i,4)/3
                gb_tt(i,3)=( 77*gb_tt(i,1)   + 400*gb_tt(i,4)
     &                       -225*gb_tt(i,5) + 48*gb_tt(i,6) )/300
                gb_tx(i,3)=( 77*gb_tx(i,1)   + 400*gb_tx(i,4)
     &                       -225*gb_tx(i,5) + 48*gb_tx(i,6) )/300
                gb_xx(i,3)=( 77*gb_xx(i,1)   + 400*gb_xx(i,4)
     &                       -225*gb_xx(i,5) + 48*gb_xx(i,6) )/300
                gb_yy(i,3)=( 77*gb_yy(i,1)   + 400*gb_yy(i,4)
     &                       -225*gb_yy(i,5) + 48*gb_yy(i,6) )/300
                psi(i,3)  =( 77*psi(i,1)     + 400*psi(i,4)
     &                       -225*psi(i,5)   + 48*psi(i,6) )/300
                gb_xy(i,3)=2*gb_xy(i,4)/3
                gb_ty(i,3)=2*gb_ty(i,4)/3
              ! 3-pt axis bcs, using 5-pt zero derivative to get non-axis pts; 
              ! set gb_yy=psi at y=0
              else if (regtype.eq.7) then  
                gb_tt(i,1)=(4*gb_tt(i,4)-gb_tt(i,7))/3
                gb_tx(i,1)=(4*gb_tx(i,4)-gb_tx(i,7))/3
                gb_xx(i,1)=(4*gb_xx(i,4)-gb_xx(i,7))/3
                psi(i,1)=(4*psi(i,4)-psi(i,7))/3
                gb_yy(i,1)=psi(i,1)
                gb_xy(i,1)=0
                gb_ty(i,1)=0
                gb_tt(i,2)=( 107*gb_tt(i,1) + 100*gb_tt(i,4)
     &                       -75*gb_tt(i,5) + 18*gb_tt(i,6) )/150
                gb_tx(i,2)=( 107*gb_tx(i,1) + 100*gb_tx(i,4)
     &                       -75*gb_tx(i,5) + 18*gb_tx(i,6) )/150
                gb_xx(i,2)=( 107*gb_xx(i,1) + 100*gb_xx(i,4)
     &                       -75*gb_xx(i,5) + 18*gb_xx(i,6) )/150
                gb_yy(i,2)=( 107*gb_yy(i,1) + 100*gb_yy(i,4)
     &                       -75*gb_yy(i,5) + 18*gb_yy(i,6) )/150
                psi(i,2)  =( 107*psi(i,1)   + 100*psi(i,4)
     &                       -75*psi(i,5)   + 18*psi(i,6) )/150
                gb_xy(i,2)=gb_xy(i,4)/3
                gb_ty(i,2)=gb_ty(i,4)/3
                gb_tt(i,3)=( 77*gb_tt(i,1)   + 400*gb_tt(i,4)
     &                       -225*gb_tt(i,5) + 48*gb_tt(i,6) )/300
                gb_tx(i,3)=( 77*gb_tx(i,1)   + 400*gb_tx(i,4)
     &                       -225*gb_tx(i,5) + 48*gb_tx(i,6) )/300
                gb_xx(i,3)=( 77*gb_xx(i,1)   + 400*gb_xx(i,4)
     &                       -225*gb_xx(i,5) + 48*gb_xx(i,6) )/300
                gb_yy(i,3)=( 77*gb_yy(i,1)   + 400*gb_yy(i,4)
     &                       -225*gb_yy(i,5) + 48*gb_yy(i,6) )/300
                psi(i,3)  =( 77*psi(i,1)     + 400*psi(i,4)
     &                       -225*psi(i,5)   + 48*psi(i,6) )/300
                gb_xy(i,3)=2*gb_xy(i,4)/3
                gb_ty(i,3)=2*gb_ty(i,4)/3
              else
                write(*,*) 'axi_reg_g error: invalid regtype'
                write(*,*) 'regtype=',regtype
                stop
              endif

           else

              gb_tt(i,1)=0 
              gb_tx(i,1)=0
              gb_xx(i,1)=0
              gb_yy(i,1)=0
              gb_xy(i,1)=0
              gb_ty(i,1)=0
              psi(i,1)=0

           end if
        end do

        return
        end

c----------------------------------------------------------------------
        subroutine axi_reg_phi(phi,chr,ex,L,x,y,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        real*8 phi(Nx,Ny),chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)

        real*8 PI,dx,dy
        parameter (PI=3.141592653589793d0)
        integer i

        integer regtype

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (abs(y(1)).gt.dy/2) return

        do i=1,Nx
           if (chr(i,1).ne.ex) then

              if (regtype.eq.1) then
                phi(i,1)=(4*phi(i,2)-phi(i,3))/3
              else if (regtype.eq.2) then
                phi(i,1)=(4*phi(i,2)-phi(i,3))/3
              else if (regtype.eq.3) then
                phi(i,1)=(4*phi(i,3)-phi(i,5))/3
                phi(i,2)  = phi(i,1)/4   + 3*phi(i,3)/2
     &                     -phi(i,4)     +   phi(i,5)/4
              else if (regtype.eq.4) then
                phi(i,1)=(4*phi(i,3)-phi(i,5))/3
                phi(i,2)=(
     &             phi(i,1)+6*phi(i,3)-4*phi(i,4)+phi(i,5)
     &             )/4
              else if (regtype.eq.5) then
                phi(i,1)=(4*phi(i,3)-phi(i,5))/3
                phi(i,2)=(
     &            25*phi(i,1)+36*phi(i,3)-16*phi(i,4)+3*phi(i,5)
     &            )/48
              else if (regtype.eq.6) then
                phi(i,1)=(4*phi(i,4)-phi(i,7))/3
                phi(i,2)  =( 107*phi(i,1)   + 100*phi(i,4)
     &                       -75*phi(i,5)   + 18*phi(i,6) )/150
                phi(i,3)  =( 77*phi(i,1)     + 400*phi(i,4)
     &                       -225*phi(i,5)   + 48*phi(i,6) )/300
              else if (regtype.eq.7) then
                phi(i,1)=(4*phi(i,4)-phi(i,7))/3
                phi(i,2)  =( 107*phi(i,1)   + 100*phi(i,4)
     &                       -75*phi(i,5)   + 18*phi(i,6) )/150
                phi(i,3)  =( 77*phi(i,1)     + 400*phi(i,4)
     &                       -225*phi(i,5)   + 48*phi(i,6) )/300
              else
                write(*,*) 'axi_reg_phi error: invalid regtype'
                write(*,*) 'regtype=',regtype
                stop
              endif
           
           else 

              phi(i,1)=0

           end if
        end do
           
        return
        end

c----------------------------------------------------------------------
        subroutine axi_reg_Hb(Hb_t,Hb_x,Hb_y,chr,ex,L,x,y,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        real*8 Hb_t(Nx,Ny),Hb_x(Nx,Ny),Hb_y(Nx,Ny)
        real*8 gb_xx(Nx,Ny),gb_xy(Nx,Ny),gb_yy(Nx,Ny)
        real*8 psi(Nx,Ny),chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)

        real*8 PI,dx,dy
        parameter (PI=3.141592653589793d0)
        integer i
 
        integer regtype

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (abs(y(1)).gt.dy/2) return

        do i=1,Nx
           if (chr(i,1).ne.ex) then

              if (regtype.eq.1) then
                Hb_t(i,1)=(4*Hb_t(i,2)-Hb_t(i,3))/3
                Hb_x(i,1)=(4*Hb_x(i,2)-Hb_x(i,3))/3
                Hb_y(i,1)=0 
              else if (regtype.eq.2) then
                Hb_t(i,1)=(4*Hb_t(i,2)-Hb_t(i,3))/3
                Hb_x(i,1)=(4*Hb_x(i,2)-Hb_x(i,3))/3
                Hb_y(i,1)=0 
              else if (regtype.eq.3) then
                Hb_t(i,1)=(4*Hb_t(i,3)-Hb_t(i,5))/3
                Hb_x(i,1)=(4*Hb_x(i,3)-Hb_x(i,5))/3
                Hb_y(i,1)=0
                Hb_t(i,2)  = Hb_t(i,1)/4   + 3*Hb_t(i,3)/2
     &                     -Hb_t(i,4)     +   Hb_t(i,5)/4
                Hb_x(i,2)  = Hb_x(i,1)/4   + 3*Hb_x(i,3)/2
     &                     -Hb_x(i,4)     +   Hb_x(i,5)/4
                Hb_y(i,2)=0.5*Hb_y(i,3)
              else if (regtype.eq.4) then
                Hb_t(i,1)=(4*Hb_t(i,3)-Hb_t(i,5))/3
                Hb_x(i,1)=(4*Hb_x(i,3)-Hb_x(i,5))/3
                Hb_y(i,1)=0
                Hb_t(i,2)=(
     &             Hb_t(i,1)+6*Hb_t(i,3)-4*Hb_t(i,4)+Hb_t(i,5)
     &             )/4
                Hb_x(i,2)=(
     &             Hb_x(i,1)+6*Hb_x(i,3)-4*Hb_x(i,4)+Hb_x(i,5)
     &             )/4
                Hb_y(i,2)=0.5*Hb_y(i,3)
              else if (regtype.eq.5) then
                Hb_t(i,1)=(4*Hb_t(i,3)-Hb_t(i,5))/3
                Hb_x(i,1)=(4*Hb_x(i,3)-Hb_x(i,5))/3
                Hb_y(i,1)=0
                Hb_t(i,2)=(
     &            25*Hb_t(i,1)+36*Hb_t(i,3)-16*Hb_t(i,4)+3*Hb_t(i,5)
     &            )/48
                Hb_x(i,2)=(
     &            25*Hb_x(i,1)+36*Hb_x(i,3)-16*Hb_x(i,4)+3*Hb_x(i,5)
     &            )/48
                Hb_y(i,2)=0.5*Hb_y(i,3)
              else if (regtype.eq.6) then
                Hb_t(i,1)=(4*Hb_t(i,4)-Hb_t(i,7))/3
                Hb_x(i,1)=(4*Hb_x(i,4)-Hb_x(i,7))/3
                Hb_y(i,1)=0
                Hb_t(i,2)  =( 107*Hb_t(i,1)   + 100*Hb_t(i,4)
     &                       -75*Hb_t(i,5)   + 18*Hb_t(i,6) )/150
                Hb_x(i,2)  =( 107*Hb_x(i,1)   + 100*Hb_x(i,4)
     &                       -75*Hb_x(i,5)   + 18*Hb_x(i,6) )/150
                Hb_y(i,2)=Hb_y(i,4)/3
                Hb_t(i,3)  =( 77*Hb_t(i,1)     + 400*Hb_t(i,4)
     &                       -225*Hb_t(i,5)   + 48*Hb_t(i,6) )/300
                Hb_x(i,3)  =( 77*Hb_x(i,1)     + 400*Hb_x(i,4)
     &                       -225*Hb_x(i,5)   + 48*Hb_x(i,6) )/300
                Hb_y(i,3)=2*Hb_y(i,4)/3
              else if (regtype.eq.7) then
                Hb_t(i,1)=(4*Hb_t(i,4)-Hb_t(i,7))/3
                Hb_x(i,1)=(4*Hb_x(i,4)-Hb_x(i,7))/3
                Hb_y(i,1)=0
                Hb_t(i,2)  =( 107*Hb_t(i,1)   + 100*Hb_t(i,4)
     &                       -75*Hb_t(i,5)   + 18*Hb_t(i,6) )/150
                Hb_x(i,2)  =( 107*Hb_x(i,1)   + 100*Hb_x(i,4)
     &                       -75*Hb_x(i,5)   + 18*Hb_x(i,6) )/150
                Hb_y(i,2)=Hb_y(i,4)/3
                Hb_t(i,3)  =( 77*Hb_t(i,1)     + 400*Hb_t(i,4)
     &                       -225*Hb_t(i,5)   + 48*Hb_t(i,6) )/300
                Hb_x(i,3)  =( 77*Hb_x(i,1)     + 400*Hb_x(i,4)
     &                       -225*Hb_x(i,5)   + 48*Hb_x(i,6) )/300
                Hb_y(i,3)=2*Hb_y(i,4)/3
              else
                write(*,*) 'axi_reg_Hb error: invalid regtype'
                write(*,*) 'regtype=',regtype
                stop
              endif

           else 

              Hb_t(i,1)=0
              Hb_x(i,1)=0
              Hb_y(i,1)=0

           end if
        end do
           
        return
        end

c----------------------------------------------------------------------
        subroutine ref_sym_g(gb_tt,gb_tx,gb_ty,
     &                       gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
     &                       L,x,y,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        integer regtype
        real*8 gb_tt(Nx,Ny),gb_tx(Nx,Ny),gb_ty(Nx,Ny)
        real*8 gb_xx(Nx,Ny),gb_xy(Nx,Ny),gb_yy(Nx,Ny)
        real*8 psi(Nx,Ny),tfunction(Nx,Ny),chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)

        integer j
        real*8 dx,dy
        real*8 PI
        parameter (PI=3.141592653589793d0)

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (abs(x(1)).gt.dx/2) return
 
        do j=1,Ny
           if (chr(1,j).ne.ex) then

              gb_tt(1,j)=(4*gb_tt(2,j)-gb_tt(3,j))/3
              gb_tx(1,j)=(4*gb_tx(2,j)-gb_tx(3,j))/3
              gb_ty(1,j)=(4*gb_ty(2,j)-gb_ty(3,j))/3
              gb_xx(1,j)=(4*gb_xx(2,j)-gb_xx(3,j))/3
              gb_xy(1,j)=(4*gb_xy(2,j)-gb_xy(3,j))/3
              gb_yy(1,j)=(4*gb_yy(2,j)-gb_yy(3,j))/3
              psi(1,j)=(4*psi(2,j)-psi(3,j))/3

           end if
        end do

        return
        end

c----------------------------------------------------------------------
        subroutine ref_sym_phi(phi,chr,ex,L,x,y,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        real*8 phi(Nx,Ny),chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)

        real*8 PI,dx,dy
        parameter (PI=3.141592653589793d0)
        integer j

        integer regtype

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (abs(x(1)).gt.dx/2) return

        do j=1,Ny
           if (chr(1,j).ne.ex) then

              phi(1,j)=(4*phi(2,j)-phi(3,j))/3

           end if
        end do

        return
        end

c----------------------------------------------------------------------
        subroutine ref_sym_Hb(Hb_t,Hb_x,Hb_y,chr,ex,L,x,y,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        real*8 Hb_t(Nx,Ny),Hb_x(Nx,Ny),Hb_y(Nx,Ny)
        real*8 gb_xx(Nx,Ny),gb_xy(Nx,Ny),gb_yy(Nx,Ny)
        real*8 psi(Nx,Ny),chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)

        real*8 PI,dx,dy
        parameter (PI=3.141592653589793d0)
        integer j

        integer regtype

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (abs(x(1)).gt.dx/2) return

        do j=1,Ny
           if (chr(1,j).ne.ex) then

              Hb_t(1,j)=(4*Hb_t(2,j)-Hb_t(3,j))/3
              Hb_x(1,j)=(4*Hb_x(2,j)-Hb_x(3,j))/3
              Hb_y(1,j)=(4*Hb_y(2,j)-Hb_y(3,j))/3

           end if
        end do

        return
        end
