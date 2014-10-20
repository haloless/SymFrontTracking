

#ifndef SDIM
#define SDIM 2
#endif


      program myft
      !
      use const_module
      use geom_module
      use ns_module
      use front_module
      use prob_module
      !
      implicit none
      double precision, allocatable :: u(:,:), v(:,:)
      double precision, allocatable :: p(:,:)
      double precision, allocatable :: ut(:,:), vt(:,:)
      double precision, allocatable :: tmp1(:,:), tmp2(:,:)
      double precision, allocatable :: uu(:,:), vv(:,:) ! node velocity
      !
      double precision, allocatable :: fx(:,:), fy(:,:)
      !
      double precision, allocatable :: x(:), y(:), xh(:), yh(:)
      double precision, allocatable :: r(:), rh(:) ! cell/face radius
      !
      integer, allocatable :: mask(:,:)
      double precision, allocatable :: color(:,:)
      double precision, allocatable :: dens(:,:), visc(:,:)
      !
      integer :: bctype(SDIM,2)
      !
      double precision, allocatable :: xf(:), yf(:)
      double precision, allocatable :: uf(:), vf(:)
      double precision, allocatable :: tx(:), ty(:)
      
      integer :: Nf
      
      !
      integer :: istep, iplot
      double precision :: time
      
      
      integer :: i,j,l
      
      !
      integer :: untin
      namelist /fortin/ &
        coordsys, Lx,Ly, nx,ny, &
        bc_xlo, bc_xhi, bc_ylo, bc_yhi, &
        gx,gy, rho1,rho2, mu1,mu2, sigma, rro, &
        dt, nstep, maxiter, maxerr, beta, &
        Nf, &
        probtype, &
        xc, yc, rad, &
        eta, omega, theta, &
        plot_int
      
      !
      integer :: untchk
      double precision :: chkvol, chkvol0
      
      ! read namelists
      untin = 9
      open(unit=untin, file='prob.fortin', form='formatted', status='old')
      read(unit=untin, nml=fortin)
      close(unit=untin)
      
      
      ! check
      if (coordsys .eq. COORDSYS_CART) then
      elseif (coordsys .eq. COORDSYS_RZ) then
      else
        print *, 'Invalid coordsys=',coordsys
        stop
      endif
      
      bctype(1,1) = bc_xlo
      bctype(1,2) = bc_xhi
      bctype(2,1) = bc_ylo
      bctype(2,2) = bc_yhi
      
      untchk = 19
      open(unit=untchk, file='chk.csv', status='replace')
      write(untchk,*) 'time,vol'

      ! setup grid
      dx = Lx / nx
      dy = Ly / ny
      
      if (bc_xhi == BCTYPE_PER) then
        iumax = nx+1
        jumax = ny+1
      else
        iumax = nx
        jumax = ny+1
      endif
      if (bc_yhi == BCTYPE_PER) then
        ivmax = nx+1
        jvmax = ny+1
      else
        ivmax = nx+1
        jvmax = ny
      endif
      
      print *, 'nx=',nx, 'ny=',ny
      print *, 'dx=',dx, 'dy=',dy
      print *, 'iumax=',iumax, 'jumax=',jumax
      print *, 'ivmax=',ivmax, 'jvmax=',jvmax
      
      ! allocate buffer
      ! staggered velocity has an extra slot for periodic BC
      allocate(u(nx+2,ny+2), v(nx+2,ny+2))      
      allocate(ut(nx+2,ny+2), vt(nx+2,ny+2))
      allocate(p(nx+2,ny+2))
      allocate(tmp1(nx+2,ny+2), tmp2(nx+2,ny+2))
      allocate(fx(nx+2,ny+2), fy(nx+2,ny+2))
      !
      allocate(uu(nx+1,ny+1), vv(nx+1,ny+1))
      !
      allocate(x(nx+2), y(ny+2), xh(nx+1), yh(ny+1))
      allocate(r(nx+2), rh(nx+1))
      !
      allocate(mask(nx+2,ny+2))
      allocate(color(nx+2,ny+2))
      allocate(dens(nx+2,ny+2), visc(nx+2,ny+2))
      !
      allocate(xf(maxnf+2), yf(maxnf+2))
      allocate(uf(maxnf+2), vf(maxnf+2))
      allocate(tx(maxnf+2), ty(maxnf+2))
      
      ! grid position
      do i = 1, nx+2
        x(i) = dx * (dble(i)-1.5d0)
      enddo
      do j = 1, ny+2
        y(j) = dy * (dble(j)-1.5d0)
      enddo
      do i = 1, nx+1
        xh(i) = dx * (dble(i)-1.0d0)
      enddo
      do j = 1, ny+1
        yh(j) = dy * (dble(j)-1.0d0)
      enddo
      ! radius
      if (coordsys == COORDSYS_CART) then
        r(:) = 1.0d0
        rh(:) = 1.0d0
      elseif (coordsys == COORDSYS_RZ) then
        r(:) = x(:)
        rh(:) = xh(:)
      endif
      
      
      ! initialize fluid state
      u(:,:) = 0.0d0
      v(:,:) = 0.0d0
      p(:,:) = 0.0d0
      
      ! initialize density and viscosity
      mask(:,:) = 1
      color(:,:) = 0.0d0
      dens(:,:) = rho1
      visc(:,:) = mu1
      if (probtype == PROBTYPE_DROPLET) then
        do i = 2, nx+1
        do j = 2, ny+1
          if ((x(i)-xc)**2 + (y(j)-yc)**2 < rad**2) then
            color(i,j) = 1.0d0
          endif
        enddo
        enddo
      elseif (probtype == PROBTYPE_RT_INSTAB) then
        do i = 2, nx+1
        do j = 2, ny+1
          if (y(j) > 0.5d0*Ly + eta*Lx*cos(omega*x(i)/Lx+theta)) then
            color(i,j) = 1.0d0
          endif
        enddo
        enddo
      else
        print *, 'Unknown probtype=',probtype
        stop
      endif
      call bcfill_cell_homo(color, nx,ny)
      !color(1,:) = color(2,:)
      !color(nx+2,:) = color(nx+1,:)
      !color(:,1) = color(:,2)
      !color(:,ny+2) = color(:,ny+1)
      call ft_color_prop(color, dens,visc, nx,ny,dx,dy)
      
      
      ! setup front elements
      uf(:) = 0.0d0
      vf(:) = 0.0d0
      if (probtype == PROBTYPE_DROPLET) then
        do l = 2, Nf+1
          xf(l) = xc + rad*sin(PI*(l-2)/(Nf-1))
          yf(l) = yc - rad*cos(PI*(l-2)/(Nf-1))
          !xf(l) = xc - rad*sin(2.0d0*PI*(l-1)/Nf)
          !yf(l) = yc + rad*cos(2.0d0*PI*(l-1)/Nf)
        enddo
        ! two wall nodes
        xf(2) = 0.0d0
        xf(Nf+1) = 0.0d0
      else if (probtype == PROBTYPE_RT_INSTAB) then
        do l = 2, Nf+1
          xf(l) = Lx/(Nf-1) * (l-2)
          yf(l) = 0.5d0*Ly + eta*Lx*cos(omega*xf(l)/Lx+theta)
          !yf(l) = Ly*0.5d0 + Lx*0.1d0*cos(xf(l))
          !yf(l) = Ly*0.5d0 + Lx*0.1d0*cos(2.0d0*PI*xf(l)/(2.0d0*Lx))
        enddo
        ! two wall nodes
        xf(2) = 0.0d0
        xf(Nf+1) = Lx
      else
        print *, 'Unknown probtype=',probtype
        stop
      endif
      ! add two ghost nodes
      call ft_bndry_pos(Nf, xf,yf)
      
      
      ! smooth properties on initialization
      if (.true.) then
        call ft_color_grad(Nf, xf,yf, mask,color, fx,fy, nx,ny, dx,dy)
        call ft_color_smooth(fx,fy, mask,color, dens,visc, r,rh,nx,ny,dx,dy)
      endif
      
      
      ! initial output
      time = 0.0d0
      istep = 0
      iplot = 0
      call vel_node(u,v, uu,vv, nx,ny, dx,dy)
      call output_csv(iplot, time, &
        color,dens,visc, p, uu,vv, &
        x,y, xh,yh, nx,ny)
      call output_front(iplot, xf, yf, uf, vf, tx, ty, Nf)
      iplot = iplot + 1
      
      call chk_color_volume(chkvol0, color, r,rh,nx,ny,dx,dy)
      write(untchk,*) time,',',chkvol0

      
      ! main loop
      do istep = 1, nstep
        print *, 'step=',istep, 'time=',time
        
		! velocity BC
        if (.false.) then
        if (.false.) then
          u(:,1) = -u(:,2)
          u(:,ny+2) = -u(:,ny+1)
          v(1,:) = -v(2,:)
          v(nx+2,:) = -v(nx+1,:)
        else ! free-slip
          u(:,1) = u(:,2)
          u(:,ny+2) = u(:,ny+1)
          v(1,:) = v(2,:)
          v(nx+2,:) = v(nx+1,:)
        endif
        else
          call bcfill_umac(u, nx,ny)
          call bcfill_vmac(v, nx,ny)
        endif
        
        
		!
        call surf_tension(Nf, xf,yf, tx,ty, fx,fy, nx,ny,dx,dy,coordsys)
		
        if (.true.) then
          call chk_avg_dens(rro, dens, r,rh,nx,ny,dx,dy)
          print *, 'rro=',rro
        endif
		
        ! velocity predictor
        ut(:,:) = u(:,:)
        vt(:,:) = v(:,:)
        
        call vel_conv(ut,vt, u,v, dens,visc, r,rh, nx,ny,dx,dy,dt)
        
        call vel_diff(ut,vt, u,v, dens,visc, r,rh, nx,ny,dx,dy,dt)
        
        call vel_jump(ut,vt, fx,fy, dens,visc, nx,ny,dx,dy,dt)
        
        !call pres_solve(ut,vt, dens, tmp1,tmp2, p, r,rh, nx,ny,dx,dy,dt)
        call bcfill_umac(ut, nx,ny)
        call bcfill_vmac(vt, nx,ny)
        call pres_solve2(ut,vt, dens, tmp1,tmp2, p, r,rh, nx,ny,dx,dy,dt)
        call bcfill_cell_homo(p, nx,ny)
		!
        call vel_corr(ut,vt, dens, u,v, p, nx,ny, dx,dy,dt)
        if (.false.) then
        if (.false.) then
          u(:,1) = -u(:,2)
          u(:,ny+2) = -u(:,ny+1)
          v(1,:) = -v(2,:)
          v(nx+2,:) = -v(nx+1,:)
        else ! free-slip
          u(:,1) = u(:,2)
          u(:,ny+2) = u(:,ny+1)
          v(1,:) = v(2,:)
          v(nx+2,:) = v(nx+1,:)
        endif
        else
          call bcfill_umac(u, nx,ny)
          call bcfill_vmac(v, nx,ny)
        endif
        
        !
        call ft_adv(Nf, xf,yf, uf,vf, u,v, nx,ny, dx,dy,dt)
        call ft_recons(Nf, xf,yf, dx,dy)
        !
        call ft_color_grad(Nf, xf,yf, mask,color, fx,fy, nx,ny,dx,dy)
        call ft_color_smooth(fx,fy, mask,color, dens,visc, r,rh, nx,ny,dx,dy)
        
        if (.false.) then
          call chk_color_volume(chkvol, color, r,rh,nx,ny,dx,dy)
          call ft_volconst(Nf,xf,yf,uf,vf, nx,ny,dx,dy,dt, chkvol,chkvol0)
          call ft_color_grad(Nf, xf,yf, mask,color, fx,fy, nx,ny,dx,dy)
          call ft_color_smooth(fx,fy, mask,color, dens,visc, r,rh, nx,ny,dx,dy)
        endif

        
        time = time + dt
        
        if (mod(istep,plot_int)==0 .or. istep==nstep) then
          call vel_node(u,v, uu,vv, nx,ny, dx,dy)
          call output_csv(iplot, time, &
            color,dens,visc, p, uu,vv, &
            x,y, xh,yh, nx,ny)
          call output_front(iplot, xf, yf, uf, vf, tx, ty, Nf);
          iplot = iplot + 1
          
          !
          call chk_color_volume(chkvol, color, r,rh,nx,ny,dx,dy)
          write(untchk,*) time,',',chkvol
          flush(untchk)
        endif
        
      enddo ! end main loop
      
      close(untchk)
      
      end program myft
      
      
      
      
      
      
      
      
      
      
      
      
      subroutine vel_node( &
        u, v, uu, vv, &
        nx, ny, dx, dy)
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy
      double precision, intent(in) :: u(nx+2,ny+2), v(nx+2,ny+2)
      double precision, intent(out) :: uu(nx+1,ny+1), vv(nx+1,ny+1)
      !
      integer :: i, j
      
      do i = 1, nx+1
      do j = 1, ny+1
        uu(i,j) = 0.5d0 * (u(i,j)+u(i,j+1))
        vv(i,j) = 0.5d0 * (v(i,j)+v(i+1,j))
      enddo
      enddo
      
      return
      end subroutine vel_node
      
      subroutine chk_color_volume(&
        chkvol, color, &
        r, rh, nx,ny,dx,dy)
      !use ns_module, only: coordsys
      implicit none
      double precision, intent(out) :: chkvol
      integer, intent(in) :: nx,ny
      double precision, intent(in) :: dx,dy
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      double precision, intent(in) :: color(nx+2,ny+2)
      !
      integer :: i,j
      
      chkvol = 0.0d0
      do i = 2, nx+1
      do j = 2, ny+1
        chkvol = chkvol + color(i,j)*r(i)*dx*dy
      enddo
      enddo
      
      return
      end subroutine chk_color_volume
      
      subroutine chk_avg_dens( &
        chkrho, dens, &
        r, rh, nx,ny,dx,dy)
      implicit none
      double precision, intent(out) :: chkrho
      integer, intent(in) :: nx,ny
      double precision, intent(in) :: dx,dy
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      double precision, intent(in) :: dens(nx+2,ny+2)
      !
      integer :: i,j
      double precision :: vsum
      
      chkrho = 0.0d0
      vsum = 0.0d0
      do i = 2, nx+1
      do j = 2, ny+1
        chkrho = chkrho + dens(i,j)*r(i)*dx*dy
        vsum = vsum + r(i)*dx*dy
      enddo
      enddo
      chkrho = chkrho / vsum
      
      return
      end subroutine chk_avg_dens
      
      

      
      
      
      
      
      






