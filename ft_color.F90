


      subroutine distrib_face( &
        xf, yf, & ! front position
        nfx, nfy, & ! front quantity
        fx, fy, &
        nx, ny, dx, dy)
      implicit none
      double precision, intent(in) :: xf, yf
      double precision, intent(in) :: nfx, nfy
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy
      double precision, intent(out) :: fx(nx+2,ny+2), fy(nx+2,ny+2)
      !
      integer :: ip, jp
      double precision :: ax, ay
      
      ! UMAC face
      ip = floor(xf/dx) + 1
      jp = floor((yf+0.5d0*dy)/dy) + 1
      ax = xf/dx - ip + 1
      ay = (yf+0.5d0*dy)/dy - jp + 1
      fx(ip,jp) = fx(ip,jp) + (1.0d0-ax)*(1.0d0-ay)*nfx/(dx*dy)
      fx(ip+1,jp) = fx(ip+1,jp) + ax*(1.0d0-ay)*nfx/(dx*dy)
      fx(ip,jp+1) = fx(ip,jp+1) + (1.0d0-ax)*ay*nfx/(dx*dy)
      fx(ip+1,jp+1) = fx(ip+1,jp+1) + ax*ay*nfx/(dx*dy)
      
      ! VMAC face
      ip = floor((xf+0.5d0*dx)/dx) + 1
      jp = floor(yf/dy) + 1
      ax = (xf+0.5d0*dx)/dx - ip + 1
      ay = yf/dy - jp + 1
      fy(ip,jp) = fy(ip,jp) + (1.0d0-ax)*(1.0d0-ay)*nfy/(dx*dy)
      fy(ip+1,jp) = fy(ip+1,jp) + ax*(1.0d0-ay)*nfy/(dx*dy)
      fy(ip,jp+1) = fy(ip,jp+1) + (1.0d0-ax)*ay*nfy/(dx*dy)
      fy(ip+1,jp+1) = fy(ip+1,jp+1) + ax*ay*nfy/(dx*dy)
      
      return
      end subroutine distrib_face

      subroutine ft_color_grad( &
        Nf, xf, yf, &
        mask, color, fx, fy, &
        nx, ny, dx, dy)
      !use ns_module, only: rho1, rho2, mu1, mu2
      implicit none
      integer, intent(in) :: Nf
      double precision, intent(in) :: xf(1:Nf+2), yf(1:Nf+2)
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy
      integer, intent(out) :: mask(nx+2,ny+2)
      double precision, intent(out) :: color(nx+2,ny+2)
      double precision, intent(out) :: fx(nx+2,ny+2), fy(nx+2,ny+2)
      !
      integer :: l
      double precision :: xx,yy, nfx, nfy
      integer :: i, j
      double precision, parameter :: drho = 1.0d0 ! phase 0 & 1
      integer, parameter :: nband = 3
      integer :: ip, jp
      
      fx(:,:) = 0.0d0
      fy(:,:) = 0.0d0
      
      mask(:,:) = 1
      
      do l = 2, Nf
        ! normal vector
        !nfx = -0.5d0 * (yf(l+1)-yf(l-1)) * (rho2-rho1)
        !nfy = 0.5d0 * (xf(l+1)-xf(l-1)) * (rho2-rho1)
        !nfx = -0.5d0 * (yf(l+1)-yf(l-1)) * drho
        !nfy = 0.5d0 * (xf(l+1)-xf(l-1)) * drho
        !call distrib_face(xf(l),yf(l), nfx,nfy, fx,fy, nx,ny, dx,dy)
        xx = 0.5d0 * (xf(l+1)+xf(l))
        yy = 0.5d0 * (yf(l+1)+yf(l))
        nfx = -(yf(l+1)-yf(l)) * drho
        nfy = (xf(l+1)-xf(l)) * drho
        call distrib_face(xx,yy, nfx,nfy, fx,fy, nx,ny,dx,dy)
        
        if (.false.) then
        ip = floor(xf(l)/dx) + 1
        jp = floor(yf(l)/dy) + 1
        do i = ip-nband, ip+nband
        do j = jp-nband, jp+nband
          if (2<=i .and. i<=nx+1 .and. 2<=j .and. j<=ny+1) then
            mask(i,j) = 1
          endif
        enddo
        enddo
        endif
      enddo ! end loop front
      
      if (.true.) then ! enforce symmetry BC
        do j = 2, ny+1
          fx(1,j) = 0.0d0
          fx(nx+1,j) = 0.0d0
        enddo
        do i = 2, nx
          fx(i,2) = fx(i,2) + fx(i,1)
          fx(i,ny+1) = fx(i,ny+1) + fx(i,ny+2)
        enddo
        
        do i = 2, nx+1
          fy(i,1) = 0.0d0
          fy(i,ny+1) = 0.0d0
        enddo
        do j = 2, ny
          fy(2,j) = fy(2,j) + fy(1,j)
          fy(nx+1,j) = fy(nx+1,j) + fy(nx+2,j)
        enddo
      endif
      
      return
      end subroutine ft_color_grad
      
      subroutine ft_color_smooth( &
        fx, fy, &
        mask, color, &
        dens, visc, &
        r, rh, &
        nx, ny, dx, dy)
      use ns_module, only: maxiter,maxerr,beta
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      double precision, intent(in) :: fx(1:nx+2,1:ny+2), fy(1:nx+2,1:ny+2)
      integer, intent(in) :: mask(nx+2,ny+2)
      double precision, intent(out) :: color(nx+2,ny+2)
      double precision, intent(out) :: dens(1:nx+2,1:ny+2), visc(1:nx+2,1:ny+2)
      !
      integer :: i,j
      double precision :: rrhs(nx+2,ny+2), rdiag(nx+2,ny+2)
      double precision :: rcoefx(nx+1,ny+2), rcoefy(nx+2,ny+1)
      double precision :: rsave(nx+2,ny+2)
      integer :: iter
      double precision :: drmax, rval
      
      rrhs(:,:) = 0.0d0
      rdiag(:,:) = 0.0d0
      do i = 2, nx+1
      do j = 2, ny+1
        !rrhs(i,j) = dx*fx(i-1,j) - dx*fx(i,j) + dy*fy(i,j-1) - dy*fy(i,j)
        rrhs(i,j) = dy * (rh(i)*fx(i,j) - rh(i-1)*fx(i-1,j)) + &
          dx*r(i) * (fy(i,j)-fy(i,j-1))
        
        if (i .ne. 2) then
          rdiag(i,j) = rdiag(i,j) + dy/dx*rh(i-1)
        endif
        if (i .ne. nx+1) then
          rdiag(i,j) = rdiag(i,j) + dy/dx*rh(i)
        endif
        if (j .ne. 2) then
          rdiag(i,j) = rdiag(i,j) + dx/dy*r(i)
        endif
        if (j .ne. ny+1) then
          rdiag(i,j) = rdiag(i,j) + dx/dy*r(i)
        endif
      enddo
      enddo
      
      rcoefx(:,:) = 0.0d0
      do i = 2, nx
      do j = 2, ny+1
        rcoefx(i,j) = dy/dx * rh(i)
      enddo
      enddo
      rcoefy(:,:) = 0.0d0
      do i = 2, nx+1
      do j = 2, ny
        rcoefy(i,j) = dx/dy * r(i)
      enddo
      enddo
      
      do iter = 1, maxiter
        rsave(:,:) = color(:,:)
        drmax = 0.0d0
        
        do i = 2, nx+1
        do j = 2, ny+1
          if (mask(i,j) == 1) then
            !dens(i,j) = (1.0d0-beta)*dens(i,j) + beta* &
            !  0.25d0 * (dens(i+1,j)+dens(i-1,j)+dens(i,j+1)+dens(i,j-1) &
            !  + rrhs(i,j))
            !rval = 0.25d0 * (color(i+1,j) + color(i-1,j) + &
            !  color(i,j+1)+color(i,j-1) + rrhs(i,j))
            
            if (.true.) then ! GS iter.
            rval = rcoefx(i-1,j)*color(i-1,j) + rcoefx(i,j)*color(i+1,j) &
              + rcoefy(i,j-1)*color(i,j-1) + rcoefy(i,j)*color(i,j+1)
            else ! Jacobi iter.
            rval = rcoefx(i-1,j)*rsave(i-1,j) + rcoefx(i,j)*rsave(i+1,j) &
              + rcoefy(i,j-1)*rsave(i,j-1) + rcoefy(i,j)*rsave(i,j+1)
            endif
            rval = (rval-rrhs(i,j)) / rdiag(i,j)
            
            color(i,j) = (1.0d0-beta)*color(i,j) + beta*rval
            
            !drmax = max(drmax, abs(rsave(i,j)-dens(i,j)))
            drmax = max(drmax, abs(rsave(i,j)-color(i,j)))
          endif
        enddo
        enddo
        
        if (drmax < maxerr) then
          print *,"Color convergence, iter=", iter
          exit
        endif
      enddo
      
      ! cutoff
      do i = 2, nx+1
      do j = 2, ny+1
        color(i,j) = min(max(0.0d0,color(i,j)),1.0d0)
      enddo
      enddo
      
      ! homogeneous BC
      color(1,:) = color(2,:)
      color(nx+2,:) = color(nx+1,:)
      color(:,1) = color(:,2)
      color(:,ny+2) = color(:,ny+1)
      
      call ft_color_prop(color, dens, visc, nx,ny,dx,dy)      
      
      return
      end subroutine ft_color_smooth
      
      
      subroutine ft_color_prop( &
        color, dens, visc, &
        nx, ny, dx, dy)
      use ns_module, only: rho1, rho2, mu1, mu2
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy
      double precision, intent(in) :: color(nx+2,ny+2)
      double precision, intent(out) :: dens(nx+2,ny+2), visc(nx+2,ny+2)
      !
      integer :: i,j
      
      do i = 1, nx+2
      do j = 1, ny+2
        dens(i,j) = rho1 + (rho2-rho1)*color(i,j)
        visc(i,j) = mu1 + (mu2-mu1)*color(i,j)
      enddo
      enddo
      
      return
      end subroutine ft_color_prop
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

