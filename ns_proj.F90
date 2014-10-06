
#ifndef SDIM
#define SDIM 2
#endif


      subroutine pres_solve( &
        ustar, vstar, &
        dens, &
        tmp1, tmp2, &
        p, &
        r, rh, &
        nx, ny, dx, dy, dt)
      use ns_module, only: maxiter, maxerr, beta
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy, dt
      double precision, intent(in) :: ustar(1:nx+1,1:ny+2), vstar(1:nx+2,1:ny+1)
      double precision, intent(in) :: dens(1:nx+2,1:ny+2)
      double precision, intent(out) :: tmp1(nx+2,ny+2), tmp2(nx+2,ny+2) ! RHS and DIAG
      double precision, intent(out) :: p(1:nx+2,1:ny+2)
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      !
      integer :: i, j
      double precision, parameter :: large = 1.0d4
      double precision :: re, rw, rn, rs, rp
      integer :: iter
      double precision :: dpmax, pval
      double precision :: rtmp(1:nx+2,1:ny+2)
      double precision :: coefx(nx+1,ny+2), coefy(nx+2,ny+1)
      double precision :: psave(1:nx+2,1:ny+2)
      !
      double precision :: sol(2:nx+1,2:ny+1), rhs(2:nx+1,2:ny+1)
      double precision :: mat(0:4,2:nx+1,2:ny+1)
      
      rtmp(:,:) = dens(:,:)
      rtmp(:,1) = large
      rtmp(:,ny+2) = large
      rtmp(1,:) = large
      rtmp(nx+2,:) = large
      
      ! PPE diagonal and RHS
      tmp1(:,:) = 0.0d0
      tmp2(:,:) = 0.0d0
      do i = 2, nx+1
      do j = 2, ny+1
        tmp1(i,j) = 1.0d0/dt * ( &
          dy * (rh(i)*ustar(i,j) - rh(i-1)*ustar(i-1,j)) + &
          dx * (r(i)*vstar(i,j) - r(i)*vstar(i,j-1)))
        tmp1(i,j) = -tmp1(i,j)
        
        if (i .ne. 2) then
          rw = 0.5d0 * (dens(i,j)+dens(i-1,j))
          tmp2(i,j) = tmp2(i,j) + dy/dx*rh(i-1)/rw
        endif
        if (i .ne. nx+1) then
          re = 0.5d0 * (dens(i,j)+dens(i+1,j))
          tmp2(i,j) = tmp2(i,j) + dy/dx*rh(i)/re
        endif
        if (j .ne. 2) then
          rs = 0.5d0 * (dens(i,j)+dens(i,j-1))
          tmp2(i,j) = tmp2(i,j) + dx/dy*r(i)/rs
        endif
        if (j .ne. ny+1) then
          rn = 0.5d0 * (dens(i,j)+dens(i,j+1))
          tmp2(i,j) = tmp2(i,j) + dx/dy*r(i)/rn
        endif
        tmp2(i,j) = 1.0d0 / tmp2(i,j)
      enddo
      enddo
      
      !
      coefx(:,:) = 0.0d0
      do i = 2, nx
      do j = 2, ny+1
        rp = 0.5d0 * (dens(i,j)+dens(i+1,j))
        coefx(i,j) = dy/dx * rh(i) / rp
      enddo
      enddo
      coefy(:,:) = 0.0d0
      do i = 2, nx+1
      do j = 2, ny
        rp = 0.5d0 * (dens(i,j)+dens(i,j+1))
        coefy(i,j) = dx/dy * r(i) / rp
      enddo
      enddo
      
      ! solve
      if (.false.) then
      do iter = 1, maxiter
        psave(:,:) = p(:,:)
        dpmax = 0.0d0
        
        do i = 2, nx+1
        do j = 2, ny+1
          if (.false.) then
          re = 0.5d0 * (rtmp(i+1,j)+rtmp(i,j))
          rw = 0.5d0 * (rtmp(i-1,j)+rtmp(i,j))
          rn = 0.5d0 * (rtmp(i,j+1)+rtmp(i,j))
          rs = 0.5d0 * (rtmp(i,j-1)+rtmp(i,j))
          p(i,j) = (1.0d0-beta)*p(i,j) + beta*tmp2(i,j) * ( &
            1.0d0/dx**2 * (p(i+1,j)/re + p(i-1,j)/rw) + &
            1.0d0/dy**2 * (p(i,j+1)/rn + p(i,j-1)/rs) - &
            tmp1(i,j))
          else
            pval = coefx(i-1,j)*p(i-1,j) + coefx(i,j)*p(i+1,j) &
              + coefy(i,j-1)*p(i,j-1) + coefy(i,j)*p(i,j+1)
            pval = (pval-tmp1(i,j)) * tmp2(i,j)
            p(i,j) = (1.0d0-beta)*p(i,j) + beta*pval
          endif
          dpmax = max(dpmax, abs(psave(i,j)-p(i,j)))
        enddo
        enddo

        if (dpmax < maxerr) then
          print *, "Pressure convergence, iter=", iter
          exit
        endif
      enddo
      else
        tmp2(:,:) = 0.0d0
        
        call solve_poisson(tmp2,coefx,coefy, tmp1, p, nx,ny)
      endif
      
      !print *, p
      return
      end subroutine pres_solve
      
      subroutine pres_solve2( &
        ustar, vstar, &
        dens, &
        tmp1, tmp2, &
        p, &
        r, rh, &
        nx, ny, dx, dy, dt)
      use const_module
      use geom_module, only: bc_xlo,bc_ylo,bc_xhi,bc_yhi
      !
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy, dt
      double precision, intent(in) :: ustar(nx+2,ny+2), vstar(nx+2,ny+2)
      double precision, intent(in) :: dens(nx+2,ny+2)
      double precision, intent(out) :: tmp1(nx+2,ny+2), tmp2(nx+2,ny+2) ! RHS and DIAG
      double precision, intent(out) :: p(nx+2,ny+2)
      double precision, intent(in) :: r(nx+2), rh(nx+1)
      !
      integer :: i, j
      double precision :: rho, coef
      !
      integer :: periodic(SDIM)
      double precision :: sol(2:nx+1,2:ny+1), rhs(2:nx+1,2:ny+1)
      double precision :: mat(0:4,2:nx+1,2:ny+1)
      integer :: iref, jref
      
      sol(:,:) = 0.0d0
      rhs(:,:) = 0.0d0
      mat(:,:,:) = 0.0d0
      
      ! PPE matrix and RHS
      do i = 2, nx+1
      do j = 2, ny+1
        rhs(i,j) = -1.0d0/dt * ( &
          (rh(i)*ustar(i,j) - rh(i-1)*ustar(i-1,j)) / dx + &
          (r(i)*vstar(i,j) - r(i)*vstar(i,j-1)) / dy)
        
        if ((bc_xlo.eq.BCTYPE_PER) .or. (i.ne.2)) then
          rho = 0.5d0 * (dens(i,j)+dens(i-1,j))
          coef = 1.0d0/(dx**2) * rh(i-1) / rho
          mat(0,i,j) = -coef
          mat(4,i,j) = mat(4,i,j) + coef
        endif
        if ((bc_xhi.eq.BCTYPE_PER) .or. (i.ne.nx+1)) then
          rho = 0.5d0 * (dens(i,j)+dens(i+1,j))
          coef = 1.0d0/(dx**2) * rh(i) / rho
          mat(2,i,j) = -coef
          mat(4,i,j) = mat(4,i,j) + coef
        endif
        if ((bc_ylo.eq.BCTYPE_PER) .or. (j.ne.2)) then
          rho = 0.5d0 * (dens(i,j)+dens(i,j-1))
          coef = 1.0d0/(dy**2) * r(i) / rho
          mat(1,i,j) = -coef
          mat(4,i,j) = mat(4,i,j) + coef
        endif
        if ((bc_yhi.eq.BCTYPE_PER) .or. (j.ne.ny+1)) then
          rho = 0.5d0 * (dens(i,j)+dens(i,j+1))
          coef = 1.0d0/(dy**2) * r(i) / rho
          mat(3,i,j) = -coef
          mat(4,i,j) = mat(4,i,j) + coef
        endif
      enddo
      enddo
      
      if (.true.) then
        iref = 2; jref = 2;
        mat(4,iref,jref) = 1.0d8
        !rhs(iref,jref) = 0.0d0
      endif
      
      periodic(:) = 0
      if ((bc_xlo.eq.BCTYPE_PER) .and. (bc_xhi.eq.BCTYPE_PER)) then
        periodic(1) = nx
      endif
      if ((bc_ylo.eq.BCTYPE_PER) .and. (bc_yhi.eq.BCTYPE_PER)) then
        periodic(2) = ny
      endif
      
      
      ! solve
      call solve_poisson2(mat, rhs, sol, nx,ny, periodic)
      
      ! get solution
      do i = 2, nx+1
      do j = 2, ny+1
        p(i,j) = sol(i,j)
      enddo
      enddo
      
      return
      end subroutine pres_solve2
      
      subroutine vel_corr( &
        ustar, vstar, dens, &
        u, v, p, &
        nx, ny, dx, dy, dt)
      use geom_module, only: iumax,jvmax
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: dx, dy, dt
      double precision, intent(in) :: ustar(nx+2,ny+2), vstar(nx+2,ny+2)
      double precision, intent(in) :: dens(1:nx+2,1:ny+2)
      double precision, intent(out) :: u(nx+2,ny+2), v(nx+2,ny+2)
      double precision, intent(in) :: p(1:nx+2,1:ny+2)
      !
      integer :: i, j
      double precision :: rp
      
      !do i = 2, nx
      do i = 2, iumax
      do j = 2, ny+1
        rp = 0.5d0 * (dens(i+1,j)+dens(i,j))
        u(i,j) = ustar(i,j) - dt/dx * (p(i+1,j)-p(i,j)) / rp
      enddo
      enddo
      
      do i = 2, nx+1
      !do j = 2, ny
      do j = 2, jvmax
        rp = 0.5d0 * (dens(i,j+1)+dens(i,j))
        v(i,j) = vstar(i,j) - dt/dy * (p(i,j+1)-p(i,j)) / rp
      enddo
      enddo
      
      return
      end subroutine vel_corr
      
      


