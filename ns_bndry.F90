
#ifndef SDIM
#define SDIM 2
#endif

      subroutine bcfill_cell_homo( &
        phi, nx,ny)
      use const_module
      use geom_module, only: bc_xlo,bc_xhi,bc_ylo,bc_yhi
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(inout) :: phi(nx+2,ny+2)
      !integer, intent(in) :: bctype(SDIM,2)
      !
      if (bc_xlo == BCTYPE_PER) then
        phi(1,:) = phi(nx+1,:)
      else
        phi(1,:) = phi(2,:)
      endif
      if (bc_xhi == BCTYPE_PER) then
        phi(nx+2,:) = phi(2,:)
      else
        phi(nx+2,:) = phi(nx+1,:)
      endif
      if (bc_ylo == BCTYPE_PER) then
        phi(:,1) = phi(:,ny+1)
      else
        phi(:,1) = phi(:,2)
      endif
      if (bc_yhi == BCTYPE_PER) then
        phi(:,ny+2) = phi(:,2)
      else
        phi(:,ny+2) = phi(:,ny+1)
      endif
      return
      end subroutine bcfill_cell_homo
      
      
      subroutine bcfill_umac( &
        umac, nx,ny)
      use const_module
      use geom_module, only: bc_xlo,bc_xhi,bc_ylo,bc_yhi
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(inout) :: umac(nx+2,ny+2)
      !
      if (bc_xlo == BCTYPE_PER) then
        umac(1,2:ny+1) = umac(nx+1,2:ny+1)
      else
        umac(1,2:ny+1) = 0.0d0
      endif
      !
      if (bc_xhi == BCTYPE_PER) then
        umac(nx+2,2:ny+1) = umac(2,2:ny+1)
      else
        umac(nx+1,2:ny+1) = 0.0d0
      endif
      !
      if (bc_ylo == BCTYPE_PER) then
        umac(:,1) = umac(:,ny+1)
      else ! free-slip
        umac(:,1) = umac(:,2)
      endif
      !
      if (bc_yhi == BCTYPE_PER) then
        umac(:,ny+2) = umac(:,2)
      else ! free-slip
        umac(:,ny+2) = umac(:,ny+1)
      endif
      return
      end subroutine bcfill_umac
      
      subroutine bcfill_vmac( &
        vmac, nx,ny)
      use const_module
      use geom_module, only: bc_xlo,bc_xhi,bc_ylo,bc_yhi
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(inout) :: vmac(nx+2,ny+2)
      !
      if (bc_ylo == BCTYPE_PER) then
        vmac(2:nx+1,1) = vmac(2:nx+1,ny+1)
      else 
        vmac(2:nx+1,1) = 0.0d0
      endif
      !
      if (bc_yhi == BCTYPE_PER) then
        vmac(2:nx+1,ny+2) = vmac(2:nx+1,2)
      else
        vmac(2:nx+1,ny+1) = 0.0d0
      endif
      !
      if (bc_xlo == BCTYPE_PER) then
        vmac(1,:) = vmac(nx+1,:)
      else ! free-slip
        vmac(1,:) = vmac(2,:)
      endif
      !
      if (bc_xhi == BCTYPE_PER) then
        vmac(nx+2,:) = vmac(2,:)
      else
        vmac(nx+2,:) = vmac(nx+1,:)
      endif
      return
      end subroutine bcfill_vmac





