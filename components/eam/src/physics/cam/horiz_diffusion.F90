!---------------------------------------------------------------------
! compute horizontal turbulent diffusion
!---------------------------------------------------------------------
module horiz_diffusion

  use constituents,   only: pcnst
  use element_mod,    only: element_t     
  use dimensions_mod, only: np, npsq, nelemd, nlev
  use dof_mod,        only: UniquePoints, PutUniquePoints
  use dyn_comp,       only: dyn_export_t, dyn_import_t, TimeLevel
  use dyn_grid,       only: get_gcol_block_d
  use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
  use physics_types,  only: physics_state, physics_tend
  use spmd_dyn,       only: block_buf_nrecs, chunk_buf_nrecs
  use spmd_utils,     only: mpicom, iam
  use parallel_mod,   only: par
  use kinds,          only: real_kind, int_kind
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use perf_mod,       only: t_startf, t_stopf, t_barrierf  
  use phys_grid,      only: get_ncols_p, get_gcol_all_p, &
                            transpose_block_to_chunk, transpose_chunk_to_block,   &
                            chunk_to_block_send_pters, chunk_to_block_recv_pters, &
                            block_to_chunk_recv_pters, block_to_chunk_send_pters  
  
CONTAINS
  !===================================================================
  !===================================================================
  ! Top level subroutine for horizontal diffusion
  
  subroutine turbulent_horiz_diffusion(phys_state, dyn_out, phys_tend)

    implicit none
  
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    
    type(dyn_export_t), intent(inout)  :: dyn_out ! dynamics export
    
    ! LOCAL VARIABLES
    type(element_t), pointer :: elem(:)           ! pointer to dyn_out element array

    real (kind=real_kind), dimension(npsq,pver,nelemd)       :: T_phys  ! temp array to hold T
    real (kind=real_kind), dimension(npsq,2,pver,nelemd)     :: uv_phys ! temp array to hold u and v
    real (kind=real_kind), dimension(npsq,pver,nelemd) :: q_phys  ! temp to hold advected constituents  
    real (kind=real_kind), dimension(np,np,pver,nelemd) :: T_dyn, Q_dyn
    real (kind=real_kind), dimension(np,np,2,pver,nelemd) :: uv_dyn  
    integer                  :: cpter(pcols,0:pver)   ! offsets into chunk buffer for packing 
    integer                  :: bpter(npsq,0:pver)    ! offsets into block buffer for unpacking
    integer                  :: tsize, lchnk, icol, ilyr, ncols
    integer                  :: i, j, m, ie
    integer                  :: nphys, nphys_sq
    integer                  :: tl_f
    ! Transpose buffers
    real (kind=real_kind), allocatable, dimension(:) :: bbuffer 
    real (kind=real_kind), allocatable, dimension(:) :: cbuffer     
    
    !---------------------------------------------------------------------------   
    
    if (par%dynproc) then
      elem => dyn_out%elem
      tl_f = TimeLevel%n0
    else
      nullify(elem)
    endif
    
    nphys = np
    nphys_sq = nphys*nphys
    
    T_phys  = 0.0_r8
    uv_phys = 0.0_r8
    q_phys  = 0.0_r8    

    ! for now horizontal diffusion will only operate on T, Q, U, V
    !   (i.e. no constiuents at this point in time)
    tsize = 4

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  We need to go from physics to dynamics space
    !!!  Thus this section of code follows that p_d_coupling

    allocate( bbuffer(tsize*block_buf_nrecs) )
    allocate( cbuffer(tsize*chunk_buf_nrecs) )

    !$omp parallel do private (lchnk, ncols, cpter, i, icol, ilyr, m)   
    !- First transform physics variables into dynamics variables
    do lchnk = begchunk,endchunk
      ncols = get_ncols_p(lchnk)
      call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)
      do i = 1,ncols
        cbuffer(cpter(i,0):cpter(i,0)+2+pcnst) = 0.0_r8
      end do
      do icol = 1,ncols
        do ilyr = 1,pver
          cbuffer(cpter(icol,ilyr))   = phys_state(lchnk)%t(icol,ilyr)
          cbuffer(cpter(icol,ilyr)+1) = phys_state(lchnk)%u(icol,ilyr)
          cbuffer(cpter(icol,ilyr)+2) = phys_state(lchnk)%v(icol,ilyr)
          cbuffer(cpter(icol,ilyr)+3) = phys_state(lchnk)%q(icol,ilyr,1)
        end do ! ilyr
      end do ! icol
    end do ! lchnk

    call t_barrierf('sync_chk_to_blk', mpicom)
    call t_startf ('chunk_to_block')
    call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
    call t_stopf  ('chunk_to_block')

    if (par%dynproc) then
      !$omp parallel do private (ie, bpter, icol, ilyr, m, ncols)
      do ie = 1,nelemd
        ncols = elem(ie)%idxP%NumUniquePts
        call chunk_to_block_recv_pters(elem(ie)%GlobalID,nphys_sq,pver+1,tsize,bpter(1:nphys_sq,:))
        do icol = 1,ncols
          do ilyr = 1,pver
            T_phys  (icol,ilyr,ie)   = bbuffer(bpter(icol,ilyr))
            uv_phys (icol,1,ilyr,ie) = bbuffer(bpter(icol,ilyr)+1)
            uv_phys (icol,2,ilyr,ie) = bbuffer(bpter(icol,ilyr)+2)
            q_phys(icol,ilyr,ie) = bbuffer(bpter(icol,ilyr)+3)
          end do ! ilyr
        end do ! icol
      end do ! ie
    end if ! par%dynproc

    deallocate( bbuffer )
    deallocate( cbuffer )
    
    if (par%dynproc) then
      call t_startf('putUniquePoints')
      do ie = 1,nelemd
        ncols = elem(ie)%idxP%NumUniquePts
        ! Will need to convert this to vtheta_dp
        call putUniquePoints(elem(ie)%idxP,nlev,T_phys(1:ncols,:,ie),T_dyn(:,:,:,ie))
        call putUniquePoints(elem(ie)%idxP,2,nlev,uv_phys(1:ncols,:,:,ie),uv_dyn(:,:,:,:,ie))
        call putUniquePoints(elem(ie)%idxP,nlev,q_phys(1:ncols,:,ie),Q_dyn(:,:,:,ie))
      end do ! ie
      call t_stopf('putUniquePoints')
    end if ! par%dynproc
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !!!  Main part of the horizontal diffusion process goes here
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !!!  Transfer back from dynamics to physics space
    !!!  Section of code follows that of d_p_coupling
    
    if ( par%dynproc) then 
      
      call t_startf('UniquePoints')
      do ie = 1,nelemd
        ncols = elem(ie)%idxP%NumUniquePts
        call UniquePoints(elem(ie)%idxP,nlev,T_dyn(:,:,:,ie),T_phys(1:ncols,:,ie))
        call UniquePoints(elem(ie)%idxP,2,nlev,uv_dyn(:,:,:,:,ie),uv_phys(1:ncols,:,:,ie))
        call UniquePoints(elem(ie)%idxP,nlev,Q_dyn(:,:,:,ie), Q_phys(1:ncols,:,ie))
      end do
      call t_stopf('UniquePoints')
    
    else
    
      T_phys(:,:,:) = 0._r8
      uv_phys(:,:,:,:) = 0._r8
      q_phys(:,:,:) = 0._r8
    
    endif
    
    ! Now map back to physics
    
    tsize = 4
    allocate(bbuffer(tsize*block_buf_nrecs))
    allocate(cbuffer(tsize*chunk_buf_nrecs)) 
    
    if (par%dynproc) then

      !$omp parallel do private (ie, bpter, icol, ilyr, m, ncols)
      do ie = 1,nelemd
        call block_to_chunk_send_pters(elem(ie)%GlobalID,nphys_sq,pver+1,tsize,bpter(1:nphys_sq,:))
        ncols = elem(ie)%idxP%NumUniquePts
        do icol = 1,ncols
          do ilyr = 1,pver
            bbuffer(bpter(icol,ilyr))   = T_phys(icol,ilyr,ie)
            bbuffer(bpter(icol,ilyr)+1) = uv_phys(icol,1,ilyr,ie)
            bbuffer(bpter(icol,ilyr)+2) = uv_phys(icol,2,ilyr,ie)
            bbuffer(bpter(icol,ilyr)+3) = Q_phys(icol,ilyr,ie)
          end do ! ilyr
        end do ! icol
      end do ! ie
    
    else
      bbuffer(:) = 0._r8
    endif 
    
  end subroutine turbulent_horiz_diffusion
  
end module horiz_diffusion
