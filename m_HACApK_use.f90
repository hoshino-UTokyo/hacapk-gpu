!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 1.0.0                                             *
!                                                                     *
!   License                                                           *
!     This file is part of HACApK.                                    *
!     HACApK is a free software, you can use it under the terms       *
!     of The MIT License (MIT). See LICENSE file and User's guide     *
!     for more details.                                               *
!                                                                     *
!   ppOpen-HPC project:                                               *
!     Open Source Infrastructure for Development and Execution of     *
!     Large-Scale Scientific Applications on Post-Peta-Scale          *
!     Supercomputers with Automatic Tuning (AT).                      *
!                                                                     *
!   Sponsorship:                                                      *
!     Japan Science and Technology Agency (JST), Basic Research       *
!     Programs: CREST, Development of System Software Technologies    *
!     for post-Peta Scale High Performance Computing.                 *
!                                                                     *
!   Copyright (c) 2015 <Akihiro Ida and Takeshi Iwashita>             *
!                                                                     *
!=====================================================================*
!C**************************************************************************
!C  This file includes examples of integrating routines for H-matrices
!C  created by Akihiro Ida at Kyoto University on May 2012
!C  last modified by Akihiro Ida on Sep. 2014
!C**************************************************************************
module m_HACApK_use
 use m_HACApK_solve
 use m_HACApK_base
 implicit real*8(a-h,o-z)
contains

!*** HACApK_gensolv
 integer function HACApK_gensolv(st_leafmtxp,st_bemv,st_ctl,gmid,rhs,sol,ztol)
#ifdef _OPENACC
   use openacc
#endif
   implicit none
   include 'mpif.h'
   type(st_HACApK_leafmtxp) :: st_leafmtxp
   type(st_HACApK_leafmtx_acc) :: st_lf_acc
   type(st_HACApK_lcontrol) :: st_ctl
   type(st_HACApK_calc_entry) :: st_bemv
   real*8 :: gmid(st_bemv%nd,3),rhs(st_bemv%nd),sol(st_bemv%nd),ztol
   real*8 :: ww(st_bemv%nd),aww(st_bemv%nd),aww1(st_bemv%nd)
   integer*4, dimension(:),pointer :: ndl,ndt,nstrtl,nstrtt,kt
   integer*4, dimension(:),pointer :: idx1,idx2,idx3,orgid
   real*8, dimension(:),pointer :: a1,a2,a3
   integer*4,dimension(:),allocatable :: lnmtx
   integer*4 :: mpinr, mpilog, nrank, icomm, nthr, lrtrn, ierr, nstp
   real*8 :: st_s, st_create_hmtx, st_measure_time_ax, en_measure_time_ax,st1,en1,st2,en2,st3,en3
#ifdef _OPENACC
   integer(acc_device_kind) :: devType
   integer :: numOfDev, devNum
#endif
   
1000 format(5(a,i10)/)
2000 format(5(a,1pe15.8)/)
!   st2 = MPI_Wtime()

   mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); icomm=st_ctl%lpmd(1); nthr=st_ctl%lpmd(20)
   icomm=st_ctl%lpmd(1)

#ifdef _OPENACC
   devType = acc_get_device_type()
   call acc_init(devType)
   numOfDev = acc_get_num_devices(devType)
!   print *, mpinr, numOfDev
   call acc_set_device_num(mod(mpinr, numOfDev),devType)
!   print *, mpinr, acc_get_device_num(devType)
#endif
   call MPI_Barrier( icomm, ierr )

!   st1 = MPI_Wtime()
#ifndef _OPENACC
   lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,gmid,ztol)
#else
   allocate(lnmtx(3))
!   call acc_init(0)
   lrtrn=HACApK_generate_frame(st_leafmtxp,st_lf_acc,st_bemv,st_ctl,gmid,ztol,st_s,st_create_hmtx,lnmtx)
   !$acc data  &
!   !$acc create(st_lf_acc%a1) &
!   !$acc create(st_lf_acc%a2) &
!   !$acc create(st_lf_acc%a3) &
#if OPT==7
#else
   !$acc create(st_lf_acc%idx1) &
   !$acc create(st_lf_acc%idx2) &
#endif
   !$acc create(st_lf_acc%idx3) &
   !$acc create(st_lf_acc%kt) &
   !$acc copyin(st_lf_acc%nstrtl) &
   !$acc copyin(st_lf_acc%nstrtt) &
   !$acc copyin(st_lf_acc%ndl) &
   !$acc copyin(st_lf_acc%ndt) 
   lrtrn=HACApK_generate_acc(st_leafmtxp,st_lf_acc,st_bemv,st_ctl,gmid,ztol,st_s,st_create_hmtx,lnmtx,numOfDev)
   
   ! print *, mpinr, st_lf_acc%nlf, st_lf_acc%nlfkt
   ! stop
#endif
   call MPI_Barrier( icomm, ierr )
   ! en1 = MPI_Wtime()
   ! print *, "time gen", en1-st1
   ! st1 = MPI_Wtime()
#ifndef _OPENACC
   lrtrn=HACApK_solve(st_leafmtxp,st_bemv,st_ctl,rhs,sol,ztol)
#else
   lrtrn=HACApK_solve_acc(st_lf_acc,st_bemv,st_ctl,rhs,sol,ztol)
#endif
   call MPI_Barrier( icomm, ierr )
   ! en1 = MPI_Wtime()
   ! print *, "time solve", en1-st1
   st_measure_time_ax=MPI_Wtime()
#ifndef _OPENACC
   call HACApK_measurez_time_ax_lfmtx(st_leafmtxp,st_ctl,st_bemv%nd,nstp,lrtrn)
#else
   call HACApK_measurez_time_ax_lfmtx_acc(st_lf_acc,st_ctl,st_bemv%nd,nstp,lrtrn)
#endif
   en_measure_time_ax=MPI_Wtime()
   if(st_ctl%param(1)>0 .and. mpinr==0)  write(6,2000) 'lfmtx; time_AX_once  =',(en_measure_time_ax - st_measure_time_ax)/st_ctl%param(99)
#ifdef _OPENACC
   !$acc end data
   !$acc exit data &
#if OPT>=7
   !$acc delete(st_lf_acc%idx1) &
   !$acc delete(st_lf_acc%idx2) &
#endif
   !$acc delete(st_lf_acc%a1) &
   !$acc delete(st_lf_acc%a2) &
   !$acc delete(st_lf_acc%a3) 
#endif
   ! en2 = MPI_Wtime()
   ! print *, "time ax", en2-en1
   ! print *, "time gensolv", en2-st2
9999 continue
   HACApK_gensolv=lrtrn
 endfunction HACApK_gensolv
 
!*** HACApK_generate
 integer function HACApK_generate(st_leafmtxp,st_bemv,st_ctl,coord,ztol)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_calc_entry) :: st_bemv
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: coord(st_bemv%nd,*)
 integer*8 :: mem8,nth1_mem,imem
 integer*4 :: ierr
 integer*4,dimension(:),allocatable :: lnmtx(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,1pe15.8)/)
 
 lrtrn=0
 nofc=st_bemv%nd; nffc=1; ndim=3
 mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); icomm=st_ctl%lpmd(1); nthr=st_ctl%lpmd(20)
 st_ctl%param(71)=ztol
 
 call HACApK_chk_st_ctl(st_ctl)
 
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'***************** HACApK start ********************'
 if(st_ctl%param(1)>0)  write(mpilog,1000) 'irank=',mpinr,', nrank=',nrank
 nd=nofc*nffc
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'nd=',nd,' nofc=',nofc,' nffc=',nffc
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'nrank=',nrank,' nth=',nthr
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'param:'
 if(st_ctl%param(1)>1 .and. mpinr==0) write(*,10000) st_ctl%param(1:100)
 10000 format(10(1pe10.3)/)
 allocate(lnmtx(3))
 call MPI_Barrier( icomm, ierr )
 st_s=MPI_Wtime()

 call HACApK_generate_frame_leafmtx(st_leafmtxp,st_bemv,st_ctl,coord,lnmtx,nofc,nffc,ndim)

 call MPI_Barrier( icomm, ierr )
 st_create_hmtx=MPI_Wtime()
 st_bemv%lp61=0
 if(st_ctl%param(61)==2)then
   call HACApK_cal_matnorm(znrm2,st_bemv,st_ctl%lpmd,nd)
   call MPI_Barrier( icomm, ierr )
   call MPI_Allreduce( znrm2, znrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr );
   znrm=dsqrt(znrm)/nd
!   print*,'irank=',mpinr,'znrm2=',znrm2,' znrm=',znrm
 elseif(st_ctl%param(61)==3)then
   ndnr_s=st_ctl%lpmd(6); ndnr_e=st_ctl%lpmd(7); ndnr=st_ctl%lpmd(5)
   allocate(st_bemv%ao(nd)); st_bemv%ao(:)=0.0d0; zsqnd=sqrt(real(nd))
   do il=ndnr_s,ndnr_e
     zad=HACApK_entry_ij(il,il,st_bemv)
     st_bemv%ao(il)=1.0d0/dsqrt(zad/zsqnd)
   enddo
   call MPI_Barrier( icomm, ierr )
   call HACApK_impi_allgv(st_bemv%ao,st_ctl%lpmd,nd)
!   call MPI_Barrier( icomm, ierr )
   znrm=1.0/nd
   st_bemv%lp61=3
 else
   znrm=0.0d0
 endif
 call MPI_Barrier( icomm, ierr )
 st_cal_matnorm=MPI_Wtime()
 if(st_ctl%param(1)>0)  write(mpilog,1000) 'ndnr_s=',st_ctl%lpmd(6),', ndnr_e=',st_ctl%lpmd(7),', ndnr=',st_ctl%lpmd(5)
 if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,' ndlf_s=',st_ctl%lpmd(11),', ndlf_e=',st_ctl%lpmd(12),', ndlf=',st_leafmtxp%nlf
 lnps=nd+1; lnpe=0
 if(st_leafmtxp%nlf<1)then
   print*,'ERROR!; sub HACApK_generate; irank=',mpinr,' nlf=',st_leafmtxp%nlf
 endif
 if(st_ctl%param(10)==0) return
 call HACApK_fill_leafmtx_hyp(st_leafmtxp%st_lf,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_leafmtxp%nlf,lnps,lnpe,st_ctl%lthr)
! call HACApK_fill_leafmtx(st_leafmtxp%st_lf,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_leafmtxp%nlf,lnps,lnpe)
 call MPI_Barrier( icomm, ierr )
 if(st_ctl%param(1)>0)  write(mpilog,*) 'No. of nsmtx',lnmtx(1:3)
 ndnr_s=st_ctl%lpmd(6); ndnr_e=st_ctl%lpmd(7); ndnr=st_ctl%lpmd(5)

 st_fill_hmtx=MPI_Wtime()
 if(st_ctl%param(1)>0)  write(mpilog,2000)  'time_supermatrix             =',st_create_hmtx- st_s
 if(st_ctl%param(1)>0)  write(mpilog,2000)  'time_fill_hmtx               =',st_fill_hmtx-st_cal_matnorm
 if(st_ctl%param(1)>0)  write(mpilog,2000)  'time_construction_Hmatrix    =',st_fill_hmtx-st_s

 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_supermatrix             =',st_create_hmtx - st_s
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_fill_hmtx               =',st_fill_hmtx - st_cal_matnorm
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_construction_Hmatrix    =',st_fill_hmtx - st_s

 call MPI_Barrier( icomm, ierr )

 call HACApK_chk_leafmtx(st_leafmtxp,st_ctl,lnmtx,nd,mem8)

 ktp=0
 call HACApK_setcutthread(st_ctl%lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp)
      
 call MPI_Barrier( icomm, ierr )
! print*,'mpinr=',mpinr,lnps,lnpe
 st_ctl%lnp(mpinr+1)=lnpe-lnps
 call MPI_Barrier( icomm, ierr )
 call MPI_Allgather(lnpe-lnps,1,MPI_INTEGER,st_ctl%lnp,1, MPI_INTEGER, icomm, ierr )
 st_ctl%lsp(mpinr+1)=lnps
 call MPI_Allgather(lnps,1,MPI_INTEGER,st_ctl%lsp,1, MPI_INTEGER, icomm, ierr )
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'lnp=',st_ctl%lnp(:nrank)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'lsp=',st_ctl%lsp(:nrank)
 
 if(st_ctl%param(11)/=0) then
   call MPI_Barrier( icomm, ierr )
   call HACApK_accuracy_leafmtx(st_leafmtxp,st_bemv,st_ctl,st_ctl%lod,st_ctl%lod,st_ctl%lpmd,nofc,nffc)
 endif
9999 continue
 HACApK_generate=lrtrn
 endfunction

! #ifndef GPUMEMLIMIT      
! #define GPUMEMLIMIT 8000000000
! #endif
! #define LIMITSIZE ( GPUMEMLIMIT / 8 / 8 )
#define LIMITSIZE 800000000

 integer function HACApK_generate_frame(st_leafmtxp,st_lf_acc,st_bemv,st_ctl,coord,ztol,st_s,st_create_hmtx,lnmtx)
   include 'mpif.h'
   type(st_HACApK_leafmtxp) :: st_leafmtxp
   type(st_HACApK_leafmtx_acc) :: st_lf_acc
   type(st_HACApK_calc_entry) :: st_bemv
   type(st_HACApK_lcontrol) :: st_ctl
   real*8 :: st_s,st_create
   real*8 :: coord(st_bemv%nd,*)
   integer*8 :: mem8,nth1_mem,imem
   integer*4 :: ierr
!   integer*4,dimension(:),allocatable :: lnmtx(:)
   integer*4,dimension(:) :: lnmtx
   integer*4 :: nlf,n_low,a3size,ndt_max_dense
1000 format(5(a,i10)/)
2000 format(5(a,1pe15.8)/)
   
   lrtrn=0
   nofc=st_bemv%nd; nffc=1; ndim=3
   mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); icomm=st_ctl%lpmd(1); nthr=st_ctl%lpmd(20)
   st_ctl%param(71)=ztol
   
   call HACApK_chk_st_ctl(st_ctl)
   
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,'***************** HACApK start ********************'
   if(st_ctl%param(1)>0)  write(mpilog,1000) 'irank=',mpinr,', nrank=',nrank
   nd=nofc*nffc
   if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'nd=',nd,' nofc=',nofc,' nffc=',nffc
   if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'nrank=',nrank,' nth=',nthr
   if(st_ctl%param(1)>1 .and. mpinr==0) print*,'param:'
   if(st_ctl%param(1)>1 .and. mpinr==0) write(*,10000) st_ctl%param(1:100)
10000 format(10(1pe10.3)/)
!   allocate(lnmtx(3))
   call MPI_Barrier( icomm, ierr )
    st_s=MPI_Wtime()
   
   call HACApK_generate_frame_leafmtx(st_leafmtxp,st_bemv,st_ctl,coord,lnmtx,nofc,nffc,ndim)
   
   call MPI_Barrier( icomm, ierr )
   st_create_hmtx=MPI_Wtime()

   nlf = st_leafmtxp%nlf
   allocate(st_lf_acc%ltmtx(nlf))
   allocate(st_lf_acc%kt(nlf))
   allocate(st_lf_acc%nstrtl(nlf))
   allocate(st_lf_acc%ndl(nlf))
   allocate(st_lf_acc%nstrtt(nlf))
   allocate(st_lf_acc%ndt(nlf))
#if OPT>=7
#else
   allocate(st_lf_acc%idx1(nlf))
   allocate(st_lf_acc%idx2(nlf))
   allocate(st_lf_acc%idx3(nlf))
#endif
   allocate(st_lf_acc%orgid(nlf))

   st_lf_acc%nlf = nlf
   n_low = 0; n_dense = 0
   do ip=1, nlf
      if(st_leafmtxp%st_lf(ip)%ltmtx==1) then
         n_low = n_low + 1
      elseif(st_leafmtxp%st_lf(ip)%ltmtx==2) then
         n_dense = n_dense + 1
      else
         write(*,*) 'unknown ltmtx at ip = ', ip,__FILE__,'(',__LINE__,')'
      endif
   enddo

   st_lf_acc%nlfkt = n_low
   i = 1; j = n_low+1
   do ip=1, nlf
      if(st_leafmtxp%st_lf(ip)%ltmtx==1) then
         st_lf_acc%ltmtx(i)  = st_leafmtxp%st_lf(ip)%ltmtx
         st_lf_acc%nstrtl(i) = st_leafmtxp%st_lf(ip)%nstrtl
         st_lf_acc%ndl(i)    = st_leafmtxp%st_lf(ip)%ndl
         st_lf_acc%nstrtt(i) = st_leafmtxp%st_lf(ip)%nstrtt
         st_lf_acc%ndt(i)    = st_leafmtxp%st_lf(ip)%ndt
         st_lf_acc%orgid(i)  = ip
         i = i + 1
      elseif(st_leafmtxp%st_lf(ip)%ltmtx==2) then
         st_lf_acc%ltmtx(j)  = st_leafmtxp%st_lf(ip)%ltmtx
         st_lf_acc%nstrtl(j) = st_leafmtxp%st_lf(ip)%nstrtl
         st_lf_acc%ndl(j)    = st_leafmtxp%st_lf(ip)%ndl
         st_lf_acc%nstrtt(j) = st_leafmtxp%st_lf(ip)%nstrtt
         st_lf_acc%ndt(j)    = st_leafmtxp%st_lf(ip)%ndt
         st_lf_acc%orgid(j)  = ip
         j = j + 1
      else
         write(*,*) 'unknown ltmtx at ip = ', ip,__FILE__,'(',__LINE__,')'
      endif
   enddo


#if OPT>=3
#if OPT<6
   call quick_sort(st_lf_acc%nstrtl,st_lf_acc%ndl,st_lf_acc%nstrtt,st_lf_acc%ndt,st_lf_acc%orgid,1,n_low)
#endif
! #if OPT>=7
!    call quick_sort(st_lf_acc%nstrtl,st_lf_acc%ndl,st_lf_acc%nstrtt,st_lf_acc%ndt,st_lf_acc%orgid,n_low+1,n_dense)
! #endif
#endif

   a3size= 0
   ndt_max_dense= 0
   do ip=n_low+1, n_low+n_dense
      a3size= a3size + st_lf_acc%ndl(ip) * st_lf_acc%ndt(ip)
      ndt_max_dense= max(ndt_max_dense,st_lf_acc%ndt(ip))
   enddo
   st_lf_acc%ndt_max_dense = ndt_max_dense

!   print *, "a3size", a3size
   allocate(st_lf_acc%a1(1:LIMITSIZE))
   allocate(st_lf_acc%a2(1:LIMITSIZE))
   allocate(st_lf_acc%a3(1:a3size))
   allocate(st_lf_acc%idx3(n_dense*ndt_max_dense+1))

   HACApK_generate_frame=lrtrn

 endfunction HACApK_generate_frame

 integer function HACApK_generate_acc(st_leafmtxp,st_lf_acc,st_bemv,st_ctl,coord,ztol,st_s,st_create_hmtx,lnmtx,numOfDev)
   implicit none
   include 'mpif.h'
   type(st_HACApK_leafmtxp), intent(inout) :: st_leafmtxp
   type(st_HACApK_leafmtx_acc), intent(inout) :: st_lf_acc
   type(st_HACApK_calc_entry), intent(inout) :: st_bemv
   type(st_HACApK_lcontrol), intent(inout) :: st_ctl
   real*8, intent(in) :: coord(st_bemv%nd,*)
   real*8, intent(in) :: ztol,st_s,st_create_hmtx
   integer*4, intent(in), dimension(:) :: lnmtx
   integer*4, intent(in) :: numOfDev
   integer*8 :: mem8,nth1_mem,imem
   integer*4 :: ierr
   integer*4 :: lrtrn,nofc,nffc,ndim,mpinr,mpilog,nrank,icomm,nthr,nd,ndnr_s,ndnr_e,ndnr,il,lnps,lnpe,ktp
   real*8 :: znrm2,znrm,zsqnd,zad,st_cal_matnorm,st_fill_hmtx
1000 format(5(a,i10)/)
2000 format(5(a,1pe15.8)/)
   
   lrtrn=0
   nofc=st_bemv%nd; nffc=1; ndim=3
   mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); icomm=st_ctl%lpmd(1); nthr=st_ctl%lpmd(20)
   st_ctl%param(71)=ztol
   nd=nofc*nffc

   
   st_bemv%lp61=0
   if(st_ctl%param(61)==2)then
      call HACApK_cal_matnorm(znrm2,st_bemv,st_ctl%lpmd,nd)
      call MPI_Barrier( icomm, ierr )
      call MPI_Allreduce( znrm2, znrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr );
      znrm=dsqrt(znrm)/nd
      !   print*,'irank=',mpinr,'znrm2=',znrm2,' znrm=',znrm
   elseif(st_ctl%param(61)==3)then
      ndnr_s=st_ctl%lpmd(6); ndnr_e=st_ctl%lpmd(7); ndnr=st_ctl%lpmd(5)
      allocate(st_bemv%ao(nd)); st_bemv%ao(:)=0.0d0; zsqnd=sqrt(real(nd))
      do il=ndnr_s,ndnr_e
         zad=HACApK_entry_ij(il,il,st_bemv)
         st_bemv%ao(il)=1.0d0/dsqrt(zad/zsqnd)
      enddo
      call MPI_Barrier( icomm, ierr )
      call HACApK_impi_allgv(st_bemv%ao,st_ctl%lpmd,nd)
      !   call MPI_Barrier( icomm, ierr )
      znrm=1.0/nd
      st_bemv%lp61=3
   else
      znrm=0.0d0
   endif
   call MPI_Barrier( icomm, ierr )
   st_cal_matnorm=MPI_Wtime()
   if(st_ctl%param(1)>0)  write(mpilog,1000) 'ndnr_s=',st_ctl%lpmd(6),', ndnr_e=',st_ctl%lpmd(7),', ndnr=',st_ctl%lpmd(5)
   if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,' ndlf_s=',st_ctl%lpmd(11),', ndlf_e=',st_ctl%lpmd(12),', ndlf=',st_leafmtxp%nlf
   lnps=nd+1; lnpe=0
   if(st_leafmtxp%nlf<1)then
      print*,'ERROR!; sub HACApK_generate; irank=',mpinr,' nlf=',st_leafmtxp%nlf
   endif
   if(st_ctl%param(10)==0) return
#if OPT==6
!   call HACApK_fill_leafmtx_acc4(st_leafmtxp%st_lf,st_lf_acc,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_leafmtxp%nlf,lnps,lnpe,st_ctl%lthr,st_leafmtxp%nlfkt)
!   print *, mpinr,znrm,nd,st_leafmtxp%nlf,lnps,lnpe,st_ctl%lthr(0)
   call HACApK_fill_leafmtx_acc4(st_leafmtxp%st_lf,st_lf_acc,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_lf_acc%nlf,lnps,lnpe,st_ctl%lthr,st_lf_acc%nlfkt,numOfDev)
#elif OPT>=7
   call HACApK_fill_leafmtx_acc7(st_leafmtxp%st_lf,st_lf_acc,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_lf_acc%nlf,lnps,lnpe,st_ctl%lthr,st_lf_acc%nlfkt,numOfDev)
#else
   call HACApK_fill_leafmtx_acc3(st_leafmtxp%st_lf,st_lf_acc,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_leafmtxp%nlf,lnps,lnpe,st_ctl%lthr,st_leafmtxp%nlfkt)
#endif
   call MPI_Barrier( icomm, ierr )
   if(st_ctl%param(1)>0)  write(mpilog,*) 'No. of nsmtx',lnmtx(1:3)
   ndnr_s=st_ctl%lpmd(6); ndnr_e=st_ctl%lpmd(7); ndnr=st_ctl%lpmd(5)

   st_fill_hmtx=MPI_Wtime()
   if(st_ctl%param(1)>0)  write(mpilog,2000)  'time_supermatrix             =',st_create_hmtx- st_s
   if(st_ctl%param(1)>0)  write(mpilog,2000)  'time_fill_hmtx               =',st_fill_hmtx-st_cal_matnorm
   if(st_ctl%param(1)>0)  write(mpilog,2000)  'time_construction_Hmatrix    =',st_fill_hmtx-st_s
   
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_supermatrix             =',st_create_hmtx - st_s
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_fill_hmtx               =',st_fill_hmtx - st_cal_matnorm
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_construction_Hmatrix    =',st_fill_hmtx - st_s
   
   call MPI_Barrier( icomm, ierr )
   
!   call HACApK_chk_leafmtx(st_leafmtxp,st_ctl,lnmtx,nd,mem8)
   call HACApK_chk_leafmtx_acc(st_lf_acc,st_ctl,lnmtx,nd,mem8)
   
   ktp=0
!   call HACApK_setcutthread(st_ctl%lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp)
   
   call MPI_Barrier( icomm, ierr )
   ! print*,'mpinr=',mpinr,lnps,lnpe
   st_ctl%lnp(mpinr+1)=lnpe-lnps
   call MPI_Barrier( icomm, ierr )
   call MPI_Allgather(lnpe-lnps,1,MPI_INTEGER,st_ctl%lnp,1, MPI_INTEGER, icomm, ierr )
   st_ctl%lsp(mpinr+1)=lnps
   call MPI_Allgather(lnps,1,MPI_INTEGER,st_ctl%lsp,1, MPI_INTEGER, icomm, ierr )
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,'lnp=',st_ctl%lnp(:nrank)
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,'lsp=',st_ctl%lsp(:nrank)
   
   if(st_ctl%param(11)/=0) then
      !$acc update self(st_lf_acc%a1,st_lf_acc%a2,st_lf_acc%a3)
      call MPI_Barrier( icomm, ierr )
      call HACApK_accuracy_leafmtx_acc(st_lf_acc,st_bemv,st_ctl,st_ctl%lod,st_ctl%lod,st_ctl%lpmd,nofc,nffc)
   endif
9999 continue
   HACApK_generate_acc=lrtrn
 endfunction HACApK_generate_acc


!*** HACApK_solve
 integer function HACApK_solve(st_leafmtxp,st_bemv,st_ctl,rhs,sol,ztol)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: rhs(st_bemv%nd),sol(st_bemv%nd),ztol
 real*8,pointer :: param(:)
 real*8,dimension(:),allocatable :: u,b,www,ao
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,1pe15.8)/)

! print *, "solve start"
 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;lod => st_ctl%lod(:); param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1); nthr=lpmd(20)
 param(91)=ztol
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_solve start'
 nofc=st_bemv%nd;nffc=1;ndim=3
 nd=nofc*nffc
 if(st_ctl%param(1)>1) write(*,*) 'irank=',mpinr,' lthr=',lthr(0:nthr-1)
 allocate(u(nd),b(nd)); u(:nd)=sol(lod(:nd)); b(:nd)=rhs(lod(:nd))
 if(param(61)==3)then
!   do il=ndnr_s,ndnr_e
   do il=1,nd
     u(il)=u(il)/st_bemv%ao(lod(il))
     b(il)=b(il)*st_bemv%ao(lod(il))
   enddo
 endif
 if(param(83)>0) then
   allocate(ao(nd))
   do il=1,nd
     zzz=HACApK_entry_ij(il,il,st_bemv)
     ao(il)=1.0d0/zzz
   enddo
   
   call MPI_Barrier( icomm, ierr )
   st_measure_time_bicgstab=MPI_Wtime()
   if(param(85)==1)then
!      call HACApK_bicgstab_lfmtx(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
     call HACApK_bicgstab_lfmtx_hyp(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
   elseif(param(85)==2)then
     call HACApK_gcrm_lfmtx(st_leafmtxp,st_ctl,st_bemv,u,b,param,nd,nstp,lrtrn)
   else
   endif
   call MPI_Barrier( icomm, ierr )
   en_measure_time_bicgstab=MPI_Wtime()
   time_bicgstab = en_measure_time_bicgstab - st_measure_time_bicgstab
   if(st_ctl%param(1)>0 .and. mpinr==0)  write(6,2000)              'time_HACApK_solve  =',time_bicgstab
   if(st_ctl%param(1)>0 .and. mpinr==0 .and. nstp>1)  write(6,2000) '       time_1step  =',time_bicgstab/nstp
   allocate(www(nd))
   sol(:nd)=0.0d0; www(lod(:nd))=u(:nd); sol(:nd)=www(:nd)
   deallocate(www)
   if(param(61)==3)then
     do il=1,nd
       sol(il)=sol(il)*st_bemv%ao(il)
     enddo
   endif
 endif
9999 continue
 HACApK_solve=lrtrn
 endfunction


 integer function HACApK_solve_acc(st_lf_acc,st_bemv,st_ctl,rhs,sol,ztol)
   implicit none
   include 'mpif.h'
   type(st_HACApK_leafmtx_acc) :: st_lf_acc
   type(st_HACApK_lcontrol) :: st_ctl
   type(st_HACApK_calc_entry) :: st_bemv
   real*8 :: rhs(st_bemv%nd),sol(st_bemv%nd),ztol
   real*8,pointer :: param(:)
   real*8,dimension(:),allocatable :: u,b,www,ao
   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
   integer*4 :: mpinr,mpilog,nrank,icomm,nthr,nofc,nffc,ndim,nd,il,ierr,nstp,lrtrn
   real*8 :: zzz,st_measure_time_bicgstab,en_measure_time_bicgstab,time_bicgstab
1000 format(5(a,i10)/)
2000 format(5(a,1pe15.8)/)
   
   lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;lod => st_ctl%lod(:); param=>st_ctl%param(:)
   mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1); nthr=lpmd(20)
   param(91)=ztol
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_solve start'
   nofc=st_bemv%nd;nffc=1;ndim=3
   nd=nofc*nffc
   if(st_ctl%param(1)>1) write(*,*) 'irank=',mpinr,' lthr=',lthr(0:nthr-1)
   allocate(u(nd),b(nd)); u(:nd)=sol(lod(:nd)); b(:nd)=rhs(lod(:nd))
   if(param(61)==3)then
      !   do il=ndnr_s,ndnr_e
      do il=1,nd
         u(il)=u(il)/st_bemv%ao(lod(il))
         b(il)=b(il)*st_bemv%ao(lod(il))
      enddo
   endif
   if(param(83)>0) then
      allocate(ao(nd))
      do il=1,nd
         zzz=HACApK_entry_ij(il,il,st_bemv)
         ao(il)=1.0d0/zzz
      enddo
      
      call MPI_Barrier( icomm, ierr )
      st_measure_time_bicgstab=MPI_Wtime()
      if(param(85)==1)then
         !     call HACApK_bicgstab_lfmtx(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
         !     call HACApK_bicgstab_lfmtx_hyp(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
         !     call HACApK_bicgstab_lfmtx_acc2(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
#if GPUOPTIMIZED==0
         call HACApK_bicgstab_lfmtx_acc_unified(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
#else
         call HACApK_bicgstab_lfmtx_acc(st_lf_acc,st_ctl,u,b,param,nd,nstp,lrtrn)
!         call HACApK_bicgstab_lfmtx_cuda(st_lf_acc,st_ctl,u,b,param,nd,nstp,lrtrn)
#endif
      elseif(param(85)==2)then
         !     call HACApK_gcrm_lfmtx(st_leafmtxp,st_ctl,st_bemv,u,b,param,nd,nstp,lrtrn)
         print *, "gcrm not supported yet"
      else
      endif
      call MPI_Barrier( icomm, ierr )
      en_measure_time_bicgstab=MPI_Wtime()
      time_bicgstab = en_measure_time_bicgstab - st_measure_time_bicgstab
      if(st_ctl%param(1)>0 .and. mpinr==0)  write(6,2000)              'time_HACApK_solve  =',time_bicgstab
      if(st_ctl%param(1)>0 .and. mpinr==0 .and. nstp>1)  write(6,2000) '       time_1step  =',time_bicgstab/nstp
      allocate(www(nd))
      sol(:nd)=0.0d0; www(lod(:nd))=u(:nd); sol(:nd)=www(:nd)
      deallocate(www)
      if(param(61)==3)then
         do il=1,nd
            sol(il)=sol(il)*st_bemv%ao(il)
         enddo
      endif
   endif
9999 continue
   HACApK_solve_acc=lrtrn
 endfunction HACApK_solve_acc

endmodule m_HACApK_use
