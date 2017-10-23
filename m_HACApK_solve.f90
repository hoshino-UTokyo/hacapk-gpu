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
!C***********************************************************************
!C  This file includes routines for utilizing H-matrices, such as solving
!C  linear system with an H-matrix as the coefficient matrix and 
!C  multiplying an H-matrix and a vector,
!C  created by Akihiro Ida at Kyoto University on May 2012,
!C  last modified by Akihiro Ida on Sep 2014,
!C***********************************************************************
module m_HACApK_solve
 use m_HACApK_base
 implicit real*8(a-h,o-z)
 implicit integer*4(i-n)
contains

!***HACApK_adot_lfmtx_p
 subroutine HACApK_adot_lfmtx_p(zau,st_leafmtxp,st_ctl,zu,nd)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(nd),zu(nd)
 real*8,dimension(:),allocatable :: wws,wwr
 integer*4 :: ISTATUS(MPI_STATUS_SIZE),isct(2),irct(2)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)
 allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
 zau(:)=0.0d0
 call HACApK_adot_body_lfmtx(zau,st_leafmtxp,st_ctl,zu,nd)
 if(nrank==1) return
 wws(1:lnp(mpinr))=zau(lsp(mpinr):lsp(mpinr)+lnp(mpinr)-1)
 ncdp=mod(mpinr+1,nrank)
 ncsp=mod(mpinr+nrank-1,nrank)
! write(mpilog,1000) 'destination process=',ncdp,'; source process=',ncsp
 isct(1)=lnp(mpinr);isct(2)=lsp(mpinr); 
! irct=lnp(ncsp)
 do ic=1,nrank-1
!   idp=mod(mpinr+ic,nrank) ! rank of destination process
!   isp=mod(mpinr+nrank+ic-2,nrank) ! rank of source process
   call MPI_SENDRECV(isct,2,MPI_INTEGER,ncdp,1, &
                     irct,2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
!   write(mpilog,1000) 'ISTATUS=',ISTATUS,'; ierr=',ierr
!   write(mpilog,1000) 'ic=',ic,'; isct=',isct(1),'; irct=',irct(1),'; ivsps=',isct(2),'; ivspr=',irct(2)

   call MPI_SENDRECV(wws,isct,MPI_DOUBLE_PRECISION,ncdp,1, &
                     wwr,irct,MPI_DOUBLE_PRECISION,ncsp,1,icomm,ISTATUS,ierr)
!   write(mpilog,1000) 'ISTATUS=',ISTATUS,'; ierr=',ierr
   
   zau(irct(2):irct(2)+irct(1)-1)=zau(irct(2):irct(2)+irct(1)-1)+wwr(:irct(1))
   wws(:irct(1))=wwr(:irct(1))
   isct=irct
!   write(mpilog,1000) 'ic=',ic,'; isct=',isct
 enddo
 deallocate(wws,wwr)
 end subroutine HACApK_adot_lfmtx_p
 
!***HACApK_adot_lfmtx_hyp
 subroutine HACApK_adot_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(*),zu(*),wws(*),wwr(*)
 integer*4 :: isct(*),irct(*)
 integer*4 :: ISTATUS(MPI_STATUS_SIZE)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
! integer*4,dimension(:),allocatable :: ISTATUS
! 1000 format(5(a,i10)/)
! 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
! allocate(ISTATUS(MPI_STATUS_SIZE))
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)
 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 if(nrank>1)then
   wws(1:lnp(mpinr))=zau(lsp(mpinr):lsp(mpinr)+lnp(mpinr)-1)
   ncdp=mod(mpinr+1,nrank)
   ncsp=mod(mpinr+nrank-1,nrank)
   isct(1)=lnp(mpinr);isct(2)=lsp(mpinr); 
   do ic=1,nrank-1
     call MPI_SENDRECV(isct,2,MPI_INTEGER,ncdp,1, &
                       irct,2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
     call MPI_SENDRECV(wws,isct,MPI_DOUBLE_PRECISION,ncdp,1, &
                       wwr,irct,MPI_DOUBLE_PRECISION,ncsp,1,icomm,ISTATUS,ierr)
     zau(irct(2):irct(2)+irct(1)-1)=zau(irct(2):irct(2)+irct(1)-1)+wwr(:irct(1))
     wws(:irct(1))=wwr(:irct(1))
     isct(:2)=irct(:2)
   enddo
 endif
!$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp


 subroutine HACApK_adot_lfmtx_acc_unified(zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
   implicit none
   include 'mpif.h'
   type(st_HACApK_leafmtxp) :: st_leafmtxp
   type(st_HACApK_lcontrol) :: st_ctl
   real*8 :: zau(*),zu(*),wws(*),wwr(*)
   integer*4 :: isct(*),irct(*)
   integer,intent(in) :: nd
   integer*4 :: ISTATUS(MPI_STATUS_SIZE)
   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
   integer*4 :: mpinr,mpilog,nrank,icomm,ndnr_s,ndnr_e,ndnr,ncdp,ncsp,ic,ierr
   
   lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
   mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
   ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)
   !$acc kernels
   zau(:nd)=0.0d0
   !$acc end kernels
   call HACApK_adot_body_lfmtx_acc_unified(zau,st_leafmtxp,st_ctl,zu,nd)
   
   if(nrank>1)then
      !$acc kernels
      wws(1:lnp(mpinr))=zau(lsp(mpinr):lsp(mpinr)+lnp(mpinr)-1)
      !$acc end kernels
      ncdp=mod(mpinr+1,nrank)
      ncsp=mod(mpinr+nrank-1,nrank)
      isct(1)=lnp(mpinr);isct(2)=lsp(mpinr); 
      do ic=1,nrank-1
         call MPI_SENDRECV(isct,2,MPI_INTEGER,ncdp,1, &
              irct,2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
         call MPI_SENDRECV(wws,isct,MPI_DOUBLE_PRECISION,ncdp,1, &
              wwr,irct,MPI_DOUBLE_PRECISION,ncsp,1,icomm,ISTATUS,ierr)
         
         !$acc kernels
         zau(irct(2):irct(2)+irct(1)-1)=zau(irct(2):irct(2)+irct(1)-1)+wwr(:irct(1))
         !$acc end kernels
         !$acc kernels
         wws(:irct(1))=wwr(:irct(1))
         !$acc end kernels
         isct(:2)=irct(:2)
      enddo
   endif
   
 end subroutine HACApK_adot_lfmtx_acc_unified
 
!***HACApK_adot_lfmtx_acc
 subroutine HACApK_adot_lfmtx_acc(zau,st_lf_acc,st_ctl,zu,wws,wwr,isct,irct,nd)
   implicit none
   include 'mpif.h'
   type(st_HACApK_leafmtx_acc) :: st_lf_acc
   type(st_HACApK_lcontrol) :: st_ctl
   real*8 :: zau(*),zu(*),wws(*),wwr(*)
   integer*4 :: isct(*),irct(*)
   integer*4, intent(in) :: nd
   integer*4 :: ISTATUS(MPI_STATUS_SIZE)
   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
   integer*4 :: mpinr,mpilog,nrank,icomm,ndnr_s,ndnr_e,ndnr,wwsize,ncdp,ncsp,ic,ierr
   ! integer*4,dimension(:),allocatable :: ISTATUS
   ! 1000 format(5(a,i10)/)
   ! 2000 format(5(a,f10.4)/)
   
   lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
   ! allocate(ISTATUS(MPI_STATUS_SIZE))
   mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
   ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)
   
   wwsize = maxval(lnp(0:nrank-1))
   
   !$acc data &
   !$acc copyin(zu(1:nd)) &
   !$acc copyout(zau(1:nd)) &
   !$acc copy(wws(1:wwsize),wwr(1:wwsize))

   !$acc kernels
   zau(:nd)=0.0d0
   !$acc end kernels
   
   ! call HACApK_adot_body_lfmtx_acc(zau,st_leafmtxp,st_ctl,zu,nd,n_low,n_dense,low_or_dense)
#ifdef USECUDA
   call HACApK_adot_body_lfmtx_cuda(zau,st_lf_acc,st_ctl,zu,nd)
#else
   call HACApK_adot_body_lfmtx_acc(zau,st_lf_acc,st_ctl,zu,nd)
#endif
   
   if(nrank>1)then
      !$acc kernels
      wws(1:lnp(mpinr))=zau(lsp(mpinr):lsp(mpinr)+lnp(mpinr)-1)
      !$acc end kernels
      !$acc update host(wws(1:lnp(mpinr)))
      ncdp=mod(mpinr+1,nrank)
      ncsp=mod(mpinr+nrank-1,nrank)
      isct(1)=lnp(mpinr);isct(2)=lsp(mpinr); 

      do ic=1,nrank-1
         call MPI_SENDRECV(isct,2,MPI_INTEGER,ncdp,1, &
              irct,2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
         !       !$acc host_data use_device(wws,wwr)
         call MPI_SENDRECV(wws,isct,MPI_DOUBLE_PRECISION,ncdp,1, &
              wwr,irct,MPI_DOUBLE_PRECISION,ncsp,1,icomm,ISTATUS,ierr)
         !       !$acc end host_data

         !$acc update device(wwr(1:irct(1)))
         !$acc kernels
         zau(irct(2):irct(2)+irct(1)-1)=zau(irct(2):irct(2)+irct(1)-1)+wwr(:irct(1))
         !$acc end kernels
         !$acc kernels
         wws(:irct(1))=wwr(:irct(1))
         !$acc end kernels
         !$acc update host(wws(1:irct(1)))
         isct(:2)=irct(:2)
      enddo
   endif
   
   !$acc end data
   
 end subroutine HACApK_adot_lfmtx_acc
 
!***HACApK_adot_body_lfmtx
 RECURSIVE subroutine HACApK_adot_body_lfmtx(zau,st_leafmtxp,st_ctl,zu,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(nd),zu(nd)
 real*8,dimension(:),allocatable :: zbu
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 nlf=st_leafmtxp%nlf
 do ip=1,nlf
   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     allocate(zbu(kt)); zbu(:)=0.0d0
     do il=1,kt
       do it=1,ndt; itt=it+nstrtt-1
         zbu(il)=zbu(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo
     do il=1,kt
       do it=1,ndl; ill=it+nstrtl-1
         zau(ill)=zau(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbu(il)
       enddo
     enddo
     deallocate(zbu)
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
         zau(ill)=zau(ill)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo
   endif
 enddo
 end subroutine HACApK_adot_body_lfmtx

!***HACApK_adot_body_lfmtx_hyp
 subroutine HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(*),zu(*)
 real*8,dimension(:),allocatable :: zbut
 real*8,dimension(:),allocatable :: zaut
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),ltmp(:)
 real*8,dimension(:,:),pointer :: a1,a2
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 ! lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;ltmp(0:) => st_ctl%lthr
 lpmd => st_ctl%lpmd(:);ltmp(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 nlf=st_leafmtxp%nlf; ktmax=st_leafmtxp%ktmax
 ith = omp_get_thread_num()
 ith1 = ith+1
 nths=ltmp(ith); nthe=ltmp(ith1)-1
 allocate(zaut(nd)); zaut(:)=0.0d0
 allocate(zbut(ktmax)) 
 ls=nd; le=1
 do ip=nths,nthe
   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
   a1 => st_leafmtxp%st_lf(ip)%a1
   a2 => st_leafmtxp%st_lf(ip)%a2
   if(nstrtl<ls) ls=nstrtl; if(nstrtl+ndl-1>le) le=nstrtl+ndl-1
   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     zbut(1:kt)=0.0d0
     do il=1,kt
       do it=1,ndt; itt=it+nstrtt-1
!         zbut(il)=zbut(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
         zbut(il)=zbut(il)+a1(it,il)*zu(itt)
       enddo
     enddo
     do il=1,kt
       do it=1,ndl; ill=it+nstrtl-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)
         zaut(ill)=zaut(ill)+a2(it,il)*zbut(il)
       enddo
     enddo
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
         zaut(ill)=zaut(ill)+a1(it,il)*zu(itt)
       enddo
     enddo
   endif
 enddo
 deallocate(zbut)
 
 do il=ls,le
!$omp atomic
   zau(il)=zau(il)+zaut(il)
 enddo
 end subroutine HACApK_adot_body_lfmtx_hyp

!***HACApK_adot_body_lfmtx_acc
#ifdef USECUDA
 subroutine HACApK_adot_body_lfmtx_cuda(zau,st_lf_acc,st_ctl,zu,nd)
   use m_HACApK_solve_cuda
   implicit none
   type(st_HACApK_leafmtx_acc) :: st_lf_acc
   type(st_HACApK_lcontrol) :: st_ctl
   real*8 :: zau(*),zu(*)
   integer*4,intent(in) :: nd
   real*8,dimension(:),allocatable :: zbut
   real*8,dimension(:),allocatable :: zaut
   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),ltmp(:)
   real*8,dimension(:),pointer :: a1,a2,a3
   real*8 :: zbutt, zautt
   integer*4 :: ii,il,it,itt,ill,ilil, idx1,idx2,idx3
   integer*4 :: nlf, n_low, n_dense
   integer*4, dimension(:),pointer :: ndlp,ndtp,nstrtlp,nstrttp,ktp
   integer*4, dimension(:,:),pointer :: subidt,subidl
   integer*4, dimension(:),pointer :: subidd
   integer*4, dimension(:),pointer :: ktpre
   real*8, dimension(:),pointer :: zbu_tmp
   integer*4 :: zbu_len
   logical,save :: flag = .true.
#if OPT>=7
   integer*4, dimension(:,:),pointer :: idxp1,idxp2
   integer*4, dimension(:),pointer :: idxp3
#else
   integer*4, dimension(:),pointer :: idxp1,idxp2,idxp3
#endif
!   real*8 :: zbuttt(32)
   nlf = st_lf_acc%nlf
   n_low = st_lf_acc%nlfkt
   n_dense = st_lf_acc%nlf - st_lf_acc%nlfkt

   ndlp    => st_lf_acc%ndl
   ndtp    => st_lf_acc%ndt
   nstrtlp => st_lf_acc%nstrtl
   nstrttp => st_lf_acc%nstrtt
   ktp     => st_lf_acc%kt
   a1      => st_lf_acc%a1
   a2      => st_lf_acc%a2
   a3      => st_lf_acc%a3
   idxp1    => st_lf_acc%idx1
   idxp2    => st_lf_acc%idx2
   idxp3    => st_lf_acc%idx3
   subidt   => st_lf_acc%subidt
   subidl   => st_lf_acc%subidl
   subidd   => st_lf_acc%subidd
   ktpre    => st_lf_acc%ktpre
   !$acc data &
   !$acc copy(zau(1:nd)) &
   !$acc copyin(zu(1:nd),ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,subidt,subidl,subidd,ktpre)

   allocate(zbu_tmp(st_lf_acc%lenkt))
   !$acc data create(zbu_tmp)
   !$acc kernels
   zbu_tmp(:) = 0
   !$acc end kernels
   
   !$acc host_data use_device(zau,zu,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,subidt,subidl,subidd,ktpre,zbu_tmp)
#if OPT>=7
   call HACApK_adot_body_lfmtx_cuda_wrapper2                                          &
        (zau,zu,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
        nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp)
#else
   call HACApK_adot_body_lfmtx_cuda_wrapper                                           &
        (zau,zu,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
        nd,nlf,n_low,n_dense)
#endif
   !$acc end host_data

   !$acc end data
   deallocate(zbu_tmp)
   !$acc end data
 end subroutine HACApK_adot_body_lfmtx_cuda

#else
 subroutine HACApK_adot_body_lfmtx_acc_blocked(zau,st_lf_acc,st_ctl,zu,nd)
   implicit none
   type(st_HACApK_leafmtx_acc) :: st_lf_acc
   type(st_HACApK_lcontrol) :: st_ctl
   real*8 :: zau(*),zu(*)
   integer*4, intent(in) :: nd
   real*8,dimension(:),allocatable :: zbut
   real*8,dimension(:),allocatable :: zaut
   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),ltmp(:)
   real*8,dimension(:),pointer :: a1,a2,a3
   real*8 :: zbutt, zautt
   integer*4 :: ii,il,it,itt,ill,ilil, idx1,idx2,idx3
   integer*4 :: n_low, n_dense
   integer*4, dimension(:),pointer :: ndlp,ndtp,nstrtlp,nstrttp,ktp
   integer*4, dimension(:),pointer :: idxp1,idxp2,idxp3
   integer*4 :: ndl,ndt,ns,nstrtl,nstrtt,kt
   real*8 :: zbuttt(32)
   n_low = st_lf_acc%nlfkt
   n_dense = st_lf_acc%nlf - st_lf_acc%nlfkt
   
   ndlp    => st_lf_acc%ndl
   ndtp    => st_lf_acc%ndt
   nstrtlp => st_lf_acc%nstrtl
   nstrttp => st_lf_acc%nstrtt
   ktp     => st_lf_acc%kt
   a1      => st_lf_acc%a1
   a2      => st_lf_acc%a2
   a3      => st_lf_acc%a3
   idxp1    => st_lf_acc%idx1
   idxp2    => st_lf_acc%idx2
   idxp3    => st_lf_acc%idx3
   !$acc data &
   !$acc copy(zau(1:nd)) &
   !$acc copyin(zu(1:nd),ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3)
   
   !$acc kernels
   !$acc loop independent gang private(zbuttt)
   do ii=1,n_low
      ndl   =ndlp(ii)   ; ndt   =ndtp(ii)   ; ns=ndl*ndt
      nstrtl=nstrtlp(ii); nstrtt=nstrttp(ii)
      idx1 = idxp1(ii)  ; idx2 = idxp2(ii)
      kt=ktp(ii)
      
      do ilil=1,kt,32
         do il=ilil,min(kt,ilil+31)
            zbutt = 0.0d0
            !$acc loop vector independent reduction(+:zbutt)
            do it=1,ndt
               itt=it+nstrtt-1
               zbutt=zbutt+a1(it+(il-1)*ndt+idx1-1)*zu(itt)
            enddo
            zbuttt(il-ilil+1) = zbutt
         enddo
         do il=ilil,min(kt,ilil+31)
            !$acc loop vector independent
            do it=1,ndl
               ill=it+nstrtl-1 
               !$acc atomic
               zau(ill)=zau(ill)+a2(it+(il-1)*ndl+idx2-1)*zbuttt(il-ilil+1)
            enddo
         enddo
      enddo
      
   enddo
   !$acc end kernels
   
   !$acc kernels
   !$acc loop independent gang
   do ii=n_low+1,n_low+n_dense
      ndl   =ndlp(ii)   ; ndt   =ndtp(ii)   ; ns=ndl*ndt
      nstrtl=nstrtlp(ii); nstrtt=nstrttp(ii)
      idx3 = idxp3(ii)
      !$acc loop independent
      do il=1,ndl
         ill=il+nstrtl-1
         zautt = 0.0d0
         !$acc loop independent reduction(+:zautt)
         do it=1,ndt
            itt=it+nstrtt-1
            zautt = zautt + a3(it+(il-1)*ndt+idx3-1)*zu(itt)
         enddo
         !$acc atomic
         zau(ill)=zau(ill)+zautt
      enddo
   enddo
   !$acc end kernels
   
   !$acc end data
   
 end subroutine HACApK_adot_body_lfmtx_acc_blocked
 

 subroutine HACApK_adot_body_lfmtx_acc(zau,st_lf_acc,st_ctl,zu,nd)
   implicit none
   type(st_HACApK_leafmtx_acc) :: st_lf_acc
   type(st_HACApK_lcontrol) :: st_ctl
   real*8 :: zau(*),zu(*)
   integer*4, intent(in) :: nd
   real*8,dimension(:),allocatable :: zbut
   real*8,dimension(:),allocatable :: zaut
   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),ltmp(:)
   real*8,dimension(:),pointer :: a1,a2,a3
   real*8 :: zbutt, zautt
   integer*4 :: ii, idx1,idx2,idx3
   integer*4 :: n_low, n_dense
   integer*4, dimension(:),pointer :: ndlp,ndtp,nstrtlp,nstrttp,ktp
   integer*4, dimension(:),pointer :: idxp1,idxp2,idxp3
   integer*4 :: il,it,itt,ill
   integer*4 :: ndl,ndt,ns,nstrtl,nstrtt,kt
   n_low = st_lf_acc%nlfkt
   n_dense = st_lf_acc%nlf - st_lf_acc%nlfkt
   
   ndlp    => st_lf_acc%ndl
   ndtp    => st_lf_acc%ndt
   nstrtlp => st_lf_acc%nstrtl
   nstrttp => st_lf_acc%nstrtt
   ktp     => st_lf_acc%kt
   a1      => st_lf_acc%a1
   a2      => st_lf_acc%a2
   a3      => st_lf_acc%a3
   idxp1    => st_lf_acc%idx1
   idxp2    => st_lf_acc%idx2
   idxp3    => st_lf_acc%idx3
   !$acc data &
   !$acc copy(zau(1:nd)) &
   !$acc copyin(zu(1:nd),ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3)
   
   !$acc kernels
   !$acc loop independent gang
   do ii=1,n_low
      ndl   =ndlp(ii)   ; ndt   =ndtp(ii)   ; ns=ndl*ndt
      nstrtl=nstrtlp(ii); nstrtt=nstrttp(ii)
      idx1 = idxp1(ii)  ; idx2 = idxp2(ii)
      kt=ktp(ii)
      
      do il=1,kt
         zbutt = 0.0d0
         !$acc loop vector independent reduction(+:zbutt)
         do it=1,ndt
            itt=it+nstrtt-1
            zbutt=zbutt+a1(it+(il-1)*ndt+idx1-1)*zu(itt)
         enddo
         !$acc loop vector independent
         do it=1,ndl
            ill=it+nstrtl-1 
            !$acc atomic
            zau(ill)=zau(ill)+a2(it+(il-1)*ndl+idx2-1)*zbutt
         enddo
      enddo
      
   enddo
   !$acc end kernels
   
   !$acc kernels
   !$acc loop independent gang
   do ii=n_low+1,n_low+n_dense
      ndl   =ndlp(ii)   ; ndt   =ndtp(ii)   ; ns=ndl*ndt
      nstrtl=nstrtlp(ii); nstrtt=nstrttp(ii)
      idx3 = idxp3(ii)
      !$acc loop independent
      do il=1,ndl
         ill=il+nstrtl-1
         zautt = 0.0d0
         !$acc loop independent reduction(+:zautt)
         do it=1,ndt
            itt=it+nstrtt-1
            zautt = zautt + a3(it+(il-1)*ndt+idx3-1)*zu(itt)
         enddo
         !$acc atomic
         zau(ill)=zau(ill)+zautt
      enddo
   enddo
   !$acc end kernels
   
   !$acc end data
   
 end subroutine HACApK_adot_body_lfmtx_acc
#endif

 subroutine HACApK_adot_body_lfmtx_acc_unified(zau,st_leafmtxp,st_ctl,zu,nd)
   implicit none
   type(st_HACApK_leafmtxp) :: st_leafmtxp
   type(st_HACApK_lcontrol) :: st_ctl
   real*8 :: zau(*),zu(*)
   integer*4, intent(in) :: nd
   real*8,dimension(:),allocatable :: zbut
   real*8,dimension(:),allocatable :: zaut
   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),ltmp(:)
   real*8,dimension(:,:),pointer :: a1,a2
   real*8 :: zbutt, zautt
   integer*4 :: il,it,itt,ill,ip
   integer*4 :: ndl,ndt,ns,nstrtl,nstrtt,kt
   ! 1000 format(5(a,i10)/)
   ! 2000 format(5(a,f10.4)/)
   
   ! lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;ltmp(0:) => st_ctl%lthr
   ! lpmd => st_ctl%lpmd(:);ltmp(0:) => st_ctl%lthr
   
   
   !$acc data copyin(zau(1:nd),zu(1:nd))

   !$acc kernels
   !$acc loop independent gang
   do ip=1,st_leafmtxp%nlf
      ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
      nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
      kt=st_leafmtxp%st_lf(ip)%kt
      a1 => st_leafmtxp%st_lf(ip)%a1
      a2 => st_leafmtxp%st_lf(ip)%a2

      if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
         kt=st_leafmtxp%st_lf(ip)%kt
      
         do il=1,kt
            zbutt = 0.0d0
            !$acc loop vector independent reduction(+:zbutt)
            do it=1,ndt
               itt=it+nstrtt-1
               zbutt=zbutt+a1(it,il)*zu(itt)
            enddo
         !$acc loop vector independent
            do it=1,ndl
               ill=it+nstrtl-1 
               !$acc atomic
               zau(ill)=zau(ill)+a2(it,il)*zbutt
            enddo
         enddo
      elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
         !$acc loop independent
         do il=1,ndl
            ill=il+nstrtl-1
            zautt = 0.0d0
            !$acc loop independent reduction(+:zautt)
            do it=1,ndt
               itt=it+nstrtt-1
               zautt = zautt + a1(it,il)*zu(itt)
            enddo
            !$acc atomic
            zau(ill)=zau(ill)+zautt
         enddo
      end if
   end do
   !$acc end kernels

   !$acc end data
 
 end subroutine HACApK_adot_body_lfmtx_acc_unified

 
!***HACApK_adotsub_lfmtx_p
 subroutine HACApK_adotsub_lfmtx_p(zr,st_leafmtxp,st_ctl,zu,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zu(nd),zr(nd)
 real*8,dimension(:),allocatable :: zau
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;
 allocate(zau(nd))
 call HACApK_adot_lfmtx_p(zau,st_leafmtxp,st_ctl,zu,nd)
 zr(1:nd)=zr(1:nd)-zau(1:nd)
 deallocate(zau)
 end subroutine HACApK_adotsub_lfmtx_p
 
!***HACApK_adotsub_lfmtx_hyp
 subroutine HACApK_adotsub_lfmtx_hyp(zr,zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zr(*),zau(*),zu(*),wws(*),wwr(*)
 integer*4 :: isct(*),irct(*)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 call HACApK_adot_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp workshare
 zr(1:nd)=zr(1:nd)-zau(1:nd)
!$omp end workshare
 end subroutine HACApK_adotsub_lfmtx_hyp

!***HACApK_adotsub_lfmtx_acc
 subroutine HACApK_adotsub_lfmtx_acc_unified(zr,zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
   implicit none
   type(st_HACApK_leafmtxp) :: st_leafmtxp
   type(st_HACApK_lcontrol) :: st_ctl
   real*8 :: zr(*),zau(*),zu(*),wws(*),wwr(*)
   integer*4 :: isct(*),irct(*)
   integer*4, intent(in) :: nd
   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
1000 format(5(a,i10)/)
2000 format(5(a,f10.4)/)
   
   lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
   call HACApK_adot_lfmtx_acc_unified(zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
   !$acc kernels
   zr(1:nd)=zr(1:nd)-zau(1:nd)
   !$acc end kernels
 end subroutine HACApK_adotsub_lfmtx_acc_unified

 subroutine HACApK_adotsub_lfmtx_acc(zr,zau,st_lf_acc,st_ctl,zu,wws,wwr,isct,irct,nd)
   implicit none
   type(st_HACApK_leafmtx_acc) :: st_lf_acc
   type(st_HACApK_lcontrol) :: st_ctl
   real*8 :: zr(*),zau(*),zu(*),wws(*),wwr(*)
   integer*4 :: isct(*),irct(*)
   integer*4,intent(in) :: nd
   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
1000 format(5(a,i10)/)
2000 format(5(a,f10.4)/)
   
   lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
   call HACApK_adot_lfmtx_acc(zau,st_lf_acc,st_ctl,zu,wws,wwr,isct,irct,nd)
   !$acc data copyin(zau) copy(zr)
   !$acc kernels
   zr(1:nd)=zr(1:nd)-zau(1:nd)
   !$acc end kernels
   !$acc end data
 end subroutine HACApK_adotsub_lfmtx_acc

 
!***HACApK_bicgstab_lfmtx
 subroutine HACApK_bicgstab_lfmtx(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: u(nd),b(nd)
 real*8 :: param(*)
 real*8,dimension(:),allocatable :: zr,zshdw,zp,zt,zkp,zakp,zkt,zakt
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
   call MPI_Barrier( icomm, ierr )
   st_measure_time=MPI_Wtime()
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_bicgstab_lfmtx start'
 mstep=param(83)
 eps=param(91)
 allocate(zr(nd),zshdw(nd),zp(nd),zt(nd),zkp(nd),zakp(nd),zkt(nd),zakt(nd))
 zp(1:nd)=0.0d0; zakp(1:nd)=0.0d0
 alpha = 0.0;  beta = 0.0;  zeta = 0.0;
 zz=HACApK_dotp_d(nd, b, b); bnorm=dsqrt(zz);
 zr(:nd)=b(:nd)
 call HACApK_adotsub_lfmtx_p(zr,st_leafmtxp,st_ctl,u,nd)
 zshdw(:nd)=zr(:nd)
 zrnorm=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'Original relative residual norm =',zrnorm/bnorm
 if(zrnorm/bnorm<eps) return
! mstep=1
 do in=1,mstep
   zp(:nd) =zr(:nd)+beta*(zp(:nd)-zeta*zakp(:nd))
   zkp(:nd)=zp(:nd)
   call HACApK_adot_lfmtx_p(zakp,st_leafmtxp,st_ctl,zkp,nd)
! exit
   znom=HACApK_dotp_d(nd,zshdw,zr); zden=HACApK_dotp_d(nd,zshdw,zakp);
   alpha=znom/zden; znomold=znom;
   zt(:nd)=zr(:nd)-alpha*zakp(:nd)
   zkt(:nd)=zt(:nd)
   call HACApK_adot_lfmtx_p(zakt,st_leafmtxp,st_ctl,zkt,nd)
   znom=HACApK_dotp_d(nd,zakt,zt); zden=HACApK_dotp_d(nd,zakt,zakt);
   zeta=znom/zden;
   u(:nd)=u(:nd)+alpha*zkp(:nd)+zeta*zkt(:nd)
   zr(:nd)=zt(:nd)-zeta*zakt(:nd)
   beta=alpha/zeta*HACApK_dotp_d(nd,zshdw,zr)/znomold;
   zrnorm=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm)
   call MPI_Barrier( icomm, ierr )
   en_measure_time=MPI_Wtime()
   time = en_measure_time - st_measure_time
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,in,time,log10(zrnorm/bnorm)
   if(zrnorm/bnorm<eps) exit
 enddo
end subroutine HACApK_bicgstab_lfmtx

!***HACApK_bicgstab_lfmtx_hyp
 subroutine HACApK_bicgstab_lfmtx_hyp(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: u(nd),b(nd)
 real*8 :: param(*)
 real*8,dimension(:),allocatable :: zr,zshdw,zp,zt,zkp,zakp,zkt,zakt
 real*8,dimension(:),allocatable :: wws,wwr
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 integer*4 :: isct(2),irct(2)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)
 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
   call MPI_Barrier( icomm, ierr )
   st_measure_time=MPI_Wtime()
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_bicgstab_lfmtx_hyp start'
 mstep=param(83)
 eps=param(91)
 allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
 allocate(zr(nd),zshdw(nd),zp(nd),zt(nd),zkp(nd),zakp(nd),zkt(nd),zakt(nd))
 alpha = 0.0;  beta = 0.0;  zeta = 0.0;
 zz=HACApK_dotp_d(nd, b, b); bnorm=dsqrt(zz);
!$omp parallel
!$omp workshare
 zp(1:nd)=0.0d0; zakp(1:nd)=0.0d0
 zr(:nd)=b(:nd)
!$omp end workshare
 call HACApK_adotsub_lfmtx_hyp(zr,zshdw,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp workshare
 zshdw(:nd)=zr(:nd)
!$omp end workshare
!$omp single
 zrnorm=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm)
 if(mpinr==0) print*,'Original relative residual norm =',zrnorm/bnorm
!$omp end single
 do in=1,mstep
   if(zrnorm/bnorm<eps) exit
!$omp workshare
   zp(:nd) =zr(:nd)+beta*(zp(:nd)-zeta*zakp(:nd))
   zkp(:nd)=zp(:nd)
!$omp end workshare
   call HACApK_adot_lfmtx_hyp(zakp,st_leafmtxp,st_ctl,zkp,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp single
   znom=HACApK_dotp_d(nd,zshdw,zr); zden=HACApK_dotp_d(nd,zshdw,zakp);
   alpha=znom/zden; znomold=znom;
!$omp end single
!$omp workshare
   zt(:nd)=zr(:nd)-alpha*zakp(:nd)
   zkt(:nd)=zt(:nd)
!$omp end workshare
   call HACApK_adot_lfmtx_hyp(zakt,st_leafmtxp,st_ctl,zkt,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp single
   znom=HACApK_dotp_d(nd,zakt,zt); zden=HACApK_dotp_d(nd,zakt,zakt);
   zeta=znom/zden;
!$omp end single
!$omp workshare
   u(:nd)=u(:nd)+alpha*zkp(:nd)+zeta*zkt(:nd)
   zr(:nd)=zt(:nd)-zeta*zakt(:nd)
!$omp end workshare
!$omp single
   beta=alpha/zeta*HACApK_dotp_d(nd,zshdw,zr)/znomold;
   zrnorm=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm)
   nstp=in
   call MPI_Barrier( icomm, ierr )
   en_measure_time=MPI_Wtime()
   time = en_measure_time - st_measure_time
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,in,time,log10(zrnorm/bnorm)
!$omp end single
 enddo
!$omp end parallel
end subroutine HACApK_bicgstab_lfmtx_hyp

!***HACApK_bicgstab_lfmtx_acc
subroutine HACApK_bicgstab_lfmtx_acc_unified(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
  implicit none
  include 'mpif.h'
  type(st_HACApK_leafmtxp) :: st_leafmtxp
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: u(nd),b(nd)
  real*8 :: param(*)
  integer*4, intent(in) :: nd
  integer*4, intent(out) :: nstp,lrtrn
  real*8,dimension(:),allocatable :: zr,zshdw,zp,zt,zkp,zakp,zkt,zakt
  real*8,dimension(:),allocatable :: wws,wwr
  integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
  integer*4 :: isct(2),irct(2)
  integer*4 :: mpinr,mpilog,nrank,icomm,mstep,ierr,nlf,in
  real*8 :: zz,eps,alpha,beta,zeta,znom,zden,znomold,bnorm,zrnorm
  real*8 :: st_measure_time,en_measure_time,time
  interface
     real*8 function HACApK_dotp_d_acc(nd,za,zb)
       real*8 :: za(nd),zb(nd)
       real*8 :: sum
       integer :: i
     end function HACApK_dotp_d_acc
  end interface

1000 format(5(a,i10)/)
2000 format(5(a,f10.4)/)
  lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
  mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
  call MPI_Barrier( icomm, ierr )
  st_measure_time=MPI_Wtime()
  if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_bicgstab_lfmtx_acc_unified start'
  mstep=param(83)
  eps=param(91)
  allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
  allocate(zr(nd),zshdw(nd),zp(nd),zt(nd),zkp(nd),zakp(nd),zkt(nd),zakt(nd))
  alpha = 0.0;  beta = 0.0;  zeta = 0.0;
  
  zz=HACApK_dotp_d_acc(nd, b, b); bnorm=dsqrt(zz);
  !$acc kernels
  zp(1:nd)=0.0d0; zakp(1:nd)=0.0d0
  zr(:nd)=b(:nd)
  !$acc end kernels
  call HACApK_adotsub_lfmtx_acc_unified(zr,zshdw,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
  !$acc kernels
  zshdw(:nd)=zr(:nd)
  !$acc end kernels
  zrnorm=HACApK_dotp_d_acc(nd,zr,zr); zrnorm=dsqrt(zrnorm)
  if(mpinr==0) print*,'Original relative residual norm =',zrnorm/bnorm
  do in=1,mstep
     if(zrnorm/bnorm<eps) exit
     !$acc kernels
     zp(:nd) =zr(:nd)+beta*(zp(:nd)-zeta*zakp(:nd))
     zkp(:nd)=zp(:nd)
     !$acc end kernels
     call HACApK_adot_lfmtx_acc_unified(zakp,st_leafmtxp,st_ctl,zkp,wws,wwr,isct,irct,nd)
     znom=HACApK_dotp_d_acc(nd,zshdw,zr); zden=HACApK_dotp_d_acc(nd,zshdw,zakp);
     alpha=znom/zden; znomold=znom;
     !$acc kernels
     zt(:nd)=zr(:nd)-alpha*zakp(:nd)
     zkt(:nd)=zt(:nd)
     !$acc end kernels
     call HACApK_adot_lfmtx_acc_unified(zakt,st_leafmtxp,st_ctl,zkt,wws,wwr,isct,irct,nd)
     znom=HACApK_dotp_d_acc(nd,zakt,zt); zden=HACApK_dotp_d_acc(nd,zakt,zakt);
     zeta=znom/zden;
     !$acc kernels
     u(:nd)=u(:nd)+alpha*zkp(:nd)+zeta*zkt(:nd)
     zr(:nd)=zt(:nd)-zeta*zakt(:nd)
     !$acc end kernels
     beta=alpha/zeta*HACApK_dotp_d_acc(nd,zshdw,zr)/znomold;
     zrnorm=HACApK_dotp_d_acc(nd,zr,zr); zrnorm=dsqrt(zrnorm)
     nstp=in
     call MPI_Barrier( icomm, ierr )
     en_measure_time=MPI_Wtime()
     time = en_measure_time - st_measure_time
     if(st_ctl%param(1)>0 .and. mpinr==0) print*,in,time,log10(zrnorm/bnorm)
  enddo
end subroutine HACApK_bicgstab_lfmtx_acc_unified

!***HACApK_bicgstab_lfmtx_acc
subroutine HACApK_bicgstab_lfmtx_acc(st_lf_acc,st_ctl,u,b,param,nd,nstp,lrtrn)
  use cudafor
  implicit none
  include 'mpif.h'
  ! type(st_HACApK_leafmtxp) :: st_leafmtxp
  type(st_HACApK_leafmtx_acc) :: st_lf_acc
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: u(nd),b(nd)
  real*8 :: param(*)
  integer*4 :: nd,nstp,lrtrn
  real*8,dimension(:),allocatable :: zr,zshdw,zp,zt,zkp,zakp,zkt,zakt
  real*8,dimension(:),allocatable :: wws,wwr
  integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
  integer*4 :: isct(2),irct(2)
  integer*4 :: mpinr,mpilog,nrank,icomm,mstep,ierr,nlf,in
  real*8 :: zz,eps,alpha,beta,zeta,znom,zden,znomold,bnorm,zrnorm
  real*8 :: st_measure_time,en_measure_time,time
  interface
     real*8 function HACApK_dotp_d_acc(nd,za,zb)
       real*8 :: za(nd),zb(nd)
       real*8 :: sum
       integer :: i
     end function HACApK_dotp_d_acc
  end interface


1000 format(5(a,i10)/)
2000 format(5(a,f10.4)/)
  lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
  mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
  call MPI_Barrier( icomm, ierr )
  st_measure_time=MPI_Wtime()
  if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_bicgstab_lfmtx_acc start'
  mstep=param(83)
  eps=param(91)
  allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
  allocate(zr(nd),zshdw(nd),zp(nd),zt(nd),zkp(nd),zakp(nd),zkt(nd),zakt(nd))
  alpha = 0.0;  beta = 0.0;  zeta = 0.0;

  
  !$acc data &
  !$acc create(zr,zshdw,zp,zt,zkp,zakp,zkt,zakt,wws,wwr) &
  !$acc copy(u,b)
  
  zz=HACApK_dotp_d_acc(nd, b, b); bnorm=dsqrt(zz);
  !$acc kernels
  zp(1:nd)=0.0d0; zakp(1:nd)=0.0d0
  zr(:nd)=b(:nd)
  !$acc end kernels
  ! call HACApK_adotsub_lfmtx_acc(zr,zshdw,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
  call HACApK_adotsub_lfmtx_acc(zr,zshdw,st_lf_acc,st_ctl,u,wws,wwr,isct,irct,nd)
  !$acc kernels
  zshdw(:nd)=zr(:nd)
  !$acc end kernels
  zrnorm=HACApK_dotp_d_acc(nd,zr,zr); zrnorm=dsqrt(zrnorm)
  if(mpinr==0) print*,'Original relative residual norm =',zrnorm/bnorm
  do in=1,mstep
     if(zrnorm/bnorm<eps) exit

     !$acc kernels
     zp(:nd) =zr(:nd)+beta*(zp(:nd)-zeta*zakp(:nd))
     zkp(:nd)=zp(:nd)
     !$acc end kernels
     !   call HACApK_adot_lfmtx_acc(zakp,st_leafmtxp,st_ctl,zkp,wws,wwr,isct,irct,nd)
     call HACApK_adot_lfmtx_acc(zakp,st_lf_acc,st_ctl,zkp,wws,wwr,isct,irct,nd)
     znom=HACApK_dotp_d_acc(nd,zshdw,zr); zden=HACApK_dotp_d_acc(nd,zshdw,zakp);
     alpha=znom/zden; znomold=znom;
     !$acc kernels
     zt(:nd)=zr(:nd)-alpha*zakp(:nd)
     zkt(:nd)=zt(:nd)
     !$acc end kernels
     !   call HACApK_adot_lfmtx_acc(zakt,st_leafmtxp,st_ctl,zkt,wws,wwr,isct,irct,nd)
     call HACApK_adot_lfmtx_acc(zakt,st_lf_acc,st_ctl,zkt,wws,wwr,isct,irct,nd)
     znom=HACApK_dotp_d_acc(nd,zakt,zt); zden=HACApK_dotp_d_acc(nd,zakt,zakt);
     zeta=znom/zden;
     !$acc kernels
     u(:nd)=u(:nd)+alpha*zkp(:nd)+zeta*zkt(:nd)
     zr(:nd)=zt(:nd)-zeta*zakt(:nd)
     !$acc end kernels
     beta=alpha/zeta*HACApK_dotp_d_acc(nd,zshdw,zr)/znomold;
     zrnorm=HACApK_dotp_d_acc(nd,zr,zr); zrnorm=dsqrt(zrnorm)
 
     nstp=in
     call MPI_Barrier( icomm, ierr )
     en_measure_time=MPI_Wtime()
     time = en_measure_time - st_measure_time
     if(st_ctl%param(1)>0 .and. mpinr==0) print*,"bcg",in,time,log10(zrnorm/bnorm)
  enddo
  !$acc end data
end subroutine HACApK_bicgstab_lfmtx_acc

#if OPT>=7
subroutine HACApK_bicgstab_lfmtx_cuda(st_lf_acc,st_ctl,u,b,param,nd,nstp,lrtrn)
  use cudafor
  use m_HACApK_solve_cuda
  implicit none
  include 'mpif.h'
  type(st_HACApK_leafmtx_acc) :: st_lf_acc
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: u(nd),b(nd)
  real*8 :: param(*)
  integer*4 :: nd,nstp,lrtrn
  real*8,dimension(:),allocatable :: zr,zshdw,zp,zt,zkp,zakp,zkt,zakt
  real*8,dimension(:),allocatable :: wws,wwr
  real*8,dimension(:),allocatable :: zbu_tmp,znom_tmp,zden_tmp
  integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
  integer*4 :: isct(2),irct(2)
  integer*4,allocatable,dimension(:,:) :: isrct,isrct_d
  real*8,dimension(:,:),allocatable :: ww,ww_d
  integer*4 :: mpinr,mpilog,nrank,icomm,mstep,ierr,nlf,in
!  real*8 :: zz,eps,alpha,beta,zeta,znom,zden,znomold,bnorm,zrnorm
  real*8 :: zz,eps,bnorm,zrnorm
  real*8,dimension(1) :: alpha,beta,zeta,znom,zden,znomold
  real*8 :: st_measure_time,en_measure_time,time
  integer(kind=cuda_stream_kind) :: stream1, stream2
  integer :: istat
  integer*4 :: n_low, n_dense
  real*8,dimension(:),pointer :: a1,a2,a3
  integer*4, dimension(:),pointer :: ndlp,ndtp,nstrtlp,nstrttp,ktp
  integer*4, dimension(:,:),pointer :: idxp1,idxp2
  integer*4, dimension(:),pointer :: idxp3
  integer*4, dimension(:,:),pointer :: subidt,subidl
  integer*4, dimension(:),pointer :: subidd
  integer*4, dimension(:),pointer :: ktpre
  integer*4 :: ic,ncdp,ncsp
  integer*4 :: ISTATUS(MPI_STATUS_SIZE)
  interface
     real*8 function HACApK_dotp_d_acc(nd,za,zb)
       real*8 :: za(nd),zb(nd)
       real*8 :: sum
       integer :: i
     end function HACApK_dotp_d_acc
  end interface

1000 format(5(a,i10)/)
2000 format(5(a,f10.4)/)
  lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
  mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
  call MPI_Barrier( icomm, ierr )
  st_measure_time=MPI_Wtime()
  if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_bicgstab_lfmtx_acc start'
  mstep=param(83)
  eps=param(91)
  allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
  allocate(zr(nd),zshdw(nd),zp(nd),zt(nd),zkp(nd),zakp(nd),zkt(nd),zakt(nd))
  alpha = 0.0;  beta = 0.0;  zeta = 0.0;

  allocate(isrct(2,0:nrank-1))
  allocate(isrct_d(2,0:nrank-1))
  isrct(1,mpinr)=lnp(mpinr)
  isrct(2,mpinr)=lsp(mpinr)
  allocate(ww(maxval(lnp(0:nrank-1)),0:nrank-1))
  allocate(ww_d(maxval(lnp(0:nrank-1)),0:nrank-1))
  ww = 0.0d0; ww_d = 0.0d0
  ! do ic=1,nrank-1
  !    ncdp=mod(mpinr+ic,nrank)
  !    ncsp=mod(mpinr+nrank-ic,nrank)
  !    call MPI_SENDRECV(isrct(1,mpinr),2,MPI_INTEGER,ncdp,1, &
  !         isrct(1,ncsp),2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
  ! enddo
  isrct(1,:)=lnp(:)
  isrct(2,:)=lsp(:)

  ! allocate(isrct(2,nrank))
  ! allocate(isrct_d(2,nrank))
  ! isrct(1,mpinr+1)=lnp(mpinr+1)
  ! isrct(2,mpinr+1)=lsp(mpinr+1)
  ! allocate(ww(maxval(lnp(0:nrank-1)),nrank))
  ! allocate(ww_d(maxval(lnp(0:nrank-1)),nrank))
  ! ww = 0.0d0; ww_d = 0.0d0
  ! do ic=1,nrank-1
  !    ncdp=mod(mpinr+ic,nrank)
  !    ncsp=mod(mpinr+nrank+ic-2,nrank)
  !    call MPI_SENDRECV(isrct(1,mpinr+1),2,MPI_INTEGER,ncdp,1, &
  !         isrct(1,ncsp+1),2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
  ! enddo
  isrct_d(:,:) = isrct(:,:)

  istat = cudaStreamCreate(stream1)

  nlf = st_lf_acc%nlf
  n_low = st_lf_acc%nlfkt
  n_dense = st_lf_acc%nlf - st_lf_acc%nlfkt
  
  ndlp    => st_lf_acc%ndl
  ndtp    => st_lf_acc%ndt
  nstrtlp => st_lf_acc%nstrtl
  nstrttp => st_lf_acc%nstrtt
  ktp     => st_lf_acc%kt
  a1      => st_lf_acc%a1
  a2      => st_lf_acc%a2
  a3      => st_lf_acc%a3
  idxp1   => st_lf_acc%idx1
  idxp2   => st_lf_acc%idx2
  idxp3   => st_lf_acc%idx3
  subidt  => st_lf_acc%subidt
  subidl  => st_lf_acc%subidl
  subidd  => st_lf_acc%subidd
  ktpre   => st_lf_acc%ktpre

  allocate(zbu_tmp(st_lf_acc%lenkt))
  allocate(znom_tmp(nd))
  allocate(zden_tmp(nd))
  
  !$acc data &
  !$acc create(zr,zshdw,zp,zt,zkp,zakp,zkt,zakt,zbu_tmp,znom_tmp,zden_tmp) &
  !$acc copyin(ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3) &
  !$acc copyin(idxp1,idxp2,idxp3,subidt,subidl,subidd,ktpre) &
  !$acc copyin(alpha,beta,zeta,znomold,isrct_d,ww_d) &
  !$acc copy(u,b)
  
  zz=HACApK_dotp_d_acc(nd, b, b); bnorm=dsqrt(zz);
  !$acc kernels
  zp(1:nd)=0.0d0; zakp(1:nd)=0.0d0
  zr(:nd)=b(:nd)
  !$acc end kernels
  call HACApK_adotsub_lfmtx_acc(zr,zshdw,st_lf_acc,st_ctl,u,wws,wwr,isct,irct,nd)
  !$acc kernels
  zshdw(:nd)=zr(:nd)
  !$acc end kernels
  zrnorm=HACApK_dotp_d_acc(nd,zr,zr); zrnorm=dsqrt(zrnorm)
  if(mpinr==0) print*,'Original relative residual norm =',zrnorm/bnorm

  !$acc wait
  !$acc host_data use_device(ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,subidt,subidl,subidd,ktpre,zbu_tmp,zr,zshdw,zp,zt,zkp,zakp,zkt,zakt,alpha,beta,zeta,znomold,u,znom_tmp,zden_tmp,isrct_d,ww_d)

  do in=1,mstep

     if(zrnorm/bnorm<eps) exit 

     ! zp(:nd) =zr(:nd)+beta*(zp(:nd)-zeta*zakp(:nd))
     ! zkp(:nd)=zp(:nd)
     call zp_zkp_cuda(zp,zkp,zr,zakp,beta,zeta,nd,stream1)

     call HACApK_adot_body_lfmtx_cuda_wrapper3                                          &
          (zakp,zkp,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
          nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp,st_ctl,wws,wwr,stream1)
     ! call HACApK_adot_body_lfmtx_cuda_wrapper4                                          &
     !      (zakp,zkp,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
     !      nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp,st_ctl,wws,wwr,stream1, &
     !      isrct,isrct_d,ww,ww_d)

     ! znom=HACApK_dotp_d_acc(nd,zshdw,zr); zden=HACApK_dotp_d_acc(nd,zshdw,zakp);
     ! alpha=znom/zden; znomold=znom;
     ! zt(:nd)=zr(:nd)-alpha*zakp(:nd)
     ! zkt(:nd)=zt(:nd)
     call zt_zkt_cuda(zt,zkt,zr,zakp,zshdw,nd,znomold,alpha,znom_tmp,zden_tmp,stream1)

     call HACApK_adot_body_lfmtx_cuda_wrapper3                                          &
          (zakt,zkt,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
          nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp,st_ctl,wws,wwr,stream1)
     ! call HACApK_adot_body_lfmtx_cuda_wrapper4                                          &
     !      (zakt,zkt,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
     !      nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp,st_ctl,wws,wwr,stream1, &
     !      isrct,isrct_d,ww,ww_d)

     ! znom=HACApK_dotp_d_acc(nd,zakt,zt); zden=HACApK_dotp_d_acc(nd,zakt,zakt);
     ! zeta=znom/zden;
     ! u(:nd)=u(:nd)+alpha*zkp(:nd)+zeta*zkt(:nd)
     ! zr(:nd)=zt(:nd)-zeta*zakt(:nd)
     ! beta=alpha/zeta*HACApK_dotp_d_acc(nd,zshdw,zr)/znomold;
     ! zrnorm=HACApK_dotp_d_acc(nd,zr,zr); zrnorm=dsqrt(zrnorm)
     call u_zr_cuda(u,zr,zkp,zkt,zt,zakt,zshdw,alpha,beta,zeta,znomold,zrnorm,nd,znom_tmp,zden_tmp,stream1)

     nstp=in
     call MPI_Barrier( icomm, ierr )
     en_measure_time=MPI_Wtime()
     time = en_measure_time - st_measure_time
     if(st_ctl%param(1)>0 .and. mpinr==0) print*,"bcg",in,time,log10(zrnorm/bnorm)
  enddo
  !$acc end host_data

  !$acc end data
  deallocate(zbu_tmp)
  deallocate(znom_tmp)
  deallocate(zden_tmp)
  deallocate(isrct)
  deallocate(ww,ww_d)
  istat = cudaStreamDestroy(stream1)

end subroutine HACApK_bicgstab_lfmtx_cuda


subroutine HACApK_bicgstab_lfmtx_cuda2(st_lf_acc,st_ctl,u,b,param,nd,nstp,lrtrn)
  use cudafor
  use m_HACApK_solve_cuda
  implicit none
  include 'mpif.h'
  type(st_HACApK_leafmtx_acc) :: st_lf_acc
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: u(nd),b(nd)
  real*8 :: param(*)
  integer*4 :: nd,nstp,lrtrn
  real*8,dimension(:),allocatable :: zr,zshdw,zp,zt,zkp,zakp,zkt,zakt
  real*8,dimension(:),allocatable :: wws,wwr
  real*8,dimension(:),allocatable :: zbu_tmp,znom_tmp,zden_tmp
  integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
  integer*4 :: isct(2),irct(2)
  integer*4,allocatable,dimension(:,:) :: isrct,isrct_d
  real*8,dimension(:,:),allocatable :: ww,ww_d
  integer*4 :: mpinr,mpilog,nrank,icomm,mstep,ierr,nlf,in
!  real*8 :: zz,eps,alpha,beta,zeta,znom,zden,znomold,bnorm,zrnorm
  real*8 :: zz,eps,bnorm,zrnorm
  real*8,dimension(1) :: alpha,beta,zeta,znom,zden,znomold
  real*8 :: st_measure_time,en_measure_time,time
  integer(kind=cuda_stream_kind) :: stream1, stream2
  integer :: istat
  integer*4 :: n_low, n_dense
  real*8,dimension(:),pointer :: a1,a2,a3
  integer*4, dimension(:),pointer :: ndlp,ndtp,nstrtlp,nstrttp,ktp
  integer*4, dimension(:,:),pointer :: idxp1,idxp2
  integer*4, dimension(:),pointer :: idxp3
  integer*4, dimension(:,:),pointer :: subidt,subidl
  integer*4, dimension(:),pointer :: subidd
  integer*4, dimension(:),pointer :: ktpre
  integer*4 :: ic,ncdp,ncsp
  integer*4 :: ISTATUS(MPI_STATUS_SIZE)
  interface
     real*8 function HACApK_dotp_d_acc(nd,za,zb)
       real*8 :: za(nd),zb(nd)
       real*8 :: sum
       integer :: i
     end function HACApK_dotp_d_acc
  end interface

1000 format(5(a,i10)/)
2000 format(5(a,f10.4)/)
  lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
  mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
  call MPI_Barrier( icomm, ierr )
  st_measure_time=MPI_Wtime()
  if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_bicgstab_lfmtx_acc start'
  mstep=param(83)
  eps=param(91)
  allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
  allocate(zr(nd),zshdw(nd),zp(nd),zt(nd),zkp(nd),zakp(nd),zkt(nd),zakt(nd))
  alpha = 0.0;  beta = 0.0;  zeta = 0.0;

  allocate(isrct(2,0:nrank-1))
  allocate(isrct_d(2,0:nrank-1))
  isrct(1,mpinr)=lnp(mpinr)
  isrct(2,mpinr)=lsp(mpinr)
  allocate(ww(maxval(lnp(0:nrank-1)),0:nrank-1))
  allocate(ww_d(maxval(lnp(0:nrank-1)),0:nrank-1))
  ww = 0.0d0; ww_d = 0.0d0
  ! do ic=1,nrank-1
  !    ncdp=mod(mpinr+ic,nrank)
  !    ncsp=mod(mpinr+nrank-ic,nrank)
  !    call MPI_SENDRECV(isrct(1,mpinr),2,MPI_INTEGER,ncdp,1, &
  !         isrct(1,ncsp),2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
  ! enddo
  isrct(1,:)=lnp(:)
  isrct(2,:)=lsp(:)

  ! allocate(isrct(2,nrank))
  ! allocate(isrct_d(2,nrank))
  ! isrct(1,mpinr+1)=lnp(mpinr+1)
  ! isrct(2,mpinr+1)=lsp(mpinr+1)
  ! allocate(ww(maxval(lnp(0:nrank-1)),nrank))
  ! allocate(ww_d(maxval(lnp(0:nrank-1)),nrank))
  ! ww = 0.0d0; ww_d = 0.0d0
  ! do ic=1,nrank-1
  !    ncdp=mod(mpinr+ic,nrank)
  !    ncsp=mod(mpinr+nrank+ic-2,nrank)
  !    call MPI_SENDRECV(isrct(1,mpinr+1),2,MPI_INTEGER,ncdp,1, &
  !         isrct(1,ncsp+1),2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
  ! enddo
  isrct_d(:,:) = isrct(:,:)

  istat = cudaStreamCreate(stream1)

  nlf = st_lf_acc%nlf
  n_low = st_lf_acc%nlfkt
  n_dense = st_lf_acc%nlf - st_lf_acc%nlfkt
  
  ndlp    => st_lf_acc%ndl
  ndtp    => st_lf_acc%ndt
  nstrtlp => st_lf_acc%nstrtl
  nstrttp => st_lf_acc%nstrtt
  ktp     => st_lf_acc%kt
  a1      => st_lf_acc%a1
  a2      => st_lf_acc%a2
  a3      => st_lf_acc%a3
  idxp1   => st_lf_acc%idx1
  idxp2   => st_lf_acc%idx2
  idxp3   => st_lf_acc%idx3
  subidt  => st_lf_acc%subidt
  subidl  => st_lf_acc%subidl
  subidd  => st_lf_acc%subidd
  ktpre   => st_lf_acc%ktpre

  allocate(zbu_tmp(st_lf_acc%lenkt))
  allocate(znom_tmp(nd))
  allocate(zden_tmp(nd))
  
  !$acc data &
  !$acc create(zr,zshdw,zp,zt,zkp,zakp,zkt,zakt,zbu_tmp,znom_tmp,zden_tmp) &
  !$acc copyin(ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3) &
  !$acc copyin(idxp1,idxp2,idxp3,subidt,subidl,subidd,ktpre) &
  !$acc copyin(alpha,beta,zeta,znomold,isrct_d,ww_d) &
  !$acc copy(u,b)
  
  zz=HACApK_dotp_d_acc(nd, b, b); bnorm=dsqrt(zz);
  !$acc kernels
  zp(1:nd)=0.0d0; zakp(1:nd)=0.0d0
  zr(:nd)=b(:nd)
  !$acc end kernels
  call HACApK_adotsub_lfmtx_acc(zr,zshdw,st_lf_acc,st_ctl,u,wws,wwr,isct,irct,nd)
  !$acc kernels
  zshdw(:nd)=zr(:nd)
  !$acc end kernels
  zrnorm=HACApK_dotp_d_acc(nd,zr,zr); zrnorm=dsqrt(zrnorm)
  if(mpinr==0) print*,'Original relative residual norm =',zrnorm/bnorm

  !$acc wait
  !$acc host_data use_device(ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,subidt,subidl,subidd,ktpre,zbu_tmp,zr,zshdw,zp,zt,zkp,zakp,zkt,zakt,alpha,beta,zeta,znomold,u,znom_tmp,zden_tmp,isrct_d,ww_d)

  do in=1,mstep

#if 0
     if(zrnorm/bnorm<eps) exit 
#endif

     ! zp(:nd) =zr(:nd)+beta*(zp(:nd)-zeta*zakp(:nd))
     ! zkp(:nd)=zp(:nd)
     call zp_zkp_cuda(zp,zkp,zr,zakp,beta,zeta,nd,stream1)

     call HACApK_adot_body_lfmtx_cuda_wrapper3                                          &
          (zakp,zkp,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
          nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp,st_ctl,wws,wwr,stream1)
     ! call HACApK_adot_body_lfmtx_cuda_wrapper4                                          &
     !      (zakp,zkp,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
     !      nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp,st_ctl,wws,wwr,stream1, &
     !      isrct,isrct_d,ww,ww_d)

     ! znom=HACApK_dotp_d_acc(nd,zshdw,zr); zden=HACApK_dotp_d_acc(nd,zshdw,zakp);
     ! alpha=znom/zden; znomold=znom;
     ! zt(:nd)=zr(:nd)-alpha*zakp(:nd)
     ! zkt(:nd)=zt(:nd)
     call zt_zkt_cuda(zt,zkt,zr,zakp,zshdw,nd,znomold,alpha,znom_tmp,zden_tmp,stream1)
!     call zt_zkt_cuda2(zt,zkt,zr,zakp,zshdw,nd,znomold,alpha,znom_tmp,zden_tmp,stream1)

     call HACApK_adot_body_lfmtx_cuda_wrapper3                                          &
          (zakt,zkt,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
          nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp,st_ctl,wws,wwr,stream1)
     ! call HACApK_adot_body_lfmtx_cuda_wrapper4                                          &
     !      (zakt,zkt,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
     !      nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp,st_ctl,wws,wwr,stream1, &
     !      isrct,isrct_d,ww,ww_d)

     ! znom=HACApK_dotp_d_acc(nd,zakt,zt); zden=HACApK_dotp_d_acc(nd,zakt,zakt);
     ! zeta=znom/zden;
     ! u(:nd)=u(:nd)+alpha*zkp(:nd)+zeta*zkt(:nd)
     ! zr(:nd)=zt(:nd)-zeta*zakt(:nd)
     ! beta=alpha/zeta*HACApK_dotp_d_acc(nd,zshdw,zr)/znomold;
     ! zrnorm=HACApK_dotp_d_acc(nd,zr,zr); zrnorm=dsqrt(zrnorm)
     call u_zr_cuda(u,zr,zkp,zkt,zt,zakt,zshdw,alpha,beta,zeta,znomold,zrnorm,nd,znom_tmp,zden_tmp,stream1)
!     call u_zr_cuda2(u,zr,zkp,zkt,zt,zakt,zshdw,alpha,beta,zeta,znomold,zrnorm,nd,znom_tmp,zden_tmp,stream1)

     nstp=in
     call MPI_Barrier( icomm, ierr )
     en_measure_time=MPI_Wtime()
     time = en_measure_time - st_measure_time
     if(st_ctl%param(1)>0 .and. mpinr==0) print*,"bcg",in,time,log10(zrnorm/bnorm)
  enddo
  !$acc end host_data

  !$acc end data
  deallocate(zbu_tmp)
  deallocate(znom_tmp)
  deallocate(zden_tmp)
  deallocate(isrct)
  deallocate(ww,ww_d)
  istat = cudaStreamDestroy(stream1)

end subroutine HACApK_bicgstab_lfmtx_cuda2
#endif

! real*8 function HACApK_bicgstab_loop_body_cuda(st_lf_acc,st_ctl,zp,zr,zakp,zkp,wws,wwr,isct,irct,nd,zshdw,zt,zkt,zakt,u)
!   implicit none
!   type(st_HACApK_leafmtx_acc) :: st_lf_acc
!   type(st_HACApK_lcontrol) :: st_ctl
!   real*8,dimension(:) :: zp,zr,zakp,zkp,zshdw,zt,zkt,zakt,u
!   real*8,dimension(:) :: wws,wwr
!   integer*4 :: nd
!   integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
!   integer*4 :: isct(2),irct(2)
!   integer*4 :: mpinr,mpilog,nrank,icomm,mstep,ierr,nlf,in
!   real*8 :: zz,eps,alpha,beta,zeta,znom,zden,znomold,bnorm,zrnorm

!    use m_HACApK_solve_cuda
!    implicit none
!    type(st_HACApK_leafmtx_acc) :: st_lf_acc
!    type(st_HACApK_lcontrol) :: st_ctl
!    real*8 :: zau(*),zu(*)
!    integer*4,intent(in) :: nd
!    real*8,dimension(:),allocatable :: zbut
!    real*8,dimension(:),allocatable :: zaut
!    integer*4,pointer :: lpmd(:),lnp(:),lsp(:),ltmp(:)
!    real*8,dimension(:),pointer :: a1,a2,a3
!    real*8 :: zbutt, zautt
!    integer*4 :: ii,il,it,itt,ill,ilil, idx1,idx2,idx3
!    integer*4 :: nlf, n_low, n_dense
!    integer*4, dimension(:),pointer :: ndlp,ndtp,nstrtlp,nstrttp,ktp
! #if OPT>=7
!    integer*4, dimension(:,:),pointer :: idxp1,idxp2
!    integer*4, dimension(:),pointer :: idxp3
!    integer*4, dimension(:,:),pointer :: subidt,subidl
!    integer*4, dimension(:),pointer :: subidd
!    integer*4, dimension(:),pointer :: ktpre
!    real*8, dimension(:),pointer :: zbu_tmp
!    integer*4 :: zbu_len
!    logical,save :: flag = .true.
! #else
!    integer*4, dimension(:),pointer :: idxp1,idxp2,idxp3
! #endif
! !   real*8 :: zbuttt(32)
!    nlf = st_lf_acc%nlf
!    n_low = st_lf_acc%nlfkt
!    n_dense = st_lf_acc%nlf - st_lf_acc%nlfkt

!    ndlp    => st_lf_acc%ndl
!    ndtp    => st_lf_acc%ndt
!    nstrtlp => st_lf_acc%nstrtl
!    nstrttp => st_lf_acc%nstrtt
!    ktp     => st_lf_acc%kt
!    a1      => st_lf_acc%a1
!    a2      => st_lf_acc%a2
!    a3      => st_lf_acc%a3
!    idxp1    => st_lf_acc%idx1
!    idxp2    => st_lf_acc%idx2
!    idxp3    => st_lf_acc%idx3
!    subidt   => st_lf_acc%subidt
!    subidl   => st_lf_acc%subidl
!    subidd   => st_lf_acc%subidd
!    ktpre    => st_lf_acc%ktpre
!    !$acc data &
!    !$acc copy(zau(1:nd)) &
!    !$acc copyin(zu(1:nd),ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,subidt,subidl,subidd,ktpre)

!    allocate(zbu_tmp(st_lf_acc%lenkt))
!    !$acc data create(zbu_tmp)
!    !$acc kernels
!    zbu_tmp(:) = 0
!    !$acc end kernels
   
!    !$acc host_data use_device(zau,zu,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,subidt,subidl,subidd,ktpre,zbu_tmp)
! #if OPT>=7
!    call HACApK_adot_body_lfmtx_cuda_wrapper2                                          &
!         (zau,zu,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
!         nd,nlf,n_low,n_dense,subidt,subidl,subidd,ktpre,zbu_tmp)
! #else
!    call HACApK_adot_body_lfmtx_cuda_wrapper                                           &
!         (zau,zu,ndlp,ndtp,nstrtlp,nstrttp,ktp,a1,a2,a3,idxp1,idxp2,idxp3,             &
!         nd,nlf,n_low,n_dense)
! #endif
!    !$acc end host_data

!    !$acc end data
!    deallocate(zbu_tmp)
!    !$acc end data


!   !$acc kernels
!   zp(:nd) =zr(:nd)+beta*(zp(:nd)-zeta*zakp(:nd))
!   zkp(:nd)=zp(:nd)
!   !$acc end kernels
!   !   call HACApK_adot_lfmtx_acc(zakp,st_leafmtxp,st_ctl,zkp,wws,wwr,isct,irct,nd)
!   call HACApK_adot_lfmtx_acc(zakp,st_lf_acc,st_ctl,zkp,wws,wwr,isct,irct,nd)
!   znom=HACApK_dotp_d_acc(nd,zshdw,zr); zden=HACApK_dotp_d_acc(nd,zshdw,zakp);
!   alpha=znom/zden; znomold=znom;
!   !$acc kernels
!   zt(:nd)=zr(:nd)-alpha*zakp(:nd)
!   zkt(:nd)=zt(:nd)
!   !$acc end kernels
!   !   call HACApK_adot_lfmtx_acc(zakt,st_leafmtxp,st_ctl,zkt,wws,wwr,isct,irct,nd)
!   call HACApK_adot_lfmtx_acc(zakt,st_lf_acc,st_ctl,zkt,wws,wwr,isct,irct,nd)
!   znom=HACApK_dotp_d_acc(nd,zakt,zt); zden=HACApK_dotp_d_acc(nd,zakt,zakt);
!   zeta=znom/zden;
!   !$acc kernels
!   u(:nd)=u(:nd)+alpha*zkp(:nd)+zeta*zkt(:nd)
!   zr(:nd)=zt(:nd)-zeta*zakt(:nd)
!   !$acc end kernels
!   beta=alpha/zeta*HACApK_dotp_d_acc(nd,zshdw,zr)/znomold;
!   zrnorm=HACApK_dotp_d_acc(nd,zr,zr); zrnorm=dsqrt(zrnorm)

! end subroutine HACApK_bicgstab_lfmtx_acc


!***HACApK_gcrm_lfmtx
 subroutine HACApK_gcrm_lfmtx(st_leafmtxp,st_ctl,st_bemv,u,b,param,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: u(nd),b(nd)
 real*8 :: param(*)
 real*8,dimension(:),allocatable :: zr,zar,capap
 real*8,dimension(:,:),allocatable,target :: zp,zap
 real*8,pointer :: zq(:)
 real*8,dimension(:),allocatable :: wws,wwr
 integer*4 :: isct(2),irct(2)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
   call MPI_Barrier( icomm, ierr )
   st_measure_time=MPI_Wtime()
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'gcr_lfmtx_hyp start'
 mstep=param(83)
 mreset=param(87)
 eps=param(91)
 allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
 allocate(zr(nd),zar(nd),zp(nd,mreset),zap(nd,mreset),capap(mreset))
 alpha = 0.0
 zz=HACApK_dotp_d(nd, b, b); bnorm=dsqrt(zz);
 call HACApK_adot_lfmtx_hyp(zar,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
 zr(:nd)=b(:nd)-zar(:nd)
 zp(:nd,1)=zr(:nd)
 zrnorm2=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm2)
   call MPI_Barrier( icomm, ierr )
   en_measure_time=MPI_Wtime()
   time = en_measure_time - st_measure_time
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,0,time,log10(zrnorm/bnorm)
 if(zrnorm/bnorm<eps) return
 call HACApK_adot_lfmtx_hyp(zap(:nd,1),st_leafmtxp,st_ctl,zp(:nd,1),wws,wwr,isct,irct,nd)
 do in=1,mstep
   ik=mod(in-1,mreset)+1
   zq=>zap(:nd,ik)
   znom=HACApK_dotp_d(nd,zq,zr); capap(ik)=HACApK_dotp_d(nd,zq,zq)
   alpha=znom/capap(ik)
   u(:nd)=u(:nd)+alpha*zp(:nd,ik)
   zr(:nd)=zr(:nd)-alpha*zq(:nd)
   zrnomold=zrnorm2
   zrnorm2=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm2)
   call MPI_Barrier( icomm, ierr )
   en_measure_time=MPI_Wtime()
   time = en_measure_time - st_measure_time
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,in,time,log10(zrnorm/bnorm)
   if(zrnorm/bnorm<eps .or. in==mstep) exit
   call HACApK_adot_lfmtx_hyp(zar,st_leafmtxp,st_ctl,zr,wws,wwr,isct,irct,nd)
   ikn=mod(in,mreset)+1
   zp(:nd,ikn)=zr(:nd)
   zap(:nd,ikn)=zar(:nd)
   do il=1,ik
     zq=>zap(:nd,il)
     znom=HACApK_dotp_d(nd,zq,zar)
     beta=-znom/capap(il)
     zp(:nd,ikn) =zp(:nd,ikn)+beta*zp(:nd,il)
     zap(:nd,ikn)=zap(:nd,ikn)+beta*zq(:nd)
   enddo
 enddo
 nstp=in
end subroutine

!***HACApK_measurez_time_ax_lfmtx
 subroutine HACApK_measurez_time_ax_lfmtx(st_leafmtxp,st_ctl,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8,dimension(:),allocatable :: wws,wwr,u,b
 integer*4 :: isct(2),irct(2)
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr; param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 mstep=param(99)
 allocate(u(nd),b(nd),wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
 !$omp parallel private(il)
 do il=1,mstep
    u(:)=1.0; b(:)=1.0
    call HACApK_adot_lfmtx_hyp(u,st_leafmtxp,st_ctl,b,wws,wwr,isct,irct,nd)
 enddo
 !$omp end parallel
 deallocate(wws,wwr)
end subroutine HACApK_measurez_time_ax_lfmtx

 subroutine HACApK_measurez_time_ax_lfmtx_acc_unified(st_leafmtxp,st_ctl,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8,dimension(:),allocatable :: wws,wwr,u,b
 integer*4 :: isct(2),irct(2)
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr; param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 mstep=param(99)
 allocate(u(nd),b(nd),wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
 do il=1,mstep
    !$acc kernels
    u(:)=1.0
    !$acc end kernels
    !$acc kernels
    b(:)=1.0
    !$acc end kernels
    call HACApK_adot_lfmtx_acc_unified(u,st_leafmtxp,st_ctl,b,wws,wwr,isct,irct,nd)
 enddo
 deallocate(wws,wwr)
end subroutine HACApK_measurez_time_ax_lfmtx_acc_unified

 subroutine HACApK_measurez_time_ax_lfmtx_acc(st_lf_acc,st_ctl,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtx_acc) :: st_lf_acc
 type(st_HACApK_lcontrol) :: st_ctl
 real*8,dimension(:),allocatable :: wws,wwr,u,b
 integer*4 :: isct(2),irct(2)
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr; param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 mstep=param(99)
 allocate(u(nd),b(nd),wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
 !$acc data create(u,b)
 do il=1,mstep
    !$acc kernels
    u(:)=1.0
    !$acc end kernels
    !$acc kernels
    b(:)=1.0
    !$acc end kernels
    call HACApK_adot_lfmtx_acc(u,st_lf_acc,st_ctl,b,wws,wwr,isct,irct,nd)
 enddo
 !$acc end data
 deallocate(wws,wwr)
end subroutine HACApK_measurez_time_ax_lfmtx_acc

!***HACApK_adot_pmt_lfmtx_p
 integer function HACApK_adot_pmt_lfmtx_p(st_leafmtxp,st_bemv,st_ctl,aww,ww)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: ww(st_bemv%nd),aww(st_bemv%nd)
 real*8,dimension(:),allocatable :: u,au
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lrtrn=0
 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;lod => st_ctl%lod(:)
 mpinr=st_ctl%lpmd(3); icomm=st_ctl%lpmd(1); nd=st_bemv%nd
 allocate(u(nd),au(nd)); u(:nd)=ww(st_ctl%lod(:nd))
 call MPI_Barrier( icomm, ierr )
 call HACApK_adot_lfmtx_p(au,st_leafmtxp,st_ctl,u,nd)
 aww(st_ctl%lod(:nd))=au(:nd)
 HACApK_adot_pmt_lfmtx_p=lrtrn
end function HACApK_adot_pmt_lfmtx_p

!***HACApK_adot_pmt_lfmtx_hyp
 integer function HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp,st_bemv,st_ctl,aww,ww)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: ww(*),aww(*)
 real*8,dimension(:),allocatable :: u,au,wws,wwr
 integer*4,dimension(:),allocatable :: isct,irct
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lrtrn=0
 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;lod => st_ctl%lod(:)
 mpinr=st_ctl%lpmd(3); icomm=st_ctl%lpmd(1); nd=st_bemv%nd; nrank=st_ctl%lpmd(2)
 allocate(u(nd),au(nd),isct(2),irct(2)); u(:nd)=ww(st_ctl%lod(:nd))
 allocate(wws(maxval(st_ctl%lnp(:nrank))),wwr(maxval(st_ctl%lnp(:nrank))))
 call MPI_Barrier( icomm, ierr )
!$omp parallel
!$omp barrier
 call HACApK_adot_lfmtx_hyp(au,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp end parallel
 call MPI_Barrier( icomm, ierr )
 aww(st_ctl%lod(:nd))=au(:nd)
 HACApK_adot_pmt_lfmtx_hyp=lrtrn
end function HACApK_adot_pmt_lfmtx_hyp

 
endmodule m_HACApK_solve

