cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
c
c
c    $Date$
c    $Revision$
c
c
c     ROTATION VIA PROJECTION  (FORTRAN 77 and 90 VERSIONS).
c
c     Multipole expansions for complex-valued functions.
c
c     Requires FFT and Associated Legendre Function Libraries.
c
c     User-callable f77 routine is rotviaproj2. 
c     User-callable f90 routine is rotviaproj2f90.
c     The other routines are used internally.
c
c***********************************************************************
      subroutine rotviaproj2(beta,alpha,nterms,m1,m2,mpole,lmp,
     $     marray2,lmpn,w,lw,lused)
c***********************************************************************
c       Purpose:
c
c	Fast and stable algorithm for applying rotation operator about
c	the z-axis determined by angle alpha plus the y-axis determined
c	by angle beta. After rotation, the expansion pole is moved to
c	location (beta, alpha) in spherical coordinates (theta, phi).
c
c       The method is based on computing the induced potential and
c       its theta-derivative on the rotated equator
c       for each order (first index). The coefficients of  the rotated
c       expansion can then be obtained by FFT and projection.
c
c       There is some loss in speed over using recurrence relations 
c       but it is stable to all orders whereas the recurrence schemes 
c       are not.
c
c       If the rotation operator is to be used multiple times, and
c       memory is available, one can precompute and store the 
c       multipliers used in evalall (see below). This has not yet been
c       implemented.
c
c       Our definition of complex spherical harmonics is
c
c       Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                       Pnm(cos theta) e^(im phi), 
c       Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                       Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c       Note that we do not include the Condon-Shortley phase (-1)^m
c
C---------------------------------------------------------------------
c       INPUT:
c
c       beta:  the rotation angle about the y-axis.
c       alpha:  the rotation angle about the z-axis.
c       nterms: order of multipole expansion
c
c       m1    : max m index for first expansion   
c               NOT IMPLEMENTED but integer argument must be supplied.
c       m2    : max m index for second expansion  
c               NOT IMPLEMENTED but integer argument must be supplied.
c
c               That is, parameters m1,m2 are currently ignored.
c
C       mpole   coefficients of original multiple expansion
C       lmp     leading dim for mpole (must exceed nterms)
C       lmpn    leading dim for marray2 (must exceed nterms)
c       w     :  work array 
c       lw    :  length of work array 
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray2  coefficients of rotated expansion.
c       lused    amount of workspace used.
c       ier      error return flag
c                0 successful execution
c                1 insufficient memory
c
C---------------------------------------------------------------------
c
c
c
      implicit real *8 (a-h,o-z)
      integer nquad
      real *8 w(lw)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
c
c      nquad = 2*nterms+2
c        write(*,*) nquad
      nquad = next235_cproj((2*nterms+2)*1.0d0)
c        write(*,*) nquad
      ictheta = 1
      istheta = ictheta+nquad
      icphi = istheta+nquad
      isphi = icphi+nquad
      iynm = isphi+nquad
      iynmd = iynm + (nterms+1)**2
      irat1 = iynmd + (nterms+1)**2
      irat2 = irat1 + (nterms+1)**2
      iuval = irat2 + (nterms+1)**2
      iuder = iuval + 2*nquad*(nterms+1)
      iephi = iuder + 2*nquad*(nterms+1)
      iwsave = iephi + 2*(2*nterms+1)
      iavec = iwsave + 4*nquad+20
      ibvec = iavec + 2*nquad
      lused = ibvec + 2*nquad
      if (lused.gt.lw) stop
c
      call rotviaproj20(beta,alpha,nquad,nterms,m1,m2,mpole,lmp,
     1           marray2,lmpn,w(ictheta),w(istheta),
     1           w(icphi),w(isphi),w(iynm),w(iynmd),
     1           w(irat1),w(irat2),w(iuval),w(iuder),
     1           w(iephi),w(iwsave),w(iavec),w(ibvec))
      return
      end
c
c
c***********************************************************************
      subroutine rotviaproj20(beta,alpha,nquad,nterms,m1,m2,mpole,lmp,
     1           marray2,lmpn,cthetas,sthetas,cphis,sphis,ynm,ynmd,
     1           rat1,rat2,uval,uder,ephis,wsave,avec,bvec)
c***********************************************************************
C
c       INPUT:
c
c       beta:  the rotation angle about the y-axis.
c       alpha:  the rotation angle about the z-axis.
c       nquad:  number of quadrature points on equator
c       nterms: order of multipole expansion
c
c       m1    : max m index for first expansion   
c               NOT IMPLEMENTED but integer argument must be supplied.
c       m2    : max m index for second expansion
c               NOT IMPLEMENTED but integer argument must be supplied.
c
c       mpole:   coefficients of original multiple expansion
c       lmp:     leading dimension of mpole
c       lmpn:    leading dimension of output array marray2
c       cthetas: workspace of dimension nquad
c       sthetas: workspace of dimension nquad
c       cphis:   workspace of dimension nquad
c       sphis:   workspace of dimension nquad
c       ynm:     workspace for spherical harmonics
c       ynmd:    workspace for theta derivative of spherical harmonics
c       rat1:    workspace of same dimension as ynm for precomputation
c       rat2:    workspace of same dimension as ynm for precomputation
c       uval:    workspace 
c       uder:    workspace 
c       ephis:   workspace for exp(i m phi)
c       wsave:   workspace 
c       avec:    workspace 
c       bvec:    workspace 
c
c       OUTPUT:
c
c       marray2  coefficients of rotated expansion.
c---------------------------------------------------------------------
c
      implicit real *8 (a-h,o-z)
      integer nquad,nterms
      real *8 cthetas(nquad),cphis(nquad)
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
      complex *16 avec(nquad)
      complex *16 bvec(nquad)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
      complex *16 uder(nquad,0:nterms),uval(nquad,0:nterms)
      complex *16 ephis(-nterms:nterms)
      real *8 wsave(4*nquad+20)
c
c     Algorithm:
c     1) get locations of quadrature nodes
c     2) evaluate u and du/dtheta
c     3) project onto spherical harmonics.
c
      call getmeridian(beta,nquad,cthetas,sthetas,cphis,sphis)    
      call evalall2a(beta,alpha,nquad,cthetas,sthetas,cphis,sphis,mpole,
     2           lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
      call projectonynm2(nquad,uval,uder,ynm,ynmd,marray2,lmpn,nterms,
     2           m2,wsave,avec,bvec,rat1,rat2)
      return
      end
C
C
C***********************************************************************
      subroutine evalall0a(beta,alpha,nquad,cthetas,sthetas,cphis,sphis,
     1           mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     beta    : angle of rotation about y-axis.
C     alpha   : angle of rotation about z-axis.
C     nquad    : number of target point son unit sphere
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     mpole    : original multipole expansion
C     nterms   : order of multipole expansion
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     ephis    : work array for exp(i m phi) values
C     rat1     : work array for accelerating ynm calculation.
C     rat2     : work array for accelerating ynm calculation.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     uval(i,j) : contribution to potential 
C                 of multipole terms of order j at ith quad node.
C     uder(i,j) : contributions to theta derivative of potential
C                 of multipole terms of order j at ith quad node.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer ndeg,morder, nquad
      real *8 cthetas(nquad),cphis(nquad)
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 ephi1,ephis(-nterms:nterms)
      complex *16 uder(nquad,0:nterms),uval(nquad,0:nterms)
      complex *16 uv,utheta,uphi,ztmp1,ztmp2,ztsum,ztdif
      complex *16 ux,uy,uz,ima
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
C
      data ima/(0.0d0,1.0d0)/
      pi = 4.0d0*datan(1.0d0)
C
      cbeta = cos(beta)
      sbeta = -sin(beta)
      call ylgndrini(nterms,rat1,rat2)
c
      do jj=1,nquad/2
	 ctheta = cthetas(jj)
	 stheta = sthetas(jj)
	 cphi = cphis(jj)
	 sphi = sphis(jj)
         dir1 = -sbeta
         dir2 = 0
         dir3 = cbeta
         tang1 = cphi*ctheta
         tang2 = sphi*ctheta
         tang3 = -stheta
         proj2 = tang1*dir1 + tang2*dir2 + tang3*dir3
         tang1 = -sphi
         tang2 = cphi
         tang3 = 0
         proj1 = tang1*dir1 + tang2*dir2 + tang3*dir3
	 call ylgndru2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)
         ephi1 = dcmplx(cphis(jj),sphis(jj)) *exp(-ima*alpha)
	 ephis(1) = ephi1
	 ephis(-1) = dconjg(ephi1)
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	    ephis(-i) = dconjg(ephis(i))
	 enddo
c
	 do ndeg = 0,nterms
	    uv=0
	    utheta=0
	    uphi=0
	    do morder = 1,ndeg
               ztmp1=ephis(morder)*mpole(ndeg,morder)
               ztmp2=ephis(-morder)*mpole(ndeg,-morder)
	       ztsum=ztmp1+ztmp2
	       ztdif=ztmp1-ztmp2
	       uv=uv+ynm(ndeg,morder)*ztsum
	       utheta=utheta+ynmd(ndeg,morder)*ztsum
	       uphi=uphi-ynm(ndeg,morder)*morder*ztdif
	    enddo
	    uv=stheta*uv+ynm(ndeg,0)*mpole(ndeg,0)
	    utheta=utheta+ynmd(ndeg,0)*stheta*mpole(ndeg,0)
c
c       ... apply the periodizing operator
c
            uval(jj,ndeg) = uv
            uder(jj,ndeg) = (utheta*proj2+uphi*ima*proj1)
            if( mod(ndeg,2) .eq. 0 ) then
            uval(jj+nquad/2,ndeg) = +uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = -uder(jj,ndeg)
            endif
            if( mod(ndeg,2) .eq. 1 ) then
            uval(jj+nquad/2,ndeg) = -uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = +uder(jj,ndeg)
            endif
c
c       ... alternative form of the periodizing operator
c
cc            uval(jj,ndeg) = 2*uv
cc            uder(jj,ndeg) = 2*(utheta*proj2+uphi*ima*proj1)
cc            uval(jj+nquad/2,ndeg) = 0
cc            uder(jj+nquad/2,ndeg) = 0
	 enddo
      enddo
      return
      end
C
C
C
C
C
C
C***********************************************************************
      subroutine evalall1a(beta,alpha,nquad,cthetas,sthetas,cphis,sphis,
     1           mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     beta    : angle of rotation about y-axis.
C     alpha   : angle of rotation about z-axis.
C     nquad    : number of target point son unit sphere
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     mpole    : original multipole expansion
C     nterms   : order of multipole expansion
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     ephis    : work array for exp(i m phi) values
C     rat1     : work array for accelerating ynm calculation.
C     rat2     : work array for accelerating ynm calculation.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     uval(i,j) : contribution to potential 
C                 of multipole terms of order j at ith quad node.
C     uder(i,j) : contributions to theta derivative of potential
C                 of multipole terms of order j at ith quad node.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer ndeg,morder, nquad
      real *8 cthetas(nquad),cphis(nquad)
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 ephi1,ephis(-nterms:nterms)
      complex *16 uder(nquad,0:nterms),uval(nquad,0:nterms)
      complex *16 uv,utheta,uphi,ztmp1,ztmp2,ztsum,ztdif
      complex *16 ux,uy,uz,ima
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
C
      data ima/(0.0d0,1.0d0)/
      pi = 4.0d0*datan(1.0d0)
C
      cbeta = cos(beta)
      sbeta = -sin(beta)
      call ylgndrini(nterms,rat1,rat2)
c
      do jj=1,nquad/2
	 ctheta = cthetas(jj)
	 stheta = sthetas(jj)
	 cphi = cphis(jj)
	 sphi = sphis(jj)
ccc         write(*,*) ctheta,cphi,sphi
         dir1 = -sbeta
         dir2 = 0
         dir3 = cbeta
         tang1 = cphi*ctheta
         tang2 = sphi*ctheta
         tang3 = -stheta
         proj2 = tang1*dir1 + tang2*dir2 + tang3*dir3
         tang1 = -sphi
         tang2 = cphi
         tang3 = 0
         proj1 = tang1*dir1 + tang2*dir2 + tang3*dir3
	 call ylgndru2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)
         ephi1 = dcmplx(cphis(jj),sphis(jj)) *exp(-ima*alpha)
	 ephis(1) = ephi1
	 ephis(-1) = dconjg(ephi1)
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	    ephis(-i) = dconjg(ephis(i))
	 enddo
ccc         call prin2('ephis=*',ephis(-nterms),2*(2*nterms+1))
c
	 do ndeg = 0,nterms
	    uv=0
	    utheta=0
	    uphi=0
	    do morder = 1,ndeg
               ztmp1=ephis(morder)*mpole(ndeg,morder)
               ztmp2=ephis(-morder)*mpole(ndeg,-morder)
	       ztsum=ztmp1+ztmp2
	       ztdif=ztmp1-ztmp2
	       uv=uv+ynm(ndeg,morder)*ztsum
	       utheta=utheta+ynmd(ndeg,morder)*ztsum
	       uphi=uphi-ynm(ndeg,morder)*morder*ztdif
	    enddo
	    uv=stheta*uv+ynm(ndeg,0)*mpole(ndeg,0)
	    utheta=utheta+ynmd(ndeg,0)*stheta*mpole(ndeg,0)
c
c       ... apply the periodizing operator
c
            uval(jj,ndeg) = uv
            uder(jj,ndeg) = (utheta*proj2+uphi*ima*proj1)
            if( mod(ndeg,2) .eq. 0 ) then
            uval(jj+nquad/2,ndeg) = +uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = -uder(jj,ndeg)
            endif
            if( mod(ndeg,2) .eq. 1 ) then
            uval(jj+nquad/2,ndeg) = -uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = +uder(jj,ndeg)
            endif
c
c       ... alternative form of the periodizing operator
c
cc            uval(jj,ndeg) = 2*uv
cc            uder(jj,ndeg) = 2*(utheta*proj2+uphi*ima*proj1)
cc            uval(jj+nquad/2,ndeg) = 0
cc            uder(jj+nquad/2,ndeg) = 0
	 enddo
      enddo
      return
      end
C
C
C***********************************************************************
      subroutine evalall2a(beta,alpha,nquad,cthetas,sthetas,cphis,sphis,
     1           mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     beta    : angle of rotation about y-axis.
C     alpha   : angle of rotation about z-axis.
C     nquad    : number of target point son unit sphere
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     mpole    : original multipole expansion
C     nterms   : order of multipole expansion
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     ephis    : work array for exp(i m phi) values
C     rat1     : work array for accelerating ynm calculation.
C     rat2     : work array for accelerating ynm calculation.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     uval(i,j) : contribution to potential 
C                 of multipole terms of order j at ith quad node.
C     uder(i,j) : contributions to theta derivative of potential
C                 of multipole terms of order j at ith quad node.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer ndeg,morder, nquad
      real *8 cthetas(nquad),cphis(nquad)
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 ephi1,ephis(-nterms:nterms)
      complex *16 uder(nquad,0:nterms),uval(nquad,0:nterms)
      complex *16 uv,utheta,uphi,ztmp1,ztmp2,ztsum,ztdif
      complex *16 ux,uy,uz,ima
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
C
      data ima/(0.0d0,1.0d0)/
      pi = 4.0d0*datan(1.0d0)
C
      cbeta = cos(beta)
      sbeta = -sin(beta)
      call ylgndrini(nterms,rat1,rat2)
c
ccc        write(*,*) nquad, nquad/2, nquad/4
      nquad2=nquad/2
      if( mod(nquad2,2) .eq. 0 ) nquad4=nquad2/2+1
      if( mod(nquad2,2) .eq. 1 ) nquad4=nquad2/2+1
c
      do jj=1,nquad4
	 ctheta = cthetas(jj)
	 stheta = sthetas(jj)
	 cphi = cphis(jj)
	 sphi = sphis(jj)
ccc         write(*,*) ctheta,cphi,sphi
         dir1 = -sbeta
         dir2 = 0
         dir3 = cbeta
         tang1 = cphi*ctheta
         tang2 = sphi*ctheta
         tang3 = -stheta
         proj2 = tang1*dir1 + tang2*dir2 + tang3*dir3
         tang1 = -sphi
         tang2 = cphi
         tang3 = 0
         proj1 = tang1*dir1 + tang2*dir2 + tang3*dir3
	 call ylgndru2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)
         ephi1 = dcmplx(cphis(jj),sphis(jj)) *exp(-ima*alpha)
	 ephis(1) = ephi1
	 ephis(-1) = dconjg(ephi1)
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	    ephis(-i) = dconjg(ephis(i))
	 enddo
ccc         call prin2('ephis=*',ephis(-nterms),2*(2*nterms+1))
c
	 do ndeg = 0,nterms
	    uv=0
	    utheta=0
	    uphi=0
	    do morder = 1,ndeg
               ztmp1=ephis(morder)*mpole(ndeg,morder)
               ztmp2=ephis(-morder)*mpole(ndeg,-morder)
	       ztsum=ztmp1+ztmp2
	       ztdif=ztmp1-ztmp2
	       uv=uv+ynm(ndeg,morder)*ztsum
	       utheta=utheta+ynmd(ndeg,morder)*ztsum
	       uphi=uphi-ynm(ndeg,morder)*morder*ztdif
	    enddo
	    uv=stheta*uv+ynm(ndeg,0)*mpole(ndeg,0)
	    utheta=utheta+ynmd(ndeg,0)*stheta*mpole(ndeg,0)
c
c       ... apply the periodizing operator
c
            uval(jj,ndeg) = uv
            uder(jj,ndeg) = (utheta*proj2+uphi*ima*proj1)
            if( mod(ndeg,2) .eq. 0 ) then
            uval(jj+nquad/2,ndeg) = +uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = -uder(jj,ndeg)
            endif
            if( mod(ndeg,2) .eq. 1 ) then
            uval(jj+nquad/2,ndeg) = -uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = +uder(jj,ndeg)
            endif
c
c       ... alternative form of the periodizing operator
c
cc            uval(jj,ndeg) = 2*uv
cc            uder(jj,ndeg) = 2*(utheta*proj2+uphi*ima*proj1)
cc            uval(jj+nquad/2,ndeg) = 0
cc            uder(jj+nquad/2,ndeg) = 0
	 enddo

         if_reflect=1
         if( jj .eq. 1 ) if_reflect=0
         if( jj .eq. nquad4 .and. mod(nquad2,2) .eq. 0 ) if_reflect=0
c
         if( if_reflect .eq. 1 ) then
         jjr=nquad2 - jj + 2
ccc         write(*,*) jj, jjr
         call ylgndr2pm_opt(nterms,ynm,ynmd)

	 ctheta = -cthetas(jj)
	 stheta = sthetas(jj)
	 cphi = -cphis(jj)
	 sphi = sphis(jj)
ccc         write(*,*) ctheta,cphi,sphi
         dir1 = -sbeta
         dir2 = 0
         dir3 = cbeta
         tang1 = cphi*ctheta
         tang2 = sphi*ctheta
         tang3 = -stheta
         proj2 = tang1*dir1 + tang2*dir2 + tang3*dir3
         tang1 = -sphi
         tang2 = cphi
         tang3 = 0
         proj1 = tang1*dir1 + tang2*dir2 + tang3*dir3

         ephi1 = dcmplx(cphi,sphi) *exp(-ima*alpha)
	 ephis(1) = ephi1
	 ephis(-1) = dconjg(ephi1)
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	    ephis(-i) = dconjg(ephis(i))
	 enddo
ccc         call prin2('ephis=*',ephis(-nterms),2*(2*nterms+1))
c
	 do ndeg = 0,nterms
	    uv=0
	    utheta=0
	    uphi=0
	    do morder = 1,ndeg
               ztmp1=ephis(morder)*mpole(ndeg,morder)
               ztmp2=ephis(-morder)*mpole(ndeg,-morder)
	       ztsum=ztmp1+ztmp2
	       ztdif=ztmp1-ztmp2
	       uv=uv+ynm(ndeg,morder)*ztsum
	       utheta=utheta+ynmd(ndeg,morder)*ztsum
	       uphi=uphi-ynm(ndeg,morder)*morder*ztdif
	    enddo
	    uv=stheta*uv+ynm(ndeg,0)*mpole(ndeg,0)
	    utheta=utheta+ynmd(ndeg,0)*stheta*mpole(ndeg,0)
c
c       ... apply the periodizing operator
c
            uval(jjr,ndeg) = uv
            uder(jjr,ndeg) = (utheta*proj2+uphi*ima*proj1)
            if( mod(ndeg,2) .eq. 0 ) then
            uval(jjr+nquad/2,ndeg) = +uval(jjr,ndeg)
            uder(jjr+nquad/2,ndeg) = -uder(jjr,ndeg)
            endif
            if( mod(ndeg,2) .eq. 1 ) then
            uval(jjr+nquad/2,ndeg) = -uval(jjr,ndeg)
            uder(jjr+nquad/2,ndeg) = +uder(jjr,ndeg)
            endif
c
c       ... alternative form of the periodizing operator
c
cc            uval(jjr,ndeg) = 2*uv
cc            uder(jjr,ndeg) = 2*(utheta*proj2+uphi*ima*proj1)
cc            uval(jjr+nquad/2,ndeg) = 0
cc            uder(jjr+nquad/2,ndeg) = 0
	 enddo
         endif
      enddo
      return
      end
C
C
C
C
C
c***********************************************************************
      subroutine rotviaproj2f90(beta,alpha,nterms,m1,m2,mpole,lmp,
     1           marray2,lmpn)
c***********************************************************************
c       Purpose:
c
c	Fast and stable algorithm for applying rotation operator about
c	the z-axis determined by angle alpha plus the y-axis determined
c	by angle beta. After rotation, the expansion pole is moved to
c	location (beta, alpha) in spherical coordinates (theta, phi).
c
c       The method is based on computing the induced potential and
c       its theta-derivative on the rotated equator
c       for each order (first index). The coefficients of  the rotated
c       expansion can then be obtained by FFT and projection.
c
c       There is some loss in speed over using recurrence relations 
c       but it is stable to all orders whereas the recurrence schemes 
c       are not.
c
c       If the rotation operator is to be used multiple times, and
c       memory is available, one can precompute and store the 
c       multipliers used in evalall (see below). This has not yet been
c       implemented.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       beta:  the rotation angle about the y-axis.
c       alpha:  the rotation angle about the z-axis.
c       nterms: order of multipole expansion
C       mpole   coefficients of original multiple expansion
C       lmp     leading dim for mpole (must exceed nterms)
C       lmpn    leading dim for marray2 (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray2  coefficients of rotated expansion.
c
C---------------------------------------------------------------------
c
c
c
      implicit none
      integer nquad,ier,m1,m2,nterms,lmp,lmpn, next235_cproj
      integer ictheta,istheta,icphi,isphi,iynm,iynmd,irat1,irat2
      integer iuval,iuder,iephi,iwsave,iavec,ibvec,lused
      real *8 beta,alpha
      real *8, allocatable :: w(:)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
      complex *16, allocatable :: cw(:)
c
c      nquad = 2*nterms+2
c        write(*,*) nquad
      nquad = next235_cproj((2*nterms+2)*1.0d0)
c        write(*,*) nquad
      ictheta = 1
      istheta = ictheta+nquad
      icphi = istheta+nquad
      isphi = icphi+nquad
      iynm = isphi+nquad
      iynmd = iynm + (nterms+1)**2
      irat1 = iynmd + (nterms+1)**2
      irat2 = irat1 + (nterms+1)**2
      iwsave = irat2 + (nterms+1)**2
      lused = iwsave + 4*nquad+20
      allocate (w(lused), stat=ier)
      iuval = 1
      iuder = iuval + nquad*(nterms+1)
      iephi = iuder + nquad*(nterms+1)
      iavec = iephi + (2*nterms+1)
      ibvec = iavec + 2*nquad
      lused = ibvec + 2*nquad
      allocate (cw(lused), stat=ier)
      if (ier.ne.0) then 
         write(6,*) ' alloc failure in rotviaproj2f90'
         stop
      endif
c
      call rotviaproj20(beta,alpha,nquad,nterms,nterms,nterms,
     1           mpole,lmp,marray2,lmpn,w(ictheta),w(istheta),
     1           w(icphi),w(isphi),w(iynm),w(iynmd),
     1           w(irat1),w(irat2),cw(iuval),cw(iuder),
     1           cw(iephi),w(iwsave),cw(iavec),cw(ibvec))
      deallocate(w)
      return
      end
c
c
