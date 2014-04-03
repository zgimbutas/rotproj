        implicit real *8 (a-h,o-z)
c
        real *8 sqc(20 000)
        real *8 rd1(10 000 000)
        real *8 rd2(10 000 000)
        complex *16 mpole(10 000 000)
        complex *16 marray(10 000 000)
        dimension w(10 000 000)
        dimension errs(100 000)
c
        lw=10 000 000
c
c       initilialize print routines to unit numbers
c       6 and 13. This prints to screen and to fort.13.
c
        call prini(6,13)
c
c       set expansion length nterms
c
        nterms=100
c
c       create mutipole expansion
c
        call mpoleinit(nterms,mpole)
c
c
c       set rotation angle alpha
c
        done=1
        pi=4*atan(done)
        alpha=pi/3
c
        ntimes=10
c
c       carry out rotation
c
        call prinf('proj_test*',i,0)
	call prinf('nterms=*',nterms,1)
	call prin2('angle=*',alpha,1)
c
        t1=second()
        do k=1,ntimes
        call rotviaproj(alpha,nterms,nterms,nterms,mpole,nterms,
     $     marray,nterms,w,lw,lused)
        enddo
        t2=second()
        call prin2('time/call=*',(t2-t1)/ntimes,1)
        call prin2('speed=*',ntimes/(t2-t1),1)
c
ccc        call prinm1(marray,nterms,nterms2)
c
c       check precision by comparing norm of original and rotated
c       expansions.
c
        call mpolecheck(nterms,mpole,marray,errs)
        call prin2('errors=*',errs,nterms+1)
c
        stop
        end
c
c
c
c
c       
        subroutine mpoleinit(nterms,mpole)
        implicit real *8 (a-h,o-z)
        complex *16 mpole(0:nterms,-nterms:nterms)
c       
c       initialize multipole expansion
c
        itype=3
c
        if( itype .eq. 1 ) then
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        mpole(n,0)=1
        enddo
        endif
c
        if( itype .eq. 2 ) then
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=n+m
        enddo
        mpole(n,0)=1
        enddo
        endif
c
        if( itype .eq. 3 ) then
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=1/sqrt(2*n+1.0d0)
        enddo
        enddo
        endif
c
        return
        end
c
c
c
c
c
        subroutine mpolecheck(nterms,mpole,marray,errs)
        implicit real *8 (a-h,o-z)
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 marray(0:nterms,-nterms:nterms)
        real *8 errs(0:nterms)
c       
c       check that norm is preserved under rotation
c
        do n=0,nterms
        d1=0
        d2=0
        do m=-n,n
        d1=d1+mpole(n,m)*dconjg(mpole(n,m))
        d2=d2+marray(n,m)*dconjg(marray(n,m))
        enddo
        d1=sqrt(d1/(2*n+1.0d0))
        d2=sqrt(d2/(2*n+1.0d0))
        errs(n)=d1-d2
        enddo
c
        return
        end
c
c
c
C
C
C
C
        SUBROUTINE PRINI(IP1,IQ1)
        CHARACTER *1 MES(1), AA(1)
         save
        REAL *4 A(1)
        REAL *8 A2(1)
ccc        REAL *16 A4(1)
        INTEGER *4 IA(1)
        INTEGER *2 IA2(1)
        IP=IP1
        IQ=IQ1
        RETURN

C
C
C
C
C
        ENTRY PRIN(MES,A,N)
        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
         RETURN
C
C
C
C
        ENTRY PRIN2(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C
C
C
C
ccc        ENTRY PRINQ(MES,A4,N)
ccc        CALL MESSPR(MES,IP,IQ)
ccc        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
ccc        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
ccc 1500 FORMAT(6(2X,e11.5))
ccc        RETURN
C
C
C
C
        ENTRY PRINF(MES,IA,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C
C
C
C
        ENTRY PRINF2(MES,IA2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY PRINA(MES,AA,N)
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END
c
c
c
c
c
        SUBROUTINE MESSPR(MES,IP,IQ)
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END
C
C
C
C
C
        SUBROUTINE ZTIME(I)
        J=1
        J=7-I+J
CCCC    I=MRUN(J)
        RETURN
        END

C
C
      SUBROUTINE PRINM(MPOLE,NTERMS)
      implicit real *8 (a-h,o-z)
      real *8 MPOLE0(0:NTERMS,0:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,-nterms:NTERMS)
      INTEGER NTERMS
C
C     print out coefficients of multipole expansion
C
1000  FORMAT(6E12.5)
1001  FORMAT(/)
      DO 100 L = 0,NTERMS
         WRITE(6,1000)(MPOLE(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
100   CONTINUE
        return
C
C
C
C
      ENTRY PRINM0(MPOLE0,NTERMS)
      DO 200 L = 0,NTERMS
         WRITE(6,1000)(MPOLE0(L,M),M=0,L)
         WRITE(13,1000)(MPOLE0(L,M),M=0,L)
         WRITE(6,1001)
         WRITE(13,1001)
200   CONTINUE
c
      RETURN
c
c
c
      ENTRY PRINM1(MPOLE,NTERMS,NTERMS2)
      DO 300 L = 0,NTERMS2
         WRITE(6,1000)(MPOLE(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
300   CONTINUE
c
      RETURN
        end
c
c
