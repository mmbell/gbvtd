C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE GOODCIRC(AZ,VRERR,IRAYPROC,IGAP,FLAG,numcoeff)

C  This subroutine checks whether there is an unacceptable data gap in
C  the ring.  Flag is returned .TRUE. if and only if the ring is good.

      IMPLICIT NONE
      INTEGER IRAYPROC,I,J,JJ,IGAP(5),ISUM
      REAL AZ(IRAYPROC),vrerr(irayproc),pi,dtor
      LOGICAL DEGREE_SECTOR(360),FLAG
      INTEGER NUMCOEFF,SAVENUMCOEFF,scale

C      print *,"in subroutine goodcirc"
      NUMCOEFF = 9
      savenumcoeff=numcoeff
      pi=acos(-1.)
      dtor=pi/180.
c      scale = 360/irayproc
      DO 10 J=1,360
         DEGREE_SECTOR(J)=.FALSE.
10    CONTINUE

      DO 20 I=1,IRAYPROC
         IF(VRERR(I).LE.90. .and. vrerr(i).gt.-90.)THEN
            J=IFIX(AZ(I)/dtor)+1
            if(j.gt.360) then
              j=j-360
              print *, 'j > 360, j=',j
            endif
            DEGREE_SECTOR(J)=.TRUE.
c            print *, i, j, az(i)/dtor, vrerr(i), degree_sector(j)
         ENDIF
20    CONTINUE
 
      do i=1,5
        IF(IGAP(I).GE.360)THEN
           FLAG=.TRUE.
           RETURN
        ENDIF
      enddo

      ISUM=0
      DO 30 JJ=1,720
         J=MOD(JJ-1,360)+1
         IF(DEGREE_SECTOR(J))THEN
            ISUM=0
            IF(JJ.GE.360)THEN
               FLAG=.TRUE.
C               PRINT *,'NUMCOEFF',NUMCOEFF
               RETURN
            ENDIF
         ELSE
            ISUM=ISUM+1
            IF(ISUM.GT.IGAP(1))THEN
               if(isum.gt.igap(2)) then
                  if(isum.gt.igap(3)) then
                    if(isum.gt.igap(4)) then
                       FLAG=.FALSE.
                       RETURN
                    else
                       savenumcoeff=3
                    endif
                 else
                    savenumcoeff=5
                 endif
              else
                 savenumcoeff=7
              ENDIF
            else
               savenumcoeff=9
            endif
         ENDIF
         IF(NUMCOEFF .GT. SAVENUMCOEFF) THEN
            NUMCOEFF = SAVENUMCOEFF
         ENDIF
c         print *,'jj =',jj,'isum =',isum,'numcoeff',numcoeff
30    CONTINUE

      END

*****************************************************************
*   FINAL1, This is for 9 coefficients (wave 0-3)
*****************************************************************

        subroutine final1(aa,vtcoef,vrc0,V_Mcos,xc,yc,ip,numcoeff,
     X                    flag1,flagwrite)

        implicit none
       INCLUDE 'gbvtd_parameter.f'
        real vrc0,V_Mcos,xc,yc,pi,R_T,cp,sp,a0,b0,a1,b1,a2,b2,a3,b3
        real amvt,a4,b4,vt0,vt1,vt2,vt3
        integer ip,numcoeff,flag1,flagwrite
	real vtcoef(maxcoeffs),aa(maxcoeffs)

	pi=acos(-1.)
        R_T=sqrt(xc*xc+yc*yc)
        cp=sqrt(R_T*R_T-ip*ip)/R_T  !cos(\alpha_max)
        sp=float(ip)/R_T            !sin(\alpha_max)
        A0=0.
        B0=0.
        A1=0.
        B1=0.
        A2=0.
        B2=0.
        A3=0.
        B3=0.
        A4=0.
        B4=0.
        if(numcoeff.eq.3 .or. sp .gt. 0.8) then
	  A0=aa(1)
	  B0=0.
	  A1=aa(3)
	  B1=aa(2)
        elseif(numcoeff.eq.5 .or. sp .gt. 0.6) then
	  A0=aa(1)
	  B0=0.
	  A1=aa(3)
	  B1=aa(2)
	  A2=aa(5)
	  B2=aa(4)
        elseif(numcoeff.eq.7 .or. sp .gt. 0.4) then
	  A0=aa(1)
	  B0=0.
	  A1=aa(3)
	  B1=aa(2)
	  A2=aa(5)
	  B2=aa(4)
	  A3=aa(7)
	  B3=aa(6)
        elseif(numcoeff.eq.9) then
	  A0=aa(1)
	  B0=0.
	  A1=aa(3)
	  B1=aa(2)
	  A2=aa(5)
	  B2=aa(4)
	  A3=aa(7)
	  B3=aa(6)
          A4=aa(9)
          B4=aa(8)
        else
          print *, 'too few points'
        endif
        if(flag1.eq.1 .or. flag1.eq.3) then       !VDAD
	  vtcoef(2)=0.                                !VTS0
	  vtcoef(1)=-B1-B3                            !VTC0
	  vtcoef(3)=-2.*(B2+B4)                       !VTC1
	  vtcoef(4)=2.*(A2+A4)                        !VTS1
	  vtcoef(5)=-2.*B3                           !VTC2
	  vtcoef(6)=2.*A3                             !VTS2
          vtcoef(7)=-2.*B4                            !VTC3
          vtcoef(8)=2.*A4                             !VTS3
	  if (sp.ge.0.99)then
	    V_Mcos=-999.
	    vrc0=-999.
	  else
	    vrc0=((A1+A3)*R_T-(a0+a2)*float(ip))/cp
	    V_Mcos=((A0+A2)*R_T-(a1+a3)*float(ip))/cp
	    amvt=vtcoef(1)
	  endif
        elseif(flag1.eq.2 .or. flag1.eq.4) then    !GBVTD
	  vtcoef(2)=0.                                !VTS0
C -2.*sqrt(1-(A0+A2)**2.)*sp has been removed from the following equation
C  because VM is unknown
CC old code before changes at 1/12/98
cc	  vtcoef(1)=-B1-B3                            !VTC0
cc	  vtcoef(4)=A2-A0+(A0+A2)*cp+2.*A4            !VTS1
CC modified by wcl on 1/12/98
cc	  vtcoef(1)=-B1+B3                            !VTC0
cc modified again by wcl on 12/19/00, the 1/12/98 correction on B3 was
cc   WRONG                            !VTC0
	  vtcoef(1)=-B1-B3                            !VTC0
	  vtcoef(4)=A2-A0+A4+(A0+A2+A4)*cp            !VTS1
	  vtcoef(3)=-2.*(B2+B4)                       !VTC1
	  vtcoef(5)=-2.*B3                            !VTC2
	  vtcoef(6)=2.*A3                             !VTS2
          vtcoef(7)=-2.*B4                            !VTC3
          vtcoef(8)=2.*A4                             !VTS3
	  if (sp.ge.0.99)then
	    V_Mcos=-999.
	    vrc0=-999.
	  else
	    vrc0=A1+A3
	    V_Mcos=A0+A2+A4
	    amvt=vtcoef(1)
            if(abs(vrc0).gt.100.) vrc0=-999.0
            if(abs(V_Mcos).gt.100.) V_Mcos=-999.0 
	  endif
        else
          print *, "flag1 out of range", flag1
          stop
        endif
        vt1=sqrt(vtcoef(4)**2+vtcoef(3)**2)
        vt2=sqrt(vtcoef(6)**2+vtcoef(5)**2)
        vt3=sqrt(vtcoef(8)**2+vtcoef(7)**2)
        vt0=vtcoef(1)
        if (vt0.lt.vt1.or.vt0.lt.vt2.or.vt0.lt.vt3) then
           numcoeff=1
           print *,'VT Error! ',vtcoef
        endif
        if(flagwrite.eq.1) then    !write coefficients
        write(*,'(1x,a8,i5,a4,2x,a15,f8.2,2x,a17,f8.2))')'RADIUS= ',ip, 
     X         ' KM ','sin(alphamax)= ', sp, 'cos(alphamax)=', cp
        write(*,'(7f8.2)')aa 
        write(*,'(2(1x,a7,f6.1))')'VTC0=',vtcoef(1), ' VTS0 =',vtcoef(2)
        write(*,'(2(1x,a7,f6.1))')'VTC1=',vtcoef(3), ' VTS1 =',vtcoef(4)
        write(*,'(2(1x,a7,f6.1))')'VTC2=',vtcoef(5), ' VTS2 =',vtcoef(6)
        write(*,'(2(1x,a7,f6.1))')'VTC3=',vtcoef(7), ' VTS3 =',vtcoef(8)
        write(*,'(2(1x,a7,f6.1))')'VRC0=',vrc0, ' MLP  =',V_Mcos
        endif
	return
	end

*************************************************************
*  CURFIT
*************************************************************

         SUBROUTINE CURFIT(C,DDZ,A,thetamax)
         implicit none
         integer i,thetamax
         real pi,th
	 REAl C(7),A(thetamax),DDZ(thetamax)
	 PI=ACOS(-1.)
	 DO I=1,thetamax
	 TH=(I-1)*2.*PI/89.
	 A(I)=DDZ(I)
	 IF (INT(DDZ(I)).EQ.99)THEN
	 A(I)=C(1)+C(3)*COS(TH)+C(2)*SIN(TH)
         ENDIF
         ENDDO
	 RETURN
	 END

*************************************************************
*   GOODDATA 
************************************************************* 
	SUBROUTINE GOODDATA(DDZ,TH,ADDZ,ATH,N,thetamax,special)
        implicit none
        integer k,n,thetamax
	REAL DDZ(thetamax),ATH(thetamax),TH(thetamax),ADDZ(thetamax)
        real special
	N=0
	DO K=1,thetamax
c	  IF (INT(DDZ(K)).NE.99)THEN
          if(abs(ddz(k)-special).gt.0.001)then
	    N=N+1
            ATH(N)=TH(K)
            ADDZ(N)=DDZ(K)
c            print *,N,ath(N)*180/3.1415926,addz(N)
          ENDIF
	ENDDO
        RETURN
	END
***********************************************************
*   THETA2SI
***********************************************************

      SUBROUTINE theta2si(ip,th,xc,yc,si,theta_t)

      implicit none
      real th,xc,yc,xx,yy,xy,pi,afi,si,theta_t
      integer ip

      pi=acos(-1.)
      xx=xc+float(ip)*cos(th+theta_t) 
      yy=yc+float(ip)*sin(th+theta_t) 
      afi=atan2(yy,xx)-theta_t  
      si=th-afi
      call fixangle(si)
c      print *, ip, theta_t,afi/pi*180., th/pi*180., si/pi*180.
      return 
      end
	
***********************************************************
*   SI2THETA th is thetaprime, not theta in math coordinate
***********************************************************

      SUBROUTINE si2theta(ip,si,xc,yc,th,theta_t)

      implicit none
      real th,xc,yc,rt,pi,afi,si,theta_t
      integer ip

      pi=acos(-1.)
      rt=sqrt(xc*xc+yc*yc)
      afi=asin(float(ip)/rt*sin(si))
      th=si+afi
      call fixangle(th)
c      print *, ip, theta_t,afi/pi*180., th/pi*180., si/pi*180.
      return 
      end
	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C STD_DVI_NEW  This routine computes the standard deviation between
C           the data points and the fitted curve, for any number of points
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccCCcc

      Subroutine std_dvi_new(A,C,X,DVD,phi,numcoeff,sd,xc,yc,flag,k,
     X                   theta_t,dummy,NN)

        implicit none
        include 'gbvtd_parameter.f'
        integer i,in,flag,numcoeff,k,j,NN
	REAL dummy,DVD(NN),phi(NN),pi,xc,yc,theta_t
        real A(NN),X(NN),C(maxcoeffs),theta,scale,sd

        SD=0.
	IN=0
        pi=acos(-1.)
        scale=360./float(thetamax)
	DO I=1,NN
          if(flag.eq.2 .or. flag.eq.4) then !GBVTD,needs to convert si to theta
            call si2theta(k,phi(i),xc,yc,theta,theta_t)
            X(i)=theta/pi*180. !transform x from deg to grid index for ploting
C                              !4. is delta theta
          endif
          A(i)=c(1)
          do j=1,numcoeff/2
            print *, i,j,a(i),dvd(i),sd, in,numcoeff
             A(i)=A(i)+c(2*j+1)*cos(float(j)*phi(i))
     *                 +c(2*j)*sin(float(j)*phi(i))
          enddo
c           A(I)=C(1)+C(3)*COS(TH)+C(2)*SIN(TH)+C(5)*COS(2.*TH)
c     *	    +C(4)*SIN(2.*TH)+C(7)*COS(3.*TH)+C(6)*SIN(3.*TH)
	  IF (abs(DVD(I)-dummy).gt.0.001)THEN
	    SD=SD+(A(I)-DVD(I))**2.
	    IN=IN+1
c            print *, i,j,a(i),dvd(i),sd, in
	  ENDIF
        ENDDO
        if(IN.gt.0) then
          SD=SQRT(SD/IN)
c          print *, in,sd
        else
          SD=dummy
        endif
        return
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C STD_DVI  This routine computes the standard deviation between
C           the data points and the fitted curve
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccCCcc

      Subroutine std_dvi(A,C,X,DVD,phi,numcoeff,sd,xc,yc,flag,k,
     X                   theta_t,dummy)

        implicit none
        include 'gbvtd_parameter.f'
	REAL dummy,DVD(thetamax),phi(thetamax),sd,pi,xc,yc,theta_t
        real A(thetamax),X(thetamax),C(maxcoeffs),theta,scale
        integer i,in,flag,numcoeff,k,j

	SD=0.
	IN=0
        pi=acos(-1.)
        scale=360./float(thetamax)
	DO I=1,thetamax
          if(flag.eq.2 .or. flag.eq.4) then !GBVTD,needs to convert si to theta
            call si2theta(k,phi(i),xc,yc,theta,theta_t)
            X(i)=theta/pi*180./scale !transform x from deg to grid index for ploting
C                              !4. is delta theta
          endif
          A(i)=c(1)
          do j=1,numcoeff/2
             A(i)=A(i)+c(2*j+1)*cos(float(j)*phi(i))
     *                 +c(2*j)*sin(float(j)*phi(i))
          enddo
c           A(I)=C(1)+C(3)*COS(TH)+C(2)*SIN(TH)+C(5)*COS(2.*TH)
c     *	    +C(4)*SIN(2.*TH)+C(7)*COS(3.*TH)+C(6)*SIN(3.*TH)
	  IF (abs(DVD(I)-dummy).gt.0.001)THEN
	    SD=SD+(A(I)-DVD(I))**2.
	    IN=IN+1
c            print *, i,a(i),dvd(i),sd, in
	  ENDIF
        ENDDO
        if(IN.gt.0) then
          SD=SQRT(SD/IN)
c          print *, in,sd
        else
           SD=dummy
        endif
        return
        end

*************************************************************************
* Fixangle, put angles back between 0 and 360 degree
*************************************************************************
       Subroutine fixangle(angle)
       implicit none
       real angle, pi

       pi=acos(-1.)
       do while (abs(angle-pi).gt. pi)
          if(angle.gt.0.) then
             angle=angle-2.*pi
          else
             angle=angle+2.*pi
          endif
       enddo
       return
       end
