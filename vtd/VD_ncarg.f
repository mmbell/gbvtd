      Program VD

C ****** OBJECT:TO GET THE VORTEX STRUCTURE USING VDAD METHOD ******      
C ****** NAME: VD.f                                           ******     
C ****** COMPLETED DATE: 1995/11/1                            ******      
C ****** AUTHOR:CHANG BAU_LUNG, WEN-CHAU LEE
C ****** All angles (theta and si) are in GBVTD coordinate,
C ****** i.e., not in math coordiate. 0 deg is TA see lee et al. (1997)
C ****** So, when plotting the results, you need to rotate the
C ****** coordinate systems by theta_t (atan2(yc,xc)).
C ****** VD is the Doppler velocity
C ****** VA is VT plus the component of VM
      implicit none
      INCLUDE 'gbvtd_parameter.f'
      integer N,KZ,IZ,K,nn,i,iz1,numcoeff,kc,l,index,rflag
      real pi,xc(20),yc(20),dummy,vrcoef,vm,s,theta_t,rad_lat,rad_long
      CHARACTER FNAME*80,infile*80,idealfile*80
      character prefix*20
      REAL THETAZ(thetamax),d,rt
      REAL DZ(nx,ny),DDZ(thetamax),XDD(thetamax),A(maxrays)
      Real A1(thetamax),X1(thetamax),max_pos,max_neg
      real YDD(thetamax),THETAV(maxrays),dddz(rmax,thetamax)
      real meanvt(rmax,zmax),meanvr(rmax,zmax),x(maxrays)
      real meandz(rmax,zmax),sd,meansd(rmax,zmax),vdif(rmax,thetamax)
      real vdifp(rmax,thetamax),xc1,yc1
      real xdddz(thetamax), ydddz(thetamax),VT(rmax,thetamax)
      real tempdz(thetamax),vr(rmax,thetamax),va(rmax,thetamax)
      REAL XS(maxcoeffs,maxrays),YS(maxrays),WEIGHT(maxrays)
      real CC(maxcoeffs),STAND_ERR(maxcoeffs),theta
      REAL TH(maxrays),AZ(maxrays),VTCOEF(maxcoeffs),dist
      real V_D(nx,ny),DVD(maxrays),DDVD(rmax,thetamax),dvd1(thetamax)
      real ori_lat,ori_long,temp,si,dtor,theta1,ddz1(thetamax)
      real asym,ampvt1,xdist, ydist,angle,alpha,power,wave,rotate
       real phasevt1,ampvt2,phasevt2,ampvt3,phasevt3,vorticity
       real ampvr1,phasevr1,ampvr2,phasevr2,ampvr3,phasevr3
       real rad_v,vrmean,u,meanu,v,special,rg,meanflow
       real meanangle,meanv,r_t,vta,ar,zone,vt1d(thetamax)
       character outfiledz*80,outfileve*80,inputfile*80
c      character outfilevt*80, outfilerad*80
      INTEGER NDATA,numgap,gap(5),flag1,j,dj,flag2,num_ray
      integer iori,jori,n1,numloc,ring_start, ring_end,j1,n2,n3,xmax
      integer*2 id(510),sf,af,id1(510),id2(510),num_x, num_y
      integer*2 field1,field2,bid(510)
      integer nargs,iargc
      LOGICAL FLAG
      data rad_lat, rad_long /25.0694,121.2083/
      data sf, af /100,64/
       common /ideal/vorticity,ar,meanflow,meanangle,phasevt1,ampvt1,
     X               phasevt2,ampvt2,phasevt3,ampvt3,phasevr1,ampvr1,
     X               phasevr2,ampvr2,phasevr3,ampvr3,vrmean,alpha,wave,
     X               power,rotate
      CALL OPNGKS
      PI=ACOS(-1.)
      DUMMY=-999.
      dtor=pi/180.
      do i=1,510
        id(i)=-999
      enddo

      nargs = iargc()
       IF ( nargs.GT.0 ) THEN
	  Call getarg(1, infile)
       ELSE
          print *, 'Please enter input file name:'
          read *, infile
       ENDIF
      OPEN(12,FILE=infile,STATUS='OLD')
      read(12,111)flag1
      read(12,111)flag2
      read(12,114)field1,field2
      read(12,112)prefix
      read(12,111)IZ1
      read(12,111)KZ        !number of vertical levels
      read(12,111)numgap
      do i=1,numgap
        read(12,111)gap(i)
      enddo
      read(12,111)ring_start
      read(12,111)ring_end
      read(12,111)rflag
      read(12,113)zone
113   format(f8.2)
111   format(i3)
112   format(a20)
 114  format(a2,/,a2)
c      call intswap(field1,field1)
c      call intswap(field2,field2)
C     IZ:total number of tropical cyclone volumes
      DO IZ=1,IZ1
        READ(12,'(a80)',END=999)FNAME,idealfile
        if(flag2.eq.2 .and. (flag1.eq.3 .or. flag1.eq.4)) then
c          print *, 'please enter ideal input filename:
c          read *, inputfile
          OPEN(21,file=idealfile,status='old',form='formatted',
     +               err=998)
c        READ(21,'(A80/A80/A80/A80/f6.1,1x,f6.1)') outfiledz,
c     +          outfileve,outfilevt,outfilerad,xc,yc
          READ(21,'(/,/,/,/,f6.1,1x,f6.1)')xc1,yc1
          read(21,'(f8.1,/,f5.1)') vorticity,ar
c          read(21,'(f8.1,/,f5.1)') alpha,wave
c          read(21,'(f8.1,/,f5.1)') power,rotate
          read(21,'(f8.1,/,f8.1)')meanflow,meanangle
          read(21,'(f8.1,/,f8.1)')phasevt1,ampvt1
          read(21,'(f8.1,/,f8.1)')phasevt2,ampvt2
          read(21,'(f8.1,/,f8.1)')phasevt3,ampvt3
          read(21,'(f8.1)')vrmean
          read(21,'(f8.1,/,f8.1)')phasevr1,ampvr1
          read(21,'(f8.1,/,f8.1)')phasevr2,ampvr2
          read(21,'(f8.1,/,f8.1)')phasevr3,ampvr3
          alpha = 0
          wave = 0
          power = 1
          rotate = 0
        endif
        call openfile(1,numloc,iz,prefix)
        OPEN(10,FILE=fname,STATUS='UNKNOWN',FORM='FORMATTED')
        read(10,116)id1
        do N=1,kz
          READ(12,'(i3,1x,F7.2,1x,F7.2)')n3,XC(N),YC(N)     !read center
          if(n3.ne.N) then
             print *,'center mismatch with level'
             stop
          endif
          read(10,117)N1
117       format(5x,i2)
          write(14,988)xc(N),yc(N)
988       format('xc= ',f8.2, 'yc= ',f8.2)
          write(14,989)N
989       format('height=', i5)
          r_t=sqrt(xc(N)*xc(N)+yc(N)*yc(N))
          if(r_t.le.float(ring_end)) ring_end=int(r_t*1.0)
          iori=id1(309)/sf
          jori=id1(310)/sf
          num_x=id1(162)
          num_y=id1(167)
          theta_t=atan2(yc(n), xc(n))
          call fixangle(theta_t)
          if(N1.ne.N) then
             print *,'level does not match'
             stop
          endif
C ********
C The last parameter of subroutine READDAA is either 1 or 2
C   1: Velocity field of VDAD
C   2: Velocity field of GBVTD or the reflectivity field
C ********
C ******  READ Reflectivity and Doppler Velocity Field ******
          CALL READDATA(DZ,V_D,10,N,XC(N),YC(N),id1,iori,jori,2,
     X         field1,field2)
          if(N.eq.1) then
            call header(id,id1,xc(N),yc(N),sf,af,thetamax,flag1,
     X                  ring_end,kz)
c Change header numbers to print fieldnames correctly
            do i=176,299
               if(id(i).ge.256) then
                  call intswap(id(i),bid(i))
               else
                  bid(i)=id(i)
               endif
            enddo
            if(flag1.gt.2) then
               call header2(id,id2,sf,af,thetamax)
               write(22,116)id2
            endif
            WRITE(15,116)ID
116         format(10i8)
          endif
          write(15,17)N
          if(flag1.gt.2) write(22,17)N
17        format('level',i2)
C reset all data points to -999.0
          do k=1,rmax
            meansd(k,N)=dummy
            meanvt(k,N)=dummy
            meanvr(k,N)=dummy
            do j=1,thetamax
               dddz(k,j)=dummy
               ddvd(k,j)=dummy
               vt(k,j)=dummy
               vr(k,j)=dummy
            enddo
          enddo
C radius of GBVTD rings from ring_start to ring_end km every 1 km
          do K=1,ring_end   
C ******  GRID POINT    ******
C ******  Calculate grid point location for VDAD or GBVTD
C ******  the grid point is in corresponding si coordinate for GBVTD
C ******  theta is angle array (si for GBVTD and theta for VDAD)
            CALL GRID(XC(N),YC(N),XDD,YDD,FLOAT(K),THETAV,theta_t,flag1,
     X                    thetamax)
            call grid(xc(N),yc(N),xdddz,ydddz,float(k),thetaZ,theta_t,1,
     X                    thetamax)
C ******  GET GRID POINT DATA FROM Doppler Velocity Field ******
            CALL INTE(DZ,DUMMY,XDDdz,YDDdz,thetamax,DDZ,2.0,0.5,
     X                iori,jori,id1(162),id1(167))
ccc            CALL INTE(V_D,DUMMY,XDD,YDD,thetamax,DVD,2.0,0.5,iori,
ccc     X                jori,id1(162),id1(167))

            if(flag1.eq.1 .or. flag1.eq.3) then !vdad
              print *, 'int'
              do i=1,thetamax
                print *, dvd(i),xdd(i),ydd(i),iori,jori,d,r_t
                d=sqrt((xdd(i))**2.+(ydd(i))**2.)
                dvd(i) = dvd(i)*d/r_t
              enddo
            endif
            if(flag2.eq.2 .and. (flag1.eq.3 .or. flag1.eq.4)) then  !use analytic wind
              if(k.eq.1) then
                write(17,988)xc(N),yc(N)
                write(17,989)N
              endif
              call analytic(V_D,dummy,xdd,ydd,thetamax,dvd1,
     X                iori,jori,xc(N),yc(N),k,ddz1,flag1)
c                do i=1,thetamax
c                  print *, k,i,xdd(i),xdddz(i),ydd(i),ydddz(i),
c     X                   dvd(i),dvd1(i)
c                enddo
              do i=1,thetamax
                DVD(i)=dvd1(i)
                DDZ(i)=DDZ1(i)
                thetav(i)=360.*(float(i-1))/float(thetamax)*dtor
                print *, k, i, dvd(i), ddz(i), thetav(i)
              enddo
              num_ray=thetamax
            elseif(k.ge.ring_start) then
              num_ray=0
c     i,j are CAPPI indices
              do i=1,num_x
                do j=1,num_y
c     xdist,ydist are distance relative to tc center
                  xdist=i-(xc(n)+iori)
                  ydist=j-(yc(n)+jori)
                  dist=sqrt(xdist**2+ydist**2)
                  if(dist.ge.(k-zone).and.dist.le.(k+zone).and.
     X              v_d(i,j).gt.-90.)then
                    angle=atan2(ydist, xdist)
                    call fixangle(angle)
c                   angle=amod((2.5*pi-angle),2.*pi)-theta_t
                    angle=angle-theta_t
                    call fixangle(angle)
                    if(flag1.eq.2 .or. flag1.eq.4) then  !gbvtd
                      call theta2si(k,angle,xc(n),yc(n),theta,theta_t)
c                     print *, i,j,angle/dtor,theta_t/dtor,theta/dtor,
c     X                 v_d(i,j)
                    else
                      theta=angle
                    endif
                    num_ray=num_ray+1
c     array of angles/velocities in irregular psi, before fourier analysis
                    thetav(num_ray)=theta
                    DVD(num_ray)=V_D(i,j)
c                    print *, i,j,k,num_ray,angle/dtor,
c     X                 thetav(num_ray)/dtor,dvd(num_ray), dist
                  endif
                enddo
              enddo
           endif
          if(k.ge.ring_start) then
            print *, num_ray
            if (num_ray.gt.maxrays) then
              print *, num_ray, maxrays
              stop
            endif
            call maxvel(dvd,num_ray,max_pos,max_neg)
            write(22,265) n,k,max_pos,max_neg
 265        format(1x,2i4,2f8.2)
            call goodcirc(thetaV,DVD,num_ray,GAP,FLAG,numcoeff)
            if(flag) then  !if numcoeff >=3, otherwise, go to next ring
C         GET GOOD DATA NUMBER
C         DVD is the original array, TH is the cleaned up angle array
            CALL GOODDATA(DVD,THETAV,AZ,TH,NN,num_ray,dummy)
C ******  LEASTSQURRE CURVE FIT AND FIND COE. OF FOURIER sqrt.o  ******
            CALL CALCULATE_XY(NN,numcoeff,TH,AZ,XS,YS,NDATA,WEIGHT)
            CALL LLS(numcoeff,NDATA,NDATA,XS,YS,WEIGHT,S,STAND_ERR,
     X                 CC,FLAG)
c           if(flag1.eq.1 .or. flag1.eq.3) then !VDAD
              CALL FINAL1(CC,VTCOEF,VRCOEF,VM,XC(N),YC(N),K,numcoeff,
     X                      flag1,rflag)
            call wrapup(dddz,ddvd,vt,ddz,dummy,vtcoef,flag1,
     *             NN,vrcoef,vm,dtor,k,N,vr,theta_t,vt1d,dvd,
     *             xc(n),yc(n))
C calculate standard deviation
C X is the theta array converted from si (has NN point in deg)
            call std_dvi_new(A,CC,X,AZ,TH,numcoeff,sd,xc(N),yc(N),
     X                     flag1,k,theta_t,dummy,NN)
            if(abs(vtcoef(1)).le.0.001) then
              meansd(k,N)=dummy
              meanvr(k,N)=dummy
              meanvt(k,N)=dummy
            else
              meansd(k,N)=sd
              meanvr(k,N)=vrcoef
              meanvt(k,N)=vtcoef(1)
            endif
C X1 is the theta array converted from si (has thetamax point in deg)
            call std_dvi(A1,CC,X1,AZ,thetaz,numcoeff,sd,xc(N),yc(N),
     X                     flag1,k,theta_t,dummy)
c            print *, k,l,n,numcoeff,theta_t*180./pi, VTCOEF
C write out points file for plotting
            if(k.eq.20 .and. N.eq.2) then
               write(18,180)n,k,NN,CC
               do j=1,thetamax
                  write(18,181)j,thetaz(j)/pi*180.,X1(j)*4.,a1(j)
               enddo
               do j=1,NN
                  write(18,181)j,th(j)/pi*180.,x(j),az(j)
               enddo
            endif
            if(k.eq.35 .and. N.eq.2) then
               write(19,180)n,k,NN,CC
               do j=1,thetamax
                  write(19,181)j,thetaz(j)/pi*180.,X1(j)*4.,a1(j)
               enddo
               do j=1,NN
                  write(19,181)j,th(j)/pi*180.,x(j),az(j)
               enddo
            endif
            if(k.eq.50 .and. N.eq.2) then
               write(20,180)n,k,NN,CC
               do j=1,thetamax
                  write(20,181)j,thetaz(j)/pi*180.,X1(j)*4.,a1(j)
               enddo
               do j=1,NN
                  write(20,181)j,th(j)/pi*180.,x(j),az(j)
               enddo
            endif
181         format(1x,i4,3f8.2)
180         format(1x,i3,i3,i5,/,13f9.2)
c ****  Add mean cross tropical cyclone flow into VT
            call absolute(VT,vr,X,VM,xc(N),yc(N),va,dummy,theta_t,
     X                    k,flag1)
C ****	PLOT  FOURIER COMPONENT  ****
C             IF (K.eq.30.AND.N.EQ.4)THEN
	    CALL PLOTDATA(K,N,IZ,DVD,CC,XC(N),YC(N),flag1,theta_t,
     X                          dummy,a1,X,x1,sd,numcoeff,thetaZ,NN)
            call frame
C	        ENDIF
           else  !  ifnumcoeff < 3
             do j=1,thetamax
               VT(k,j)=dummy
               VA(k,j)=dummy
               VR(k,j)=dummy
               DDVD(k,j)=dummy
             enddo
             do j=1,8 
                vtcoef(j)=dummy
             enddo
             cc(1)=dummy
             vrcoef=dummy
             vm=dummy
             meansd(k,n)=dummy
             j1=-999
             call wrapupdz(ddz, dddz, dummy,theta_t,k)
           endif  !if numcoeff>=3
            write(14,980)k,vrcoef,VM,cc(1),(vtcoef(j),j=1,8),
     X                   meansd(k,n),j1
980         format(i4, 12f6.1,i4)
            do j=1, maxcoeffs
              cc(j)=dummy
              vtcoef(j)=dummy
            enddo
            vrcoef=dummy
           else
             do j=1,thetamax
               VT(k,j)=dummy
               VA(k,j)=dummy
               VR(k,j)=dummy
               DDVD(k,j)=dummy
             enddo
             call wrapupdz(ddz, dddz, dummy,theta_t,k)
           endif
            temp=0.
            index=0
            do j=1,thetamax
              if(dddz(k,j).gt.-100.) then
                temp=temp+dddz(k,j)
                index=index+1
              endif
            enddo
            if (index.gt.0) then
               meandz(k,N)=temp/float(index)
            else
               meandz(k,N)=dummy
            endif
          ENDDO
          do j=1,thetamax
            write(15,18)j
18          format('azimuth',i3)
            write(15,19)bid(176),bid(177),bid(178),bid(179),
     X           (dddz(k,j),k=1,ring_end)
19          format(4a2,/,(8E10.3))
            write(15,19)bid(181),bid(182),bid(183),bid(184),
     X             (VT(k,j),k=1,ring_end)
            write(15,19)bid(186),bid(187),bid(188),bid(189),
     X             (VR(k,j),k=1,ring_end)
            write(15,19)bid(191),bid(192),bid(193),bid(194),
     X             (DDVD(k,j),k=1,ring_end)
            write(15,19)bid(196),bid(197),bid(198),bid(199),
     X             (VA(k,j),k=1,ring_end)
            if(flag1.gt.2) then
              write(22,18) j
              do k=1,ring_end
                 if(abs(VA(k,j)).gt.0.001 .and. 
     X                  abs(dddz(k,j)).le.100. .and. 
     X                  abs(dddz(k,j)).gt.0.001 .and. 
     X                  abs(VA(k,j)).le.100.) then
                    VDIF(k,j)=VA(k,j)-dddz(k,j)
Cc modified 1/12/98 to compute differences in percentage ----------
                    VDIFP(k,j)=(VA(k,j)-dddz(k,j))/dddz(k,j)*100.
CC----------------
                 else
                    VDIF(k,j)=-999.0
                    VDIFP(k,j)=-999.0
                 endif
              enddo
              write(15,19)bid(201),bid(202),bid(203),bid(204),
     X             (VDIF(k,j),k=1,ring_end)
              write(22,19)id2(176),id2(177),id2(178),id2(179),
     X             (VDIFP(k,j),k=1,ring_end)
            endif
          enddo
        enddo
        close(10)
        call header1(id,af,sf,kz,ring_end)
c Change header numbers to print fieldnames correctly
        do i=176,299
           if(id(i).ge.256) then
              call intswap(id(i),bid(i))
           else
              bid(i)=id(i)
           endif
        enddo        
        WRITE(16,116)ID
        N2=1
        write(16,17)N2
        do j=1,zmax
          write(16,18)j
          write(16,19)bid(176),bid(177),bid(178),bid(179),
     X         (meandz(k,j),k=1,ring_end)
          write(16,19)bid(181),bid(182),bid(183),bid(184),
     X         (meanVT(k,j),k=1,ring_end)
          write(16,19)bid(186),bid(187),bid(188),bid(189),
     X         (meanVR(k,j),k=1,ring_end)
          write(16,19)bid(191),bid(192),bid(193),bid(194),
     X         (meansd(k,j),k=1,ring_end)
        enddo
        call std(meansd,rmax,kz,ring_end,kz,dummy)
        CLOSE(15)
        close(16)
        close(14)     
      ENDDO
999   CONTINUE
      CALL FRAME
      CALL CLSGKS
      STOP
998       print *, 'cannot open ideal input file'
      END

************************************************************************
*  std
************************************************************************

      subroutine std(meansd,rmax,zmax,ring_end, kz,dummy)

      implicit none
      integer rmax,zmax,i,j,levelnum,volumenum,ring_end,kz
      real meansd(rmax,zmax),levelmeansd(20),volumemeansd,dummy
      real levelsd,volumesd

      volumenum=0
      volumemeansd=0.
      volumesd=0.
      do j=1,kz
         levelnum=0
         levelmeansd(j)=0.
         levelsd=0.
         do i=1,ring_end
            if(abs(meansd(i,j)-dummy).gt.0.0001) then
               levelmeansd(j)=levelmeansd(j)+meansd(i,j)
               levelnum=levelnum+1
            endif
         enddo
         if(levelnum.gt.0) then
            volumenum=volumenum+1
            levelmeansd(j)=levelmeansd(j)/float(levelnum)
            volumemeansd=volumemeansd+levelmeansd(j)
         endif
         do i=1,rmax
            if(abs(meansd(i,j)-dummy).gt.0.0001) then
               levelsd=levelsd+(meansd(i,j)-levelmeansd(j))**2
            endif
         enddo
         levelsd=sqrt(levelsd/float(levelnum))            
         write(14, 15) j, levelnum, levelmeansd(j), levelsd
      enddo
      volumemeansd=volumemeansd/float(volumenum)
      do j=1,zmax
         if(abs(levelmeansd(j)-dummy).gt.0.0001) then
            volumesd=volumesd+(levelmeansd(j)-volumemeansd)**2
         endif
      enddo
      volumesd=sqrt(volumesd/float(volumenum))
      write(14,16) volumenum, volumemeansd, volumesd
15    format(i3, i5, 2f9.2)
16    format(i6, 2f9.2)
      return
      end 
   

************************************************************************
* absolute VA
************************************************************************

      Subroutine absolute(VT,vr,X,VM,xc,yc,va,dummy,theta_t,i,flag1)
      implicit none
      INCLUDE 'gbvtd_parameter.f'
      real vt(rmax,thetamax),va(rmax,thetamax),X(thetamax)
      real xc,yc,theta_t,pi,factor,dummy,vm,th,temp,u,v
      real vr(rmax,thetamax)
      integer i,j,flag1

      pi=acos(-1.)
c      VM=10.
      do j=1,thetamax
        if((vt(i,j)-dummy).gt.0.001)then
         temp=float(j-1)*360./float(thetamax)*pi/180.
c         factor=-vm*cos(temp+theta_t)
         u=-vt(i,j)*sin(temp)+vm*cos(theta_t)
     X                 +vr(i,j)*cos(temp)
         v=vt(i,j)*cos(temp)+vm*sin(theta_t)
     X                 +vr(i,j)*sin(temp)
         va(i,j)=sqrt(u*u+v*v)
        endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c         factor=vm*cos(temp+theta_t)
c         if(abs(vt(i,j)-dummy).gt.0.001) then
c            va(i,j)=vt(i,j)+factor
c         else
c            va(i,j)=dummy
c         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(vt(i,j).gt.60.) then
         print *,va(i,j),vt(i,j),u,v,vm*cos(theta_t),
     X          vm*sin(theta_t),(temp)/pi*180.,theta_t/pi*180.
         endif
      enddo
      return
      end

      
************************************************************************
* CAL_VT
************************************************************************

      Subroutine wrapup(dddz,ddvd,vt,ddz,dummy,vtcoef,flag1,
     *    NN,vrcoef,vm,dtor,k,N,vr,theta_t,vt1d,dvd,xc,yc)

      implicit none
      INCLUDE 'gbvtd_parameter.f'
      integer j,k,dj,i,flag1,NN,j1,N
      real dummy, temp, si, theta_t,xc,yc,theta1,vrcoef
      real dddz(rmax,thetamax), ddvd(rmax,thetamax),vt(rmax,thetamax)
      real tempdz(thetamax),dvd(thetamax),vtcoef(maxcoeffs)
      real ddz(thetamax), vr(rmax,thetamax),theta(thetamax)
      real dtor,vm,vt1d(thetamax)

c Comment out following for tornado data MMB
c      IF (NN.LE.60)THEN
c         IF (NN.LE.40)THEN
c            DO I=1,6
c               VTCOEF(I)=dummy
c	    ENDDO
c         ENDIF
c         VRCOEF=dummy
c         VM=dummy
c      ENDIF

      do i=1,thetamax
         tempdz(i)=ddz(i)
          vr(k,i)=vrcoef
      enddo
      dj=int(theta_t/dtor/360.*float(thetamax))
      DO I=1,6
         if(vtcoef(i).le.-90.) vtcoef(i)=0.
      ENDDO
      if(vrcoef .ge. 90.) vrcoef=dummy
      if(vm .le. -90.) vm=0.
      j1=0          
c     j is in theta coordinates
      do j=1,thetamax
         temp=dummy
         theta1=(j-1)*360./float(thetamax)*dtor
C                         !theta1 is azimuth angle for VDAD. For GBVTD,
C                         !theta1 needs to be converted to si for GBVTD 
         if(flag1.eq.2 .or. flag1.eq.4) then !GBVTD
            call theta2si(k,theta1,xc,yc,si,theta_t)
            theta1=si
         endif
         if(abs(vtcoef(1)-dummy).gt. 0.001) then
            temp=vtcoef(1)
            if(abs(vtcoef(3)-dummy).gt.0.001 .and.
     X           abs(vtcoef(4)-dummy) .gt. 0.001) then
               temp=temp+vtcoef(3)*cos(theta1)
     X             +vtcoef(4)*sin(theta1)
               if(abs(vtcoef(5)-dummy).gt.0.001 .and.
     X                abs(vtcoef(6)-dummy) .gt. 0.001) then
                  temp=temp+vtcoef(5)*cos(theta1*2.)
     X                    +vtcoef(6)*sin(theta1*2.)
                  if(abs(vtcoef(7)-dummy).gt.0.001 .and.
     X                abs(vtcoef(8)-dummy) .gt. 0.001) then
                     temp=temp+vtcoef(7)*cos(theta1*3.)
     X                    +vtcoef(8)*sin(theta1*3.)
                  endif
               endif
            endif
         endif
c         print *,n,k,j,temp,vm,vrcoef
         if(j+dj .gt. thetamax)then
            dddz(k,j+dj-thetamax)=tempdz(j)
            ddvd(k,j+dj-thetamax)=DVD(j)
            vt(k,j+dj-thetamax)=temp
            vt1d(j+dj-thetamax)=temp
         elseif(j+dj .lt. 1) then
            dddz(k,j+dj+thetamax)=tempdz(j)
            ddvd(k,j+dj+thetamax)=DVD(j)
            vt(k,j+dj+thetamax)=temp
            vt1d(j+dj+thetamax)=temp
         else
            dddz(k,j+dj)=tempdz(j)
            ddvd(k,j+dj)=dvd(j)
            vt(k,j+dj)=temp
            vt1d(j+dj)=temp
         endif
      enddo
      return
      end

*********************************************************************
      Subroutine wrapupdz(ddz,dddz,dummy,theta_t,k)

      implicit none
      INCLUDE 'gbvtd_parameter.f'
      integer j,k,dj,i,flag1,NN,j1,N
      real dummy, temp, si, theta_t,xc,yc,theta1,vrcoef
      real dddz(rmax,thetamax)
      real tempdz(thetamax)
      real ddz(thetamax),dtor

      dtor=acos(-1.)/180.
      do i=1,thetamax
         tempdz(i)=ddz(i)
      enddo
      dj=int(theta_t/dtor/360.*float(thetamax))
      j1=0          
      do j=1,thetamax
         temp=dummy
         if(j+dj .gt. thetamax)then
            dddz(k,j+dj-thetamax)=tempdz(j)
         elseif(j+dj .lt. 1) then
            dddz(k,j+dj+thetamax)=tempdz(j)
         else
            dddz(k,j+dj)=tempdz(j)
         endif
      enddo
      return
      end

************************************************************************
* OPENFILE
************************************************************************

      Subroutine openfile(l,numloc,kc,prefix)

      integer l, numloc,kc,iz,lnblnk
      CHARACTER FNAME*80,STR*2,XFNAME*30
      character prefix*20
      iz =kc-1
            if(iz.lt.10) then
              str='0'
              write(str,'(i1)')IZ
            else
              write(str,'(i2)')IZ
            endif
            xfname=prefix(:lnblnk(prefix))//'f'//STR
            if(L.lt.10) then
              str='0'
              write(str,'(i1)')L
            else
              write(str,'(i2)')L
            endif
            xfname=xfname(:lnblnk(xfname))//'c'
            xfname=xfname(:lnblnk(xfname))//STR
            OPEN(15,FILE=xfname,STATUS='UNKNOWN',FORM='FORMATTED')
            if(iz.lt.10) then
              str='0'
              write(str,'(i1)')IZ
            else
              write(str,'(i2)')IZ
            endif
            xfname=prefix(:lnblnk(prefix))//'f'//STR
            if(L.lt.10) then
              str='0'
              write(str,'(i1)')L
            else
              write(str,'(i2)')L
            endif
            xfname=xfname(:lnblnk(xfname))//'z'
            xfname=xfname(:lnblnk(xfname))//STR
            OPEN(16,FILE=xfname,STATUS='UNKNOWN',FORM='FORMATTED')
            if(iz.lt.10) then
              str='0'
              write(str,'(i1)')IZ
            else
              write(str,'(i2)')IZ
            endif
            xfname=prefix(:lnblnk(prefix))//'f'//STR
            if(L.lt.10) then
              str='0'
              write(str,'(i1)')L
            else
              write(str,'(i2)')L
            endif
            xfname=prefix(:lnblnk(prefix))//STR
            xfname=xfname(:lnblnk(xfname))//'.out'
            OPEN(14,FILE=xfname,STATUS='unknown',FORM='formatted')
            xfname=xfname(:lnblnk(xfname)-3)//'ana'
            open(17,file=xfname,status='unknown',form='formatted')
            xfname=xfname(:lnblnk(xfname)-3)//'dif'
            open(22,file=xfname,status='unknown',form='formatted')
            xfname=xfname(:lnblnk(xfname)-3)//'point1'
            open(18,file=xfname,status='unknown',form='formatted')
            xfname=xfname(:lnblnk(xfname)-1)//'2'
            open(19,file=xfname,status='unknown',form='formatted')
            xfname=xfname(:lnblnk(xfname)-1)//'3'
            open(20,file=xfname,status='unknown',form='formatted')
       return
      end

********************************************************************
      subroutine header(id,id1,xc,yc,sf,af,thetamax,flag1,xmax,kz)

      implicit none
      
      real xc,yc,temp,rad_lat,rad_long,ori_lat,ori_long
      integer thetamax,flag1,xmax,kz
      integer*2 id(510),sf,af,i,id1(510)
c            call latlong(rad_lat,rad_long,xc,yc,ori_lat,ori_long)
c            id(33)=int(ori_lat)
c            temp=(ori_lat-float(id(33)))*60.
c            id(34)=int(temp)
c            id(35)=(int(temp*60.-float(id(34)))*60.)*sf
c            id(36)=int(ori_long)
c            temp=(ori_long-float(id(36)))*60.
c            id(37)=int(temp)
c            id(38)=(int(temp*60.-float(id(37)))*60.)*sf
      do i=1,510
        id(i)=id1(i)
      enddo
      id(309)=-xc*sf
      id(310)=-yc*sf
      call str2int('PO',id(16))
      call str2int('LA',id(17))
      id(68)=100
      id(69)=64
      id(160)=1*sf
      id(161)=xmax*sf
      id(162)=xmax
      id(163)=1000
      id(165)=0*af
      id(166)=356*af
      id(167)=thetamax
      id(168)=4*af
C------------------------------------------------------------------
C commented by wcl 8/9/97, The levels in polar coordinate should
C    agree with the levels in CAPPI, therefore, no need to set
C    new value here.
cc mbell 04/04/03 Changed to allow grid2ps to work if less levels
cc are processed than in the cappi file
C------------------------------------------------------------------
c      id(170)=1000
      if (id(172).gt.kz) then
c Reduce it
         id(171)=kz*1000
         id(172)=kz
      endif
      id(173)=1000
      id(175)=5
      call str2int('DZ',id(176))
      call str2int('  ',id(177))
      call str2int('  ',id(178))
      call str2int('  ',id(179))
      id(180)=1
      call str2int('VT',id(181))
      call str2int('  ',id(182))
      call str2int('  ',id(183))
      call str2int('  ',id(184))
      id(185)=1
      call str2int('VR',id(186))
      call str2int('  ',id(187))
      call str2int('  ',id(188))
      call str2int('  ',id(189))
      id(190)=1
      call str2int('VD',id(191))
      call str2int('  ',id(192))
      call str2int('  ',id(193))
      call str2int('  ',id(194))
      id(195)=1
      call str2int('VA',id(196))
      call str2int('  ',id(197))
      call str2int('  ',id(198))
      call str2int('  ',id(199))
      id(200)=1
      if(flag1.gt.2) then
        call str2int('VD',id(201))
        call str2int('IF',id(202))
        call str2int('  ',id(203))
        call str2int('  ',id(204))
        id(205)=1
        id(175)=6
      endif
      id(303)=1 !number of radars
c      call str2int('CA',id(306))
c      call str2int('A ',id(307))

      return
      end
***********************************************************************      
      subroutine header1(id,af,sf,kz,xmax)

      implicit none
      
      integer*2 af,sf
      integer zmax,kz,xmax
      integer*2 id(510),i

      call str2int('CA',id(16))
      call str2int('RT',id(17))
      id(160)=1*sf
      id(161)=xmax*sf
      id(162)=xmax
      id(163)=1000
      id(165)=1*sf
      id(166)=kz*sf
      id(167)=kz
      id(168)=1000
      id(170)=1000
      id(171)=1000
      id(172)=1
      id(173)=1000
      id(175)=4
      call str2int('MD',id(176))
      call str2int('Z ',id(177))
      call str2int('  ',id(178))
      call str2int('  ',id(179))
      id(180)=1
      call str2int('MV',id(181))
      call str2int('T ',id(182))
      call str2int('  ',id(183))
      call str2int('  ',id(184))
      id(185)=1
      call str2int('MV',id(186))
      call str2int('R ',id(187))
      call str2int('  ',id(188))
      call str2int('  ',id(189))
      id(190)=1
      call str2int('MS',id(191))
      call str2int('D ',id(192))
      call str2int('  ',id(193))
      call str2int('  ',id(194))
      id(195)=1
      id(303)=1 !number of radars
c      call str2int('CA',id(306))
c      call str2int('A ',id(307))

      return
      end
      
********************************************************************
      subroutine header2(id,id2,sf,af,thetamax)

      implicit none
      
      real xc,yc,temp
      integer thetamax,flag1
      integer*2 id(510),sf,af,i,id2(510)
      do i=1,510
        id2(i)=id(i)
      enddo
      call str2int('CA',id2(16))
      call str2int('RT',id2(17))
      id2(175)=1
      call str2int('VD',id2(176))
      call str2int('IF',id2(177))
      call str2int('P ',id2(178))
      call str2int('  ',id2(179))
      id2(165)=1*sf
      if(thetamax.le.320) then 
        id2(166)=thetamax*sf
        id2(167)=thetamax
      else
        print *, 'thetamax too large'
        stop
      endif
      id2(168)=1000
      id2(303)=1 !number of radars
      return
      end
          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Calculate Typhoon center lat and long based on the radar lat and 
C   long and the relative coordinate XC and YC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine latlong(rad_lat,rad_long,xc,yc,ori_lat,ori_long)

      implicit none

      real rad, rad_lat, rad_long, ori_lat, ori_long, xc, yc
      real pi, rad_laty, rad_longx,radr 
      RAD=6371.25
      PI=ACOS(-1.)
C Lat and Long of the CAA radar
      RAD_LAT=25.0694
      RAD_LONG=121.2083
      rad_laty=yc/rad*180./pi    !calculate delta lat from the center y coord.
      radr=rad*cos((rad_lat+rad_laty)*pi/180.) !calculate radius at the typhoon
C                                      !latitude
      rad_longx=xc/radr*180./pi   !calculate delta long from the center x coord.
      ori_lat=rad_lat+rad_laty           !typhoon lat.
      ori_long=rad_long+rad_longx           !typhoon long.

      return
      end


*************************************************************
*   DATA
*************************************************************

      subroutine data(dddz,thetamax)
      implicit none
      integer kk,nn1,nn2,index,i,k,thetamax
      real dddz(thetamax)
      nn1=0
      nn2=0
      do kk=1,thetamax
      if (int(dddz(kk)).ne.99) then
      nn1=kk
      goto 12
      endif
      enddo
12    do kk=thetamax,1,-1
      if (int(dddz(kk)).ne.99) then
      INDEX=1
      nn2=kk
      goto 13
      endif
      enddo
13    if (nn1.eq.0)then
      INDEX=0
      dddz(i)=0.
      goto 30
      endif
      do i=nn1+1,nn2-1
      if (int(dddz(i)).eq.99) then
      do k=i+1,nn2
      if (int(dddz(k)).ne.99) then
      dddz(i)=dddz(i-1)+(dddz(k)-dddz(i-1))/(k-i+1)
      goto 10
      endif
      enddo
      endif
10    enddo
30    return
      end


***********************************************************
*   PLOTDATA
***********************************************************

	SUBROUTINE PLOTDATA(K,N,IZ,DVD,C,XC,YC,flag,theta_t,
     X                      dummy,A1,X,X1,sd,numcoeff,theta,NN)

        implicit none
        include 'gbvtd_parameter.f'
        real pi,th,sd,w0,w1,w2,w3,w4,si,xc,yc,theta_t
        integer i,in,iz,k,n,flag,numcoeff,NN
	CHARACTER STR0*16,STR1*14,STR2*14,STR3*3,S(5)*4
	CHARACTER*16  SW0,SW1,SW2,SW3,SW4
	REAL DVD(NN),dummy,theta(thetamax)
        real A1(thetamax),X1(thetamax),X(NN),C(maxcoeffs)
C	DATA S /'1041','1101','1121','1141','2:IGL:P'/
	DATA S /'0917','0947','0602','0632','2:IGL:P'/
	CALL SETUSV('LW',2500)
	IF (IZ.eq.1)CALL SET(0.05,0.45,0.75,0.95,1.,360.,-80.,80.,1)
	IF (IZ.eq.2)CALL SET(0.5,0.9,0.75,0.95,1.,90.,-80.,80.,1)
	IF (IZ.eq.3)CALL SET(0.05,0.45,0.45,0.65,1.,90.,-80.,80.,1)
	IF (IZ.eq.4)CALL SET(0.5,0.9,0.45,0.65,1.,90.,-80.,80.,1)
	PI=ACOS(-1.)
        do i=1,NN
          CALL PLCHHQ(X(I),DVD(I),':IRL:5',0.020,0.,0.)
        enddo
        do i=1,thetamax
          if(flag.eq.2 .or. flag.eq.4) then     !gbvtd
            theta(i)=x1(i)*360./float(thetamax)
          else
            theta(i)=theta(i)*180./pi/360.*float(thetamax)
          endif
c          theta(i)=theta(i)*180./pi/360.*float(thetamax)
        enddo
        W0=C(1)
	W1=sqrt(C(2)**2.+C(3)**2.)
	W2=sqrt(C(4)**2.+C(5)**2.)
	W3=sqrt(C(6)**2.+C(7)**2.)
        W4=sqrt(C(8)**2.+C(9)**2.)
	WRITE(SW0,'(A7,F5.1,A4)')'wave 0=',W0,' m/s '
	WRITE(SW1,'(A7,F5.1,A4)')'wave 1=',W1,' m/s '
	WRITE(SW2,'(A7,F5.1,A4)')'wave 2=',W2,' m/s '
	WRITE(SW3,'(A7,F5.1,A4)')'wave 3=',W3,' m/s '
	WRITE(SW4,'(A7,F5.1,A4)')'wave 3=',W4,' m/s '
	IF (SD.LT.1.)THEN
	WRITE(STR0,'(A7,A3,F2.1,A4)')' MSE  =','  0',SD,' m/s '
	ELSE
	WRITE(STR0,'(A7,F5.1,A4)')' MSE  =',SD,' m/s '
	ENDIF
	WRITE(STR2,'(A6,I4,A4)')'RAD.=',K,' Km  '
	WRITE(STR1,'(A6,I4,A4)')'HEI.=',N,' Km  '
        CALL CURVE(theta,A1,thetamax)
	CALL PERIM(0,4,0,8)
	CALL SETUSV('LW',1500)
c        CALL LINE(1.,0.,90.,0.)
        CALL PLCHHQ(25.,35.,STR0,0.008,0.,0.)
        CALL PLCHHQ(25.,28.,SW0,0.008,0.,0.)
        CALL PLCHHQ(25.,21.,SW1,0.008,0.,0.)
        CALL PLCHHQ(25.,14.,SW2,0.008,0.,0.)
        CALL PLCHHQ(25.,7., SW3,0.008,0.,0.)
        CALL PLCHHQ(25.,1., SW4,0.008,0.,0.)
        CALL PLCHHQ(75.,-25.,S(IZ),0.010,0.,0.)
	CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
c	call pcsetc('FC','/')
c	CALL PLCHHQ(0.22,0.7,'TYPHOON GLADYS 1994.09.01',0.013,0.,0.)
	CALL PLCHHQ(0.22,0.7,'TYPHOON ALEX 1987.07.27',0.013,0.,0.)
	CALL PLCHHQ(0.82,0.7,'CURVE FIT',0.010,0.,0.)
        CALL PLCHHQ(0.50,0.7,STR1,0.010,0.,0.)
        CALL PLCHHQ(0.68,0.7,STR2,0.010,0.,0.)
	DO I=1,9,2
	WRITE(STR3,'(I3)')(I-5)*15
	IF (IZ.EQ.1)CALL PLCHHQ(0.03,0.75+(I-1)*0.2/8.,STR3,0.006,0.,0.)
	IF (IZ.EQ.2)CALL PLCHHQ(0.03,0.45+(I-1)*0.2/8.,STR3,0.006,0.,0.)
	IF (IZ.EQ.3)CALL PLCHHQ(0.48,0.75+(I-1)*0.2/8.,STR3,0.006,0.,0.)
	IF (IZ.EQ.4)CALL PLCHHQ(0.48,0.45+(I-1)*0.2/8.,STR3,0.006,0.,0.)
	ENDDO
	IF (IZ.EQ.1)CALL PLCHHQ(0.01,0.85,'m/s',0.006,0.,0.)
	IF (IZ.EQ.2)CALL PLCHHQ(0.01,0.55,'m/s',0.006,0.,0.)
	IF (IZ.EQ.1)CALL PLCHHQ(0.05,0.74,' 0 ',0.006,0.,0.)
	IF (IZ.EQ.2)CALL PLCHHQ(0.05,0.44,' 0 ',0.006,0.,0.)
	IF (IZ.EQ.3)CALL PLCHHQ(0.5,0.74,' 0 ',0.006,0.,0.)
	IF (IZ.EQ.4)CALL PLCHHQ(0.5,0.44,' 0 ',0.006,0.,0.)
	IF (IZ.EQ.1)CALL PLCHHQ(0.25,0.74,':PGL:P',0.006,0.,0.)
	IF (IZ.EQ.2)CALL PLCHHQ(0.25,0.44,':PGL:P',0.006,0.,0.)
	IF (IZ.EQ.3)CALL PLCHHQ(0.7,0.74,':PGL:P',0.006,0.,0.)
	IF (IZ.EQ.4)CALL PLCHHQ(0.7,0.44,':PGL:P',0.006,0.,0.)
	IF (IZ.EQ.1)CALL PLCHHQ(0.45,0.74,'2:PGL:P',0.006,0.,0.)
	IF (IZ.EQ.2)CALL PLCHHQ(0.45,0.44,'2:PGL:P',0.006,0.,0.)
	IF (IZ.EQ.3)CALL PLCHHQ(0.9,0.74,'2:PGL:P',0.006,0.,0.)
	IF (IZ.EQ.4)CALL PLCHHQ(0.9,0.44,'2:PGL:P',0.006,0.,0.)
	RETURN
	END

*****************************************************************
*   READDATA
*****************************************************************

       subroutine readdata(dz,ve,unit,ik,xc,yc,id,iori,jori,flag,
     X     field1,field2)

       implicit none
       integer ik,i,j,n,nx,ny,kz,k,iori,jori,ix,iy,j1
       real xc,yc,pi,special,d,R_T
       integer flag,unit,indexdz,indexve,index
       real dz(241,241),ve(241,241),verel(241,241)
       integer*2 id1(510),id(510)
       integer*2 field1,field2
      
       special=-999.
       pi=acos(-1.)
       R_T=sqrt(xc*xc+yc*yc)
        ix=id(162)
        iy=id(167)

       do j=1,id(175)   !id(175):total number of fields in the ascii file
          index=176+(j-1)*5
          if(id(index).eq.field1) indexdz=j
          if(id(index).eq.field2) indexve=j
       enddo 
       
       do j=1,iy
         read(unit,18)j1
c           print *, j1,j
         if(j1.ne.j) then
           print *, j1,j,' second index not match'
           stop
         endif
18       format(7x,i3)
         do k = 1, id(175)
           if(k.eq.indexdz) then
             read(unit,19)id(176),id(177),id(178),id(179),
     X             (dz(i,j),i=1,ix)
           elseif(k.eq.indexve) then
             read(unit,19)id(181),id(182),id(183),id(184),
     X             (ve(i,j),i=1,ix)
           else
             read(unit,19)id(186),id(187),id(188),id(189),
     X             (verel(i,j),i=1,ix)
           endif
         enddo
       enddo
19       format(4a2,/,(8E10.3))
c       DO j = 1, iy
c         DO i = 1, ix
c           IF (abs(ve(i,j)-special).gt.0.001) then
c             if(flag.eq.1 .or. flag.eq.3) then !vdad
c               d=sqrt((i-iori)**2.+(j-jori)**2.)
c               ve(i,j) = ve(i,j)*d/R_T
c             elseif(flag.ne.2 .and. flag.ne.4) then ! not gbvtd or reflectivity
c               print *, 'flag out of range'
c             endif
c           endif
c         END DO
c       END DO
       return
       end

***************************************************************
*   INTE
***************************************************************

       SUBROUTINE INTE(dz,dummy,xd,yd,nx,gdz,rang,gama,iori,jori,
     X                 xmax,ymax)

       implicit none
       real rang, gama,dummy,pi,wt,det0,diff,totdiff
       integer nx, idummy, ik1, ik2, ij1, ij2, i
       integer iori,jori
       real dz(241,241),xd(nx),yd(nx)
       real gdz(nx),temp1,temp2,ave1,ave2,ave3
       integer*2 xmax,ymax

C       real tme(nx)
       PI=ACOS(-1.)
       idummy=dummy
       wt=-4.*gama/(rang*rang)
       DET0=PI/210.
       do 100 i=1,nx
         ik1=int(xd(i))+iori
         ik2=int(xd(i))+1+iori
         ij1=int(yd(i))+jori
         ij2=int(yd(i))+1+jori
         if (ik1.lt.1)goto 99
         if (ij1.lt.1)goto 99
         if (ik2.gt.xmax)goto 99
         if (ij2.gt.ymax)goto 99
C-----------------------------------------------------
C Commented by wcl 8/9/97 because closest point interpolation
C   produced a bias on the pattern which is not tolerant by GBVTD analysis
C   So, change it to bi-linear interpolation.
C--------------------------------------------------------------
c         if ((dz(ik1,ij1)-dummy).gt.0.0001) then 
c           gdz(i)=dz(ik1,ij1)
c         else
c           gdz(i)=dummy
c         endif
         temp1=dz(ik1,ij1)
         temp2=dz(ik1,ij2)
         totdiff=float(ij2-ij1)
         diff=yd(i)-int(yd(i))
         call ave(temp1,temp2,diff,totdiff,ave1,dummy)
         temp1=dz(ik2,ij1)
         temp2=dz(ik2,ij2)
         call ave(temp1,temp2,diff,totdiff,ave2,dummy)
         totdiff=float(ik2-ik1)
         diff=xd(i)-int(xd(i))
         call ave(ave1,ave2,diff,totdiff,ave3,dummy)
         gdz(i)=ave3
c         print *,i,ave1,ave2,ave3,gdz(i)
         goto 100
99       gdz(i)=dummy
100    continue
       return
       end

********************************************************************
*ANALYTIC
********************************************************************

      subroutine analytic(dz,dummy,xd,yd,nx,gdz,iori,jori,
     X                    xc,yc,k,ddz,flag)

       implicit none
       real dummy,pi
       integer nx, i
       integer iori,jori,k,flag
       real disturb
       real dz(241,241),xd(nx),yd(nx),rt,xg,yg,theta
       real gdz(nx),radius,ddz(nx),xc,yc,asymvt,ampvt1
       real phasevt1,dtor,ampvt2,phasevt2,ampvt3,phasevt3,a,vorticity
       real ampvr1,phasevr1,ampvr2,phasevr2,ampvr3,phasevr3
       real rad_v,vrmean,u,meanu,v,special,rg,meanflow
       real meanangle,meanv,r_t,vta,angle,theta_t,asymvr
       real alpha,s,power,rotate,theta_v,vtas,rad_vs
       character outfiledz*80,outfileve*80,inputfile*80
       common /ideal/vorticity,a,meanflow,meanangle,phasevt1,ampvt1,
     X               phasevt2,ampvt2,phasevt3,ampvt3,phasevr1,ampvr1,
     X               phasevr2,ampvr2,phasevr3,ampvr3,vrmean,alpha,s,
     X               power,rotate

      
       PI=ACOS(-1.)
       dtor=pi/180.
       rt=sqrt(xc*xc+yc*yc)
       theta_t=atan2(yc,xc)/dtor
       print *, theta_t
       meanu=meanflow*sin(meanangle*dtor-pi)
       meanv=meanflow*cos(meanangle*dtor-pi)
       print *, 'ana'
c Computing velocity 
        do i=1,nx
c xd,yd are x,y of each ring point in the regular psi coordinates
          xg=(xd(i))
          yg=(yd(i))
          theta_v=atan2(yg-yc,xg-xc) !Math angle in vortex coords.
          angle=360./float(nx)*float(i-1)*dtor
C For GBVTD, angle is the si coord., the following if statement computes
C the corresponding math angle so the analytical value can be assigned.
          if(flag.eq.2 .or. flag.eq.4) then
            call si2theta(k,angle,xc,yc,angle,theta_t)
          endif
          theta=amod((angle/dtor+theta_t),360.)*dtor
          radius=float(k)
          asymvt=ampvt1*cos(theta-(phasevt1*dtor))
     X              +ampvt2*cos(2*(theta-(phasevt2*dtor)))
     X              +ampvt3*cos(3*(theta-(phasevt3*dtor)))
c          print *, i,j, theta/dtor, asym,phase1,phase2,phase3
          asymvr=ampvr1*cos(theta-(phasevr1*dtor))
     X              +ampvr2*cos(2*(theta-(phasevr2*dtor)))
     X              +ampvr3*cos(3*(theta-(phasevr3*dtor)))
          if(radius.le.a) then  !solid rotation
c            VTa=0.5*vorticity*radius*(1.+asymvt)-(vorticity/2*alpha*
c     X        ((radius**(s-1))/(a**(s-1)))*cos(s*theta+rotate*dtor))
c             rad_v=0.1*vrmean*sqrt((a-radius)*radius)*(1.+asymvr)-
c     X           (alpha*vorticity/2)*(radius**(s-1))/
c     X           (a**(s-1))*sin(s*theta+rotate*dtor)
            VTa=0.5*vorticity*radius*(1.+asymvt)
            rad_v=0.1*vrmean*sqrt((a-radius)*radius)*(1.+asymvr)
            if(s.gt.0.5) then
               disturb=-(vorticity*alpha/2.)
               vtas=disturb*cos(s*(theta+rotate*dtor))
               rad_vs=disturb*sin(s*(theta+rotate*dtor))
c               vtas=-(vorticity/2.)*(alpha*(radius**(s-1.))/  
c     X          (a**(s-1.))*cos(s*theta+rotate*dtor))
c               rad_vs=-(alpha*vorticity/2.)*(radius**(s-1.))/
c     X           (a**(s-1.))*sin(s*theta+rotate*dtor)
c               print *, "in asy", s, s*theta/dtor, cos(s*theta),vta,
c     x                  vtas, rad_v, rad_vs

               VTa=VTas+VTa
               rad_v=rad_vs+rad_v
C Uncomment these two lines to have just a disturbance
c               vta=vtas
c               rad_v=rad_vs
           endif
          else
c            VTa=0.5*vorticity*a*((a/radius)**power)*(1.+asymvt)+
c     X            (vorticity/2)*(alpha*(a**(s+1))/(radius**(s+1))
c     X            *cos(s*theta+rotate*dtor))
c            rad_v=-vrmean*sqrt((radius-a))*a/radius*(1.+asymvr)-
c     X           (vorticity*alpha/2)*((a**(s+1))/
c     X           (radius**(s+1)))*sin(s*theta+rotate*dtor)
            VTa=0.5*vorticity*a*((a/radius)**power)*(1.+asymvt)
            rad_v=-vrmean*sqrt((radius-a))*a/radius*(1.+asymvr)
            if(s.gt.0.5) then
              VTas=(vorticity/2.)*(alpha*(a**(s+1.))/(radius**(s+1.))
     X            *cos(s*(theta+rotate*dtor)))
              rad_vs=-(vorticity*alpha/2.)*((a**(s+1.))/
     X           (radius**(s+1.)))*sin(s*(theta+rotate*dtor))
               print *, "in asy", s
               VTa=VTas+VTa
               rad_v=rad_vs+rad_v
C Uncomment these two lines to have just a disturbance
c               vta=vtas
c               rad_v=rad_vs
            endif
          endif
c-- writeout summary analytic TC
          if(i.eq.1) then
            write(17,141)k,meanflow,meanangle,vta/(1.+asymvt),
     X               ampvt1,phasevt1, ampvt2,phasevt2,ampvt3,
     X               phasevt3,rad_v/(1.+asymvr)
 141        format(i4,10f6.1)
          endif
          if(radius.ge. 2.) then 
             u=meanu-vta*sin(theta)+rad_v*cos(theta)
             v=meanv+vta*cos(theta)+rad_v*sin(theta)
            rg=sqrt(xg*xg+yg*yg)
            if(rg.ge.0.01) then
              gdz(i)=(u*xg+v*yg)/rg
            else
              gdz(i)=dummy
            endif
          else
            gdz(i)=dummy
          endif
c Reflectivity has to be computed in psi, since it is not recomputed elsewhere
          angle=360./float(nx)*float(i-1)*dtor
          theta=amod((angle/dtor+theta_t),360.)*dtor
          radius=float(k)
          asymvt=ampvt1*cos(theta-(phasevt1*dtor))
     X              +ampvt2*cos(2*(theta-(phasevt2*dtor)))
     X              +ampvt3*cos(3*(theta-(phasevt3*dtor)))
c          print *, i,j, theta/dtor, asym,phase1,phase2,phase3
          asymvr=ampvr1*cos(theta-(phasevr1*dtor))
     X              +ampvr2*cos(2*(theta-(phasevr2*dtor)))
     X              +ampvr3*cos(3*(theta-(phasevr3*dtor)))
          if(radius.le.a) then  !solid rotation
c            VTa=0.5*vorticity*radius*(1.+asymvt)-(vorticity/2*alpha*
c     X        ((radius**(s-1))/(a**(s-1)))*cos(s*theta+rotate*dtor))
c             rad_v=0.1*vrmean*sqrt((a-radius)*radius)*(1.+asymvr)-
c     X           (alpha*vorticity/2)*(radius**(s-1))/
c     X           (a**(s-1))*sin(s*theta+rotate*dtor)
            VTa=0.5*vorticity*radius*(1.+asymvt)
            rad_v=0.1*vrmean*sqrt((a-radius)*radius)*(1.+asymvr)
            if(s.gt.0.5) then
               disturb=-(vorticity*alpha/2.)
               vtas=disturb*cos(s*(theta+rotate*dtor))
               rad_vs=disturb*sin(s*(theta+rotate*dtor))
c               vtas=-(vorticity/2.)*(alpha*(radius**(s-1.))/  
c     X          (a**(s-1.))*cos(s*theta+rotate*dtor))
c               rad_vs=-(alpha*vorticity/2.)*(radius**(s-1.))/
c     X           (a**(s-1.))*sin(s*theta+rotate*dtor)
c               print *, "in asy", s, s*theta/dtor, cos(s*theta),vta,
c     x                  vtas, rad_v, rad_vs

               VTa=VTas+VTa
               rad_v=rad_vs+rad_v
C Uncomment these two lines to have just a disturbance
c               vta=vtas
c               rad_v=rad_vs
           endif
          else
c            VTa=0.5*vorticity*a*((a/radius)**power)*(1.+asymvt)+
c     X            (vorticity/2)*(alpha*(a**(s+1))/(radius**(s+1))
c     X            *cos(s*theta+rotate*dtor))
c            rad_v=-vrmean*sqrt((radius-a))*a/radius*(1.+asymvr)-
c     X           (vorticity*alpha/2)*((a**(s+1))/
c     X           (radius**(s+1)))*sin(s*theta+rotate*dtor)
            VTa=0.5*vorticity*a*((a/radius)**power)*(1.+asymvt)
            rad_v=-vrmean*sqrt((radius-a))*a/radius*(1.+asymvr)
            if(s.gt.0.5) then
              VTas=(vorticity/2.)*(alpha*(a**(s+1.))/(radius**(s+1.))
     X            *cos(s*(theta+rotate*dtor)))
              rad_vs=-(vorticity*alpha/2.)*((a**(s+1.))/
     X           (radius**(s+1.)))*sin(s*(theta+rotate*dtor))
               print *, "in asy", s
               VTa=VTas+VTa
               rad_v=rad_vs+rad_v
C Uncomment these two lines to have just a disturbance
c               vta=vtas
c               rad_v=rad_vs
            endif
          endif
          if(radius.ge. 2.) then 
             u=meanu-vta*sin(theta)+rad_v*cos(theta)
             v=meanv+vta*cos(theta)+rad_v*sin(theta)
             ddz(i)=sqrt(u*u+v*v)
          else
             ddz(i)=-999.
          endif
           IF (abs(gdz(i)-dummy).gt.0.001) then
c             print *, angle/dtor,theta/dtor,vta,u,v,gdz(i)
             if(flag.eq.1 .or. flag.eq.3) then !vdad
               gdz(i) = gdz(i)*rg/RT
             elseif(flag.ne.2 .and. flag.ne.4) then ! not gbvtd or reflectivity
               print *, 'flag out of range'
             endif
           else
              print *, gdz(i), xd(i),yd(i),iori,jori,rg,rt
           endif
c           print *, gdz(i),ddz(i),xd(i),yd(i),xg,yg,rg,rt
           print *, gdz(i), vta, rad_v, theta/dtor, u, v, radius
         enddo
          return
          end

*****************************************************************
*   AVE Perform linear interpolation between two points
*****************************************************************

      subroutine ave(temp1,temp2,diff,totdiff,temp,dummy)
 
      implicit none
      real diff,totdiff,temp1,temp2,temp,dummy

      if(abs(temp1-dummy).ge.0.01) then
        if(abs(temp2-dummy).ge.0.01) then
          if(abs(totdiff).le.0.0001) then
             print *,'totdiff=0.'
             stop
          endif
          temp=(temp1*(totdiff-diff)+temp2*diff)/totdiff
        else
          temp=temp1
        endif
      elseif(abs(temp2-dummy).ge.0.01) then
        temp=temp2
      else
        temp=dummy
      endif  
      return
      end    


******************************************************************
*   GRID
******************************************************************
        subroutine grid(xc,yc,xdd,ydd,rmeter,tth,theta_t,status,
     X                  thetamax)

        implicit none
        real pi,rt,r,rmeter,xc,yc,theta_t,alpha
        integer i,status,thetamax
	real xdd(thetamax),ydd(thetamax) 
	real theta,tth(thetamax)

	pi=acos(-1.)
	rt=sqrt(xc**2+yc**2)
	r=rmeter
c        theta_t=atan2(yc,xc) !commented by wcl 98/10/26
        theta_t=amod(atan2(yc,xc)+2.*pi,2.*pi)  !put it back to 0-2pi
C A total of 90 points around a ring
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C tth(i) is the angle thetaprime (VDAD) or si (gbvtd) in the 
C hurricane coordinate
C theta is the math. angle counterclockwise from east to a point E
C in the gbvtd formulation, thetaprime=alpha+si (see Lee et al. 1997)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
	do i=1,thetamax 
           tth(i)=(i-1)*(360./float(thetamax))*pi/180.
           if(status.eq.1 .or. status.eq.3) then ! VDAD
             theta=amod((tth(i)+theta_t),2.*pi) 
           elseif(status.eq.2 .or. status.eq.4) then !gbvtd, calculate theta for a given si
             alpha=asin(r/rt*sin(tth(i)))  !law of sine
             theta=amod((tth(i)+alpha+theta_t),2.*pi) 
           endif
           xdd(i)=xc+r*cos(theta)
           ydd(i)=yc+r*sin(theta)
        enddo
        return
        end

      subroutine str2int(str,final_int)

      implicit none
      character*2 str
      integer*2 int1,int2,final_int

      int1 = iachar(str(1:1))
      int2 = iachar(str(2:2))
c      print *,'int:',int1,' ',int2
      final_int = int1*256+int2
      return
      end

      subroutine intswap(oldint,newint)
      implicit none
      integer*2 oldint,int1,int2,newint, intfix
      intfix = 256
      int1 = mod(oldint,intfix)
      int2 = (oldint - int1)/intfix
      newint = int1*intfix + int2
      return
      end

      subroutine maxvel(doppvel,numvel,max_pos,max_neg)
      implicit none
      integer numvel,i
      real doppvel(numvel), max_pos, max_neg
      
      max_pos = 1.;
      max_neg = -1.;
      do i=1,numvel
         if (doppvel(i).lt.max_neg) then
            max_neg = doppvel(i)
         else if (doppvel(i).gt.max_pos) then
            max_pos = doppvel(i)
         endif
      enddo

      if (max_neg.eq.-1.) then
         max_neg=-999;
      endif

      if (max_pos.eq.1.) then
         max_pos=-999;
      endif

      return
      end

