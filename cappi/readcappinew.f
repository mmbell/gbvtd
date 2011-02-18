       program READ_CAPPI_NEW
C----------------------------------------------------------------------------#
C   #1 READ CAA ROW DATA(unformat type in hp9000)
C   #2 FULL DOMAIN(241Km*241Km, grid interval 2Km) CAPPI
C   #3 USE "input.readcappi" AS PARAMETER INPUT FILE
C   #4 "TEMP.CAPPI" IS A TEMPFILE, DON'T KILL AT RUNNING PROGRAM
C   #5 MAIN PROGRAM CREATED BY                              FENG.  1992.3.30.
C----------------------------------------------------------------------------#
       implicit none
       include 'readcappi_para.h'
       CHARACTER infiledz*100,infilevr*100,outfile1*100,inputfile*100
       character*100 swpfile(maxnel)
       character test*2,code*2
       integer ier,lev,i,j,irad,jrad,k,ii,index
       Integer input_format,indx,parm_num,jj
       real diff,dz(maxnr,maxnth)
       real tempdz(nx,ny),tempve(nx,ny),tempvrel(nx,ny)
       real dbzmin,dbzmax,velmin,velmax,vel,deg
       real temp
       real tempsp,dzcaa(120,420)
       common /radar/irad,jrad,altrad
       data vel,deg /5.5,0.0/
       integer nargs,iargc

c**    different objection analysis for dBZ or Vr ------------------- 
c**    inflence radius = 2 km  

***** INITILIZE ARRAY ID ***** 
       do i=1,510
         id(i)=-999
       enddo

***** READ INPUT FILE ***** 

       nargs = iargc()
       IF ( nargs.GT.0 ) THEN
	  Call getarg(1, inputfile)
       ELSE
          print *, 'Please enter input file name:'
          read *, inputfile
       ENDIF
       OPEN(12,file=inputfile,status='old',form='formatted',
     +               err=998)
 100   CONTINUE
       READ(12,'(f4.1,1x,f4.1/i2/A2,6x,A2)',end=999)
     +          xc,yc,index,flag(1),flag(2)
       read(12,75)xmin,xmax,xgridsp
       read(12,75)ymin,ymax,ygridsp
       read(12,75)zmin,zmax,zgridsp
       read(12,75)dbzmin,dbzmax
       read(12,75)velmin,velmax
75     format(3f8.1)
       read(12,75)dummy,special,tempsp
       read(12,76)sf,cf
       read(12,77)irad,jrad,altrad
       read(12,79)adv_speed,adv_angle
       read(12,78)infiledz
       read(12,78)infilevr
       read(12,78)outfile1
76     format(3I8)
77     format(2I8,f8.3)
78     format(A100)
 79    format(f8.3,f8.3)
       ix=nint((xmax-xmin)/xgridsp+1.)
       iy=nint((ymax-ymin)/ygridsp+1.)
       iz=nint((zmax-zmin)/zgridsp+1.)
       adv_angle=mod(450.-adv_angle,360.)

       do i=1,iz
         zt(i)=zmin+(i-1)*zgridsp
       enddo

       if(index.eq.6)then
          tndo_flag=1
       else
          tndo_flag=0
       endif

       call assignid(xmin,xmax,
     X           ymin,ymax,zmin,zmax,irad,jrad,altrad)
            

***** SET UP GRID SIZE ***** 
       call grid(xc,yc,xd,yd,xgridsp,ygridsp,xmin,ymin)

***** READ INPUT FOR DZ*****
        code=flag(1)
	IF ((index.EQ.3).or.(index.eq.6)) THEN
           call readwsr88d(dz3d,swpfile,code)
        ELSEIF(index.eq.2) THEN
           call rankine(infiledz,dz3d,code)
        ELSEIF(index.eq.1 .or. index.eq.4) then
           call readcaa(infiledz,dz3d,ier,code)
           dummy=99.0
c Should not activate the following 3 lines for header information purposes
c        ELSEIF(index.eq.4) then
c           call readwenchau(infiledz,dz3d,1,nr,nth,nel)
c           dummy=-999.0
        ELSE
           print *,'UNKNOWN INPUT FORMAT'
           STOP
        endif

        do lev=1,nel
          DO j=1,nth(lev)
	    DO i=1,nr(lev)
              dz(i,j)=dz3d(i,j,lev)
            ENDDO
	  ENDDO
           
          IF (flag(1).eq.'dz' .or. flag(1).eq.'DZ' .or. flag(1).eq. 
     X         'ZZ' .or. flag(1).eq.'DM') THEN
            DO j=1,nth(lev)
              DO i=1,nr(lev)
                call dbzZe(dz(i,j),1,dummy) !CONVERT ZE TO LINEAR SCALE
	      ENDDO
            ENDDO
          ENDIF
          call fill(dz,dummy,maxnr,maxnth,nr(lev),nth(lev))
          DO j=1,nth(lev)
            DO i=1,nr(lev)
              dz3d(i,j,lev)=dz(i,j)
            ENDDO
          ENDDO
        enddo

        PRINT*,'DOING CAPPI.....'
        call intebilinear(dz3d,dzcappi)

        print *,'after intebilinear for dbz'

***** READ INPUT FOR VR*****
        code=flag(2)
	IF ((index.EQ.3).or.(index.eq.6)) THEN
           call readwsr88d(dz3d,swpfile,code)
        ELSEIF(index.eq.2) THEN
           call rankine(infilevr,dz3d,code)
        ELSEIF(index.eq.1) then
           call readcaa(infilevr,dz3d,ier,code)
           dummy=99.0
        ELSEIF(index.eq.4) then
           call readwenchau(infilevr,dz3d,2,nr,nth,nel)
           dummy=-999.0
        ELSE
           print *,'UNKNOWN INPUT FORMAT'
           STOP
        endif
        do lev=1,nel
          DO j=1,nth(lev)
            if(index.eq.4) then
              jj=mod((525-j),420)+1
            else
              jj=j
            endif
	    DO i=1,nr(lev)
              dz(i,jj)=dz3d(i,j,lev)
            ENDDO
	  ENDDO
           
          IF (flag(2).eq.'dz' .or. flag(2).eq.'DZ' .or. flag(2) 
     X        .eq. 'ZZ' .or. flag(2).eq.'DM') THEN
            DO j=1,nth(lev)
              DO i=1,nr(lev)
                call dbzZe(dz(i,j),1,dummy) !CONVERT ZE TO LINEAR SCALE
	      ENDDO
            ENDDO
          ENDIF
          call fill(dz,dummy,maxnr,maxnth,nr(lev),nth(lev))
          DO j=1,nth(lev)
            DO i=1,nr(lev)
              dz3d(i,j,lev)=dz(i,j)
            ENDDO
          ENDDO
          call relative(dz,vel,deg,nr(lev),nth(lev),dummy)
          do i=1,nr(lev)
            do j=1,nth(lev)
              vrel3d(i,j,lev)=dz(i,j)
            enddo
          enddo
        enddo

        PRINT*,'DOING CAPPI.....'
        call intebilinear(dz3d,vecappi)
        call intebilinear(vrel3d,vrelcappi)

***** SMOOTH THE DATA *****
       do k = 1,iz
          do j=1,iy
            do i=1,ix
              tempdz(i,j)=dzcappi(i,j,k)
              if(index.eq.4) then
                if(abs(tempdz(i,j)-99.0).le.0.001)tempdz(i,j)=special
              endif
              tempve(i,j)=vecappi(i,j,k)
              tempvrel(i,j)=vrelcappi(i,j,k)
            enddo
          enddo
          call smooth2d(tempdz,dummy,nx,ny,ix,iy)
          call smooth2d(tempve,dummy,nx,ny,ix,iy)
          call smooth2d(tempvrel,dummy,nx,ny,ix,iy)

          do j=1,iy
            do i=1,ix
              dzcappi(i,j,k)=tempdz(i,j)
              vecappi(i,j,k)=tempve(i,j)
              vrelcappi(i,j,k)=tempvrel(i,j)
            enddo
          enddo
       enddo

c***** WRITE CAPPI *****
       PRINT*,'WRITE CAPPI..... WRITE TO :',outfile1
       CALL writecappi(outfile1,dzcappi,vecappi,vrelcappi)
c       PRINT*,'File deleted !(TEMP.CAPPI)'
c       CLOSE(20)
c       PRINT*,'GOTO 100(doing again)'
c       GO TO 100
       stop
 996   print*, 'unable to open sweeplist'
       stop
 997   PRINT*,'FILE ERROR: TEMP.CAPPI    FILE ACCESS ERROR !'
       STOP
 998   PRINT*,'FILE OPEN ERROR : input.readcappi NOT FOUND'
       CLOSE(12)
C       CLOSE(20)
       STOP
 999   END


*****************************************************************
*   READDATA
*****************************************************************

       subroutine readwenchau(infilevr,temp,flag,nr,nth,nel)

       implicit none
       integer maxnr,maxnth,maxflds,maxnel
       parameter(maxnr=2000,maxnth=500,maxflds=15,maxnel=20)
       integer ik,i,j,n,nx,ny,iz,k,ix,iy,j1
       real xc,yc,pi,special,d,R_T
       integer flag,unit
       character infilevr*100
       real ve(maxnr,maxnth,maxnel),temp(maxnr,maxnth,maxnel)
       real dz(maxnr,maxnth,maxnel)
       integer nr(maxnel),nth(maxnel),nel
       integer*2 id(510)
      
       special=-999.
       unit=10
       pi=acos(-1.)
       OPEN(unit,FILE=infilevr,STATUS='OLD',FORM='FORMATTED')
       read(unit,116)id
116    format(10i8)
       ix=id(162)
       iy=id(167)
       iz=id(172)
       nel=iz
       do k=1,iz
         nr(k)=ix
         nth(k)=iy
         read(unit,117)ik
117      format(5x,i2)
         do j=1,iy
           read(unit,18)j1
c             print *, j1,j
           if(j1.ne.j) then
             print *, j1,j,' second index not match'
             stop
           endif
18         format(7x,i3)
           if(id(175).ge.1) then
             read(unit,19)id(176),id(177),id(178),id(179),
     X          (dz(i,j,k),i=1,ix)
           endif
           if(id(175).ge.2) then
             read(unit,19)id(181),id(182),id(183),id(184),
     X        (ve(i,j,k),i=1,ix)
           endif
           if(flag.eq.1) then
             do i=1,ix
               temp(i,j,k)=dz(i,j,k)
             enddo
           else
             do i=1,ix
               temp(i,j,k)=ve(i,j,k)
             enddo
           endif
         enddo
       enddo
19     format(4a2,/,(8E10.3))
       close(unit)
       return
       end


C--------------------------------------------------------------------
      subroutine relative(temp,vel,deg,nr,nth,dummy)

      implicit none
c      include 'readcappi_para.h'
      integer i,j,nr,nth
      real temp(nr,nth),vel,deg,pi,azn,deltaaz,comp,dummy,special


      pi=acos(-1.)
      azn=nth/(2.*pi)
      deltaaz=1/azn
      do j=1,nth
        comp=vel*cos(deg/180.*pi-deltaaz*(real(j-1)))
        do i=1,nr
          if(abs(temp(i,j)-dummy).gt.0.0001)then 
            if(temp(i,j).gt.2000. .or. temp(i,j) .le. -2000.) then
              print *, i,j,temp(i,j),comp,deg,deltaaz*(real(j-1))
            endif
            temp(i,j)=temp(i,j)-comp
          endif
        enddo
      enddo
      return
      end


***** ASSIGN ID ***** 
       subroutine assignid(xmin,xmax,
     X           ymin,ymax,zmin,zmax,irad,jrad,altrad)

       implicit none
       include 'readcappi_para.h'
       integer irad,jrad
c       real altrad

       if(tndo_flag.eq.1) then !Need to scale the data for center-finding
          xmin=xmin*10
          xmax=xmax*10
          ymin=ymin*10
          ymax=ymax*10
          zmin=zmin*10
          zmax=zmax*10
          xgridsp=xgridsp*10
          ygridsp=ygridsp*10
          zgridsp=zgridsp*10
       endif

       id(40)=90
       id(68)=sf
       id(69)=cf
       id(160)=int(xmin*sf)
       id(161)=int(xmax*sf)
       id(162)=ix
       id(163)=xgridsp*1000. !convert xgridsp from km to m
       id(164)=1
       id(165)=int(ymin*sf)
       id(166)=int(ymax*sf)
       id(167)=iy
       id(168)=ygridsp*1000. !convert ygridsp from km to m
       id(169)=2
       id(170)=int(zmin*1000.) !convert zmin from km to m
       id(171)=int(zmax*1000.) !convert zmax from km to m
       id(172)=iz
       id(173)=zgridsp*1000. !convert zgridsp from km to m
       id(174)=3
       id(175)=3
c       id(176)='DZ'
       read(flag(1), 1)id(176)
1      format(a2)     
       print *, id(176)
       call str2int('  ',id(177))
       call str2int('  ',id(178))
       call str2int('  ',id(179))
       id(180)=1
c       id(181)='VV'
       read(flag(2),1)id(181)
       print *, id(181)
       call str2int('  ',id(182))
       call str2int('  ',id(183))
       call str2int('  ',id(184))
       id(185)=1
       call str2int('VR',id(186))
       print *, id(186)
       call str2int('EL',id(187))
       call str2int('  ',id(188))
       call str2int('  ',id(189))
       id(190)=1
       id(309)=irad*sf
       id(310)=jrad*sf
       id(311)=altrad*1000

       if(tndo_flag.eq.1) then !Unscale the data for rest of program
          xmin=xmin/10
          xmax=xmax/10
          ymin=ymin/10
          ymax=ymax/10
          zmin=zmin/10
          zmax=zmax/10
          xgridsp=xgridsp/10
          ygridsp=ygridsp/10
          zgridsp=zgridsp/10
       endif

       return
       end

       subroutine readcaa(file1,dz3d,ier,code)

       include 'readcappi_para.h'
       CHARACTER FILE1*100,code*2
       INTEGER*2 IYER,IMON,IDAY,IHOR,IMIN,ISEC,KK
       INTEGER*2 IDEXZ,IDEXW,ISLPZ,ISLPW,ISRDZ,ISRDW
       INTEGER*2 IDZ(120),LO1,LO2,LA1,LA2,IELEV(20)
       INTEGER*2 ITIME(420),IZ1(120,420),ITIME1(420)
       real rad_lat,rad_long
       real tim(120,420)
       integer lnblnk,i,numgate,numray
       data numgate,numray /120,420/
       data rad_lat, rad_long /25.0694,121.2083/

       common /header/ iyer,imon,iday,ihor,imin,isec
c** read one level CAA doppler radar dBZ(zn*.dat) or Vr(wn*.dat)
c** file(input)          : data file name
c** lev(input)           : data's level(its elevation angle)
c** elev(output)         : elevation angle(degrees)
c** tim(420)(output)     : azimuthal data's time(unit: sec)
c** dz(120,420)(output)  : dBZ data(no data =dummy)
c** tim(120,420)(output) : time  data(no data =dummy)
c** dummy(input)         : no data value
c** ier(output)          : ier=0 is normal , ier=1 is error

c********************* Converting radar lat and long into deg,min,sec
       id(33)=int(rad_lat)
       temp=(rad_lat-real(id(33)))*60.
       id(34)=int(temp)
       id(35)=(int(temp*60.-real(id(34)))*60.)*sf
       id(36)=int(rad_long)
       temp=(rad_long-real(id(36)))*60.
       id(37)=int(temp)
       id(38)=(int(temp*60.-real(id(37)))*60.)*sf
       ier=0
       GO TO 9
8      PRINT*,' ',file1,' NOT FOUND'
       ier=1
       return
9      CONTINUE
       file1=file1(:lnblnk(file1))
c       OPEN(10,FILE=FILE1,STATUS='OLD',FORM='UNFORMATTED',
c     +      READONLY,ERR=8)
       OPEN(10,FILE=FILE1,STATUS='OLD',FORM='UNFORMATTED',ERR=8)
       READ(10)IYER,IMON,IDAY,IHOR,IMIN,ISEC,LO1,LO2,LA1,LA2,
     *          ISLPZ,ISRDZ,IDEXZ,ISLPW,ISRDW,IDEXW,KK
       idate=iyer*10000+imon*100+iday
       ntim=ihor*10000+imin*100+isec
       nel=kk
c       if( (lev.lt.1).or.(lev.gt.kk) )then
c	print*,'level is error ',lev,' max level= ',kk
c	ier=1
c	return
c       endif
       RLO=LO1+LO2/10000.
       RLA=LA1+LA2/10000.
       SLOPZ=ISLPZ/1000.
       SORDZ=ISRDZ/100.
       SLOPW=ISLPW/1000.
       SORDW=ISRDW/100.
       call str2int('CA',id(13))
       call str2int('A ',id(14))
       call str2int('  ',id(15))
       call str2int('PP',id(16))
       call str2int('I ',id(17))
       call str2int('  ',id(18))
       id(21)=iyer
       id(22)=imon
       id(23)=iday
       id(24)=ihor
       id(25)=imin
       id(26)=isec
       id(33)=int(la1)
       id(34)=int(la2/10000.*60.)
       id(35)=int((la2-(int(la2/100)*100))/10000.*3600.*sf)
       id(36)=int(lo1)
       id(37)=int(lo2/10000.*60.)
       id(38)=int((lo2-(int(lo2/100)*100))/10000.*3600.*sf)
       id(67)=-999
       id(116)=iyer
       id(117)=imon
       id(118)=iday
       id(119)=ihor
       id(120)=imin
       id(121)=isec

       READ(10)IELEV
       if(kk .gt. maxnel) then
          print *,'maxlev exceeds maxnel in CAAREAD'
          stop
       endif
       DO 666 K=1,kk
         nr(k)=numgate
         nth(k)=numray
	 startang(k)=0.
	 stopang(k)=360.
         if(nth(k).eq.0) stop
	 beamsep(k)=360.0/nth(k)
         gatesp(k)=1. 
         elev(k)=ielev(k)/100.
         DO 630 J=1,nth(k)
           angle(j,k)=beamsep(k)*(j-1)
           READ(10)ITIME(J),IDZ
           DO 20 I=1,nr(k)
             IZ1(I,J)=IDZ(I)
             tim(i,j)=itime(j)/10.
 20        CONTINUE
 630     CONTINUE
         IF(K.EQ.kk-1) THEN
           DO 40 J=1,nth(k)
             ITIME1(J)=ITIME(J)
 40        CONTINUE
         ENDIF
         DO 60 J=1,nth(k)
           DO 50 I=1,nr(k)
             IF(lev.EQ.kk) THEN
               IF(ITIME1(J).EQ.ITIME(J)) THEN
	         dz3d(i,j,k)=dummy
                 GO TO 60
               ENDIF
             ENDIF
             IF( (IZ1(I,J).EQ.0).or.(IZ1(i,j).eq.255) ) THEN
	       dz3d(i,j,k)=dummy
             else
               dz3d(I,J,k)=(IZ1(I,J)-IDEXZ)*SLOPZ+SORDZ
             endif
             if(dz3d(i,j,k).gt.70.) dz3d(i,j,k)=dummy 
 50        continue
 60      continue
 666   continue
       CLOSE(10)
       return
       END
C--------------------------------------------------------------------
      subroutine rankine(file1,dz3d,code)

      implicit none
      include 'readcappi_para.h'
      integer i,j,lev,k,numr,numth,numz,kk
      character file1*100, code*2
      real rad_lat,rad_long,temp
       data rad_lat, rad_long /25.0694,121.2083/

      elev(1)=0.0
c********************* Converting radar lat and long into deg,min,sec
       id(33)=int(rad_lat)
       temp=(rad_lat-real(id(33)))*60.
       id(34)=int(temp)
       id(35)=(int(temp*60.-real(id(34)))*60.)*sf
       id(36)=int(rad_long)
       temp=(rad_long-real(id(36)))*60.
       id(37)=int(temp)
       id(38)=(int(temp*60.-real(id(37)))*60.)*sf
      OPEN(101,FILE=FILE1,STATUS='OLD',FORM='FORMATTED',ERR=8)
      read(101,103) numr,numth,nel
      do i=1,nel
        nr(i)=numr
        nth(i)=numth
	startang(i)=0.
	stopang(i)=360.
        if(nth(i).eq.0) then
           print *, 'Number of ray is 0.'
           stop
        endif
	beamsep(i)=360.0/real(nth(i))
        do j=1,nth(i)
           angle(j,i)=startang(i)+real(j-1)*beamsep(i)
        enddo
        gatesp(i)=1. 
      enddo
      do kk=1,nel
        do i=1,nr(kk)
          read(101,102)code,k,(dz3d(i,j,kk),j=1,nth(kk))
          print *,code,k
        enddo
      enddo
      return
 102    format(a2,i4,/(8E10.3))
 103    format(3I8)
 8    print *, 'unable to open input file: ', file1
      return
      end
*******************************************************************
      SUBROUTINE READWSR88D(dz3d,swpfile,test)

      implicit none
       include 'structures.h'
       include 'readcappi_para.h'
       character*100 swpfile(maxnel)
       integer ier,lev,i,j,irad,jrad,k,ii,index,counter
       Integer input_format,indx1,indx2,parm_num
       integer newcount,angleflag(maxnth,maxnel)
       real tempangle(maxnth,maxnel),tempdz(maxnr,maxnth,maxnel)
       integer temphour(maxnth,maxnel),tempminute(maxnth,maxnel)
       integer tempsecond(maxnth,maxnel),dt
       character*2  test
       character*100 test1
       real diff,temp,tempsp,testang
        Type (VOLD_INFO) vold
        Type (RADD_INFO) radd
        Type (CELV_INFO) celv
        Type (CFAC_INFO) cfac
        Type (PARM_INFO) parm(maxflds)
        Type (SWIB_INFO) swib
        Type (RYIB_INFO) ryib(maxnth)
        Type (ASIB_INFO) asib(maxnth)
        Type (RDAT_INFO) rdat(maxnth)
	Type (PARM_NAME) pname

c EQUIVALENCE doesnt work
c        EQUIVALENCE(pname%int,pname%name)

C          PROCESS DZ
C	   READ SWEEP FILE UNTIL LIST IS EXHAUSTED
           lev=1
36         if(test.eq.'DZ' .or. test.eq.'dz' .or. test.eq.'ZZ'
     X        .or. test.eq.'DM') then
             Read (12,*,end=37),swpfile(lev)
	   Call sweepread(swpfile(lev),vold,radd,celv,cfac,parm,swib,
     $                    ryib,asib,rdat,id(176))
           else
             test1=swpfile(lev)
             if(test1(1:3).eq.'end') go to 37
	   Call sweepread(swpfile(lev),vold,radd,celv,cfac,parm,swib,
     $                    ryib,asib,rdat,id(181))
           endif

	   print*, 'processing ',swpfile(lev)
C          ****************** assign lat,long,alt into cedric header, only
C          ****************** do it on the first sweep
           if(lev.eq.1) then
             call str2int(radd%rad_name(1:2),id(13))
             call str2int(radd%rad_name(3:4),id(14))
             call str2int(radd%rad_name(5:6),id(15))
c             id(13)=radd%rad_name(1:2)
c             id(14)=radd%rad_name(3:4)
c             id(15)=radd%rad_name(5:6)
             id(33)=int(radd%radar_lat)
             temp=(radd%radar_lat-real(id(33)))*60.
             id(34)=int(temp)
             id(35)=(int((temp-real(id(34)))*60.))*sf
             id(36)=int(radd%radar_lon)
             temp=(radd%radar_lon-real(id(36)))*60.
             id(37)=int(temp)
             id(38)=(int((temp-real(id(37)))*60.))*sf
             id(39)=int(radd%radar_alt*1000.)
             id(116)=vold%year
             id(117)=vold%mon
             id(118)=vold%day
             id(119)=vold%hour
             id(120)=vold%min
             id(121)=vold%sec
           endif
           id(122)=vold%year
           id(123)=vold%mon
           id(124)=vold%day
           id(125)=vold%hour
           id(126)=vold%min
           id(127)=vold%sec
                        
	   nr(lev)=celv%total_gates
           if(nr(lev).gt.maxnr) then
              print *,'maxnr too small!!!'
              stop
           endif
	   nth(lev)=swib%num_rays
           if(nth(lev).gt.maxnth) then
              print *,'maxnth too small!!!'
              stop
           endif
	   nel=lev
	   startang(lev)=swib%start_ang
	   stopang(lev)=swib%stop_ang
           if(nth(lev).eq.0) stop
C--  Watch out, this is a big source of error especially in PPI scans
C--  Modified by wcl 4-27-00 attempt to handle PPI data by checking the
C--    direction (clockwise or counterclockwise) of the scan then compute
C--    the beamwidth
C--------------------------The original code
c	   beamsep(lev)=amod((stopang(lev)+360.0-startang(lev)),360.)
c     $                  /nth(lev) 
c           if(beamsep(lev).le.0.1) then
c             print *, 'beamwidth(',lev,') too small', beamsep(lev)
c  	     beamsep(lev)=(amod((stopang(lev)+360.0-startang(lev)),360.)
c     $                  +360.)/nth(lev)
c             print *, 'new beamwidth=', lev, beamsep(lev)
c           endif 
C--------------------------End the original code

           if (radd%scan_mode .eq. 1) then   !PPI mode
             testang=amod(ryib(2)%azimuth+360., 360.) + 360. -
     $               amod(ryib(1)%azimuth+360., 360.)
             if(testang .gt. 0. .and. testang .lt. 180.) then !clockwise
	        beamsep(lev)=amod((stopang(lev)+360.0-startang(lev)),
     $                  360.)/nth(lev)
                counter=1 
             else   !counterclockwise
                beamsep(lev)=amod((startang(lev)+360.0-stopang(lev)),
     $                  360.)/nth(lev)
                counter=2
             endif
           endif
                


C	   GET THE PARAMETER NUMBER
	   DO i=1,radd%num_param_desc
              if(test.eq.'DZ' .or. test.eq.'dz' .or. test.eq.'ZZ'
     X             .or. test.eq.'DM') then
                 pname%int(1)=id(176)
                 pname%int(2)=id(177)
                 pname%int(3)=id(178)
                 pname%int(4)=id(179)
              else
                 pname%int(1)=id(181)
                 pname%int(2)=id(182)
                 pname%int(3)=id(183)
                 pname%int(4)=id(184)
              endif
              pname%name(1:1) = achar(pname%int(1))
              pname%name(2:2) = achar(pname%int(2))
              pname%name(3:3) = achar(pname%int(3))
              pname%name(4:4) = achar(pname%int(4))
              indx1=index(pname%name,' ')-1
	      indx2 = index(parm(i)%parm_name,CHAR(0))-1
	      IF (indx1.EQ.indx2) THEN
	         IF (parm(i)%parm_name(1:indx2).EQ.
     $               pname%name(1:indx1)) THEN
		    parm_num=i
	         ENDIF
	      ENDIF
	   ENDDO
c           if(test.eq.'VR' .or. test.eq.'vr') then
c	     DO j=1,nth(lev)
c	      DO i=1,nr(lev)
c		 print *, i,j,rdat(j)%data(i)
c	      ENDDO
c	     ENDDO
c           endif

C          THROW OUT RAYS WHERE THE ACTUAL ELEVATION ANGLE (RYIB) IS
C          +/- .25 DEGREE DIFFERENCE THAN THE FIXED ELEVATION 
C          IF WE ARE NOT USING DOW TORNADO DATA

	   DO j=1,nth(lev)
             if(tndo_flag.eq.0) then
	      diff=abs(swib%fixed_ang-ryib(j)%elevation)
	      IF (diff.GT..25) THEN
                 print *,j,swib%fixed_ang, ryib(j)%elevation
	         DO i=1,nr(lev)
		    rdat(j)%data(i)=parm(parm_num)%baddata_flag
		 ENDDO
	      ENDIF
             endif
C-- WCL, 98/11/5, fill azimuth angle array for each beam to solve 
C--   missing beam problem
     
              angle(j,lev)=amod(ryib(j)%azimuth+360.
     X            +cfac%c_azimuth, 360.)
              hour(j,lev)=ryib(j)%hour
              minute(j,lev)=ryib(j)%min
              second(j,lev)=ryib(j)%sec
     	   ENDDO

C          STORE DATA IN 3 DIM ARRAY
	   DO j=1,nth(lev)
	      DO i=1,nr(lev)
C                 Uncomment this if block if using MM wavelength data
C                 since the bad flag is 0

c                 if(abs(rdat(j)%data(i)) .le. 0.001) then
c                    dz3d(i,j,lev)=-999.0
c                 else
     		    dz3d(i,j,lev)=rdat(j)%data(i)
c                 endif
	      ENDDO
	   ENDDO

c         Find multiple beams
           do j=2,nth(lev)
              if (angle(j,lev).eq.angle(j-1,lev)) then
                 !Flag angle
                 angleflag(j-1,lev)=1
              else
                 !Good angle
                 angleflag(j-1,lev)=0
              endif
           enddo

c         Map old beams to new
           newcount=1
           do j=1,nth(lev)
              if (angleflag(j,lev).eq.0) then
                 tempangle(newcount,lev)=angle(j,lev)
                 temphour(newcount,lev)=hour(j,lev)
                 tempminute(newcount,lev)=minute(j,lev)
                 tempsecond(newcount,lev)=second(j,lev)
                 do i=1,nr(lev)
                  tempdz(i,newcount,lev)=dz3d(i,j,lev)
                 enddo
                 newcount=newcount+1
              endif
           enddo
           newcount=newcount-1

           do j=1,newcount
              angle(j,lev)=tempangle(j,lev)
              hour(j,lev)=temphour(j,lev)
              minute(j,lev)=tempminute(j,lev)
              second(j,lev)=tempsecond(j,lev)
              if (j.gt.1) then
                 dt = hour(j,lev)*3600+minute(j,lev)*60+second(j,lev)
     X           -hour(j-1,lev)*3600-minute(j-1,lev)*60-second(j-1,lev)
                 if(dt.gt.5) then
c Problem with times and possible advection, go back and edit! MMB 10/02
                    print *,'Time difference ',dt,' at ',j,lev
                    stop
                 endif
              endif
             do i=1,nr(lev)
              dz3d(i,j,lev)=tempdz(i,j,lev)
             enddo
           enddo
           nth(lev)=newcount

C          GATE SPACING (IN KM)
           gatesp(lev)=(celv%gate_spacing(2)-celv%gate_spacing(1))/1000.

C          RANGE TO FIRST GATE
           range_first=celv%gate_spacing(1)/1000

C          ELEVATION ANGLE
	   elev(lev)=swib%fixed_ang
           if((test.eq.'dz' .or. test.eq.'DZ' .or. test.eq.'ZZ'
     X        .or. test.eq.'DM') .and. lev .le. 2) then
             Read (12,*,end=37),swpfile(lev)
           endif
c           if(lev.eq.iz) return
	   lev=lev+1
           go to 36
37	   continue  ! SWEEP LIST EXHAUSTED!
           swpfile(lev)='end'

      return
      end

**************************************************
      subroutine fill(temp,dummy,maxnr,maxnth,nr,nth)

      integer i,j,k,maxnr,maxnth,nr,nth
      real temp(maxnr,maxnth),temp1,dummy
    
      do j=2,nth-1
        do i = 2,nr-1
          k=0
          temp1=0.
          if(abs(temp(i,j)-dummy).lt.0.001) then
            if(abs(temp(i-1,j)-dummy).gt.0.001) then
              temp1=temp1+temp(i-1,j)
              k=k+1
            endif
            if(abs(temp(i+1,j)-dummy).gt.0.001) then
              temp1=temp1+temp(i+1,j)
              k=k+1
            endif
            if(abs(temp(i,j-1)-dummy).gt.0.001) then 
              temp1=temp1+temp(i,j-1)
              k=k+1
            endif
            if(abs(temp(i,j+1)-dummy).gt.0.001) then
              temp1=temp1+temp(i,j+1)
              k=k+1
            endif
          endif
          if(k.gt.2) temp(i,j)=temp1/real(k)
        enddo
      enddo
      return
      end
**************************************************
      subroutine smooth2d(temp,dummy,nx,ny,ix,iy)

      integer i,j,k,nx,ny,ix,iy
      real temp(nx,ny),temp1,dummy
    
      do j=2,iy-1
        do i = 2,ix-1
          k=0
          temp1=0.
          if(abs(temp(i,j)-dummy).gt.0.001) then
            temp1=temp1+4*temp(i,j)
            k=k+4
          endif
          if(abs(temp(i-1,j)-dummy).gt.0.001) then
            temp1=temp1+temp(i-1,j)
            k=k+1
          endif
          if(abs(temp(i+1,j)-dummy).gt.0.001) then
            temp1=temp1+temp(i+1,j)
            k=k+1
          endif
          if(abs(temp(i,j-1)-dummy).gt.0.001) then 
            temp1=temp1+temp(i,j-1)
            k=k+1
          endif
          if(abs(temp(i,j+1)-dummy).gt.0.001) then
            temp1=temp1+temp(i,j+1)
            k=k+1
          endif
          if(k.gt.2) temp(i,j)=temp1/real(k)
        enddo
      enddo
      return
      end
**************************************************
      Subroutine intebilinear(dz3d,dzcappi)

      implicit none
      include 'readcappi_para.h'
      integer i,j,k,ngate1,ngate2,nh,nazimuth,irad,jrad,nazimuth_1
      integer nazplus1_1,oobflag,dirindex(nel)
      real hrange,range,el,el1,el2,el_deg1,el_deg2
      real az_deg,deltaaz,deltaaz_1,dra1,dra2
      real pi,daz,del,temp1,temp2,ave1,ave2,ave3,el_deg,daz_1
      real deltael,ave4,ave5,dx,dy,deltaz,azimuth,diffaz
      real sepi,sepi1,rprime,earth
      real x1,y1,x2,y2,range1,range2,az_deg1,az_deg2
      integer nazplus1,d,deltatime
c      data range_first /1.0/
       common /radar/irad,jrad,altrad

c     Get directions
      do d=1,nel
         call direction(d,dirindex(d))
      enddo

      pi=acos(-1.)
      rprime=(4./3.*6374.)
      x1=0
      y1=0
      x2=0
      y2=0
      do k=1,iz
        deltaz=(zt(k)-altrad)
        do j=1,iy
          dy=(j-jrad)*ygridsp
          do i=1,ix
            dzcappi(i,j,k)=dummy
            dx=(i-irad)*xgridsp
            hrange=sqrt(dx*dx+dy*dy)
            range=sqrt(hrange*hrange+deltaz*deltaz)
c--take into account earth curvature, wcl 6/13/00, from Rinehart
            earth=sqrt(hrange*hrange+rprime*rprime)-rprime
            el=atan2(deltaz-earth,hrange)
            el_deg=el/pi*180
C-comment for 1 ppi only
            call findz(el_deg,elev,nh,nel)
c            nh=1
            if(nh.ne.0) then
              if((gatesp(nh).lt.0.01) .and.(nh.ne.(nel+1))) then
                print *, i,j,k
                stop
              endif
              azimuth=mod(2.5*pi-atan2(dy,dx),2.*pi) !to meteorological coord.
              az_deg=azimuth*(180./pi)

              if (nh.lt.nel+1) then
                 call search(az_deg, nh, nazimuth,sepi,
     X                sepi1,dirindex(nh),oobflag)
c Now we know az,lev so we can advect origin
               if (oobflag.eq.0) then
                  call calc_time(nazimuth,nh,deltatime)
                  x2= adv_speed/1000*deltatime*cos(adv_angle*pi/180.)
                  y2= adv_speed/1000*deltatime*sin(adv_angle*pi/180.)
               endif
            else
               call search(az_deg, nel, nazimuth,sepi,
     X              sepi1,dirindex(nel),oobflag)
               if (oobflag.eq.0) then
                  call calc_time(nazimuth,nel,deltatime)
                  x2= adv_speed/1000*deltatime*cos(adv_angle*pi/180.)
                  y2= adv_speed/1000*deltatime*sin(adv_angle*pi/180.)
               endif
            endif

            if(nh.gt.1) then
               call search(az_deg, nh-1, nazimuth_1,
     X              sepi,sepi1,dirindex(nh-1),oobflag)
               if (oobflag.eq.0) then
                  call calc_time(nazimuth_1,nh-1,deltatime)
                  x1= adv_speed/1000*deltatime*cos(adv_angle*pi/180.)
                  y1= adv_speed/1000*deltatime*sin(adv_angle*pi/180.)
               endif
            endif

            hrange = sqrt((dx-x1)**2+(dy-y1)**2)
            earth=sqrt(hrange*hrange+rprime*rprime)-rprime
            el1=atan2(deltaz-earth,hrange)
            el_deg1=el1/pi*180
            range1 = sqrt(hrange*hrange+deltaz*deltaz)
            azimuth=mod(2.5*pi-atan2(dy-y1,dx-x1),2.*pi) !to meteorological coord.
            az_deg1=azimuth*(180./pi)

            hrange = sqrt((dx-x2)**2+(dy-y2)**2)
            earth=sqrt(hrange*hrange+rprime*rprime)-rprime
            el2=atan2(deltaz-earth,hrange)
            el_deg2=el2/pi*180
            range2 = sqrt(hrange*hrange+deltaz*deltaz)
            azimuth=mod(2.5*pi-atan2(dy-y2,dx-x2),2.*pi) !to meteorological coord.
            az_deg2=azimuth*(180./pi)
            
            el_deg=(el_deg1+el_deg2)/2

C---- take into account range to the first gate  (wcl 4/29/00)
c              ngate=int((range-range_first)/gatesp(nh))
C---- take into account changing gate spacing within a volume
C---- For Josh Wurman !!!!! (wcl 2/13/01)
              if(nh.ne.1) then
                ngate1=int((range1-range_first)/gatesp(nh-1))
                if(((range1-range_first)/gatesp(nh-1)).gt.0
     X               .and.ngate1.lt.1) then
                   ngate1=1
                endif
                dra1=range1-(ngate1*gatesp(nh-1)+range_first)
              endif
              if(nh.ne.nel+1) then
                ngate2=int((range2-range_first)/gatesp(nh))
                dra2=range2-(ngate2*gatesp(nh)+range_first)
              endif
c              if((ngate1.gt.1 .and. ngate1.lt.nr(nh-1)).and.
c     X           (ngate2.gt.1 .and. ngate2.lt.nr(nh))) then !process ngate1 within range
c                dra=range-(ngate*gatesp(nh)+range_first)
c                azimuth=mod(2.5*pi-atan2(dy,dx),2.*pi) !to meteorological coord.
c	        az_deg=azimuth*(180./pi)
                if (nh.lt.nel+1) then
                   call search(az_deg2, nh, nazimuth,sepi,
     X                  sepi1,dirindex(nh),oobflag)
                else
                   call search(az_deg2, nel, nazimuth,sepi,
     X                  sepi1,dirindex(nel),oobflag)
                endif


                if (oobflag.eq.0) then   !Go ahead and interpolate
                   if(nh.lt.nel+1) then
                      nazplus1=mod(nazimuth,nth(nh))+1 !put azimuth back to the range
                      if(nazplus1.gt.nth(nh)) print *,i,j,k
                   else
                      nazplus1=mod(nazimuth,nth(nel))+1 !put azimuth back to the range
                      if(nazplus1.gt.nth(nel)) print *,i,j,k
                   endif
                deltaaz=sepi+sepi1
                daz=sepi
                if(sepi.lt.0. .or. sepi1.lt.0.) then
                  print *, i,j,k,sepi,sepi1
                  stop
                endif
c azsweep and azsweep_1 are necessary because the starting angles
c of wsr88d are different for each sweep, this is a more generalized
c case
                if(nh.gt.1 .and. nh.lt.nel+1) then
c                   if (nh.eq.2) print *,'5', az_deg1, dx, dy
                  call search(az_deg1, nh-1, nazimuth_1,
     X              sepi,sepi1,dirindex(nh-1),oobflag)
                  if (oobflag.eq.0) then !Keep going

                  deltaaz_1=sepi+sepi1
                  nazplus1_1=mod(nazimuth_1,nth(nh-1))+1 !put azimuth back to the range
                  daz_1=sepi
                  del=el_deg-elev(nh-1)
                  deltael=elev(nh)-elev(nh-1)
c average across azimuth on gate ngate and elev nh-1
                  temp1=dz3d(ngate1,nazimuth_1,nh-1)
                  temp2=dz3d(ngate1,nazplus1_1,nh-1)
                  call ave(temp1,temp2,daz_1,deltaaz_1,ave1,dummy)
c average across azimuth on gate ngate and elev nh
                  temp1=dz3d(ngate2,nazimuth,nh)
                  temp2=dz3d(ngate2,nazplus1,nh)
                  call ave(temp1,temp2,daz,deltaaz,ave2,dummy)
c average across elevation on gate ngate1 and elev nh-1 and ngate2 on nh
                  call ave(ave1,ave2,del,deltael,ave3,dummy)
c average across azimuth on gate ngate+1 and elev nh-1
                  temp1=dz3d(ngate1+1,nazimuth_1,nh-1)
                  temp2=dz3d(ngate1+1,nazplus1_1,nh-1)
                  call ave(temp1,temp2,daz_1,deltaaz_1,ave1,dummy)
c average across azimuth on gate ngate+1 and elev nh
                  temp1=dz3d(ngate2+1,nazimuth,nh)
                  temp2=dz3d(ngate2+1,nazplus1,nh)
                  call ave(temp1,temp2,daz,deltaaz,ave2,dummy)
c average across elevation on gate ngate+1 and elev nh-1 and nh
                  call ave(ave1,ave2,del,deltael,ave4,dummy)
c average across range 
                  call ave(ave3,ave4,dra1,gatesp(nh),ave5,dummy)
c assign value to the grid
                  dzcappi(i,j,k)=ave5
                endif
               endif    !oobflag end
                if(nh.eq.1 .and. (range2*(el2-(elev(1)/180.*pi)))
     X                 .lt.zgridsp) then
c average across azimuth on gate ngate and elev nh
                  temp1=dz3d(ngate2,nazimuth,nh)
                  temp2=dz3d(ngate2,nazplus1,nh)
                  call ave(temp1,temp2,daz,deltaaz,ave1,dummy)
c                 print *, ngate, nazimuth, temp1,temp2,ave1,i,j,k
c average across azimuth on gate ngate+1 and elev nh
                  temp1=dz3d(ngate2+1,nazimuth,nh)
                  temp2=dz3d(ngate2+1,nazplus1,nh)
                  call ave(temp1,temp2,daz,deltaaz,ave2,dummy)
c average across range 
                  call ave(ave1,ave2,dra2,gatesp(nh),ave5,dummy)
c                 print *, ngate+1, nazplus1, temp1,temp2,ave2,ave5
c assign value to the grid
                  dzcappi(i,j,k)=ave5
                endif
C when grid is above the highest sweep
                if(nh.eq. nel+1 .and. (range1*(el1-(elev(nel)/180.*pi)))
     X                 .lt.zgridsp) then
c average across azimuth on gate ngate and elev nh
                  temp1=dz3d(ngate1,nazimuth,nel)
                  temp2=dz3d(ngate1,nazplus1,nel)
                  call ave(temp1,temp2,daz,deltaaz,ave1,dummy)
c average across azimuth on gate ngate+1 and elev nh
                  temp1=dz3d(ngate1+1,nazimuth,nel)
                  temp2=dz3d(ngate1+1,nazplus1,nel)
                  call ave(temp1,temp2,daz,deltaaz,ave2,dummy)
c average across range 
                  call ave(ave1,ave2,dra1,gatesp,ave5,dummy)
c assign value to the grid
                  dzcappi(i,j,k)=ave5
                endif
              endif
            endif
c            endif   ! This ends the oobflag condition
          enddo
c          print *, i,j,k,(dzcappi(i,j,k),i=1,ix)
        enddo
      enddo
      return
      end
******************************************************************
      subroutine direction(nh,index)

      implicit none
      include 'readcappi_para.h'
      integer nh,index,test
      real sep

      test=nth(nh)/2
       sep=angle(test,nh)-angle(test-10,nh)
c      sep=angle(2,nh)-angle(1,nh)
      if(sep.gt.0.) then
        if(sep.gt.180.) then ! counterclockwise
          index=2        ! counterclockwise, cross 0
        else     !clockwise
          index=1
        endif
      elseif(sep.le.-180.) then !clockwise, cross 0      
        index=1
      else !counterclockwise
        index=2
      endif
      return
      end
         
******************************************************************
      subroutine search(az_deg, nh,nazimuth,sepi,sepi1,dindex,oobflag)

      implicit none
      include 'readcappi_para.h'
      real az_deg, sepi,sepi1,oobtest1,oobtest2
      integer nh, nazimuth, i, dindex,oobflag


      if (az_deg.eq.-999) then
         oobflag=1
         return
      endif
      if (nh.eq.1) then
c Stop here!
      endif

      oobtest1 = az_deg-angle(nth(nh),nh)
      oobtest2 = az_deg-angle(1,nh)
      if((abs(angle(nth(nh),nh)-angle(1,nh)).le.7.).or.
     X     (abs(angle(nth(nh),nh)-angle(1,nh)).gt.352.)) then
         oobflag=0
      elseif((angle(1,nh).gt.210.).and.(angle(nth(nh),nh).lt.150.)) then
         if(oobtest1.gt.0. .and. oobtest2.lt.0.) then
            oobflag=1
            return
         else
            oobflag=0
         endif
      elseif((angle(nth(nh),nh).gt.210.).and.(angle(1,nh).lt.150.)) then
         if(oobtest1.lt.0. .and. oobtest2.gt.0.) then
            oobflag=1
            return
         else
            oobflag=0
         endif
      elseif(oobtest1.gt.0. .and. oobtest2.gt.0.) then !Out of bounds
         oobflag=1
         return
c         angle(i,nh)=angle((nth(nh)/2),nh)
c         angle(i+1,nh)=angle((nth(nh)/2)+1,nh)
      elseif (oobtest1.lt.0. .and. oobtest2.lt.0.) then !Also oob
         oobflag=1
         return
      else
         oobflag=0
      endif

      if(dindex .eq. 1) then !clockwise
        do i = 1, nth(nh)
          if(angle(i,nh).lt.angle(i+1,nh)) then  !normal
            if(angle(i,nh).le.az_deg .and. angle(i+1,nh).gt.az_deg) then
              nazimuth=i
              sepi1=angle(i+1,nh)-az_deg
              sepi=az_deg-angle(i,nh)
c              print *, 'normal, clockwise'
              return
            endif
          elseif((angle(i,nh).gt.angle(i+1,nh)).and.
     X           (angle(i+1,nh).ge.2.0)) then !Problem
c Ray is retrograding, and we'll just skip it. 04/13/01 MB
c              print *, 'Retrograde!(cw) ',angle(i,nh),' ',angle(i+1,nh)
c             oobflag=1
c             return
          elseif((az_deg.lt.180.).and.(angle(i+1,nh) .ge. 
     X           az_deg)) then  !crossing 0
              nazimuth=i
              sepi1=angle(i+1,nh)-az_deg
              sepi=(360.-angle(i,nh))+az_deg
c              print *, 'cross zero(cw) <', angle(i+1,nh),' ',az_deg
              return
          elseif(az_deg .gt. 180. .and. angle(i,nh).le.az_deg) then
              nazimuth=i
              sepi=az_deg-angle(i,nh)
              sepi1=360.-az_deg+angle(i+1,nh)
c              print *, 'cross zero(cw) >', angle(i,nh),' ',az_deg
              return
          elseif(i.eq.nth(nh)) then
             if(angle(i,nh).le.az_deg .and. angle(1,nh).gt.az_deg) then
                nazimuth=i
                sepi1=angle(1,nh)-az_deg
                sepi=az_deg-angle(i,nh)
                return
             endif
          endif
        enddo
      else             !counterclockwise
        do i = 1, nth(nh)
          if(angle(i,nh).gt.angle(i+1,nh)) then  !normal
            if(angle(i+1,nh).le.az_deg .and. angle(i,nh).gt.az_deg) then
              nazimuth=i
              sepi1=az_deg-angle(i+1,nh)
              sepi=angle(i,nh)-az_deg
c              print *, 'normal, counterclockwise'
              return
            endif
          elseif((angle(i,nh).lt.angle(i+1,nh)).and.
     X           (angle(i,nh).ge.1.0)) then !Problem
c Ray is retrograding, and we'll just skip it. 04/13/01 MB
c              print *, 'Retrograde!(ccw) ',angle(i,nh),' ',angle(i+1,nh)
c               oobflag=1
c               return
          elseif((az_deg.gt.180. .and.angle(i+1,nh) .lt. az_deg)) then
              nazimuth=i
              sepi1=az_deg-angle(i+1,nh)
              sepi=360.-az_deg+angle(i,nh)
c              print *, 'cross zero <', angle(i+1,nh),' ',az_deg
              return
          elseif(az_deg.lt.180. .and. angle(i,nh) .gt. az_deg) then 
              nazimuth=i
              sepi1=360.-angle(i+1,nh)+az_deg
              sepi=angle(i,nh)-az_deg
c              print *, 'cross zero >', angle(i,nh),' ',az_deg
              return
          elseif(i.eq.nth(nh)) then
             if(angle(i,nh).ge.az_deg .and. angle(1,nh).lt.az_deg) then
                nazimuth=i
                sepi1=angle(1,nh)-az_deg
                sepi=az_deg-angle(i,nh)
                return
             endif
          endif
        enddo
      endif
      print *, 'Something is wrong in Search.' 
      print *, az_deg,startang(nh),nh,dindex
      oobflag=1
      return
c      stop
      end
**************************************************
      subroutine ave(temp1,temp2,diff,totdiff,temp,dummy)
 
      implicit none
      real diff,totdiff,temp1,temp2,temp,dummy

      if(abs(temp1-dummy).ge.0.01) then
        if(abs(temp2-dummy).ge.0.01) then
          if(abs(totdiff).lt.0.00001) then
             print *,'ave'
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
**************************************************
      subroutine findz(el,elev,nh,nel)

      implicit none
      integer l,nel,nh
      real el,elev(nel),el_tol

      nh=0
      do l=1,nel
        if(el.le.elev(l)) then
           nh=l
           return
        endif
      enddo
      nh=nel+1
      return

c      do l=1,(nel-1)
c        el_tol=(elev(l+1)-elev(l))/2
c        if(abs(el-elev(l)).le.el_tol) then
c          nh=l
c          return
c        endif
c      enddo 
c      el_tol=(elev(nel)-elev(nel-1))/2
c      if (abs(el-elev(nel)).le.el_tol)then
c         nh=nel
c         return
c      endif
c      return 

       end    
**************************************************
      subroutine calc_time(az,nh,deltatime)

      implicit none
      include 'readcappi_para.h'
      integer az,nh,deltatime
      integer oldtime,newtime

      oldtime = hour(1,1)*3600+minute(1,1)*60+second(1,1)

      if (hour(az,nh).lt.hour(1,1)) then
c Change days
         newtime = (hour(az,nh)+25)*3600+minute(az,nh)*60+second(az,nh)
      else
         newtime = hour(az,nh)*3600+minute(az,nh)*60+second(az,nh)         
      endif
      deltatime = newtime-oldtime

      return
      end

**************************************************
       SUBROUTINE writecappi(outfile,dzcappi,vecappi,vrelcappi)
       implicit none
       include 'readcappi_para.h'
       CHARACTER outfile*100,time*12
       integer i,j,k,lev
       real dat(nx,ny),vedat(nx,ny),vreldat(nx,ny)

       print*,'output file name(CAPPI 15Km)=',outfile
c       OPEN(25,FILE=outfile,STATUS='unknown',FORM='unFORMATTED',
c     +      ERR=90)
       OPEN(15,FILE=outfile,STATUS='unknown',FORM='FORMATTED',
     +      ERR=90)

c       write(25)id
       print*,'time,ix,iy,lev ',time,ix,iy,lev

c       write(25)(int(zt(k)),k=1,iz)
       do k=1,iz
         if(k.eq.1) WRITE(15,116)ID
116      format(10i8)
         write(15,17)k
17       format('level',i2)
         print*,'write level= ',k
         do j=1,iy
           do i=1,ix
             dat(i,j)=dzcappi(i,j,k)
             if (dat(i,j).le.0.) then
                dat(i,j)=dummy
             endif
             call dbzZe(dat(i,j),2,dummy) !convert reflectivity to log scale
             vedat(i,j)=vecappi(i,j,k)
             vreldat(i,j)=vrelcappi(i,j,k)
             if(abs(dat(i,j)-dummy).le.0.01) dat(i,j)=special
C             if(abs(vedat(i,j)-dummy).le.0.01 .or. abs(vedat(i,j))
c     X           .le. 0.01) vedat(i,j)=special
C             if(abs(vreldat(i,j)-dummy).le.0.01 .or. abs(vreldat(i,j))
C     X           .le. 0.01) vreldat(i,j)=special
             if(abs(vedat(i,j)-dummy).le.0.01) vedat(i,j)=special
             if(abs(vreldat(i,j)-dummy).le.0.01) vreldat(i,j)=special
            enddo
c           write(*,*)i,j,k,(vedat(i,j),i=1,ix)
           write(15,18)j
18         format('azimuth',i3)
           write(15,19)id(176),id(177),id(178),id(179),
     X               (dat(i,j),i=1,ix)
           write(15,19)id(181),id(182),id(183),id(184),
     X               (vedat(i,j),i=1,ix)
           write(15,19)id(186),id(187),id(188),id(189),
     X               (vreldat(i,j),i=1,ix)
19         format(4a2,/,(8E10.3))
         enddo
c         write(25)((dat(i,j),i=1,ix),j=1,iy)
c Commented by wcl on 11/11/96
c       write(25)tmh
       enddo
       close(25)
 90    RETURN
       END



*********************************************************************
* DBZZE converts reflectivity factor to reflectivity or vise versa
*   in=1 log to linear scale, in=others linear to log scale
*********************************************************************       
        subroutine dbzZe(dz,in,dummy) 

        implicit none
        integer in
	real dz,dummy

            if (abs(dz-dummy).le.0.001)return
            if (in.eq.1)then
              dz=10.**(dz/10.) 
            else
	      dz=10.*alog10(dz)
	    endif
        return
	end
**************************************************
       subroutine grid(xc,yc,xd,yd,xdet,ydet,xmin,ymin)

       implicit none
       include 'readcappi_para.h' 
       integer i,j
       real pi,xdet,ydet

       pi=acos(-1.)
       do 10 i=1,nx
       do 10 j=1,ny
       xd(i,j)=xmin+xc+(i-1)*xdet+0.01
       yd(i,j)=ymin+yc+(j-1)*ydet+0.01
 10    continue
       return
       end
**************************************************
       subroutine writcaa(no,file,lev,nx,ny,dummy,irr)
       CHARACTER FILE*100
       INTEGER*2 IYER,IMON,IDAY,IHOR,IMIN,ISEC,KK
       INTEGER*2 IDEXZ,IDEXW,ISLPZ,ISLPW,ISRDZ,ISRDW
       INTEGER*2 LO1,LO2,LA1,LA2,IELEV(20),iev(20)
       integer lnblnk
       irr=1
       file=file(:lnblnk(file))
       open(10,file=file,status='old',form='unformatted',err=8)
       GO TO 9
8      PRINT*,' ',file,' NOT FOUND',file
       close(10)
       return
9      CONTINUE
       READ(10)IYER,IMON,IDAY,IHOR,IMIN,ISEC,LO1,LO2,LA1,LA2,
     *          ISLPZ,ISRDZ,IDEXZ,ISLPW,ISRDW,IDEXW,KK
       iy=iyer
       im=imon
       id=iday
       ih=ihor
       imn=imin
       is=isec
       lev=kk
       READ(10)IELEV
       do 10 k=1,lev
       iev(k)=ielev(k)*64./100.
 10    continue
       write(no)iy,im,id,ih,imn,is,nx,ny,dummy,lev
       write(no)(iev(k),k=1,lev)
       print*,'time,nx,lev ',iy,im,id,ih,imn,is,nx,ny,lev
       irr=0
       close(10)
       return
       end
**************************************************
       subroutine wrtdn(no,dz,tme,ht,nx,ny)
       real dz(nx,ny),ht(nx,ny),tme(nx,ny)
       write(no)dz
       write(no)ht
       write(no)tme
       return
       end
**************************************************
       subroutine caaread(file1,lev,elev,dz,ier,dummy)

       CHARACTER FILE1*100
       INTEGER*2 IYER,IMON,IDAY,IHOR,IMIN,ISEC,KK
       INTEGER*2 IDEXZ,IDEXW,ISLPZ,ISLPW,ISRDZ,ISRDW
       INTEGER*2 IDZ(120),LO1,LO2,LA1,LA2,IELEV(20)
       INTEGER*2 ITIME(420),IZ1(120,420),ITIME1(420)
       real dz(120,420),tim(120,420)
       integer*2 id(510),sf,cf
       integer lnblnk
       common /header1/ sf,cf,id
       common /header/ iyer,imon,iday,ihor,imin,isec
c** read one level CAA doppler radar dBZ(zn*.dat) or Vr(wn*.dat)
c** file(input)          : data file name
c** lev(input)           : data's level(its elevation angle)
c** elev(output)         : elevation angle(degrees)
c** tim(420)(output)     : azimuthal data's time(unit: sec)
c** dz(120,420)(output)  : dBZ data(no data =dummy)
c** tim(120,420)(output) : time  data(no data =dummy)
c** dummy(input)         : no data value
c** ier(output)          : ier=0 is normal , ier=1 is error
       ier=0
       GO TO 9
8      PRINT*,' ',file1,' NOT FOUND'
       ier=1
       return
9      CONTINUE
       file1=file1(:lnblnk(file1))
c       OPEN(10,FILE=FILE1,STATUS='OLD',FORM='UNFORMATTED',
c     +      READONLY,ERR=8)
       OPEN(10,FILE=FILE1,STATUS='OLD',FORM='UNFORMATTED',ERR=8)
       READ(10)IYER,IMON,IDAY,IHOR,IMIN,ISEC,LO1,LO2,LA1,LA2,
     *          ISLPZ,ISRDZ,IDEXZ,ISLPW,ISRDW,IDEXW,KK
       idate=iyer*10000+imon*100+iday
       ntim=ihor*10000+imin*100+isec
       if( (lev.lt.1).or.(lev.gt.kk) )then
	print*,'level is error ',lev,' max level= ',kk
	ier=1
	return
       endif
       RLO=LO1+LO2/10000.
       RLA=LA1+LA2/10000.
       SLOPZ=ISLPZ/1000.
       SORDZ=ISRDZ/100.
       SLOPW=ISLPW/1000.
       SORDW=ISRDW/100.
       call str2int('CA',id(13))
       call str2int('A ',id(14))
       call str2int('  ',id(15))
       call str2int('PP',id(16))
       call str2int('I ',id(17))
       call str2int('  ',id(18))
       id(21)=iyer
       id(22)=imon
       id(23)=iday
       id(24)=ihor
       id(25)=imin
       id(26)=isec
       id(33)=int(la1)
       id(34)=int(la2/10000.*60.)
       id(35)=int((la2-(int(la2/100.)*100))/10000.*3600.*sf)
       id(36)=int(lo1)
       id(37)=int(lo2/10000.*60.)
       id(38)=int((lo2-(int(lo2/100.)*100))/10000.*3600.*sf)
       id(67)=-999
       id(116)=iyer
       id(117)=imon
       id(118)=iday
       id(119)=ihor
       id(120)=imin
       id(121)=isec

       READ(10)IELEV
       elev=ielev(lev)/100.
       DO 666 K=1,lev-1
       DO 630 J=1,420
       READ(10)ITIME(J),IDZ
 630   CONTINUE
       IF(K.EQ.kk-1) THEN
        DO 40 J=1,420
        ITIME1(J)=ITIME(J)
 40     CONTINUE
       ENDIF
 666   continue
       DO 30 J=1,420
       READ(10)ITIME(J),IDZ
       DO 20 I=1,120
       IZ1(I,J)=IDZ(I)
       tim(i,j)=itime(j)/10.
 20    CONTINUE
 30    CONTINUE
       DO 50 J=1,420
       DO 50 I=1,120
       IF(lev.EQ.kk) THEN
        IF(ITIME1(J).EQ.ITIME(J)) THEN
	 dz(i,j)=dummy
         GO TO 50
        ENDIF
       ENDIF
        IF( (IZ1(I,J).EQ.0).or.(IZ1(i,j).eq.255) ) THEN
	 dz(i,j)=dummy
        else
         dz(I,J)=(IZ1(I,J)-IDEXZ)*SLOPZ+SORDZ
        endif
        if(dz(i,j).gt.70.) dz(i,j)=dummy 
 50    CONTINUE
       CLOSE(10)
       return
       END
**************************************************
       SUBROUTINE intecaa(dz,tim,elev,dummy,xd,yd,nx,ny,gdz,tme,
     *                    ht,rang,gama,dex)
       real dz(120,420),tim(120,420),xd(nx,ny),yd(nx,ny)
       real ht(nx,ny)
       real gdz(nx,ny),tme(nx,ny)
       real wkdz(1000),wkran(1000),wkht(1000),wktme(1000)
c** This subroutine is interpolation CAA dopller radar dBZ or Vr by 
c**   barns mehtods. 
c**   dz(120,420)(input)  : some one elevation angle PPI dBZ data
c**   elev(input)         : elevation angle (degree)
c**   dummy(input)        : no data value
c**   xd(nx,ny)(input)    : x-grid coordinate
c**   yd(nx,ny)(input)    : y-grid coordinate
c**   gdz(nx,ny)(output)  : output dBZ data
c**   ht(nx,ny)(output)   : the height of analysis data
c**   tme(nx,ny)(output)  : the tme of analysis data
c**   dex(input)          : grid interval(km)
c**   deg(input)          : the angle(degree) between new x-axis and W-E axis
c**   rang(input)         : inflence radius(km) of cressman or barns scheme
c**   gama(input)         : parameter of barns cheme,its range(0.1 to 1.0) 

       PI=ACOS(-1.)
       if(rang*rang .lt. 0.0001) stop
       wt=-4.*gama/(rang*rang)
       DET0=PI/210.
       do 110 i=1,nx
         do 100 j=1,ny
           rad=sqrt(xd(i,j)**2+yd(i,j)**2)
           if(abs(cos(pi*elev/180.)).lt.0.000001) stop
           ran=rad/cos(pi*elev/180.)
           IR=ran
           if( (ir.le.1).or.(ir.ge.119) ) then
             gdz(i,j)=dummy
  	     go to 100  !outside the radar range, get another point
           endif
           DER=ran-IR
c           IF(XD(i,j).GE.0.) THEN
c             ANG=ACOS(yd(i,j)/RAD)
c           ELSE
c             ANG=2.*PI-ACOS(yd(i,j)/RAD)
c           ENDIF
           ANG=atan2(yd(i,j),xd(i,j))
           if(abs(det0).lt.0.00001) stop
           JANG=ANG/DET0+1.
           J1=JANG
           J2=J1+1
           IF(J1.LT.1)THEN
             J1=420-J1
           ENDIF
           IF(J2.GT.420)J2=J2-420
           if(abs(det0).lt.0.00001) stop
           DET=ANG/DET0+1-J1
           ik1=ir-rang-1
           ik2=ir+rang+1
           if(ik2.ge.119)ik2=119
           if(ik1.le.5)ik1=5
           if(abs(det0*ran).lt.0.00001) stop
           jth=rang/(ran*det0)+1
           if(jth.gt.10)jth=10
           if(jth.lt.1)jth=1
           ij1=j1-jth
           if(ij1.lt.1)ij1=420+ij1
           ij2=ij1+2*jth+1
           nn=0
           do 21 i1=ik1,ik2,1
             do 20 jn=ij1,ij2,1
               jn1=jn
               if(jn1.gt.420)jn1=jn1-420
               if( int(dz(i1,jn1)).eq.dummy)go to 20
               az=(jn1-1)*det0
               xx=real(i1)*sin(az)*cos(pi*elev/180.)-xd(i,j)
               yy=real(i1)*cos(az)*cos(pi*elev/180.)-yd(i,j)
               if( (xx.gt.rang).or.(yy.gt.rang) )go to 20
               dist=sqrt(xx*xx+yy*yy)
               if(dist.gt.rang)go to 20
               nn=nn+1
               wkht(nn)=i1*sin(pi*elev/180.)+i1*i1/17000.
               wkdz(nn)=dz(i1,jn1)
               wktme(nn)=tim(i1,jn1)
               wkran(nn)=dist
 20          continue
 21        continue
           if(nn.lt.2) then
	     gdz(i,j)=dummy
             go to 100
           endif
           wtsum=0.
           dzsum=0.
           htsum=0.
           tmsum=0.
           do 30 kk=1,nn
             if(wkran(kk)**2 .lt. 0.00001) stop
             wtsum=wtsum+1./(wkran(kk)**2)
             dzsum=dzsum+wkdz(kk)/(wkran(kk)**2)
             htsum=htsum+wkht(kk)/(wkran(kk)**2)
             tmsum=tmsum+wktme(kk)/(wkran(kk)**2)
 30        continue
           if (abs(wtsum).lt.0.00001) stop
           outdz=dzsum/wtsum
           outht=htsum/wtsum
           outtm=tmsum/wtsum
           do 40 kk=1,nn
             wkdz(kk)=wkdz(kk)-outdz
             wkht(kk)=wkht(kk)-outht
             wktme(kk)=wktme(kk)-outtm
 40        continue
           wtsum=0.
           dzsum=0.
           vesum=0.
           htsum=0.
           tmsum=0.
           do 50 kk=1,nn
             disqt=wkran(kk)**2
             weigt=exp(disqt*wt)
             wtsum=wtsum+weigt
             dzsum=dzsum+weigt*wkdz(kk)
             htsum=htsum+weigt*wkht(kk)
             tmsum=tmsum+weigt*wktme(kk)
 50        continue
           if(wtsum.lt.0.0000001)then
	     gdz(i,j)=outdz
	     ht(i,j)=outht
	     tme(i,j)=outtm
	     go to 100
           endif
           gdz(i,j)=outdz+dzsum/wtsum
           ht(i,j)=outht+htsum/wtsum
           tme(i,j)=outtm+tmsum/wtsum
 100     continue
 110   continue
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
      integer*2 oldint,int1,int2,newint
      int1 = mod(oldint,256)
      int2 = (oldint - int1)/256
      newint = int1*256 + int2
      return
      end

