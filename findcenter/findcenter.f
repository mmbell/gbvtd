      PROGRAM findcenter
      implicit none
      include 'gbvtd_parameter.f'
      INTEGER 	IMAX,JMAX,KMAX,i,j
      real SX,SY,SZ,XZ,YZ,ZZ,ox(20),oy(20),Eps,x_len(20),y_len(20)
      real dz(nx,ny),ve(nx,ny),hgt,radii,value,ox_inc(20),oy_inc(20)
      real mwind(50),mean_wind_new
      Integer flag1,flag2,iz1,numgap,gap(5),kc,k,k1,k2,n
      integer ring_start,ring_end,It_max,rflag,iii,ii,ring
      Integer*2 id(510),ni,n_vol
      real vert_x, vert_y, vert_x_mean, vert_y_mean,a,b,ring_width
      complex vert_sum,vert_mean,vert_good,ai
      real gooddone,sd,mean_std,mean_std_new
      real vert_x_mean_new,vert_y_mean_new,vert_std_new,wind,mean_wind
      real ox1,oy1,radii_temp,vert_std,xctr(50),yctr(50),zone
      integer radi_beg, radi_end, radi_inc,num_x,num_y
      integer nargs,iargc
      integer*2 field1,field2
      logical P_Max        
      character*80 fname(20),infile,outfile_sum(20),outfile_ind(20)
c
      common /gbvtd_p/ flag1,flag2,gap,id,ring_start,ring_end,rflag,
     X                 zone,num_x, num_y
      common /CRITR/Eps,It_Max   
      common /MaxMin/P_Max      
      common /goodverts/ k1,k2,radi_beg,radi_end
      common /std/ sd

      imax=nx
      jmax=ny
      kmax=nz
c
C read parameter file
      nargs=iargc()
      if(nargs.gt.0) then
         call getarg(1,infile)
      else
         print *, 'Please enter input file name:'
         read *, infile
      endif
      OPEN(12,FILE=infile,STATUS='OLD')
      read(12,111)flag1
      read(12,111)flag2
      read(12,115)field1,field2
      read(12,111)K1        !lower z-level
      read(12,111)K2        !upper z-level
      read(12,111)numgap
      do i=1,numgap
        read(12,111)gap(i)
      enddo
      read(12,111)radi_beg
      read(12,111)radi_end
      read(12,113)zone
      read(12,111)rflag     !write GBVTD results or not
111   format(i3)
112   format(a20)
      read(12,111)P_Max
      read(12,113)radii     !radius of influence for circle of SIMPLEX
      read(12,113)Eps       !converging criteria (initial value 0.001)
      read(12,111)It_Max    !Max number of iterations
      read(12,111)n_vol
113   format(f8.2, f8.2)
 115  format(a2,/,a2)
      call intswap(field1,field1)
      call intswap(field2,field2)
      do ni=1,n_vol
        READ(12,'(a80)')FNAME(ni)
        read(12,113)ox(ni),oy(ni)
        read(12,113)ox_inc(ni), oy_inc(ni)
        read(12,113)x_len(ni),y_len(ni)
        read(12,'(a80)') outfile_ind(ni)
        read(12,'(a80)') outfile_sum(ni)
        OPEN(60,FILE=outfile_ind(ni),STATUS='unknown',FORM='FORMATTED')
        OPEN(70,FILE=outfile_sum(ni),STATUS='unknown',FORM='FORMATTED')

        OPEN(10,FILE=fname(ni),STATUS='UNKNOWN',FORM='FORMATTED')
c 
        DO ring=radi_beg, radi_end
           ring_start=ring
           ring_end=ring
          write(70,114) ring_start, zone
          write(60,114) ring_start, zone
  114     format(i3,'-',f8.2)
c        print *, 'Eps=',Eps,'  It_Max=',It_Max          
          DO K=k1,k2 
            write(1,'("Processing data from level #",i2/)') k
            CALL READDATA(DZ,VE,10,K,id,imax,jmax,2,field1,field2)
            num_x=id(162)
            num_y=id(167)

c altitude
            HGT=float(id(170)+(k-1)*id(173))/1000.
            write(60,61) HGT
61          format('Height= ', f8.2)

c 
c   xctr,yctr : position of first guess center location 
c               in 3-D Doppler arrays use olat,olon 
c   rad : influence radius for circle that circumscribes the simplex
c         typically use 4 km
c   data : array of wd/ws values form one level in the 3-D Doppler array
c   imax,jmax : data array dimensions 
c   sx,sy : resolution of the data array
c   radi (2) : inner and outer radial limits for computation of TANW
c 
            write(*,'(" * Finding center on level #",i2/)') k 
            write(*,'("* Finding center at",f4.1," km")') hgt
            vert_x=0.
            vert_y=0.
            iii=0
            ii=0
            vert_std=0.
            vert_std_new=0.
            vert_x_mean=0.
            vert_y_mean=0.
            vert_x_mean_new=0.
            vert_y_mean_new=0.
            wind=0.
            mean_wind=0.
            mean_wind_new=0.
            mean_std=0.
            mean_std_new=0.
            do i=1,int(x_len(ni)/ox_inc(ni))+1
              do j=1,int(y_len(ni)/oy_inc(ni))+1
                radii_temp=radii  !reset radius of influence to initial value
                ox1=ox(ni)+(i-1)*ox_inc(ni)
                oy1=oy(ni)+(j-1)*oy_inc(ni)
                print *, 'processing initial center guess at', ox1,oy1
                call SIMCON(ox1,oy1,radii_temp,dz,ve,imax,jmax,sx,sy,
     X                 value,wind)
                if (wind .le. 100.) then
                  iii=iii+1
                  xctr(iii)=ox1
                  yctr(iii)=oy1
                  mwind(iii)=wind
                  vert_x=vert_x+xctr(iii)
                  vert_y=vert_y+yctr(iii)
                  mean_wind=mean_wind+mwind(iii)
                endif
              enddo 
            enddo
            vert_x_mean=vert_x/float(iii)
            vert_y_mean=vert_y/float(iii)
            mean_wind=mean_wind/float(iii)
            do i=1,iii
              vert_std=vert_std+(xctr(i)-vert_x_mean)**2
     X               +(yctr(i)-vert_y_mean)**2
              mean_std=mean_std+(mwind(i)-mean_wind)**2
            enddo
            vert_std=sqrt(vert_std/float(iii))
            mean_std=sqrt(mean_std/float(iii))
            do i=1,iii
              if(sqrt((xctr(i)-vert_x_mean)**2
     X              +(yctr(i)-vert_y_mean)**2).lt. vert_std) then
c     X              .and.abs(mwind(i)-mean_wind).lt.mean_std) then
                ii=ii+1
                xctr(ii)=xctr(i)
                yctr(ii)=yctr(i)
                mwind(ii)=mwind(i)
                vert_x_mean_new=vert_x_mean_new+xctr(ii)
                vert_y_mean_new=vert_y_mean_new+yctr(ii)
                mean_wind_new=mean_wind_new+mwind(ii)
              endif
            enddo
            if(ii.gt.0) then
               vert_x_mean_new=vert_x_mean_new/float(ii)
               vert_y_mean_new=vert_y_mean_new/float(ii)
               mean_wind_new=mean_wind_new/float(ii)
               do i=1,ii
                  vert_std_new=vert_std_new+(xctr(i)-vert_x_mean_new)**2
     X                 +(yctr(i)-vert_y_mean_new)**2
                  mean_std_new=mean_std_new+(mwind(i)-mean_wind_new)**2
               enddo
               vert_std_new=sqrt(vert_std_new/float(ii))
               mean_std_new=sqrt(mean_std_new/float(ii))
            else
               vert_x_mean_new=-999.
               vert_y_mean_new=-999.
               mean_wind_new=0.
               vert_std_new=99.
               mean_std_new=99.
            endif
            write(60,161) vert_x_mean, vert_y_mean, vert_std,iii,
     X                 vert_x_mean_new, vert_y_mean_new,vert_std_new, ii
161         format(' Mean X:',f8.3, ' Mean Y:', f8.3, ' STD:', f8.3,
     X      ' #=',i5,/,' Mean X New:', f8.3, ' Mean Y New:', f8.3, 
     X       ' STD_NEW:' ,f8.3, ' #=',i5)
            write(60,162) mean_wind,mean_std,mean_wind_new,mean_std_new
 162        format('Mean VT:',f8.3,16x,' STD:', f8.3,/,
     X           'Mean VT New:', f8.3,20x,' STD_NEW:' ,f8.3)
            write(*,161) vert_x_mean_new, vert_y_mean_new, vert_std,iii,
     X               vert_x_mean_new, vert_y_mean_new,vert_std_new, ii  
            WRITE (1,'("Finished level:",I2)') k
            write(70,163)k,vert_x_mean_new,vert_y_mean_new,vert_std_new,
     X                  ii, iii,mean_wind_new
 163        format(2x, i2, 5x, f9.3, 4x, f9.3, 4x, f6.3, 3x,i2,'/',i2
     X           ,4x,f9.2)
            
c Add good array to check final centers 03/02/01 MB
c            call goodvert(vert_x_mean_new,vert_y_mean_new,
c     X                    dz,ve,ring,k)
               
          ENDDO
        ENDDO
c Add new loops to process stddev
c        do ring=ring_start,ring_end
c           do k=k1,k2
c              gooddone=-999.
c              call goodvert(gooddone,gooddone,ring,k)
c              print *,'good:',vert_good
c           enddo
c        enddo
        close(10)
        close(60)
        close(70)
      enddo
      stop       
c 
c *********************** 
C Error Messages: 
c
998   print *, 'cannot open ideal input file'
   
      stop
      END 

      subroutine goodvert(vx,vy,dz,ve,ring,k)

      implicit none
      include 'gbvtd_parameter.f'
      real vx,vy,rleg,dz(nx,ny),ve(nx,ny),origsd,percent
      real sd,goodsd,xrsd,xlsd,yusd,ydsd,goodwind
      complex*8 good,good_xr,good_xl,good_yu,good_yd,best
      integer ring,k,k1,k2,radi_beg,radi_end,sdflag,first
      integer count
      common /goodverts/ k1,k2,radi_beg,radi_end
      common /std/ sd
      real x_good(radi_beg:radi_end,k1:k2)
      real y_good(radi_beg:radi_end,k1:k2)     
      complex ai

c Add good array to check final centers 03/02/01 MB
      if((vx.ne.0.0).or.(vy.ne.0.0)) then
         x_good(ring,k)=vx
         y_good(ring,k)=vy
         sdflag=1
         ai=(0,1.0)
         best=(-999.,-999.)
         first=0
         count=0
         do while ((best.eq.(-999.,-999.)).and.(count.le.100))
            good=x_good(ring,k)+ai*y_good(ring,k)
            good_xr=(x_good(ring,k)+.1)+ai*y_good(ring,k)
            good_xl=(x_good(ring,k)-.1)+ai*y_good(ring,k)
            good_yu=x_good(ring,k)+ai*(y_good(ring,k)+.1)
            good_yd=x_good(ring,k)+ai*(y_good(ring,k)-.1)
            call gbvtd(good,goodwind,dz,ve,sdflag)   
            goodsd=sd
            if (first.eq.0) then
               origsd=goodsd
               first=1
            endif
            call gbvtd(good_xr,rleg,dz,ve,sdflag)   
            xrsd=sd
            call gbvtd(good_xl,rleg,dz,ve,sdflag)   
            xlsd=sd
            call gbvtd(good_yu,rleg,dz,ve,sdflag)   
            yusd=sd
            call gbvtd(good_yd,rleg,dz,ve,sdflag)   
            ydsd=sd
            print *,'good=',goodsd,' r=',xrsd,' l=',xlsd
            print *,'u=',yusd,' d=',ydsd

            if ((xrsd.lt.goodsd).and.(xrsd.lt.xlsd).and.
     X         (xrsd.lt.yusd).and.(xrsd.lt.ydsd)) then
               x_good(ring,k)=x_good(ring,k)+.1
            elseif ((xlsd.lt.goodsd).and.(xlsd.lt.xrsd).and.
     X             (xlsd.lt.yusd).and.(xlsd.lt.ydsd)) then
               x_good(ring,k)=x_good(ring,k)-.1
            elseif ((yusd.lt.goodsd).and.(yusd.lt.xrsd).and.
     X             (yusd.lt.xlsd).and.(yusd.lt.ydsd)) then
               y_good(ring,k)=y_good(ring,k)+.1
            elseif ((ydsd.lt.goodsd).and.(ydsd.lt.xrsd).and.
     X             (ydsd.lt.xlsd).and.(ydsd.lt.yusd)) then
               y_good(ring,k)=y_good(ring,k)-.1
            else
               if(((x_good(ring,k)-vx)**2
     X              +(y_good(ring,k)-vy)**2).lt.1.) then
                  best=good
                  print *,'best center= ',best
                  percent=(origsd-goodsd)/origsd*100
                  vx = REAL(best)
                  vy = AIMAG(best)
                  write(70,164)vx,vy,goodsd,percent,goodwind
               else
                  best=vx+ai*vy
                  call gbvtd(best,goodwind,dz,ve,sdflag)
                  percent=0
                  write(70,164)vx,vy,goodsd,percent,goodwind
               endif
            endif
            count=count+1
         enddo
         if(count.ge.100) then
            best=vx+ai*vy
            call gbvtd(best,goodwind,dz,ve,sdflag)
            percent=0
            write(70,164)vx,vy,goodsd,percent,goodwind
         endif
      endif
 164  format ('Adjusted-',f9.3,4x,f9.3,4x,f6.3,3x,f5.2,'%',3x,f9.2)
      end



*****************************************************************
*   READDATA
*****************************************************************

       subroutine readdata(dz,ve,unit,ik,id,ix,iy,flag,
     X     field1,field2)
       implicit none
       include 'gbvtd_parameter.f'
       integer ik,i,j,ix,iy,j1,N1,k,n
       real special
       integer flag,unit,index,indexdz,indexve
       real dz(nx,ny),ve(nx,ny),verel(nx,ny)
       integer*2 id(510)
       integer*2 field1,field2
c       data field1,field2 /'MD','EV'/

       special=-999.
       read(unit,116)id
116    format(10i8)
        ix=id(162)
        iy=id(167)

       do j=1,id(175)   !id(175):total number of fields in the ascii file
          index=176+(j-1)*5
          if(id(index).eq.field1) indexdz=j
          if(id(index).eq.field2) indexve=j
       enddo

       do k=1,ik
         read(unit,117)N1
117      format(5x,i2)
         if(n1.ne.k) then
            print *, 'Height index not match', N1,k
            stop
         endif
         do j=1,iy
           read(unit,18)j1
c           print *, j1,j
           if(j1.ne.j) then
             print *, j1,j,' second index does not match'
             stop
           endif
18         format(7x,i3)
           do n=1,id(175)
           if(n.eq.indexdz) then
             read(unit,19)id(176),id(177),id(178),id(179),
     X           (dz(i,j),i=1,ix)
c           write(*,19)id(176),id(177),id(178),id(179),(dz(i,j),i=1,ix)
           elseif(n.eq.indexve) then
             read(unit,19)id(181),id(182),id(183),id(184),
     X           (ve(i,j),i=1,ix)
c           write(*,19)id(181),id(182),id(183),id(184),(ve(i,j),i=1,ix)
           else
             read(unit,19)id(186),id(187),id(188),id(189),
     X           (verel(i,j),i=1,ix)
           endif
           enddo
         enddo
       enddo
19       format(4a2,/,(8E10.3))
       return
       end

c **********************
      subroutine SIMCON (xctr,yctr,rad,dz,ve,imax,jmax,sx,sy,value,
     X  wind)
      real dz(imax,jmax), ve(imax,jmax),radi(2),rad_temp
      complex*8 Vert(3),Rvert,Xvert,Cvert,Svert(3) 
      complex*8 temp,temp1
      real LegAr(3),Rleg,Xleg,Cleg,RCoef,XCoef,CCoef,SCoef,wind
      logical Done,DOIT                                       !Mod 04SEP86
      integer iwrit,i       
      character*3 CMxn
      integer sdflag
                                                          
      data RCoef,XCoef,CCoef,SCoef/-1.0,-2.0,0.5,0.5/, iwrit/1/ 
c 
c   xctr,yctr : position of first guess center location 
c               in 3-D Doppler arrays use olat,olon 
c   rad : influence radius for circle that circumscribes the simplex
c         typically use 4 km
c   DZ,VE : array of wd/ws values form one level in the 3-D Doppler array
c   imax,jmax : data array dimensions 
c   sx,sy : resolution of the data array
c   radi (2) : inner and outer radial limits for computation of TANW
c 
      sdflag=0
      xctr_ini=xctr
      yctr_ini=yctr
      call SETUP(xctr,yctr,rad,dz,ve,imax,jmax,sx,sy,radi, 
     #        Vert,LegAr,Done,I_Bst,I_Wst,Kount,CMxn) 
c 
      do while (.NOT. Done)                                                     
        Kount = Kount + 1                                                       
        print *, 'Kcount, VERT,LegAr', kount,vert,LegAr 
c       write (iwrit,10)Kount 
        call SIZE(Vert,Rvert,RCoef,I_Bst,I_Wst) 
c        call ALSAND(data,imax,jmax,sx,sy,radi,Rvert,Rleg)
        call gbvtd(Rvert,Rleg,dz,ve,sdflag)
        print *, 'Rvert,Rleg', rvert,rleg
        if (DOIT(Rleg,LegAr(I_Bst))) then                     !Mod 04SEP86      
          call SIZE(Vert,Xvert,XCoef,I_Bst,I_Wst)       
c          call ALSAND(data,imax,jmax,sx,sy,radi,Xvert,Xleg) 
          call gbvtd(Xvert,Xleg,dz,ve,sdflag) 
          print *,'Xvert,Xleg',xvert,xleg
          if (DOIT(Xleg,LegAr(I_Bst))) then                   !Mod 04SEP86      
            call STASH(Vert,LegAr,Xvert,Xleg,I_Wst)     
          else                  
            call STASH(Vert,LegAr,Rvert,Rleg,I_Wst)     
          endif                 
        else
          if (DOIT(Rleg,LegAr(I_Wst))) then                   !Mod 04SEP86      
            call STASH(Vert,LegAr,Rvert,Rleg,I_Wst)     
          else                  
            call SIZE(Vert,Cvert,CCoef,I_Bst,I_Wst)     
c            call ALSAND(data,imax,jmax,sx,sy,radi,Cvert,Cleg)       
            call gbvtd(Cvert,Cleg,dz,ve,sdflag)
            print *, 'Cvert,Cleg', cvert,cleg 
            if(DOIT(Cleg,LegAr(I_Wst)).AND.DOIT(Cleg,Rleg))then  !Mod 04SEP86   
              call STASH(Vert,LegAr,Cvert,Cleg,I_Wst)   
            else                
              call SHRINK(Vert,Svert,SCoef,I_Bst,I_Wst) 
c              call ALSAND(data,imax,jmax,sx,sy,radi,Svert(2),LegAr(2))
              call gbvtd(Svert(2),LegAr(2),dz,ve,sdflag)
              print *,'Svert(2),LegAr(2)',svert(2),legar(2) 
              Vert(2) = Svert(2)
c              call ALSAND(data,imax,jmax,sx,sy,radi,Svert(I_Wst), 
c     @                     LegAr(I_Wst))
              call gbvtd(Svert(I_Wst),LegAr(I_Wst),dz,ve,sdflag)
          print *,'Svert(I_Wst),LegAr(I_Wst)',Svert(I_Wst),LegAr(I_Wst) 
              Vert(I_Wst) = Svert(I_Wst)    
            endif               
          endif                 
        endif                   
        call RNKVRT(Vert,LegAr) 
        call CONTST(Kount,Done,Vert)        
c       write(iwrit,11)Kount,LegAr(I_Bst),Vert(I_Bst) 
        value=LegAr(I_bst)
        xctr=real(Vert(I_bst))
        yctr=aimag(Vert(I_bst))
        print *,'xctr,yctr',xctr,yctr 
        temp1=vert(1)-vert(2)
        rad_temp=cabs(temp1)
        if(rad_temp.le.rad) rad=rad_temp
        temp1=vert(1)-vert(3)
        rad_temp=cabs(temp1)
        if(rad_temp.le.rad) rad=rad_temp
        temp1=vert(2)-vert(3)
        rad_temp=cabs(temp1)
        if(rad_temp.le.rad) rad=rad_temp 
      enddo 
      wind=LegAr(I_Bst)
      write(iwrit,10)Kount
      write(iwrit,13) CMxn,LegAr(I_Bst),Vert(I_Bst) 
      write(*,10)Kount
      write(*,13) CMxn,LegAr(I_Bst),Vert(I_Bst) 
10    format(' - - Iteration #',i4)         
c11   format(' Iter:',i4,';   Max. Value: ',g12.5,' at ', 2F10.5)       
13    format('     Reporting ',a3,' of ',g12.5,' at ',2f10.5//)
      write(60,61) xctr_ini,yctr_ini,Kount,Vert(I_Bst),rad, 
     X     LegAr(I_Bst)
61    format(2f8.2,i5,2f8.3,3x,f8.3,f10.5)  
      return
      end 
      
      subroutine SETUP(xctr,yctr,rad,dz,ve,imax,jmax,sx,sy,radi, 
     #        Vert,LegAr,Done,I_Bst,I_Wst,Kount,CMxn)   
      real dz(imax,jmax), ve(imax,jmax), radi(2)
c      common/CRITR/Eps,It_Max   
      common /MaxMin/P_Max      
      complex*8 Vert(3),ai
      character*3 CMxn          
      real LegAr(3),xctr,yctr,rad,sqr32 
      logical Done,P_Max        
      data sqr32/0.866025/           !squareroot of 3 divided by 2
c 
      Done = .FALSE.            
c      P_Max = .TRUE.            
c      P_Max = .FALSE.            
      if (P_Max)then            
              I_Bst = 3         
              I_Wst = 1         
              CMxn = 'Max'      
      else
              I_Bst = 1         
              I_Wst = 3         
              CMxn = 'Min'      
      endif 
      Kount = 0                 
c 
c  initial vertices are at the points of an equilateral triangle
c   circumscribed in a circle with radius "rad" 
c 
      ai=(0,1.0)
      Vert(1) = xctr + ai*(yctr+rad)
      Vert(2) = xctr + sqr32*rad + ai*(yctr-.5*rad) 
      Vert(3) = xctr - sqr32*rad + ai*(yctr-.5*rad) 
c      Eps = 0.001               
c      It_Max = 60               
      do i = 1,3                
c        call ALSAND(data,imax,jmax,sx,sy,radi,Vert(i),LegAr(i))
        call gbvtd(Vert(i),LegAr(i),dz,ve,sdflag) 
      enddo 
      call RNKVRT(Vert,LegAr)   
      return
      end 
      
      subroutine RNKVRT(Vert,LegAr)         
      complex*8 Vert(3),Tvert   
      real LegAr(3),Tleg      
      do i = 1,3                
        do j = 1,2              
          if (LegAr(j) .GT. LegAr(j+1)) then
            Tvert = Vert(j)     
            Tleg = LegAr(j)     
            Vert(j) = Vert(j+1) 
            LegAr(j) = LegAr(j+1) 
            Vert(j+1) = Tvert   
            LegAr(j+1) = Tleg   
          endif                 
        enddo                   
      enddo 
c     do i = 1,3                
c       write(1,10)i,Vert(i),LegAr(i) 
c     enddo 
c10   format(' RNKVRT: Vertex# ',i2,' Vert = ',2f10.5,' LegAr = 'g12.5)         
      return
      end 
      
      subroutine STASH(Vert,LegAr,Nvert,Nleg,I_Wst)     
      complex*8 Vert(3),Nvert   
      real LegAr(3),Nleg      
c 
c     write(iwrit,9)Vert(I_Wst),LegAr(I_Wst)
      Vert(I_Wst) = Nvert       
      LegAr(I_Wst) = Nleg       
c     write(iwrit,10)Vert(I_Wst),LegAr(I_Wst) 
c9    format(' STASH: Old Vert(I_Wst) and LegAr(I_Wst) ',2f10.5,g12.5)          
c10   format(' STASH: New Vert(I_Wst) and LegAr(I_Wst) ',2f10.5,g12.5)          
      return
      end 
      
      subroutine SIZE(Vert,Nvert,Coef,I_Bst,I_Wst)      
      complex*8 Vert(3),Nvert,Bvert         
      Bvert = 0.5 * (Vert(I_Bst) + Vert(2)) 
      Nvert = Bvert + Coef*(Vert(I_Wst)-Bvert)          
c     write(iwrit,10)Coef,Vert(I_Wst),Nvert 
c10   format(' SIZE: Coef =',f10.5,' Old =',2f10.5,' New =',2f10.5) 
      return
      end 
      
      subroutine SHRINK(Vert,Svert,SCoef,I_Bst,I_Wst)   
      complex*8 Vert(3),Svert(3)
      Svert(2)=Vert(I_Bst)+SCoef*(Vert(2)-Vert(I_Bst))  
c     write(iwrit,10)2,Vert(2),Svert(2) 
      Svert(I_Wst)=Vert(I_Bst)+SCoef*(Vert(I_Wst)-Vert(I_Bst))      
c     write(iwrit,10)I_Wst,Vert(I_Wst),Svert(I_Wst) 
c10   format(' SHRINK: Vertex#',i2,' Old =',2f10.5,' New =',2f10.5) 
      return
      end 
      
      subroutine CONTST(Kount,Done,Vert)    
      complex*8 Vert(3)         
      Logical Done              
      common/CRITR/Eps,It_Max   
      SMax = 0.0                
      do i = 1,3                
        j = i + 1               
        if (j .GT. 3) j = j - 3 
        SS = CABS(Vert(i)-Vert(j))          
        SMax = amax1(SS,SMax)   
c       write(iwrit,10)i,j,SS 
      enddo 
      Done = (SMax .LT. Eps).OR.(Kount.GE.It_Max)       
c10   format('CONTST: Side between ',i2,' and ',i2,' has length ',F10.5)        
      return
      end 
      
      logical function DOIT(ALn,ALo)        
      common /MaxMin/P_Max      
      logical P_Max             
      if(P_Max) then            
              DOIT = ALn .GT. ALo 
      else
              DOIT = ALn .LT. ALo 
      endif 
      return
      end 


********************************************************************
*  Subroutine GBVTD
********************************************************************

      subroutine gbvtd(Vert,Rleg,dz,vd,sdflag)

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
      integer K,nn,i,numcoeff,n
      real pi,xc,yc,dummy,vrcoef,vm,s,theta_t,rad_lat,rad_long
      REAL THETAZ(thetamax),d
      REAL DZ(nx,ny),DDZ(thetamax),XDD(thetamax),A(maxrays)
      real x(thetamax),sd
      real YDD(thetamax),THETAV(maxrays),dddz(rmax,thetamax)
      real xdddz(thetamax), ydddz(thetamax),VT(rmax,thetamax)
      REAL XS(maxcoeffs,maxrays),YS(maxrays),WEIGHT(maxrays)
      real CC(maxcoeffs),STAND_ERR(maxcoeffs),theta,dist
      REAL TH(maxrays),AZ(maxrays),VTCOEF(maxcoeffs),xdist,ydist
      real VD(nx,ny),DVD(maxrays),DDVD(rmax,thetamax),dvd1(thetamax)
      real temp,ddz1(thetamax)
      real ampvt1,angle,zone
       real phasevt1,dtor,ampvt2,phasevt2,ampvt3,phasevt3,vorticity
       real ampvr1,phasevr1,ampvr2,phasevr2,ampvr3,phasevr3
       real meanflow,radi(2),vrmean,xctemp,yctemp
       real meanangle,r_t,ar,rleg
       character inputfile*80,outfiledz,outfileve
      real vt0max,vt1max,vr0max
      INTEGER NDATA,gap(5),flag1,j,flag2,rflag,num_ray
      integer iori,jori,ring_start, ring_end,iii,num_x,num_y
      integer sdflag
      integer*2 id(510),sf,af,id1(510)
      complex*8 Vert         
      LOGICAL FLAG,P_Max
      data rad_lat, rad_long /25.0694,121.2083/
      data sf, af /100,64/
       common /ideal/vorticity,ar,meanflow,meanangle,phasevt1,ampvt1,
     X               phasevt2,ampvt2,phasevt3,ampvt3,phasevr1,ampvr1,
     X               phasevr2,ampvr2,phasevr3,ampvr3
      common /gbvtd_p/ flag1,flag2,gap,id,ring_start,ring_end,rflag,
     X                 zone,num_x,num_y
      common /MaxMin/ P_Max
      common /std/ sd
c      CALL OPNGKS
      print *, maxrays
      PI=ACOS(-1.)
      DUMMY=-999.
      dtor=pi/180.
      vt0max=.0
      vt1max=.0
      vr0max=.0
c          CALL READDATA(DZ,VD,10,N,id,nx,ny,2)
      iori=id(309)/sf
      jori=id(310)/sf
      iii=0
      if(flag2.eq.2 .and. (flag1.eq.3 .or. flag1.eq.4)) then
c             print *, 'please enter ideal input filename:'
c             read *, inputfile
             inputfile='/scr/science9/wlee/gbvtd/ideal/rankine.sym.inp'
             OPEN(1200,file=inputfile,status='old',form='formatted',
     +               err=998)
             READ(1200,'(A80/A80/f7.2,1x,f7.2)') outfiledz,outfileve,
     +          xctemp,yctemp
             read(1200,'(f8.1,/,f5.1)') vorticity,ar
             read(1200,'(f8.1,/,f8.1)')meanflow,meanangle
             read(1200,'(f8.1,/,f8.1)')phasevt1,ampvt1
             read(1200,'(f8.1,/,f8.1)')phasevt2,ampvt2
             read(1200,'(f8.1,/,f8.1)')phasevt3,ampvt3
             read(1200,'(f8.1)')vrmean
             read(1200,'(f8.1,/,f8.1)')phasevr1,ampvr1
             read(1200,'(f8.1,/,f8.1)')phasevr2,ampvr2
             read(1200,'(f8.1,/,f8.1)')phasevr3,ampvr3
             close(1200)
c             print *,xctemp,yctemp,vorticity,ar
      endif
             OPEN(13,file='result_coeff',status='unknown',
     +            form='formatted')

      xc = REAL(Vert)           
      yc = AIMAG(Vert)
      if(rflag.eq.1) then
        print *, 'Process GBVTD on center:', Vert
      endif
      
cc      do k=1,ring_end
cc        do j=1,thetamax
cc          dddz(k,j)=dummy
cc          ddvd(k,j)=dummy
cc          vt(k,j)=dummy
cc        enddo
cc      enddo
 
C radius of GBVTD rings from ring_start to ring_end km every 1 km
      do K=ring_start,ring_end   
C ******  GRID POINT    ******
C ******  Calculate grid point location for VDAD or GBVTD
C ******  the grid point is in corresponding si coordinate for GBVTD
C ******  theta is angle array (si for GBVTD and theta for VDAD)
cc         CALL GRID(XC,YC,XDD,YDD,FLOAT(K),THETAV,theta_t,flag1,
cc     X                    thetamax)
cc        call grid(xc,yc,xdddz,ydddz,float(k),thetaZ,theta_t,1,
cc     X                    thetamax)
C ******  GET GRID POINT DATA FROM Doppler Velocity Field ******
cc         CALL INTE(DZ,DUMMY,XDDdz,YDDdz,thetamax,DDZ,2.0,0.5,
cc     X                iori,jori,id(162),id(167))
cc         CALL INTE(VD,DUMMY,XDD,YDD,thetamax,DVD,2.0,0.5,iori,
cc     X                jori,id(162),id(167))
cc         if(flag1.eq.1 .or. flag1.eq.3) then !vdad
cc           r_t=sqrt(xc*xc+yc*yc)
ccc                print *, 'int'
cc           do i=1,thetamax
ccc                  print *, dvd(i),xdd(i),ydd(i),iori,jori,d,r_t
cc             d=sqrt((xdd(i))**2.+(ydd(i))**2.)
cc             dvd(i) = dvd(i)*d/r_t
cc           enddo
cc         endif
cc         if(flag2.eq.2) then
cc           if(flag1.eq.3 .or. flag1.eq.4) then  !use analytic wind
ccc             call analytic(VD,dummy,xdd,ydd,thetamax,dvd1,2.0,0.5,
ccc     X                iori,jori,xc,yc,k,ddz1,flag1)
cc             call analytic(VD,dummy,xdd,ydd,thetamax,dvd1,2.0,0.5,
cc     X                iori,jori,xctemp,yctemp,k,ddz1,flag1)
cc             do i=1,thetamax
c c              DVD(i)=dvd1(i)
cc               DDZ(i)=DDZ1(i)
cc            enddo
cc           endif
cc         endif
            num_ray=0
            do i=1,num_x
              do j=1,num_y
                xdist=i-(xc+iori)
                ydist=j-(yc+jori)
                dist=sqrt(xdist**2+ydist**2)
                if(dist.ge.(k-zone).and.dist.le.(k+zone))then
                  angle=atan2(ydist, xdist)
                  call fixangle(angle)
c                 angle=amod((2.5*pi-angle),2.*pi)-theta_t
                  angle=angle-theta_t
                  call fixangle(angle)
                  if(flag1.eq.2 .or. flag1.eq.4) then  !gbvtd
                    call theta2si(k,angle,xc,yc,theta,theta_t)
c                 print *, i,j,angle/dtor,theta_t/dtor,theta/dtor,
c     X                 vd(i,j)
                  else
                    theta=angle
                  endif
                  num_ray=num_ray+1
                  thetav(num_ray)=theta
                  DVD(num_ray)=VD(i,j)
c                  print *, i,j,k,num_ray,angle/dtor,
c     X                 thetav(num_ray)/dtor,dvd(num_ray), dist
                endif
              enddo
            enddo
         call goodcirc(thetaV,DVD,num_ray,GAP,FLAG,numcoeff)
c Modified 1/11/01 to test result of flag
         if (flag) then
C         GET GOOD DATA NUMBER
C         DVD is the original array, TH is the cleaned up angle array
         CALL GOODDATA(DVD,THETAV,AZ,TH,NN,num_ray,dummy)
C ******  LEASTSQURRE CURVE FIT AND FIND COE. OF FOURIER sqrt.o  ******
         CALL CALCULATE_XY(NN,numcoeff,TH,AZ,XS,YS,NDATA,WEIGHT)
         CALL LLS(numcoeff,NDATA,NDATA,XS,YS,WEIGHT,S,STAND_ERR,
     X                 CC,FLAG)
         CALL FINAL1(CC,VTCOEF,VRCOEF,VM,XC,YC,K,numcoeff,flag1,rflag)

C Added std_dev call to test final center, 03/02/01 MB
C X is the theta array converted from si (has NN point in deg)
           if (sdflag.eq.1) then
            call std_dvi(A,CC,X,AZ,TH,numcoeff,sd,xc,yc,
     X                     flag1,k,theta_t,dummy)
           endif   

         else
            xc=-999.0
            yc=-999.0
            do n=1,8
               vtcoef(n)=0.0
            enddo
            vrcoef=0.0
            VM=0.0
         endif

c-- modified 9/30/98 to use P_Max to determine the minimization/maximization
C--  criteria. if P_Max is .true., then, Rleg is the mean tangential wind.
C--  If P_Max is .false., then, Rleg is the amplitude of the wavenumber 1.
         write(13, 1202) xc,yc,(vtcoef(n),n=1,8),vrcoef,VM
1202     format(3f8.3, 9f6.1)
c         if(P_Max.eq. .true.) then  !maximize mean tangential wind 
c           if(vtcoef(1).gt.0.) then
c             vt0max=vt0max+vtcoef(1)
c             iii=iii+1
cc            print *,k,vtcoef(1),ddz(1)
c           endif
c         elseif(P_Max.eq. .false.) then !minimize wavenumber 1 amplitude
c           temp=sqrt(vtcoef(3)*vtcoef(3)+vtcoef(4)*vtcoef(4))
c           if(temp.gt.0.0001) then
c             vt1max=vt1max+temp
c             iii=iii+1
c           endif 
c         endif
      ENDDO
c     if(iii.gt.0) then
        if(P_Max.eqv. .true.) then
          Rleg=vtcoef(1)
        else
          Rleg=sqrt(vtcoef(3)*vtcoef(3)+vtcoef(4)*vtcoef(4))
        endif
c      else
c        Rleg=0.
c      endif
999   CONTINUE
      rewind(10)
      return
998       print *, 'cannot open ideal input file'
      END

      
************************************************************************
* OPENFILE
************************************************************************

      Subroutine openfile(l,numloc,kc,prefix)

      integer l, numloc,kc,iz,lnblnk
      CHARACTER FNAME*80,STR*2,XFNAME*30
      character prefix*20
            if(iz.lt.10) then
              str='0'
              write(str,'(i1)')IZ
            else
              write(str,'(i2)')IZ
            endif
            xfname=prefix(:len_trim(prefix))//'f'//STR
            if(L.lt.10) then
              str='0'
              write(str,'(i1)')L
            else
              write(str,'(i2)')L
            endif
            xfname=xfname(:len_trim(xfname))//'c'
            xfname=xfname(:len_trim(xfname))//STR
            OPEN(15,FILE=xfname,STATUS='UNKNOWN',FORM='FORMATTED')
            if(iz.lt.10) then
              str='0'
              write(str,'(i1)')IZ
            else
              write(str,'(i2)')IZ
            endif
            xfname=prefix(:len_trim(prefix))//'f'//STR
            if(L.lt.10) then
              str='0'
              write(str,'(i1)')L
            else
              write(str,'(i2)')L
            endif
            xfname=xfname(:len_trim(xfname))//'z'
            xfname=xfname(:len_trim(xfname))//STR
            OPEN(16,FILE=xfname,STATUS='UNKNOWN',FORM='FORMATTED')
            if(iz.lt.10) then
              str='0'
              write(str,'(i1)')IZ
            else
              write(str,'(i2)')IZ
            endif
            xfname=prefix(:len_trim(prefix))//'f'//STR
            if(L.lt.10) then
              str='0'
              write(str,'(i1)')L
            else
              write(str,'(i2)')L
            endif
            xfname=prefix(:len_trim(prefix))//STR
            xfname=xfname(:len_trim(xfname))//'.out'
            OPEN(14,FILE=xfname,STATUS='UNKNOWN',FORM='FORMATTED')
            xfname=xfname(:len_trim(xfname)-3)//'point1'
            open(171,file=xfname,status='unknown',form='formatted')
            xfname=xfname(:len_trim(xfname)-1)//'2'
            open(172,file=xfname,status='unknown',form='formatted')
            xfname=xfname(:len_trim(xfname)-1)//'3'
            open(173,file=xfname,status='unknown',form='formatted')
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

**************************************************************
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

      subroutine analytic(dz,dummy,xd,yd,nx,gdz,rang,gama,iori,jori,
     X                 xc,yc,k,ddz,flag)

       implicit none
       real rang, gama,dummy,pi
       integer nx, i
       integer iori,jori,k,flag
       real dz(241,241),xd(nx),yd(nx),rt,xg,yg,theta
       real gdz(nx),radius,ddz(nx),xc,yc,asymvt,ampvt1
       real phasevt1,dtor,ampvt2,phasevt2,ampvt3,phasevt3,a,vorticity
       real ampvr1,phasevr1,ampvr2,phasevr2,ampvr3,phasevr3
       real rad_v,vrmean,u,meanu,v,special,rg,meanflow
       real meanangle,meanv,r_t,vta,angle,theta_t,asymvr
       character outfiledz*80,outfileve*80,inputfile*80
       common /ideal/vorticity,a,meanflow,meanangle,phasevt1,ampvt1,
     X               phasevt2,ampvt2,phasevt3,ampvt3,phasevr1,ampvr1,
     X               phasevr2,ampvr2,phasevr3,ampvr3

      
       PI=ACOS(-1.)
       dtor=pi/180.
       rt=sqrt(xc*xc+yc*yc)
       theta_t=atan2(yc,xc)/dtor
c       print *, theta_t
       meanu=meanflow*sin(meanangle*dtor-pi)
       meanv=meanflow*cos(meanangle*dtor-pi)
c       print *, 'ana'
        do i=1,nx
          xg=(xd(i))
          yg=(yd(i))
          angle=360./float(nx)*float(i-1)*dtor
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
            VTa=0.5*vorticity*radius*(1.+asymvt)
c            rad_v=0.5*radius
            rad_v=-vrmean*(1.+asymvr)
          else
            VTa=0.5*vorticity*a*a/radius*(1.+asymvt)
c            rad_v=-0.5*a*a/radius
            rad_v=vrmean*(1.+asymvr)
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
          ddz(i)=sqrt(u*u+v*v)
           IF (abs(gdz(i)-dummy).gt.0.001) then
c             print *, angle/dtor,theta/dtor,vta,u,v,gdz(i)
c              print *, gdz(i), xd(i),yd(i),iori,jori,rg,rt
             if(flag.eq.1 .or. flag.eq.3) then !vdad
               gdz(i) = gdz(i)*rg/RT
             elseif(flag.ne.2 .and. flag.ne.4) then ! not gbvtd or reflectivity
               print *, 'flag out of range'
             endif
           endif
c           print *, gdz(i),xd(i),yd(i),xg,yg,rg,rt
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
        theta_t=atan2(yc,xc)
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

      subroutine intswap(oldint,newint)
      implicit none
      integer*2 intfix
      integer*2 oldint,int1,int2,newint
      intfix = 256
      int1 = mod(oldint,intfix)
      int2 = (oldint - int1)/intfix
      newint = int1*intfix + int2
      return
      end
















