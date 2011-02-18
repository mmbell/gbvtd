       integer nx,ny,nz,nel,ix,iy,iz
       integer maxnr,maxnth,maxflds,maxnel
c       parameter(nx=81,ny=91,nz=1)
       parameter(nx=241,ny=241,nz=20)
       parameter(maxnr=1000,maxnth=500,maxflds=10,maxnel=16)
c       parameter(maxnr=1000,maxnth=1000,maxflds=16,maxnel=20)
       real dz3d(maxnr,maxnth,maxnel),dzcappi(nx,ny,nz),dummy,special
       integer minute(maxnth,maxnel),second(maxnth,maxnel)
       integer hour(maxnth,maxnel)
       real adv_speed, adv_angle
       real zmax,zmin,xmax,xmin,ymax,ymin,vecappi(nx,ny,nz)
       real xgridsp,ygridsp,zgridsp,zt(nz),vrelcappi(nx,ny,nz)
       real vrel3d(maxnr,maxnth,maxnel),angle(maxnth,maxnel)
       real startang(maxnel),beamsep(maxnel),stopang(maxnel)
       real altrad,range_first
       real xc, yc, xd(nx,ny), yd(nx,ny),elev(maxnel)
       real gatesp(maxnel)
       Integer nr(maxnel),nth(maxnel),tndo_flag
       integer*2 id(510),sf,cf
       character*2 flag(2)
       common /xyz/ ix,iy,iz,zt,xgridsp,ygridsp,zgridsp,angle
       common /header1/ sf,cf,id
       common /constant1/ dummy,special
       common /flag/ flag
       common /para/ startang,stopang,beamsep,gatesp,elev,nr,nth,nel
	 common /adv/ hour,minute,second,adv_speed,adv_angle
       common /indices/ range_first,tndo_flag
