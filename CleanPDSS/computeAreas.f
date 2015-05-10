
C  computeAreas.f vers 1.0  program to test the getLayerAreas()  subroutine
C  compile  gfortran  -g -w -fbounds-check computeAreas.f  -o computeAreas
      SAVE
C  Program expects a single input line  such as:
C   11  3010.240  0.00500  0.00000  0.00500    0.0   1402.0  102.   72.  101.   73.
C  from which the last 4 items are the diameters.


      character*300  inp300,jnp300   ! used to read in lines as chars.
      character*1    inpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)

      integer nlayers, ierr, nlayerEdges
      real  Areas(100),deltStr(100),deltSts(100)
      real  xstrain(100),ystress(100),energy(100)

      real*8 dtime,dtime2,ev1(100),ev2(100)
      real dummy1,dummy2,dummy3,press,pos,xYod,xYid,xXod,xXid


      read(5,*)nlayerEdges, dtime,dummy1,dummy2,dummy3,press,pos,
     &         xYod,xYid,xXod,xXid
      nlayers= nlayerEdges-1
C      xYod=102.
C      xYid=72.
C      xXod=101.
C      xXid=73.
      call getLayerAreas(nlayers,xYod,xYid,xXod,xXid,Areas,ierr)

C     Read in the stuff in a very rigid manner.  No extras!
C     read in the blank line
      read(5,"(a1)")inpone(1) ! 2nd line is always empty
      read(5,"(a1)")inpone(1) ! 3rd line, 1st char should be a "#"
C     Fetch the data
      do 20 i=1,nlayerEdges
        read(5,*)ilayEdge,deltStr(i),deltSts(i),ev1(i),ev2(i),idummy,
     &           xstrain(i),ystress(i),dtime2,dE,energy(i)
        if(ilayEdge .ne. i)then
          write(0,*)"#ERROR:computArea: InputFile: ",
     &        "layerEdge Nos. messed up"
          stop
        endif
        if(ev2(i) .ne. dtime)then
          write(0,*)"#ERROR:computArea: InputFile: ",
     &        "event Nos. messed up"
          stop
        endif
   20 continue
      
      write(6,40)
   40 format("#Layer     event2   aveStr  avSts "
     &       ,"avDeltStr avDeltSts Area avEnergy EnergyXarea ")
      do 60 i=1,nlayers
        aveStr=(xstrain(i+1)+xstrain(i))/2.0
        aveSts=(ystress(i+1)+ystress(i))/2.0
        aveDeltStr= (deltStr(i+1)+deltStr(i))/2.0
        aveDeltSts= (deltSts(i+1)+deltSts(i))/2.0
        aveEnergy=  (energy(i+1)+energy(i))/2.0
        energyXarea=aveEnergy * Areas(i)
C        write(6,50)i,ev2(i),xstrain(i),ystress(i),
C     &        xstrain(i+1),ystress(i+1),aveStr,aveSts,
        write(6,50)i,ev2(i),aveStr,aveSts,
     &        aveDeltStr,aveDeltSts,
     &        Areas(i),aveEnergy,energyXarea
   50   format(1x,i3,1x,f13.3, 2(1x,f8.5,1x,f6.1),
     &         1x,f6.1,1x,f6.2,1x,f7.1)
   60 continue

      stop
      end


      SUBROUTINE getLayerAreas(nlayers,xYod,xYid,xXod,xXid,Areas,iret)
      SAVE
C     Compute the area of each layer.  Assume circle and that xYod,xYid
C     represent the  outside and inside diameters of circular tube.
C     i.e. no adjustment for tube flattening.

C     Ref.1: http://planetcalc.com/1421/
C     Ref.2:   mathworld.wolfram.com/CircularSegment.html
C     (in Ref.2  it is not clear what is meant by " cos^-1  " 
C      thus Ref.1 is followed here)

      real*4  Areas(100)
      real*4  xYod,xYid,xXod,xXid
      integer*4 nlayers, iret

      iret=0   ! error return indicator
      maxlayers=100


      if(nlayers .gt. maxlayers)then
        write(0,90)nlayers,maxlayers
        write(6,90)nlayers,maxlayers
   90   format("#ERROR:pdStressStrain:getLayerAreas:  Too many layers"/
     &         "# requested: req.= ",i4,"     max= ",i4/
     #         "# recompile subroutine getLayerAreas")
        iret=1
        stop   ! for now just simply stop here.
      endif


C     We will work our way from top to bottom of the cross section.
C     There are   nlayers.

      pi=3.1415927
      deltaH=xYod/float(nlayers)
      R= xYod/2.0   ! radius outside
      Ri=xYid/2.0   !       inside
      wallThick= (R-Ri)
      totalBigArea= pi*R**2
      totalInArea= pi*Ri**2
      totalWallArea=  totalBigArea-totalInArea
      sumArea=0.0
      write(6,80)xYod,xYid, xXod,xXid,deltaH,R,Ri,wallThick,
     &      totalBigArea,totalInArea,totalWallArea
   80 format("#getLayerAreas #xYod= ",f7.1/
     &       "#getLayerAreas #xYid= ",f7.1/
     &       "#getLayerAreas #xXod= ",f7.1/
     &       "#getLayerAreas #xXid= ",f7.1/
     &       "#getLayerAreas #deltaH= ",f7.1/
     &       "#getLayerAreas #R=  ",f7.1/
     &       "#getLayerAreas #Ri= ",f7.1/
     &       "#getLayerAreas #wallThick= ",f7.1/
     &       "#getLayerAreas #totalBigArea= ",f7.1/
     &       "#getLayerAreas #totalInArea= ",f7.1/
     &       "#getLayerAreas #totalWallArea= ",f7.1/
     &      )

C     We could just do half  and map to other half,  but in future we
C     may have a non-symmetrical section,  thus do in  nlayer steps

      h=0.0

      do 1000  i=1,nlayers 
      h=h+deltaH
C     for this cord compute the angle between lines that go from
C     circle center to cord ends.
      alpha=2.0*acos(1.0-(h/R) )
C     Now get the cord segment area.  eq. 1  of  Ref.2 (see above)
      areaBig= 0.5 *(R**2) * (alpha- sin(alpha) )

      hinside=h-wallThick
      write(6,94)h,hinside,alpha,areaBig
   94 format("#getLayerAreas #h= ",f7.1," hinside= ",f7.1,
     &       " alpha= ",e14.7," areaBig= ",f10.1)
C     check if we are inside the hollow portion
      if(hinside .gt. 0.0 .and. hinside .lt. xYid)then  
C       Yes, inside hollow: compute the area of the inner hole segment
        beta=2.0 * acos(1.0- (hinside/Ri) )
        areaIn= 0.5 * (Ri**2) * (beta- sin(beta) )
C       We can subtract the two areas
        xNetArea= areaBig-areaIn
        write(6,95)beta,areaIn,xNetArea
   95   format("#getLayerAreas #beta= ",e14.7," areaIn= ",f10.1,
     &         "  xNetArea= ",f10.1,"  of segment")
        goto 800
      endif

C     No,  we are either within the top wall   or the bottom wall
      if(hinside .le. 0.0)then  ! inside top wall
C        area is just  areaBig
         xNetArea=areaBig
         goto 800
      endif

C     inside bottom wall?
      if(hinside .ge. xYid)then  ! yes
         xNetArea=areaBig - totalInArea
         goto 800
      endif

C     We should never get to this point in logic.  So just in case:
      write(0,450)i,h,hinside,xYod,xYid
      write(6,450)i,h,hinside,xYod,xYid
  450 format("#ERROR:pdStressStrain:getLayerAreas:Logic error:"/
     &       "# layer= ",i4," h= ",f7.2," hinside= ",f7.2,
     &       "  Yod= ",f7.2,"   Yid= ",f7.2 )
      stop

C     Area computation complete for full size of tube. Get Layer area:
  800 continue
C     We have the total area of the inside and outside tube circles given
C     this value  of h
C     Now we need to get the area of just this layer  by subtracting the
C     previous computed  sumArea
      Areas(i)=xNetArea - sumArea
      write(6,810)i,Areas(i),xNetArea,sumArea
  810 format("#getLayerAreas #i= ",i3,"  Areas(i)= ",f10.1,
     &       "   xNetArea= ",f10.1,"  Previous sumArea= ",f10.1)
      sumArea=xNetArea

 1000 continue

C     we are done with all layers
      return
      end

