C  pdStressStrain.f   vers. 1.0  Plate Edge Crack Prop.  FAC nov 21 2013
      SAVE
C  Push-Down List cyclic deformation material memory program.
C Compile:  gfortran  -g -w -fbounds-check pdStressStrain.f  -o pdStressStrain
C Usage:   pdStressStrain  emMag  ebMag  <loadHistory >outputFile
C           where emMag and ebMag are scale factors for em and eb.

C  The program is made availble to help students envision what happens
C  to the stress-strain path during variable amplitude loading
C  Program based on program "RCROK" Conle 1979 PhD thesis pg.128 and its
C  2012/3 derivitives  for BS7910 simulation of crack prop. with memory.

C  Copyright (C) 2015  Al Conle
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation; either version 2 of the license, or (at
C your option) any later version.
C  This  file is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTA-
C BILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
C License for more details.
C  You should have received a copy of the GNU General PUblic License along
C with this program; if not, write to the Free Software Foundation, Inc.,
C 59 Temple Place -Suite 330, Boston, MA 02111-1307, USA. Try also their
C web site: http://www.gnu.org/copyleft/gpl.html
C Note that some subroutines included in this work are from GNU GPL licenced
C program:  http://fde.uwaterloo.ca/Fde/Calcs/saefcalc3.html
C
C vers. 1.0  Fix bug near 2350. Bad compute of Smax for SWaT  Mar 5 2015
C vers. 0.9  Introduce emMag  and  ebMag  magnification factors Feb2 2015
C Fork: pdStressStrain.f  is a fork from plateEdgeFlaw.f  Jan23 2015
C            Remove all crack stuff.  Only run the basic memory
C            model in STRAIN control (no Neuber).  Strain history input is
C            similar to push down list crack model input: time em  eb
C                time em  eb  P  Pos  Yod Yid Xod Xid
C            The code for Neuber transform  P --> e  is commented out but
C            left in the listing for future use.
C Version changes in plateEdgeFlaw.f:
C vers. 3.10 Replace getPeakLoads() s/r  to remove small cycles Oct 26 2013
C            and align the end and begining of points in the history block.
C            Also introduce (but not yet use) a cycle repetition factor.
C            This latter requires a small change in hist. plot of makereport5
C Fork  plateEdgeFlaw.f  from plateLongSurfFlaw.f  Aug. 31 2013
C       Edge is very similar except that we have a/W  rather than a/B. 
C       Also no weld feature code for Mkm and Mkb.  Also fw=M=1
C vers. 3.06   Other catch-up version fixes:
C       3.03 change "a" accumulator to  REAL*8
C            Add output of crack info to a binary file.      
C            Binary file read for  FAD  post processing.
C       2.30 Add Pm  and Pb  stress to the #crk= line printout
C       2.24 Activate SAVELEVEL code to reduce output.  
C vers. 2.23 In s/r getPeakLoads() print out "Filtered" when no changes.
C            to allow makereport1's grep  to function properly  Feb28 2013
C vers. 2.21 Correct: material file name read error.  Jan 6 2013
C Fork to  plateEdgeFlaw.f  from plateWeldflaw.f   Dec.31 2012 vers. 2.2
C       Remove S/R: readMmMb00, readMmMb90, fwSurfFlaw, readMkmMkb00,
C                   readMkmMkb90, getMmMb00a90, getMkmMkb00a90
C       Add routine getLongMmMb()

C vers. 2.1  Correction: multipy Y(sigma) by sqrt(pi*a)  Dec.29 2012
C vers. 2.0  Discretize lobj00, ldo00 lobj90,ldo90, dld00,dld90
C vers. 1.1  Divide damage or dadn  by 2.0 to make it per 1/2 cycle.
C            Add output for loads to be rainflow counted and used for
C            crack initiation analysis.
C vers. 1.00 runs ok.  
C vers. 0.91 Remove the old unused functions
C vers. 0.9 Fork: plateWeldflaw.f and plateWeldflaw+ss.f
C           In this version we have eliminated the local stress-strain stuff.
C           Thus it does the memory thing using deltaK00 and deltaK90 with no
C           tracking of the local ss.  Nov.17 2012
C           Local ss stuff is commented out with "Css"
C vers. 0.8 Create two parallel push-down lists. One for 00 and one for 90 deg.
C           crack tip.  Each has its own deltaK, stress, strain etc.
C           Each is P.D. list counted using its deltaK.
C vers. 0.6 replace Mkm and Mkb read s/r's. Add getMm and getMk s/r's.
C              (extract from testgetMmMkm.f test prog.)
C       Add fw  read and interpolate s/r. subr_getfw.f -> fwSurfFlaw()
C       Add peak pick  s/r. getPeakLoads()
C       Add getStress2Strain(), getLoad2StressStrain()
C       De-activate old functions from thesis version: XKP(),SMITH(),FLD(),DET()
C vers. 0.5 places read and interpolation for Mkm and Mkb into S/Rs. 
C           also adds read and interp. for Mm  and Mb


C In general:
C  1. Counts using Nominal "Load" == etotali =  em+eb
C  2. Computes Epsilon from Load.
C  3. Computes Sigma from Epsilon 
C  4. Computed damage as crack length  ( not in this version)

C  Material behaviour is stabilized.  No cyclic mean stress relaxation
C  or cyclic hardening or softening is performed in this version.  
C  Also no "disappearing" memory of prior deformation due to crack prop. is done.

C  Explanation of primary variables used=
C    tlimL, climL = 1 dimensional arrays of vectors containing
C       the push-ddown list loads.  TLIML is the tensile, CLIML is compr.
C  tlstr, clstr = the p.d. list associated Tens and Compr. strains
C  tlsts, clsts =               "          "               stresses
C  tldam, cldam =               "          "               1/2 cycle fat. damage

C  nptt,nptc =  push-down list pointers.
C  lorg, eorg, sorg  The load, strain, stress origin of a given 1/2 cycle
C  lobj, eobj, sobj         "           "      end point        "    "
C  dld, de,ds  =     the change in "    "  during a 1/2 cycle
C  totdam  =  total damage. In this program is crack length
C  ef = monotonic fracture strain
C  maxpd = maximum no of points possible in PD list
C  maxpoints = max no. of pts in Load history
C  maxnio = max no. reversals at which output is requested


C     Save the fitted Strain-Life-Stress  digital curves here:
C     Check dimensions in the various subroutines too    !!!!!!!!!!!!
      logical debugMat
      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain,Emod
C     Use these to transform to ALOG  to speed up interpolations
      real logStrainAmp(250),logLifecycles(250),logStressAmp(250),
     &     logSwtAmp(250),logSelasticAmp(250)
      Common/MaterialA/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 logStrainAmp,logLifecycles,logStressAmp,
     &                 logSwtAmp,logSelasticAmp,
     &                 ndata,FractureStress,FractureStrain,Emod


      real tstrain(250),tstress(250)
C          tstrain & tstress are used only for loop plotting and passing
c          from s/r getloop


Crack      logical logiExist  !used in opening  fadInput.rand   file

C     These are the push-down list members:
      real*4 tlimL90(4000),climL90(4000)
     & ,tlstr90(4000), clstr90(4000),tlsts90(4000),clsts90(4000)
     & ,tldam90(4000),cldam90(4000)
      integer itlnrev90(4000),iclnrev90(4000)
C     The damage (crack length) counters need to be real*8 because we
C     may be trying to add a very small crack increment to a large crack length.
      real*8 tltotCrk90(4000),cltotCrk90(4000) !p.d. list damage storage.
      real*8 dtotdam90, daminc90, damold90 !these are total damage counters
C      The reason for pushing the tltotCrk and tlnrev numbers are for future
C      disappearing memory due to crack extension. 
C

      character*300  inp300,jnp300,Cwebpage ! used to read in lines as chars.
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)

      character*10 name1, name2, stressunits
      character*10 strainunits, lifeunits, sorttype
      character*30 names30(10), firstfield, ctail
      character*30 cfiletype, Cdatatype
      character*80 argv, matfile!, fwfile, kprimefile, dadnfile
      character*80 ctID,probType,histfile!,dadnType,dadnTableFile

C     Use these to read the command line arguments emMag  ebMag
      integer*4 iargc,nargc

c     Load history storage.  If changing dimensions, also change same 
C                             in S/R getPeakLoads() if used.
C     Note!: The following will actually be  STRAINS :
C     Stress membrane, bending, and total storage:
      real*8 dtime(5000)
      real*4 em(5000),eb(5000),etot(5000),press(5000),pos(5000),
     *       Yod(5000),Yid(5000),Xod(5000),Xid(5000)
      integer*4 iloadflag(5000) ! = 1 if ok, =0 not used, or intValue if repeat
      real*4 etotmax,etotmin,etotwindow
      real*8 dxtime,doldtime deltaTime dinsertTime !used to add tiny 1/2 cycs.
      integer*4 nloads
      logical debugLoads
      common /LOADS/ dtime,em,eb,etot,press,pos,Yod,Yid,Xod,Xid,
     &         iloadflag,nloads,debugLoads,etotmax,etotmin,etotwindow


      real*4 ldo90,lobj90

c     Not really used here, but left in code:
c     Reversals at which PDList output occurs
      integer*4 nio(23)/2,3,4,5,8,70,100,200,500,1000,2000,5000,
     &    10000,20000,50000,100000,200000,300000,500000,1000000,
     &    5000000,10000000,100000000/
c     After last number in nio() put out data every nio(maxnio)
c     Note that as the crack propagates away from the Notch stress 
c     raiser that the PDLists may be VERY large: = maxpd
c     Thus we will need to place a limit on the no. PDlist printed out.

      logical prin
      logical stabl
      common /STAB/ stabl

      pi=3.1415927

      debugMat=.false.     !use to debug the Neuber and stress-strain portion
      debugLoads=.false.   ! use to debug load read and peakpick section.

C     Storage area limits. Change this if you change above real* and int*
      maxmatdata=1500 ! Stress,strain,life store max

      maxloads=5000   ! max load storage
      maxpd=4000     ! Push down list storage max
      m=maxpd
      maxnio=23       ! Max store for Rev data output
      nioLast=nio(maxnio)

C     Test if -Wuninitilized  warning works in gfortran:
C      i=ifix(x)
C     Works!
    
C      Initilize some variables.
C      Origin and zero
      eorg=0.0
      sorg=0.0
      nact=0

      nptc90=0
      nptt90=0
      do 20 i=1,maxpd
        tlimL90(i)=0.    !peak in tensile nominal load or deltaK
        climL90(i)=0.    !peak in compr.  nominal load or deltaK
Css        tlstr90(i)=0. !peak in tens. strain
Css        clstr90(i)=0. !  "     comp.   "
Css        tlsts90(i)=0. !peak in tens. stress
Css        clsts90(i)=0. !  "     comp.   "
        tldam90(i)=0.    !damage or delta a for that ramp
        cldam90(i)=0.
        tltotCrk90(i)=0. ! total damage at that point
        cltotCrk90(i)=0.
        itlnrev90(i)=0.   ! rev. at which the above occured.
        iclnrev90(i)=0.
   20 continue

      smallStrain=0.00005 !size of artificial reversal

C---------------------------  Run time input data------------------
  184 continue
      write(6,185)
      write(0,185)
  185 format("# pdStressStrain.f vers. 1.0"/
     & "#Usage: pdStressStrain  emScale ebScale <histfile  >outfile"/)

      nargc = iargc()
C     Note that in HP fortran the arg numbers must be changed by -1
C     and that iargc probably includes the "pdStressStrain" as an arg.
      if( nargc .ne. 2)then
        write(0,*)" pdStressStrain:  usage ERROR"
        write(6,*)" pdStressStrain:  usage ERROR"
        write(0,186)
        write(6,186)
  186  format(
     &   "#Usage:      pdStressStrain emScale ebScale <histfile ",
     &   " >outputFile"//
     &   "# Where *Scale is the multiplier applied to all history",
     &   " points. "/
     &   "# NOTE!: The above factors will be applied to the Strains"/
     &   "# history effectively AFTER the factors "/
     &   "# in the ""pdss.env"" file are applied "
     &   "# Stopping now."
     &   )
        stop
      endif

C       The first arg is the history multiplication factor
        jvect=1
        call getarg(jvect,argv)
        read(argv,*,err= 178)emMag
        write(6,*)"#emMag= ",emMag
        if(emMag .eq. 0.0)go to 178

        jvect=jvect+1
        call getarg(jvect,argv)
        read(argv,*,err=179)ebMag
        write(6,*)"#ebMag= ",ebMag
C       Both emMag and ebMag have been read in successfully
        go to 190

C     Bad arguments in command line
  178 write(0,*)"# ERROR: bad scaleValue emMag argument=",argv
      stop
  179 write(0,*)"# ERROR: bad scaleValue ebMag argument= ",argv
      stop

  190 continue
      write(6,191)
      write(0,191)
  191 format("#Opening pdss.env file...")
C     Initilize the things to be read to checkable items
C     to make sure they have been entered.
      probType= " "  !not really used
C     The calling script will read the material file name from the tube
C     segment's history file and place the
C         results into the variable matfile:
      matfile= " "  ! if this does not get changed it will cause error.
      isegment= 0   ! segment number will be read into this.

C     The history file will be read from standard input (device 5)
      histfile= "stdin"
      xmagfactorm= -1.0e20    ! Set no data read flag
      xmagfactorb= -1.0e20
      xmeanAddm=   -1.0e20  ! Hopefully no one will ever use this :(
      xmeanAddb=   -1.0e20  !   (unlikely since it is a strain)
      
      maxHistReps= -1
      blockSkip= -9999.
      isavelevel= 0    !this is the default level if none is specified



C -----------   Open and read in the pdss.env  problem environment file

C     In this file all lines should begin with a #  or are blank lines
      open(unit=10,file="pdss.env")
      ninput=0
  800 continue
c     Loop back to here for next input line.

      read(10,"(a300)",end=380)inp300
      ninput=ninput+1
Cdebug      write(0,*)" read input line ",ninput

C     Check for blank line
      if(inp300.eq." ")then
C        write(6,"(a1)")" "
        go to 800
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 805 i=1,300
           if(inpone(i).eq." ") go to 805
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 803 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  803        continue
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 820
          else
C         The first non blank is not a #
          write(0,*)"# Skipped garbage in *.env file: line no.= ",ninput
          write(0,*)"# Text is= ",inp300
          write(6,*)"# Skipped garbage in *.env file: line no.= ",ninput
          write(6,*)"# Text is= ",inp300
          go to 800
          endif
  805    continue


      endif
  820 continue
C     First char was a #
C     pdprop.env file input line has a # in 1st col.  
C     See if it has a keyword.

      if(inpone(1).ne."#")then
C       something is bad in program
        write(0,*)" ERROR 820, sorry prog. messed up call ? admin"
        stop
      endif

C     Ok, its a nice comment.  Figure out if its a special tag.
      read(inp300,*)firstfield


C     Type is not used in this program. left here for future.
      if(firstfield .eq."#TYPE=" .or.
     &   firstfield .eq."#Type=" .or.
     &   firstfield .eq."#type=" )then
         read(inp300,*) firstfield, probType
         write(jnp300,833)probType
  833    format("#found: #TYPE= ",a80)
         call WR(6,jnp300)
         go to 800
      endif


      if(firstfield .eq."#MAXREPS=" .or.
     &   firstfield .eq."#Maxreps=" .or.
     &   firstfield .eq."#maxreps=" )then
         read(inp300,*) firstfield, maxHistReps
         write(6,850)maxHistReps
  850    format("#found: #MAXREPS= ",i10)
         go to 800
      endif

      if(firstfield .eq."#MAGFACTOR_m=" .or.
     &   firstfield .eq."#Magfactor_m=" .or.
     &   firstfield .eq."#magfactor_m=" )then
         read(inp300,*) firstfield, xmagfactorm
         write(6,852)xmagfactorm
  852    format("#found: #MAGFACTOR_m= ",e14.7)
         go to 800
      endif

      if(firstfield .eq."#MAGFACTOR_b=" .or.
     &   firstfield .eq."#Magfactor_b=" .or.
     &   firstfield .eq."#magfactor_b=" )then
         read(inp300,*) firstfield, xmagfactorb
         write(6,854)xmagfactorb
  854    format("#found: #MAGFACTOR_b= ",e14.7)
         go to 800
      endif

C     It is generally very confusing if a user applies these.
      if(firstfield .eq."#MEANADD_m=" .or.
     &   firstfield .eq."#Meanadd_m=" .or.
     &   firstfield .eq."#meanadd_m=" )then
         read(inp300,*) firstfield, xmeanAddm
         write(6,855)xmeanAddm
  855    format("#found: #MEANADD_m= ",e14.7)
         go to 800
      endif

      if(firstfield .eq."#MEANADD_b=" .or.
     &   firstfield .eq."#Meanadd_b=" .or.
     &   firstfield .eq."#meanadd_b=" )then
         read(inp300,*) firstfield, xmeanAddb
         write(6,856)xmeanAddb
  856    format("#found: #MEANADD_b= ",e14.7)
         go to 800
      endif
C

C     Not really used in pdStressStrain.f    future...
      if(firstfield .eq."#SAVELEVEL=" .or.
     &   firstfield .eq."#Savelevel=" .or.
     &   firstfield .eq."#savelevel=" )then
         read(inp300,*) firstfield, isavelevel
         write(6,894)isavelevel
  894    format("#found: #SAVELEVEL= ",i2)
         go to 800
      endif


  350 continue
C     None of the above?  Then it must be a plain old
C     comment line.   Write it out too.
C       Count backwards and see where the last char is
        do 360 i=1,300
          j=300-(i-1)
          if(inpone(j).ne." ")then
C           found last char
            lastloc=j
            go to 362
          endif
  360   continue

  362   continue
      write(6,"(300a1)")(inpone(i),i=1,lastloc)
C     Go read another line
      go to 800

C     End of pdprop.env file reached
  380 continue
C     Check if critical items have been read in.
      close(unit=10)
      istop=0
      if(probType .eq. " ")then
        write(0,*)"#ERROR: #Type= not found."
        write(6,*)"#ERROR: #Type= not found."
        istop=1
      endif
C 
      if(xmagfactorm .eq. -1.0e20)then
        write(0,*)"ERROR:  #MAGFACTOR_m=  not found"
        write(6,*)"ERROR:  #MAGFACTOR_m=  not found"
        istop=1
      endif
      if(xmagfactorb .eq. -1.0e20)then
        write(0,*)"ERROR:  #MAGFACTOR_b=  not found"
        write(6,*)"ERROR:  #MAGFACTOR_b=  not found"
        istop=1
      endif
      if(xmeanAddm .eq. -1.0e20)then
        write(0,*)"ERROR:  #MEANDADD_m=  not found"
        write(6,*)"ERROR:  #MEANDADD_m=  not found"
        istop=1
      endif
      if(xmeanAddb .eq. -1.0e20)then
        write(0,*)"ERROR:  #MEANDADD_b=  not found"
        write(6,*)"ERROR:  #MEANDADD_b=  not found"
        istop=1
      endif


      if(istop.eq. 1)then
         write(0,*)"# Stopping..."
         write(6,*)"# Stopping..."
         stop
      endif





      izout=1 ! used in fatigue runs to control output cyc. no. Not used here.

C     initilize the functions
      nrev=0
c
C------------------Get Load history file------------------------------
      iupdown=0   ! 0=flat  1=going into tension,  -1= going into compression

C     Expected format:
C     #Comment 1
C     #ID= XYZ5134     #Tube Identification
C     #Material=   mergedSAE950X.html    # material fatigue data file for this segment.
C     #SEGMENT=  1030   # position of segment measured from ??? Reel Zero ?
C     #Comment line(s) 
C     #secs?  strain strain  mpa       ft?     mm    mm   mm   mm
C     #time    em     eb       P       Pos    Yod    Yid  Xod  Xid
C     0         0     0        0      -1030    102   72  101   73    #comments...
C     598.      0     0        0      -1030    102   72  101   73    #comments...
C     601.      0     0.0080   0      -1030    102   72  101   73    #comments...
C     ....  etc


C     Membrane and Bending strains will be altered by the #MAGFACTORm  etc lines
C     in the pdss.env file and the calling args.

      write(0,901)
      write(6,901)
  901 format("# Getting Strain history file from std.input ")
C      open(unit=10,file=histfile)

      ninput=0
      nloads=0
  900 continue     !Loop back to here for next input line.
      read(5,"(a300)",end=950)inp300
C      call WR(6,inp300)  ! debug
      ninput=ninput+1

      if(inp300.eq." ")then    ! Check for blank line
C        write(6,"(a1)")" "
        go to 900
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 905 i=1,300
           if(inpone(i).eq." ") go to 905
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 903 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  903        continue
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 930
           else
C            The first non blank is not a #
             go to 910
           endif
  905      continue
           go to 900  !assume garbage in line.

  910      continue
C        1st non blank is not a #  Check if its a number.
         if(inpone(loc).ne."+" .and.
     &      inpone(loc).ne."-" .and.
     &      inpone(loc).ne."." .and.
     &      inpone(loc).ne."0" .and.
     &      inpone(loc).ne."1" .and.
     &      inpone(loc).ne."2" .and.
     &      inpone(loc).ne."3" .and.
     &      inpone(loc).ne."4" .and.
     &      inpone(loc).ne."5" .and.
     &      inpone(loc).ne."6" .and.
     &      inpone(loc).ne."7" .and.
     &      inpone(loc).ne."8" .and.
     &      inpone(loc).ne."9"  )then
C           It must be a letter, not a number. Time to bomb out.
            write(0,915)ninput
            write(6,915)ninput
  915       format(" ERROR stdin: input line no. ",I5,
     &      " not a # and not a number. Fix data in strain hist. file")
            stop
         endif
C        Ok, its a number ----------------------
         read(inp300,*,err=925)dxtime,xem,xeb,xpress,
     &     xpos,xYod,xYid,xXod,xXid  ! dxtime is real*8

         nloads=nloads+1
         if(nloads .gt. maxloads)then
           write(0,916)nloads
           write(6,916)nloads
  916      format(" Error too many data points:",i5,
     &      "recompile pdStressStrain.f or reduce strain history data")
           stop
         endif
           xetot=xem+xeb
C        Adjust to mag. and mean shift in *.env file
         xem= (xem*xmagfactorm + xmeanAddm )*emMag
         xeb= (xeb*xmagfactorb + xmeanAddb )*ebMag
C        Find the total strain. also for Rainflow counting
         xetot=xem+xeb
         if(nloads .eq. 1)then  ! 1st point, set to old and set direction.
           etotmax=xetot ! Set here for max and min scan
           etotmin=xetot
           if(xetot .ne. 0.) then ! establish an initial direction
             if(xetot .lt. 0.) iupdown=-1
             if(xetot .gt. 0.) iupdown=1
C            otherwise leave iupdown=0  its initial value. i.e. flat history.
           endif
           goto 919  ! save the items read in.
         endif

C        If nloads  > 1  then see if we have a direction change (i.e. a reversal)
         if(xetot .eq. etot(nloads-1) )then   ! same total strain as before
C          just write this set into the pile.  iupdown does not need change
           goto 919 !save items
         endif
         
C          --------------- load going up---------------
         if(xetot .gt. etot(nloads-1) )then ! going higher in tension.
          if(iupdown .eq.0) then ! Was flat change to Tens. and save
             iupdown= 1
             goto 919
          endif
          if(iupdown .eq. -1)then ! Was in Comp.  Change to Tens and save
             iupdown= 1
             goto 919
          endif
          if(iupdown .eq. +1)then ! Was already going up.  This is extension!
C           Insert a very small unload reversal just after the previous rev.
            j=nloads-1
            doldtime=dtime(j)  ! all time stuff is real*8
            deltaTime=dxtime-doldtime
            dinsertTime=doldtime + deltaTime/20.0
            emInsert=em(j)
            ebInsert=eb(j) - smallStrain  ! note the minus sign
            etotInsert= emInsert+ebInsert
C           C write out this small reversal in strain
            write(6,917)dinsertTime,emInsert,ebInsert,etotInsert,
     &           press(j),pos(j),Yod(j),Yid(j),Xod(j),Xid(j)
  917       format("#history.adj ",f15.3,1x,3(f8.5,1x),f6.1,1x,f8.1,
     &             4(1x,f5.0), " #an insert")
C           nloads has already been incremented after data read line above.
C           Store the small cycle items
            dtime(nloads)=dinsertTime
            em(nloads)=  emInsert
            eb(nloads)=  ebInsert
            etot(nloads)=etotInsert
            press(nloads)=press(j)
            pos(nloads)=pos(j)
            Yod(nloads)=Yod(j)
            Yid(nloads)=Yid(j)
            Xod(nloads)=Xod(j)
            Xid(nloads)=Xid(j)
C           iupdown will remain = +1
C           Save the new (extended) point's data
            nloads=nloads+1  ! need to increase. Old value was used above.
            goto 919
          endif
         endif  !  for check if xetot > etot(nloads-1)


         if(xetot .lt.etot(nloads-1))then  ! We are going down next.
           if(iupdown .eq. 0) then ! Was flat change to Comp. and save
             iupdown= -1
             goto 919
          endif
          if(iupdown .eq. 1)then ! Was in Tens.  Change to Comp and save
             iupdown= -1
             goto 919
          endif
          if(iupdown .eq. -1)then ! Was already going DOWN.  This is extension!
C           Insert a very small unload reversal
            j=nloads-1
            doldtime=dtime(j)
            deltaTime=dxtime-doldtime
            dinsertTime=doldtime + deltaTime/20.0
            emInsert=em(j)
            ebInsert=eb(j) + smallStrain   ! note the plus sign
            etotInsert= emInsert+ebInsert
C           C write out this small reversal in strain
            write(6,917)dinsertTime,emInsert,ebInsert,etotInsert,
     &           press(j),pos(j),Yod(j),Yid(j),Xod(j),Xid(j)
C           nloads has already been incremented after data read line above.
            dtime(nloads)=dinsertTime
            em(nloads)=  emInsert
            eb(nloads)=  ebInsert
            etot(nloads)=etotInsert
            press(nloads)=press(j)
            pos(nloads)=pos(j)
            Yod(nloads)=Yod(j)
            Yid(nloads)=Yid(j)
            Xod(nloads)=Xod(j)
            Xid(nloads)=Xid(j)
C           iupdown will remain = -1
C           Save the new points data
            nloads=nloads+1  ! need to increase. Old value was used above.
            goto 919
          endif
        endif  ! end of  xetot < etot(nloads-1) if

C       We should never get to here.  Throw an error.
        write(6,*)"#ERROR: Processing history logic: ninput, nloads=",
     &            ninput,nloads
        call WR(6,inp300)
        write(0,*)"#ERROR: Processing history logic: ninput, nloads=",
     &            ninput,nloads
        call WR(0,inp300)
        stop

  919   continue !Now save the point that came in at the read statement.
        if(xetot .gt. etotmax) etotmax=xetot ! check for max and min change
        if(xetot .lt. etotmin) etotmin=xetot
        dtime(nloads)=dxtime
        em(nloads)= xem
        eb(nloads)= xeb
        etot(nloads)=xetot
        press(nloads)=xpress
        pos(nloads)=xpos
        Yod(nloads)=xYod
        Yid(nloads)=xYid
        Xod(nloads)=xXod
        Xid(nloads)=xXid
        write(6,920)dxtime,xem,xeb,xetot,xpress,xpos,xYod,xYid,xXod,xXid
  920       format("#history.adj ",f15.3,1x,3(f8.5,1x),f6.1,1x,f8.1,
     &             4(1x,f5.0) )
        goto 900   ! fetch the next input line




  925    write(0,926)nloads,inp300
         write(6,926)nloads,inp300
  926    format("#READ ERROR: history file: data line no.= ",i9/
     &   "#LINE= ",a300/ "# Stopping now.")
         stop

      endif !  end of if  first letter not #


  930 continue
C        Yes, first char was a #
C        See it is a special tag=  line
         read(inp300,*)firstfield
C        We are looking for  #ID=   #Material=   #SEGMENT=
        if(firstfield .eq."#ID=" .or.
     &    firstfield .eq."#Id=" .or.
     &    firstfield .eq."#id=" )then
          read(inp300,*) firstfield,ctID
          write(jnp300,931)ctID  !this is just to strip extra chars
  931     format("#ID= ",a80)
          call WR(6,jnp300)
          go to 900
        endif

        if(firstfield .eq."#MATERIAL=" .or.
     &    firstfield .eq. "#Material=" .or.
     &    firstfield .eq. "#material=" )then
          read(inp300,*) firstfield,matfile
          write(jnp300,932)matfile  !this is just to strip extra chars
  932     format("#MATERIAL= ",a80)
          call WR(6,jnp300)
          go to 900
        endif

        if(firstfield .eq."#SEGMENT=" .or.
     &    firstfield .eq."#Segment=" .or.
     &    firstfield .eq."#segment=" )then
          read(inp300,*) firstfield,isegment
          write(6,933)isegment
  933     format("#SEGMENT= ",i10)
          go to 900
        endif

        go to 900


  950    continue    ! end of file comes here
C        Check that we have the stuff we need
         istop=0
         if(nloads.eq.0)then
            write(6,941)
            write(0,941)
  941       format("#Error: histfile: contains no data !")
            istop=1
         endif
         if(matfile .eq. " ")then
           write(6,942)
           write(0,942)
  942      format("#Error: histfile: missing #MATERIAL= ")
           istop=1
         endif

         if(ctID .eq. " ")then
           write(6,943)
           write(0,943)
  943      format("#Error: histfile: missing #ID= ")
           istop=1
         endif

         if(isegment .eq. 0 )then
           write(6,944)
           write(0,944)
  944      format("#Error: histfile: missing #SEGMENT= integer .")
           istop=1
         endif
         if(istop .eq. 1) stop

C        All data is in. Close file
C         close(unit=10)
         write(6,952)nloads,etotmax,etotmin
         write(0,952)nloads,etotmax,etotmin
  952    format("#history #Data input completed. nloads= ",i5/
     &   "#history #etotMax= ",f7.5/"#history #etotMin= ",f7.5/
     &   "#history #Where etot = emembrane + ebending"//)
         


C------------- Read in the Material Fitted Fatigue file Table-------------------

C   The material fatigue file name is specified in the  strain history file.
C     matfile has been read in from the strain history file above.
      istrainCtrl= 1  ! Set up for a strain control simulation.
C         This means that   SelasticAmpl() = StrainAmpl()   thus no Neuber.
      call readStrainLifeStress(matfile,istrainCtrl, iret)





C  Future code for long histories
CC----------------------binary output file-------------------------------------
CC     Open the binary (direct access) record file for data output
CC     Get a filename for the standard input stream (doesnt work)
CC      inquire(unit=5, NAME=stdinFileName)
C      write(0,460)
C      write(6,460)
C  460 format(/"# Opening random access output file:  fadInput.rand ...")
C      inquire(file="fadInput.rand", EXIST= logiExist)
C      if(logiExist)then
C        write(0,*)"#ERROR: an old copy of file  fadInput.rand  exists."
C        write(0,*)"#   You need to rename it or remove it before we"
C        write(0,*)"#   can run the simulation.  Stopping now..."
C
C        write(6,*)"#ERROR: an old copy of file  fadInput.rand  exists."
C        write(6,*)"#   You need to rename it or remove it before we"
C        write(6,*)"#   can run the simulation.  Stopping now..."
C        stop
C      endif
C
CC     Ok,  previous copies do not exist. open it.
C      open(unit=60, file="fadInput.rand", access="direct",
C     &     form= "unformatted", status= "new", recl= 36 )
C      write(0,462)
C      write(6,462)
C  462 format("#Random Access output file: fadInput.rand   opened.")
C      nrecord=1  !this will be incremented by 1 upon 1st use.
CC                 The 1st record is written at the end of test.
C


C---------------------- Start Cycling--------------------------------------------
C     Load or nominal stress (in this case = Strain ) is used to model 
C     material memory. 
C     Program flow and logic are pretty much the same as concepts of
C     original program  rcrock.f from thesis.
      nblk=1
      nrev=0  !total rev counter
      nact=0  !tracks which rev we are on in a block


      lobj90=0. ! target of half cycle. Nominal Stress
      ldo90=0.  ! origin of half cycle. Nominal Stress
      eobj90=0.
      sobj90=0.
C     eo90 and so90 origins will be set below.

      StrOut=0.0 ! Used to track stress-strain output
      StsOut=0.0

      iupdown90=+1


C      dtotdam90=azero  !crack length "a" for 90deg depth crack point
      dtotdam90=0. ! in pdStressStrain.f this is just Miner's rule damage.
      damold90=0.
      daminc90=0.
C     note that there are two damage accumulators:  totdam90  and dtotdam90
C     They both contain the same numbers but dtotdam90   is REAL*8  in size
C     in order to better accumulate very small damage numbers  when the 
C     total damage is large.

C      if(.not.lactivatefw)xfw=1.0  ! See if FiniteWidth Corr. is shutdown
C       xfw=1.0    ! See BS7910 2005  M.3.5


 3000 continue  ! =================Top of Cycling loop====================
      ldo90=lobj90
      eo90=eobj90
      so90=sobj90
      
      nact=nact+1
      if(nptt90.eq.maxpd .or. nptc90.eq.maxpd )then
         write(0,118)maxpd,nptt90,nptc90
         write(6,118)maxpd,nptt90,nptc90
 118     format("#Error: one or more of the PushDown list counters >",
     &   "maxpd= ",i6," : nptt90=",i6," nptc90=",i6
     &   )
         stop
      endif

cccccccccccccccccccccccccccccccccccccccccccccc  Block Repeat
C     The user must make certain that end of block is  compatible with 
C     the begin of block.  The program will start at 0 stress(load) but
C     the history does not need to include this start point. The last
C     point in the history and the first and 2nd points MUST form reversals.
C    E.g:     2nd_last_pt    last_pt   1st_pt   2nd_pt
C                 +100        -90       90       -90       is allowed
C                 +100        +50       0         +50      is NOT allowed

 3001 continue
CC     Check if no change.  At times nact will point to items with no change
CC     in the strain values.  This would occur for example if there was a
CC     internal pressure test at begin with no change in em or eb:
C      if(nrev .eq. 1 .and. em(nact) .eq. 0.0 .and. 
C     &   eb(nact) .eq. 0.0)then
C         write(6,3002)dtime(nact),em(nact),eb(nact),etot(nact)
C     &             press(nact),pos(nact)
C 3002    format("#nullrev: ",f13.3,3(1x,f7.5),1x,f5.1,1x,
C         nact=nact+1  !zero range, skip this half cycle
C         goto 3001
C      endif
C
      if(nact .gt. nloads)then
C       End of block, increment the blk count and reset nact=1
        nblk=nblk+1  
C       Check for MAX. history repeats
        if(nblk .gt. maxHistReps)then
          if(totdam90 .gt. 0.0)then
           xreps=1.0/totdam90
           write(0,130) nblk,totdam90,xreps,nrev
           write(6,130) nblk,totdam90,xreps,nrev
  130      format("# Max no. of History Reps. Reached: nblk= ",i10/
     &      "#TOTDAM90= ",E14.7," allowed Repeats= ",f10.1," nrev= ",i9)
          else
           write(0,132)nblk,totdam90,nrev
           write(6,132)nblk,totdam90,nrev
  132      format("# Max no. of History Reps. Reached: nblk= ",i10/
     &     "#TOTDAM90= ",E14.7," allowed Repeats=      inf   nrev= ",i9)
          endif
          goto 9000  !stop
        endif
        nact=1 ! reset to repeat history
      endif
      

C       daminc and damold s  are used to skip a block if damage is same
C       Mostly for crack prop version.
        daminc90=dtotdam90 - damold90
        damold90=dtotdam90

      totdam90=sngl(dtotdam90)

      lobj90= em(nact)+eb(nact)
    
      if(lobj90.eq.ldo90 )then
C       Skip, its a nothing reversal.  It should probably be a delta
C       check.  (but what if one is zero, but the other is not?)
C        write(0,3007)nrev,nblk,nact
        write(6,3007)nrev,nblk,nact
 3007     format("#Warning: Same target as previous lobj90=ldo90 "
     &    ," at nrev=",i9," nblk=",i9," nact=",i9/
     &     "# Skipping this stress...")
        write(6,3008)dtime(nact),em(nact),eb(nact),etot(nact),
     &       press(nact),pos(nact),lobj90,eobj90,sobj90,totdam90
 3008   format(f13.3,1x,3(f7.5,1x),f5.1,1x,f7.0,1x,
     &         f7.5,1x,f7.5,1x,f6.1,1x,e14.7" #nullrev")
C       Although we have a null rev. we should print out the stress-strain
C       status for proper plotting on a time scale.
        write(6,3009)StrOut,StsOut,dtime(nact)
 3009   format("#plotloops ",f7.5,1x,f6.1,1x,f13.3)
        go to 3000
      endif

C     It is a valid reversal.  Start the ramp.
C     Crack length variables may have to be real*8  because adding
C     a very small da increment may not show up if crack is big.
      nrev=nrev+1
      if(nrev .eq. 1)then
        write(6,*)"#plotloops ",StrOut,StsOut,dtime(nact)
      endif
      totdam90=sngl(dtotdam90)
      if(isavelevel .gt. 0)then
        write(6,120)dtime(nact),em(nact),eb(nact),etot(nact),
     &     press(nact),pos(nact),lobj90,eobj90,sobj90,totdam90
  120   format(f13.3,1x,3(f7.5,1x),f5.1,1x,f7.0,1x,
     &         f7.5,1x,f7.5,1x,f6.1,1x,e14.7," #Target" )
      endif

CC     We always write to the binary out file. RecSize is 12*4 =48
C        nrecord=nrecord+1
C        write(60,rec=nrecord)nrev,totdam90,nblk,nact,
C     &             lobj90,xMm90,xMb90, 
C     &             stsMembrane,stsBending


C     Often too much output.  Turn this off if isavelevel <3
      if(isavelevel .lt. 3) go to 3050


C     Debug:  dump the push-down lists:
      write(6,121) nrev,totdam90,nblk,nact
  121 format("#PD Before REV= ",i9,"  SWT damage= ",E14.7
     &  ,3X,"NBLK=",i6," NACT= ",I6
     & /"#PD: ",4x,"CLDAM90",3x,"CLIML90" )
      do 3039 ipr=1,nptc90
 3039 write(6,122) ipr,cldam90(ipr),climL90(ipr)
  122 format("#PD:",1x,i2,1x,E12.5,2x,F9.2)
   
      write(6,123)
  123 format("#PD:",5x,"TLDAM90",4x,"TLIML90")
      do 3045 ipr=1,nptt90
      write(6,122) ipr,tldam90(ipr),tlimL90(ipr)
 3045 continue


 3050 continue
C ----------------------------------------------------------------
C
      dld90=abs(lobj90-ldo90)
Cdebug      write(6,*)"#stsm(),stsb(),nact=",stsm(nact),stsb(nact),nact
Cdebug      write(6,*)"#dld90,lobj90,ldo90=",dld90,lobj90,ldo90,nrev
C      if(dld90.le. 0.05)then   !check for no actual half cycle



C     Going up or down ?   -----------------Crack Direction 90 deg------------
      if(lobj90 .lt. ldo90) go to 2100

C************ Going UP,  Tensile Direction **************************
C     Are we on the monotonic curve?
 1100 if(nptt90 .eq. 0) go to 1500
c     No. Has a exceedence  occurred?
 1110 if(lobj90 .gt. tlimL90(nptt90)) go to 1120

c     Is this Same as previous amplitude? e.g.: Constant ampl. test
 1130 if(lobj90 .eq. tlimL90(nptt90) .and. nptt90 .ne. 1) go to 1135
      go to 1350

c     Is a loop being closed without the monotonic?
 1120 if(nptt90 .ne. 1) go to 1400

c     Is closure in connection  with the monotonic?
 1150 if(nptc90 .eq. 2) go to 1170

c     Check for impossible
 1160 if(nptc90 .eq. 1) go to 1165
      idump=1160
      go to 9999

C     Same load level as previous level.
 1135 continue
Css      dam90=SMITH( (lobj-ldo),tlsts(nptt90),clsts(nptc90),tlstr(nptt90),
Css     &            clstr(nptc90),totdam,nrev)
C     lobj90 and ldo90 exist. They are the strain values. Use them
      dld90=lobj90-ldo90 !this is total Strain range in pdStressStrain.f
      if(dld90.le.0.)then
        write(6,*)"#Error: 1135:dld90,lobj90,ldo90 = ",
     &            dld90,lobj90,ldo90,nrev
      endif
C     Get the stress,strain, damage stuff from the previous rev's storage:
      xStsmax=tlsts90(nptt90)
      xStrmax=tlstr90(nptt90)
      xStsmin=clsts90(nptc90)
      xStrmin=clstr90(nptc90)
      dam90=tldam90(nptt90)
C     get the stress-strain path of this half cycle
      xStsAmpl=(xStsmax-xStsmin)/2.0
      call getloop(xStsAmpl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 1137 i=1,nreturn
        xpoint=eo90+tstrain(i)*2.0
        ypoint=so90+tstress(i)*2.0
        if(xpoint .le.StrOut .or. ypoint .le.StsOut)goto 1137
        write(6,1136)xpoint,ypoint,dtime(nact)
 1136   format("#plotloops ",f8.5,1x,f6.1,1x,f13.3)
 1137 continue
      StrOut=xpoint
      StsOut=ypoint

      nptc90=nptc90-1
      dtotdam90=dtotdam90+dble(dam90)
c     The closed loop is erased and the point of rev. is already there
C     in the push-down list.
      eobj90=tlstr90(nptt90)
      sobj90=tlsts90(nptt90)
      tldam90(nptt90)   =dam90
      tltotCrk90(nptt90)=dtotdam90
      itlnrev90(nptt90) =nrev
      go to 5000  ! Ramp is done go to 00 code


C     An unmatched 1/2 cycle is being forgotten and a return to the 
c     monotonic curve is occurring. Count the unmatched 1/2 cycle.
C     The monotonic is usually not counted for damage.
 1165 continue
Css      dam90=SMITH ( (tlimL90(1)-climL90(1)), tlsts90(1), clsts90(1),
Css     &         tlstr90(1),clstr90(1),totdam90,nrev)
      dld90=tlimL90(1)-climL90(1) ! Strain Rg of unmatched 1/2 cyc
      if(dld90.eq.0.)then
        write(6,*)"#Error: 1165:dld90=0 tlimL90(1),climL90(1) = ",
     &            tlimL90(1),climL90(1),nrev
      endif
      xStsmax=tlsts90(1)
      xStrmax=tlstr90(1)
      xStsmin=clsts90(1)
      xStrmin=clstr90(1)
      Smea= xStsmax * (xStrmax-xStrmin)/2.0
      call getSwat2Life(Smea,Swatlife,iexit)
      write(6,*)"#getSwatl:",Smea,Swatlife
C     Infinite life will be a -ve number
      if(Swatlife .gt.0)then
C       this is just a 1/2 cycle, thus multiply by 2
        dam90=2.0/Swatlife
        dtotdam90=dtotdam90+dble(dam90)
      else
C       Inf. life,  damage is zero.
        dam90=0.
      endif
C     get or get the finish of this half cycle
      xStsAmpl=(xStsmax-xStsmin)/2.0
      call getloop(xStsAmpl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output. Proceed from there.
      do 1167 i=1,nreturn
        xpoint=eo90+tstrain(i)*2.0
        ypoint=so90+tstress(i)*2.0
        if(xpoint .le.StrOut .or. ypoint .le.StsOut)goto 1167
        write(6,1136)xpoint,ypoint,dtime(nact)
 1167 continue
      StrOut=xpoint
      StsOut=ypoint
      nptc90=0
      nptt90=0
      go to 1500  !go back on monotonic. P.D. lists are empty.

C     A loop is being closed and a return to the tensile monotonic 
c     curve is occurring.  Damage is the same as the other half of
c     the cycle being closed.
 1170 continue
Css      dam90=SMITH( (tlimL90(1)-climL90(2)),
Css     &    tlsts90(1),clsts90(2),
Css     &    tlstr90(1),clstr90(2), totdam,nrev)
      xStsmax=tlsts90(1)
      xStrmax=tlstr90(1)
      xStsmin=clsts90(2)
      xStrmin=clstr90(2)
      dld90=tlimL90(1)-climL90(2) !"load" is same as strain in this case
      if(dld90.eq.0.)then
        write(6,*)"#Error: 1170:dld90=0 tlimL90(1),climL90(2)= ",
     &            tlimL90(1),climL90(2),nrev
      endif
      Smea= xStsmax * (xStrmax-xStrmin)/2.0
      call getSwat2Life(Smea,Swatlife,iexit)
      write(6,*)"#getSwatl:",Smea,Swatlife,dtime(nact)
C     Infinite life will be a -ve number
      if(Swatlife .gt.0)then
C       this is just a 1/2 cycle, thus multiply by 2
        dam90=2.0/Swatlife
        dtotdam90=dtotdam90+dble(dam90)
      else
C       Inf. life,  damage is zero.
        dam90=0.
      endif
      xStsAmpl=(xStsmax-xStsmin)/2.0
      call getloop(xStsAmpl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 1172 i=1,nreturn
        xpoint=eo90+tstrain(i)*2.0
        ypoint=so90+tstress(i)*2.0
        if(xpoint .le.StrOut .or. ypoint .le.StsOut)goto 1172
        write(6,1136)xpoint,ypoint,dtime(nact)
 1172 continue
      StrOut=xpoint
      StsOut=ypoint

c     Eliminate the closed loop.
      nptc90=0
      nptt90=0
      go to 1500  ! go back to monotonic

C     A new entry in the P.D. list is being made. We are at a rev. point
 1350 continue
C     compute the co-ordinates
      dld90=abs(lobj90-climL90(nptc90))
Css      de90=FLD(dld,nrev) !old code for ref
Css      ds90=DET(de,nrev)
Css      eobj90=eo90+de90
Css      sobj90=so90+ds90
Css      dam90=SMITH( dld90,sobj90,so90,eobj90,eo90,totdam90,nrev)
      if(dld90.eq.0.)then
        write(6,*)"#Error: 1350:dld90=0 lobj90,climL90(nptc90)= ",
     &            lobj90,climL90(nptc90),nrev
      endif
      dld90ampl=dld90/2.0
      call getLoad2StressStrain(dld90ampl,ds90ampl,de90ampl,iret)
      eobj90=eo90+(de90ampl*2.0)
      sobj90=so90+(ds90ampl*2.0)
      Smea= sobj90 * (eobj90-eo90)/2.0
      call getSwat2Life(Smea,Swatlife,iexit)
      write(6,*)"#getSwatl:",Smea,Swatlife,dtime(nact)
C     Infinite life will be a -ve number
      if(Swatlife .gt.0)then
C       this is just a 1/2 cycle, thus multiply by 2
        dam90=2.0/Swatlife
        dtotdam90=dtotdam90+dble(dam90)
      else
C       Inf. life,  damage is zero.
        dam90=0.
      endif
      call getloop(ds90ampl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 1352 i=1,nreturn
        xpoint=eo90+tstrain(i)*2.0
        ypoint=so90+tstress(i)*2.0
        if(xpoint .le.StrOut .or. ypoint .le.StsOut)goto 1352
        write(6,1136)xpoint,ypoint,dtime(nact)
 1352 continue
      StrOut=xpoint
      StsOut=ypoint

c     Enter the rev in the P.D. list
      nptt90=nptt90+1
      tlimL90(nptt90)=lobj90
      tlstr90(nptt90)=eobj90
      tlsts90(nptt90)=sobj90
      tldam90(nptt90)   =dam90
      tltotCrk90(nptt90)=dtotdam90
      itlnrev90(nptt90) =nrev
      go to 5000  ! Ramp is done go to 00 code

C     Deformation is occurring on the monotonic curve again.
 1500 continue
      dld90=abs(lobj90)
      if(dld90.eq.0.)then
        write(6,*)"#Error: 1500:dld90=0 lobj90= ",
     &            lobj90,nrev
      endif
      dld90ampl=dld90
      call getLoad2StressStrain(dld90ampl,ds90ampl,de90ampl,iret)
      eobj90=de90ampl
      sobj90=ds90ampl
      eo90=0.
      so90=0.
c     Subtract damage of previous use of monotonic curve (=0 in cyc )
C      dtotdam90=dtotdam90-dble(tldam90(1) )
C     Decision Jan31 2015: damage is not counted on the monotonic curve
      dam90=0.

C     Fill in the stress-strain path extension.
      call getloop(ds90ampl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 1502 i=1,nreturn
        xpoint=eo90+tstrain(i) ! on monot. we do not multiply by 2
        ypoint=so90+tstress(i)
        if(xpoint .le.StrOut .or. ypoint .le.StsOut)goto 1502
        write(6,1136)xpoint,ypoint,dtime(nact)
 1502 continue
      StrOut=eobj90 
      StsOut=sobj90 
      write(6,*)"#plotloops #StrOut, StsOut= ",StrOut,StsOut

c     Set the appropriate push-down list arrays for rev on monotonic.
C     We are on the tensile monotonic.
      nptt90=1
      nptc90=1
      climL90(1)=-abs(lobj90)
      tlimL90(1)=abs(lobj90)
      tlstr90(1) =abs(eobj90)
      clstr90(1) =-abs(eobj90)
      tlsts90(1) =abs(sobj90)
      clsts90(1) =-abs(sobj90)

      tldam90(1)   =dam90 ! monot. damage is zero
      tltotCrk90(1)=dtotdam90
      itlnrev90(1) =nrev
      cldam90(1)   =dam90
      cltotCrk90(1)=dtotdam90
      iclnrev90(1) =nrev
      go to 5000 ! half cycle is done, go to 00 direction crack code



C     A loop is being closed, count this remaining halfcycle of the
c     loop and then eliminate this closing loop. Then continue on.
 1400 continue
      dld90=tlimL90(nptt90)-climL90(nptc90)
Css      dam90=SMITH( (tlimL90(nptt90)-climL90(nptc90)), 
Css     &      tlsts90(nptt90),clsts90(nptc90),
Css     &      tlstr90(nptt90),clstr90(nptc90),totdam90,nrev)
      if(dld90.eq.0.)then
       write(6,*)"#Error:1400:dld90=0 tlimL90(nptt90),climL90(nptc90)=",
     &            tlimL90(nptt90),climL90(nptc90),nrev
      endif
      xStsmax= tlsts90(nptt90)
      xStrmax= tlstr90(nptt90)
      xStsmin= clsts90(nptc90)
      xStrmin= clstr90(nptc90)
      Smea= xStsmax * (xStrmax-xStrmin)/2.0
      call getSwat2Life(Smea,Swatlife,iexit)
      write(6,*)"#getSwatl:",Smea,Swatlife,dtime(nact)
C     Infinite life will be a -ve number
      if(Swatlife .gt.0)then
C       this is just a 1/2 cycle, thus multiply by 2
        dam90=2.0/Swatlife
        dtotdam90=dtotdam90+dble(dam90)
      else
C       Inf. life,  damage is zero.
        dam90=0.
      endif

      xStsAmpl=(xStsmax-xStsmin)/2.0
      call getloop(xStsAmpl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 1402 i=1,nreturn
        xpoint=eo90+tstrain(i)*2.0
        ypoint=so90+tstress(i)*2.0
        if(xpoint .le.StrOut .or. ypoint .le.StsOut)goto 1402
        write(6,1136)xpoint,ypoint,dtime(nact)
 1402 continue
      StrOut=xpoint
      StsOut=ypoint

      nptc90=nptc90-1
C     Subtract the old half cycle damage of the Re-incurred
c     stress-strain path.
      dtotdam90=dtotdam90-dble(tldam90(nptt90) )
      nptt90=nptt90-1

c     Reset the local stress-strain origin to begining of older loop.
      eo90=clstr90(nptc90)
      so90=clsts90(nptc90)
      ldo90=climL90(nptc90)
      go to 1100 ! contiue onwards in tens. direction



C*****  Going Down, Compressive Direction ****************************

C     Are we on the monotonic curve?
 2100 if(nptc90 .eq. 0) go to 2500

c     No. Has a exceedence  occurred?
 2110 if(lobj90 .lt. climL90(nptc90)) go to 2120

c     Is this Same as previous amplitude?
 2130 if(lobj90 .eq. climL90(nptc90) .and. nptc90 .ne. 1) go to 2135
      go to 2350

c     Is a loop being closed without the monotonic?
 2120 if(nptc90 .ne. 1) go to 2400

c     Is closure in connection  with the monotonic?
 2150 if(nptt90 .eq. 2) go to 2170

c     Check for impossible
 2160 if(nptt90 .eq. 1) go to 2165
      idump=2160
      go to 9999

C     Same load level as in previous cycle's level.
 2135 continue
Css      dam90=SMITH( (ldo90-lobj90),tlsts90(nptt90),clsts90(nptc90),
Css     &      tlstr90(nptt90),clstr90(nptc90),totdam90,nrev)
C     lobj90 and ldo90 exist. They are the strain values. Use them
      dld90=ldo90-lobj90 !this is total Strain range in pdStressStrain.f
      if(dld90.eq.0.)then
       write(6,*)"#Error:2135:dld90=0 ldo90,lobj90=",
     &            ldo90,lobj90,nrev
      endif
C     Get the stress,strain, damage stuff from the previous rev's storage:
      xStsmax=tlsts90(nptt90)
      xStrmax=tlstr90(nptt90)
      xStsmin=clsts90(nptc90)
      xStrmin=clstr90(nptc90)
      dam90=  cldam90(nptc90)
C     get the stress-strain path of this half cycle
      xStsAmpl=(xStsmax-xStsmin)/2.0
      call getloop(xStsAmpl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 2137 i=1,nreturn
        xpoint=eo90-(tstrain(i)*2.0)
        ypoint=so90-(tstress(i)*2.0)
        if(xpoint .ge.StrOut .or. ypoint .ge.StsOut)goto 2137
        write(6,2136)xpoint,ypoint,dtime(nact)
 2136   format("#plotloops ",f8.5,1x,f6.1,1x,f13.3)
 2137 continue
      StrOut=xpoint
      StsOut=ypoint

      dtotdam90=dtotdam90+dble(dam90)
      nptt90=nptt90-1
c     The closed loop is erased and the point of rev. is already there
      eobj90=clstr90(nptc90)
      sobj90=clsts90(nptc90)
      cldam90(nptc90)   =dam90
      cltotCrk90(nptc90)=dtotdam90
      iclnrev90(nptc90) =nrev
      go to 5000  ! Ramp is done go to 00 code

C     An unmatched 1/2 cycle is being forgotten and a return to the
c     monotonic curve is occurring. Count the unmatched 1/2 cycle.
C     The monotonic is not counted for damage.
 2165 continue
Css      dam90=SMITH( (tlimL90(1)-climL90(1)),tlsts90(1),clsts90(1),
Css     &        tlstr90(1),clstr90(1),totdam90,nrev)
      dld90=tlimL90(1)-climL90(1)
      if(dld90.eq.0.)then
       write(6,*)"#Error:2165:dld90=0 tlimL90(1),climL90(1)=",
     &            tlimL90(1),climL90(1),nrev
      endif
      xStsmax=tlsts90(1)
      xStrmax=tlstr90(1)
      xStsmin=clsts90(1)
      xStrmin=clstr90(1)
      Smea= xStsmax * (xStrmax-xStrmin)/2.0
      write(6,*)"#dbug:xStsmax,xStrmax,xStrmin",xStsmax,xStrmax,xStrmin
      call getSwat2Life(Smea,Swatlife,iexit)
      write(6,*)"#getSwatl:",Smea,Swatlife
C     Infinite life will be a -ve number
      if(Swatlife .gt.0)then
C       this is just a 1/2 cycle, thus multiply by 2
        dam90=2.0/Swatlife
        dtotdam90=dtotdam90+dble(dam90)
      else
C       Inf. life,  damage is zero.
        dam90=0.
      endif
C     get or get the finish of this half cycle
      xStsAmpl=(xStsmax-xStsmin)/2.0
      call getloop(xStsAmpl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output. Proceed from there.
      do 2167 i=1,nreturn
        xpoint=eo90-(tstrain(i)*2.0)
        ypoint=so90-(tstress(i)*2.0)
        if(xpoint .ge.StrOut .or. ypoint .ge.StsOut)goto 2167
        write(6,2136)xpoint,ypoint,dtime(nact)
 2167 continue
      StrOut=xpoint
      StsOut=ypoint
      nptc90=0

      nptc90=0
      nptt90=0
      go to 2500 !go back on monotonic. P.D. lists are empty.

C     A loop is being closed and a return to the monotonic
c     curve is occurring.  Damage is the same as the other half of
c     the cycle being closed.
 2170 continue
Css      dam90=SMITH( (tlimL90(2)-climL90(1)),tlsts90(2),
Css     &       clsts90(1),tlstr90(2),clstr90(1),totdam90,nrev)
      dld90=tlimL90(2)-climL90(1)
      if(dld90.eq.0.)then
       write(6,*)"#Error:2170:dld90=0 tlimL90(2),climL90(1)=",
     &            tlimL90(2),climL90(1),nrev
      endif
      dld90=tlimL90(2)-climL90(1) !"load" is same as strain in this case
      if(dld90.eq.0.)then
        write(6,*)"#Error: 1170:dld90=0 tlimL90(1),climL90(2)= ",
     &            tlimL90(1),climL90(2),nrev
      endif
      xStsmax=tlsts90(2)
      xStrmax=tlstr90(2)
      xStsmin=clsts90(1)
      xStrmin=clstr90(1)
      Smea= xStsmax * (xStrmax-xStrmin)/2.0
      write(6,*)"#dbug:xStsmax,xStrmax,xStrmin",xStsmax,xStrmax,xStrmin
      call getSwat2Life(Smea,Swatlife,iexit)
      write(6,*)"#getSwatl:",Smea,Swatlife,dtime(nact)
C     Infinite life will be a -ve number
      if(Swatlife .gt.0)then
C       this is just a 1/2 cycle, thus multiply by 2
        dam90=2.0/Swatlife
        dtotdam90=dtotdam90+dble(dam90)
      else
C       Inf. life,  damage is zero.
        dam90=0.
      endif
      xStsAmpl=(xStsmax-xStsmin)/2.0
      call getloop(xStsAmpl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 2172 i=1,nreturn
        xpoint=eo90-(tstrain(i)*2.0)
        ypoint=so90-(tstress(i)*2.0)
        if(xpoint .ge.StrOut .or. ypoint .ge.StsOut)goto 2172
        write(6,2136)xpoint,ypoint,dtime(nact)
 2172 continue
      StrOut=xpoint
      StsOut=ypoint

C     Eliminate the closed loop.
      nptc90=0
      nptt90=0
      go to 2500 !go back to monotonic


C     A new entry in the P.D. list is being made. We are at a rev. point
 2350 continue
      dld90=abs(lobj90-tlimL90(nptt90))
Css      de90=FLD(dld90,nrev) !old code for ref
Css      ds90=DET(de90,nrev)
Css      eobj90=eo90-de90
Css      sobj90=so90-ds90
Css      dam90=SMITH(dld90,so90,sobj90,eo90,eobj90,totdam90,nrev)
      if(dld90.eq.0.)then
       write(6,*)"#Error:2350:dld90=0 lobj90,tlimL90(nptt90)=",
     &            lobj90,tlimL90(nptt90),nrev
      endif
      dld90ampl=dld90/2.0
      call getLoad2StressStrain(dld90ampl,ds90ampl,de90ampl,iret)
      eobj90=eo90-(de90ampl*2.0)
      sobj90=so90-(ds90ampl*2.0)
C     Fix bug here Smax was not  sobj90   March 5 2015
      Smea= so90 * de90ampl
      write(6,*)"#dbug:so90,de90ampl: ",so90,de90ampl
      call getSwat2Life(Smea,Swatlife,iexit)
      write(6,*)"#getSwatl:",Smea,Swatlife,dtime(nact)
C     Infinite life will be a -ve number
      if(Swatlife .gt.0.)then
C       this is just a 1/2 cycle, thus multiply by 2
        dam90=2.0/Swatlife 
        dtotdam90=dtotdam90+dble(dam90) 
      else
C       Inf. life,  damage is zero.
        dam90=0.
      endif 
      call getloop(ds90ampl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 2352 i=1,nreturn
        xpoint=eo90-(tstrain(i)*2.0)
        ypoint=so90-(tstress(i)*2.0)
        if(xpoint .ge.StrOut .or. ypoint .ge.StsOut)goto 2352
        write(6,2136)xpoint,ypoint,dtime(nact)
 2352 continue
      StrOut=xpoint
      StsOut=ypoint

C     Enter the rev in the P.D. list
      nptc90=nptc90+1
      climL90(nptc90)=lobj90
      clstr90(nptc90)=eobj90
      clsts90(nptc90)=sobj90
      cldam90(nptc90)   =dam90
      cltotCrk90(nptc90)=dtotdam90
      iclnrev90(nptc90) =nrev
      go to 5000  ! Ramp is done go to 00 code


C     Deformation is occurring on the monotonic curve again. --------
 2500 continue
      dld90=abs(lobj90)
      if(dld90.eq.0.)then
        write(6,*)"#Error: 2500:dld90=0 lobj90= ",
     &            lobj90,nrev
      endif
      dld90ampl=dld90
      call getLoad2StressStrain(dld90ampl,ds90ampl,de90ampl,iret)
      eobj90= - de90ampl
      sobj90= - ds90ampl
      eo90=0.
      so90=0.
c     Subtract damage of previous use of monotonic curve (=0 in cyc )
C      dtotdam90=dtotdam90 - dble(tldam90(1) )
C     Decision Jan31 2015: damage is not counted on the monotonic curve
      dam90=0.

C     Fill in the stress-strain path extension.
      call getloop(ds90ampl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 2502 i=1,nreturn
        xpoint=eo90-(tstrain(i)) ! on monot. we do not multiply by 2
        ypoint=so90-(tstress(i))
        if(xpoint .ge.StrOut .or. ypoint .ge.StsOut)goto 2502
        write(6,2136)xpoint,ypoint,dtime(nact)
 2502 continue
      StrOut=eobj90
      StsOut=sobj90

c     Set the appropriate push-down list arrays for rev on monotonic.
      nptt90=1
      nptc90=1
      climL90(1)=-abs(lobj90)
      tlimL90(1)=abs(lobj90)
      tlstr90(1) =abs(eobj90)
      clstr90(1) =-abs(eobj90)
      tlsts90(1) =abs(sobj90)
      clsts90(1) =-abs(sobj90)

      tldam90(1)   =dam90
      tltotCrk90(1)=dtotdam90
      itlnrev90(1) =nrev
      cldam90(1)   =dam90
      cltotCrk90(1)=dtotdam90
      iclnrev90(1) =nrev
      go to 5000    ! half cycle is done, go to  00 code

C     A loop is being closed, count the remaining half of the
c     loop and then eliminate the loop.
 2400 continue
Css      dam90=SMITH( (tlimL90(nptt90)-climL90(nptc90)), 
Css     &       tlsts90(nptt90),clsts90(nptc90),
Css     &       tlstr90(nptt90),clstr90(nptc90),totdam90,nrev)
      dld90=tlimL90(nptt90)-climL90(nptc90)
      if(dld90.eq.0.)then
       write(6,*)"#Error:2400:dld90=0 tlimL90(nptt90),climL90(nptc90)=",
     &            tlimL90(nptt90)-climL90(nptc90),nrev
      endif
      xStsmax= tlsts90(nptt90)
      xStrmax= tlstr90(nptt90)
      xStsmin= clsts90(nptc90)
      xStrmin= clstr90(nptc90)
      Smea= xStsmax * (xStrmax-xStrmin)/2.0
      write(6,*)"#dbug:xStsmax,xStrmax,xStrmin",xStsmax,xStrmax,xStrmin
      call getSwat2Life(Smea,Swatlife,iexit)
      write(6,*)"#getSwatl:",Smea,Swatlife,dtime(nact)
C     Infinite life will be a -ve number
      if(Swatlife .gt.0)then
C       this is just a 1/2 cycle, thus multiply by 2
        dam90=2.0/Swatlife
        dtotdam90=dtotdam90+dble(dam90)
      else
C       Inf. life,  damage is zero.
        dam90=0.
      endif

      xStsAmpl=(xStsmax-xStsmin)/2.0
      call getloop(xStsAmpl,tstrain,tstress,nreturn,iret)
C     StrOut and StsOut are the points of previous loop output.
C     Proceed from there.
      do 2402 i=1,nreturn
        xpoint=eo90-(tstrain(i)*2.0)
        ypoint=so90-(tstress(i)*2.0)
        if(xpoint .ge.StrOut .or. ypoint .ge.StsOut)goto 2402
        write(6,2136)xpoint,ypoint,dtime(nact)
 2402 continue
      StrOut=xpoint
      StsOut=ypoint

      nptt90=nptt90-1
C     Subtract the old half cycle damage of the Re-incurred
c     stress-strain path.
      dtotdam90=dtotdam90-dble(cldam90(nptc90) )
      nptc90=nptc90-1

c     Reset the local stress-strain origin to begining of older loop.
      so90=tlsts90(nptt90)
      eo90=tlstr90(nptt90)
      ldo90=tlimL90(nptt90)
      go to 2100     !    continue the ramp in comp. direction.



 5000 continue
      go to 3000    !go back to begin of overall program cycle loop

 9000 continue
C   future binary output file code:
CC      write out the -ve number of the last rec.  that was written.
CC      All the other variables in this rec are endof test values.
C       ndummy=-nrecord  !make -ve to make it unique.
C       xdummy=0.
C       write(60,rec=1 )ndummy,xdummy,nblk,nact,
C     &            lobj90,xMm90,xMb90,
C     &            stsMembrane,stsBending
C      close(unit=60)

       write(0,9050)nrev,totdam90,nblk,nact,nrecord
       write(6,9050)nrev,totdam90,nblk,nact,nrecord
 9050  format("#Last: nrev= ",i10," a= ",e14.7,
     &        "   nblk= ",i10," nact= ",i10," nrecord= ",i10/)
      stop

 9999 write(6,9998) idump
 9998 format("# *** Error: near statement label no. ",i5/)
      stop
      end           !  END of MAINLINE



C==============================================================
      SUBROUTINE WR(UN,INP1)
C S/R TO WRITE A LINE BUT CUT THE TRAILING BLANKS
      CHARACTER*1 INP1(300)
      INTEGER UN
      LONG=300
C FIND LAST NON BLANK CHAR
      DO 10 I=LONG,1,-1
        IF(INP1(I).EQ.' ')GO TO 10
C       NO? THIS IS IT
        N=I
        GO TO 20
   10 CONTINUE

C COMPLETELY BLANK LINE
      WRITE(UN,11)
   11 FORMAT(' ')
      RETURN

   20 WRITE(UN,24)(INP1(J),J=1,N)
   24 FORMAT(300A1)
      RETURN
      END

C==========================================================================

      SUBROUTINE readStrainLifeStress(fname,istraincontrol,iret)
      SAVE
C     fname          Contains material fitted StrainLifeStress filename
C     istrainControl  =0  The input history will be Elastic FEA type stress
C                     =1  The input history will be strain.
C                         (i.e. no Neuber correction for plasticity)
C     iret= 0  No error,      =1 error encountered.

C Linux Compile if converted to a mainline:  
C   gfortran -g -w -fbounds-check readStrainLifeStress.f  -o readStrainLifeStress
C   
C    Opens the filename in fname and creates two common areas: 
C          /MaterialA/ contains the data from file in"fname" plus a few more
C                     for low speed interpolations such as in saefcalc2.f
C         !!!     Common area is different from  saefcalc2.f version.

C---------------------------------------------------------------------------
C  Copyright (C) 2015 A.Conle added features.

C Fork as subroutine from : saefcalc2.f  #version= 2.1 
C  Copyright (C) 2004 SAE Fatigue Design and Evaluation Committee
C  This program is free software; you can redistribute it and/or
C  modify it under the terms of the GNU General Public License as
C  published by the Free Software Foundation; either version 2 of the
C  license, or (at your option) any later version.
C
C  This program is distributed in the hope tha it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General PUblic License for more details.
C
C  You should have received a copy of the GNU General PUblic License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA
C  Try also their web site: http://www.gnu.org/copyleft/gpl.html
C---------------------------------------------------------------------------


C Assumptions and other items of interest:
C  Interpolation, on curves such as stress-life or strain-life is done in
C   log-log space.  Generally anything that plots "better" on log-log
C   co-ords is also interpolated there.

C The material input file is assumed to be a "fitted" file i.e.: the computer
C  will assume that a line could be drawn through the successive input data
C  points and that this line would make a reasonable curve, reasonalble 
C  enough to allow one to compute values between the points. 
C  It is expected that the tag 
C                #DataType= fitted
C  will appear in the input stream.  If not the program should error out.
C  If you are changing this routine please consider that sooner or later
C  someone will send a different type of file by mistake.
C  It is also expected that the tag:
C                #FileType= strain_life
C  will appear in the input stream.  But, in the intrests of possibly 
C  allowing other types of data files (e.g.: plastic materials), I'm not
C  certain this rule will always be enforced.

C Other important identifiers expected in the input stream along with
C defaults assumed and others allowed:
C
C   Tag        Example      Default  Others_possible......................
C  #name=  SAE1010          Unknown  anything that starts with a letter

C  #E=     29000.           None(error out)

C  #Stress_Units= MPA       MPA      ksi,   psi

C  #Strain_Units= strain    strain   microstrain

C  #Life_Units= Reversals   Reversals  Cycles

C  #Su=  119.9              none
C  #Sy=   92.0              none
C  #BHN= 243                none

C Future feature?:  Allow material file tag
C #Sort= life (default)   or  strain,  stress
C and sort the digital curves accordingly.  Life sort should have biggest
C first, Stress & strain sorts should have smallest first.
C This would accomodate stress life curves that have shorter lives at low
C stress levels. (i.e.: Cup towards left, or Cup right etc)
C Cup or Cap type data is not allowed.

C The first part of the program was adapted from saedigcurve.f.  It reads
C in the information from the material file, and decodes the input arguments.


C-----------------------------------------------------------------------------
C Eg. of SAE standard form fatigue data file:

C   ______ first column in file
C  |
C  v
C
C  # SAE Exchange File Format.   
C  # data collected from ASTM E606 axial fatigue test data.
C  # Note: The data below is not a real file
C   
C  #DataType= fitted
C  
C  #NAME= SAE1045
C  #NAME= SAE350X
C  #NAME= SAE050X
C  #Ford= 34   
C  #Stress_Units= KSI
C  #Strain_Units= strain
C  #Life_Units= reversals
C  #E=  30000.
C  #Su= 89.
C  #Sy= 50.
C  #%RA= 85.
C  #BHN= 325
C  #WebPage= http://fde.uwaterloo.ca/Fde/Materials/Steels/ASTM-A588C/g40.21-50A_non_os.html
C   
C  # Total Strain   2Nf  Stress  Mean   Plastic Strain   Initial
C  #    Amp               Amp   Stress      Amp        Elastic Mod.
C  0.0125          180   279.   .0       0.0030        30100.    #Fitted_point
C  0.0095          490   253.   .0       0.0011        29400.    #Fitted_point
C  0.0090          950   229.   .0       0.0007        29800.    #Fitted_point
C  0.0075         2260   220.   .0       0.0002        30050.    #Fitted_point
C  0.0050        38000.  149.   .0       0.0           29900.    #Fitted_point
C  0.0040       770000   119.   .0       0             30700.    #Fitted_point
C  #
C-----------------------------------------------------------------------------


      CHARACTER*1   INP1(80)
      CHARACTER*5   INP5(16),JNP5(16)
      CHARACTER*10  INP10(8),JNP10(56)
      character*11  clifeout(10)

      CHARACTER*80  JNP80(7),INP80, INP80temp 
      EQUIVALENCE  (INP1(1),INP5(1),INP10(1),INP80),
     &             (JNP5(1),JNP10(1),JNP80(1))

      character*300  inp300,jnp300,Cwebpage
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)

      character*10 name1, name2, stressunits
      character*10 strainunits, lifeunits, sorttype
      character*30 names30(10), firstfield, ctail
      character*30 cfiletype, Cdatatype
      character*80 argv, fname

      integer uno
      integer*4 iargc, argc

      integer istrainControl !  =0  or  =1


C     Save the fitted Strain-Life-Stress  digital curves here:
C     Check dimensions in the various subroutines too    !!!!!!!!!!!!
      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain,Emod
C     Use these to transform to ALOG  to speed up interpolations
      real logStrainAmp(250),logLifecycles(250),logStressAmp(250),
     &     logSwtAmp(250),logSelasticAmp(250)
      Common/MaterialA/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 logStrainAmp,logLifecycles,logStressAmp,
     &                 logSwtAmp,logSelasticAmp,
     &                 ndata,FractureStress,FractureStrain,Emod


C     tstrain & tstress are used only for loop plotting and passing
c     stuff back from s/r getloop
C      real tstrain(250),tstress(250)


C     Future Code.  Not used now.
C     Storage for higher speed interp
C      real StressAmpX(1500), StrainAmpX(1500), SnominalAmpX(1500),
C     &     LifecyclesX(1500),SwtAmpX(1500)
C      integer*4 nmatdata,maxmatdata
CC      real FractureStress,FractureStrain,Emod
C      logical debugMat
C      common/XCALC4/ StressAmpX,StrainAmpX,SnominalAmpX,LifecyclesX,
C     &      snominalInterval,nmatdata,maxmatdata,debugMat
C
C      maxmatdata= 1500


Cdidntwork to s/rs      PARAMETER ( IDIM=250 )
      idimSS=250
C         = max dimension of digital curve stuff



      iret=0
      write(6,50)
   50 format("#version= 1.0  readStrainLifeStress.f starts...")


      XMPAS=6.894759
C     Set some default values in case user forgets
      stressunits=" "
      strainunits=" "
      lifeunits="reversals"
      EMOD=0.
      Cdatatype=" "
      cfiletype=" "
      Sult=0.
      Syield=0.
      percentRA=0.
      Cwebpage=" "
      FractureStress=0.
      FractureStrain=0.
      nbrinell=0
      sorttype="life"


C        call getarg(jvect,fname)
C        write(6,*)"#Opening material file= ",fname," as unit 10"
C        write(6,*)"#matfile= ",fname," as unit 10"
CC          #matfile=   is the string used to get filename in 
CC                      makeInitReport script after this runs.
        open(unit=10,file=fname)


  610   continue



C---------------------------- Read in Fitted material file from stdin-------------
  800 continue
c     Loop back to here for next input line.

C     Input lines may either be 1. data lines,  2. #Comment lines or 3. blank
C        It will be hard to distinguish the real comment from the junk comment.
C        Leave it up to the user to edit the comment section.

      read(10,"(a300)",end=980)inp300
      ninput=ninput+1
Cdebug      write(0,*)" read input line ",ninput
   
C     Check for blank line
      if(inp300.eq." ")then
C        write(6,"(a1)")" "
        go to 800
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 805 i=1,300
           if(inpone(i).eq." ") go to 805
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 803 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  803        continue
Cdebug             write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 880
           else
C           The first non blank is not a #
             go to 810
           endif
  805    continue

      
  810    continue
C        1st non blank is not a #.  Check if its a number.
         if(inpone(loc).ne."+" .and.
     &      inpone(loc).ne."-" .and.
     &      inpone(loc).ne."." .and.
     &      inpone(loc).ne."0" .and.
     &      inpone(loc).ne."1" .and.
     &      inpone(loc).ne."2" .and.
     &      inpone(loc).ne."3" .and.
     &      inpone(loc).ne."4" .and.
     &      inpone(loc).ne."5" .and.
     &      inpone(loc).ne."6" .and.
     &      inpone(loc).ne."7" .and.
     &      inpone(loc).ne."8" .and.
     &      inpone(loc).ne."9"  )then
C           It must be a letter, not a number. Time to bomb out.
            write(0,815)ninput
            write(6,815)ninput
  815       format(" ERROR pdStressStrain: input line no. ",I5,
     &      " not a # and not a number. Edit fitted input file")
            stop
         endif


C        Ok, it should be like this one:
C        0.0125          180   279.   .0       0.0030        30100.    #specimen comment
C        #Total Strain   2Nf  Stress  Mean   Plastic Strain   Initial
C        ##   Amp               Amp   Stress      Amp        Elastic Mod.
C        Note that the trailing comment field "#specimen comment" may or may not be there.

C        We should be able to read the first 5 fields as real numbers.  Since 2Nf may have
C        more than 7 signif. digits, it needs to be read in as a double precision.  It may
C        be an integer.  The value can be changed later to whatever the local format is.

C      Since this is all "Fitted" data expressing a curve, we only need the first 3 items.

         ndata=ndata+1
         if(ndata.gt.idimSS)then
           write(0,816)
           write(6,816)
  816      format("# Error too many data points in strain-life file:",
     &            i5," recompile s/r readStrainLifeStress: See idimSS ")
           stop
         endif
         read(inp300,*)StrainAmp(ndata), Lifecycles(ndata), 
     &                 StressAmp(ndata)
Cdebug         write(6,*)StrainAmp(ndata), Longlife(ndata), 
Cdebug      &            StressAmp(ndata)

C          If we have not crashed by here, the data from the line has been read in. 
C          Now figure out if there is comment at the end of line.
C          In saefcalc2 we don't really need to do this, but the code is here for
C          some possible future use. ? :)
C          Brute force it. Hunt for a #

           do 820 i=1,300
             if(inpone(i).eq."#")then
C               we found it in col i
Cdebug                write(0,*)"found # trailer in data line"
                loc=i
                go to 822
             endif
  820      continue
C          If we got to here, the trailing field is empty
           go to 855

  822      continue
C          Now count backwards and see where the last char is 
           do 830 i=1,300
             j=300-(i-1)
             if(inpone(j).ne." ")then
C              found last char
               lastloc=j
               go to 832
             endif
  830      continue

  832      continue
C          It is possible that the #field is an html tag, so it may be real long
C          Move the stuff to begin of a new field for decoding
           j=0
           do 840 i=loc,lastloc
             j=j+1
             jnpone(j)=inpone(i)
  840      continue
Cdebug           write(0,*)"Trailer # on data: ",(jnpone(i),i=1,j)

C          Check to see if any special tags are in this field
C          If the field is blank, then what?
           read(jnp300,*)ctail
C           if(ctail.eq."#Runout" .or.
C     &        ctail.eq."#RUNOUT" .or.
C     &        ctail.eq."#runout" )then
C            The data point was a runout
C            In the local company's format, runouts are -ve nos.
C            Change it to whatever you like !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C             longlife=-longlife
C           endif
C          #knifeEdge and #outsideGage  are also known tags
C           but we dont do anything with them right now.
C           So just write out, and say which point.
Cdebug           write(6,850)ndata,(jnpone(i),i=1,j)
  850      format(" Pt No.", I5, " : ", 300a1)

C          Now save the data in the local program fields
  855    continue


C       End of data line processing, write it out later, when all are in.
        go to 800
       endif

C     Input line has a # in first col, check it just in case
  880 continue

      if(inpone(1).ne."#")then
C       something is bad in program
       write(0,*)" ERROR 880, sorry prog. messed up call ? admin"
        stop
      endif
C     Ok, its a nice comment.  Figure out if its a special tag.
      read(inp300,*)firstfield
      if(firstfield .eq."#FileType=" .or.
     &   firstfield .eq."#filetype=" .or.
     &   firstfield .eq."#FILETYPE=" )then
         read(inp300,*) firstfield, cfiletype
         if(cfiletype.ne."strain_life")then
C          We have the wrong kind of sae file here folks
           write(0,*)" WARNING: wrong type of SAE std file. ",
     &               " Not #FileType= strain_life"
         endif
C        filetype is ok.  Write it out as a comment
         go to 950
      endif

      if(firstfield .eq."#DataType=" .or.
     &   firstfield .eq."#DATATYPE=" .or.
     &   firstfield .eq."#datatype=" )then
         read(inp300,*) firstfield, Cdatatype
         if(Cdatatype.ne."fitted")then
C          We have the wrong kind of sae file here folks
           write(0,*)" ERROR, wrong type of SAE std file. ",
     &               " Not #DataType= fitted"
           stop
         endif
C        DataType is ok.  Write it out as a comment
         go to 950
      endif

      if(firstfield .eq."#NAME=" .or.
     &   firstfield .eq."#Name=" .or.
     &   firstfield .eq."#name=" )then
         numnames=numnames+1
         read(inp300,*) firstfield, names30(numnames)
         go to 950
      endif

      if(firstfield .eq."#UNITS=" .or.
     &   firstfield .eq."#units=" .or.
     &   firstfield .eq."#Units=" .or.
     &   firstfield .eq."#Stress_units=" .or.
     &   firstfield .eq."#STRESS_UNITS=" .or.
     &   firstfield .eq."#Stress_Units=" .or.
     &   firstfield .eq."#stress_units=" )then
         read(inp300,*) firstfield, stressunits
         go to 950
      endif

      if(firstfield .eq."#STRAIN_UNITS=" .or.
     &   firstfield .eq."#Strain_units=" .or.
     &   firstfield .eq."#STRAIN_units=" .or.
     &   firstfield .eq."#Strain_Units=" .or.
     &   firstfield .eq."#strain_units=" )then
         read(inp300,*) firstfield, strainunits
         go to 950
      endif

      if(firstfield .eq."#LIFEUNITS=" .or.
     &   firstfield .eq."#Life_units=" .or.
     &   firstfield .eq."#LIFE_UNITS=" .or.
     &   firstfield .eq."#Life_Units=" .or.
     &   firstfield .eq."#life_units=" )then
         read(inp300,*) firstfield, lifeunits
         go to 950
      endif
 
      if(firstfield .eq."#Su=" .or.
     &   firstfield .eq."#SU=" )then
         read(inp300,*) firstfield, Sult
         go to 950
      endif
       if(firstfield .eq."#Sy=" .or.
     &    firstfield .eq."#SY=" )then
          read(inp300,*) firstfield, Syield
          go to 950
      endif
       if(firstfield .eq."#E=" .or. firstfield .eq."#e=" .or.
     &    firstfield.eq."#EMOD=" .or. firstfield .eq."#Emod=".or.
     &    firstfield .eq."#emod=" .or.
     &    firstfield .eq."#MODULUS=" .or.
     &    firstfield .eq."#Modulus=" .or.
     &    firstfield .eq."#modulus="  )then
          read(inp300,*) firstfield, EMOD
          go to 950
      endif

       if(firstfield .eq."#%RA=" .or.
     &    firstfield .eq."#%Ra=" )then 
          read(inp300,*) firstfield, percentRA
          go to 950 
      endif 
  
      if(firstfield .eq."#WebPage=" .or.
     &   firstfield .eq."#Webpage=" .or.
     &   firstfield .eq."#WEBPAGE=" .or.
     &   firstfield .eq."#webpage=" )then
         read(inp300,*) firstfield, Cwebpage
         go to 950
      endif

      if(firstfield .eq."#FractureStress=" .or.
     &   firstfield .eq."#fracturestress=" .or.
     &   firstfield .eq."#FRACTURESTRESS=" .or.
     &   firstfield .eq."#Fracturestress=" )then
         read(inp300,*) firstfield, FractureStress
         go to 950
      endif

      if(firstfield .eq."#FractureStrain=" .or.
     &   firstfield .eq."#Fracturestrain=" .or.
     &   firstfield .eq."#fracturestrain=" .or.
     &   firstfield .eq."#FRACTURESTRAIN=" )then
         read(inp300,*) firstfield, FractureStrain
         go to 950
      endif


      if(firstfield .eq."#BHN=" .or.
     &   firstfield .eq."#bhn=" .or.
     &   firstfield .eq."#Bhn=" .or.
     &   firstfield .eq."#HBN=" .or.
     &   firstfield .eq."#HB=" )then
         read(inp300,*) firstfield, xnbrinell
         nbrinell=IFIX(xnbrinell)
         go to 950
      endif

C       Look for things we should skip over in the output
        write(INP80temp,939)
  939   format("#Here is the bottom part of the ",
     &         "html  graph/calc wrapper-----------")
        if(INP300.eq. INP80temp)go to 800

        if(INP300.eq.
     &"#</textarea>"
     &  )go to 800

        if(INP300.eq.
     &"#</pre></DL></FORM></body></html>"
     &  )go to 800





  950 continue
C     Write out the comment line
      j=6  !unit no.
      call WR(j,inp300)
C     Go read another line
      go to 800

  980 continue
C     All input lines have been read.  The comments were
C     put out along the way.  Its time to dump out the data
C     in whatever format the local machine wants it.  This bit
C     is probably site specific.

C     See if some of the critical values are missing:
      istop=0
      if(EMOD.eq.0.)then
        write(0,*)" ERROR #EMOD= missing value"
        write(6,*)" ERROR #EMOD= missing value"
        istop=1
      endif
      if(Cdatatype.eq." ")then
        write(0,*)" ERROR #DataType= missing value"
        write(6,*)" ERROR #DataType= missing value"
        istop=1
      endif
      if(Sult.eq.0.)then
        write(0,*)" WARNING: #Sult= missing value"
        write(6,*)" WARNING: #Sult= missing value"
        write(6,*)" Cannot compute Goodman Damage."
        istop=1
      endif
      if(stressunits .eq. " ")then
        write(0,*)" ERROR #stressunits= missing value"
        write(6,*)" ERROR #stressunits= missing value"
        istop=1
      endif
      if(strainunits .eq. " ")then
        write(0,*)" ERROR #strainunits= missing value"
        write(6,*)" ERROR #strainunits= missing value"
        istop=1
      endif
      if(istop .eq. 1) stop ! ERROR MISSING VALUES.  STOPPING !

      write(6,*)"#CHECK THESE :"
      write(6,*)"#Stress_units=",stressunits," Strain_units=",
     &          strainunits," Life_units=",lifeunits
      write(6,*)"#EMOD=",EMOD, " Sult=",Sult," Syield=",Syield,
     &          " %RA=",percentRA," Fracture_Stress=",FractureStress,
     &          " Fracture_Strain=",FractureStrain
      write(6,*)"#BHN=",nbrinell

      
      if(stressunits .eq."mpa" .or.
     &   stressunits .eq."MPA" .or.
     &   stressunits .eq."MPa" .or.
     &   stressunits .eq."Mpa" )then
         stressunits="MPa"
         go to 1010
      endif
      
      if(stressunits .eq."ksi" .or.
     &   stressunits .eq."KSI" .or.
     &   stressunits .eq."Ksi" )then
C        We really do not need to change units?? What if loads
C        are in mpa and material in ksi?? hm, Ok make everything mpa
         EMOD=EMOD*XMPAS
         Sult=Sult*XMPAS
         Syield=Syield*XMPAS
         FractureStress=FractureStress*XMPAS
         do 983 i=1,ndata
           StressAmp(i)=StressAmp(i)*XMPAS
  983    continue
         write(6,*)"#Material file ksi -> MPa."
         stressunits="MPa"
         go to 1010
      endif
      
      if(stressunits .eq."psi" .or.
     &   stressunits .eq."PSI" .or.
     &   stressunits .eq."Psi" )then
         xtemp =XMPAS/1000.
         EMOD=EMOD*xtemp
         Sult=Sult*xtemp
         Syield=Syield*xtemp
         FractureStress=FractureStress*xtemp
         do 985 i=1,ndata
           StressAmp(i)=StressAmp(i)*xtemp
  985    continue
         write(6,*)"#Material file: psi -> MPa."
         stressunits="MPa"
      endif
 1010 continue

      if(strainunits .eq."Microstrain" .or.
     &   strainunits .eq."MICROSTRAIN" .or.
     &   strainunits .eq."MicroStrain" )then
         FractureStrain=FractureStrain/1.0E6
         do 1020 i=1,ndata
           StrainAmp(i)=StrainAmp(i)/1.0E+6
 1020   continue
         write(6,*)"#Material file Microstrain -> strain"
         strainunits="strain"
      endif

      if(lifeunits .eq."reversals" .or.
     &   lifeunits .eq."REVERSALS" .or.
     &   lifeunits .eq."Reversals" )then
         do 1030 i=1,ndata
           Lifecycles(i)=Lifecycles(i)/2.
 1030   continue
         write(6,*)"#Material file: Reversals -> Cycles"
         lifeunits="cycles"
      endif



C     Now sort the stress,strain, life values  ------------------
      do 1090 i=1,ndata
        ifind=0
        biglife=0.

        do 1080 j=i,ndata
         if(biglife.gt.Lifecycles(j)) go to 1080
C        Nope, its bigger
         biglife=Lifecycles(j)
         ifind=j
 1080  continue

Cdebug       write(6,*)" Sorting ",ifind, " into ",i
C      biglife is now the biggest & is at location j
C      Save it in a temp location
       tempstrain=StrainAmp(ifind)
       tempstress=StressAmp(ifind)
       templife=Lifecycles(ifind)

C      Now move the stuff pointed to by the big loop
C      into "biglife's" location
       StrainAmp(ifind)  =StrainAmp(i)
       StressAmp(ifind)  =StressAmp(i)
       Lifecycles(ifind)=Lifecycles(i)

C      Put the temp stuff into the "i" location
       StrainAmp(i)=tempstrain
       StressAmp(i)=tempstress
       Lifecycles(i)=templife
 1090 continue
      do 1092 i=1,ndata
C        Look for FractureStress & Strain  value at N=0.5
         if(Lifecycles(i) .eq. 0.5)then
C          We have a good approx of  FractureStress & Strain
           FractureStress=StressAmp(i)
           write(6,*)"#Took SigmaF_primed value =",
     &                FractureStress," at N=0.5"
           FractureStrain=StrainAmp(i)
           write(6,*)"#Took Frac.Strain value =",
     &                FractureStrain," at N=0.5"
         endif
Cdebug      write(6,*)StrainAmp(i),StressAmp(i),Lifecycles(i)
 1092 continue
 
      if(FractureStress.eq. 0.)then
C       If no N=0.5 occured we need to extrapolate for FractureStress
        xslope=
     &   ( alog10(StressAmp(ndata-1 )) -alog10(StressAmp(ndata )) ) /
     &   ( alog10(Lifecycles(ndata-1)) -alog10(Lifecycles(ndata)) )
        xlogSigfp= alog10( StressAmp(ndata-1)) -
     &    xslope*( alog10(Lifecycles(ndata-1)) -alog10(0.5) )
        FractureStress=10.0**xlogSigfp
        write(6,*)"# Extrapolated to FractureStress =",
     &             FractureStress
        xslope=
     &   ( alog10(StrainAmp(ndata-1)) -alog10(StrainAmp(ndata)) ) /
     &   ( alog10(Lifecycles(ndata-1)) -alog10(Lifecycles(ndata)) )
        xlogS= alog10(StrainAmp(ndata-1)) -
     &    xslope*( alog10(Lifecycles(ndata-1)) -alog10(0.5) )
        FractureStrain=10.0**xlogS
        write(6,*)"# Extrapolated to FractureStrain =",
     &             FractureStrain
      endif


C     All is well with input data. Fill in the other "columns"
C     for the material file matrix (some folks might 
C     call this a spreadsheet).
Cbugfix  Oct24/04: change to #xcalc1 tag :
C
C#xcalc1 Strain_Amp     Cycles  Stress_Amp  Elas_Str_Amp  Plas_Str_Amp  Smax*Str_Amp  Snominal_Amp
C#xcalc1  0.88485         0.5  795.0  0.00390  0.88095  0.7034255E+03 11967.0
C#xcalc1  0.00914      2500.0  359.2  0.00176  0.00738  0.3283243E+01   817.6
C#xcalc1  0.00665      5000.0  336.5  0.00165  0.00500  0.2237487E+01   674.9
      write(6,1096)
 1096 format("#xcalc1 Strain_Amp     Cycles  Stress_Amp",
     &  "  Elas_Str_Amp ",
     &  " Plas_Str_Amp  Smax*Str_Amp  Snominal_Amp")
      do 1100 i=1,ndata
C       Put out in reverse order to storage
        j=ndata+1-i
        ElstrainAmp(j)=StressAmp(j)/EMOD
        PlstrainAmp(j)=StrainAmp(j)-ElstrainAmp(j)
        SigmaxStrainAmp(j)=StressAmp(j)*StrainAmp(j)
        if(istrainControl .eq. 1)then
C         Its not going to be Neuber,  but just strain history inputs.
          SelasticAmp(j)=StrainAmp(j)
        else
C         It will be a Neuber controlled input history. e.g. FEAelastic stress
          SelasticAmp(j)=SQRT( StressAmp(j)*StrainAmp(j)*EMOD)
        endif
C        save a log copy for making XCALC4 stuff
         logStrainAmp(j)= ALOG10(StrainAmp(j) )
         logLifecycles(j)=ALOG10(Lifecycles(j) )
         logStressAmp(j)= ALOG10(StressAmp(j) )
         logSwtAmp(j)=    ALOG10(SigmaxStrainAmp(j) )
         logSelasticAmp(j)=ALOG10(SelasticAmp(j) )

      write(6,1098) StrainAmp(j),Lifecycles(j),
     &   StressAmp(j),ElstrainAmp(j),PlstrainAmp(j),
     &   SigmaxStrainAmp(j),SelasticAmp(j)
Cbugfix  Oct24/04: change to #xcalc1 tag :
 1098 format("#xcalc1 ",f8.5,1x,f11.1,1x,f6.1,
     &       1x,f8.5,1x,f8.5,1x,E14.7,1x,f7.1)
 1100 continue

      return
      end
   
C   PROGRAMMING PROBLEM.  It would actually be "faster"  for long 
C   rev by rev histories to make seperate  equal interval files:
C         1. for  StsFEA  --> Stress-Strain
C         2. for  Strain --> Stress, life
C         3. for  SWT --> life
C         4. for  Stress --> life

C  This would be ok,  if we have lots of cycles.  If we have only a few
C  it would take more time to prep these files than to do the ALOG10 
C  interpolations of only the original fitted points  as done in the
C  saefcalc2.f program.

C  Thus for the CoilTubing exercise, we are going to abandon this 
C  program for now.  If each segment of the tube needs to be calculated
C  we will need to re-consider this decision.







C==============================================================
      SUBROUTINE getSwat2Life(Smea,Cycles,iexit)
      SAVE
C     Given SmithWatsonTopper product, interpolate Cycles
C     Smea= Smax * totStrainAmpl    (input)
C     Cycles and iexit              (output)

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain
C     Use these to transform to ALOG  to speed up interpolations
      real logStrainAmp(250),logLifecycles(250),logStressAmp(250),
     &     logSwtAmp(250),logSelasticAmp(250)
      Common/MaterialA/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 logStrainAmp,logLifecycles,logStressAmp,
     &                 logSwtAmp,logSelasticAmp,
     &                 ndata,FractureStress,FractureStrain,Emod


      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the Smea is below the first SigmaxStrainAmp point, then it is below
C       the fatigue limit (the first life value in the table).

C      check if data is below the fatigue limit
      if(Smea .lt. SigmaxStrainAmp(1))then
C       Yes, its below. Set it equal to - fat_limit
        Cycles= -Lifecycles(1)  !  Note the negative sign !
        return
      endif

Cbug fixed here in May 19 2004:
C     Check if data is above FractureStress
      if(Smea .ge. SigmaxStrainAmp(ndata))then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Cycles=0.5
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the Smea value
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Smea - SigmaxStrainAmp(idat) ) 100,200,300 !computed goto
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.


  100   continue
C       The Smea value is less than the point on the curve, we have arrived
C       Interpolate.  Since this is Smea vs Life, use log-log interpl.
C        xslope=
C     &     (ALOG10(SigmaxStrainAmp(idat))
C     &         -ALOG10(SigmaxStrainAmp(idat-1))  )
C     &    /(ALOG10(Lifecycles(idat))
C     &         -ALOG10(Lifecycles(idat-1)) )
C        xClog=  (ALOG10(Smea)-ALOG10(SigmaxStrainAmp(idat-1)) )
C     &            /xslope   +ALOG10(Lifecycles(idat-1))
C        Cycles= 10.0**xClog

C       Using pre-calculated log values for speed:
        xslope= (logSwtAmp(idat) -logSwtAmp(idat-1) ) /
     &          (logLifecycles(idat) -logLifecycles(idat-1) )
        xClog= (ALOG10(Smea)-logSwtAmp(idat-1)) / xslope 
     &          +logLifecycles(idat-1)
        Cycles= 10.0**xClog

        return

  200   continue
C       The stress is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the Smea-life curve, and
C       thus the next point would have the same Smea value. Check
Cbugfix   May2005 :
        im=idat
  201   im=im+1
        if(Smea .eq. SigmaxStrainAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger Smea, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Cycles=Lifecycles(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        go to 1000

 1000   continue
        write(0,*)" ERROR:pdStressStrain:getSwat2Life.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Cycles=0.5
        return
        end





C---------------------------------------------------------
      SUBROUTINE getStress2Life(Stress,Cycles,iexit)
      SAVE
C     Given Stress Ampl., interpolate Cycles

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain
C     Use these to transform to ALOG  to speed up interpolations
      real logStrainAmp(250),logLifecycles(250),logStressAmp(250),
     &     logSwtAmp(250),logSelasticAmp(250)
      Common/MaterialA/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 logStrainAmp,logLifecycles,logStressAmp,
     &                 logSwtAmp,logSelasticAmp,
     &                 ndata,FractureStress,FractureStrain,Emod


      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the stress is below the first stress point, then it is below
C       the fatigue limit (the first life value in the table).

C      check if data is below the fatigue limit
      if(Stress .lt. StressAmp(1))then
C       Yes, its below. Set it equal to - fat_limit
        Cycles= -Lifecycles(1)
        return
      endif

C     Check if data is above FractureStress
      if(Stress .ge.FractureStress)then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Cycles=0.5
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Stress - StressAmp(idat) ) 100,200,300
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.


  100   continue
C       Our stress is less than the point on the curve, we have arrived
C       Interpolate.  Since this is Stress vs Life, use log-log interpl.
C        xslope=
C     &     (ALOG10(StressAmp(idat)) -ALOG10(StressAmp(idat-1))  )
C     &    /(ALOG10(Lifecycles(idat))-ALOG10(Lifecycles(idat-1)) )
C        xClog=  (ALOG10(Stress)-ALOG10(StressAmp(idat-1)) ) /xslope
C     &           +ALOG10(Lifecycles(idat-1))
C        Cycles= 10.0**xClog

C       Use pre-calculated log values for speed.
        xslope= (logStressAmp(idat)-logStressAmp(idat-1) ) /
     &          (logLifecycles(idat)-logLifeCycles(idat-1) )
        xClog= (ALOG10(Stress)-logStressAmp(idat-1) ) /xslope
     &          + logLifecycles(idat-1)
        Cycles= 10.0**xClog
        return

  200   continue
C       The stress is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the stress-life curve, and
C       thus the next point would have the same stress value. Check
        im=idat
  201   im=im+1
        if(Stress .eq. StressAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger stress, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Cycles=Lifecycles(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        go to 1000

 1000   continue
        write(0,*)" ERROR:pdStressStrain:getStress2Life.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Cycles=0.5
        return
        end




C---------------------------------------------------------
      SUBROUTINE getStrain2Life(Strain,Cycles,iexit)
      SAVE
C     Given Strain Ampl., interpolate Cycles

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain
C     Use these to transform to ALOG  to speed up interpolations
      real logStrainAmp(250),logLifecycles(250),logStressAmp(250),
     &     logSwtAmp(250),logSelasticAmp(250)
      Common/MaterialA/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 logStrainAmp,logLifecycles,logStressAmp,
     &                 logSwtAmp,logSelasticAmp,
     &                 ndata,FractureStress,FractureStrain,Emod

      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the strain is below the first strain point, then it is below
C       the fatigue limit (the first life value in the table).

C      check if data is below the fatigue limit
      if(Strain .lt. StrainAmp(1))then
C       Yes, its below. Set it equal to - fat_limit
        Cycles= -Lifecycles(1)
        return
      endif

C     Check if data is above FractureStrain
      if(Strain .ge.FractureStrain)then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Cycles=0.5
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Strain - StrainAmp(idat) ) 100,200,300
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.


  100   continue
C       Our stress is less than the point on the curve, we have arrived
C       Interpolate.  Since this is Strain vs Life, use log-log interpl.
C        xslope=
C     &     (ALOG10(StrainAmp(idat)) -ALOG10(StrainAmp(idat-1))  )
C     &    /(ALOG10(Lifecycles(idat))-ALOG10(Lifecycles(idat-1)) )
C        xClog=  (ALOG10(Strain)-ALOG10(StrainAmp(idat-1)) )/xslope
C     &          +ALOG10(Lifecycles(idat-1))
C        Cycles= 10.0**xClog

C       Use pre-computed log values for speed
        xslope= ( logStrainAmp(idat)-logStrainAmp(idat-1) ) /
     &          ( logLifecycles(idat)-logLifecycles(idat-1) )
        xClog= ( ALOG10(Strain) - logStrainAmp(idat-1) ) /xslope
     &          + logLifecycles(idat-1)
        Cycles= 10.0**xClog

        return

  200   continue
C       The Strain is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the Strain-life curve, and
C       thus the next point would have the same Strain value. Check
        im=idat
  201   im=im+1
        if(Strain .eq. StrainAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger Strain, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Cycles=Lifecycles(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going

 1000   continue
         write(0,*)" ERROR:pdStressStrain:getStrain2Life.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Cycles=0.5
        return
        end

           
C========================================================================
      SUBROUTINE getStress2Strain(Stress,Strain,iexit)
      SAVE
C     Given Stress Ampl., interpolate Strain

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain
C     Use these to transform to ALOG  to speed up interpolations
      real logStrainAmp(250),logLifecycles(250),logStressAmp(250),
     &     logSwtAmp(250),logSelasticAmp(250)
      Common/MaterialA/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 logStrainAmp,logLifecycles,logStressAmp,
     &                 logSwtAmp,logSelasticAmp,
     &                 ndata,FractureStress,FractureStrain,Emod


      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the stress is below the first stress point, then it is below
C       the fatigue limit (the first life value in the table).

      if(Stress .lt. StressAmp(1))then
C       Yes, its below. Assume straight line from (0,0)
        Strain=(Stress/StressAmp(1))*StrainAmp(1) 
C       This can cause some problems when not equal to E modulus,
C       but should be ok. The data is what it is.
        return
      endif

C     Check if data is above Sigfp
      if(Stress .ge. Fracturestress)then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Strain=Fracturestrain
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.
      do 1000 idat=1,ndata
        if(Stress - StressAmp(idat) ) 100,200,300  !computed goto
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.


  100   continue
C       Our stress is less than the point on the curve, we have arrived
C       Interpolate.  Try log-log interpl.  If we ever get -ve stresses
C       or strains we will definitely have to go to linear interpl.
C        xslope=
C     &     (ALOG10(StressAmp(idat)) -ALOG10(StressAmp(idat-1)))
C     &    /(ALOG10(StrainAmp(idat))-ALOG10(StrainAmp(idat-1)) )
C        xClog=  (ALOG10(Stress)-ALOG10(StressAmp(idat-1)) )/xslope
C     &          +ALOG10(StrainAmp(idat-1))
C        Strain= 10.0**xClog

C       Use pre-computed log values for speed:
        xslope=( logStressAmp(idat) - logStressAmp(idat-1) ) /
     &         ( logStrainAmp(idat) - logStrainAmp(idat-1) )
        xClog= (ALOG10(Stress)-logStressAmp(idat-1) ) /xslope
     &          +logStrainAmp(idat-1)
        Strain= 10.0**xClog
        return

  200   continue
C       The stress is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the stress-strain curve, and
C       thus the next point would have the same stress value. Check
Cbugfix   May2005 : im didnt increment
        im=idat
  201   im=im+1
        if(Stress .eq. StressAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger stress, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Strain=StrainAmp(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        go to 1000
 1000   continue

        write(0,*)" ERROR:pdStressStrain:getStress2Strain.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Strain=Fracturestrain
        return
        end




C==========================================================================
      SUBROUTINE getLoad2StressStrain(Sneuber,Stress,Strain,
     &           iexit)
      SAVE
C     Given Load_Amp (could be FEA Stress Ampl.), interpolate Strain &Stress
C     The saefcalc2.f and saefcalc3.f versions use log-log interp.
C     The Crack Prop versions use linear interp. 
C     of ~1000 equally space points for speed.

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain
C     Use these to transform to ALOG  to speed up interpolations
      real logStrainAmp(250),logLifecycles(250),logStressAmp(250),
     &     logSwtAmp(250),logSelasticAmp(250)
      Common/MaterialA/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 logStrainAmp,logLifecycles,logStressAmp,
     &                 logSwtAmp,logSelasticAmp,
     &                 ndata,FractureStress,FractureStrain,Emod


      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the "load", actually = Elastic_local_Stress, is below 
C       the first SelasticAmp point, then it is below
C       the fatigue limit (the first life value in the table).

      if(Sneuber .lt. SelasticAmp(1))then
C       Yes, its below. Assume straight line from (0,0)
        Strain=(Sneuber/SelasticAmp(1))*StrainAmp(1)
        Stress=(Sneuber/SelasticAmp(1))*StressAmp(1)
        return
      endif

C     Check if data is above Sigfp
C     Assumes: ElasMod= First Stress/ 1st strain
C     If these first points are not on the "elastic" line we could have
C     a problem here, but it shouldnt be too bad (?)
      E=StressAmp(1)/StrainAmp(1)
      if(Sneuber .ge.
     &        SQRT(Fracturestress*Fracturestrain*E))then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Strain=Fracturestrain
        Stress=Fracturestress
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Sneuber - SelasticAmp(idat) ) 100,200,300
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.

  100   continue
C       Our stress is less than the point on the curve, we have arrived
C       Interpolate.  Try log-log interpl.  If we ever get -ve stresses
C       or strains we will definitely have to go to linear interpl.
C        xslope=
C     &     (ALOG10(SelasticAmp(idat)) -ALOG10(SelasticAmp(idat-1)))
C     &    /(ALOG10(StressAmp(idat))-ALOG10(StressAmp(idat-1)) )
C        xClog=  (ALOG10(Sneuber)-ALOG10(SelasticAmp(idat-1)) )/xslope
C     &          +ALOG10(StressAmp(idat-1))
C        Stress= 10.0**xClog
C        xslope=
C     &     (ALOG10(SelasticAmp(idat)) -ALOG10(SelasticAmp(idat-1)))
C     &    /(ALOG10(StrainAmp(idat))-ALOG10(StrainAmp(idat-1)) )
C        xClog=  (ALOG10(Sneuber)-ALOG10(SelasticAmp(idat-1)) )/xslope
C     &          +ALOG10(StrainAmp(idat-1))
C        Strain= 10.0**xClog

        xslope=
     &     (logSelasticAmp(idat) -logSelasticAmp(idat-1) )
     &    /(logStressAmp(idat)   -logStressAmp(idat-1)   ) 
        xClog=  (ALOG10(Sneuber)-logSelasticAmp(idat-1) )/xslope
     &          + logStressAmp(idat-1)
        Stress= 10.0**xClog
        xslope=
     &     (logSelasticAmp(idat) -logSelasticAmp(idat-1) )
     &    /(logStrainAmp(idat)   -logStrainAmp(idat-1)   )
        xClog=  (ALOG10(Sneuber)-logSelasticAmp(idat-1)  ) /xslope
     &          +logStrainAmp(idat-1)
        Strain= 10.0**xClog

        return

Cx  100   continue
CxC       The stress is less than the point on the curve, we have arrived.
CxC       Interpolate.  Try log-log interpl.  If we ever get -ve stresses
CxC       or strains we will definitely have to go to linear interpl.
Cx
Cx        fraction=(ALOG10(Sneuber)-ALOG10(SelasticAmp(idat-1))  )/
Cx     &    (ALOG10(SelasticAmp(idat))-ALOG10(SelasticAmp(idat-1)))
Cx        Clog=ALOG10(StressAmp(idat-1)) + fraction *
Cx     &       (ALOG10(StressAmp(idat))-ALOG10(StressAmp(idat-1)))
Cx        Stress= 10.0**Clog

CxC       Since the Sneuber,Strain, Stress occur as triple points, we
CxC       can use the same fraction for Strain interpolation:
Cx        Clog=ALOG10(StrainAmp(idat-1)) + fraction *
Cx     &       (ALOG10(StrainAmp(idat))-ALOG10(StrainAmp(idat-1)))
Cx        Strain= 10.0**Clog
Cx        return

  200   continue
C       Sneuber is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the stress-strain curve, and
C       thus the next point would have the same stress value. 
C       Check :
Cbugfix   May2005 : im didnt increment.
        im=idat
  201   im=im+1
        if(Sneuber .eq. SelasticAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger number, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Strain=StrainAmp(im-1)
        Stress=StressAmp(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        go to 1000

 1000   continue

        write(0,*)" ERROR:pdStressStrain:getload2stressStrain:",
     &            " endofDOloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Strain=Fracturestrain
        Stress=Fracturestress
        return
        end


C=======================================================================
      SUBROUTINE getloop(Stress,Straintemp,Stresstemp,npts,iexit)
      SAVE
C     Given Stress Ampl., get all stress-strain points from 0 to Stress
C     Routine used in saefcalc2.f   Not used in  pdStressStrain.f 

C     These are the returned value storage:
      real Stresstemp(250),Straintemp(250)


      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain
C     Use these to transform to ALOG  to speed up interpolations
      real logStrainAmp(250),logLifecycles(250),logStressAmp(250),
     &     logSwtAmp(250),logSelasticAmp(250)
      Common/MaterialA/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 logStrainAmp,logLifecycles,logStressAmp,
     &                 logSwtAmp,logSelasticAmp,
     &                 ndata,FractureStress,FractureStrain,Emod

      iexit=0

      Stresstemp(1)=0.
      Straintemp(1)=0
      npts=1
C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the stress is below the first stress point, then it is below
C       the fatigue limit (the first life value in the table).

      if(Stress .lt. StressAmp(1))then
C       Yes, its below. Assume straight line from (0,0)
        npts=npts+1
        Straintemp(npts)=(Stress/StressAmp(1))*StrainAmp(1)
        Stresstemp(npts)=Stress
        go to 9000
      endif

C     Check if data is above Sigfp
      if(Stress .ge. Fracturestress)then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        do 50 i=1,ndata
          npts=npts+1
          Straintemp(npts)=StrainAmp(i)
          Stresstemp(npts)=StressAmp(i)
   50   continue
        npts=npts+1
        Straintemp(npts)=Fracturestrain
        Stresstemp(npts)=Fracturestress
        iexit=1
        go to 9000
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Stress - StressAmp(idat) ) 100,200,300
C         In-elegant but fast:        -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.

  100   continue
C       The stress is less than the point on the curve, we have arrived
C       Interpolate.  Try log-log interpl.  If we ever get -ve stresses
C       or strains we will definitely have to go to linear interpl. but
C       then we will need many more points on the strain life curve (?)
C       (perhaps its same result ?)
C        xslope=
C     &     (ALOG10(StressAmp(idat)) -ALOG10(StressAmp(idat-1)))
C     &    /(ALOG10(StrainAmp(idat))-ALOG10(StrainAmp(idat-1)) )
C        xClog=  (ALOG10(Stress)-ALOG10(StressAmp(idat-1)) )/xslope
C     &          +ALOG10(StrainAmp(idat-1))

C       Use pre-calculated log values for speed:
        xslope=( logStressAmp(idat) -logStressAmp(idat-1)  )
     &        /( logStrainAmp(idat) -logStrainAmp(idat-1)  )
        xClog= (ALOG10(Stress) -logStressAmp(idat-1) ) /xslope
     &          +logStrainAmp(idat-1)

        npts=npts+1
        Straintemp(npts)= 10.0**xClog

C        deltS=StressAmp(idat)-StressAmp(idat-1)
C        fraction=(Stress-StressAmp(idat-1))/deltS
C        Stresstemp(npts)=StressAmp(idat-1) + fraction*deltS
        Stresstemp(npts)=Stress
        go to 9000

  200   continue
C       The stress is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the stress-strain curve, and
C       thus the next point would have the same stress value. Check
Cbugfix   May2005 : im didnt increment
        im=idat
  201   im=im+1
        if(Stress .eq. StressAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger stress, use the last
C       one that was equal as a match
        npts=npts+1
        Straintemp(npts)=StrainAmp(im-1)
        Stresstemp(npts)=Stress
        go to 9000
  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        npts=npts+1
        Straintemp(npts)=StrainAmp(idat)
        Stresstemp(npts)=StressAmp(idat)

C      end of hunting loop
 1000  continue

        write(0,*)" ERROR:pdStressStrain:getloop.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        npts=npts+1
        Straintemp(npts)=Fracturestrain
        Stresstemp(npts)=Fracturestress
        go to 9000

 9000   continue
C        write(6,*)"#debugPlotLoop for Stress Amp: ",Stress
C        write(6,*)"#debugPlotLoop   Amplitudes        Ampl.*2"
C        do 9005 i=1,npts
C        write(6,9003)Straintemp(i),Stresstemp(i),
C     &   Straintemp(i)*2.0,Stresstemp(i)*2.0
C 9003   format("#debugPlotLoop ",1x,f7.5,1x,f6.1,4x,f7.5,1x,f6.1)
 9005   continue
        end


