C deltaEvents.f  v1.1  remove lines from a file before and after two events
      SAVE
C  Compile:  gfortran  -g -w -fbounds-check deltaEvents.f  -o deltaEvents
C  Usage eg.:   deltaEvents  601.000  1600.320   <infile >outfile
C         Deletes all data before event 601.000  and after ev. 1600.320
C         Also computes the accumulated stress-strain energy in sequence.

C  Max no. of chars per line is 300
C  Comment lines between the two events are left in.

C  Copyright (C) 2012  Al Conle
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

C vers 1.1 Add energy computation as stress-strain locus moves

C Items below are :
C    Strain  Stress     Event_time
C Input file e.g.:
C   0.00000    0.0         0.000
C   0.00000    0.0       598.000
C   0.0000000       0.0000000       601.00000000000000
C   0.00131  254.6       601.000
C   0.00131  254.6       601.000
C   0.00145  271.9       601.000
C   0.00158  284.3       601.000
C   0.00175  296.1       601.000
C   0.00204  310.8       601.000
C   0.00234  321.5       601.000
C   0.00273  332.1       601.000
C   0.00346  346.0       601.000
C   0.00424  356.6       601.000
C   0.00532  367.4       601.000
C   0.00739  382.1       601.000
C   0.00800  385.4       601.000
C   #StrOut, StsOut=   8.00000038E-03   385.44385
C   0.00800  385.4       999.000
C   0.00800  385.4      1000.100
C   0.00800  385.4      1298.000
C   0.00800  385.4      1300.020
C   0.00800  385.4      1598.000
C   0.00795  375.7      1598.116
C   0.00800  385.4      1600.320
C   0.00900  390.5      1600.320
C   #StrOut, StsOut=   8.99999868E-03   390.46323
C   0.00638 -118.7      2048.000
C   0.00638 -118.7      2048.000

      character*300  inp300 ! used to read in lines as chars.
      character*1    inpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)

      character*80 argv
      integer*4 iargc,nargc

C     Time must be in  double precision
      real*8  eventTime1,eventTime2,xtime,dtime

      real  energy,energy1,energy2,totalenergy

      write(6,165)
      write(0,165)
  165 format("# deltaEvents.f vers. 1.1"/
     & "#Usage e.g.: deltaEvents 602.000  2048.000  <infile  >outfile"/
     & "#        ")

      nargc = iargc()
C     Note that in HP fortran the arg numbers must be changed by -1
C     and that iargc probably includes the "deltaEvents" as an arg.
      if( nargc .ne. 2)then
        write(0,174)
        write(6,174)
  174   format("# deltaEvents:  usage ERROR"/
     &  "# Usage e.g.: deltaEvents ev1 ev2 <infile  >outfile")
        stop
      endif

C       The first arg is the no. of points removed between saved point
        jvect=1
        call getarg(jvect,argv)
        read(argv,*,err= 176)eventTime1
        write(0,*)"#EventTime1= ",eventTime1
        write(6,*)"#EventTime1= ",eventTime1
        go to 178
C     Bad arguments in command line
  176 continue
        write(0,174)
        write(6,174)
      stop

  178 continue
        jvect=jvect+1
        call getarg(jvect,argv)
        read(argv,*,err= 179)eventTime2
        write(0,*)"#EventTime2= ",eventTime2
        write(6,*)"#EventTime2= ",eventTime2
        go to 180
C     Bad arguments in command line
  179 continue
        write(0,174)
        write(6,174)
      stop




  180 continue


  200 continue
C     Get the first data line with eventTime1
      
      ninput=0
      ndelete=0
      nsavedata=0

  700 continue
c     Loop back to here for next input line.
      read(5,"(a300)",end=750)inp300
      ninput=ninput+1
Cdebug      write(0,*)" read input line ",ninput

C     Check for blank line
      if(inp300.eq." ")then
C        write(6,"(a1)")" "
        go to 700
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 705 i=1,300
           if(inpone(i).eq." ") go to 705
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 703 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  703        continue
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 730
           else
C            The first non blank is not a #
             go to 710
           endif
  705      continue
  710      continue
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
            write(0,715)ninput
            write(6,715)ninput
  715       format("# ERROR : input line no. ",I5,
     &      " not a # and not a number. Fix input file. ")
            call WR(0,inpone)
            stop
         endif

C        Ok, its a number  Fetch in the  Strain  Stress and  eventTime
         read(inp300,*)xstrain,ystress,xtime
         ndata=ndata+1
         if(ndata .eq. 1)then !set up for energy calc
           oldstrain=xstrain
           oldstress=ystress
           energy=0.
           totalenergy=0.
         endif

         if(xtime .lt.eventTime1)goto 700  ! not there yet, keep reading
C        Yes, it is above eventTime1
         if(xtime .gt. eventTime2)goto 750  ! we are done. Stop

C        Ok, it is in our desired even time interval. Save it.
         nsavedata=nsavedata+1
         if(ndata .gt. 1)then
C          energy segment is average stress * delta strain
C          Note that it is possible to a -ve energy here.  e.g. by
C          unloading one does not need to apply energy.  It is a spring.
C          Energy comes out, and is subtracted from the total.
C          Energy is put into the material by either plastic straining
C          in the tensile or in the compressive direction.  Thus is one
C          performs a small loop inside a bigger path  this loop too is
C          input into shear on slip planes = heat.
C            Plan:
C          1. any unloading towards zero gets its energy subtracted from
C             the total.
C          2. Any loading, above or below zero that causes plastic strain
C             or elastic strain is added to the sum.  When loading reverses
C             again the elastic stored energy is subtracted from the total
C            but the previous plastic strain energy always stays in the sum.

C          The Stress-strain points in the files define line segments.
C          We need to check the two above rules for each segment.

           if(ystress .eq. 0. .and. oldstress .eq. 0.)then
             energy=0.
             goto 720
           endif

           if(ystress .ge. 0. .and. oldstress .ge. 0.)then
C            line segment is above zero
             energy=(((ystress+oldstress)/2.0)  * (xstrain-oldstrain) )
C            Note that the sign of (xstrain-oldstrain) will show up or down
C            and energy wil be +ve or -ve accordingly.
             totalenergy=totalenergy + energy
             goto 720
           endif

           if(ystress .le. 0. .and. oldstress .le.0.)then
C            line segment is below zero
             avedstress= abs(ystress+oldstress)/2.0 !average stress
             dstrain= (xstrain - oldstrain)  !-ve means  down slope.
C            If downslope Add energy,  if upslope  subtract
C            Thus multiply by  -1     (See notes if confused).
             energy= -(avedstress * dstrain)
             totalenergy=totalenergy + energy
             goto 720
           endif

C          The neither of above, then two points straddle  the zero load line
C          Compute strain at zero crossing
           if(ystress .gt. 0.)then  !we are going upwards
             dstrain=(xstrain-oldstrain) ! +ve no.
             dstress=(ystress-oldstress) ! +ve no.
             frac=(0.-oldstress)/dstress ! oldstress is -ve, thus frac +ve
             zstrain=oldstrain + frac*dstrain
C            upwards implies that below zero stress, energy is unloading -ve
C            Above zero is loading and  +ve energy. Thus
             energy1=-( (0.0-oldstress)*(zstrain-oldstrain) ) / 2.0
C            Yes, I know,  its just for clarity.
             totalenergy=totalenergy+energy1  ! energy is -ve

             energy2= ( (ystress-0.0)*(xstrain-zstrain) ) /2.0
C            since we are going up this gets added
             totalenergy=totalenergy+energy2
             energy=energy1+energy2
             goto 720
           endif

C          We must be going downwards
           if(ystress .gt. oldstress)then ! error check.
             write(0,721)oldstrain,oldstress,xstrain,ystress
             write(6,721)oldstrain,oldstress,xstrain,ystress
  721        format("#ERROR: deltaevents: logic in energy compute:"/
     &       "#oldstrain,oldstress,xstrain,ystress= ",2(f7.5,1x,f6.1) )
             stop
           endif
C          Ok, going downwards.  Find strain at zero load first.
           dstrain=xstrain-oldstrain !-ve number
           frac=abs(0.0-oldstress)/(abs(ystress-oldstress) )
           zstrain=(oldstrain+(frac*dstrain) ) !where dstrain is -ve
C          Since we are going downwards the above zero load energy is
C          unloading energy. It will be subtracted from total.
           energy1= ( (oldstress -0.)*(zstrain-oldstrain) ) /2.0 ! -ve str rg.
C          Due to the strain fraction energy1 is a -ve no.
           totalenergy=totalenergy+energy1 ! note addition of a -ve no.

           energy2= (abs(ystress-0.))*(abs(xstrain-zstrain)) ! a +ve no.
C          we are loading up in compression, thus +ve energy.
           totalenergy=totalenergy+energy2  ! add the energy
           energy=energy1+energy2
           goto 720

         endif!  end of:  if(ndata .gt. 1)then

  720    continue
         write(6,724)xstrain,ystress,xtime,energy,totalenergy
         oldstress=ystress
         oldstrain=xstrain
         
  724    format(f7.5,1x,f6.1,1x,f13.3,1x,e14.7,1x,e14.7)
C         call WR(6,inp300)
         go to 700 ! go get more input

         
        endif
  730    continue
C        Simply write out the comment line
         if(nsavedata .ne. 0) then  ! we are in time interval, save comm.
           call WR(6,inp300)
         endif
         go to 700

  750    continue   ! end of file comes here
C      All data is in. 
       write(0,755)ninput,ndata,nsavedata
  755  format("#deltaEvents Done. Tot. input lines incl. comments= ",i8/
     &  "#of which ",i8," were data lines."/
     &  "# Event Interval data Lines saved= ",i8)
       stop
       end



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
   24 FORMAT(200A1)
      RETURN
      END

