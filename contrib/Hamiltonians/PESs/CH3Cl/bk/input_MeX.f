C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE LES(X,Y,NUM,WT,PAR,NPAR,MA,SEED,GTYPE,NSET,CHKPNT,RSTRT,
C%      OLD,ITER,MAXIT,MICIT,NSEL,SPREAD,MUT,DIFPER,FREEZE)
C%
C% Reads all parameters and data from genetic input file. The input
C% is controled by so called cards. The fall into two types: One are
C% the system specific cards (set in keys.incl) and the others are 
C% the general cards for genetic listet later in the header. Every
C% card ends with a colon. It does not matter if its written in capital
C% or in small letters. Lines with an exclamation mark are ignort. They
C% are like comments. A line can be wraped by setting an & at the and
C% of the line. There is an example for an input in  example_input.in.
C%
C% Attention: Interface changed (7.4.16), new variable hybrid
C% Attention: Interface changed (28.6.16), new variable so
C%
C% Input cards:
C% CHK:      Name of checkpointfile,<meake new or use old one>
C% COORD:    Coordinate for plot
C% DATA:     Data for Fit
C% DIFPER:   Minimum difference of selected parent parameter sets [%]
C% FREEZE:   No fit will be done
C% GSPREAD:  General spread for each parameter
C% GTYPE:    integer that selects random number generator [0-5]
C% MAXIT:    Number of macro iterations
C% MICIT:    Maxiumum number of micro iterations
C% MUT:      percentage number of how many mutations occur in parameters
C% NSET:     Number of Parameter sets that should be generated
C% NPOINTS:  Number of points per dataset
C% RANDOM:   Random parameters will be used
C% SEED:     Seed for random.f
C% SEL:      percentege of parents that should have babies
C% SETS:     Number of datasets
C% ACTSTATE: Active Cards for States. The number of active Cards must
C%           be equal to the number of states. 0: no fit, 1: fit
C% HYBRID:   CIs are going to be fitted to
C% SO:       If set then the Spin orbit fit will be done
C%
C% Recomandation for naming Parameter Cards:
C% They must be shorter than 12 Characters and this Prefixes should be
C% used:
C% N number of parameters 
C% P parameter values 
C% A active parameters (0 or 1)
C% S Spread of parameters 
C% 
C% Variables of the Function:
C% x:      vector of coordinates for every point (double[])
C% y:      vector of data for every point (double[])
C% num:    number of read datapoints (int)
C% wt:     vector of weights for every point (double[])
C% par_:    vector of parameters (doube[])
C% npar:   number of read parameters (int)
C% ma:     ative cards (int[])
C% seed:   seed for random.f (int)
C% gtype:  selects random number generator [0-5] (int)
C% nset:   number of read datasets (int)
C% chkpnt: name of checkpoint file (char)
C% rstrt:  restart variable (0 = none, 1 = new, 2 = restart) (int)
C% old:    (double)
C% iter:   Sarting point after checkpoint (int)
C% maxit:  number of macro iterations (int)
C% micit:  maximum number of micro iterations (int)
C% nsel:   number of parents selected for having babies (int)
C% spread: spread vector for parameters (double[])
C% mut:    percentage number of how many mutations occur in parameters
C%         (double)
C% difper: minimum difference of selected parent parameter sets [%]
C%         (double)
C% freeze: fit or not (logical)
C% hybrid: 1.d0: CIs are going to be fitted to, 0.d0 diabatisazion
C%         per Ansatz (double)
C% so:     if .true. so fit will be done
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine les(x,y,num,wt,par_,npar,ma,seed,gtype,nset,
     $     chkpnt,rstrt,old,iter,maxit,micit,nsel,spread,mut,difper,
     $     freeze,hybrid,so)
      implicit none

      include 'states.incl'    !some System informations
      include 'params.incl' !some array sizes
      include 'common.incl' !pst, keys and so on
      include 'morse.incl'  !???
      include 'plot.incl'   !Some variables for write (plotfiles)
      include 'damp.incl'   !Variables for dampig function

c     Functions
      integer iargc

c     Variables
      integer mmax !maximum file length in lines
      integer sl   !maximum linelength in characters
      parameter (mmax=2500000,sl=1000)

      integer num !number of datapoints and index for reading them
      integer ii  !index for ndata array
      integer jj  !number of read datapoints

      double precision x(qn,ntot*mx) !read coordinates
      double precision y(ntot*mx)    !read energy values and ci values

      integer wtact(nstat) !active crads for states

      double precision par_(px*nx) !parametervector
      integer npar                !total number of given parameters
      integer nset                !number of parametersets
      integer ma(px)              !active cards for each parameter
      double precision spread(px) !size of manipulation of parametersets
      double precision gspread    !Genereral spread for each parameter

      double precision mut    !percentege number of mutation
      integer nsel            !number of parents having babies
      double precision sel    !percentege of parents having babies
      double precision difper !minimum difference of selected parent

      double precision wt(ntot*mx) !vector of weights for each point
      double precision hybrid      !ci fit or not

      double precision old !???

      integer lauf !number of lines without comments

      integer i,j,k !running indeces
      
      integer rstrt !restart variable (0 = none, 1 = new, 2 = restart)

      integer maxit !number of macro iterations (with mutation of
                    !the former beat parameters)
      integer micit !maximum number of iterations in one marq run
                    !(micro iterations)

      integer seed !Seed for random.f

      integer gtype !selects random number generator

      integer pcount !number of read parameters for each block

      character str*(sl)      !line string (contains one line)
      character in(mmax)*(sl) !input string (contains file)
      integer st              !counts the number of read lines

      character datnam*30   !name of input file
      character chkpnt*(80) !name of checkpoint file
      integer iter          !Sarting point after checkpoint
      logical chk           !exists the checkpoint file
      logical lchk          !was a checkpoint file read
      logical new           !Will a new checkpoint be created
      character ntst*3      !for reading new marker of chkpoint

      character reader*10   !='read' for calling chkfile str argument

      logical lrnd !Should Random parameters be used

      logical dbg !More Output for debugging

      logical freeze !fit or not

      logical so !spin orbit fit or not

      logical signswap !sign change at start

      include 'keys.incl' !System keys for parameters

      dbg=.false.
c      dbg=.true.

c     initialization of par_
      do i=1,px*nx
         par_(i)=0.d0
      enddo

!-------------------------------------------------------------------
!.....read input from disk:
      j=iargc()
      if (j.gt.0) then
        call getarg(1,datnam)
      else
        write (6,100)
        read (5,*) datnam
      endif
100   format ('input file name?')

      open (10,file=datnam)

      lauf=0
      do i=1,mmax
        read (10,'(a750)',end=20) str
        if (str(1:3).eq.'---') goto 20          !end of input mark
         call capital(in,str,lauf,mmax,sl)
      enddo

20    close (10)
!-------------------------------------------------------------------

      st=0

!..   read 'RANDOM:' card:
      lrnd=.false.
      do i=1,lauf
        if (in(i)(1:7).eq.'RANDOM:') then
          if (i+1.gt.st) st=i+1
	  lrnd=.true.
	  write (6,70)
	  goto 69
	endif
      enddo
69    continue
70    format('Random parameters will be chosen')

!-------------------------------------------------------------------
!..   read random seed:
      seed=8236475
      do i=1,lauf
        if (in(i)(1:5).eq.'SEED:') then
	  read (in(i)(6:sl),*) seed
          if (i+1.gt.st) st=i+1
	  goto 75
	endif
      enddo
      write (6,76) seed
75    continue
76    format('No random seed specified; seed set to',i9)

!-------------------------------------------------------------------
!..   read random generator type GTYPE:
      gtype=1
      do i=1,lauf
        if (in(i)(1:6).eq.'GTYPE:') then
	  read (in(i)(7:sl),*) gtype
          if (i+1.gt.st) st=i+1
	  if ((gtype.lt.0).or.(gtype.gt.5)) then
	    write (6,87) gtype
	    stop
	  endif
	  goto 85
	endif
      enddo
      write (6,86) gtype
85    continue
86    format('No random generator type specified; gtype set to',i2)
87    format('Error: GTYPE=',i1,' out of allowed range! [0-5]')

!-------------------------------------------------------------------
!..   read number of parameter sets:
      nset =0
      do i=1,lauf
        if (in(i)(1:5).eq.'NSET:') then
            read (in(i)(6:sl),*) nset
        endif
      enddo
      if (nset.le.0) then
        write(6,'(''Error: no NSET card found'')')
        stop
      endif  

!-------------------------------------------------------------------
!..   we select the sel*100 percent best of the input parameter sets:
      sel=0.15d0
      do i=1,lauf
        if (in(i)(1:4).eq.'SEL:') then
	  read (in(i)(5:sl),*) sel
          if (i+1.gt.st) st=i+1
	  goto 200
	endif
      enddo
      write (6,201) sel*100.
200   continue
201   format('No SEL specified, default of',f5.1,'% sets will be selecte
     $d as parents')
      nsel=nint(nset*sel)
      if(nsel.eq.0) nsel=1

!-------------------------------------------------------------------
c     read DIFPER: card to select minimum percentage of difference between
c     parent parameter sets:
      difper=0.05d0
      do i=1,lauf
        if (in(i)(1:7).eq.'DIFPER:') then
          read (in(i)(8:sl),*) difper
          if (i+1.gt.st) st=i+1
          goto 205
        endif
      enddo
      write(6,206) difper
205   continue
206   format('No DIFPER specified, default difference set to',f5.1)

!-------------------------------------------------------------------
!..   read maximum number of iterations:
      maxit=100
      do i=1,lauf
        if (in(i)(1:6).eq.'MAXIT:') then
	  read (in(i)(7:sl),*) maxit
          if (i+1.gt.st) st=i+1
	  goto 210
	endif
      enddo
      write (6,211) maxit
210   continue
211   format('No MAXIT specified; default number of iterations is',i5)
      old=1.e+5

!-------------------------------------------------------------------
!..   read maximum number of micro iterations:
      micit=1000
      do i=1,lauf
        if (in(i)(1:6).eq.'MICIT:') then
          read (in(i)(7:sl),*) micit
          if (i+1.gt.st) st=i+1
          goto 215
        endif
      enddo
      write (6,216) micit
215   continue
216   format('No MICIT specified; default number of iterations is',i5)

!-------------------------------------------------------------------
!..   read general spread for each parameter:
      gspread=1.d2
      do i=1,lauf
        if (in(i)(1:8).eq.'GSPREAD:') then
          read(in(i)(9:sl),*) gspread
          goto 333
        endif  
      enddo
333   continue
      write(6,'(''GSPREAD set to'',g12.1)') gspread


!-------------------------------------------------------------------
!..   read general mutation parameter:
      mut=0.d0
      do i=1,lauf
        if (in(i)(1:4).eq.'MUT:') then
          read(in(i)(5:sl),*) mut
          mut=mut/1.d2
          goto 343
        endif  
      enddo
343   continue
      write(6,'(''MUTATION set to'',g12.1)') mut

!..   refers to the number of data sets
      do i=1,lauf
        if (in(i)(1:5).eq.'SETS:') then
            read (in(i)(6:sl),*) sets
        endif
      enddo          
      if (dbg) write(6,*) 'no. of data sets:', sets
      if (sets.gt.sx) then
        write(6,'(''Error: Too many data sets! Increase sx'')')
	stop
      endif	
            
!------------------------------------------------------------------------
c     read parameters for dynamic weights and split reference and 
c     correction scheme:
      wtref=1.d20
      wtrho=1.d20
      e0ref=1.d20
      e0rho=1.d20
      do i=1,lauf
        if (in(i)(1:6).eq.'WTREF:') then
            read (in(i)(7:sl),*) wtref
        endif
        if (in(i)(1:6).eq.'WTRHO:') then
            read (in(i)(7:sl),*) wtrho
        endif
        if (in(i)(1:6).eq.'E0REF:') then
            read (in(i)(7:sl),*) e0ref
        endif
        if (in(i)(1:6).eq.'E0RHO:') then
            read (in(i)(7:sl),*) e0rho
        endif
      enddo

!------------------------------------------------------------------------
!..    read hybrid card:
      hybrid=0.d0
      do i=1,lauf
        if (in(i)(1:7).eq.'HYBRID:') then
           hybrid=1.d0
        endif
      enddo

!------------------------------------------------------------------------
!..    read so card:
      so=.false.
      do i=1,lauf
        if (in(i)(1:3).eq.'SO:') then
           so=.true.
        endif
      enddo

!------------------------------------------------------------------------
!..    read sign card:
      signswap=.false.
      do i=1,lauf
        if (in(i)(1:9).eq.'SIGNSWAP:') then
           signswap=.true.
        endif
      enddo
!------------------------------------------------------------------------

!..    read freeze card:
      freeze=.false.
      do i=1,lauf
        if (in(i)(1:7).eq.'FREEZE:') then
           freeze=.true.
        endif
      enddo

!--------------------------------------------------------------
!.....read number of parameters for two sets and coupling:
      st=0
      npar=0
      
      !Default if card is not set
      do j=1, npblocks
         pst(2,j)=0
      enddo
      !reading cards for the nuber of parameters'
      do i=1,lauf
         do j=1, npblocks
            if (in(i)(1:nkey(1,j)).eq.key(1,j)) then
               if(dbg) write(6,*) key(1,j),' read'
               read (in(i)(nkey(1,j)+1:sl),*) pst(2,j)
               st=st+1
            endif
         enddo
      enddo
      if(dbg) write(6,*) 
     &     'Ready with reading numbers of Parameters'

!..   compute total number of parameters:
      do i=1, npblocks
         npar=npar + pst(2,i)
      enddo


!..   determine start and end points of parameter blocks:
      pst(1,1)=1                ! 1 = start of block
      do i=2,npblocks
         pst(1,i)= pst(1,i-1)+pst(2,i-1)
      enddo

!.....initialize parameters:
      do i=1, npar
         par_(i)=0.d0
      enddo

      write (6,220) npar
 220  format('Total number of parameters:', i4)
      parn=npar




!----------------------------------------------------------------
!.....read parameters:
      do i=1, lauf
         do j=1, npblocks
            if(in(i)(1:nkey(2,j)).eq.key(2,j)) then
               if(dbg) write(6,*) key(2,j),' read'
               if (pst(2,j).gt.0) then
                  call longreal(in(i),sl,nkey(2,j),par_,pst(1,j)
     &                 ,pst(2,j),pcount)
                  if (pcount.lt.pst(2,j)) then
                     write(6,*) 'Error reading', key(2,j)
                     stop
                  endif
               endif
            endif
         enddo
      enddo  

      if(dbg) then
         do i=1,npar
            write(6,'(''(P(''i4'')=''f12.5)') i,par_(i)
         enddo
      endif
  
      write (*,*) 'Parameters read'

!-------------------------------------------------------------------
!....read active parameters (ma) (1 or -1 = optimize, 0 = don't optimize)
!....default is ma(i)=1, all are optimized

      do i=1,npar
         ma(i)=2
      enddo

      st=0
      do i=1, lauf
         do j=1, npblocks
            if (in(i)(1:nkey(3,j)).eq.key(3,j)) then
               if (pst(2,j).gt.0) then
                  call longint(in(i),sl,nkey(3,j),ma,pst(1,j),pst(2,j),
     $                 pcount)
                  if (pcount.lt.pst(2,j)) then
                     write(6,*) 'Error reading ',
     &                    key(3,j), pcount, pst(2,j)
                     stop
                  endif
               endif
            endif
         enddo
      enddo  

      
!     check if everything is ok
      j=0
      do i=1, npar
         if (ma(i).eq.2) then 
            write (*,*) 'not all act set', i
            stop
         endif
      enddo


!-------------------------------------------------------------------
!....read spreads for each parameter block:

!..   setting the spreads for each parameter:
      do i=1,npar
        spread(i)=gspread
      enddo

!----------------------------------------------
!     read specific spreads for parameter blocks

      do i=1, lauf
         do j=1, npblocks
            if (in(i)(1:nkey(4,j)).eq.key(4,j)) then
               if (pst(2,j).gt.0) then
                  call longreal(in(i),sl,nkey(4,j),spread,pst(1,j),
     &                 pst(2,j),pcount)
                  if (pcount.lt.pst(2,j)) then
                     write(6,*) 'Error reading', key(4,j)
                     stop
                  endif
               endif
            endif
         enddo
      enddo

!****************************************************************************
!.....check for a checkpoint file. If a checkpoint exists, use that data
!     and ignore the selections done before, except when the file is 
!     labeled by ",new":
      iter=1
      lchk=.false.
      new=.false.
      do i=1,lauf
        if (in(i)(1:4).eq.'CHK:') then
          if (i+1.gt.st) st=i+1
	  lchk=.true.
	  do j=1,sl-3
	     if (in(i)(j:j).eq.',') then
	      ntst=in(i)(j+1:j+3)
	      call lcap(ntst,3)
              if (ntst.eq.'NEW') then
	        read (in(i)(5:j-1),*) chkpnt
	        new=.true.
	        goto 500
	      endif
	     endif
	  enddo
500       continue
          if (.not.new) read(in(i)(5:sl),*) chkpnt
	  goto 501
        endif
      enddo
501   continue

      chk=.false.
      inquire (file=chkpnt, exist=chk)
505   format('Warning: No checkpoint file found, trying to run anyway...
     $ (creating new checkpoint!)')
506   Format('Creating new checkpoint file on request')
      rstrt=0
      if (lchk) then
        if (.not.chk) then
          write(6,505)
	  rstrt=1
	else if (new) then
	  write(6,506)
	  rstrt=1
        else
           reader='read'
           call chkfile(chkpnt,par_,npar,ma,seed,gtype,nset,iter,
     $               old,reader)
          rstrt=2
        endif
      endif	

            
      write(*,*) 'test ok'
!---------------------------------------------------------------
!..   read number of data points for each block.
!     this helps to create plot files later!



      write(*,*) 'before NPOINTS:'

!..   reads number of data points per set
      do i=1,lauf
        if (in(i)(1:8).eq.'NPOINTS:') then
            if (sets.gt.0) then
	      print*, 'before entering longint'
              call longint(in(i),sl,8,ndata,1,sets,pcount)
	      print*, 'after longint'
              if (pcount.lt.sets) then
                write(6,*) 'Error reading ', 'NPOINTS:'
                write(6,*) 'index is: ', i
                stop
              endif

            else
               write (*,*) 'sets not set. stop.' 
               stop
            endif 
        elseif (in(i)(1:6).eq.'COORD:') then
            if (sets.gt.0) then
              call longint(in(i),sl,6,coord,1,sets,pcount)
              if (pcount.lt.sets) then
                write(6,*) 'Error reading ', 'COORD:'
                stop
              endif
            else
               write (*,*) 'sets not set. stop.' 
               stop
            endif
        endif
      enddo	

!----------------------------------------------------------------
!.. read keyword to check if symmetrization is necessary
!..   0: no transformation
!..   1: transform from symmetrized primitive to symm. morse
!..   2: transform from true primitive to symmetrized morse
      do i=1,lauf
         if (in(i)(1:10).eq.'SYMMCOORD:') then
              call longint(in(i),sl,10,symm,1,sets,pcount)
              if (pcount.lt.sets) then
                write(6,*) 'Error reading ', 'SYMMCOORD:'
                stop
              endif

!            read (in(i)(11:sl),*) (symm(k), k=1,sets) 
            goto 676          
         endif
      enddo
676   continue
      k=0
      do i=1,sets
        if (symm(i).eq.0) k=k+1
        if (symm(i).eq.1) k=k+1
        if (symm(i).eq.2) k=k+1
        if (symm(i).eq.3) k=k+1
      enddo
      if (k.ne. sets) then
        write(6,*) 'Error in SYMMCOORD: non valid switches'
        stop
      endif

!----------------------------------------------------------------
!..   read parameters for coordinate transformation:
      write(6,*) 'reading Morse parameters now...'
      chk=.false.
      nmorse=0
      do i=1,lauf
        if (in(i)(1:7).eq.'NMORSE:') then
c     number of different Morse coordinates
	    read (in(i)(8:sl),*) nmorse
	    chk=.true.
            st = i+1 
        endif
      enddo 
      if (nmorse.gt.mxmrse) then
        write(6,'(''Error: too many Morse coordinates!'')')
	stop
      endif	

      do i=1,mxmrse
        alpha(i)=0.d0
	re(i)=0.d0
      enddo	

      if (chk) then
        do i=st,lauf
          if (in(i)(1:6).eq.'ALPHA:') then
	      read (in(i)(7:sl),*) (alpha(j), j=1,nmorse)
	      chk=.false.
          endif
        enddo
	if (chk) then
	  write(6,'(''Error: No alpha for Morse specified!'')')
	  stop
	endif
	chk=.true.
        do i=st,lauf
          if (in(i)(1:3).eq.'RE:') then
	      read (in(i)(4:sl),*) (re(j), j=1,nmorse)
	      chk=.false.
          endif
        enddo 
	if (chk) then
	  write(6,'(''Error: No re for Morse specified!'')')
	  stop
	endif
      endif
      write(6,'(''Morse parameters:'')')
      write(6,'(''alpha:'',12f12.8)') (alpha(i), i=1,nmorse)
      write(6,'(''re:   '',12f12.8)') (re(i), i=1,nmorse)


!***************************************************************************
!  start reading data block:
!***************************************************************************

!     find the data point in the file
      do i=1,lauf
        if (in(i)(1:5).eq.'DATA:') then
            st = i+1 
        endif
      enddo 
      write(6,*) 'start reading data points now...'

!***************************************************************************
c     read active cards for states
       do i=1,lauf
        if (in(i)(1:9).eq.'ACTSTATE:') then
           call longint(in(i),sl,9,wtact,1,nstat,pcount)
           if (pcount.lt.nstat) then
              write(6,*) 'ERROR: Not enough active cards for the states'
              stop
           endif
c     loop over all weights in steps of nstat
           do j=1,ntot*mx,nstat
c     loop over all states
              do k=1,nstat
                 if(wtact(k).eq.0) then
                    wt(j+k-1)=0.d0
                 else if(wtact(k).eq.1) then 
                    wt(j+k-1)=1.d0
                 else
                    write(6,*) 'ERROR: Invalid active card for State'
                 endif
              enddo
           enddo
        endif
      enddo 

!***************************************************************************
!.....read x,y value pairs and respective weights
      ii=1
      jj=ndata(ii)*(nstat+nn) !changed !!!
      num=0
      do i=st,lauf
        if (in(i)(1:3).eq.'WT:') then
          read (in(i)(4:sl),*) wt(num)   !weight has to follow value!
        else
          num=num+1
          read (in(i),*) y(num), (x(j,num), j=1,qn)
!R	  wt(num)=1.d0
	  if (num.ge.jj) then
	    write(6,234) ii
	    ii=ii+1
	    write(6,'(i8,'' data points read'')') jj
	    jj=jj+ndata(ii)*(nstat+nn) !changed !!!
!	    write(6,'(a120)') in(i)
!	    write(6,*)
	  endif  
        endif
      enddo
234   format('Set',i4,':')

1010  format('Error: Number of data points on input inconsistent with ND
     $ATA card!')
1011  format('NDATA information:',i8)
1012  format('Number of points: ',i8)

!      STOP 'TEST1'

!      call checkingparameter(pst)

      end subroutine

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE CAPITAL(IN,STR,LAUF,MMAX,SL)
C%
C% Changes all Characters of a line to Capital Letters, removes lines
C% beginning vith '!' and cuts leading spaces.
C%
C% variables:
C% in:   input file
C% str:  current line of inputfile
C% lauf: number of current line
C% mmax: maximum number of input lines
C% sl:   maximum line length
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine capital(in,str,lauf,mmax,sl)
      implicit none

      integer i !running index

      integer mmax           !maximum number of lines
      integer lauf           !current line
      character in(mmax)*(*) !whole input string

      integer j         !real start of the line
      integer sl        !Length of one line
      character str*(*) !manipulated line (one line of input)

c     ignor empty lines
      if (str.eq.'') return

c     search for first non spacing character
      j=0
      do i=1,sl
        if (str(i:i).ne.' ') then
          j=i-1
          goto 10
        endif
      enddo
c     Write important part of the line on str
10    do i=1,sl-j
        str(i:i)=str(i+j:i+j)
      enddo
c     fill the remaining characters with spaces
      do i=sl-j+1,sl
        str(i:i)=' '
      enddo

c     Ignor lines with ! at the beginning
      if (str(1:1).eq.'!') return

c     Turns every letter into capital letter
      lauf=lauf+1
      do i=1,sl
        in(lauf)(i:i)=str(i:i) 
        if (str(i:i).eq.'a') in(lauf)(i:i)='A'
        if (str(i:i).eq.'b') in(lauf)(i:i)='B'
        if (str(i:i).eq.'c') in(lauf)(i:i)='C'
        if (str(i:i).eq.'d') in(lauf)(i:i)='D'
        if (str(i:i).eq.'e') in(lauf)(i:i)='E'
        if (str(i:i).eq.'f') in(lauf)(i:i)='F'
        if (str(i:i).eq.'g') in(lauf)(i:i)='G'
        if (str(i:i).eq.'h') in(lauf)(i:i)='H'
        if (str(i:i).eq.'i') in(lauf)(i:i)='I'
        if (str(i:i).eq.'j') in(lauf)(i:i)='J'
        if (str(i:i).eq.'k') in(lauf)(i:i)='K'
        if (str(i:i).eq.'l') in(lauf)(i:i)='L'
        if (str(i:i).eq.'m') in(lauf)(i:i)='M'
        if (str(i:i).eq.'n') in(lauf)(i:i)='N'
        if (str(i:i).eq.'o') in(lauf)(i:i)='O'
        if (str(i:i).eq.'p') in(lauf)(i:i)='P'
        if (str(i:i).eq.'q') in(lauf)(i:i)='Q'
        if (str(i:i).eq.'r') in(lauf)(i:i)='R'
        if (str(i:i).eq.'s') in(lauf)(i:i)='S'
        if (str(i:i).eq.'t') in(lauf)(i:i)='T'
        if (str(i:i).eq.'u') in(lauf)(i:i)='U'
        if (str(i:i).eq.'v') in(lauf)(i:i)='V'
        if (str(i:i).eq.'w') in(lauf)(i:i)='W'
        if (str(i:i).eq.'x') in(lauf)(i:i)='X'
        if (str(i:i).eq.'y') in(lauf)(i:i)='Y'
        if (str(i:i).eq.'z') in(lauf)(i:i)='Z'
        if (i.le.3) cycle
	if (in(lauf)(i-3:i).eq.'CHK:') then
	  in(lauf)(i+1:sl)=str(i+1:sl)
	  return
	endif
      enddo

      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE LCAP(STR,N)
C%
C% Changes all letters of a line to capital latters.
C%
C% variables:
C% str: line which should be manipulated
C% n:   length of the line
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine lcap(str,n)
      implicit none
      integer i, n
      character str*(*), dum*750

      dum=''
      do i=1,n
        dum(i:i)=str(i:i)
        if (str(i:i).eq.'a') dum(i:i)='A'
        if (str(i:i).eq.'b') dum(i:i)='B'
        if (str(i:i).eq.'c') dum(i:i)='C'
        if (str(i:i).eq.'d') dum(i:i)='D'
        if (str(i:i).eq.'e') dum(i:i)='E'
        if (str(i:i).eq.'f') dum(i:i)='F'
        if (str(i:i).eq.'g') dum(i:i)='G'
        if (str(i:i).eq.'h') dum(i:i)='H'
        if (str(i:i).eq.'i') dum(i:i)='I'
        if (str(i:i).eq.'j') dum(i:i)='J'
        if (str(i:i).eq.'k') dum(i:i)='K'
        if (str(i:i).eq.'l') dum(i:i)='L'
        if (str(i:i).eq.'m') dum(i:i)='M'
        if (str(i:i).eq.'n') dum(i:i)='N'
        if (str(i:i).eq.'o') dum(i:i)='O'
        if (str(i:i).eq.'p') dum(i:i)='P'
        if (str(i:i).eq.'q') dum(i:i)='Q'
        if (str(i:i).eq.'r') dum(i:i)='R'
        if (str(i:i).eq.'s') dum(i:i)='S'
        if (str(i:i).eq.'t') dum(i:i)='T'
        if (str(i:i).eq.'u') dum(i:i)='U'
        if (str(i:i).eq.'v') dum(i:i)='V'
        if (str(i:i).eq.'w') dum(i:i)='W'
        if (str(i:i).eq.'x') dum(i:i)='X'
        if (str(i:i).eq.'y') dum(i:i)='Y'
        if (str(i:i).eq.'z') dum(i:i)='Z'
      enddo
      str(1:n)=dum(1:n)

      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% FUNCTION CLEN(STR,SL)
C%
C% Function to test how many entries are on one line.
C%
C% variables
C% str: line, which should be checked
C% sl:  length of the line
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function clen(str,sl)
      implicit none
      integer clen, i, j, sl
      character str*(sl)

      clen=0
      j=0
      do i=sl,1,-1
        if ((str(i:i).ne.' ').and.(j.eq.0)) then
	  clen=clen+1
	  j=1
	endif
	if (str(i:i).eq.' ') j=0
      enddo

      end

!-------------------------------------------------------------------------
!  subroutine to handle the checkpoint file:
      subroutine chkfile(chkpnt,par_,npar,ma,seed,gtype,nset,iter,
     $               old,str)
      implicit none
      integer i, j, npar, ma(npar), seed, gtype, nset, iter
      double precision par_(npar), old
      character chkpnt*80, str*10

      if (str(1:4).eq.'read') then
        open (10, file=chkpnt, form='unformatted')
        read (10) npar, seed, gtype, iter, old
	read (10) (ma(i), i=1,npar)
	do i=1,nset
	  read (10) (par_(j), j=(i-1)*npar+1,i*npar)
	enddo
	close (10)
      else if (str(1:5).eq.'write') then
        open (10, file=chkpnt, form='unformatted')
        write (10) npar, seed, gtype, iter, old
	write (10) (ma(i), i=1,npar)
	do i=1,nset
	  write (10) (par_(j), j=(i-1)*npar+1,i*npar)
	enddo
	close (10)
      endif

      end



!-------------------------------------------------------------------------
!..   routine to read multiple-line input for parameters et!.
      subroutine longreal(in,sl,nkey,par_,pst1,pst2,pcount)
      implicit none
      integer i, k, ii, jj
      integer line, sta, last, pcount
      integer sl        ! defined length of character strings
      integer nkey      ! length of keyword
      integer pst1,pst2 ! starting address_ and length of parameter block
      real*8 par_(*)     ! parameter set
      character*(sl) in(*) ! current input line
      integer clen
      logical cont ! flag whether or not current line will be continued


!..   set starting point for reading for first line such that keyword 
c     is skipped:
      sta=nkey+1

!..   loop over lines:
      pcount=0                   ! counting number of read parameters
      ii=pst1                   ! beginning of parameter block
      do line=1,10000

!     print*, in(line)(1:130)
!     print*,'ii:', ii

c     first check whether current line will be continued:
         cont=.false.
         do i=sl,1,-1
            if (in(line)(i:i).ne.' ') then
               last=i
               if (in(line)(i:i).eq.'&') then 
                  cont=.true.
                  last=i-1
               endif
               goto 100
            endif  
         enddo
 100     continue	


c     now read entries from current line

         jj=clen(in(line),sl)
         if (cont) jj=jj-1      ! remove '&' at end of line
         if (line.eq.1) jj=jj-1 ! remove keyword from beginning of line
!     print*, 'jj:', jj
         pcount=pcount+jj
         read(in(line)(sta:last),*) (par_(k), k=ii, ii+jj-1)
         read(in(line)(sta:last),*) (par_(k), k=ii, min(ii+jj-1,ii+pst2))
!     print*, 'end:', ii+jj-1
         ii=ii+jj
         sta=1

         if (.not.cont) return

      enddo

      end

!------------------------------------------------------------------------
c     routine to read multiple-line input for parameters etc.
      subroutine longint(in,sl,nkey,par_,pst1,pst2,pcount)
      implicit none
      integer i, k, ii, jj
      integer line, sta, last, pcount
      integer sl        ! defined length of character strings
      integer nkey      ! length of keyword
      integer pst1,pst2 ! starting address_ and length of parameter block
      integer par_(*) ! parameter set
      character*(sl) in(*) ! current input line
      integer clen

      logical cont  ! flag whether or not current line will be continued

!      print*, 'longint entered'
!      print*, in(1)

c     set starting point for reading for first line such that 
c     keyword is skipped:
      sta=nkey+1
!      print*, 'sta:', sta

c     loop over lines:
      pcount=0                   ! counting number of read parameters
      ii=pst1                    ! beginning of parameter block
!      print*, 'ii:', ii
      do line=1,10000

c     first check whether current line will be continued:
        cont=.false.
        do i=sl,1,-1
	  if (in(line)(i:i).ne.' ') then
	    last=i
            if (in(line)(i:i).eq.'&') then 
	      cont=.true.
	      last=i-1
	    endif
	    goto 100
	  endif  
        enddo
100     continue	


!..     now read entries from current line

	jj=clen(in(line),last)
!	if (cont) jj=jj-1 ! remove '&' at end of line
	if (line.eq.1) jj=jj-1 ! remove keyword from beginning of line
!	print*, 'number of entries jj:', jj
	pcount=pcount+jj
!        write(6,*) 'minimum number of entries:'
!        write(6,*) 'min(ii+jj-1,pst2)', min(ii+jj-1,pst2)
!	write(6,'(a120)') in(i)
        read(in(line)(sta:last),*) (par_(k), k=ii, min(ii+jj-1,ii+pst2))
!        write(6,'(9x,10i4,x,10i4)')(par_(k), k=ii, min(ii+jj-1,ii+pst2))
	ii=ii+jj  ! moving to next section of parameter block
	sta=1

	if (.not.cont) return

      enddo

      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE CONVENTION(Y,NUM)
C%
C% Checks and sets some conventions for genetic. The sign of each
C% biggest ci coefficient is set positiv
C%
C% Variables:
C% y:   datavalues
C% num: number of datavalues
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine convention(y,num)
      implicit none

      include 'states.incl'

      integer num             !number of datapoints
      double precision y(num) !datavalues

      integer i,j

      logical dbg !more outpu for debugging
      
      dbg=.true.
      dbg=.false.

      do i=1,num,ntot
         if(dbg) write(6,'(''Data set: ''i6)') (i/ntot)+1
         do j=nstat,ntot-1,nstat
            if(dbg) write(6,'(''read CI: ''3f9.4)')
     $           y((i+j):(i+j+nstat-1))
            call correctPhase(y(i+j),nstat)
            if(dbg) write(6,'(''corrected CI: ''3f9.4)')
     $           y(i+j:i+j+nstat-1)
         enddo
      enddo

      do i=1,num,ntot
         if(dbg) write(6,'(''Data set: ''i6)') (i/ntot)+1
         call sortInput(y(i+nstat:i+nstat*(nstat+1)),nstat)
      enddo

      end subroutine

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE SORTINPUT(V,N)
C%
C% Sets the molpro phase convention.
C%
C% variables:
C% v: ci-vector
C% n: length of the vector
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine sortInput(v,nstat)

      integer nstat                     !nnumber of states
      double precision v(nstat,nstat)   !ci-vector
      double precision rot(nstat,nstat) !sorting matrix

      double precision result(nstat,nstat)
      
      integer i,j,k !running indeces

c     declare sorting matrix
      do i=1,nstat
         do j=1,nstat
            rot(i,j)=0.d0
            result(i,j)=0.d0
         enddo
      enddo
      rot(1,2)=1.d0
      rot(2,3)=1.d0
      rot(3,1)=1.d0
      rot(4,4)=1.d0
      rot(5,5)=1.d0
      rot(6,6)=1.d0

c     multipli cis with it
      do i=1,nstat
         do j=1,nstat
            do k=1,nstat
               result(i,j)=result(i,j)+rot(k,i)*v(k,j)
            enddo
         enddo
      enddo
      
      do i=1,nstat
         do j=1,nstat
            v(i,j)=result(i,j)
         enddo
      enddo
      
      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE CORRECTPHASE(V,N)
C%
C% Sets the molpro phase convention.
C%
C% variables:
C% v: ci-vector
C% n: length of the vector
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine correctPhase(v,n)
      implicit none

      integer n              !length of the vector
      double precision v(n)  !ci-vector

      integer i !running index
      
      double precision biggest !biggest ci-value
      double precision sig     !sign of the biggest ci-value

      double precision sign !internal function      

c     searching biggest ci
      biggest=v(1)
      do i=2,n
         if(abs(biggest).lt.abs(v(i))) biggest=v(i)
      enddo

c     determing right sign
      sig=1.d0
      sig = sign(sig,biggest)

c     setting the right signs
      do i=1,n
         v(i)=sig*v(i)
      enddo

      end subroutine

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine calcSOData(num,y,npar,ma,pst,p,q,somat)
      implicit none

      include 'params.incl' !this includes the size definitions
      include 'states.incl'    !this includes number of states and modes

      double complex somat(2*nstat,2*nstat,mx) !so matrecies
      integer num
      integer numso
      double precision y(2*ntot*mx)
      double precision yso(2*nstat*mx)

      integer npar
      integer ma(npar)
      integer pst(2,np)
      double precision v(ntot)
      double precision q(qn*ntot*mx)
      double precision p(px*nx)
      logical skip
      double complex z(2*nstat,2*nstat)
      double precision inter

      double complex re

      integer i,j

      double precision cm
      parameter (cm=219474.69d0)

      skip=.false.
      re=(1.d0,0.d0)

      numso=(num/ntot)*2*nstat

      do i=1,num/ntot
         call adia(v,i,q((i-1)*qn*ntot+1),p,npar,pst,skip,
     &        y((i-1)*ntot)+1) 
         inter = v(1)
         do j=0,nstat-1
            v(nstat-j)=(v(nstat-j)-v(1))*cm
         enddo
         do j=1,nstat
            somat(j,j,i) = somat(j,j,i)+ re*v(j)
            somat(j+nstat,j+nstat,i) = somat(j+nstat,j+nstat,i)
     &           + re*v(j)
         enddo
         call hermdiag(6,6,somat(1,1,i)
     &        ,yso((i-1)*2*nstat+1),0,z)
         do j=1,2*nstat
            yso((i-1)*2*nstat+j)=yso((i-1)*2*nstat+j)/cm+inter
         enddo
      enddo

      num=numso
      do i=1,num
         y(i)=yso(i)
      enddo

      end
