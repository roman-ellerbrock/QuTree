cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Random generator From Marius Lewerenz
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gautrg(gran,nran,iseed,iout)
c
c      vectorized portable unit standard deviation gaussian random
c      number generator using the box-mueller transformation method.
c      this method is faster than the central limit method, when the
c      uniform random numbers are comparatively expensive.
c      version working on a square with sine and cosine functions. the
c      same sequence is produced independent of the sequence of calls
c      if no calls to other vranf based prn generators are made.
c
c      gran  : vector of length nran returning the random numbers
c      nran  : number of desired random numbers, with nran=0 and iseed
c              not 0 only generator initialization is done. no action
c              when both are zero.
c      iseed : if not 0, integer to start generator. use 0 to continue
c              a previously used/initialized random sequence, unchanged
c              on output.
c      iout  : unit number for messages.
c
c      times for the generation of 10**6 prn's:
c      rs6000/350 41 mhz    2.31 s 
c      r8000      75 mhz    1.43 s
c      axp21064  200 mhz    1.20 s
c      axp21164  250 mhz    0.572 s
c      t3e-900   450 mhz    0.52 s
c      cray-t90  450 mhz    0.0728 s   1.28*10**7 prn/s
c
c      subroutines called: vranf, r1mach
c      m. lewerenz 6/may/90, modified 17/jun/91, mar/95
c
      implicit real*8 (a-h,o-z)
      parameter (two=2.d0,twom=-two,one=1.d0,npoly=11)
      dimension gran(nran),ccf(npoly),scf(npoly)
      save isave,gsave,tiny,twopi,pi4,ccf,scf
      data isave/-1/
C
C      POLYNOMIAL FROM CHEBYSHEV APPROXIMATION ON [ 0.000, 0.790]
C      FOR COS(X) WITH ABSOLUTE ERROR LESS THAN 0.2220E-14
C
      DATA ccf/ 0.9999999999999986D+00, 0.6612476846390664D-13,
     #         -0.4999999999989523D+00,-0.5434658088910759D-10,
     #          0.4166666737609693D-01,-0.4648977428396692D-08,
     #         -0.1388871052129944D-02,-0.4228394738587799D-07,
     #          0.2486361945804866D-04,-0.5317743184071916D-07,
     #         -0.2539224676809412D-06/
C
C      POLYNOMIAL FROM CHEBYSHEV APPROXIMATION ON [ 0.000, 0.790]
C      FOR SIN(X) WITH ABSOLUTE ERROR LESS THAN 0.2220E-14
C
      DATA scf/-0.9443414574112906D-15, 0.1000000000000244D+01,
     #         -0.1224236196202217D-10,-0.1666666664242968D+00,
     #         -0.2471495821870120D-08, 0.8333348067492644D-02,
     #         -0.5482536616811601D-07,-0.1982815612039858D-03,
     #         -0.2017619095413939D-06, 0.2948964053761139D-05,
     #         -0.1051448397925916D-06/
c
      if(isave.lt.0) then
        isave=0
        tiny=r1mach(0)**2
        pi4=atan(one)
        twopi=pi4*8
      end if
      if(iseed.ne.0) call vranf(gran,0,iseed,iout)
c
      if(nran.gt.0) then
        newran=nran-isave
        if(isave.eq.1) gran(1)=gsave
        call vranf(gran(isave+1),newran,0,iout)
        do 100 i=1,newran-1,2
        fac=sqrt(twom*log(gran(isave+i)+tiny))
        x=pi4*gran(isave+i+1)
        cx=(((((((((ccf(11)*x+ccf(10))*x+ccf(9))*x+ccf(8))*x
     #                       +ccf(7))*x+ccf(6))*x+ccf(5))*x
     #                       +ccf(4))*x+ccf(3))*x+ccf(2))*x+ccf(1)
        sx=(((((((((scf(11)*x+scf(10))*x+scf(9))*x+scf(8))*x
     #                       +scf(7))*x+scf(6))*x+scf(5))*x
     #                       +scf(4))*x+scf(3))*x+scf(2))*x+scf(1)
        sxi=(two*cx)*sx
        cxi=(two*cx)*cx-one
        sxi=(two*cxi)*sxi
        cxi=(two*cxi)*cxi-one
        sxi=(two*cxi)*sxi
        cxi=(two*cxi)*cxi-one
        gran(isave+i)=fac*sxi
        gran(isave+i+1)=fac*cxi
  100   continue
c
        if(mod(newran,2).eq.1) then
          call vranf(xran,1,0,iout)
          fac=sqrt(twom*log(gran(nran)+tiny))
          trig=twopi*xran
          gran(nran)=fac*sin(trig)
          gsave=fac*cos(trig)
          isave=1
        else
          isave=0
        end if
      end if
      return
      end
C
C----------------------------------------------------------------------
C
      FUNCTION R1MACH (IDUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ONE=1.D0,TWO=2.D0,HALF=0.5D0)
      SAVE ICALL,EPS
      DATA ICALL,EPS/0,ONE/
C ---------------------------------------------------------------------
C THIS ROUTINE COMPUTES THE UNIT ROUNDOFF OF THE MACHINE.
C THIS IS DEFINED AS THE SMALLEST POSITIVE MACHINE NUMBER
C U SUCH THAT 1.0 + U .NE. 1.0E0
C ---------------------------------------------------------------------
      IF(ICALL.EQ.0) THEN
        ICALL=1
        U = ONE
  10    U = U*HALF
        COMP = ONE + U
        IF(COMP .NE. ONE) GO TO 10
        EPS = U*TWO
      END IF
      R1MACH = EPS
      RETURN
      END
c
c----------------------------------------------------------------------
c
      subroutine vcopy(v1,v2,nv)
c
c      copies vector v1 into v2.
c      subroutines called : none                           m. lewerenz 
c
      implicit real*8 (a-h,o-z)
      dimension v1(nv),v2(nv)
c
      if(nv.gt.0) then
        do 10 i=1,nv
        v2(i)=v1(i)
   10   continue
      end if
      return
      end

c-----------------------------------------------------------------------
c-------------------- ranf/vranf uniform random package ----------------
c-----------------------------------------------------------------------
c
      subroutine vranf(ranv,nran,iseed,iout)
c
c      machine independent portable uniform random number generator for
c      the interval [0,1) based on a floating point subtractive lagged
c      fibonacci method similar to the feedback shift register method
c      proposed by kirkpatrick/stoll. subtractive (or additive) version
c      gives much better prn quality than the original xor operation.
c      an additive variant is given in v r a n f a.
c
c      v r a n f  and r a n f use the same method and can be used to
c      work together on the same random sequence. either of them can be
c      used for initialization. the state of the generator can be saved
c      or retrieved for restart with routine r s a v e f.
c ----------------------------------------------------------------------
c
c      ranv  : vector of length nran for random numbers; output
c      nran  : number of desired random numbers. nran=0 and iseed.ne.0
c              -> only generator initialization. no action when both 
c              are zero.; input
c      iseed : if not 0, integer to start generator. use 0 to continue
c              a previously used/initialized random sequence, unchanged
c              on output; input
c      iout  : unit number for error messages. silent for iout.le.0
c
c ----------------------------------------------------------------------
c      method:
c      x(k+np)=x(k)-x(k+np-nq), initial x array generated by xuinit.
c      ieee standard requires real*8 to have at least 48 mantissa bits.
c      with nbit=48 this generator is entirely machine independent
c      and will always give the same random sequence. you can change
c      the period of the generator by setting a different nbit value or
c      changing np and nq appropriately.
c
c ----------------------------------------------------------------------
c      this is a floating implementation of a generator described in:
C      M.H. KALOS, P.A. WHITLOCK, MONTE CARLO METHODS, APPENDIX,
C                  WILEY 1986
C      D.W. HEERMANN, COMPUTER SIMULATION METHODS, 2ND ED.,SPRINGER 199
C                  APPENDIX A1
c      d. stauffer, f.w. hehl, v. winkelmann, j.g. zabolitzky,
c                  computer simulation & computer algebra, section 2.2
c
c      original references:
c      s. kirkpatrick, e.p. stoll, j. comput. linearizedLeaves_. 40, 517 (1981)
c      r.c. tausworthe, random numbers generated by linear recurrence
c                  modulo 2, math. comp. 19, 201 (1965)
c      t.g. lewis, w.h. payne, generalized feedback shift register
c                  pseudorandom number algorithm, j. acm 20, 456 (1973)
c
c ----------------------------------------------------------------------
c      other (np,nq) values: (17,5), (250,103), (521,158), (1279,418),
c                            (2281,715), (4423,1393), (1279,1063)
c          ref.: Bhanot et al., linearizedLeaves_. rev b 33, 7841 (1986);
c                Zierler, inf. control 15, 67 (1961)
c ----------------------------------------------------------------------
c      alternative additive formulation bypassing if statements:
c      temp=x(k)+x(k+np-nq)
c      x(k)=temp-float(int(temp))
c
c      alternative subtractive formulation bypassing if statements:
c      temp=x(k)-x(k+np-nq)
c      x(k)=(temp+one)-float(int(temp+one))
c ----------------------------------------------------------------------
c
c      timing in s for 10**7 random numbers:
c
c      machine           mhz   unroll(4)     no unrolling
c      ibm rs6000/350     41      2.67            3.9
c      dec axp3800       200      1.17
c      dec axp3600       150      1.35
c      dec alpha 21164   250      0.45
c      t3e alpha 21164   450      0.40
c      sun ultra1        ???      1.38
c      ibm 3090vf         58      ---        2.47(vec), 4.57(sc)
c      cray t90          450      0.31      ranf() takes 0.13 s
c      
c      unrolling to depth 6 gives a slight speed increase.
c ----------------------------------------------------------------------
c      subroutines called: xuinit,errprt    m. lewerenz may/91 & nov/93
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
      parameter (nratio=np/nq,nexec=4,mroll=4,zero=0.d0,one=1.d0)
      dimension ranv(nran)
      common /doctrl/ nroll
c
c      table initialization by xuinit
c
      if(iseed.ne.0) then
        call xuinit(x,np,nq,0,nexec,iseed,init,last,iout)
      end if
c
c      fibonacci generator updates elements of x in a cyclic fashion
c      and copies them into ranv in blocks of max. length np.
c      loop split into chunks of max. length nq to avoid recurrence.
c      unrolling improves performance on superscalar machines.
c
      if(nran.gt.0) then
        if(init.ne.0) then
          j=0
          left=nran
   10     continue
          if(nroll.gt.1) then
            loop=mod((min(nq,left+last)-last),mroll)
          else
            loop=min(nq,left+last)-last
          end if
convex, cray, and ibm directives
c$dir no_recurrence
cdir$ ivdep
cibmdir ignore recrdeps
          do 500 i=last+1,last+loop
          x1=x(i)-x(i+np-nq)
          if(x1.lt.zero) x1=x1+one
          x(i)=x1
          j=j+1
          ranv(j)=x1
  500     continue
          if(nroll.gt.1) then
            do 501 i=last+loop+1,min(nq,left+last),mroll
            x1=x(i)-x(i+np-nq)
            x2=x(i+1)-x(i+1+np-nq)
            x3=x(i+2)-x(i+2+np-nq)
            x4=x(i+3)-x(i+3+np-nq)
            if(x1.lt.zero) x1=x1+one
            if(x2.lt.zero) x2=x2+one
            if(x3.lt.zero) x3=x3+one
            if(x4.lt.zero) x4=x4+one
            x(i)=x1
            x(i+1)=x2
            x(i+2)=x3
            x(i+3)=x4
            ranv(j+1)=x1
            ranv(j+2)=x2
            ranv(j+3)=x3
            ranv(j+4)=x4
            j=j+4
  501       continue
          end if
c
          if(last.lt.nratio*nq) then
            do 650 k=1,nratio-1
            limit=min((k+1)*nq,left+last)
            if(nroll.gt.1) then
              loop=mod((limit-max(k*nq,last)),mroll)
            else
              loop=limit-max(k*nq,last)
            end if
convex, cray, and ibm directives
c$dir no_recurrence
cdir$ ivdep
cibmdir ignore recrdeps
            do 600 i=max(k*nq,last)+1,max(k*nq,last)+loop
            x1=x(i)-x(i-nq)
            if(x1.lt.zero) x1=x1+one
            x(i)=x1
            j=j+1
            ranv(j)=x1
  600       continue
            if(nroll.gt.1) then
              do 601 i=max(k*nq,last)+loop+1,limit,mroll
              x1=x(i)-x(i-nq)
              x2=x(i+1)-x(i+1-nq)
              x3=x(i+2)-x(i+2-nq)
              x4=x(i+3)-x(i+3-nq)
              if(x1.lt.zero) x1=x1+one
              if(x2.lt.zero) x2=x2+one
              if(x3.lt.zero) x3=x3+one
              if(x4.lt.zero) x4=x4+one
              x(i)=x1
              x(i+1)=x2
              x(i+2)=x3
              x(i+3)=x4
              ranv(j+1)=x1
              ranv(j+2)=x2
              ranv(j+3)=x3
              ranv(j+4)=x4
              j=j+4
  601         continue
            end if
  650       continue
          end if
c
          limit=min(np,left+last)
          if(nroll.gt.1) then
            loop=mod((limit-max(nratio*nq,last)),mroll)
          else
            loop=limit-max(nratio*nq,last)
          end if
convex, cray, and ibm directives
c$dir no_recurrence
cdir$ ivdep
cibmdir ignore recrdeps
          do 700 i=max(nratio*nq,last)+1,max(nratio*nq,last)+loop
          x1=x(i)-x(i-nq)
          if(x1.lt.zero) x1=x1+one
          x(i)=x1
          j=j+1
          ranv(j)=x1
  700     continue
          if(nroll.gt.1) then
            do 701 i=max(nratio*nq,last)+loop+1,limit,mroll
            x1=x(i)-x(i-nq)
            x2=x(i+1)-x(i+1-nq)
            x3=x(i+2)-x(i+2-nq)
            x4=x(i+3)-x(i+3-nq)
            if(x1.lt.zero) x1=x1+one
            if(x2.lt.zero) x2=x2+one
            if(x3.lt.zero) x3=x3+one
            if(x4.lt.zero) x4=x4+one
            x(i)=x1
            x(i+1)=x2
            x(i+2)=x3
            x(i+3)=x4
            ranv(j+1)=x1
            ranv(j+2)=x2
            ranv(j+3)=x3
            ranv(j+4)=x4
            j=j+4
  701       continue
          end if
c
          last=mod(limit,np)
          left=nran-j
          if(left.gt.0) goto 10
        else
          call errprt(iout,'vranf','incorrect initialization',-1)
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine xuinit(x,np,nq,mode,nexec,iseed,init,last,iout)
c
c      initializes a (np,nq) lagged fibonacci generator table
c      with random bits generated by a congruential generator using
c      l'ecuyers decomposition. ref.: bratley et al. p. 214
c
c      x     : vector of length np for initial random number table;
c              output
c      np,nq : parameters p and q of feedback shift register generator;
c              input
c      mode  : operation for lfg:
c              mode=<0 -> subtractive generator, mode=1 additive
c      nexec : number of warm up_ cycles for the table. nexec*nbit*np
c              random numbers are generated and discarded; input
c      iseed : integer seed for congruential generator generating
c              the bits of the initial table entries; input
c      init  : returns updated seed of congruential generator.
c              0 if table was not initialized, > 0 if ok; output
c      last  : pointer to the last used number in the table; output
c      iout  : unit number for messages, silent for iout.le.0; input
c      subroutines called : errprt           m. lewerenz mar/93, mar/98
c
      implicit real*8 (a-h,o-z)
      parameter (ia=40692,ib=52774,ic=3791,ip=2147483399)
      parameter (zero=0.d0,one=1.d0,half=0.5d0,iphalf=ip/2,nbit=48)
      logical high
      dimension x(np)
c
      if(nq.ge.np.or.iseed.eq.0) then
        call errprt(iout,'xuinit','illegal seed parameter(s)',-1)
      else
c
c      set table to zero and exercise the bit generator a little
c
        ix=iabs(iseed)
        if(ix.ne.0) then
          do i=1,np
            x(i)=zero
            k1=ix/ib
            ix=ia*(ix-k1*ib)-k1*ic
            if(ix.lt.0) ix=ix+ip
          end do
c
c      assemble np numbers with mantissa length nbit from random bits
c      'high' toggle compensates for bias from odd ip
c
          high=.true.
          do i=1,np
            add=half
            do j=1,nbit
              k1=ix/ib
              ix=ia*(ix-k1*ib)-k1*ic
              if(ix.lt.0) ix=ix+ip
              if(high) then
                if(ix.ge.iphalf) x(i)=x(i)+add
                high=.false.
              else
                if(ix.gt.iphalf) x(i)=x(i)+add
                high=.true.
              end if
              add=add*half
            end do
          end do
          if(nexec.gt.0) call xuwarm(x,np,nq,mode,nbit*nexec,iout)
        end if
        init=ix
        last=0
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine xuwarm(x,np,nq,mode,nexec,iout)
c
c      warms up_ (p,q) lagged fibonacci generators (lfg) by nexec rounds
c
c      x     : vector of length np for initial random number table;
c              output
c      np,nq : parameters p and q of feedback shift register generator;
c              input
c      mode  : operation for lfg:
c              mode=<0 -> subtractive generator, mode=1 additive
c      nexec : number of warm up_ cycles for the table. nexec*nbit*np
c              random numbers are generated and discarded; input
c      iout  : unit number for messages, silent for iout.le.0; input
c      subroutines called : errprt                   m. lewerenz mar/98
c
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0)
      dimension x(np)
c
      if(nq.ge.np.or.np.eq.0.or.nq.eq.0) then
        call errprt(iout,'xuwarm','illegal table parameter(s)',-1)
      else
c
c      exercise the generator for nexec rounds of np prn's
c      separate sections for subtractive or additive version
c
          if(mode.le.0) then
            do k=1,nexec
              do i=1,nq
                x(i)=x(i)-x(i+np-nq)
                if(x(i).lt.zero) x(i)=x(i)+one
              end do
              do i=nq+1,np
                x(i)=x(i)-x(i-nq)
                if(x(i).lt.zero) x(i)=x(i)+one
              end do
            end do
          else
            do k=1,nexec
              do i=1,nq
                x(i)=x(i)+x(i+np-nq)
                if(x(i).ge.one) x(i)=x(i)-one
              end do
              do i=nq+1,np
                x(i)=x(i)+x(i-nq)
                if(x(i).ge.one) x(i)=x(i)-one
              end do
            end do
          end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vrans(ranv,nran,iseed,iout)
c
c      floating point subtractive lagged finonacci generator for
c      uniformly distributed random numbers.
c      identical with vranf but split into two separate subroutines.
c      see there for description.
c
c      ranv  : vector of length nran for random numbers; output
c      nran  : number of desired random numbers. nran=0 and iseed.ne.0
c              -> only generator initialization. no action when both 
c              are zero.; input
c      iseed : if not 0, integer to start generator. use 0 to continue
c              a previously used/initialized random sequence, unchanged
c              on output; input
c      iout  : unit number for error messages. silent for iout.le.0
c
c      times for 10**7 prn:    sun ultra1        1.423 s
c                              cray-t90 450 mhz  0.307 s
c
c      subroutines called: xuinit, fslfg, errprt     m. lewerenz mar/98
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
      parameter (nexec=4)
      dimension ranv(nran)
c
c      table initialization by xuinit
c
      if(iseed.ne.0) then
        call xuinit(x,np,nq,0,nexec,iseed,init,last,iout)
      end if
c
c      cyclic table update and output vector generation by fslfg
c
      if(nran.gt.0) then
        if(init.ne.0) then
          call fslfg(x,np,nq,last,ranv,nran,iout)
        else
          call errprt(iout,'vrans','incorrect initialization',-1)
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vrans2(ranv,nran,a,b,iseed,iout)
c
c      floating point subtractive generator for random numbers with
c      uniform distribution on the interval [a,b).
c
c      ranv  : vector of length nran for random numbers; output
c      nran  : number of desired random numbers. nran=0 and iseed.ne.0
c              -> only generator initialization. no action when both 
c              are zero.; input
c      a,b   : interval limits for random numbers; input
c      iseed : if not 0, integer to start generator. use 0 to continue
c              a previously used/initialized random sequence, unchanged
c              on output; input
c      iout  : unit number for error messages. silent for iout.le.0
c
c      times for 10**7 prn:    sun ultra1  1.502 s
c                              cray-t90    0.355 s
c ----------------------------------------------------------------------
c      subroutines called: xuinit, fslfg2, errprt 
c                                                   m. lewerenz mar/98
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
      parameter (nexec=4)
      dimension ranv(nran)
c
c      table initialization by xuinit
c
      if(iseed.ne.0) then
        call xuinit(x,np,nq,0,nexec,iseed,init,last,iout)
      end if
c
c      cyclic table update and output vector generation by fslfg2
c
      if(nran.gt.0) then
        if(init.ne.0) then
          call fslfg2(x,np,nq,last,a,b,ranv,nran,iout)
        else
          call errprt(iout,'vrans2','incorrect initialization',-1)
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vranf2(ranv,nran,a,b,iseed,iout)
c
c      identical to vrans2, kept to resolve references to vranf2.
c      see vrans2 for description.                  m. lewerenz mar/98
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
      parameter (nexec=4)
      dimension ranv(nran)
c
c      table initialization by xuinit
c
      if(iseed.ne.0) then
        call xuinit(x,np,nq,0,nexec,iseed,init,last,iout)
      end if
c
c      cyclic table update and output vector generation by fslfg2
c
      if(nran.gt.0) then
        if(init.ne.0) then
          call fslfg2(x,np,nq,last,a,b,ranv,nran,iout)
        else
          call errprt(iout,'vranf2','incorrect initialization',-1)
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vranfx(ranv,nran,iseed,iout)
c
c      special version of vranf eliminating the store into the table
c      during the main computation. this version cannot work together
c      with ranf. 0 =< x < 1
c
c      ranv  : vector of length nran returning the random numbers
c      nran  : number of desired random numbers. with nran=0 and iseed
c              not 0 only generator initialization is done. no action
c              when both are zero.
c      iseed : if not 0, integer to start generator. use 0 to continue
c              a previously used/initialized random sequence, unchanged
c              on output
c      iout  : unit number for error messages. silent for iout.le.0
c
c      method and references : see comments in vranf
c      subroutines called: xuinit               m. lewerenz 19/jan/93
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
      parameter (nexec=4,zero=0.d0,one=1.d0)
      dimension ranv(nran)
c
c ---------- initialization of x-table and warm up_ ----------
c
      if(iseed.ne.0) then
        call xuinit(x,np,nq,0,nexec,iseed,init,last,iout)
      end if
c
c -------- generation of vector of uniform deviates ----------
c
      if(nran.gt.0) then
        if(init.ne.0) then
c
c      expand x-vector into ranv vector, the first loop is
c      not recurrent, recurrence in second loop broken by
c      explicit stripmining with vector lengths of max. nq
c
          do i=1,min(nran,nq)
            xx=x(i)-x(i+np-nq)
            if(xx.lt.zero) xx=xx+one
            ranv(i)=xx
          end do
c
          if(nran.gt.nq) then
            istart=nq+1
  150       continue
cdir$ ivdep
            do i=istart,min(np,istart+nq-1)
              xx=x(i)-ranv(i-nq)
              if(xx.lt.zero) xx=xx+one
              ranv(i)=xx
            end do
            istart=istart+nq
            if(istart.le.np) goto 150
          end if
c
c      main computation within ranv vector; recurrence broken 
c      by explicit stripmining with vector lengths of max. nq
c
          if(nran.gt.np) then
            istart=np+1
  250       continue
cdir$ ivdep
            do i=istart,min(nran,istart+nq-1)
              xx=ranv(i-np)-ranv(i-nq)
              if(xx.lt.zero) xx=xx+one
              ranv(i)=xx
            end do
            istart=istart+nq
            if(istart.le.nran) goto 250
c
c      shift the np most recent prn's back into x vector
c
            do i=1,np
              x(i)=ranv(nran-np+i)
            end do
          else
c   this may look like a recurrence to some compilers
cdir$ ivdep
            do i=1,np-nran
              x(i)=x(i+nran)
            end do
            do i=np-nran+1,np
              x(i)=ranv(i+nran-np)
            end do
          end if
        else
          call errprt(iout,'vranfx','incorrect initialization',-1)
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine fslfg(x,np,nq,last,ranv,nran,iout)
c
c      generic floating point subtractive lagged fibonacci generator.
c      performs internal table update and copies random numbers into
c      output vector. seed table must have been initialized separately
c      e.g. by x u i n i t.
c ----------------------------------------------------------------------
c      alternative additive formulation bypassing if statements:
c      temp=x(k)+x(k+np-nq)
c      x(k)=temp-float(int(temp))
c
c      alternative subtractive formulation bypassing if statements:
c      temp=x(k)-x(k+np-nq)
c      x(k)=(temp+one)-float(int(temp+one))
c ----------------------------------------------------------------------
c      subroutines called: errprt            m. lewerenz dec/96, mar/98
c
      implicit real*8 (a-h,o-z)
      parameter (mroll=4,zero=0.d0,one=1.d0)
      dimension ranv(nran),x(np)
      common /doctrl/ nroll
c
      if(np.le.0.or.nq.le.0.or.nq.ge.np.or.last.ge.np.or.last.lt.0) then
        call errprt(iout,'fslfg','invalid table parameters',-1)
      else
        if(nran.gt.0) then
          nratio=np/nq
          j=0
          left=nran
c
c ------------------- straight loop version first ----------------------
c --------------------- best for vector machines -----------------------
c
          if(nroll.le.1) then
   10       continue
cdir$ ivdep
            do i=last+1,min(nq,left+last)
              x1=x(i)-x(i+np-nq)
              if(x1.lt.zero) x1=x1+one
              x(i)=x1
              j=j+1
              ranv(j)=x1
            end do
c
            if(last.lt.nratio*nq) then
              do k=1,nratio-1
cdir$ ivdep
                do i=max(k*nq,last)+1,min((k+1)*nq,left+last)
                  x1=x(i)-x(i-nq)
                  if(x1.lt.zero) x1=x1+one
                  x(i)=x1
                  j=j+1
                  ranv(j)=x1
                end do
              end do
            end if
            limit=min(np,left+last)
cdir$ ivdep
            do i=max(nratio*nq,last)+1,limit
              x1=x(i)-x(i-nq)
              if(x1.lt.zero) x1=x1+one
              x(i)=x1
              j=j+1
              ranv(j)=x1
            end do
c
            last=mod(limit,np)
            left=nran-j
            if(left.gt.0) goto 10
c
c --------------------- unrolled version of loops ----------------------
c ---------------------- best for risc machines ------------------------
c
          else
   20       continue
            loop=mod((min(nq,left+last)-last),mroll)
            do i=last+1,last+loop
              x1=x(i)-x(i+np-nq)
              if(x1.lt.zero) x1=x1+one
              x(i)=x1
              j=j+1
              ranv(j)=x1
            end do
            do i=last+loop+1,min(nq,left+last),mroll
              x1=x(i)-x(i+np-nq)
              x2=x(i+1)-x(i+1+np-nq)
              x3=x(i+2)-x(i+2+np-nq)
              x4=x(i+3)-x(i+3+np-nq)
              if(x1.lt.zero) x1=x1+one
              if(x2.lt.zero) x2=x2+one
              if(x3.lt.zero) x3=x3+one
              if(x4.lt.zero) x4=x4+one
              x(i)=x1
              x(i+1)=x2
              x(i+2)=x3
              x(i+3)=x4
              ranv(j+1)=x1
              ranv(j+2)=x2
              ranv(j+3)=x3
              ranv(j+4)=x4
              j=j+4
            end do
c
            if(last.lt.nratio*nq) then
              do k=1,nratio-1
                limit=min((k+1)*nq,left+last)
                loop=mod((limit-max(k*nq,last)),mroll)
                do i=max(k*nq,last)+1,max(k*nq,last)+loop
                  x1=x(i)-x(i-nq)
                  if(x1.lt.zero) x1=x1+one
                  x(i)=x1
                  j=j+1
                  ranv(j)=x1
                end do
                do i=max(k*nq,last)+loop+1,limit,mroll
                  x1=x(i)-x(i-nq)
                  x2=x(i+1)-x(i+1-nq)
                  x3=x(i+2)-x(i+2-nq)
                  x4=x(i+3)-x(i+3-nq)
                  if(x1.lt.zero) x1=x1+one
                  if(x2.lt.zero) x2=x2+one
                  if(x3.lt.zero) x3=x3+one
                  if(x4.lt.zero) x4=x4+one
                  x(i)=x1
                  x(i+1)=x2
                  x(i+2)=x3
                  x(i+3)=x4
                  ranv(j+1)=x1
                  ranv(j+2)=x2
                  ranv(j+3)=x3
                  ranv(j+4)=x4
                  j=j+4
                end do
              end do
            end if
c
            limit=min(np,left+last)
            loop=mod((limit-max(nratio*nq,last)),mroll)
            do i=max(nratio*nq,last)+1,max(nratio*nq,last)+loop
              x1=x(i)-x(i-nq)
              if(x1.lt.zero) x1=x1+one
              x(i)=x1
              j=j+1
              ranv(j)=x1
            end do
            do i=max(nratio*nq,last)+loop+1,limit,mroll
              x1=x(i)-x(i-nq)
              x2=x(i+1)-x(i+1-nq)
              x3=x(i+2)-x(i+2-nq)
              x4=x(i+3)-x(i+3-nq)
              if(x1.lt.zero) x1=x1+one
              if(x2.lt.zero) x2=x2+one
              if(x3.lt.zero) x3=x3+one
              if(x4.lt.zero) x4=x4+one
              x(i)=x1
              x(i+1)=x2
              x(i+2)=x3
              x(i+3)=x4
              ranv(j+1)=x1
              ranv(j+2)=x2
              ranv(j+3)=x3
              ranv(j+4)=x4
              j=j+4
            end do
c
            last=mod(limit,np)
            left=nran-j
            if(left.gt.0) goto 20
          end if
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine fslfg2(x,np,nq,last,a,b,ranv,nran,iout)
c
c      generic floating point subtractive lagged fibonacci generator on
c      the interval [a,b]. works on the same internal table as f s l f g
c      and allows mixed usage.
c      performs internal table update and copies random numbers into
c      output vector. seed table must have been initialized separately
c      e.g. by x u i n i t.
c      subroutines called: errprt            m. lewerenz dec/96, mar/98
c
      implicit real*8 (a-h,o-z)
      parameter (mroll=4,zero=0.d0,one=1.d0)
      dimension ranv(nran),x(np)
      common /doctrl/ nroll
c
      if(np.le.0.or.nq.le.0.or.nq.ge.np.or.last.ge.np.or.last.lt.0) then
        call errprt(iout,'fslfg2','invalid table parameters',-1)
      else
        if(nran.gt.0) then
          nratio=np/nq
          j=0
          left=nran
          range=b-a
c
c ------------------- straight loop version first ----------------------
c --------------------- best for vector machines -----------------------
c
          if(nroll.le.1) then
   10       continue
cdir$ ivdep
            do i=last+1,min(nq,left+last)
              x1=x(i)-x(i+np-nq)
              if(x1.lt.zero) x1=x1+one
              x(i)=x1
              j=j+1
              ranv(j)=range*x1+a
            end do
c
            if(last.lt.nratio*nq) then
              do k=1,nratio-1
cdir$ ivdep
                do i=max(k*nq,last)+1,min((k+1)*nq,left+last)
                  x1=x(i)-x(i-nq)
                  if(x1.lt.zero) x1=x1+one
                  x(i)=x1
                  j=j+1
                  ranv(j)=range*x1+a
                end do
              end do
            end if
            limit=min(np,left+last)
cdir$ ivdep
            do i=max(nratio*nq,last)+1,limit
              x1=x(i)-x(i-nq)
              if(x1.lt.zero) x1=x1+one
              x(i)=x1
              j=j+1
              ranv(j)=range*x1+a
            end do
c
            last=mod(limit,np)
            left=nran-j
            if(left.gt.0) goto 10
c
c --------------------- unrolled version of loops ----------------------
c ---------------------- best for risc machines ------------------------
c
          else
   20       continue
            loop=mod((min(nq,left+last)-last),mroll)
            do i=last+1,last+loop
              x1=x(i)-x(i+np-nq)
              if(x1.lt.zero) x1=x1+one
              x(i)=x1
              j=j+1
              ranv(j)=range*x1+a
            end do
            do i=last+loop+1,min(nq,left+last),mroll
              x1=x(i)-x(i+np-nq)
              x2=x(i+1)-x(i+1+np-nq)
              x3=x(i+2)-x(i+2+np-nq)
              x4=x(i+3)-x(i+3+np-nq)
              if(x1.lt.zero) x1=x1+one
              if(x2.lt.zero) x2=x2+one
              if(x3.lt.zero) x3=x3+one
              if(x4.lt.zero) x4=x4+one
              x(i)=x1
              x(i+1)=x2
              x(i+2)=x3
              x(i+3)=x4
              ranv(j+1)=range*x1+a
              ranv(j+2)=range*x2+a
              ranv(j+3)=range*x3+a
              ranv(j+4)=range*x4+a
              j=j+4
            end do
c
            if(last.lt.nratio*nq) then
              do k=1,nratio-1
                limit=min((k+1)*nq,left+last)
                loop=mod((limit-max(k*nq,last)),mroll)
                do i=max(k*nq,last)+1,max(k*nq,last)+loop
                  x1=x(i)-x(i-nq)
                  if(x1.lt.zero) x1=x1+one
                  x(i)=x1
                  j=j+1
                  ranv(j)=range*x1+a
                end do
                do i=max(k*nq,last)+loop+1,limit,mroll
                  x1=x(i)-x(i-nq)
                  x2=x(i+1)-x(i+1-nq)
                  x3=x(i+2)-x(i+2-nq)
                  x4=x(i+3)-x(i+3-nq)
                  if(x1.lt.zero) x1=x1+one
                  if(x2.lt.zero) x2=x2+one
                  if(x3.lt.zero) x3=x3+one
                  if(x4.lt.zero) x4=x4+one
                  x(i)=x1
                  x(i+1)=x2
                  x(i+2)=x3
                  x(i+3)=x4
                  ranv(j+1)=range*x1+a
                  ranv(j+2)=range*x2+a
                  ranv(j+3)=range*x3+a
                  ranv(j+4)=range*x4+a
                  j=j+4
                end do
              end do
            end if
c
            limit=min(np,left+last)
            loop=mod((limit-max(nratio*nq,last)),mroll)
            do i=max(nratio*nq,last)+1,max(nratio*nq,last)+loop
              x1=x(i)-x(i-nq)
              if(x1.lt.zero) x1=x1+one
              x(i)=x1
              j=j+1
              ranv(j)=range*x1+a
            end do
            do i=max(nratio*nq,last)+loop+1,limit,mroll
              x1=x(i)-x(i-nq)
              x2=x(i+1)-x(i+1-nq)
              x3=x(i+2)-x(i+2-nq)
              x4=x(i+3)-x(i+3-nq)
              if(x1.lt.zero) x1=x1+one
              if(x2.lt.zero) x2=x2+one
              if(x3.lt.zero) x3=x3+one
              if(x4.lt.zero) x4=x4+one
              x(i)=x1
              x(i+1)=x2
              x(i+2)=x3
              x(i+3)=x4
              ranv(j+1)=range*x1+a
              ranv(j+2)=range*x2+a
              ranv(j+3)=range*x3+a
              ranv(j+4)=range*x4+a
              j=j+4
            end do
c
            last=mod(limit,np)
            left=nran-j
            if(left.gt.0) goto 20
          end if
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      function ranf(iseed,iout)
c
c      generator for uniformly distributed pseudo random numbers using
c      the fibonacci method with 48 bit mantissa output. single output
c      version of vranf. see comments there for method. successive calls
c      generate the same random sequence as vranf. ranf and vranf can be
c      used together, working on the same random sequence, and either
c      one can be used for initialization. this generator is machine
c      independent and gives the same random sequence on any machine.
c      5.52 s for 1000000 calls on convex c210. 3.06 s on ibm-3090/300s
c
c      iseed : if not 0, generator is initialized and ranf returns a
c              real echo of iseed; no random number output. iseed is
c              unchanged.
c              if 0, ranf returns the next random number from a
c              previously used or initialized random sequence.
c      iout  : unit number for error messages. silent for iout.le.0
c      subroutines called: vranf                m. lewerenz 12/may/91
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
      parameter (zero=0.d0,one=1.d0)
c
c      table initialization by vranf
c
      if(iseed.ne.0) then
        call vranf(dummy,0,iseed,iout)
        ranf=iseed
      else
        if(init.ne.0) then
          if(last.lt.nq) then
            ranf=x(last+1)-x(last+1+np-nq)
          else
            ranf=x(last+1)-x(last+1-nq)
          end if
          if(ranf.lt.zero) ranf=ranf+one
          x(last+1)=ranf
          last=mod(last+1,np)
        else
          call errprt(iout,'ranf','incorrect initialization',-1)
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine rsavef(isave,iout)
c
c      table backup for routines ranf, vranf, and vranf2.
c      uses unit iabs(isave) to save or retrieve contents of common
c      block /xrandf/ for program restart with continuation of the
c      old random sequence. isave > 0 write, isave < 0 read.
c      iout is a message unit. np and nq must be consistent with
c      vranf and ranf!
c      subroutines called: errprt,xuinit      m. lewerenz 12/may/91
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
c
      iunit=iabs(isave)
      if(isave.gt.0) then
        write(iunit) np,nq,last
        write(iunit) x
      else if(isave.lt.0) then
        read(iunit) npp,nqq,last
        if(npp.ne.np.or.nqq.ne.nq.or.nqq.ge.npp.or.last.lt.0.
     #               or.last.gt.np) then
          call errprt(iout,'rsavef',
     #                'illegal file contents -> using default seed',0)
          iseed=123456789
          call xuinit(x,np,nq,0,4,iseed,init,last,iout)
        else
          read(iunit) x
          init=1
        end if
      else
        call errprt(iout,'rsavef','illegal unit number',-1)
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine rnsetf(xsave,nx,iflag,iout)
c
c      sets or gets state of generator table used in routines ranf,
c      vranf, and vranf2.
c
c      xsave : vector of length nx containing the state of the
c              generator. nx must be > np; input or output
c      nx    : length of vector xsave, must be at least np+2 which
c              is currently 252; input
c      iflag : 0 -> xsave returns the state of the generator;
c              else -> generator state is set to xsave; input
c      iout  : unit number for messages; silent for iout.le.0; input
c      np and nq must be consistent with settings in vranf and ranf!
c      subroutines called: errprt                    m. lewerenz jan/94
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
      dimension xsave(nx)
c
      if(nx.ge.(np+2)) then
        if(iflag.eq.0) then
          call vcopy(x,xsave,np)
          xsave(np+1)=last
          xsave(np+2)=init
        else
          call vcopy(xsave,x,np)
          last=xsave(np+1)
          init=xsave(np+2)
          if(init.eq.0) then
            call errprt(iout,'rnsetf','bad generator state',-1)
          end if
        end if
      else
        call errprt(iout,'rnsetf','xsave vector too short',-1)
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vrseed(iseed,iout)
c
c      returns current status of internal congruential generator used
c      to seed the bit tables for ranf/vranf package. useful to seed
c      other copies of ranf/vranf.
c
c      iseed : seed from current congruential generator status; output
c      iout  : unit number for messages; silent for iout.le.0; input
c      subroutines called: errprt                    m. lewerenz feb/98
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
c
      iseed=init
      if(init.eq.0) call errprt(iout,'vrseed','illegal seed',1)
      return
      end
c=======================================================================
c================== last line of ranf/vranf package ====================
c=======================================================================
c-----------------------------------------------------------------------
c-------------------------- error handling -----------------------------
c-----------------------------------------------------------------------
c
      subroutine errprt(iout,pgname,text,icode)
c
c      prints error messages from library subroutines
c
c      iout   : unit number for message output, 0-> no output; input
c      pgname : name of the subroutine calling errprt; input
c      text   : message text; input
c      icode  : severity code: 0 -> warning, < 0 -> fatal error with
c               abort, else -> error but execution continues
c      subroutines called : strlen                   m. lewerenz dec/93
c
      character pgname*(*),text*(*),header*20,tail*40
      save nerror,nwarn,icall
      common /errcnt/ maxerr,maxwrn
      data icall/0/
c
      if(icall.eq.0) then
        icall=1
        nerror=0
        nwarn=0
      end if
c
      if(icode.lt.0) then
        header='  *** fatal error,'
        tail=', execution aborted ***'
      else if(icode.eq.0) then
        header='  *** warning,'
        tail=' ***'
        nwarn=nwarn+1
      else
        header='  *** error,'
        tail=', return without action ***'
        nerror=nerror+1
      end if
c
c      write the message on unit iout
c
      if(iout.gt.0) then
        call strlen(pgname,lname,iout)
        call strlen(text,ltext,iout)
        call strlen(header,lhead,iout)
        call strlen(tail,ltail,iout)
        write(iout,'(/6a/)') header(1:lhead),' ',text(1:ltext),' in ',
     #                      pgname(1:lname),tail(1:ltail)
        call flush(iout)
      end if
c
      jcode=icode
      if(maxerr.gt.0.and.nerror.ge.maxerr) then
        if(iout.gt.0) write(iout,'(/a)')
     #  '  *** maximum number of errors exceeded, program stopped *** '
        jcode=-1
      end if
      if(maxwrn.gt.0.and.nwarn.ge.maxwrn) then
        if(iout.gt.0) write(iout,'(/a)')
     #  '  *** maximum number of warnings exceeded, program stopped ***'
        jcode=-1
      end if
      if(iout.gt.0) call flush(iout)
c
      if(jcode.lt.0) stop
      return
c
c ---------------------------------------------------------------------
c      error report, returns current number of errors and warning 
c
      entry errnum(nerr,nwrn)
      nerr=nerror
      nwrn=nwarn
      return
      end
C
C----------------------------------------------------------------------
C
      SUBROUTINE STRLEN(STRING,LS,IOUT)
C
C      DETERMINES LENGTH OF STRING
C
      CHARACTER STRING*(*)
C
      LS=LEN(STRING)
   10 IF(LS.GT.0.AND.STRING(LS:LS).EQ.' ') THEN
        LS=LS-1
        GOTO 10
      END IF
      IF(LS.EQ.0) call errprt(iout,'strlen','empty string',0)
      RETURN
      END
c
c-----------------------------------------------------------------------
c
      subroutine flush(iunit)
c
c      dummy to resolve calls to the buffer flushing routine
c      available on dec axp and sgi machines
c
      return
      end
