      program ed
c     Constants
      integer   NMAX, PMAX, DMAX
      parameter (NMAX = 32, PMAX = 128, DMAX = 6)

c     Local variables
      integer        np, nd, ngen, imut, irep, ielite, ivrb, k, ip, ig,
     +               ip1, ip2, new, newtot
c
      real           ph(NMAX,2), oldph(NMAX,PMAX), newph(NMAX,PMAX)
c
      integer        gn1(NMAX*DMAX), gn2(NMAX*DMAX)

      np = 10
      n  = 1
      nd = 6

      call rninit(12345)

c     Compute initial (random but bounded) phenotypes
      do 1 ip=1,np
         do 2 k=1,n
            oldph(k,ip)=urand()
    2    continue
    1 continue

      ip1 = 1
      
      call encode(n,nd,oldph(1,ip1),gn1)
      write(*,*) 'phen= ', oldph(1,ip1),'gene= ', gn1(1:nd)
c     
      call decode(n,nd,gn1,ph(1,1))

c http://www.hao.ucar.edu/modeling/pikaia/breeding.php

      n  = 2
      nd = 8
      ip1 = 1

      oldph(1,1) = 0.123456789
      oldph(2,1) = 0.987654321
      call encode(n,nd,oldph(1,ip1),gn1)
      write(*,*) 'phen= ', oldph(1,ip1),'gene= ', gn1(1:nd*n)



      end
c*********************************************************************
      function urand()
c=====================================================================
c     Return the next pseudo-random deviate from a sequence which is
c     uniformly distributed in the interval [0,1]
c
c     Uses the function ran0, the "minimal standard" random number
c     generator of Park and Miller (Comm. ACM 31, 1192-1201, Oct 1988;
c     Comm. ACM 36 No. 7, 105-110, July 1993).
c=====================================================================
      implicit none
c
c     Input - none
c
c     Output
      real     urand
c
c     Local
      integer  iseed
      real     ran0
      external ran0
c
c     Common block to make iseed visible to rninit (and to save
c     it between calls)
      common /rnseed/ iseed
c
      urand = ran0( iseed )
      return
      end
c*********************************************************************
      subroutine rninit( seed )
c=====================================================================
c     Initialize random number generator urand with given seed
c=====================================================================
      implicit none
c
c     Input
      integer seed
c
c     Output - none
c
c     Local
      integer iseed
c
c     Common block to communicate with urand
      common /rnseed/ iseed
c
c     Set the seed value
      iseed = seed
      if(iseed.le.0) iseed=123456
      return
      end
c*********************************************************************
      function ran0( seed )
c=====================================================================
c     "Minimal standard" pseudo-random number generator of Park and
c     Miller.  Returns a uniform random deviate r s.t. 0 < r < 1.0.
c     Set seed to any non-zero integer value to initialize a sequence,
c     then do not change seed between calls for successive deviates
c     in the sequence.
c
c     References:
c        Park, S. and Miller, K., "Random Number Generators: Good Ones
c           are Hard to Find", Comm. ACM 31, 1192-1201 (Oct. 1988)
c        Park, S. and Miller, K., in "Remarks on Choosing and Imple-
c           menting Random Number Generators", Comm. ACM 36 No. 7,
c           105-110 (July 1993)
c=====================================================================
c *** Declaration section ***
c
      implicit none
c
c     Input/Output:
      integer seed
c
c     Output:
      real ran0
c
c     Constants:
      integer A,M,Q,R
      parameter (A=48271,M=2147483647,Q=44488,R=3399)
      real SCALE,EPS,RNMX
      parameter (SCALE=1./M,EPS=1.2e-7,RNMX=1.-EPS)
c
c     Local:
      integer j
c
c *** Executable section ***
c
      j = seed/Q
      seed = A*(seed-j*Q)-R*j
      if (seed .lt. 0) seed = seed+M
      ran0 = min(seed*SCALE,RNMX)
      return
      end
c**********************************************************************
      subroutine encode(n,nd,ph,gn)
c======================================================================
c     encode phenotype parameters into integer genotype
c     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
c======================================================================
c
      implicit none
c
c     Inputs:
      integer   n, nd
      real      ph(n)
c
c     Output:
      integer   gn(n*nd)
c
c     Local:
      integer   ip, i, j, ii
      real      z
c
      write(*,*) 'ph= ', ph

      z=10.**nd
      write(*,*) 'z= ', z
      ii=0
      do 1 i=1,n
         ip=int(ph(i)*z)
         write(*,*) 'ip= ', ip
         do 2 j=nd,1,-1
            gn(ii+j)=mod(ip,10)
            ip=ip/10
            write(*,*) ii, j, mod(ip,10), ip
    2   continue
        ii=ii+nd
    1 continue
 
      return
      end
 
c**********************************************************************
      subroutine decode(n,nd,gn,ph)
c======================================================================
c     decode genotype into phenotype parameters
c     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
c======================================================================
c
      implicit none
c
c     Inputs:
      integer   n, nd, gn(n*nd)
c
c     Output:
      real      ph(n)
c
c     Local:
      integer   ip, i, j, ii
      real      z
c
      z=10.**(-nd)
      ii=0
      do 1 i=1,n
         ip=0
         do 2 j=1,nd
            ip=10*ip+gn(ii+j)
    2    continue
         ph(i)=ip*z
         ii=ii+nd
    1 continue
 
      return
      end

