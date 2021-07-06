      parameter (lamax=30,nchan=3500,nchane2e=10,nchanop=1)
      parameter (nnmax=401,lnabmax=lamax,ltmax=30,ncmax=nnmax)
      parameter (lmax=400+ltmax,kmax=401,lcoul=lmax)
C  Note that the MAXR parameter must be set the same as NSIZEK7 in paratom.f
      PARAMETER (maxr=20000,larged=kmax*nchanop*2*500,largee=larged)
      parameter (maxrc = maxr)
      parameter(ispline=100)
c$$$  PARAMETER (maxr=3000,larged=kmax*nchan*lamax*1100,largee=larged)
C  The following parameters are necessary in Dima's helium code
C  NCMAX is also used, but is the number of configurations for each
C  symmetry of a helium target state. This is the max size of the of
C  the diagonalization matrix.
      parameter (lomax = lamax)
C  LOMAX is the max l in any single particle function.
      parameter (komax = nnmax, maxfac=100)
C  KOMAX is the max number of single particle functions for a particular l.
C  MAXFAC is the max number of factorials calculated.
      parameter(KNM = nchan/2)
C  KNM is the max number of helium states in the CC calculation
      parameter (nspmax = 3000)
C  NSPMAX is the max number of single particle functions <= lomax * komax.
      parameter (nmaxr = maxr)
      parameter (nsetm = 2)
C  NSETM allows for different exponential fall offs within a Laguerre basis
      parameter (nqmax = kmax)
      parameter (nchmax = nchan)
      parameter (ncmCI  = 1200, nspmCI = 1200, nicmax = 50)

* Alisher's addenda
     
