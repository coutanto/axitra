c @(#) parameter.F	AXITRA 4.15     9/22/97 4
c **********************************************************************
c *                                                                    *
c *                  PARAMETER                                         *
c *                                                                    *
c *    nsp: nbre max source
c *    nrp: nbre max stat
c *    ncp: nbre max couche
c *    ikmin: nbre min d iterations
c *    nkmax: (si fct de Bessel stockees) dimension tab fct Bessel
c **********************************************************************

      implicit real*8 (a-b,d-h,o-z)
      implicit complex*16 (c)

      integer*2   in1,in2,in3,out,out2
      parameter (in1=10,in2=11,in3=12,out=13,out2=14)
c
c nsp= max number of source
c nrp= max number of receiver
c ncp= max number of layer
      integer   nsp,nrp,ncp,nrsp,ikmin
      parameter (nsp=5,nrp=4,ncp=2,nrsp=nrp*nsp,ikmin=100)

c explim= exponent lower limit
c needs to be adjuest when running with single precision float number
c in single precision, -30 is a good value.
c don't forget to change cdexp to cexp
      real*8	explim,elim
c     parameter (explim=-600.,elim=1.d-300)
      parameter (explim=-300.,elim=1.d-300)

c bessel functions are stored up to nkmax wavenumber iterations.
c the larger, the faster, but needs more memory
      integer   nkmax
      parameter (nkmax=2000)

c convergence relative error
      real*8   rerr
      parameter (rerr=1.e-5)

c reference frequency for attenuation
      real*8   fref
      parameter (fref=1.)
