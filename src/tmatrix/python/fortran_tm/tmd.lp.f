C   This is NOT the original code by the authors listed below!
C   Modified 2009-2013 by Jussi Leinonen (jsleinonen@gmail.com)
C   to be compatible with the Python extension interface.
C   The requirement of non-profit use only has been dropped from
C   this code with the permission of M. I. Mishchenko; thus the
C   license has been made open source compatible.

C   New release including the LAPACK matrix inversion procedure.
C   We thank Cory Davis (University of Edinburgh) for pointing
C   out the possibility of replacing the proprietary NAG matrix 
C   inversion routine by the public-domain LAPACK equivalent.

C   CALCULATION OF LIGHT SCATTERING BY POLYDISPERSE, RANDOMLY          
C   ORIENTED PARTICLES OF IDENTICAL AXIALLY SYMMETRIC SHAPE      

C   This version of the code uses DOUBLE PRECISION variables
C   and must be used along with the accompanying files tmd.par.f      
C   and lpd.f.    
                                                                       
C   Last update 08/06/2005 
                                                                       
C   The code has been developed by Michael Mishchenko at the NASA      
C   Goddard Institute for Space Studies, New York. This research
C   was funded by the NASA Radiation Sciences Program.
                                                                       
C   The code can be used without limitations in any not-for-       
C   profit scientific research.  We only request that in any      
C   publication using the code the source of the code be acknowledged  
C   and relevant references (see below) be made.                                   
                                                                       
C   This version of the code is applicable to spheroids,               
C   Chebyshev particles, and finite circular cylinders.                
                                                                       
C   The computational method is based on the Watermsn's T-matrix      
C   approach and is described in detail in the following papers:           
C                                                                      
C   1.  M. I. Mishchenko, Light scattering by randomly oriented        
C       axially symmetric particles, J. Opt. Soc. Am. A,               
C       vol. 8, 871-882 (1991).                                        
C                                                                      
C   2.  M. I. Mishchenko, Light scattering by size-shape               
C       distributions of randomly oriented axially symmetric           
C       particles of a size comparable to a wavelength,                
C       Appl. Opt., vol. 32, 4652-4666 (1993).                         
C                                                                      
C   3.  M. I. Mishchenko and L. D. Travis, T-matrix computations       
C       of light scattering by large spheroidal particles,             
C       Opt. Commun., vol. 109, 16-21 (1994).                          
C                                                                      
C   4.  M. I. Mishchenko, L. D. Travis, and A. Macke, Scattering       
C       of light by polydisperse, randomly oriented, finite            
C       circular cylinders, Appl. Opt., vol. 35, 4927-4940 (1996).     
C                                                                      
C   5.  D. J. Wielaard, M. I. Mishchenko, A. Macke, and B. E. Carlson, 
C       Improved T-matrix computations for large, nonabsorbing and    
C       weakly absorbing nonspherical particles and comparison         
C       with geometrical optics approximation, Appl. Opt., vol. 36,    
C       4305-4313 (1997).                                                  
C                                                                      
C   A general review of the T-matrix approach can be found in          
C                                                                      
C   6.  M. I. Mishchenko, L. D. Travis, and D. W. Mackowski,           
C       T-matrix computations of light scattering by nonspherical      
C       particles: a review, J. Quant. Spectrosc. Radiat.             
C       Transfer, vol. 55, 535-575 (1996).                             
C                                                                      
C   The following paper provides a detailed user guide to the          
C   T-matrix code:                                                     
C                                                                      
C   7.  M. I. Mishchenko and L. D. Travis, Capabilities and            
C       limitations of a current FORTRAN implementation of the         
C       T-matrix method for randomly oriented, rotationally            
C       symmetric scatterers, J. Quant. Spectrosc. Radiat. Transfer,   
C       vol. 60, 309-324 (1998).                                       
C
C   These papers are available in the .pdf format at the web site 
C
C   http://www.giss.nasa.gov/~crmim/publications/
C
C   or in hardcopy upon request from Michael Mishchenko    
C   Please e-mail your request to crmim@giss.nasa.gov.      
C
C   A comprehensive book "Scattering, Absorption, and Emission of  
C   Light by Small Particles" (Cambridge University Press, Cambridge,
C   2002) is also available in the .pdf format at the web site 
C
C   http://www.giss.nasa.gov/~crmim/books.html
 
C   Analytical averaging over particle orientations (Ref. 1) makes     
C   this method the fastest exact technique currently available.       
C   The use of an automatic convergence procedure                      
C   (Ref. 2) makes the code convenient in massive computations.        
C   Ref. 4 describes features specific for finite cylinders as         
C   particles with sharp rectangular edges.  Ref. 5 describes further      
C   numerical improvements.                                               
                                                                       
C   The use of extended precision variables (Ref. 3) can                     
C   significantly increase the maximal convergent equivalent-sphere
C   size parameter and make it greater than 200 (depending on
C   refractive index and aspect ratio). The extended-precision code
C   is also available. However, the use of extended precision varibales     
C   results in a greater consumption of CPU time.                  
C   On IBM RISC workstations, that code is approximately               
C   five times slower than this double-precision code.  The            
C   CPU time difference between the double-precision and extended-     
C   precision codes can be larger on supercomputers.                   
                                                                       
C   This is the first part of the full T-matrix code.  The second part,    
C   lpd.f, is completely independent of the firsti part. It contains no    
C   T-matrix-specific subroutines and can be compiled separately.
C   The second part of the code replaces the previously implemented 
C   standard matrix inversion scheme based on Gaussian elimination 
C   by a scheme based on the LU factorization technique.     
C   As described in Ref. 5 above, the use of the LU factorization is          
C   especially beneficial for nonabsorbing or weakly absorbing particles.     
C   In this code we use the LAPACK implementation of the LU factorization
C   scheme. LAPACK stands for Linear Algebra PACKage. The latter is 
C   publicly available at the following internet site:
C
C   http://www.netlib.org/lapack/      
                                                                      
C   INPUT PARAMETERS:                                                  
C                                                                      
C      RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C      RAT.NE.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius                      
C      NDISTR specifies the distribution of equivalent-sphere radii   
C      NDISTR = 1 - modified gamma distribution                        
C           [Eq. (40) of Ref. 7]                                       
C               AXI=alpha                                              
C               B=r_c                                                  
C               GAM=gamma                                              
C      NDISTR = 2 - log-normal distribution                            
C           [Eq. 41) of Ref. 7]                                        
C               AXI=r_g                                                
C               B=[ln(sigma_g)]**2                                    
C      NDISTR = 3 - power law distribution                             
C           [Eq. (42) of Ref. 7]                                       
C                AXI=r_eff (effective radius)                   
C                B=v_eff (effective variance)                
C                Parameters R1 and R2 (see below) are calculated       
C                automatically for given AXI and B
C      NDISTR = 4 - gamma distribution                                 
C           [Eq. (39) of Ref. 7]                                       
C                AXI=a                                                 
C                B=b                                                   
C      NDISTR = 5 - modified power law distribution
C         [Eq. (24) in M. I. Mishchenko et al.,
C         Bidirectional reflectance of flat,
C         optically thick particulate laters: an efficient radiative
C         transfer solution and applications to snow and soil surfaces,
C         J. Quant. Spectrosc. Radiat. Transfer, Vol. 63, 409-432 (1999)].
C                B=alpha
C                                                                      
C      The code computes NPNAX size distributions of the same type     
C      and with the same values of B and GAM in one run.               
C      The parameter AXI varies from AXMAX to AXMAX/NPNAX in steps of  
C      AXMAX/NPNAX.  To compute a single size distribution, use        
C      NPNAX=1 and AXMAX equal to AXI of this size distribution.       
C                                                                      
C      R1 and R2 - minimum and maximum equivalent-sphere radii         
C           in the size distribution. They are calculated automatically              
C           for the power law distribution with given AXI and B   
C           but must be specified for other distributions        
C           after the lines           
C                                                                      
C             DO 600 IAX=1,NPNAX                            
C                AXI=AXMAX-DAX*DFLOAT(IAX-1)                
C                                                                      
C           in the main program.                                 
C           For the modified power law distribution (NDISTR=5), the
C           minimum radius is 0, R2 is the maximum radius,
C           and R1 is the intermediate radius at which the
C           n(r)=const dependence is replaced by the power law
C           dependence.
C                                                                      
C      NKMAX.LE.988 is such that NKMAX+2 is the                        
C           number of Gaussian quadrature points used in               
C           integrating over the size distribution for particles with
C           AXI=AXMAX.  For particles with AXI=AXMAX-AXMAX/NPNAX,      
C           AXMAX-2*AXMAX/NPNAX, etc. the number of Gaussian points    
C           linearly decreases.                                       
C           For the modified power law distribution, the number
C           of integration points on the interval [0,R1] is also
C           equal to NKMAX.
C                                                                      
C      LAM - wavelength of light                                       
C      MRR and MRI - real and imaginary parts of the refractive        
C                  index (MRI.GE.0)   
C      EPS and NP - specify the shape of the particles.                
C             For spheroids NP=-1 and EPS is the ratio of the          
C                 horizontal to rotational axes.  EPS is larger than   
C                 1 for oblate spheroids and smaller than 1 for       
C                 prolate spheroids.                                   
C             For cylinders NP=-2 and EPS is the ratio of the          
C                 diameter to the length.                              
C             For Chebyshev particles NP must be positive and 
C                 is the degree of the Chebyshev polynomial, while     
C                 EPS is the deformation parameter                     
C                 [Eq. (33) of Ref. 7].      
C      DDELT - accuracy of the computations                            
C      NPNA - number of equidistant scattering angles (from 0      
C             to 180 deg) for which the scattering matrix is           
C             calculated.                                              
C      NDGS - parameter controlling the number of division points      
C             in computing integrals over the particle surface.        
C             For compact particles, the recommended value is 2.       
C             For highly aspherical particles larger values (3, 4,...) 
C             may be necessary to obtain convergence.                  
C             The code does not check convergence over this parameter. 
C             Therefore, control comparisons of results obtained with  
C             different NDGS-values are recommended.                   
                                                                       

C   OUTPUT PARAMETERS:                                                 
C                                                                      
C      REFF and VEFF - effective radius and effective variance of      
C          the size distribution as defined by Eqs. (43)-(45) of
C          Ref. 7.       
C      CEXT - extinction cross section per particle                    
C      CSCA - scattering cross section per particle                    
C      W - single scattering albedo                                    
C      <cos> - asymmetry parameter of the phase function               
C      ALPHA1,...,BETA2 - coefficients appearing in the expansions
C          of the elements of the scattering matrix in
C          generalized spherical functions
C          [Eqs. (11)-(16) of Ref. 7].
C      F11,...,F44 - elements of the normalized scattering matrix [as      
C          defined by Eqs. (1)-(3) of Ref. 7] versus scattering angle

C   Note that LAM, r_c, r_g, r_eff, a, R1, and R2 must 
C   be given in the same units of length, and that 
C   the dimension of CEXT and CSCA is that of LAM squared (e.g., if        
C   LAM and AXI are given in microns, then CEXT and CSCA are         
C   calculated in square microns).    
                                                                       
C   The physical correctness of the computed results is tested using   
C   the general inequalities derived by van der Mee and Hovenier,      
C   Astron. Astrophys., vol. 228, 559-568 (1990).  Although            
C   the message that the test of van der Mee and Hovenier is satisfied 
C   does not guarantee that the results are absolutely correct,        
C   the message that the test is not satisfied can mean that something 
C   is wrong.                                                          
                                                                       
C   The convergence of the T-matrix method for particles with          
C   different sizes, refractive indices, and aspect ratios can be      
C   dramatically different.  Usually, large sizes and large aspect     
C   ratios cause problems.  The user of this code                      
C   should first experiment with different input parameters in          
C   order to get an idea of the range of applicability of this         
C   technique.  Sometimes decreasing the aspect ratio                  
C   from 3 to 2 can increase the maximum convergent equivalent-        
C   sphere size parameter by a factor of several (Ref. 7).                      
C   The CPU time required rapidly increases with increasing ratio      
C   radius/wavelength and/or with increasing particle asphericity.     
C   This should be taken into account in planning massive computations.
C   Using an optimizing compiler on IBM RISC workstations saves        
C   about 70% of CPU time.                                             
                                                                       
C   Execution can be automatically terminated if dimensions of certain 
C   arrays are not big enough or if the convergence procedure decides  
C   that the accuracy of double precision variables is insufficient    
C   to obtain a converged T-matrix solution for given particles.       
C   In all cases, a message appears explaining the cause of termination. 
                                                                       
C   The message                                                        
C        "WARNING:  W IS GREATER THAN 1"                               
C   means that the single-scattering albedo exceeds the maximum        
C   possible value 1.  If W is greater than 1 by more than             
C   DDELT, this message can be an indication of numerical              
C   instability caused by extreme values of particle parameters.       
                                                                       
C   The message "WARNING: NGAUSS=NPNG1" means that convergence over    
C   the parameter NG (see Ref. 2) cannot be obtained for the NPNG1     
C   value specified in the PARAMETER statement in the file tmd.par.f. 
C   Often this is not a serious problem, especially for compact         
C   particles.                                                         
                                                                       
C   Larger and/or more aspherical particles may require larger         
C   values of the parameters NPN1, NPN4, and NPNG1 in the file
C   tmd.par.f.  It is recommended to keep NPN1=NPN4+25 and
C   NPNG1=3*NPN1.  Note that the memory requirement increases    
C   as the third power of NPN4. If the memory of                  
C   a computer is too small to accomodate the code in its current    
C   setting, the parameters NPN1, NPN4, and NPNG1 should be
C   decreased. However, this will decrease the maximum size parameter  
C   that can be handled by the code.                                   
                                                                       
C   In some cases any increases of NPN1 will not make the T-matrix     
C   computations convergent.  This means that the particle is just     
C   too "bad" (extreme size parameter and/or extreme aspect ratio      
C   and/or extreme refractive index; see Ref. 7).              
C   The main program contains several PRINT statements which are       
C   currently commentd out.  If uncommented, these statements will     
C   produce numbers which show the convergence rate and can be         
C   used to determine whether T-matrix computations for given particle 
C   parameters will converge at all.   
                                                                       
C   Some of the common blocks are used to save memory rather than      
C   to transfer data.  Therefore, if a compiler produces a warning     
C   message that the lengths of a common block are different in        
C   different subroutines, this is not a real problem.                 
                                                                       
C   The recommended value of DDELT is 0.001.  For bigger values,       
C   false convergence can be obtained.                                 
                                                                       
C   In computations for spheres use EPS=1.000001 instead of EPS=1.     
C   The use of EPS=1 can cause overflows in some rare cases.           
                                                                       
C   To calculate a monodisperse particle, use the options              
C        NPNAX=1                                                       
C        AXMAX=R                                                       
C        B=1D-1   
C        NKMAX=-1                                                      
C        NDISTR=4                                                      
C        ...                                                           
C        DO 600 IAX=1,NPNAX                                            
C           AXI=AXMAX-DAX*DFLOAT(IAX-1)                                
C           R1=0.9999999*AXI                                           
C           R2=1.0000001*AXI                                          
C        ...                                                           
C   where R is the equivalent-sphere radius.                           
                                                                       
C   It is recommended to use the power law rather than the             
C   gamma size distribution, because in this case convergent solution  
C   can be obtained for larger REFF and VEFF assuming the same
C   maximal R2 (Mishchenko and Travis, Appl. Opt., vol. 33, 7206-7225,
C   1994).           
                                                                       
C   For some compilers, DACOS must be raplaced by DARCOS and DASIN     
C   by DARSIN.                                                         
                                                                       
C   If many different size distributions are computed and the          
C   refractive index is fixed, then another approach can be more       
C   efficient than running this code many times.  Specifically,        
C   scattering results should be computed for monodisperse particles   
C   with sizes ranging from essentially zero to some maximum value     
C   with a small step size (say, 0.02 microns).  These results         
C   should be stored on disk and can be used along with spline         
C   interpolation to compute scattering by particles with intermediate 
C   sizes.  Scattering patterns for monodisperse nonspherical          
C   particles in random orientation are (much) smoother than for 
C   monodisperse spheres, and spline interpolation usually gives good 
C   results. In this way, averaging over any size distribution is a
C   matter of seconds.  For more on size averaging, see Refs. 2 and 4. 
                                                                       
C   We would highly appreciate informing me of any problems encountered 
C   with this code.  Please send your message to the following         
C   e-mail address:  CRMIM@GISS.NASA.GOV.                              

C   WHILE THE COMPUTER PROGRAM HAS BEEN TESTED FOR A VARIETY OF CASES,
C   IT IS NOT INCONCEIVABLE THAT IT CONTAINS UNDETECTED ERRORS. ALSO,
C   INPUT PARAMETERS CAN BE USED WHICH ARE OUTSIDE THE ENVELOPE OF
C   VALUES FOR WHICH RESULTS ARE COMPUTED ACCURATELY. FOR THIS REASON,
C   THE AUTHORS AND THEIR ORGANIZATION DISCLAIM ALL LIABILITY FOR
C   ANY DAMAGES THAT MAY RESULT FROM THE USE OF THE PROGRAM. 
                                                                       
      SUBROUTINE CALCRAND(RAT,LAM,MRR,MRI,EPS,NP,DDELT,NDGS,
     &     NPNAX,AXMAX,B,GAM,NDISTR,NKMAX,NPNA,NCOEFF,
     &     REFF,VEFF,CEXTIN,CSCAT,WALB,ASYMM, F)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'tmd.par.f'
      REAL*8 REFF,VEFF,CEXTIN,CSCAT,WALB,ASYMM
      REAL*8 ALPHA(NPL,4),BETA(NPL,2),F(NPNA,4,4)
Cf2py intent(in) rat
Cf2py intent(in) lam
Cf2py intent(in) mrr
Cf2py intent(in) mri
Cf2py intent(in) eps
Cf2py intent(in) np
Cf2py intent(in) ddelt
Cf2py intent(in) ndgs
Cf2py intent(in) npnax
Cf2py intent(in) axmax
Cf2py intent(in) b
Cf2py intent(in) gam
Cf2py intent(in) ndistr
Cf2py intent(in) nkmax
Cf2py intent(in) npna
Cf2py intent(in) ncoeff
Cf2py intent(out) reff
Cf2py intent(out) veff
Cf2py intent(out) cextin
Cf2py intent(out) cscat
Cf2py intent(out) walb
Cf2py intent(out) asymm
Cf2py intent(out) alpha
Cf2py intent(out) beta
Cf2py intent(out) f

      REAL*8  LAM,MRR,MRI,X(NPNG2),W(NPNG2),S(NPNG2),SS(NPNG2),
     *        AN(NPN1),R(NPNG2),DR(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      REAL*8  XG(1000),WG(1000),TR1(NPN2,NPN2),TI1(NPN2,NPN2),
     &        ALPH1(NPL),ALPH2(NPL),ALPH3(NPL),ALPH4(NPL),BET1(NPL),
     &        BET2(NPL),XG1(2000),WG1(2000),
     &        AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BE1(NPL),BE2(NPL)
      REAL*4
     &     RT11(NPN6,NPN4,NPN4),RT12(NPN6,NPN4,NPN4),
     &     RT21(NPN6,NPN4,NPN4),RT22(NPN6,NPN4,NPN4),
     &     IT11(NPN6,NPN4,NPN4),IT12(NPN6,NPN4,NPN4),
     &     IT21(NPN6,NPN4,NPN4),IT22(NPN6,NPN4,NPN4)
 
 
      COMMON /CT/ TR1,TI1
      COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
      P=DACOS(-1D0)
 
C  OPEN FILES *******************************************************
 
      OPEN (6,FILE='tmatr.test')
      OPEN (10,FILE='tmatr.write')
 
C  INPUT DATA ********************************************************
 
C      RAT=0.5 D0
C      NDISTR=3
C      AXMAX=1D0
C      NPNAX=2
C      B=0.1D0
C      GAM=0.5D0
C      NKMAX=5
C      EPS=2D0
C      NP=-1
C      LAM=0.5D0
C      MRR=1.53 d0
C      MRI=0.008D0
C      DDELT=0.001D0
C      NPNA=19
C      NDGS=2

      NCHECK=0
      IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK=1
      IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1
      WRITE (6,5454) NCHECK
 5454 FORMAT ('NCHECK=',I1)
      DAX=AXMAX/NPNAX
      IF (ABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
      if (ABS(RAT-1D0).GT.1D-8.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
      IF (ABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
C     PRINT 8000, RAT
 8000 FORMAT ('RAT=',F8.6)
      IF(NP.EQ.-1.AND.EPS.GE.1D0) PRINT 7000,EPS
      IF(NP.EQ.-1.AND.EPS.LT.1D0) PRINT 7001,EPS
      IF(NP.GE.0) PRINT 7100,NP,EPS
      IF(NP.EQ.-2.AND.EPS.GE.1D0) PRINT 7150,EPS
      IF(NP.EQ.-2.AND.EPS.LT.1D0) PRINT 7151,EPS
      PRINT 7400, LAM,MRR,MRI
      PRINT 7200, DDELT
 7000 FORMAT('RANDOMLY ORIENTED OBLATE SPHEROIDS, A/B=',F11.7)
 7001 FORMAT('RANDOMLY ORIENTED PROLATE SPHEROIDS, A/B=',F11.7)
 7100 FORMAT('RANDOMLY ORIENTED CHEBYSHEV PARTICLES, T',
     &       I1,'(',F5.2,')')
 7150 FORMAT('RANDOMLY ORIENTED OBLATE CYLINDERS, D/L=',F11.7)
 7151 FORMAT('RANDOMLY ORIENTED PROLATE CYLINDERS, D/L=',F11.7)
 7200 FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',D8.2)
 7400 FORMAT('LAM=',F10.6,3X,'MRR=',D10.4,3X,'MRI=',D10.4)
      DDELT=0.1D0*DDELT
      DO 600 IAX=1,NPNAX
         AXI=AXMAX-DAX*DFLOAT(IAX-1)
         R1=0.89031D0*AXI
         R2=1.56538D0*AXI
         NK=INT(AXI*NKMAX/AXMAX+2)
         IF (NK.GT.1000) PRINT 8001,NK
         IF (NK.GT.1000) STOP
         IF (NDISTR.EQ.3) CALL POWER (AXI,B,R1,R2)
 8001    FORMAT ('NK=',I4,' I.E., IS GREATER THAN 1000. ',
     &           'EXECUTION TERMINATED.')
         CALL GAUSS (NK,0,0,XG,WG)
         Z1=(R2-R1)*0.5D0
         Z2=(R1+R2)*0.5D0
         Z3=R1*0.5D0
         IF (NDISTR.EQ.5) GO TO 3
         DO I=1,NK
            XG1(I)=Z1*XG(I)+Z2
            WG1(I)=WG(I)*Z1
         ENDDO
         GO TO 4
    3    DO I=1,NK
            XG1(I)=Z3*XG(I)+Z3
            WG1(I)=WG(I)*Z3
         ENDDO
         DO I=NK+1,2*NK
            II=I-NK
            XG1(I)=Z1*XG(II)+Z2
            WG1(I)=WG(II)*Z1
         ENDDO
         NK=NK*2
    4    CALL DISTRB (NK,XG1,WG1,NDISTR,AXI,B,GAM,R1,R2,
     &                REFF,VEFF,P)
         PRINT 8002,R1,R2
 8002    FORMAT('R1=',F10.6,'   R2=',F10.6)
         IF (ABS(RAT-1D0).LE.1D-6) PRINT 8003, REFF,VEFF
         IF (ABS(RAT-1D0).GT.1D-6) PRINT 8004, REFF,VEFF
 8003    FORMAT('EQUAL-VOLUME-SPHERE REFF=',F8.4,'   VEFF=',F7.4)
 8004    FORMAT('EQUAL-SURFACE-AREA-SPHERE REFF=',F8.4,
     &          '   VEFF=',F7.4)
         PRINT 7250,NK
 7250    FORMAT('NUMBER OF GAUSSIAN QUADRATURE POINTS ',
     &          'IN SIZE AVERAGING =',I4)
         DO I=1,NPL
            ALPH1(I)=0D0
            ALPH2(I)=0D0
            ALPH3(I)=0D0
            ALPH4(I)=0D0
            BET1(I)=0D0
            BET2(I)=0D0
            ALPHA(I,1)=0D0
            ALPHA(I,2)=0D0
            ALPHA(I,3)=0D0
            ALPHA(I,4)=0D0
            BETA(I,1)=0D0
            BETA(I,2)=0D0
         ENDDO      
         CSCAT=0D0
         CEXTIN=0D0
         L1MAX=0
         DO 500 INK=1,NK
            I=NK-INK+1
            A=RAT*XG1(I)
            XEV=2D0*P*A/LAM
            IXXX=XEV+4.05D0*XEV**0.333333D0
            INM1=MAX0(4,IXXX)
            IF (INM1.GE.NPN1) PRINT 7333, NPN1
            IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,  
     &       '.  EXECUTION TERMINATED')
            QEXT1=0D0
            QSCA1=0D0
            DO 50 NMA=INM1,NPN1
               NMAX=NMA
               MMAX=1
               NGAUSS=NMAX*NDGS
               IF (NGAUSS.GT.NPNG1) PRINT 7340, NGAUSS
               IF (NGAUSS.GT.NPNG1) STOP
 7340          FORMAT('NGAUSS =',I3,' I.E. IS GREATER THAN NPNG1.',
     &                '  EXECUTION TERMINATED')
 7334          FORMAT(' NMAX =', I3,'  DC2=',D8.2,'   DC1=',D8.2)
 7335 FORMAT('                              NMAX1 =', I3,'  DC2=',D8.2,
     &      '  DC1=',D8.2)
               CALL CONST(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
               CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,
     &                   DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &                      DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0D0
               QSCA=0D0
               DO N=1,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                          +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
               ENDDO
               DSCA=ABS((QSCA1-QSCA)/QSCA)
               DEXT=ABS((QEXT1-QEXT)/QEXT)
C              PRINT 7334, NMAX,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               NMIN=DFLOAT(NMAX)/2D0+1D0
               DO 10 N=NMIN,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  DQSCA=DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                      +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  DQEXT=(TR1NN+TR1NN1)*DN1
                  DQSCA=ABS(DQSCA/QSCA)
                  DQEXT=ABS(DQEXT/QEXT)
                  NMAX1=N
                  IF (DQSCA.LE.DDELT.AND.DQEXT.LE.DDELT) GO TO 12
   10          CONTINUE
   12          CONTINUE
c              PRINT 7335, NMAX1,DQSCA,DQEXT
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 55
               IF (NMA.EQ.NPN1) PRINT 7333, NPN1
               IF (NMA.EQ.NPN1) STOP      
   50       CONTINUE
   55       NNNGGG=NGAUSS+1
            IF (NGAUSS.EQ.NPNG1) PRINT 7336
            MMAX=NMAX1
            DO 150 NGAUS=NNNGGG,NPNG1
               NGAUSS=NGAUS
               NGGG=2*NGAUSS
 7336          FORMAT('WARNING: NGAUSS=NPNG1')
 7337          FORMAT(' NG=',I3,'  DC2=',D8.2,'   DC1=',D8.2)
               CALL CONST(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
               CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,
     &                   DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &                      DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0D0
               QSCA=0D0
               DO 104 N=1,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                          +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
  104          CONTINUE
               DSCA=ABS((QSCA1-QSCA)/QSCA)
               DEXT=ABS((QEXT1-QEXT)/QEXT)
c              PRINT 7337, NGGG,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 155
               IF (NGAUS.EQ.NPNG1) PRINT 7336
  150       CONTINUE
  155       CONTINUE
            QSCA=0D0
            QEXT=0D0
            NNM=NMAX*2
            DO 204 N=1,NNM
               QEXT=QEXT+TR1(N,N)
  204       CONTINUE
            IF (NMAX1.GT.NPN4) PRINT 7550, NMAX1
 7550       FORMAT ('NMAX1 = ',I3, ', i.e. greater than NPN4.',
     &              ' Execution terminated')
            IF (NMAX1.GT.NPN4) STOP              
            DO 213 N2=1,NMAX1
               NN2=N2+NMAX
               DO 213 N1=1,NMAX1
                  NN1=N1+NMAX
                  ZZ1=TR1(N1,N2)
                  RT11(1,N1,N2)=ZZ1
                  ZZ2=TI1(N1,N2)
                  IT11(1,N1,N2)=ZZ2
                  ZZ3=TR1(N1,NN2)
                  RT12(1,N1,N2)=ZZ3
                  ZZ4=TI1(N1,NN2)
                  IT12(1,N1,N2)=ZZ4
                  ZZ5=TR1(NN1,N2)
                  RT21(1,N1,N2)=ZZ5
                  ZZ6=TI1(NN1,N2)
                  IT21(1,N1,N2)=ZZ6
                  ZZ7=TR1(NN1,NN2)
                  RT22(1,N1,N2)=ZZ7
                  ZZ8=TI1(NN1,NN2)
                  IT22(1,N1,N2)=ZZ8
                  QSCA=QSCA+ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &                 +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8
  213       CONTINUE
C           PRINT 7800,0,ABS(QEXT),QSCA,NMAX
            DO 220 M=1,NMAX1
               CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &                    DDR,DRR,DRI,NMAX,NCHECK)
               NM=NMAX-M+1
               NM1=NMAX1-M+1
               M1=M+1
               QSC=0D0
               DO 214 N2=1,NM1
                  NN2=N2+M-1
                  N22=N2+NM
                  DO 214 N1=1,NM1
                     NN1=N1+M-1
                     N11=N1+NM
                     ZZ1=TR1(N1,N2)
                     RT11(M1,NN1,NN2)=ZZ1
                     ZZ2=TI1(N1,N2)
                     IT11(M1,NN1,NN2)=ZZ2
                     ZZ3=TR1(N1,N22)
                     RT12(M1,NN1,NN2)=ZZ3
                     ZZ4=TI1(N1,N22)
                     IT12(M1,NN1,NN2)=ZZ4
                     ZZ5=TR1(N11,N2)
                     RT21(M1,NN1,NN2)=ZZ5
                     ZZ6=TI1(N11,N2)
                     IT21(M1,NN1,NN2)=ZZ6
                     ZZ7=TR1(N11,N22)
                     RT22(M1,NN1,NN2)=ZZ7
                     ZZ8=TI1(N11,N22)
                     IT22(M1,NN1,NN2)=ZZ8
                     QSC=QSC+(ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &                       +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8)*2D0
  214          CONTINUE
               NNM=2*NM
               QXT=0D0
               DO 215 N=1,NNM
                  QXT=QXT+TR1(N,N)*2D0
  215          CONTINUE
               QSCA=QSCA+QSC
               QEXT=QEXT+QXT
C              PRINT 7800,M,ABS(QXT),QSC,NMAX
 7800          FORMAT(' m=',I3,'  qxt=',d12.6,'  qsc=',d12.6,
     &                '  nmax=',I3)
  220       CONTINUE
            COEFF1=LAM*LAM*0.5D0/P
            CSCA=QSCA*COEFF1
            CEXT=-QEXT*COEFF1
c           PRINT 7880, NMAX,NMAX1
 7880       FORMAT ('nmax=',I3,'   nmax1=',I3)
            CALL GSP (NMAX1,CSCA,LAM,AL1,AL2,AL3,AL4,BE1,BE2,LMAX)
            L1M=LMAX+1
            L1MAX=MAX(L1MAX,L1M)
            WGII=WG1(I)
            WGI=WGII*CSCA
            DO 250 L1=1,L1M
               ALPH1(L1)=ALPH1(L1)+AL1(L1)*WGI
               ALPH2(L1)=ALPH2(L1)+AL2(L1)*WGI
               ALPH3(L1)=ALPH3(L1)+AL3(L1)*WGI
               ALPH4(L1)=ALPH4(L1)+AL4(L1)*WGI
               BET1(L1)=BET1(L1)+BE1(L1)*WGI
               BET2(L1)=BET2(L1)+BE2(L1)*WGI
  250       CONTINUE
            CSCAT=CSCAT+WGI
            CEXTIN=CEXTIN+CEXT*WGII
C           PRINT 6070, I,NMAX,NMAX1,NGAUSS
 6070       FORMAT(4I6)
  500    CONTINUE
         DO 510 L1=1,L1MAX
            ALPH1(L1)=ALPH1(L1)/CSCAT
            ALPH2(L1)=ALPH2(L1)/CSCAT
            ALPH3(L1)=ALPH3(L1)/CSCAT
            ALPH4(L1)=ALPH4(L1)/CSCAT
            BET1(L1)=BET1(L1)/CSCAT
            BET2(L1)=BET2(L1)/CSCAT
            ALPHA(L1,1)=ALPH1(L1)
            ALPHA(L1,2)=ALPH2(L1)
            ALPHA(L1,3)=ALPH3(L1)
            ALPHA(L1,4)=ALPH4(L1)
            BETA(L1,1)=BET1(L1)
            BETA(L1,2)=BET2(L1)
  510    CONTINUE
         WALB=CSCAT/CEXTIN
         CALL HOVENR(L1MAX,ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2)
         ASYMM=ALPH1(2)/3D0
         PRINT 9100,CEXTIN,CSCAT,WALB,ASYMM
 9100    FORMAT('CEXT=',D12.6,2X,'CSCA=',D12.6,2X,
     &          2X,'W=',D12.6,2X,'<COS>=',D12.6)
         IF (WALB.GT.1D0) PRINT 9111
 9111    FORMAT ('WARNING: W IS GREATER THAN 1')
         WRITE (10,580) WALB,L1MAX
         DO L=1,L1MAX
            WRITE (10,575) ALPH1(L),ALPH2(L),ALPH3(L),ALPH4(L),
     &                     BET1(L),BET2(L)
         ENDDO   
  575    FORMAT(6D14.7)
  580    FORMAT(D14.8,I8)
         LMAX=L1MAX-1
         CALL MATR (ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,LMAX,NPNA,F)
  600 CONTINUE
      ITIME=MCLOCK()
      TIME=DFLOAT(ITIME)/6000D0
      PRINT 1001,TIME
 1001 FORMAT (' time =',F8.2,' min')
      RETURN
      END
 
C**********************************************************************
 
C**********************************************************************
 
C**********************************************************************
 
C**********************************************************************
 
C**********************************************************************
 
C**********************************************************************
 
C**********************************************************************
 
C**********************************************************************

C**********************************************************************
C                                                                     *
C   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       *
C   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  *
C   BY USING BACKWARD RECURSION. PARAMETR NNMAX DETERMINES NUMERICAL  *
C   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))                *
C                                                                     *
C**********************************************************************
 
C**********************************************************************

C**********************************************************************
 
C**********************************************************************
 
C********************************************************************
C                                                                   *
C   CALCULATION OF THE EXPANSION COEFFICIENTS FOR (I,Q,U,V) -       *
C   REPRESENTATION.                                                 *
C                                                                   *
C   INPUT PARAMETERS:                                               *
C                                                                   *
C      LAM - WAVELENGTH OF LIGHT                                    *
C      CSCA - SCATTERING CROSS SECTION                              *
C      TR AND TI - ELEMENTS OF THE T-MATRIX. TRANSFERRED THROUGH    *
C                  COMMON /CTM/                                     *
C      NMAX - DIMENSION OF T(M)-MATRICES                            *
C                                                                   *
C   OUTPUT INFORTMATION:                                            *
C                                                                   *
C      ALF1,...,ALF4,BET1,BET2 - EXPANSION COEFFICIENTS             *
C      LMAX - NUMBER OF COEFFICIENTS MINUS 1                        *
C                                                                   *
C********************************************************************
 
      SUBROUTINE GSP(NMAX,CSCA,LAM,ALF1,ALF2,ALF3,ALF4,BET1,BET2,LMAX)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-B,D-H,O-Z),COMPLEX*16 (C)
      REAL*8 LAM,SSIGN(900)
      REAL*8  CSCA,SSI(NPL),SSJ(NPN1),
     &        ALF1(NPL),ALF2(NPL),ALF3(NPL),
     &        ALF4(NPL),BET1(NPL),BET2(NPL),
     &        TR1(NPL1,NPN4),TR2(NPL1,NPN4),
     &        TI1(NPL1,NPN4),TI2(NPL1,NPN4),
     &        G1(NPL1,NPN6),G2(NPL1,NPN6),
     &        AR1(NPN4),AR2(NPN4),AI1(NPN4),AI2(NPN4),
     &        FR(NPN4,NPN4),FI(NPN4,NPN4),FF(NPN4,NPN4)
      REAL*4 B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4),
     &       B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4),
     &       D1(NPL1,NPN4,NPN4),D2(NPL1,NPN4,NPN4),
     &       D3(NPL1,NPN4,NPN4),D4(NPL1,NPN4,NPN4),
     &       D5R(NPL1,NPN4,NPN4),D5I(NPL1,NPN4,NPN4),
     &       PLUS1(NPN6*NPN4*NPN4*8)         
      REAL*4
     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
      COMPLEX*16 CIM(NPN1)
 
      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
      COMMON /CBESS/ B1R,B1I,B2R,B2I    
      COMMON /SS/ SSIGN
      EQUIVALENCE ( PLUS1(1),TR11(1,1,1) )
      EQUIVALENCE (D1(1,1,1),PLUS1(1)),        
     &            (D2(1,1,1),PLUS1(NPL1*NPN4*NPN4+1)),
     &            (D3(1,1,1),PLUS1(NPL1*NPN4*NPN4*2+1)),
     &            (D4(1,1,1),PLUS1(NPL1*NPN4*NPN4*3+1)), 
     &            (D5R(1,1,1),PLUS1(NPL1*NPN4*NPN4*4+1)) 
 
      CALL FACT
      CALL SIGNUM
      LMAX=2*NMAX
      L1MAX=LMAX+1
      CI=(0D0,1D0)
      CIM(1)=CI
      DO 2 I=2,NMAX
         CIM(I)=CIM(I-1)*CI
    2 CONTINUE
      SSI(1)=1D0
      DO 3 I=1,LMAX
         I1=I+1
         SI=DFLOAT(2*I+1)
         SSI(I1)=SI
         IF(I.LE.NMAX) SSJ(I)=DSQRT(SI)
    3 CONTINUE
      CI=-CI
      DO 5 I=1,NMAX
         SI=SSJ(I)
         CCI=CIM(I)
         DO 4 J=1,NMAX
            SJ=1D0/SSJ(J)
            CCJ=CIM(J)*SJ/CCI
            FR(J,I)=CCJ
            FI(J,I)=CCJ*CI
            FF(J,I)=SI*SJ
    4    CONTINUE
    5 CONTINUE
      NMAX1=NMAX+1
 
C *****  CALCULATION OF THE ARRAYS B1 AND B2  *****
 
      K1=1
      K2=0
      K3=0
      K4=1
      K5=1
      K6=2
 
C     PRINT 3300, B1,B2
 3300 FORMAT (' B1 AND B2')
      DO 100 N=1,NMAX
 
C *****  CALCULATION OF THE ARRAYS T1 AND T2  *****
 
 
         DO 10 NN=1,NMAX
            M1MAX=MIN0(N,NN)+1
            DO 6 M1=1,M1MAX
               M=M1-1
               L1=NPN6+M
               TT1=TR11(M1,N,NN)
               TT2=TR12(M1,N,NN)
               TT3=TR21(M1,N,NN)
               TT4=TR22(M1,N,NN)
               TT5=TI11(M1,N,NN)
               TT6=TI12(M1,N,NN)
               TT7=TI21(M1,N,NN)
               TT8=TI22(M1,N,NN)
               T1=TT1+TT2
               T2=TT3+TT4
               T3=TT5+TT6
               T4=TT7+TT8
               TR1(L1,NN)=T1+T2
               TR2(L1,NN)=T1-T2
               TI1(L1,NN)=T3+T4
               TI2(L1,NN)=T3-T4
               IF(M.EQ.0) GO TO 6
               L1=NPN6-M
               T1=TT1-TT2
               T2=TT3-TT4
               T3=TT5-TT6
               T4=TT7-TT8
               TR1(L1,NN)=T1-T2
               TR2(L1,NN)=T1+T2
               TI1(L1,NN)=T3-T4
               TI2(L1,NN)=T3+T4
    6       CONTINUE
   10    CONTINUE
 
C  *****  END OF THE CALCULATION OF THE ARRAYS T1 AND T2  *****
 
         NN1MAX=NMAX1+N
         DO 40 NN1=1,NN1MAX
            N1=NN1-1
 
C  *****  CALCULATION OF THE ARRAYS A1 AND A2  *****
 
            CALL CCG(N,N1,NMAX,K1,K2,G1)
            NNMAX=MIN0(NMAX,N1+N)
            NNMIN=MAX0(1,IABS(N-N1))
            KN=N+NN1
            DO 15 NN=NNMIN,NNMAX
               NNN=NN+1
               SIG=SSIGN(KN+NN)
               M1MAX=MIN0(N,NN)+NPN6
               AAR1=0D0
               AAR2=0D0
               AAI1=0D0
               AAI2=0D0
               DO 13 M1=NPN6,M1MAX
                  M=M1-NPN6
                  SSS=G1(M1,NNN)
                  RR1=TR1(M1,NN)
                  RI1=TI1(M1,NN)
                  RR2=TR2(M1,NN)
                  RI2=TI2(M1,NN)
                  IF(M.EQ.0) GO TO 12
                  M2=NPN6-M
                  RR1=RR1+TR1(M2,NN)*SIG
                  RI1=RI1+TI1(M2,NN)*SIG
                  RR2=RR2+TR2(M2,NN)*SIG
                  RI2=RI2+TI2(M2,NN)*SIG
   12             AAR1=AAR1+SSS*RR1
                  AAI1=AAI1+SSS*RI1
                  AAR2=AAR2+SSS*RR2
                  AAI2=AAI2+SSS*RI2
   13          CONTINUE
               XR=FR(NN,N)
               XI=FI(NN,N)
               AR1(NN)=AAR1*XR-AAI1*XI
               AI1(NN)=AAR1*XI+AAI1*XR
               AR2(NN)=AAR2*XR-AAI2*XI
               AI2(NN)=AAR2*XI+AAI2*XR
   15       CONTINUE
 
C  *****  END OF THE CALCULATION OF THE ARRAYS A1 AND A2 ****
 
            CALL CCG(N,N1,NMAX,K3,K4,G2)
            M1=MAX0(-N1+1,-N)
            M2=MIN0(N1+1,N)
            M1MAX=M2+NPN6
            M1MIN=M1+NPN6
            DO 30 M1=M1MIN,M1MAX
               BBR1=0D0
               BBI1=0D0
               BBR2=0D0
               BBI2=0D0
               DO 25 NN=NNMIN,NNMAX
                  NNN=NN+1
                  SSS=G2(M1,NNN)
                  BBR1=BBR1+SSS*AR1(NN)
                  BBI1=BBI1+SSS*AI1(NN)
                  BBR2=BBR2+SSS*AR2(NN)
                  BBI2=BBI2+SSS*AI2(NN)
   25          CONTINUE
               B1R(NN1,M1,N)=BBR1
               B1I(NN1,M1,N)=BBI1
               B2R(NN1,M1,N)=BBR2
               B2I(NN1,M1,N)=BBI2
   30       CONTINUE
   40    CONTINUE
  100 CONTINUE
 
C  *****  END OF THE CALCULATION OF THE ARRAYS B1 AND B2 ****
 
C  *****  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5  *****
 
C     PRINT 3301
 3301 FORMAT(' D1, D2, ...')
      DO 200 N=1,NMAX
         DO 190 NN=1,NMAX
            M1=MIN0(N,NN)
            M1MAX=NPN6+M1
            M1MIN=NPN6-M1
            NN1MAX=NMAX1+MIN0(N,NN)
            DO 180 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD1=0D0
               DD2=0D0
               DO 150 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B1R(NN1,M1,NN)
                  X4=B1I(NN1,M1,NN)
                  X5=B2R(NN1,M1,N)
                  X6=B2I(NN1,M1,N)
                  X7=B2R(NN1,M1,NN)
                  X8=B2I(NN1,M1,NN)
                  DD1=DD1+XX*(X1*X3+X2*X4)
                  DD2=DD2+XX*(X5*X7+X6*X8)
  150          CONTINUE
               D1(M1,NN,N)=DD1
               D2(M1,NN,N)=DD2
  180       CONTINUE
            MMAX=MIN0(N,NN+2)
            MMIN=MAX0(-N,-NN+2)
            M1MAX=NPN6+MMAX
            M1MIN=NPN6+MMIN
            DO 186 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD3=0D0
               DD4=0D0
               DD5R=0D0
               DD5I=0D0
               M2=-M+2+NPN6
               DO 183 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B2R(NN1,M1,N)
                  X4=B2I(NN1,M1,N)
                  X5=B1R(NN1,M2,NN)
                  X6=B1I(NN1,M2,NN)
                  X7=B2R(NN1,M2,NN)
                  X8=B2I(NN1,M2,NN)
                  DD3=DD3+XX*(X1*X5+X2*X6)
                  DD4=DD4+XX*(X3*X7+X4*X8)
                  DD5R=DD5R+XX*(X3*X5+X4*X6)
                  DD5I=DD5I+XX*(X4*X5-X3*X6)
  183          CONTINUE
               D3(M1,NN,N)=DD3
               D4(M1,NN,N)=DD4
               D5R(M1,NN,N)=DD5R
               D5I(M1,NN,N)=DD5I
  186       CONTINUE
  190    CONTINUE
  200 CONTINUE
 
C  *****  END OF THE CALCULATION OF THE D-ARRAYS *****
 
C  *****  CALCULATION OF THE EXPANSION COEFFICIENTS *****
 
C     PRINT 3303
 3303 FORMAT (' G1, G2, ...')
 
      DK=LAM*LAM/(4D0*CSCA*DACOS(-1D0))
      DO 300 L1=1,L1MAX
         G1L=0D0
         G2L=0D0
         G3L=0D0
         G4L=0D0
         G5LR=0D0
         G5LI=0D0
         L=L1-1
         SL=SSI(L1)*DK
         DO 290 N=1,NMAX
            NNMIN=MAX0(1,IABS(N-L))
            NNMAX=MIN0(NMAX,N+L)
            IF(NNMAX.LT.NNMIN) GO TO 290
            CALL CCG(N,L,NMAX,K1,K2,G1)
            IF(L.GE.2) CALL CCG(N,L,NMAX,K5,K6,G2)
            NL=N+L
            DO 280  NN=NNMIN,NNMAX
               NNN=NN+1
               MMAX=MIN0(N,NN)
               M1MIN=NPN6-MMAX
               M1MAX=NPN6+MMAX
               SI=SSIGN(NL+NNN)
               DM1=0D0
               DM2=0D0
               DO 270 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  IF(M.GE.0) SSS1=G1(M1,NNN)
                  IF(M.LT.0) SSS1=G1(NPN6-M,NNN)*SI
                  DM1=DM1+SSS1*D1(M1,NN,N)
                  DM2=DM2+SSS1*D2(M1,NN,N)
  270          CONTINUE
               FFN=FF(NN,N)
               SSS=G1(NPN6+1,NNN)*FFN
               G1L=G1L+SSS*DM1
               G2L=G2L+SSS*DM2*SI
               IF(L.LT.2) GO TO 280
               DM3=0D0
               DM4=0D0
               DM5R=0D0
               DM5I=0D0
               MMAX=MIN0(N,NN+2)
               MMIN=MAX0(-N,-NN+2)
               M1MAX=NPN6+MMAX
               M1MIN=NPN6+MMIN
               DO 275 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  SSS1=G2(NPN6-M,NNN)
                  DM3=DM3+SSS1*D3(M1,NN,N)
                  DM4=DM4+SSS1*D4(M1,NN,N)
                  DM5R=DM5R+SSS1*D5R(M1,NN,N)
                  DM5I=DM5I+SSS1*D5I(M1,NN,N)
  275          CONTINUE
               G5LR=G5LR-SSS*DM5R
               G5LI=G5LI-SSS*DM5I
               SSS=G2(NPN4,NNN)*FFN
               G3L=G3L+SSS*DM3
               G4L=G4L+SSS*DM4*SI
  280       CONTINUE
  290    CONTINUE
         G1L=G1L*SL
         G2L=G2L*SL
         G3L=G3L*SL
         G4L=G4L*SL
         G5LR=G5LR*SL
         G5LI=G5LI*SL
         ALF1(L1)=G1L+G2L
         ALF2(L1)=G3L+G4L
         ALF3(L1)=G3L-G4L
         ALF4(L1)=G1L-G2L
         BET1(L1)=G5LR*2D0
         BET2(L1)=G5LI*2D0
         LMAX=L
         IF(ABS(G1L).LT.1D-6) GO TO 500
  300 CONTINUE
  500 CONTINUE
      RETURN
      END
 
C****************************************************************
 
C   CALCULATION OF THE QUANTITIES F(N+1)=0.5*LN(N!)
C   0.LE.N.LE.899
 
      SUBROUTINE FACT
      REAL*8 F(900)
      COMMON /FAC/ F
      F(1)=0D0
      F(2)=0D0
      DO 2 I=3,900
         I1=I-1
         F(I)=F(I1)+0.5D0*DLOG(DFLOAT(I1))
    2 CONTINUE
      RETURN
      END
 
C************************************************************
 
C   CALCULATION OF THE ARRAY SSIGN(N+1)=SIGN(N)
C   0.LE.N.LE.899
 
      SUBROUTINE SIGNUM
      REAL*8 SSIGN(900)
      COMMON /SS/ SSIGN
      SSIGN(1)=1D0
      DO 2 N=2,899 
         SSIGN(N)=-SSIGN(N-1)
    2 CONTINUE
      RETURN
      END
 
C******************************************************************
C
C   CALCULATION OF CLEBSCH-GORDAN COEFFICIENTS
C   (N,M:N1,M1/NN,MM)
C   FOR GIVEN N AND N1. M1=MM-M, INDEX MM IS FOUND FROM M AS
C   MM=M*K1+K2
C
C   INPUT PARAMETERS :  N,N1,NMAX,K1,K2
C                               N.LE.NMAX
C                               N.GE.1
C                               N1.GE.0
C                               N1.LE.N+NMAX
C   OUTPUT PARAMETERS : GG(M+NPN6,NN+1) - ARRAY OF THE CORRESPONDING
C                                       COEFFICIENTS
C                               /M/.LE.N
C                               /M1/=/M*(K1-1)+K2/.LE.N1
C                               NN.LE.MIN(N+N1,NMAX)
C                               NN.GE.MAX(/MM/,/N-N1/)
C   IF K1=1 AND K2=0, THEN 0.LE.M.LE.N
 
 
      SUBROUTINE CCG(N,N1,NMAX,K1,K2,GG)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 GG(NPL1,NPN6),CD(0:NPN5),CU(0:NPN5)
      IF(NMAX.LE.NPN4.
     &   AND.0.LE.N1.
     &   AND.N1.LE.NMAX+N.
     &   AND.N.GE.1.
     &   AND.N.LE.NMAX) GO TO 1
      PRINT 5001
      STOP
 5001 FORMAT(' ERROR IN SUBROUTINE CCG')
    1 NNF=MIN0(N+N1,NMAX)
      MIN=NPN6-N
      MF=NPN6+N
      IF(K1.EQ.1.AND.K2.EQ.0) MIN=NPN6
      DO 100 MIND=MIN,MF
         M=MIND-NPN6
         MM=M*K1+K2
         M1=MM-M
         IF(IABS(M1).GT.N1) GO TO 90
         NNL=MAX0(IABS(MM),IABS(N-N1))
         IF(NNL.GT.NNF) GO TO 90
         NNU=N+N1
         NNM=(NNU+NNL)*0.5D0
         IF (NNU.EQ.NNL) NNM=NNL
         CALL CCGIN(N,N1,M,MM,C)
         CU(NNL)=C  
         IF (NNL.EQ.NNF) GO TO 50
         C2=0D0
         C1=C
         DO 7 NN=NNL+1,MIN0(NNM,NNF)
            A=DFLOAT((NN+MM)*(NN-MM)*(N1-N+NN))
            A=A*DFLOAT((N-N1+NN)*(N+N1-NN+1)*(N+N1+NN+1))
            A=DFLOAT(4*NN*NN)/A
            A=A*DFLOAT((2*NN+1)*(2*NN-1))
            A=DSQRT(A)
            B=0.5D0*DFLOAT(M-M1)
            D=0D0
            IF(NN.EQ.1) GO TO 5
            B=DFLOAT(2*NN*(NN-1))
            B=DFLOAT((2*M-MM)*NN*(NN-1)-MM*N*(N+1)+
     &               MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN-1)*(NN-1))
            D=D*DFLOAT((2*NN-3)*(2*NN-1))
            D=DFLOAT((NN-MM-1)*(NN+MM-1)*(N1-N+NN-1))/D
            D=D*DFLOAT((N-N1+NN-1)*(N+N1-NN+2)*(N+N1+NN))
            D=DSQRT(D)
    5       C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CU(NN)=C
    7    CONTINUE
         IF (NNF.LE.NNM) GO TO 50
         CALL DIRECT(N,M,N1,M1,NNU,MM,C)
         CD(NNU)=C
         IF (NNU.EQ.NNM+1) GO TO 50
         C2=0D0
         C1=C
         DO 12 NN=NNU-1,NNM+1,-1
            A=DFLOAT((NN-MM+1)*(NN+MM+1)*(N1-N+NN+1))
            A=A*DFLOAT((N-N1+NN+1)*(N+N1-NN)*(N+N1+NN+2))
            A=DFLOAT(4*(NN+1)*(NN+1))/A
            A=A*DFLOAT((2*NN+1)*(2*NN+3))
            A=DSQRT(A)
            B=DFLOAT(2*(NN+2)*(NN+1))
            B=DFLOAT((2*M-MM)*(NN+2)*(NN+1)-MM*N*(N+1)
     &               +MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN+2)*(NN+2))
            D=D*DFLOAT((2*NN+5)*(2*NN+3))
            D=DFLOAT((NN+MM+2)*(NN-MM+2)*(N1-N+NN+2))/D
            D=D*DFLOAT((N-N1+NN+2)*(N+N1-NN-1)*(N+N1+NN+3))
            D=DSQRT(D)
            C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CD(NN)=C
   12    CONTINUE
   50    DO 9 NN=NNL,NNF
            IF (NN.LE.NNM) GG(MIND,NN+1)=CU(NN)
            IF (NN.GT.NNM) GG(MIND,NN+1)=CD(NN)
c           WRITE (6,*) N,M,N1,M1,NN,MM,GG(MIND,NN+1)
    9    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
 
C*********************************************************************
 
      SUBROUTINE DIRECT (N,M,N1,M1,NN,MM,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(900)
      COMMON /FAC/ F
      C=F(2*N+1)+F(2*N1+1)+F(N+N1+M+M1+1)+F(N+N1-M-M1+1)    
      C=C-F(2*(N+N1)+1)-F(N+M+1)-F(N-M+1)-F(N1+M1+1)-F(N1-M1+1)
      C=DEXP(C)
      RETURN
      END
 
C*********************************************************************
C
C   CALCULATION OF THE CLEBCSH-GORDAN COEFFICIENTS
C   G=(N,M:N1,MM-M/NN,MM)
C   FOR GIVEN N,N1,M,MM, WHERE NN=MAX(/MM/,/N-N1/)
C                               /M/.LE.N
C                               /MM-M/.LE.N1
C                               /MM/.LE.N+N1
 
      SUBROUTINE CCGIN(N,N1,M,MM,G)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(900),SSIGN(900)
      COMMON /SS/ SSIGN
      COMMON /FAC/ F
      M1=MM-M
      IF(N.GE.IABS(M).
     &   AND.N1.GE.IABS(M1).
     &   AND.IABS(MM).LE.(N+N1)) GO TO 1
      PRINT 5001
      STOP
 5001 FORMAT(' ERROR IN SUBROUTINE CCGIN')
    1 IF (IABS(MM).GT.IABS(N-N1)) GO TO 100
      L1=N
      L2=N1
      L3=M
      IF(N1.LE.N) GO TO 50
      K=N
      N=N1
      N1=K
      K=M
      M=M1
      M1=K
   50 N2=N*2
      M2=M*2
      N12=N1*2
      M12=M1*2
      G=SSIGN(N1+M1+1)
     & *DEXP(F(N+M+1)+F(N-M+1)+F(N12+1)+F(N2-N12+2)-F(N2+2)
     &       -F(N1+M1+1)-F(N1-M1+1)-F(N-N1+MM+1)-F(N-N1-MM+1))
      N=L1
      N1=L2
      M=L3
      RETURN
  100 A=1D0
      L1=M
      L2=MM
      IF(MM.GE.0) GO TO 150
      MM=-MM
      M=-M
      M1=-M1
      A=SSIGN(MM+N+N1+1)
  150 G=A*SSIGN(N+M+1)
     &   *DEXP(F(2*MM+2)+F(N+N1-MM+1)+F(N+M+1)+F(N1+M1+1)
     &        -F(N+N1+MM+2)-F(N-N1+MM+1)-F(-N+N1+MM+1)-F(N-M+1)
     &        -F(N1-M1+1))
      M=L1
      MM=L2
      RETURN
      END
 
C********************************************************************
 
c********************************************************************
 
C********************************************************************
 
C********************************************************************
 
C  COMPUTATION OF R1 AND R2 FOR A POWER LAW SIZE DISTRIBUTION WITH
C  EFFECTIVE RADIUS A AND EFFECTIVE VARIANCE B
 
      SUBROUTINE POWER (A,B,R1,R2)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F
      COMMON AA,BB
      AA=A
      BB=B
      AX=1D-5
      BX=A-1D-5
      R1=ZEROIN (AX,BX,F,0D0)
      R2=(1D0+B)*2D0*A-R1
      RETURN
      END
 
C***********************************************************************
 
      DOUBLE PRECISION FUNCTION ZEROIN (AX,BX,F,TOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      EPS=1D0
   10 EPS=0.5D0*EPS
      TOL1=1D0+EPS
      IF (TOL1.GT.1D0) GO TO 10
   15 A=AX
      B=BX
      FA=F(A)
      FB=F(B)
   20 C=A
      FC=FA
      D=B-A
      E=D
   30 IF (ABS(FC).GE.ABS(FB)) GO TO 40
   35 A=B
      B=C
      C=A
      FA=FB
      FB=FC
      FC=FA
   40 TOL1=2D0*EPS*ABS(B)+0.5D0*TOL
      XM=0.5D0*(C-B)
      IF (ABS(XM).LE.TOL1) GO TO 90
   44 IF (FB.EQ.0D0) GO TO 90
   45 IF (ABS(E).LT.TOL1) GO TO 70
   46 IF (ABS(FA).LE.ABS(FB)) GO TO 70
   47 IF (A.NE.C) GO TO 50
   48 S=FB/FA
      P=2D0*XM*S
      Q=1D0-S
      GO TO 60
   50 Q=FA/FC
      R=FB/FC
      S=FB/FA
      P=S*(2D0*XM*Q*(Q-R)-(B-A)*(R-1D0))
      Q=(Q-1D0)*(R-1D0)*(S-1D0)
   60 IF (P.GT.0D0) Q=-Q
      P=ABS(P)
      IF ((2D0*P).GE.(3D0*XM*Q-ABS(TOL1*Q))) GO TO 70
   64 IF (P.GE.ABS(0.5D0*E*Q)) GO TO 70
   65 E=D
      D=P/Q
      GO TO 80
   70 D=XM
      E=D
   80 A=B
      FA=FB
      IF (ABS(D).GT.TOL1) B=B+D
      IF (ABS(D).LE.TOL1) B=B+DSIGN(TOL1,XM)
      FB=F(B)
      IF ((FB*(FC/ABS(FC))).GT.0D0) GO TO 20
   85 GO TO 30
   90 ZEROIN=B
      RETURN
      END
 
C***********************************************************************
 
      DOUBLE PRECISION FUNCTION F(R1)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON A,B
      R2=(1D0+B)*2D0*A-R1
      F=(R2-R1)/DLOG(R2/R1)-A
      RETURN
      END
 
C**********************************************************************
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
C    N - NUMBER OF POINTS                                             *
C    Z - DIVISION POINTS                                              *
C    W - WEIGHTS                                                      *
C**********************************************************************
 
C****************************************************************
 
      SUBROUTINE DISTRB (NNK,YY,WY,NDISTR,AA,BB,GAM,R1,R2,REFF,                 
     &                   VEFF,PI)                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 YY(NNK),WY(NNK)                                                    
      IF (NDISTR.EQ.2) GO TO 100                                                
      IF (NDISTR.EQ.3) GO TO 200                                                
      IF (NDISTR.EQ.4) GO TO 300                                                
      IF (NDISTR.EQ.5) GO TO 360
      PRINT 1001,AA,BB,GAM                                                      
 1001 FORMAT('MODIFIED GAMMA DISTRIBUTION, alpha=',F6.4,'  r_c=',               
     &  F6.4,'  gamma=',F6.4)                                                   
      A2=AA/GAM                                                                 
      DB=1D0/BB
      DO 50 I=1,NNK                                                             
         X=YY(I)                                                             
         Y=X**AA                                                                
         X=X*DB
         Y=Y*DEXP(-A2*(X**GAM))                                                 
         WY(I)=WY(I)*Y                                                       
   50 CONTINUE                                                                  
      GO TO 400                                                                 
  100 PRINT 1002,AA,BB                                                          
 1002 FORMAT('LOG-NORMAL DISTRIBUTION, r_g=',F8.4,                
     &       '  [ln(sigma_g)]**2=', F6.4)          
      DA=1D0/AA                                                                 
      DO 150 I=1,NNK                                                            
         X=YY(I)                                                                
         Y=DLOG(X*DA)                                                          
         Y=DEXP(-Y*Y*0.5D0/BB)/X                                             
         WY(I)=WY(I)*Y                                                          
  150 CONTINUE                                                                  
      GO TO 400                                                                 
  200 PRINT 1003                                                                
 1003 FORMAT('POWER LAW DISTRIBUTION OF HANSEN & TRAVIS 1974')                 
      DO 250 I=1,NNK                                                            
         X=YY(I)                                                                
         WY(I)=WY(I)/(X*X*X)                                                 
  250 CONTINUE                                                                  
      GO TO 400                                                                 
  300 PRINT 1004,AA,BB                                                          
 1004 FORMAT ('GAMMA DISTRIBUTION,  a=',F6.3,'  b=',F6.4)
      B2=(1D0-3D0*BB)/BB                                                        
      DAB=1D0/(AA*BB)                                                          
      DO 350 I=1,NNK                                                            
         X=YY(I)                                                                
         X=(X**B2)*DEXP(-X*DAB)                                                 
         WY(I)=WY(I)*X                                                       
  350 CONTINUE                                                                  
      GO TO 400                                                                 
  360 PRINT 1005,BB
 1005 FORMAT ('MODIFIED POWER LAW DISTRIBUTION,  alpha=',D10.4)
      DO 370 I=1,NNK
         X=YY(I)
         IF (X.LE.R1) WY(I)=WY(I)
         IF (X.GT.R1) WY(I)=WY(I)*(X/R1)**BB
  370 CONTINUE
  400 CONTINUE                                                                  
      SUM=0D0
      DO 450 I=1,NNK
         SUM=SUM+WY(I)
  450 CONTINUE
      SUM=1D0/SUM
      DO 500 I=1,NNK
         WY(I)=WY(I)*SUM
  500 CONTINUE
      G=0D0
      DO 550 I=1,NNK
         X=YY(I)
         G=G+X*X*WY(I)
  550 CONTINUE
      REFF=0D0
      DO 600 I=1,NNK
         X=YY(I)
         REFF=REFF+X*X*X*WY(I)
  600 CONTINUE
      REFF=REFF/G
      VEFF=0D0
      DO 650 I=1,NNK
         X=YY(I)
         XI=X-REFF
         VEFF=VEFF+XI*XI*X*X*WY(I)
  650 CONTINUE
      VEFF=VEFF/(G*REFF*REFF)
      RETURN                                                                    
      END                                                                       
 
C*************************************************************
 
      SUBROUTINE HOVENR(L1,A1,A2,A3,A4,B1,B2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
      DO 100 L=1,L1
         KONTR=1
         LL=L-1
         DL=DFLOAT(LL)*2D0+1D0
         DDL=DL*0.48D0
         AA1=A1(L)
         AA2=A2(L)
         AA3=A3(L)
         AA4=A4(L)
         BB1=B1(L)
         BB2=B2(L)
         IF(LL.GE.1.AND.ABS(AA1).GE.DL) KONTR=2
         IF(ABS(AA2).GE.DL) KONTR=2
         IF(ABS(AA3).GE.DL) KONTR=2
         IF(ABS(AA4).GE.DL) KONTR=2
         IF(ABS(BB1).GE.DDL) KONTR=2
         IF(ABS(BB2).GE.DDL) KONTR=2
         IF(KONTR.EQ.2) PRINT 3000,LL
         C=-0.1D0
         DO 50 I=1,11
            C=C+0.1D0
            CC=C*C
            C1=CC*BB2*BB2
            C2=C*AA4
            C3=C*AA3
            IF((DL-C*AA1)*(DL-C*AA2)-CC*BB1*BB1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL-C3)+C1.LE.-1D-4) KONTR=2
            IF((DL+C2)*(DL-C3)-C1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL+C3)-C1.LE.-1D-4) KONTR=2
            IF(KONTR.EQ.2) PRINT 4000,LL,C
   50    CONTINUE
  100 CONTINUE
      IF(KONTR.EQ.1) PRINT 2000
 2000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS SATISFIED')
 3000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3)
 4000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3,
     & '   A=',D9.2)
      RETURN
      END
 
C****************************************************************
 
C    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
C    COEFFICIENTS
 
C    A1,...,B2 - EXPANSION COEFFICIENTS
C    LMAX - NUMBER OF COEFFICIENTS MINUS 1
C    N - NUMBER OF SCATTERING ANGLES
C        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
C        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES
 
      SUBROUTINE MATR(A1,A2,A3,A4,B1,B2,LMAX,NPNA,F)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A1(NPL),A2(NPL),A3(NPL),A4(NPL),B1(NPL),B2(NPL),F(NPNA,4,4)
      N=NPNA
      DN=1D0/DFLOAT(N-1)
      DA=DACOS(-1D0)*DN
      DB=180D0*DN
      L1MAX=LMAX+1
      PRINT 1000
 1000 FORMAT(' ')
      PRINT 1001
 1001 FORMAT(' ',2X,'S',6X,'ALPHA1',6X,'ALPHA2',6X,'ALPHA3',
     &       6X,'ALPHA4',7X,'BETA1',7X,'BETA2')
      DO 10 L1=1,L1MAX
         L=L1-1
         PRINT 1002,L,A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
   10 CONTINUE
 1002 FORMAT(' ',I3,6F12.5)
      TB=-DB
      TAA=-DA
      PRINT 1000
      PRINT 1003
 1003 FORMAT(' ',5X,'<',8X,'F11',8X,'F22',8X,'F33',
     & 8X,'F44',8X,'F12',8X,'F34')
      D6=DSQRT(6D0)*0.25D0
      DO 500 I1=1,N
         TAA=TAA+DA
         TB=TB+DB
         U=DCOS(TAA)
         F11=0D0
         F2=0D0
         F3=0D0
         F44=0D0
         F12=0D0
         F34=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11=F11+A1(L1)*PP1
            F44=F44+A4(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12=F12+B1(L1)*PP4
            F34=F34+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DFLOAT(L*L1)*U
            PL3=DFLOAT(L1*(L*L-4))
            PL4=1D0/DFLOAT(L*(L1*L1-4))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DFLOAT(L*L-4))*P4)/DSQRT(DFLOAT(L1*L1-4))
            P4=PP4
            PP4=P
  400    CONTINUE
         F22=(F2+F3)*0.5D0
         F33=(F2-F3)*0.5D0
C        F22=F22/F11
C        F33=F33/F11
C        F44=F44/F11
C        F12=-F12/F11
C        F34=F34/F11
         PRINT 1004,TB,F11,F22,F33,F44,F12,F34
         F(I1,1,1)=F11
         F(I1,1,2)=F12
         F(I1,1,3)=0
         F(I1,1,4)=0
         F(I1,2,1)=F12
         F(I1,2,2)=F22
         F(I1,2,3)=0
         F(I1,2,4)=0
         F(I1,3,1)=0
         F(I1,3,2)=0
         F(I1,3,3)=F33
         F(I1,3,4)=F34
         F(I1,4,1)=0
         F(I1,4,2)=0
         F(I1,4,3)=-F34
         F(I1,4,4)=F44
  500 CONTINUE
      PRINT 1000 
 1004 FORMAT(' ',F6.2,6F11.4)
      RETURN
      END
