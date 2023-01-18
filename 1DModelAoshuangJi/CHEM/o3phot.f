      SUBROUTINE O3PHOT
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comO3BLOK.inc'
C
C   THIS SUBROUTINE CALCULATES TEMPERATURE COEFFICIENTS USED TO FIND
C   THE O(1D) QUANTUM YIELD IN O3 PHOTOLYSIS.  (SEE JPL, 1983.)
C
      DO 1 I=1,NZ
      TAU = T(I) - 230.
      TAU2 = TAU * TAU
      TAU3 = TAU2 * TAU
      AT(I) =   0.332 + 2.565E-4*TAU + 1.152E-5*TAU2 + 2.313E-8*TAU3
      BT(I) = - 0.575 + 5.590E-3*TAU - 1.439E-5*TAU2 - 3.270E-8*TAU3
      CT(I) =   0.466 + 8.883E-4*TAU - 3.546E-5*TAU2 + 3.519E-7*TAU3
   1  ALM0(I) = 308.2 + 4.4871E-2*TAU + 6.938E-5*TAU2 - 2.5452E-6*TAU3
C

c      print *, ' BT(I) = ', BT(I) same 
c      print *, 'TAU =', TAU same
c      print *, 'TAU2 =', TAU2 same 
c      print *, 'TAU3 = ', TAU3 same 
c      print *, 'AT(I) = ', AT(I) 
c      print *, ' CT(I) = ', CT(I) 
c      print *, 'ALM0(I) = ', ALM0(I) 
      RETURN
      END
