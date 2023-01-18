
      SUBROUTINE RATES

       INCLUDE '../INCLUDECHEM/parNZ.inc'
       INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
       INCLUDE '../INCLUDECHEM/parNR.inc'
       INCLUDE '../INCLUDECHEM/parNSP_NSP1_NSP2.inc'
       INCLUDE '../INCLUDECHEM/parNMAX.inc'

      DIMENSION A3(NZ),A4(NZ),A5(NZ),A6(NZ),A10(NZ),A11(NZ),A12(NZ),
     2  A13(NZ),A14(NZ),A15(NZ),A16(NZ),A17(NZ),A18(NZ),A19(NZ),
     3  A20(NZ),A21(NZ),A22(NZ),A30(NZ),A31(NZ),A32(NZ),A41(NZ),
     4  A43(NZ),A46(NZ),A47(NZ),A48(NZ)
      
      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comCBLOK.inc'     
      INCLUDE '../INCLUDECHEM/comGBLOK.inc'
      INCLUDE '../INCLUDECHEM/comRBLOK.inc'



C
C ***** TEMPERATURE-DEPENDENT RATE COEFFICIENTS *****
C
 
      DO 1 I=1,NZ
      PATM = DEN(I)*1.38E-16*T(I)/1.013E6
      A3(I) = 3.44E-13*((T(I)/298.)**2.67)*EXP(-3160./T(I))
      A4(I) = 7.7E-12*EXP(-2100./T(I))
      A5(I) = 1.4E-10*EXP(-470./T(I))
      A10(I) = 2.4E-11*EXP(110./T(I))
      A12(I) = 1.7E-12*EXP(-940./T(I))
      A13(I) = 2.7E-11*EXP(225./T(I))
      A14(I) = 1.97E-16*((T(I)/298.)**4.57)*EXP(695./T(I))
      A15(I) = 2.3E-13*EXP(600./T(I)) + 1.7E-33*EXP(1000./T(I))
     2  *DEN(I)
      A16(I) = 2.9E-12*EXP(-160./T(I))
      A17(I) = 5.21E-35*EXP(900./T(I))*DEN(I)
      A19(I) = 8.0E-12 * EXP(-2060./T(I))
      A20(I) = 6.2E-14*((T(I)/298.)**2.60)*EXP(945./T(I))
      A21(I) = 1.8E-11*EXP(110./T(I))
      A22(I) = 3.2E-11*EXP(70./T(I))
      A30(I) = 1.5E-13 * (1. + 0.6*PATM)
      A31(I) = 1.7E-33*EXP(-1510./T(I))*DEN(I)
      A32(I) = 5.29E-34*EXP(-370./T(I))*DEN(I)
      A41(I) = 2.14E-12*((T(I)/298.)**1.62)*EXP(-1090./T(I))
      A43(I) = 6.04E-33*(T(I)/298.)*DEN(I)
      A46(I) = 4.38E-30*((T(I)/298.)**-2.00)*DEN(I)
      A48(I) = 3.4E-11*EXP(-1600./T(I))
   1  CONTINUE

C
C ***** THREE-BODY COEFFICIENTS *****
      DO 3 I=1,NZ
      TT = T(I)
      DN = DEN(I)
C     AR(J,I) = TBDY(K0,KI,N,M,T,DEN)
C
C  HO2
      AR(6,I) = TBDY(4.11E-32,7.51E-11,1.10,0.,TT,DN)
   
C
C  O3
      AR(18,I) = TBDY(5.70E-34,1.E-10,2.6,0.,TT,DN)
!      AR(18,I) = 6.0e-34*(TT/300.)**(-2.6)*DN
C
C  H2O2
      AR(47,I) = TBDY(6.9E-31,2.5E-11,0.8,0.,TT,DN)
C
C   CH3O2
      AR(67,I) = TBDY(4.5E-31,1.8E-12,3.,1.7,TT,DN)
C
C   NO2
      AR(84,I) = TBDY(9.E-32,3.E-11,1.5,0.,TT,DN)
C
C   HNO2
      AR(86,I) = TBDY(7.52E-31,3.3E-11,2.4,0.3,TT,DN)
C
C   HNO3
      AR(88,I) = TBDY(3.28E-30,2.7E-11,3.75,0.,TT,DN)
C
C   HO2NO2
      AR(91,I) = TBDY(1.8E-31,4.7E-12,3.2,0.,TT,DN)
C
C   NOCL
      AR(114,I) = TBDY(9.E-32,1.E-10,1.6,0.,TT,DN)
C
C   CLONO
      AR(115,I) = TBDY(1.3E-30,1.E-10,2.,1.,TT,DN)
C
C   CLO2
      AR(117,I) = TBDY(1.44E-33,1.E-10,3.9,0.,TT,DN)
c
C   CLONO2
      AR(122,I) = TBDY(1.64E-31,7.00E-11,3.4,0.,TT,DN)
C
C   CL2O2
      AR(141,I) = TBDY(2.05E-32,1.00E-11,4.00,0.,TT,DN)
C
C   NO3
      AR(150,I) = TBDY(1.31E-31,2.3E-11,1.50,0.24,TT,DN)
C      AR(104,I) = TBDY(1.31E-31,2.3E-11,1.50,0.24,TT,DN)
C
C   N2O5
C      AR(106,I) = TBDY(3.70E-30,1.9E-12,4.10,0.20,TT,DN)
      AR(152,I) = TBDY(3.70E-30,1.9E-12,4.10,0.20,TT,DN)
C
C   HSO3
C      AR(117,I) = TBDY(4.62E-31,1.31E-12,3.90,0.70,TT,DN)   
      AR(163,I) = TBDY(4.62E-31,1.31E-12,3.90,0.70,TT,DN) 
c   N2O
C      AR(172,I) = TBDY(3.5E-37,1.E-10,0.6,0.,TT,DN)  
      AR(218,I) = TBDY(3.5E-37,1.E-10,0.6,0.,TT,DN)            

   3  CONTINUE
C
C ***** FILL UP RATE MATRIX *****
      DO 4 I=1,NZ
      AR(1,I) = 2.62E-10*EXP(65./T(I))
      AR(2,I) = 1.10E-10
      AR(3,I) = A3(I)
      AR(4,I) = A4(I)
      AR(5,I) = A5(I)
      AR(7,I) = 5.6E-12
      AR(8,I) = 2.4E-12
      AR(9,I) = 7.2E-11
      AR(10,I) = A10(I)
      AR(11,I) = 5.00E-11
      AR(12,I) = A12(I)
      AR(13,I) = A13(I)
      AR(14,I) = A14(I)
      AR(15,I) = A15(I)
      AR(16,I) = A16(I)
      AR(17,I) = A17(I)
      AR(19,I) = A19(I)
      AR(20,I) = A20(I)
      AR(21,I) = A21(I)
      AR(22,I) = A22(I)
      AR(30,I) = A30(I)
      AR(31,I) = A31(I)
      AR(32,I) = A32(I)
      AR(33,I) = 1.50E-10
      AR(34,I) = 0.      ! zero out for avoiding dependency at high O2 run: HCO + HCO --> H2CO + CO
      AR(35,I) = 1.69E-10
      AR(36,I) = 5.00E-11
      AR(37,I) = 5.00E-11
      AR(40,I) = 1.0E-2
      AR(41,I) = A41(I)
      AR(43,I) = A43(I)
      AR(44,I) = 3.5E-12 * EXP(140./T(I))
      AR(45,I) = 8.20E-12*EXP(40./T(I))
      AR(46,I) = A46(I)
      AR(48,I) = A48(I)
      AR(49,I) = 1.4E-12 * EXP(-2000./T(I))
      AR(58,I) = 4.16E-13*((T(I)/298.)**2.18)*EXP(-1230./T(I))
      AR(59,I) = 1.3E-10
      AR(60,I) = 7.51E-12
      AR(61,I) = 5.9E-11
      AR(62,I) = 6.64E-14
      AR(63,I) = 8.8E-12
      AR(64,I) = 5.E-15
      AR(65,I) = 7.1E-12*EXP(-5051./T(I))
      AR(66,I) = 6.64E-14
      AR(68,I) = 2.59E-13*((T(I)/298.)**-0.53)*EXP(-5440./T(I))
      AR(69,I) = 1.1E-10
      AR(70,I) = 5.4E-12 * EXP(-220./T(I))
      AR(71,I) = 3.8E-13 * EXP(780./T(I))
      AR(72,I) = 7.40E-13 * EXP(-520./T(I))
      AR(73,I) = 2.8E-12 * EXP(285./T(I))
      AR(74,I) = 3.9E-14*EXP(-900./T(I))
      AR(75,I) = 1.E-14
   4  CONTINUE
c      print *, ' AR(1) =', AR(1,I)
C
      DO 5 I=1,NZ
      AR(76,I) = 3.2E-13
      AR(77,I) = 6.7E-11
      AR(78,I) = 4.9E-11
      AR(79,I) = 1.5E-11*EXP(-3600./T(I))
      AR(80,I) = 2.01E-16
      AR(81,I) = 4.70E-11
      AR(82,I) = 2.09E-11*EXP(100./T(I))
      AR(83,I) = 1.4E-12 * EXP(-1310./T(I))
      AR(85,I) = 3.60E-12 * EXP(270./T(I))
      AR(87,I) = 5.50E-12 * EXP(190./T(I))
      AR(89,I) = 1.47E-10
      AK0 = 7.2E-15*EXP(785./T(I))
      AK2 = 4.1E-16*EXP(1440./T(I))
      AK3M = 1.9E-33*EXP(725./T(I))*DEN(I)
      AR(90,I) = AK0 + AK3M/(1. + AK3M/AK2)
      AR(92,I) = 1.90E-12 * EXP(270./T(I))
      AR(93,I) = 7.8E-11 * EXP(-3400./T(I))
      AR(94,I) = AR(91,I)/(2.33E-27*EXP(10870./T(I)))
      AR(96,I) = 1.9E-12 * EXP(190./T(I))
      AR(97,I) = AR(11,I)
      AR(98,I) = 1.40E-13*EXP(-2470./T(I))
      AR(99,I) = 4.5E-14*EXP(-1260./T(I))
      AR(100,I) = 1.E-11

c in may 7 2019, JL reassign order to the rxn after rip off the cl
c-AJ 08/30/22 ADD CL BACK
      AR(102,I) = 1.80E-11*EXP(110./T(I))
      AR(103,I) = 2.E-11
c in May 7, 2019 JL rip off all Cl rxn. the rxn number is the old one.
c-AJ 08/30/22 ADD CL BACK
      AR(104,I) = 2.4E-12*EXP(-1250./T(I))
      AR(105,I) = 2.80E-11*EXP(-250./T(I))
      AR(106,I) = 3.90E-11*EXP(-2310./T(I))
      AR(107,I) = 8.24E-13*((T(I)/298.)**2.49)*EXP(-610./T(I))
      AR(108,I) = 3.2E-11*EXP(-1250./T(I)) 
      AR(109,I) = 8.20E-11*EXP(-35./T(I))
      AR(110,I) = 1.1E-11*EXP(-980./T(I))
      AR(111,I) = 1.8E-11*EXP(170./T(I))
      AR(112,I) = 6.30E-11*EXP(-570./T(I))
      AR(113,I) = 6.51E-12*EXP(135./T(I))
      AR(116,I) = 5.8E-11*EXP(100./T(I))
      AR(118,I) = 2.3E-10
      AR(119,I) = 1.2E-11
      AR(120,I) = 2.5E-11*EXP(110./T(I))
      AR(121,I) = 6.2E-12*EXP(295./T(I))
      AR(123,I) = 4.8E-13*EXP(700./T(I))
      AR(124,I) = 1.1E-11*EXP(120./T(I))
      AR(125,I) = 1.70E-12*EXP(-230./T(I))
      AR(126,I) = 3.E-12*EXP(-500./T(I))
      AR(127,I) = 1.2E-12*EXP(-330./T(I))
      AR(128,I) = 1.E-11*EXP(-3300./T(I))
      AR(129,I) = 1.70E-13
      AR(130,I) = 4.5E-12*EXP(-900./T(I))
      AR(131,I) = 3.6E-12*EXP(-1200./T(I))
      AR(139,I) = AR(117,I)/(2.43E-25*EXP(2979./T(I)))
C      AR(103,I) = 1.91E-12
      AR(140,I) = 1.91E-12
   5  CONTINUE
C
      DO 6 I=1,NZ
      AR(143,I) = AR(141,I)/(3.E-27*EXP(8450./T(I)))
      AR(144,I) = AR(117,I)/(5.7E-25*EXP(2500./T(I)))
      AR(145,I) = 2.40E-11
      AR(146,I) = 2.5E-12*EXP(-130./T(I))
      AR(147,I) = 4.E-13
      AR(148,I) = 2.4E-12
      AR(149,I) = 2.4E-12*EXP(-960./T(I))

C      AR(108,I) = TBDY(1.3E-19,5.7E-14,0.,0.,TT,DN)
C      AR(109,I) = 2.E-19
      AR(154,I) = TBDY(1.3E-19,5.7E-14,0.,0.,TT,DN)
      AR(155,I) = 2.E-19
C  AR(155,I) should be viewed as a tuning parameter. The gas phase upper 
C  limit is 2.e-21. A value 100 times higher than this yields an effective
C  first-order rate of ~2.e-5 in the lower stratosphere.
 6    CONTINUE
C-AP
C-AP Adding sulfur rate constants
C
C-AP

c-AJ 08/30/22 ADD CL BACK
      DO 7 I=1,NZ
      AR(159,I) = 1.60E-13 * EXP(-2280./T(I))
      AR(160,I) = 2.8E-11
      AR(161,I) = 6.0E-31 * DEN(I)
      AR(162,I) = 8.6E-11
      AR(164,I) = 3.4E-32 * EXP(-1130./T(I)) * DEN(I)
      AR(165,I) = 6.0E-15
      AR(166,I) = 1.3E-12 *EXP(-330./T(I))
      AR(167,I) = 1.0E-11
      AR(168,I) = 1.0E-11
      AR(169,I) = 1.0E-11
      AR(170,I) = 6.10E-12 * EXP(-80./T(I))
      AR(171,I) = 3.66E-12*((T(I)/298.)**1.94)*EXP(-455./T(I))
      AR(172,I) = 9.22E-12 * EXP(-1800./T(I))
      AR(173,I) = 1.6E-10
      AR(174,I) = 4.0E-19
      AR(175,I) = 3.0E-11
      AR(176,I) = 1.2E-11
      AR(177,I) = 5.0E-11
      AR(178,I) = 1.0E-11
      AR(179,I) = 2.2E-11 * EXP(120./T(I))
      AR(180,I) = 2.10E-12
      AR(181,I) = 6.6E-11
c-AJ 08/30/22 ADD CL BACK
C      AR(113,I) = 1.60E-13 * EXP(-2280./T(I))
C      AR(114,I) = 2.8E-11
C      AR(115,I) = 6.0E-31 * DEN(I)
C      AR(116,I) = 8.6E-11
C      AR(118,I) = 3.4E-32 * EXP(-1130./T(I)) * DEN(I)
C      AR(119,I) = 6.0E-15
C      AR(120,I) = 1.3E-12 *EXP(-330./T(I))
C      AR(121,I) = 1.0E-11
C      AR(122,I) = 1.0E-11
C      AR(123,I) = 1.0E-11
C      AR(124,I) = 6.10E-12 * EXP(-80./T(I))
C      AR(125,I) = 3.66E-12*((T(I)/298.)**1.94)*EXP(-455./T(I))
C      AR(126,I) = 9.22E-12 * EXP(-1800./T(I))
C      AR(127,I) = 1.6E-10
C      AR(128,I) = 4.0E-19
C      AR(129,I) = 3.0E-11
C      AR(130,I) = 1.2E-11
C      AR(131,I) = 5.0E-11
C      AR(132,I) = 1.0E-11
C      AR(133,I) = 2.2E-11 * EXP(120./T(I))
C      AR(134,I) = 2.10E-12
C      AR(135,I) = 6.6E-11

C-AP We need to zero out this reaction S +HCO ->HS +CO because both S and HCO are short-lived
C-AP      AR(182,I) = 5.0E-11
C-KK	This reaction has been removed from the primo3s rxn list.
C-KK      AR(182,I) = 1.E-99

c-AJ 08/30/22 ADD CL BACK
      AR(183,I) = 1.5E-11
      AR(184,I) = 1.5E-11
      AR(185,I) = 1.7E-11 * EXP(-800./T(I))

C      AR(137,I) = 1.5E-11
C      AR(138,I) = 1.5E-11
C      AR(139,I) = 1.7E-11 * EXP(-800./T(I))
  
   7  CONTINUE
C-AP 
C

c-AJ 08/30/22 ADD CL BACK
c-AJ 08/30/22 O1D + CO2 = CO2 + O NOW IS AR(221)
C-AJ 08/30/22 CO + O1D = CO2 NOW IS AR(222)
      DO 8 I=1,NZ
      AR(190,I) = 1.0E-12
      AR(191,I) = 1.0E-11
      AR(192,I) = 1.5E+3
      AR(193,I) = 2.2E+4
      AR(194,I) = 1.0E-16
      AR(195,I) = 4.0E-12
      AR(196,I) = 1.5E-13
      AR(197,I) = 1.13E+3
      AR(198,I) = 7.0E-14
      AR(199,I) = 1.4E-11
      AR(200,I) = 4.50E-12 * EXP(-1170./T(I))
      AR(201,I) = 0.
      AR(202,I) = 9.50E-12 * EXP(-280./T(I))
      AR(203,I) = 2.9E-11 * EXP(240./T(I))
      AR(204,I) = 1.2E-11
      AR(205,I) = 8.3E-15
      AR(206,I) = 2.0E-15
      AR(207,I) = 1.0E-20
      AR(208,I) = 0.
      AR(209,I) = AR(44,I)
      AR(210,I) = AR(6,I)
      AR(212,I) = AR(11,I)
      AR(213,I) = AR(9,I)
      AR(214,I) = AR(7,I)
      AR(215,I) = 1.E-12
      AR(216,I) = AR(13,I)
      AR(217,I) = 1.E-11
      AR(219,I) = 9.22E-14*EXP(-2.99E-3/T(I))!5.03E-7*(298./T(I))**2.16*EXP(-18701./T(I))
      AR(220,I) = 2.E-14 * EXP(-25000./T(I))
      AR(221,I) = 6.9E-11*EXP(117./T(I))
      AR(222,I) = 0.!8.0E-11 !CO + O1D = CO2
c-AJ 08/30/22 ADD CL BACK
C      AR(144,I) = 1.0E-12
C      AR(145,I) = 1.0E-11
C      AR(146,I) = 1.5E+3
C      AR(147,I) = 2.2E+4
C      AR(148,I) = 1.0E-16
C      AR(149,I) = 4.0E-12
C      AR(150,I) = 1.5E-13
C      AR(151,I) = 1.13E+3
C      AR(152,I) = 7.0E-14
C      AR(153,I) = 1.4E-11
C      AR(154,I) = 4.50E-12 * EXP(-1170./T(I))
C      AR(155,I) = 0.
C      AR(156,I) = 9.50E-12 * EXP(-280./T(I))
C      AR(157,I) = 2.9E-11 * EXP(240./T(I))
C      AR(158,I) = 1.2E-11
C      AR(159,I) = 8.3E-15
C      AR(160,I) = 2.0E-15
C      AR(161,I) = 1.0E-20
C      AR(162,I) = 0.
C      AR(163,I) = AR(44,I)
C      AR(164,I) = AR(6,I)
C      AR(166,I) = AR(11,I)
C      AR(167,I) = AR(9,I)
C      AR(168,I) = AR(7,I)
C      AR(169,I) = 1.E-12
C      AR(170,I) = AR(13,I)
C      AR(171,I) = 1.E-11
C      AR(173,I) = 9.22E-14*EXP(-2.99E-3/T(I))!5.03E-7*(298./T(I))**2.16*EXP(-18701./T(I))
C      AR(174,I) = 2.E-14 * EXP(-25000./T(I))
C      AR(175,I) = 6.9E-11*EXP(117./T(I))
C      AR(176,I) = 0.!8.0E-11 !CO + O1D = CO2
   8  CONTINUE
      


C ***** GIORGI AND CHAMEIDES RAINOUT RATE *****
      GAM15 = 8.64E+05/2.0
      GAM8 = 7.0E+06/2.0
      AV = 6.02E+23
      WL = 1.0
      R = 1.36E-22

      NH = ZTROP/DZ + 0.01
c      print *,'NH =',NH,' ZTROP =',ZTROP,' DZ = ',DZ

C  Loop over altitude
c       print *,'In altitude loop in rainout routine'
       DO 10 I=1,NH
       ZKM = Z(I)/1.E5
       TEMP = T(I)
C       print *,'I =',I,' ZKM =',ZKM,' TEMP =',T(I)
C
C  Find appropriate GAMMA
      IF (ZKM.LE.1.51) THEN
         GAMMA = GAM15
      ELSE IF (ZKM.LT.8.) THEN
         GAMMA = GAM15 + (GAM8-GAM15)*((ZKM-1.5)/6.5)
      ELSE
         GAMMA = GAM8
      END IF
C
C  Find WH2O
      IF (ZKM.LE.1.) THEN
         X = 11.35 + 0.1*ZKM
      ELSE
         X = 11.5444 - 0.085333*ZKM - 9.1111E-03*ZKM*ZKM
      END IF
      WH2O = 10.0**X
C
C  Find F(Z)
      IF (ZKM.LE.1.51) THEN
         F = 0.1
      ELSE
         F = 0.16615 - 0.04916*ZKM + 3.37451E-03*ZKM*ZKM
      END IF

c      print *,'H(J) = ', H same 

c       print *, 'R =', R, same value, but print much more        
C
C  Loop over species
c in may 7 2019, JL makes changes on here 
       DO 10 J=1,NQ
c       print *, ' J =', J 
       RKJ = (WH2O/55.)/(AV*WL*1.0E-9 + 1./(H(J)*R*TEMP))
       QJ = 1. - F + F/(GAMMA*RKJ) * (1.0 - EXP(-RKJ*GAMMA))
   10  RAINGC(J,I) = (1. - EXP(-RKJ*GAMMA))/(GAMMA*QJ)
C
c       print *, 'Nq = ', NQ 
       NH1 = NH + 1
       DO 11 I=NH1,NZ
       DO 11 J=1,NQ
   11  RAINGC(J,I) = 0.


C
C ***** RAINOUT RATE *****
       DO 2 I=1,NZ
       ZKM = Z(I)/1.E5
       RAIN(I) = 0.
       IF(ZKM.LT.10.) RAIN(I) = 2.4E-6*EXP((6. - ZKM)/2.42)
       IF(ZKM.LT.6.) RAIN(I) = 2.4E-6
    2  CONTINUE

      RETURN

      END
