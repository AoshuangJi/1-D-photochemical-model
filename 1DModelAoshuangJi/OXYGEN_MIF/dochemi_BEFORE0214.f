      SUBROUTINE DOCHEMI(FVALI,N,USOLI)
      
!      INCLUDE 'parameters_add_CO2.txt'
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
      INCLUDE '../INCLUDECHEM/parNEQ_LDA.inc'
      INCLUDE '../INCLUDECHEM/parNR.inc'
      INCLUDE '../INCLUDECHEM/parNSP_NSP1_NSP2.inc'
      INCLUDE '../INCLUDECHEM/parNMAX.inc'
      INCLUDE '../INCLUDECHEM/parNF.inc'
      INCLUDE '../INCLUDECHEM/comNBLOK.inc'  

      INCLUDE '../INCLUDECHEM/comABLOK.inc' !NEW:JTROP DEL:PMOD
      INCLUDE '../INCLUDECHEM/comBBLOK.inc' 
      INCLUDE '../INCLUDECHEM/comCBLOK.inc' !ZTROP 
      INCLUDE '../INCLUDECHEM/comFBLOK1.inc' 
      INCLUDE '../INCLUDECHEM/comGBLOK.inc'
      INCLUDE '../INCLUDECHEM/comSATBLK.inc' ! for H2OSAT
      INCLUDE '../INCLUDECHEM/comSULBLK.inc'
      INCLUDE '../INCLUDECHEM/comLTBLOK.inc' !for ZAPNO
      PARAMETER(NQI=31, NQI1=NQI+1, NSLS=38) 
      PARAMETER(NRI=353, NSPI=78, NSPI1=NSPI+1, NSPI2=NSPI+2,NMAXI=70)

      DIMENSION FVALI(NQI,NZ),XP(NZ),XL(NZ),USOLI(NQI,NZ)
      REAL ZAPO2I(JTROP),LIGHTNO2I
!      COMMON/ABLOK/EDD(NZ),DEN(NZ),DK(NQT,NZ),Z(NZ),T(NZ),G,FSCALE,
!     2  ALB,DELZ,BOVERH,DD(NQT,NZ),DL(NQT,NZ),DU(NQT,NZ),DI(NQT,NZ),
!     3  HI(NQT,NZ),DHU(NQT,NZ),DHL(NQT,NZ),PMOD
     
      COMMON/ISOTOP/AI(NRI,NZ),ILOSSI(2,NSPI,NMAXI),
     2  JCHEMI(5,NRI),NUMLI(NSPI),NUMPI(NSPI),TPI(NQI),TLI(NQI),
     3  YP(NQI,NZ),YL(NQI,NZ),SRI(NQI),TLOSSI(NQI),PHIDEPI(NQI),
     4 ISPECI(NSPI2),DIZ(NSPI2,NZ),IPRODI(NSPI,NMAXI),PSO4AER
     5 ,LBOUNDI(NQI),FLOWI(NQI),FUPI(NQI),CONI(NQI)
C
      COMMON/RATESIS/RAINGCI(NQI,NZ),DDI(NQI,NZ),DKI(NQI,NZ),
     2  DUI(NQI,NZ),DLI(NQI,NZ),DHUI(NQI,NZ),DHLI(NQI,NZ),HIZ(NQI,NZ),
     3 VDEPI(NQI),RATI(NRI)
     
!      COMMON/FBLOK/REL(NQ,NZ),MBOUND(NQT),LBOUND(NQT),PHIDEP(NQT),
!     2  TLOSS(NQT),HBUG(NQT),HBUG2(NQT),HBUG3(NQT),HCOEFF(NQT),
!     3  H2CHEM,H2SURF,H2VOLC,VEFF(NQ),COIMPACT,CH4VOLC,H2SVOLC,S2VOLC,
!     4  CO2VOLC,FEFLUX
    
!      COMMON/SULBLK/VH2O(NF,NZ),VH2SO4(NF,NZ),FTAB(NF),H2SO4S(NZ),
!     2  FSULF(NZ),CONSO4(NZ)

C      COMMON/NIBLOK/LH2CQ,LQ,LH2Q,LQH,LHOQ,LH2OQ,LO2Q,LOQO,LCQ,LCH3OQH,
C     2  LCH3QOH,LCH3OQ,LCH3QO,LN2Q,LNQ,LNOQ,LHNOQ,LHNO2Q,
C     3  LHOQNO2,LHO2NOQ,
C     4  LNO2Q,LN2O4Q,LN2QO4,LOQ,LSQ,LSOQ,LH2SO3Q,LHSQ,LCOQ,LQ1D, 
C     5  LH3CQ,LHCQ,LSO1Q,LSO3Q,LHSO2Q,LSO2Q,LIH2CO,LIO,LIH2O,LIOH,
C     6  LIHO2,LIH2O2,LIO3,LIH,LIH2,LICH4,LICO,LICH3OOH,LICH3O2,LIN2O,
C     7  LINO,LINO2,LIHNO3,LIHO2NO2,LINO3,LIN2O5,LIO2,LIH2S,LIHS,LISO,
C     8  LISO2,LIHSO,LIO1D,LICH21,LICH23,LICH3,LIH3CO,LIHCO,LIN,LIS,
C     9  LISO21,LISO23,LIHSO3,LISO3,LIN2

C   01/23/2020 JL ADD HQO HQONO2 TO THE NIBLOK

      COMMON/NIBLOK/LH2CQ,LQ,LH2Q,LQH,LHOQ,LH2OQ,LO2Q,LOQO,LCQ,LCH3OQH,
     2  LCH3QOH,LCH3OQ,LCH3QO,LN2Q,LNQ,LNOQ,LHNOQ,LHNO2Q,
     3  LHOQNO2,LHO2NOQ,LNO2Q,LN2O4Q,LN2QO4,LOQ,LSQ,
     4  LSOQ,LH2SO3Q,LHSQ,LCOQ,LHQO,LHQONO2,LQ1D, 
     5  LH3CQ,LHCQ,LSO1Q,LSO3Q,LHSO2Q,LSO2Q,LIH2CO,LIO,LIH2O,LIOH,
     6  LIHO2,LIH2O2,LIO3,LIH,LIH2,LICH4,LICO,LICH3OOH,LICH3O2,LIN2O,
     7  LINO,LINO2,LIHNO3,LIHO2NO2,LINO3,LIN2O5,LIO2,LIH2S,LIHS,LISO,
     8  LISO2,LIHSO,LICO2,LIO1D,LICH21,LICH23,LICH3,LIH3CO,LIHCO,LIN,
     9  LIS,LISO21,LISO23,LIHSO3,LISO3,LIN2
      INTEGER N

!      print *, 'CHECKK2'
!      PRINT *, (DIZ(I,1),I=30,36)
C   CHEMICAL EQUILIBRIUM SPECIES ARE DONE FIRST.  THESE MUST CON-
C   TAIN NO NONLINEARITIES (SUCH AS S4 REACTING WITH ITSELF TO FORM
C   S8) AND MUST BE DONE IN THE PROPER ORDER (I.E. IF SPECIES A
C   REACTS TO FORM B, THEN A MUST BE FOUND FIRST).  LONG-LIVED
C   SPECIES CAN BE DONE IN ANY ORDER.
C
      DO I=1,NQI
       DO J=1,NZ
        FVALI(I,J) = 0.
       ENDDO
      ENDDO 
      DO 1 I=1,NQI
      DO 1 J=1,NZ
   1  DIZ(I,J) = USOLI(I,J) * DEN(J)

C                                                         
C
C   THIS SUBROUTINE DOES THE CHEMISTRY BY CALLING CHEMPLI.  PHOTO-
C ***** SOLVE FOR THE PHOTOCHEMICAL EQUILIBRIUM SPECIES *****
C
                                                         
      DO 3 I=NQI1,NSLS                                     
      CALL CHEMPLI(XP,XL,I)                              
      DO 3 J=1,NZ                                         
   3  DIZ(I,J) = XP(J)/XL(J)

 
C
C ***** LONG-LIVED SPECIES CHEMISTRY *****
      DO 4 I=1,NQI
      CALL CHEMPLI(XP,XL,I)


C-PL IN ORDER TO UPDATE THE NUMBER DENSITY FOR LONG-LIVED SPECIES
!      DO  J=1,NZ                                         
!      DIZ(I,J) = XP(J)/XL(J)
!      END DO 
C-PL BUT WE CAN NOT UPDATE THE DIZ HERE, SINCE WE USE THIS VARIABLE 
C-PL AT THE END OF THIS SUBROUTINE TO CAL. TPI TLI.SO WHEN STEP=1,
C-PL WE ONLY UPDATE THE SHORT-LIVED SPECIES
C-AP
 
 
      DO 4 J=1,NZ
      XLJ = XL(J) + RAINGCI(I,J)
      FVALI(I,J) = XP(J)/DEN(J) - XLJ*USOLI(I,J)
      YP(I,J) = XP(J)
   4  YL(I,J) = XLJ
 
CCCCCCCCC
C-PL ADD FOR H2O QUESTION WHETHER WE NEED ! PART 06/03
C-AP ***********************************************
C   ZERO OUT H2O TERMS IN THE TROPOSPHERE AND INCLUDE LIGHTNING
C   PRODUCTION OF NO, O2, AND CO.  (MUST INCLUDE CO IN ORDER TO
C   BALANCE THE HYDROGEN BUDGET)
C-AS This is a leftover of the low O2, high CO2 version of the 
c    photochemical code, all the commented lines on this do loop
c    are not nedeed for this highO2 version of the code

      DO J=1,JTROP
      FVALI(LH2Q,J) = 0.
      SCALE = RAIN(J)/RAIN(1)
      ZAP = ZAPNO * SCALE
      FVALI(LNQ,J) = FVALI(LNQ,J) + ZAP/DEN(J)
      YP(LNQ,J) = YP(LNQ,J) + ZAP
C-PL ONLY CONSIDER N2+O2=2NO   
      FVALI(LOQ,J) = FVALI(LOQ,J) - 2.*ZAP/DEN(J)   !0.5
      YL(LOQ,J) = YL(LOQ,J) + 2.*ZAP/DEN(J) !0.5  
      ZAPO2I(J) = 2.*ZAP/DEN(J) !0.5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      FVAL(LOQ,J) = FVAL(LOQ,J) - ZAP/DEN(J)
!      YL(LOQ,J) = YL(LOQ,J) + ZAP !C-PL ACTIVE THIS TWO LINE FOR LOW O2 CAL'

c     ZAP = ZAPO2 * SCALE
c     FVAL(LO2,J) = FVAL(LO2,J) + ZAP/DEN(J)
c     YP(LO2,J) = YP(LO2,J) + ZAP
c      ZAP = (ZAPNO + 2.*ZAPO2)*SCALE
c      FVAL(LCO,J) = FVAL(LCO,J) + ZAP/DEN(J)
c      YP(LCO,J) = YP(LCO,J) + ZAP
      enddo
!QUESTION whether need this? 06/11 won't change the result too much
! so comment out. Discuss with Jim
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  VOLCANIC OUTGASSING OF SOQ
C    Note: PRSO2 = 3.5E9, distributed over lower 10 km.
      PRSO2 = 2.*3.5E9/1.E6 !C-PL 06/13 you need times two for SOQ HERE
!      PRH2S = 3.5E8/1.e6 !h2s
      JTEN = 1.E6/DELZ + 0.01
      DO 351 J=1,JTEN
      FVALI(LSOQ,J) = FVALI(LSOQ,J) + PRSO2/DEN(J)
!      FVALI(LIH2S,J)= FVALI(LIH2S,J) + PRH2S/DEN(J)
!      YP(LIH2S,J) = YP(LIH2S,J) + PRH2S
 351  YP(LSOQ,J) = YP(LSOQ,J) + PRSO2
C
C-PL ADD FOR H2O QUESTION WHETHER WE NEED HERE
CQUESTION WHETHER WE NEED ! PART 06/03
C   H2O CONDENSATION IN THE STRATOSPHERE
C   (RHCOLD IS THE ASSUMED RELATIVE HUMIDITY AT THE COLD TRAP
      RHCOLD = 0.1
      JT1 = JTROP + 1
      CONFAC = 1.6E-5
      DO 13 J=JT1,NZ
      H2OCRT = RHCOLD * H2OSAT(J)
      IF (USOLI(LH2Q,J) .LT. H2OCRT) GO TO 13
      CONDEN(J) = CONFAC * (USOLI(LH2Q,J) - H2OCRT)
      FVALI(LH2Q,J) = FVALI(LH2Q,J) - CONDEN(J)
  13  CONTINUE

!QUESTION whether need this?06/11 won't change the result too much
! so comment out. Discuss with Jim
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-PL eliminate fei O2/CO2  07/2019
!      CORRO2 =  2.*1.7E14/1.E5!2.2521E12/1.E5
      GPPI = 2.* GPPOXY/1.E5 /10. !Divide GPP into 10 layers 
      GPPIC= 2.* GPPCDE/1.E5 /10.
C REMOVED
!      DO J=1,10
!      FVALI(LOQ,J) = FVALI(LOQ,J) - GPPI/DEN(J)
!      YL(LOQ,J) = YL(LOQ,J) + GPPI/DIZ(LOQ,J)
!      END DO
C INJECTED
      DO J=1,10
      FVALI(LOQ,J) = FVALI(LOQ,J) + GPPI/DEN(J)
      YP(LOQ,J) = YP(LOQ,J) + GPPI

      FVALI(LCOQ,J) = FVALI(LCOQ,J) + GPPIC/DEN(J)
      YP(LCOQ,J) = YP(LCOQ,J) + GPPIC
!      FVALI(LCOQ,J) = FVALI(LCOQ,J) - GPPIC/DEN(J)
!      YL(LCOQ,J) = YL(LCOQ,J) + GPPIC/DIZ(LCOQ,J)
      END DO


!      FVALI(LOQ,1) = FVALI(LOQ,1) - GPPI/DEN(1)!CORRO2/DEN(1)
!      YL(LOQ,1) = YL(LOQ,1) + GPPI/DIZ(LOQ,1)! CORRO2
!      FVALI(LCOQ,1) = FVALI(LCOQ,1) - GPPI/DEN(1)
!      YL(LCOQ,1) = YL(LCOQ,1) + GPPI/DIZ(LCOQ,1)
!      CORRCO2 = -2.*1.7E14/1.5E5/USOLI(LCOQ,1)!!1.1127E12/1.E5/USOLI(LCOQ,1)
!      FVALI(LCOQ,1) = FVALI(LCOQ,1) - GPP/DEN(1)!CORRCO2/DEN(1)
!      YL(LCOQ,1) = YL(LCOQ,1) + GPP/DIZ(LCOQ,1)!CORRCO2/DEN(1)

C   H2SO4 CONDENSATION
      CONFAC = 1.6E-5
      PSO4AER = 0.0
      DO 14 J=1,NZ 
      CONSO4(J) = CONFAC * (USOLI(LH2SO3Q,J) - (H2SO4S(J)*4.))
      CONSO4(J) = AMAX1(CONSO4(J),0.)
      FVALI(LH2SO3Q,J) = FVALI(LH2SO3Q,J) - CONSO4(J)
      YL(LH2SO3Q,J) = YL(LH2SO3Q,J) + CONFAC
      IF((H2SO4S(J)*4.).GT.USOLI(LH2SO3Q,J)) THEN !C-PL WHEN OVER SATURATION
      YP(LH2SO3Q,J) = YP(LH2SO3Q,J) + 0.0 
      ELSE
      YP(LH2SO3Q,J) = YP(LH2SO3Q,J) + CONFAC*4.*H2SO4S(J)*DEN(J)
      PSO4AER = PSO4AER + CONSO4(J)*DEN(J)*1.e5
      END IF
  14  CONTINUE

C-AP
C   SPECIAL TREATMENT FOR H2O AND NO
C-PL CALCULATE TP/L(H2O) IN TROPOSPHERE AND 
      TPH2OI = 0.
      TLH2OI = 0.
     
      DO J=1,JTROP
      TPH2OI = TPH2OI + YP(LH2Q,J)*1.e5
      TLH2OI = TLH2OI + YL(LH2Q,J)*DIZ(LH2Q,J)*1.e5
      END DO


C-PL GET RID OF 218-224 IN ORDER TO GET RIGHT TP/TL(H2O)
!      DO 5 J=1,NZ
!      IF(Z(J).GT.ZTROP) GO TO 7     
!      FVALI(LH2Q,J) = 0.
!      YP(LH2Q,J) = 0.
!      YL(LH2Q,J) = 1.E-99  
!   7  CONTINUE
!   5  CONTINUE

      IF(N.LT.1) RETURN
  
      DO 8 I=1,NQI
      TPI(I) = 0.
      TLI(I) = 0.
      DO 8 J=1,NZ
      TPI(I) = TPI(I) + YP(I,J)*1.e5
      TLI(I) = TLI(I) + YL(I,J)*DIZ(I,J)*1.e5
   8  CONTINUE

C-PL CALCULATE TL(O2) BY LIGHTNING
      LIGHTNO2I = 0.
      DO J=1,JTROP
      LIGHTNO2I = LIGHTNO2I + ZAPO2I(J)*DIZ(LOQ,J)*1.e5

      END DO

      RETURN
      END
