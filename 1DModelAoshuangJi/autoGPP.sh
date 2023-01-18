#!/bin/bash
#!/bin/awk
#set echo
#set verbose


#############
# WROTE BY PENG LIU AT JUL/2019 IN PENN STATE
# THIS SCRIPT CAN HELP YOU FIND THE APPROPRIATE GPP FOR BOTH O2 AND CO2
#                 RUN WIDE RANGES MIXING RATIO FOR BOTH O2 AND CO2
#
#############


t1=1000
depo2=0.0001     #1.e-4     O2 deposition velocity
depco2=0.00375   #3.75e-3  CO2 deposition velocity
OXYR=7 # YOU GOT 7 O2 MIXING RATIO
OXYPAL=(0.0 1.0 0.3 0.1 0.03 0.01 0.003 0.001) # O2 MIXING RATIO RANGE [1PAL-0.001PAL]
CO2R=6 # YOU GOT 6 CO2 MIXING RATIO
CO2PAL=(0.0 1.0 3.0 10.0 30.0 100.0 300.0)    # CO2 MIXING RATIO RANGE [1PAL-300PAL]

io2=1
while [ $io2 -le ${OXYR} ];do
#echo $io2
#echo ${OXYPAL[io2]}
# SET O2 MIXING RATIO
OXYM_1=$(grep "FO2=" ./IO/input_atmchem.dat)
OXYM_2=${OXYM_1//FO2=}
OXYM_3=${OXYM_2//!O2 mixing ratio}
OXYM_ORI="FO2="$OXYM_3 # THE ORIGINAL O2 MIXING RATIO IN input_atmchem.dat

OXY_CODE1=$(echo "${OXYPAL[io2]}" | awk '{printf("%.5f\n",$1*0.21)}') #CALCULATE NEW O2 MIXING RATIO
OXY_CODE="FO2=       "$OXY_CODE1"      "

sed -i "s/${OXYM_ORI}/${OXY_CODE}/g" ./IO/input_atmchem.dat  #SEND O2 MIXING RATIO TO THE MODEL
echo '#########'
echo 'OPERATING O2 MIXING RATIO='${OXYPAL[io2]}' PAL'
echo '#########'
# SET CO2 MIXING RATIO TO 1PAL AT THE BEGINNING OF EACH O2 LEVEL
CO2M_1=$(grep "FCO2=" ./IO/input_atmchem.dat)
CO2M_2=${CO2M_1//FCO2=}
CO2M_3=${CO2M_2//!CO2 surface mixing ratio}
CO2M_ORI="FCO2="$CO2M_3 # THE ORIGINAL CO2 MIXING RATIO IN input_atmchem.dat

CO2_CODE1=$(echo "${CO2PAL[1]}" | awk '{printf("%.5f\n",$1*0.0004)}') #CALCULATE NEW O2 MIXING RATIO
CO2_CODE2=$(echo $CO2_CODE1 |awk '{printf("%.2e\n",$0)}') 
CO2_CODE="FCO2=      "$CO2_CODE2"     "

sed -i "s/${CO2M_ORI}/${CO2_CODE}/g" ./IO/input_atmchem.dat  #SEND O2 MIXING RATIO TO THE MODEL
#echo 'OPERATING CO2 MIXING RATIO='${CO2PAL[1]}' PAL'
#############################################################
#RUN MODEL
if [ "$io2" -gt 1 ]
then
cp ./IO/history/O2_${OXYPAL[io2-1]}PAL_CO2_${CO2PAL[1]}PAL/atm_composition.out ./IO/atm_composition.dat
fi

./runchem > isotope.dat

# cal. vdep for O2
FL=$(grep "FLOW(O2)=" ./IO/outchem.dat)
ND=$(grep "SL(O2)="   ./IO/outchem.dat)
FLU=${FL//FLOW(O2)=}
NDS=${ND//SL(O2)=}
vdep=$(echo "scale=7;$FLU / $NDS"|bc) #calculate vdep according to the output


#################################################
# WHETHER YOU NEED FIND THE APPROPRIATE GPP FOR OXYGEN
if [ $(echo "$vdep < $depo2"|bc) -eq 1 ] 
then
echo 'You already use the right GPP for O2' # NO, YOU DON'T
#cp   isotope.dat ./IO/outchem.dat ./IO/atm_composition.dat  ./IO/atm_composition.out ./IO/int.iso_rates.out.dat  ./IO/int.rates.out.dat  ./IO/input_atmchem.dat ./IO/history/O2_${OXYPAL[io2]}PAL_CO2_1PAL 
else                # YES, YOU NEED

NGPP=$(echo "scale=4;$NDS * $depo2 / 1"|bc)   #calculate modest GPP
i=1
while [ $i -le ${t1} ];do
GPP_1=$(grep "GPPOXY=" ./IO/input_atmchem.dat)
GPP_2=${GPP_1//GPPOXY=}
GPP_3=${GPP_2//!Gross primary productivity (1.16e14 cm-2s-1)}
GPP_ORI="GPPOXY="$GPP_3

NGPPG=$(echo "scale=4;$NGPP * (1.0 + (($i - 1.0) * 0.01)) / 1"|bc) #increase GPP by 1%


NGPP_1=$(echo $NGPPG |awk '{printf("%.2e\n",$0)}')  #change to scientific notation 
NGPP_2=${NGPP_1//+}

NGPP_CODE[i]="GPPOXY=    "$NGPP_2"      "
#echo "${NGPP_CODE[i]}" 
#echo "${GPP_ORI}"
sed -i "s/${GPP_ORI}/${NGPP_CODE[i]}/g" ./IO/input_atmchem.dat

#run case
./runchem > isotope.dat

# check the new output
FLX=$(grep "FLOW(O2)=" ./IO/outchem.dat)
NDX=$(grep "SL(O2)="   ./IO/outchem.dat)
FLUX=${FLX//FLOW(O2)=}
NDSX=${NDX//SL(O2)=}
vdepX=$(echo "scale=7;$FLUX / $NDSX"|bc) #calculate UPDATED vdep according to the output

if [ $(echo "$vdepX < $depo2"|bc) -eq 1 ] 
then
echo 'keep finding for O2 GPP'
else  # the right GPP is ${NGPP_CODE[i-1]}
sed -i "s/${NGPP_CODE[i]}/${NGPP_CODE[i-1]}/g" ./IO/input_atmchem.dat
#run case
./runchem > isotope.dat
echo 'YOU GOT THE RIGHT O2 GPP!!!'
#save fils
#cp   isotope.dat ./IO/outchem.dat ./IO/atm_composition.dat  ./IO/atm_composition.out ./IO/int.iso_rates.out.dat  ./IO/int.rates.out.dat  ./IO/input_atmchem.dat ./IO/history/O2_${OXYPAL[io2]}PAL_CO2_1PAL 
break
fi
  i=$((i+1))
done
fi       # END FINDING THE GPP FOR OXYGEN
#################################################################
# CO2 LEVEL LOOP
ico2=1
while [ $ico2 -le ${CO2R} ];do
# SET CO2 MIXING RATIO 
CO2M_1=$(grep "FCO2=" ./IO/input_atmchem.dat)
CO2M_2=${CO2M_1//FCO2=}
CO2M_3=${CO2M_2//!CO2 surface mixing ratio}
CO2M_ORI="FCO2="$CO2M_3 # THE ORIGINAL CO2 MIXING RATIO IN input_atmchem.dat

CO2_CODE1=$(echo "${CO2PAL[ico2]}" | awk '{printf("%.5f\n",$1*0.0004)}') #CALCULATE NEW O2 MIXING RATIO
CO2_CODE2=$(echo $CO2_CODE1 |awk '{printf("%.2e\n",$0)}') 
CO2_CODE="FCO2=      "$CO2_CODE2"     "

sed -i "s/${CO2M_ORI}/${CO2_CODE}/g" ./IO/input_atmchem.dat  #SEND O2 MIXING RATIO TO THE MODEL
echo 'OPERATING CO2 MIXING RATIO='${CO2PAL[ico2]}' PAL'
# GPP_CO2 EQUALS TO MAX VALUE AT THE BEGINNING
GPPC_1=$(grep "GPPCDE=" ./IO/input_atmchem.dat)
GPPC_2=${GPPC_1//GPPCDE=}
GPPC_3=${GPPC_2//!Gross primary productivity CARBON DIOXIDE}
GPPC_ORI="GPPCDE="$GPPC_3

NGPPCoo=1.16e14
GPPC_ORIo="GPPCDE=    "$NGPPCoo"      "
sed -i "s/${GPPC_ORI}/${GPPC_ORIo}/g" ./IO/input_atmchem.dat

#RUN MODEL
./runchem > isotope.dat
# cal. vdep for CO2
FLC=$(grep "FLOW(CO2)=" ./IO/outchem.dat)
NDC=$(grep "SL(CO2)="   ./IO/outchem.dat)
FLUC=${FLC//FLOW(CO2)=}
NDSC=${NDC//SL(CO2)=}
vdepc=$(echo "scale=7;$FLUC / $NDSC"|bc) #calculate vdep according to the output
###################################################
# WHETHER YOU NEED FIND THE APPROPRIATE GPP FOR CO2
if [ $(echo "$vdepc < $depco2"|bc) -eq 1 ] 
then
echo 'Saving files' # NO, YOU DON'T
cp   isotope.dat ./IO/outchem.dat ./IO/atm_composition.dat  ./IO/atm_composition.out ./IO/int.iso_rates.out.dat  ./IO/int.rates.out.dat  ./IO/input_atmchem.dat ./IO/history/O2_${OXYPAL[io2]}PAL_CO2_${CO2PAL[ico2]}PAL 
else                # YES, YOU NEED

NGPPC=$(echo "scale=4;$NDSC * $depco2 / 1"|bc)   #calculate modest GPP_CO2
j=1
while [ $j -le ${t1} ];do
GPPC_1=$(grep "GPPCDE=" ./IO/input_atmchem.dat)
GPPC_2=${GPPC_1//GPPCDE=}
GPPC_3=${GPPC_2//!Gross primary productivity CARBON DIOXIDE}
GPPC_ORI="GPPCDE="$GPPC_3

NGPPGC=$(echo "scale=4;$NGPPC * (1.0 - (($j - 1.0) * 0.01)) / 1"|bc) #decrease modest GPP by 1%


NGPPC_1=$(echo $NGPPGC |awk '{printf("%.2e\n",$0)}')  #change to scientific notation 
NGPPC_2=${NGPPC_1//+}

NGPPC_CODE[j]="GPPCDE=    "$NGPPC_2"      "
#echo "${NGPPC_CODE[j]}" >>logg.dat
#echo "${GPPC_ORI}">>logg.dat
sed -i "s/${GPPC_ORI}/${NGPPC_CODE[j]}/g" ./IO/input_atmchem.dat

#run case
./runchem > isotope.dat

# check the new output
FLXC=$(grep "FLOW(CO2)=" ./IO/outchem.dat)
NDXC=$(grep "SL(CO2)="   ./IO/outchem.dat)
FLUXC=${FLXC//FLOW(CO2)=}
NDSXC=${NDXC//SL(CO2)=}
vdepXC=$(echo "scale=7;$FLUXC / $NDSXC"|bc) #calculate UPDATED vdep according to the output

if [ $(echo "$vdepXC > $depco2"|bc) -eq 1 ] 
then
echo 'keep finding for CO2 GPP'
else  # the right GPP is ${NGPP_CODE[i-1]}
#sed -i "s/${NGPPC_CODE[j]}/${NGPPC_CODE[j-1]}/g" ./IO/input_atmchem.dat
#run case
./runchem > isotope.dat
echo 'YOU GOT THE RIGHT CO2 GPP!!!'
echo '        '
#save fils
cp   isotope.dat ./IO/outchem.dat ./IO/atm_composition.dat  ./IO/atm_composition.out ./IO/int.iso_rates.out.dat  ./IO/int.rates.out.dat  ./IO/input_atmchem.dat ./IO/history/O2_${OXYPAL[io2]}PAL_CO2_${CO2PAL[ico2]}PAL 
mv  ./IO/atm_composition.out ./IO/atm_composition.dat
break
fi
  j=$((j+1))
done
fi
# END FINDING THE GPP FOR CO2
###############
  ico2=$((ico2+1))
done     # END THE CO2 RANGE LOOP
################

  io2=$((io2+1))
done     # END THE O2 RANGE LOOP

echo 'I have finished all the runs'
################























