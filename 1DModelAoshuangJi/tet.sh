#!/bin/bash
#!/bin/awk
#set echo
#set verbose



t1=5
i=1
ii=1.0

while [ $i -le ${t1} ];do

  a=$(echo "scale=1;$ii - 1.0"|bc)
#echo $a
  g=$(echo "scale=4;10.0 * (1.0 + ((i - 1.0) * 0.01)) / 1"|bc) 

#echo $g



  i=$((i+1))
  ii=$(echo "scale=1;$ii + 1.0"|bc)
done

a=1.0
#cp isotope.dat ./IO/history/O2_${a}PAL_CO2_1PAL/aa.dat

#b=$(printf "%.5f" 'echo "scale=5;0.21000 / 1.000000"|bc') #CALCULATE NEW O2 MIXING RATIO
A=5
B=16
#b=$(echo "$A $B" | awk '{printf("%.5f\n",$1/$2)}')

c=4.00e-4

b=$(echo "$c" | awk '{printf("%.5f\n",$1*3.0)}')
d=$(echo $b |awk '{printf("%.2e\n",$0)}')  
#echo $d

NGPPCoo=1.10e14
GPPC_ORIo="GPPCDE=    "$NGPPCoo"      "
echo "$GPPC_ORIo" >logg.dat

