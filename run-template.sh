#!/bin/bash
#since Bash v4
LANG=en_US
#set -e
#path=$(pwd)
#echo path= $path
#/usr/local/lib/R/site-library/irace/bin/irace

problemsComplete=("DAFJS01.txt" "DAFJS02.txt" "DAFJS03.txt"
	"DAFJS04.txt" "DAFJS05.txt" "DAFJS06.txt" "DAFJS07.txt" "DAFJS08.txt"
	"DAFJS09.txt" "DAFJS10.txt" "DAFJS11.txt" "DAFJS12.txt" "DAFJS13.txt" 
	"DAFJS14.txt" "DAFJS15.txt" "DAFJS16.txt" "DAFJS17.txt" "DAFJS18.txt"
	"DAFJS19.txt" "DAFJS20.txt" "DAFJS21.txt" "DAFJS22.txt" "DAFJS23.txt" "DAFJS24.txt"
	"DAFJS25.txt" "DAFJS26.txt" "DAFJS27.txt" "DAFJS28.txt" "DAFJS29.txt" "DAFJS30.txt" "YFJS01.txt"
	"YFJS02.txt" "YFJS03.txt" "YFJS04.txt" "YFJS05.txt" "YFJS06.txt" "YFJS07.txt" "YFJS08.txt"
	"YFJS09.txt" "YFJS10.txt" "YFJS11.txt" "YFJS12.txt" "YFJS13.txt" "YFJS14.txt" "YFJS15.txt"
	"YFJS16.txt" "YFJS17.txt" "YFJS18.txt" "YFJS19.txt" "YFJS20.txt")


problemsclassic=("mfjs01.txt" "mfjs02.txt" "mfjs03.txt" "mfjs04.txt" "mfjs05.txt" "mfjs06.txt" "mfjs07.txt"
	"mfjs08.txt" "mfjs09.txt" "mfjs10.txt" "MK01.txt" "MK02.txt" "MK03.txt" "MK04.txt" "MK05.txt" "MK06.txt"
	"MK07.txt" "MK08.txt" "MK09.txt" "MK10.txt" "MK11.txt" "MK12.txt" "MK13.txt" "MK14.txt" "MK15.txt"
	"sfjs01.txt" "sfjs02.txt" "sfjs03.txt" "sfjs04.txt" "sfjs05.txt" "sfjs06.txt" "sfjs07.txt" "sfjs08.txt" "sfjs09.txt" "sfjs10.txt")

Hurink=("la01.txt"  "la06.txt"  "la11.txt"  "la16.txt"  "la21.txt"  "la26.txt"  "la31.txt"  "la36.txt"  "mt06.txt"
"la02.txt"  "la07.txt"  "la12.txt"  "la17.txt"  "la22.txt"  "la27.txt"  "la32.txt"  "la37.txt"  "mt10.txt"
"la03.txt"  "la08.txt"  "la13.txt"  "la18.txt"  "la23.txt"  "la28.txt"  "la33.txt"  "la38.txt"  "mt20.txt"
"la04.txt"  "la09.txt"  "la14.txt"  "la19.txt"  "la24.txt"  "la29.txt"  "la34.txt"  "la39.txt"
"la05.txt"  "la10.txt"  "la15.txt"  "la20.txt"  "la25.txt"  "la30.txt"  "la35.txt"  "la40.txt")

Barnes=("mt10c1.txt"  "mt10xx.txt"  "mt10xyz.txt"  "setb4x.txt"  "setb4xy.txt"  "seti5cc.txt"  "seti5xxx.txt"
		"mt10cc.txt"  "mt10xxx.txt"  "setb4c9.txt"  "setb4xx.txt"  "setb4xyz.txt"  "seti5x.txt"  "seti5xy.txt"
		"mt10x.txt"  "mt10xy.txt"  "setb4cc.txt"  "setb4xxx.txt"  "seti5c12.txt"  "seti5xx.txt"  "seti5xyz.txt")

Dauzere=("01a.txt"  "04a.txt"  "07a.txt"  "10a.txt"  "13a.txt"  "16a.txt"
		 "02a.txt"  "05a.txt"  "08a.txt"  "11a.txt"  "14a.txt"  "17a.txt"
		 "03a.txt"  "06a.txt"  "09a.txt"  "12a.txt"  "15a.txt"  "18a.txt")


problemsMini=("YFJS_4_4_7_5_0" "YFJS_4_4_7_5_1" "YFJS_4_4_7_5_2" "YFJS_4_4_7_5_3" "YFJS_4_4_7_5_4" "YFJS_4_4_7_5_5" "YFJS_4_4_7_5_6"
			"YFJS_4_4_7_5_7" "YFJS_4_4_7_5_8" "YFJS_4_4_7_5_9" "YFJS_5_4_7_5_0" "YFJS_5_4_7_5_1" "YFJS_5_4_7_5_2" "YFJS_5_4_7_5_3" "YFJS_5_4_7_5_4"
			"YFJS_5_4_7_5_5" "YFJS_5_4_7_5_6" "YFJS_5_4_7_5_7" "YFJS_5_4_7_5_8" "YFJS_5_4_7_5_9" "YFJS_6_4_7_5_0" "YFJS_6_4_7_5_1" "YFJS_6_4_7_5_2"
			"YFJS_6_4_7_5_3" "YFJS_6_4_7_5_4" "YFJS_6_4_7_5_5" "YFJS_6_4_7_5_6" "YFJS_6_4_7_5_7" "YFJS_6_4_7_5_8" "YFJS_6_4_7_5_9"
   			"DA_2_5_1" "DA_2_5_10" "DA_2_5_11" "DA_2_5_12" "DA_2_5_13" "DA_2_5_14" "DA_2_5_15" "DA_2_5_16" "DA_2_5_17"
			"DA_2_5_18" "DA_2_5_19" "DA_2_5_2" "DA_2_5_20" "DA_2_5_3" "DA_2_5_4" "DA_2_5_5" "DA_2_5_6" "DA_2_5_7" 
			"DA_2_5_8" "DA_2_5_9" "DAFJS_3_5_1" "DAFJS_3_5_10" "DAFJS_3_5_2" "DAFJS_3_5_3" "DAFJS_3_5_4" "DAFJS_3_5_5" "DAFJS_3_5_6" "DAFJS_3_5_7" "DAFJS_3_5_8" "DAFJS_3_5_9" 
			)

let j=1
for rep in $(seq 1 1 5); do
	for alpha in $(seq .1 .1 .3); do
		for inst in ${Hurink[@]}; do
			#./FJS -i "hurink/edata/$inst" -o newclassic-TS-RN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 0 -itmax 0 -tsize 9 &
			#./FJS -i "hurink/edata/$inst" -o newclassic-TS-CN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 1 -itmax 0 -tsize 5 &
			#./FJS -i "hurink/rdata/$inst" -o newclassic-TS-RN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 0 -itmax 0 -tsize 9 &
			#./FJS -i "hurink/rdata/$inst" -o newclassic-TS-CN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 1 -itmax 0 -tsize 5 &
			#./FJS -i "hurink/vdata/$inst" -o newclassic-TS-RN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 0 -itmax 0 -tsize 9 &
			#./FJS -i "hurink/vdata/$inst" -o newclassic-TS-CN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 1 -itmax 0 -tsize 5 &
			#./FJS -i "hurink/edata/$inst" -o newclassic-SA.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh SA -tol 0 -c 0 -itmax -1 --pertMin 3 --pertMax 3 --T0m 0.78 --T0p 0.79 --Tf 0.001 --deltaMin 0.82 --deltaMax 0.82 -fType fix &
			#./FJS -i "hurink/rdata/$inst" -o newclassic-SA.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh SA -tol 0 -c 0 -itmax -1 --pertMin 3 --pertMax 3 --T0m 0.78 --T0p 0.79 --Tf 0.001 --deltaMin 0.82 --deltaMax 0.82 -fType fix &
			#./FJS -i "hurink/vdata/$inst" -o newclassic-SA.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh SA -tol 0 -c 0 -itmax -1 --pertMin 3 --pertMax 3 --T0m 0.78 --T0p 0.79 --Tf 0.001 --deltaMin 0.82 --deltaMax 0.82 -fType fix &
			#./FJS -i "hurink/vdata/$inst"  -o newclassic-grasp-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 1 -itmax 1 -alphaGRASP 0.59 &
			#./FJS -i "hurink/vdata/$inst"  -o newclassic-grasp-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 0 -itmax 1 -alphaGRASP 0.38 &
			#./FJS -i "hurink/vdata/$inst"  -o newclassic-ILS-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 0 -itmax 1 -ellmin 2 -ellmax 4 &
			#./FJS -i "hurink/vdata/$inst"  -o newclassic-ILS-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 1 -itmax 1 -ellmin 1 -ellmax 3 &
			#./FJS -i "hurink/edata/$inst"  -o newclassic-grasp-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 1 -itmax 1 -alphaGRASP 0.59 &
			#./FJS -i "hurink/edata/$inst"  -o newclassic-grasp-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 0 -itmax 1 -alphaGRASP 0.38 &
			#./FJS -i "hurink/edata/$inst"  -o newclassic-ILS-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 0 -itmax 1 -ellmin 2 -ellmax 4 &
			#./FJS -i "hurink/edata/$inst"  -o newclassic-ILS-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 1 -itmax 1 -ellmin 1 -ellmax 3 &
			#./FJS -i "hurink/rdata/$inst"  -o newclassic-grasp-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 1 -itmax 1 -alphaGRASP 0.59 &
			#./FJS -i "hurink/rdata/$inst"  -o newclassic-grasp-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 0 -itmax 1 -alphaGRASP 0.38 &
			#./FJS -i "hurink/rdata/$inst"  -o newclassic-ILS-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 0 -itmax 1 -ellmin 2 -ellmax 4 &
			#./FJS -i "hurink/rdata/$inst"  -o newclassic-ILS-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 1 -itmax 1 -ellmin 1 -ellmax 3 &
			./FJS -i "hurink/rdata/$inst"  -o newclassic-Tayebi.csv -t 300 -a -$alpha -s 0 -ls TAYEBI &
			./FJS -i "hurink/edata/$inst"  -o newclassic-Tayebi.csv -t 300 -a -$alpha -s 0 -ls TAYEBI &
			./FJS -i "hurink/vdata/$inst"  -o newclassic-Tayebi.csv -t 300 -a -$alpha -s 0 -ls TAYEBI &

			#./FJS -i $inst -o classic-ECT-definitivo.csv -t 300 -a -$alpha -s 0 -ls None -lse None -he ECT -mh None -tol 0 -c 0 &
			#./FJS -i $inst -o classic-SPT-definitivo.csv -t 300 -a -$alpha -s 0 -ls None -lse None -he SPT -mh None -tol 0 -c 0 &
			#./FJS -i $inst -o classic-LocalSearch-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh None -tol 0 -c 0 &
			#./FJS -i $inst -o classic-ModeloCP.csv -t 3600 -a -$alpha -s 0 -model MODELOCP &
			#./FJS -i "hurink/edata/$inst" -o classic2-ModeloMILP.csv -t 3600 -a -$alpha -s 0 -model MODELOMILP &
			#./FJS -i "hurink/rdata/$inst" -o classic2-ModeloMILP.csv -t 3600 -a -$alpha -s 0 -model MODELOMILP &
			#./FJS -i "hurink/vdata/$inst" -o classic2-ModeloMILP.csv -t 3600 -a -$alpha -s 0 -model MODELOMILP &
			#echo "dauzere/$inst"
			j=$(( $j + 1 ));
			if [ $( expr  $j % 8) == 0 ]
			then
				wait
			fi
		done
		for inst in ${Barnes[@]}; do
			#./FJS -i "barnes/$inst" -o newclassic-TS-RN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 0 -itmax 0 -tsize 9 &
			#./FJS -i "barnes/$inst" -o newclassic-TS-CN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 1 -itmax 0 -tsize 5 &
			#./FJS -i "barnes/$inst" -o newclassic-SA.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh SA -tol 0 -c 0 -itmax -1 --pertMin 3 --pertMax 3 --T0m 0.78 --T0p 0.79 --Tf 0.001 --deltaMin 0.82 --deltaMax 0.82 -fType fix &
			#./FJS -i "barnes/$inst"  -o newclassic-grasp-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 1 -itmax 1 -alphaGRASP 0.59 &
			#./FJS -i "barnes/$inst"  -o newclassic-grasp-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 0 -itmax 1 -alphaGRASP 0.38 &
			#./FJS -i "barnes/$inst"  -o newclassic-ILS-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 0 -itmax 1 -ellmin 2 -ellmax 4 &
			#./FJS -i "barnes/$inst"  -o newclassic-ILS-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 1 -itmax 1 -ellmin 1 -ellmax 3 &
			./FJS -i "barnes/$inst"  -o newclassic-Tayebi.csv -t 300 -a -$alpha -s 0 -ls TAYEBI &
			j=$(( $j + 1 ));
			if [ $( expr  $j % 24) == 0 ]
			then
				wait
			fi
		done
		for inst in ${Dauzere[@]}; do
			#./FJS -i "dauzere/$inst" -o newclassic-TS-RN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 0 -itmax 0 -tsize 9 &
			#./FJS -i "dauzere/$inst" -o newclassic-TS-CN.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 1 -itmax 0 -tsize 5 &
			#./FJS -i "dauzere/$inst" -o newclassic-SA.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh SA -tol 0 -c 0 -itmax -1 --pertMin 3 --pertMax 3 --T0m 0.78 --T0p 0.79 --Tf 0.001 --deltaMin 0.82 --deltaMax 0.82 -fType fix &
			#./FJS -i "dauzere/$inst"  -o newclassic-grasp-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 1 -itmax 1 -alphaGRASP 0.59 &
			#./FJS -i "dauzere/$inst"  -o newclassic-grasp-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 0 -itmax 1 -alphaGRASP 0.38 &
			#./FJS -i "dauzere/$inst"  -o newclassic-ILS-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 0 -itmax 1 -ellmin 2 -ellmax 4 &
			#./FJS -i "dauzere/$inst"  -o newclassic-ILS-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 1 -itmax 1 -ellmin 1 -ellmax 3 &
			./FJS -i "dauzere/$inst"  -o newclassic-Tayebi.csv -t 300 -a -$alpha -s 0 -ls TAYEBI &
			j=$(( $j + 1 ));
			if [ $( expr  $j % 24) == 0 ]
			then
				wait
			fi
		done
	done	
done

#let j=1
#for rep in $(seq 1 1 5); do
#	for alpha in $(seq .1 .1 .3); do
#		for inst in ${problemsclassic[@]}; do
#			./FJS -i $inst  -o classic-grasp-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 1 -itmax 1 -alphaGRASP 0.59 &
#			./FJS -i $inst  -o classic-grasp-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 0 -itmax 1 -alphaGRASP 0.38 &
#			./FJS -i $inst  -o classic-SA-definitivo.csv -t 300 -a -$alpha -s 0 -ls None -lse None -he Best -mh SA -tol 0 -c 0 -itmax -1 --pertMin 3 --pertMax 3 --T0m 0.78 --T0p 0.79 --Tf 0.001 --deltaMin 0.82 --deltaMax 0.82 -fType fix &
#			./FJS -i $inst  -o classic-ILS-RN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 0 -itmax 1 -ellmin 2 -ellmax 4 &
#			./FJS -i $inst  -o classic-ILS-CN-definitivo.csv -t 300 -a -$alpha -s 0 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 1 -itmax 1 -ellmin 1 -ellmax 3 &
#			./FJS -i $inst  -o classic-Tayebi.csv -t 300 -a -$alpha -s 0 -mh TAYEBI &
#			j=$(( $j + 1 ));
#			if [ $( expr  $j % 4) == 0 ]
#			then
#				wait
#			fi
#		done
#	done	
#done

