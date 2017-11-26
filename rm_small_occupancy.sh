a=(`cat $1 | awk '{print $1}'`)
for ((i=0;i<${#a[*]};i++));do
	unset b
	((j=$i+1))
	b=(`sed -n ${j}p $1 | cut -f2-`)
	for ((p=0;p<${#b[*]};p++));do
		((c=`echo "${b[$p]} * 10 " | bc | awk '{printf "%f", $0}' | cut -f1 -d"."`))
		if (( $c >= 6 ));then
			d[$i]=${a[$i]}
			break
		fi
	done
done

d=(`echo ${d[*]}`)

rm -f big_occupancy.dat
for ((i=0;i<${#d[*]};i++));do
	grep ${d[$i]} $1 | tr '\t ' '\012' > ${i}_temp.dat
done

i=`ls *_temp.dat | sort -n`
paste $i > big_occupancy_temp.dat
((p=0))
echo Time > time_gap1_temp.dat
((n=`cat big_occupancy_temp.dat | wc -l`))
for i in {0..200..10};do
	((j=$i+10))
	echo ${i}\-${j}ns >> time_gap1_temp.dat
	((i=$i+10))
done
sed -n 1,${n}p time_gap1_temp.dat > time_gap_temp.dat
paste time_gap_temp.dat big_occupancy_temp.dat > big_occupancy.dat
rm *_temp.dat

cat > hbond.gnu << \EOF
reset
set terminal pdfcairo
set boxwidth 0.9 absolute
set style fill solid 1.00 border lt -1
set key outside right top vertical Right noreverse noenhanced autotitle nobox
set style histogram clustered gap 5 title textcolor lt -1
set datafile missing '-'
set style data histograms
set xtics border in scale 0,0 nomirror rotate by -45  autojustify
set xtics ()
set yl "Average Hydrogen Bonds Number (/ps)"
set title "Hydrogen bonds through structure water"
EOF
name=`echo $1 | cut -f1 -d"."`
((n=${#d[*]}+1))

cat >> hbond.gnu << EOF
set output "${name}.pdf"
p for [i=2:${n}] 'big_occupancy.dat' u i:xticlabels(1) ti col
set output
EOF

plot hbond.gnu
