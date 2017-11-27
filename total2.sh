dir=`pwd`
cd $dir
rm -f list
for k in ligand receptor uniq_pair;do
	aa=()
	for i in `ls -lrtx1 ${k}*ns.dat`;do
		aa+=(`cat $i | awk '{print $1}'`)
		echo $i >> list
	done
	aa=(`echo ${aa[*]} | tr ' ' '\n' | sort | uniq`)
	echo ${aa[*]} | tr ' ' '\n' > temp_0.dat
	
	((p=1))
	for j in `ls -lrtx1 ${k}*ns.dat`;do
	count=()
	echo $j >> list
		for i in ${aa[*]};do
			aa_count=`grep -w $i $j | awk '{print $2}'`
			if [ -z "$aa_count" ]; then aa_count=0; fi
			aa_count=`echo "$aa_count / 1000" | bc -l | awk '{printf "%f", $0}' | cut -c -4`
			count+=($aa_count)
		done
	echo ${count[*]} | tr ' ' '\n' > temp_${p}.dat
	((p++))
	done


	i=`ls -lrtx1 temp_*.dat`;echo $i >> list
	paste $i > ${k}_total.dat
done
rm temp_*.dat

 
