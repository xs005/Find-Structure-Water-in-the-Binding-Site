sed -i s/\ /+/g $1
cat $1 | sort | uniq -c > temp_uniq
sed -i /0+0/d temp_uniq
cat temp_uniq | awk '{print $2 " " $1}' > uniq_$1
rm temp_uniq
 
