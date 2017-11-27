import numpy as np
import os
#By Xiaoquan Sun
#This script is used to find structure water and corresponding amino acid on receptor and ligand.
#There are following arrays to store the data
ligand_wat=[];receptor_wat=[] #The water residue number of ligand or receptor
ligand_wat_uniq=[] # Water residue number which form hbond with ligand, remove repeated terms
receptor_wat_uniq=[] # Water residue number which form hbond with receptor, remove repeated terms
wat_total_uniq=[] # ligand_wat_uniq + receptor_wat_uniq and remove the repeated terms
wat_uniq_count=[] # times of elements in wat_total_uniq
stru_wat=[] # element in wat_total_uniq has wat_uniq_count >=2
frame_vs_stru_wat=[] # the residue number of structure water at each frame
frame_vs_stru_wat_num=[] # store the number of structure water at each frame, Format "frame number_of_structure_wat"
ligand_aa=[] # amino acid form hbond with structure water
receptor_aa=[] # amino acid form hbond with structure water
frame_stru_water_aa_ligand=[] # store the amino acid pair form hbond with the same structure water, Format "frame water# amino_acid#"
frame_stru_water_aa_receptor=[]
stru_wat_ligand_aa=[] # all amino acid in ligand form hbond with structure water at each frame
stru_wat_receptor_aa=[] # all amino acid in receptor form hbond with structure water at each frame
ligand_aa_total=[];ligand_aa_count=[] #  count the amino aicd times for each 1000 frame
receptor_aa_total=[]; receptor_aa_count=[] # count the amino acid times for each 1000 frame


for i in range(2,3006):
	#set the file name of each loop
	la='avg_'+str(i)+'la.dat'
	lb='avg_'+str(i)+'lb.dat'
	ra='avg_'+str(i)+'ra.dat'
	rb='avg_'+str(i)+'rb.dat'

	##ligand part
	ligand_wat=[]; ligand_aa=[]
	with open(la) as file:
		la=[line.split() for line in file]
	for j in range(1,len(la)): # the first line is title, so start from 1 rather than 0
		ligand_wat.append(la[j][0].split("@")[0].split("_")[1])
		ligand_aa.append(la[j][2].split("@")[0])#.split("_")[1])

	with open(lb) as file:
		lb=[line.split() for line in file]
	for j in range(1,len(lb)): # the first line is title, so start from 1 rather than 0
		ligand_wat.append(lb[j][2].split("@")[0].split("_")[1])
		ligand_aa.append(lb[j][0].split("@")[0])#.split("_")[1])	

	ligand_wat_uniq=np.unique(ligand_wat)
	
	###################	

	##receptor part
	receptor_wat=[]; receptor_aa=[]
	with open(ra) as file:
		ra=[line.split() for line in file]
	for j in range(1,len(ra)): # the first line is title, so start from 1 rather than 0
		receptor_wat.append(ra[j][0].split("@")[0].split("_")[1])
		receptor_aa.append(ra[j][2].split("@")[0])#.split("_")[1])

	with open(rb) as file:
		rb=[line.split() for line in file]
	for j in range(1,len(rb)): # the first line is title, so start from 1 rather than 0
		receptor_wat.append(rb[j][2].split("@")[0].split("_")[1])
		receptor_aa.append(rb[j][0].split("@")[0])#.split("_")[1])	

	receptor_wat_uniq=np.unique(receptor_wat)
	############################################

	###################	find structure water number at each frame

	wat_total_uniq, wat_uniq_count=np.unique(np.append(ligand_wat_uniq, receptor_wat_uniq), return_counts=True)
	stru_wat=[]
	for j in range(len(wat_total_uniq)):
		if wat_uniq_count[j] > 1:
			stru_wat.append(wat_total_uniq[j])

	if not stru_wat:
		stru_wat.append('0') #no structure water found add 0
		frame_vs_stru_wat.append([[i], [0]])
		frame_vs_stru_wat_num.append([[i], [0]])
	else:
		frame_vs_stru_wat.append([[i],stru_wat])
		frame_vs_stru_wat_num.append([[i],[len(stru_wat)]])
	 
	##############################################
	#########find amino acid form hbond with structure in both ligand and receptor
	temp0=[]
	for j in range(len(stru_wat)):
		temp=[]
		for k in range(len(ligand_wat)):
			if stru_wat[j] == ligand_wat[k]:
				temp0.append(ligand_aa[k])
				temp.append(ligand_aa[k])
		if not temp:
			temp0.append([0])
			temp.append([0])
		frame_stru_water_aa_ligand.append([[i],[stru_wat[j]],temp])
	stru_wat_ligand_aa.append([[i],temp0])
	
	temp0=[]
	for j in range(len(stru_wat)):
		temp=[]
		for k in range(len(receptor_wat)):
			if stru_wat[j] == receptor_wat[k]:
				temp0.append(receptor_aa[k])
				temp.append(receptor_aa[k])
		if not temp:
			temp0.append('0')
			temp.append('0')
		frame_stru_water_aa_receptor.append([[i],[stru_wat[j]],temp])
	stru_wat_receptor_aa.append([[i],temp0])	
	############################################		

np.savetxt('1.dat', frame_vs_stru_wat, fmt='%s')
np.savetxt('2.dat', frame_vs_stru_wat_num, fmt='%s')
np.savetxt('3.dat', frame_stru_water_aa_ligand, fmt='%s')
np.savetxt('4.dat', frame_stru_water_aa_receptor, fmt='%s')
np.savetxt('5.dat', stru_wat_ligand_aa, fmt='%s')
np.savetxt('6.dat', stru_wat_receptor_aa, fmt='%s')
os.system('for i in {1..6};do source ~/common/structure_water/remove_symbol.sh ${i}.dat; done') #remove all symbols
 ####################

#the hbond number of each 10ns simulation for ligand and receptor
with open('5.dat') as file:
	list5=[line.split() for line in file]
with open('6.dat') as file:
 	list6=[line.split() for line in file]

cutoff=0 #hbond bonds in 1000 frames above this cutoff will output
temp=[]
for i in range(len(list5)):
 	temp.append(list5[i][1:])
for i in range(len(temp)/1000):
 	m=i*1000; n=m+1000
 	file='ligand'+str(n/100)+'ns.dat'
 	temp0=temp[m:n]
 	temp1=[]
 	for j in range(len(temp0)):
 		temp1.extend(temp0[j])
 	aa, aa_counts=np.unique(temp1, return_counts=True)
 	temp0=[]
 	for j in range(len(aa)):
 		if aa[j] != '0' and aa_counts[j] > cutoff:
 			temp0.append([aa[j],aa_counts[j]])
 	np.savetxt(file, temp0, fmt='%s')

temp=[]
for i in range(len(list6)):
 	temp.append(list6[i][1:])
for i in range(len(temp)/1000):
 	m=i*1000; n=m+1000
 	file='receptor'+str(n/100)+'ns.dat'
 	temp0=temp[m:n]
 	temp1=[]
 	for j in range(len(temp0)):
 		temp1.extend(temp0[j])
 	aa, aa_counts=np.unique(temp1, return_counts=True)
 	temp0=[]
 	for j in range(len(aa)):
 		if aa[j] != '0' and aa_counts[j] > cutoff:
 			temp0.append([aa[j],aa_counts[j]])
 	np.savetxt(file, temp0, fmt='%s')
 ############################################################


######the following content is using output the residue pair for the same structure water
with open('3.dat') as file:
	list3=[line.split() for line in file]
with open('4.dat') as file:
	list4=[line.split() for line in file]

temp0=[]
for i in range(len(list3)):
	m=list3[i][:2]
	for j in range(len(list3[i][2:])):
		for k in range(len(list4[i][2:])):
			temp0.append([m,list3[i][2+j],list4[i][2+k]])
np.savetxt('7.dat', temp0, fmt='%s')
os.system('for i in 7;do source ~/common/structure_water/remove_symbol.sh ${i}.dat; done') 

#####################################################################################
with open('7.dat') as file:
	list7=np.array([line.split() for line in file])

temp0=[]
cutoff2=0 #the cutoff for the pair
x=np.array(list7[:,0],dtype=int)
for i in range(len(list5)/1000):
	m=i*1000; n=m+1000
	file='pair'+str(n/100)+'ns.dat'
	o=np.where(np.logical_and(x>=m,x<n))
	temp=np.array(list7[o][:,2:])
#	ncols = temp.shape[1]
#	dtype = temp.dtype.descr * ncols
#	struct = temp.view(dtype)
#	uniq, uniq_counts = np.unique(struct, return_counts=True)
#	uniq = uniq.view(temp.dtype).reshape(-1, ncols)
#	for j in range(len(uniq)):
#		if uniq[j][0] != '0' and uniq_counts[j] > cutoff2:
#			temp0.append([uniq[j],uniq_counts[j]])
	np.savetxt(file, temp, fmt='%s')

os.system('for i in `ls -lrtx1 pair*ns.dat`;do source ~/common/structure_water/remove_symbol.sh ${i}; done') 
os.system('for i in `ls -lrtx1 pair*ns.dat`;do source ~/common/structure_water/uniq_pair.sh ${i}; done') 

# #the following script is used to convert the format for gnuplot
os.system('source ~/common/structure_water/total2.sh') 

# #the following script is using gnuplot to generate the figure
os.system('source ~/common/structure_water/rm_small_occupancy.sh ligand_total.dat; mv big_occupancy.dat ligand_big_occupancy.dat') 
os.system('source ~/common/structure_water/rm_small_occupancy.sh receptor_total.dat; mv big_occupancy.dat receptor_big_occupancy.dat')
os.system('source ~/common/structure_water/rm_small_occupancy.sh uniq_pair_total.dat; mv big_occupancy.dat uniq_pair_big_occupancy.dat')
os.system('cp *_total.pdf ../')


