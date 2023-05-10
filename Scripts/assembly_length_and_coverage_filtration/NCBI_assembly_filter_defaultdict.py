
# >NODE_1_length_712560_cov_106.220214
# 3, 5


from collections import defaultdict
import sys

if len(sys.argv)!=3: 
    print("Usage:python3 NCBI_assembly_filter_defaultdict.py <input.txt> <output.txt> ")
    sys.exit()

handle = open(sys.argv[1],'r')
lines = handle.read().split('\n')

my_dict= defaultdict(list)
for line in lines:
    if line.startswith(">"):
        header = line
    else:
        my_dict[header].append(line)
handle.close()

outfile=open(sys.argv[2],'w')
for line1 in my_dict.keys():
    if line1.startswith(">"):
        (length,random, cov) = line1.split('-')[0].split('_')[3:]
        if float(length) >= 200 and float(cov) >= 5:
            outfile.write(line1+"\n"+''.join(my_dict[line1])+'\n')
outfile.close()


