# This program takes raw allele counts at SNP sites and assesses freqency 
#from each parent. In this case the Ref is B73 (spider mite sensitive) and 
# the Alt is from the spider mite resistant haplotype on chromosome 6. The
# input files were processed using the pysam module from the combined .bam
# files for sensitive and resistant bulks in crosses of B73 (sensitive) to
# B49, B75, and B96. The two input files are:
# Sen_Chr6_Allele_Counts.txt
# Res_Chr6_Allele_Counts.txt

# First, process combined sensitive bulk data

file = open('Sen_Chr6_Allele_Counts.txt','r')

outfile = open('Combined_Sen_Chr6_freq.txt','w')

outfile.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('begin', 'end', 'snp_count', 'ref_count', 'alt_count', 'total', 'freq_resistant') )

ranges = []

start_pos = 7895601
end_pos = 109684471
window_size = 5000000
step = 1000000

while (start_pos + window_size) <= end_pos:
    ranges.append([start_pos,(start_pos+window_size),0,0,0])
    start_pos += step

print(ranges)

for line in file:
    if line[0] != 'c':
        line=line.strip()
        line=line.split('\t')
        pos = int(line[1])
        ref = int(line[4])
        alt = int(line[5])
        for j in range(len(ranges)):
            if ranges[j][0] <= pos <= ranges[j][1]:
                ranges[j][2] = ranges[j][2] + 1
                ranges[j][3] = ranges[j][3] + ref
                ranges[j][4] = ranges[j][4] + alt

for k in range(len(ranges)):
    freq_R = ranges[k][4] / (ranges[k][3] + ranges[k][4])
    outfile.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\n'  %  (ranges[k][0], ranges[k][1], ranges[k][2], ranges[k][3], ranges[k][4], (ranges[k][3] + ranges[k][4]), round(freq_R,3))  )

file.close()
outfile.close()


# Second, process combined resistant bulk data

file = open('Res_Chr6_Allele_Counts.txt','r')

outfile = open('Combined_Res_Chr6_freq.txt','w')

outfile.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('begin', 'end', 'snp_count', 'ref_count', 'alt_count', 'total', 'freq_resistant') )

ranges = []

start_pos = 7895601
end_pos = 109684471
window_size = 5000000
step = 1000000

while (start_pos + window_size) <= end_pos:
    ranges.append([start_pos,(start_pos+window_size),0,0,0])
    start_pos += step

print(ranges)

for line in file:
    if line[0] != 'c':
        line=line.strip()
        line=line.split('\t')
        pos = int(line[1])
        ref = int(line[4])
        alt = int(line[5])
        for j in range(len(ranges)):
            if ranges[j][0] <= pos <= ranges[j][1]:
                ranges[j][2] = ranges[j][2] + 1
                ranges[j][3] = ranges[j][3] + ref
                ranges[j][4] = ranges[j][4] + alt

for k in range(len(ranges)):
    freq_R = ranges[k][4] / (ranges[k][3] + ranges[k][4])
    outfile.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\n'  %  (ranges[k][0], ranges[k][1], ranges[k][2], ranges[k][3], ranges[k][4], (ranges[k][3] + ranges[k][4]), round(freq_R,3))  )

file.close()
outfile.close()

# Finally, open Sen and Res output files and process to get snp counts, window
# midpoint, and ratio of sensitive to resistant (the input for plotting Fig. S5).

Sen = open('Combined_Sen_Chr6_freq.txt','r')
Res = open('Combined_Res_Chr6_freq.txt','r')
outfile = open('Input_for_ploting_5_Mb.txt','w')

order = []
processed = {} # [win_begin, wind_end] as KEY
# As VALUE create [midpoint, snps_sen, ratio_sen, snps_res, ratio_res]

for i in Sen:
    i=i.strip()
    i=i.split('\t')
    if i[0] != 'begin':
        print(i)
        order.append((int(i[0]),int(i[1])))
        midpoint = ( int(i[0]) + int(i[1])  ) / 2
        midpoint = round(midpoint,0)
        processed[int(i[0]),int(i[1])] = [midpoint,int(i[2]),float(i[-1]),0,0]
        print(processed)

for i in Res:
    i=i.strip()
    i=i.split('\t')
    if i[0] != 'begin':
        print(i)
        if (int(i[0]),int(i[1])) not in order:
            raise
        processed[int(i[0]),int(i[1])][3] = int(i[2])
        processed[int(i[0]),int(i[1])][4] = float(i[-1])
        print(processed)

outfile.write('%s\t%s\t%s\t%s\t%s\n' % ('Midpoint','SNP_count','Res_Freq_Sen_Pool','Res_Freq_Res_Pool','Freq_Difference' ))

for i in order:
    i = processed[i]
    diff = round( (i[4]-i[2]),3)
    outfile.write(   '%s\t%s\t%s\t%s\t%s\n' % (i[0],i[1],i[2],i[4],diff)   )

Sen.close()
Res.close()
outfile.close()