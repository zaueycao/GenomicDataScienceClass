import wget
#print out the sequence and count ATCG in the sequences
site_url = 'http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa'
file_name = wget.download(site_url)
print(file_name)
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
genome = readGenome('lambda_virus.fa')
print(genome[:100])
#count the number of occurences in the base
count = {'A':0, 'C':0, 'G':0,'T':0 }
for base in genome:
    count[base] += 1
print(count)

#this part is to read the Fastq dataset
site_url_2 =  'http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/SRR835775_1.first1000.fastq'
file_name_2 = wget.download(site_url_2)
print(file_name_2)
def readfastq(file_name_2):
    sequence = []
    quality = []
    with open(file_name_2,'r') as fh:
        while True: 
            fh.readline() # skip name line
            seq = fh.readline().rstrip() #read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() # read quality line
            if len(seq) == 0:
                break
            sequence.append(seq)
            quality.append(qual)
        return sequence, quality 
seq, qual = readfastq('SRR835775_1.first1000.fastq')

def phred33toQ(qual):
    return ord(qual)-33
#this finction is to have a histogram for quality of the base
def histQuality(qualities):
    hist = [0]*50
    for read in qualities:
        for phred in read:
            q = phred33toQ(phred)
            hist[q] += 1
    return hist 
h = histQuality(qual)
print(h)

#plot the histogram

import matplotlib.pyplot as plt
plt.plot(range(len(h)), h)
plt.show()