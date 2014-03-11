__author__ = 'saltikov'


def plate_header():
    """Generates the A1-H12 header"""
    header = [ ]
    beginNum = ord('A')
    endNum = ord('H')
    for number in xrange(beginNum, endNum+1):
        for i in xrange(1,13):
            header.append(chr(number)+i.__str__())
    return header


basedir = '/Users/saltikov/Documents/platereader_script/'

big_line = [ ]

outfile = basedir+'output.txt'
infile = basedir+"DeltaRNAseq031014.txt"

fileout = open(outfile, "w")


head = plate_header()
for item in head:
    fileout.write(item+"\t")
fileout.write("\n")

i = 1
for line in open(infile, "r"):
    line = line.strip()
    line = line.split("\t")
    if len(line) == 13:
        big_line.extend(line[1:13])
        i +=1
    if len(line) == 12:
        big_line.extend(line[0:12])
        i +=1
    if i == 9:
        for item in big_line:
            fileout.write(item+'\t')
        fileout.write('\n')
        big_line = [ ]
        i = 1
