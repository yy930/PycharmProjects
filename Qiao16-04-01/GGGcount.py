import re
import sys

# dictionary for barcode,umi
rowbarcodes = {'AGT': 'A', 'ACC': 'B', 'CTA': 'C', 'CAG': 'D', 'GCG': 'E', 'TTT': 'F', 'TAC': 'G', 'GGA': 'H'}

def ReverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])

colbarcodes = {ReverseComplement('CGTGAT'): '1', ReverseComplement('ACATCG'): '2', ReverseComplement('GCCTAA'): '3',
               ReverseComplement('TGGTCA'): '4', ReverseComplement('CACTGT'): '5', ReverseComplement('ATTGGC'): '6',
               ReverseComplement('GATCTG'): '7', ReverseComplement('TCAAGT'): '8', ReverseComplement('CTGATC'): '9',
               ReverseComplement('AAGCTA'): '10', ReverseComplement('GTAGCC'): '11', ReverseComplement('TACAAG'): '12'}

#transposon
TN5 = 'CTGTCTCTTATACACATCTGACGC'
TN5R = 'GCGTCAGATGTGTATAAGAGACAG'
P5 = 'AATGATACGGCGACCACCGAT'
P5R = 'ATCGGTGGTCGCCGTATCATT'
P7 = 'CAAGCAGAAGACGGCATACGA'
P7R = 'CGTATGCCGTCTTCTGCTTG'

noBarcodeCount = 0
noGGCount = 0
passCount = 0
# output=''
# undef=''
# def main(inputFileName, outputFileName, undefFileName):
#     inputFileName = sys.argv[1]
#     outputFileName = sys.argv[2]
#     undefFileName = sys.argv[3]
#     global output
#     global undef
#     output = outputFileName
#     undef = undefFileName
#     #from collections import defaultdict
#     parseFile(inputFileName)
#     #parseFile("/Users/yingy_adm/data/1-1311-1_ATCACG_L001_R1_001.fastq")

class Read:
    def __init__(self, lines):
        self.id = lines[0]
        self.sequence = lines[1]
        self.quality = lines[3]
        self._extractbarcodeumi()

    def _extractbarcodeumi(self):
        global noBarcodeCount
        global noGGCount
        self.barcode = ''
        self.umi = ''
        self.GG = 0

        for rbarc in rowbarcodes:
            if rbarc == self.sequence[0:3]: # barcode is matched
                self.barcode = rowbarcodes[rbarc]
                self.umi = self.sequence[3:8]

                if self.sequence[8:10] == 'GG': #GG matched in the right position
                    self.GG = 1
                    #remove Gs, barcode, umi
                    Gposition = re.search('G{2,}',self.sequence[8:])
                    self.sequence = self.sequence[Gposition.end()+8:]
                    self.quality = self.quality[Gposition.end()+8:]

                    # if Transposon is found, remove all the following
                    if TN5 in self.sequence:
                        TN5pos = self.sequence.find(TN5)
                        # transposon=line[TN5pos:TN5pos+len(TN5)]
                        self.sequence = self.sequence[:TN5pos]
                        self.quality = self.quality[:TN5pos]
                    elif TN5R in self.sequence:
                        TN5Rpos = self.sequence.find(TN5R)
                        # self.transposon=line[TN5pos:len(TN5R)]
                        self.sequence = self.sequence[:TN5Rpos]
                        self.quality = self.quality[:TN5Rpos]
                    # identify barcode, and remove all the following
                    for cbarc in colbarcodes:
                        cbarcpos = self.sequence.find(cbarc)
                        if cbarcpos != -1:
                            self.sequence = self.sequence[:cbarcpos]
                            self.quality = self.quality[:cbarcpos]
                            break
                break

    def writeNF(self, nf):
        nf.write(self.id)  # write id, barcode'
        nf.write(self.sequence.strip('\n'))  # remove all the \n in sequence line
        nf.write('\n')
        nf.write('+\n')
        nf.write(self.quality.strip('\n'))
        nf.write('\n')

    def writeF(self, f):
        global passCount
        header = '|'.join([self.id.strip('\n'), 'parsedreads', self.barcode, self.umi])# new id
        f.write(header)  # write id, barcode'
        f.write('\n')
        f.write(self.sequence.strip('\n'))  # remove all the \n in sequence line
        f.write('\n')
        f.write('+\n')
        f.write(self.quality.strip('\n'))
        f.write('\n')

#def parseFile(fn, output, undef):
reads = []  # list of all reads in a file 'fn'
collectedLines = []  # 4 lines each time
fn = sys.argv[1]
output = sys.argv[2]
undef = sys.argv[3]
for line in open(fn):
    if line[0] == '@':
        if len(collectedLines) == 4:
            read = Read(collectedLines)
            reads.append(read)
            collectedLines = []
    collectedLines.append(line)

#readsPerBarcode = defaultdict(list)
f = open(output, 'w')
nf = open(undef,'w')
for read in reads:
    if read.barcode == '':
        noBarcodeCount += 1
        read.writeNF(nf)
    elif read.GG == 0:
        noGGCount += 1
        read.writeNF(nf)
    else:
        #readsPerBarcode[read.barcode].append(read) # write into dict
        passCount += 1
        read.writeF(f)
# here you can write directly to a file,
f.close()
nf.close()
al= noBarcodeCount+noGGCount+passCount
print al, "reads in all."
print noBarcodeCount,  float(noBarcodeCount)*100/al, "% barcodes were not found."
print noGGCount,  float(noGGCount)*100/al,"% GG were not found"
print passCount,  float(passCount)*100/al,"% reads passed."

#if __name__ == "__main__":
    #from collections import defaultdict
    #f = open('file1.fastq', 'a')
    #parseFile("./data/1-1311-1_ATCACG_L001_R1_001.fastq")
    #parseFile("/Users/yingy_adm/data/1-1311-1_ATCACG_L001_R1_001.fastq")
    #f.close()
# parseFile('/Users/yingy_adm/data/2-1311-2_CGATGT_L001_R1_001.fastq','/Users/yingy_adm/data/fle2', '/Users/yingy_adm/data/undefile2')