import HTSeq
import sys, getopt








def swap_dna(dnastring):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
        }
    protein = []
    end = len(dnastring) - (len(dnastring) %3) - 1
    for i in range(0,end,3):
        codon = dnastring[i:i+3]
        if codon in table:
            aminoacid = table[codon]
            protein.append(aminoacid)
        else:
            protein.append("N")
    return "".join(protein)





def main(argv):
    inputfile = ""
    outputfile = ""
    printprots = False

    try:
      opts, args = getopt.getopt(argv,"hi:o:k:p",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'dna2prot4u.py -i <inputfile> REDIRECT TO OWN FILE -p'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'dna2prot4u.py -i <inputfile> -o <outputfile> -p'
            print "-p prints the protein sequences in the terminal"
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt == "-p":
            printprots = True

    return inputfile, outputfile, printprots



if __name__ == "__main__":
   inputfile, outputfile, printprots = main(sys.argv[1:])



outputting =  open(outputfile, "w")  
for read in HTSeq.FastaReader(inputfile):
    name = read.name
    DNA = read.seq
    Prot = swap_dna(DNA)
    outputting.write('>' + name + '\n')
    outputting.write(Prot + '\n')
outputting.close()
