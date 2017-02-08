import HTSeq
import sys, getopt


def frame_id(seq,key):
    '''
    Usage: frame_id(dictionary['key'])
    frame_id(sequences['>Seq1'])
    six_frames = frame_id(sequences['>Seq1'])
'''
    frames = {'+1':[],'+2':[],'+3':[],'-1':[],'-2':[],'-3':[]}
    seq_rev = rev_seq(seq)
    for j in range(0,3):
        temp = ''.join([seq[j::]])
        temp_rev = ''.join([seq_rev[j::]])
        seq_trans = swap_dna(temp)
        #print key + '--' + str(j)
        #print seq_trans
        seq_rev_trans = swap_dna(temp_rev)
        #print key + '--' + str(j) + '_trans'
        #print seq_rev_trans
        if j==0:
            frames['+1']=seq_trans
            frames['-1']=seq_rev_trans
        if j==1:
            frames['+2']=seq_trans
            frames['-2']=seq_rev_trans
        if j==2:
            frames['+3']=seq_trans
            frames['-3']=seq_rev_trans

    return frames

def rev_seq(seq):
    trans=[]
    for i in seq:
        if i=='A':
            trans.append('T')
        elif i=='C':
            trans.append('G')
        elif i=='G':
            trans.append('C')
        elif i=='T':
            trans.append('A')
        else:
            trans.append(i)
    trans=''.join(trans)
    seq_rev= trans[::-1]
    return seq_rev



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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
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
      opts, args = getopt.getopt(argv,"hi:o:p",["ifile=","ofile="])
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
    all_Frames = frame_id(DNA,name)
    for k,v in all_Frames.iteritems():
        outputting.write('>' + name + '--' + k + '\n')
        outputting.write(v + '\n')
outputting.close()
