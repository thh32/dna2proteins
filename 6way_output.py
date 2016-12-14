### Import all functions

## Function to read fasta files

def read_fasta(fastafile):
    """
    Reads a fasta file and returns a dictionary with sequence
    number as keys and sequence code as values
    """
    sequences = []
    with open(fastafile, "r") as f:
        ls = f.readlines()
        for i in ls:
             sequences.append(i.rstrip("\n"))

    seq_id = []
    for i in sequences:
        if i.startswith(">"):
            seq_id.append(i)

    seq_id_index = []
    for i in range(len(seq_id)):
        seq_id_index.append(sequences.index(seq_id[i]))

    seq_dic = {}
    for i in range(len(seq_id_index)):
        if i == (len(seq_id_index) - 1):
            seq_dic[seq_id[i]] = sequences[seq_id_index[i]+1:]
        else:
            seq_dic[seq_id[i]] = sequences[seq_id_index[i]+1:seq_id_index[i+1]]

    seq_dic_2 = {}
    for keys, values in seq_dic.items():
        seq_dic_2[keys] = "".join(values)

    return seq_dic_2

## Writes a dictionary to a fasta file

def write_fasta(dictionary, filename):
    """
    Takes a dictionary and writes it to a fasta file
    Must specify the filename when caling the function
    """

    import textwrap
    with open(filename, "w") as outfile:
        for key, value in dictionary.items():
            outfile.write(key + "\n")
            outfile.write("\n".join(textwrap.wrap(value, 60)))
            outfile.write("\n")

    print "Success! File written"

## Swaps DNA sequencs for proteins

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

## Generates the six possible frames per one sequence

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
        seq_rev_trans = swap_dna(temp_rev)

        #OUTPUT
        minimum_int = int(minimum)
        if j == 0:
            if len(seq_trans.replace('_','')) >= minimum_int:
                f1.write(key + '\n')
                f1.write(seq_trans + '\n')
            if len(seq_rev_trans.replace('_','')) >= minimum_int:
                f2.write(key + '\n')
                f2.write(seq_rev_trans + '\n')
        elif j == 1:
            if len(seq_trans.replace('_','')) >= minimum_int:
                f3.write(key + '\n')
                f3.write(seq_trans + '\n')
            if len(seq_rev_trans.replace('_','')) >= minimum_int:
                f4.write(key + '\n')
                f4.write(seq_rev_trans + '\n')
        elif j == 2:
            if len(seq_trans.replace('_','')) >= minimum_int:
                f5.write(key + '\n')
                f5.write(seq_trans + '\n')
            if len(seq_rev_trans.replace('_','')) >= minimum_int:
                f6.write(key + '\n')
                f6.write(seq_rev_trans + '\n')        
        #print key + '--' + str(j)
        #print seq_trans
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

## Required function for the previos frame_id function

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

## Generates all the frames for all the sequences

def gen_frames(dictionary):
    all_dict = {}
    for key, value in dictionary.items():
        all_dict[key] = frame_id(dictionary[key],key)

    return all_dict

## Find the open frames in the protein sequences

def oframe(amino):
    oframes = []
    for i in range(0,len(amino)):
        if amino[i]=='M':
            temp = ''.join([amino[i::]])
            oframe=temp[0:temp.find('_')+1]
            oframes.append(oframe)
    return oframes

## Finds the longest proteins in each sequence

def find_prots(dictionary):
    prots_dict = {}
    for key, value in dictionary.items():
        poss_protein = []
        for f in value:
            poss_protein += (oframe(value[f]))
            #print key, poss_protein
            c = 0
            result = ""
            for s in poss_protein:
                if len(s) > c:
                    result = s
                    c = len(s)
                else:
                    continue
            prots_dict[key] = result

    return prots_dict

## For command line Usage

import sys, getopt

def main(argv):
    inputfile = ""
    outputfile = ""
    minimum = 5
    printprots = False
    print 'Running...'
    try:
      opts, args = getopt.getopt(argv,"hi:o:m:p",["ifile=","ofile=","min="])
    except getopt.GetoptError:
        print 'dna2proteins.py -i <inputfile> -o <output pattern> -m <minimum protein length> -p'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'dna2proteins.py -i <inputfile> -o <outputfile> -p'
            print "-p prints the protein sequences in the terminal"
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-m", "--min"):
            outputfile = arg
        elif opt == "-p":
            printprots = True

    return inputfile, outputfile,minimum, printprots

if __name__ == "__main__":
   inputfile, outputfile,minimum, printprots = main(sys.argv[1:])

print 'Reading in file'
sequences = read_fasta(inputfile)

f1 = open(outputfile + '.f1.txt','w')
f2 = open(outputfile + '.f2.txt','w')
f3 = open(outputfile + '.f3.txt','w')
f4 = open(outputfile + '.f4.txt','w')
f5 = open(outputfile + '.f5.txt','w')
f6 = open(outputfile + '.f6.txt','w')


print 'Producing protein sequeces'
sequences_frames = gen_frames(sequences)

f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()

print 'Closing files'
proteins = find_prots(sequences_frames)

if printprots == True:
    for key, values in proteins.items():
        print key
        print values

print 'Finished.'

