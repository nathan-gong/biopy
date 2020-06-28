import itertools

# DATA TABLES

codon = {
    'UUU' : 'F', # Phenylalanine (Phe)
    'UUC' : 'F', 
    'UUA' : 'L', # Leucine (Leu)
    'UUG' : 'L',
    'UCU' : 'S', # Serine (Ser)
    'UCC' : 'S',
    'UCA' : 'S',
    'UCG' : 'S',
    'UAU' : 'Y', # Tyrosine (Tyr)
    'UAC' : 'Y',
    'UAA' : 'Stop', # Stop Codon (Stop)
    'UAG' : 'Stop',
    'UGU' : 'C', # Cysteine (Cys)
    'UGC' : 'C',
    'UGA' : 'Stop',
    'UGG' : 'W', # Tryptophan (Trp)
    'CUU' : 'L',
    'CUC' : 'L',
    'CUA' : 'L',
    'CUG' : 'L',
    'CCU' : 'P', # Proline (Pro)
    'CCC' : 'P',
    'CCA' : 'P',
    'CCG' : 'P',
    'CAU' : 'H', # Histidine (His)
    'CAC' : 'H',
    'CAA' : 'Q', # Glutamine (Gln)
    'CAG' : 'Q',
    'CGU' : 'R', # Arginine (Arg)
    'CGC' : 'R',
    'CGA' : 'R',
    'CGG' : 'R',
    'AUU' : 'I', # Isoleucine (Ile)
    'AUC' : 'I',
    'AUA' : 'I',
    'AUG' : 'M', # Methionine (Met)
    'ACU' : 'T', # Threonine (Thr)
    'ACC' : 'T',
    'ACA' : 'T',
    'ACG' : 'T',
    'AAU' : 'N', # Asparagine (Asn)
    'AAC' : 'N',
    'AAA' : 'K', # Lysine (Lys)
    'AAG' : 'K',
    'AGU' : 'S',
    'AGC' : 'S',
    'AGA' : 'R', 
    'AGG' : 'R',
    'GUU' : 'V', # Valine (Val)
    'GUC' : 'V',
    'GUA' : 'V',
    'GUG' : 'V',
    'GCU' : 'A', # Alanine (Ala)
    'GCC' : 'A',
    'GCA' : 'A',
    'GCG' : 'A',
    'GAU' : 'D', # Aspartic Acid (Asp)
    'GAC' : 'D',
    'GAA' : 'E', # Glutamic Acid (Glu)
    'GAG' : 'E',
    'GGU' : 'G', # Glycine (Gly)
    'GGC' : 'G',
    'GGA' : 'G',
    'GGG' : 'G'
}

reverse_codon = {
    'F' : ['UUU', 'UUC'], # Phenylalanine (Phe)
    'L' : ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], # Leucine (Leu)
    'S' : ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], # Serine (Ser)
    'Y' : ['UAU', 'UAC'], # Tyrosine (Tyr)
    'Stop' : ['UAA', 'UAG', 'UGA'], # Stop Codon (Stop)
    'C' : ['UGU', 'UGC'], # Cysteine (Cys)
    'W' : ['UGG'], # Tryptophan (Trp)
    'P' : ['CCU', 'CCC', 'CCA', 'CCG'], # Proline (Pro)
    'H' : ['CAU', 'CAC'], # Histidine (His)
    'Q' : ['CAA', 'CAG'], # Glutamine (Gln)
    'R' : ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], # Arginine (Arg)
    'I' : ['AUU', 'AUC', 'AUA'], # Isoleucine (Ile)
    'M' : ['AUG'], # Methionine (Met)
    'T' : ['ACU', 'ACC', 'ACA', 'ACG'], # Threonine (Thr)
    'N' : ['AAU', 'AAC'], # Asparagine (Asn)
    'K' : ['AAA', 'AAG'], # Lysine (Lys)
    'V' : ['GUU', 'GUC', 'GUA', 'GUG'], # Valine (Val)
    'A' : ['GCU', 'GCC', 'GCA', 'GCG'], # Alanine (Ala)
    'D' : ['GAU', 'GAC'], # Aspartic Acid (Asp)
    'E' : ['GAA', 'GAG'], # Glutamic Acid (Glu)
    'G' : ['GGU', 'GGC', 'GGA', 'GGG'] # Glycine (Gly)
}

monoisotopic_residue_mass = {
    # Daltons (Da)
    'F' : 147.06841, # Phenylalanine (Phe)
    'L' : 113.08406, # Leucine (Leu)
    'S' : 87.03203, # Serine (Ser)
    'Y' : 163.06333, # Tyrosine (Tyr)
    'C' : 103.00919, # Cysteine (Cys)
    'W' : 186.07931, # Tryptophan (Trp)
    'P' : 97.05276, # Proline (Pro)
    'H' : 137.05891, # Histidine (His)
    'Q' : 128.05858, # Glutamine (Gln)
    'R' : 156.10111, # Arginine (Arg)
    'I' : 113.08406, # Isoleucine (Ile)
    'M' : 131.04049, # Methionine (Met)
    'T' : 101.04768, # Threonine (Thr)
    'N' : 114.04293, # Asparagine (Asn)
    'K' : 128.09496, # Lysine (Lys)
    'V' : 99.06841, # Valine (Val)
    'A' : 71.03711, # Alanine (Ala)
    'D' : 115.02694, # Aspartic Acid (Asp)
    'E' : 129.04259, # Glutamic Acid (Glu)
    'G' : 57.02146 # Glycine (Gly)
}

# GETTER FUNCTIONS

def get_codons():
    '''
    Return the codon dictionary
    '''
    return codon

def get_reverse_codons():
    '''
    Return the codon dictionary with keys and values swapped
    '''
    # rev = {y:x for x,y in codon.items()}
    return reverse_codon

def get_masses():
    '''
    Return the residue mass dictionary
    '''
    return monoisotopic_residue_mass

# ALGORITHMS

def dna_nucleotides(dnaseq):
    '''
    Count the number of each nucleotide in a DNA sequence 
    '''
    dnaseq = dnaseq.lower()
    a = 0
    c = 0
    g = 0
    t = 0
    for i in dnaseq:
        if (i == "a"):
            a += 1
        elif (i == "c"):
            c += 1
        elif (i == "g"):
            g += 1
        elif (i == "t"):
            t += 1
    # print("Adenine: %d\nCytosine: %d\nGuanine: %d\nThymine: %d" % (a, c, g, t))
    # return {'A' : a, 'C' : c, 'G' : g, 'T' : t}
    return (a, c, g, t)

def rna_nucleotides(rnaseq):
    '''
    Count the number of each nucleotide in a RNA sequence 
    '''
    rnaseq = rnaseq.lower()
    a = 0
    c = 0
    g = 0
    u = 0
    for i in rnaseq:
        if (i == "a"):
            a += 1
        elif (i == "c"):
            c += 1
        elif (i == "g"):
            g += 1
        elif (i == "u"):
            u += 1
    # print("Adenine: %d\nCytosine: %d\nGuanine: %d\nUracil: %d" % (a, c, g, u))
    # return {'A' : a, 'C' : c, 'G' : g, 'U' : u}
    return (a, c, g, u)

def complementary(dnaseq):
    '''
    Get the complentary strand of a DNA sequence
    '''
    dnaseq = dnaseq.lower()
    compseq = ""
    for i in dnaseq:
        if (i == "a"):
            compseq += "T"
        elif (i == "c"):
            compseq += "G"
        elif (i == "g"):
            compseq += "C"
        elif (i == "t"):
            compseq += "A"
    return compseq

def transcribe(dnaseq):
    '''
    Get the corresponding mRNA sequence from a DNA sequence
    '''
    strand = input("Enter C (Coding/Sense) or (T) (Template/Antisense):")
    if (strand.upper() == "C"):
        dnaseq = complementary(dnaseq)
    
    dnaseq = dnaseq.lower()
    compseq = ""
    for i in dnaseq:
        if (i == "a"):
            compseq += "U"
        elif (i == "c"):
            compseq += "G"
        elif (i == "g"):
            compseq += "C"
        elif (i == "t"):
            compseq += "A"
    return compseq

def reverse_transcribe(rnaseq):
    '''
    Get the corresponding template DNA sequence from an mRNA sequence
    '''
    rnaseq = rnaseq.lower()
    tempseq = ""
    for i in rnaseq:
        if (i == "a"):
            tempseq += "T"
        elif (i == "c"):
            tempseq += "G"
        elif (i == "g"):
            tempseq += "C"
        elif (i == "u"):
            tempseq += "A"
    return tempseq

def gc_content(dnaseq):
    '''
    Calculate the GC content of a DNA sequence
    '''
    dnaseq = dnaseq.lower()
    gc = 0
    basepairs = 0
    for i in dnaseq:
        if (i == "g" or i == "c"):
            gc += 1
        basepairs += 1
    gc_pct = round((gc / basepairs) * 100,2)
    return gc_pct

def hamming(seq1, seq2):
    '''
    Find the Hamming distance between two DNA sequences
    '''
    s1 = min(seq1, seq2)
    s2 = max(seq1, seq2)
    diff = len(s2[len(s1):])
    for i in range(len(s1)):
        if (s1[i] != s2[i]):
            diff += 1
    return diff

def translate(rnaseq):
    '''
    Get the corresponding protein sequence from an mRNA sequence
    '''
    rnaseq = rnaseq.upper()
    protein = ""
    start = 0
    for i in range(len(rnaseq)):
        if (rnaseq[i:i+3] == "AUG"):
            start = i
        else:
            pass
    rnaseq = rnaseq[start:]

    while (len(rnaseq) > 2):
        first_codon = rnaseq[:3]
        if (codon[first_codon] == 'Stop'):
            break
        else:
            protein += codon[first_codon]
        rnaseq = rnaseq[3:]

    return protein

def reverse_translate(protein):
    '''
    Get a list of all possible RNA transcripts from a protein sequence
    '''
    protein = protein.upper()
    codons = []
    for i in protein:
        codons.append(reverse_codon[i])
    
    rnaseqs = list(itertools.product(*codons))
    for i in rnaseqs:
        "".join(i)

    return rnaseqs

def motif(dnaseq, motif):
    '''
    Return a list of the location(s) of the motif in the DNA sequence
    '''
    dnaseq = dnaseq.lower()
    motif = motif.lower()
    loc = []
    for i in range(len(dnaseq)):
        if (dnaseq[i:i + len(motif)] == motif):
            loc += [i]
            # print(i)
        else:
            pass
    return loc
    
def shared_motif(arg1, arg2, *args):
    '''
    Return the longest shared motif between 2+ DNA sequences
    '''
    motif = shared_motif2(arg1, arg2)
    for arg in args:
        motif = shared_motif2(motif, arg)
    return motif 

def shared_motif2(seq1, seq2):
    '''
    Return the longest shared motif between two DNA sequences
    '''
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    l1 = len(seq1)
    l2 = len(seq2)
    suff = [[0 for i in range(l2 + 1)] for j in range(l1 + 1)]
    length = 0
    row = 0
    col = 0

    for i in range(l1 + 1):
        for j in range(l2 + 1):
            if (i == 0 or j == 0):
                suff[i][j] = 0
            elif (seq1[i - 1] == seq2[j - 1]):
                suff[i][j] = suff[i - 1][j - 1] + 1
                if (length < suff[i][j]):
                    length = suff[i][j]
                    row = i
                    col = j
            else:
                suff[i][j] = 0
    
    if (length == 0):
        motif = []
    else:
        motif =  [' '] * length

        while (suff[row][col] != 0):
            length -= 1
            motif[length] = seq1[row - 1]
            row -= 1
            col -= 1
    
    return ''.join(motif)

def protein_mass(protein):
    '''
    Calculate the total weight (Da) of a protein
    '''
    protein = protein.upper()
    weight = 0
    for i in protein:
        weight += monoisotopic_residue_mass[i]
    return weight

def palindrome(dnaseq):
    '''
    Determine if a DNA sequence is palindromic
    '''
    dnaseq = dnaseq.lower()
    for i in range(len(dnaseq)):
        if (dnaseq[i] == "a"):
            if (dnaseq[-1 * (i + 1)] != "a"):
                return False
        elif (dnaseq[i] == "c"):
            if (dnaseq[-1 * (i + 1)] != "c"):
                return False
        elif (dnaseq[i] == "g"):
            if (dnaseq[-1 * (i + 1)] != "g"):
                return False
        elif (dnaseq[i] == "t"):
            if (dnaseq[-1 * (i + 1)] != "t"):
                return False
    return True

def reverse_palindrome(dnaseq):
    '''
    Determine if a DNA sequence is equal to its reverse complement
    '''
    dnaseq = dnaseq.lower()
    for i in range(len(dnaseq)):
        if (dnaseq[i] == "a"):
            if (dnaseq[-1 * (i + 1)] != "t"):
                return False
        elif (dnaseq[i] == "c"):
            if (dnaseq[-1 * (i + 1)] != "g"):
                return False
        elif (dnaseq[i] == "g"):
            if (dnaseq[-1 * (i + 1)] != "c"):
                return False
        elif (dnaseq[i] == "t"):
            if (dnaseq[-1 * (i + 1)] != "a"):
                return False
    return True

def splice(rnaseq, *args):
    '''
    Return the translated protein after splicing introns from the pre-mRNA sequence
    '''
    rnaseq = rnaseq.lower()
    for arg in args:
        for i in range(len(rnaseq)):
            if (rnaseq[i:i + len(arg)] == arg):
                rnaseq = rnaseq[:i] + rnaseq[i + len(arg):]
            else:
                pass

    return translate(rnaseq)

def superstring(*args):
    '''
    Return the shortest superstring containing all given strings
    A superstring serves as a candidate chromosome by parsimony
    '''
    seq = []
    for arg in args:
        seq.append(arg)
    l = len(seq)

    # Compute overlap
    ol = [[0] * l for i in range(l)]
    for i, x in enumerate(seq):
        for j, y in enumerate(seq):
            if (i != j):
                for k in range(min(len(x), len(y)), -1, -1):
                    if (x.endswith(y[:k])):
                        ol[i][j] = k
                        break
    
    # Calculate mask
    dp = [[0] * l for i in range(1 << l)]
    parent = [[None] * l for i in range(1 << l)]
    for i in range(1, 1 << l):
        for j in range(l):
            if ((i >> j) & 1):
                pmask = i ^ (1 << j)
                if (pmask == 0):
                    continue
                for k in range(l):
                    if ((pmask >> k) & 1):
                        v = dp[pmask][k] + ol[k][j]
                        if (v > dp[i][j]):
                            dp[i][j] = v
                            parent[i][j] = k
    
    # Reconstruct superstring
    perm = []
    mask = (1 << l) - 1
    i = max(range(l), key = dp[-1].__getitem__)
    while i is not None:
        perm.append(i)
        mask, i = mask ^ (1 << i), parent[mask][i]
    
    perm = perm[::-1]
    seen = [False] * l
    for j in perm:
        seen[j] = True
    perm.extend([j for j in range(l) if not seen[j]])

    sup = [seq[perm[0]]]
    for j in range(1, len(perm)):
        over = ol[perm[j - 1]][perm[j]]
        sup.append(seq[perm[j]][over:])
    
    return "".join(sup)

def spliced_motif(dnaseq, motif):
    '''
    Return a list of the first occurring indices of a motif as a subsequence of a DNA sequence
    '''
    dnaseq = dnaseq.lower()
    motif = motif.lower()
    ind = [0] * len(motif)
    for i in range(len(motif)):
        for j in range(len(dnaseq)):
            if (motif[i] == dnaseq[j] and j > max(ind)):
                ind[i] = j
                break
    return ind

def point_substitute(pos, bp, *args):
    '''
    Induce a point substitution at the given position in each DNA sequence
    '''
    seqs = []
    for arg in args:
        arg = arg.upper()
        arg = arg[:pos] + bp.upper() + arg[pos + 1:]
        seqs.append(arg)
    return seqs
    
def edit_distance(seq1, seq2):
    '''
    Calculate the edit distance (number of point mutations) between two sequences
    Point mutations: 
    '''
    l1 = len(seq1)
    l2 = len(seq2)
    dist = 0

    if (l1 > l2):
        dist = l1 - l2
        seq1[:dist]
    elif (l2 > l1):
        dist = l2 - l1
        seq2[:dist]

    for i in range(l1):
        if (seq1[i] != seq2[i]):
            dist += 1

    return dist
    





        








