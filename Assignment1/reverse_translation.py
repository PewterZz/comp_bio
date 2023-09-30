def transcription():
    """RNA transcription"""
    dna = input("Enter DNA sequence: ")
    
    rna = dna.replace("T", "U")
    
    print("RNA sequence:", rna)
    
    
def validate_dna(dna_seq):
    seqm = dna_seq.upper()
    print(seqm)
    valid = seqm.count("A") + seqm.count("T") + seqm.count("C") + seqm.count("G")
    
    if valid == len(seqm):
        print("VALID")
        return True
    
    else:
        print("INVALID")
        return False

def calfreq(seq):
    dictionary = {}
    
    for i in seq.upper():
        if i in dictionary:
            dictionary[i] += 1
        else:
            dictionary[i] = 1
            
    for element in dictionary:
        dictionary[element] = dictionary[element] / len(seq) * 100
    return dictionary

codon_table = {
    "UUU":"F","CUU":"L","AUU":"I","GUU":"V",
    "UUC":"F","CUC":"L","AUC":"I","GUC":"V",
    "UUA":"L","CUA":"L","AUA":"I","GUA":"V",
    "UUG":"L","CUG":"L","AUG":"M","GUG":"V",
    "UCU":"S","CCU":"P","ACU":"T","GCU":"A",
    "UCC":"S","CCC":"P","ACC":"T","GCC":"A",
    "UCA":"S","CCA":"P","ACA":"T","GCA":"A",
    "UCG":"S","CCG":"P","ACG":"T","GCG":"A",
    "UAU":"Y","CAU":"H","AAU":"N","GAU":"D",
    "UAC":"Y","CAC":"H","AAC":"N","GAC":"D",
    "UAA":"STOP","CAA":"Q","AAA":"K","GAA":"E",
    "UAG":"STOP","CAG":"Q","AAG":"K","GAG":"E",
    "UGU":"C","CGU":"R","AGU":"S","GGU":"G",
    "UGC":"C","CGC":"R","AGC":"S","GGC":"G",
    "UGA":"STOP","CGA":"R","AGA":"R","GGA":"G",
    "UGG":"W","CGG":"R","AGG":"R","GGG":"G"
}

def translation():
    """RNA translation"""
    dictionary = {}
    rna = input("Enter RNA sequence: ")
    
    codons = [rna[i:i+3] for i in range(0, len(rna), 3)]
    protein = ""
    for codon in codons:
        if codon in codon_table:
            protein += codon_table[codon]
            if codon in dictionary:
                dictionary[codon] += 1
            else:
                dictionary[codon] = 1 
        else:
            protein += "X"
    
    print("Protein sequence:", protein)
    print(dictionary)
    
def reverse_translation():
    """Reverse translation"""
    dictionary = {}
    protein = input("Enter protein sequence: ")
    
    codons = [protein[i:i+1] for i in range(0, len(protein), 1)]
    rna = ""
    for codon in codons:
        for key, value in codon_table.items():
            if value == codon:
                rna += key
                if key in dictionary:
                    dictionary[key] += 1
                else:
                    dictionary[key] = 1
                break
    print("RNA sequence:", rna)
    print(dictionary)
    
if __name__ == '__main__':
    reverse_translation()