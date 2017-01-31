# -*- coding: utf-8 -*-
"""
TODO: YOUR HEADER COMMENT HERE

@author: Kyle Combes

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
	I added unit tests for both of the other potential nucleotides to verify all were complemented properly.

            >>> get_complement('A')
            'T'
            >>> get_complement('C')
            'G'
            >>> get_complement('G')
            'C'
            >>> get_complement('T')
            'A'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    else:
        return 'C'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
	These unit tests should be good. All the function is doing is calling get_complement() and reversing the result.
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    res = ''
    for i in range(len(dna)-1,-1,-1):
        res += get_complement(dna[i])
    return res

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        return the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    end = len(dna)
    if end < 4:
        return '' #
    for i in range(3,len(dna)-2,3):
        try:
            if dna[i] == 'T':
                if dna[i+1] == 'A':
                    if dna[i+2] == 'A' or dna[i+2] == 'G':
                        end = i
                        break
                if dna[i+1] == 'G':
                    if dna[i+2] == 'A':
                        end = i
                        break
        except IndexError:
            print('Trying to index string %s using index %i' % (dna, i))
    return dna[0:end]


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    start_index = 0
    dna_len = len(dna)
    orfs = list()
    i = 0
    while i < dna_len - 3:
        # Find the index of the start codon
        if dna[i:i+3] == 'ATG':
            # Read ORF
            orf = rest_of_ORF(dna[i:])
            orfs.append(orf)
            i += len(orf) + 3
        else:
            i += 3
    return orfs

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    res = list()
    for start_index in range(0,3):
        res.extend(find_all_ORFs_oneframe(dna[start_index:]))
    return res


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    res = find_all_ORFs(dna)
    bkwd = get_reverse_complement(dna)
    res.extend(find_all_ORFs(bkwd))
    return res


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orf_list = find_all_ORFs_both_strands(dna)
    max_len_i = 0
    max_len = 0
    for i in range(0,len(orf_list)):
        if len(orf_list[i]) > max_len:
            max_len = len(orf_list[i])
            max_len_i = i
    return orf_list[max_len_i]

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    max_len = 0
    for i in range(0,num_trials):
       sdna = shuffle_string(dna)
       max_orf_len = len(longest_ORF(sdna))
       if (max_orf_len > max_len):
           max_len = max_orf_len
    return max_len

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    proteins = ''
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if len(codon) < 3: break
        protein = aa_table[codon]
        proteins += protein
    return proteins


def gene_finder(filename):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    from load import load_seq
    dna = load_seq(filename)
    threshold = longest_ORF_noncoding(dna, 100)
    orfs = find_all_ORFs_both_strands(dna)
    amino_acids = list()
    for i in range(0,len(orfs)):
        orf = orfs[i]
        if len(orf) > threshold:
            aa = coding_strand_to_AA(orf)
            amino_acids.append(aa)
    return amino_acids


if __name__ == "__main__":
    import doctest
    import argparse
    #doctest.testmod()
    #doctest.run_docstring_examples(coding_strand_to_AA, globals())
    parser = argparse.ArgumentParser(description='Accept passing DNA filename')
    parser.add_argument('filename', help='The text file to parse for DNA code')
    args = parser.parse_args()
    res = gene_finder(args.filename)
    print('Found %i amino acid codings' % len(res))
    f = open('result.txt', 'w')
    for line in res:
        f.write('%s\n' % line)
    f.close()
