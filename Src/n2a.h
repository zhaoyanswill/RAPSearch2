#ifndef __N2A_H
#define __N2A_H

const int TOTCODON = 64;
const char STOP_AA = '.';
const char STOP_AA1 = 's';
const char UNKNOWN_AA = 'Z';
const char UNKNOWN_AA0 = 'X';
const char UNKNOWN_AA1 = 'x';
const char ATYPICAL[] = ".sZXx";


const char* nt[TOTCODON] =
{
    "TTT",
    "TTC",
    "TTA",
    "TTG",
    "TCT",
    "TCC",
    "TCA",
    "TCG",
    "TAT",
    "TAC",
    "TAA",
    "TAG",
    "TGT",
    "TGC",
    "TGA",
    "TGG",
    "CTT",
    "CTC",
    "CTA",
    "CTG",
    "CCT",
    "CCC",
    "CCA", 
    "CCG",
    "CAT",
    "CAC",
    "CAA",
    "CAG",
    "CGT",
    "CGC",
    "CGA",
    "CGG",
    "ATT",
    "ATC",
    "ATA",
    "ATG",
    "ACT",
    "ACC",
    "ACA",
    "ACG",
    "AAT",
    "AAC",
    "AAA",
    "AAG",
    "AGT",
    "AGC",
    "AGA",
    "AGG",
    "GTT",
    "GTC",
    "GTA",
    "GTG",
    "GCT",
    "GCC",
    "GCA",
    "GCG",
    "GAT",
    "GAC",
    "GAA",
    "GAG",
    "GGT",
    "GGC",
    "GGA",
    "GGG",
};


char aa[TOTCODON] =
{
    'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', STOP_AA, STOP_AA, 'C', 'C', STOP_AA, 'W', 'L', 'L', 'L', 'L', 'P', 'P','P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', 'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', 'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'
};


#endif
