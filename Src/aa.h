int tot_codon = 64;
int tot_aa = 20;


//    A  R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   @
//    0  1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19

//gbmr.10
char gbmr10s[500][20] = {"G", "D", "N", "AEFIKLMQRVW", "Y", "H", "C", "T", "S", "P"};
//group                0    1    2    3                        4    5    6    7    8   9
char gbmr10r[] = "GDNAYHCTSP"; //representative residue for each group
char gbmr10[500] = {3, 3, 2, 1, 6, 3, 3, 0, 5, 3, 3, 3, 3, 3, 9, 8, 7, 3, 4, 3};
//               A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V

//dayhoff.6
char dayhoff6s[500][20] = {"AGPST", "C", "DENQ", "FWY", "HKR", "ILMV"};
//group                     0        1    2       3      4      5
char dayhoff6r[] = "ACDFHI";
char dayhoff6[500] = {0, 4, 2, 2, 1, 2, 2, 0, 4, 5, 5, 4, 5, 3, 0, 0, 0, 3, 3, 5};
//                    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V

//murphy based on BLOSUM50, Protein Engineering, Vol. 13, No. 3, 149-152, 2000
//murphy.5
char murphy5s[500][20] = {"LVIMC", "ASGTP", "FYW", "EDNQ", "KRH"};
//group                     0        1       2      3      4
char murphy5r[] = "LAFEK";
char murphy5[500] = {1, 4, 3, 3, 0, 3, 3, 1, 4, 0, 0, 4, 0, 2, 1, 1, 1, 2, 2, 0};
//                   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V

//murphy.10
char murphy10s[500][20] = {"A", "KR", "EDNQ", "C", "G", "H", "ILVM", "FYW", "P", "S/T" };
//group                     0    1      2          3    4    5    6          7        8    9
char murphy10r[] = "AKECGHIFPS";
char murphy10[500] = {0, 1, 2, 2, 3, 2, 2, 4, 5, 6, 6, 1, 6, 7, 8, 9, 9, 7, 7, 6};
//                   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V

//murphy.9
char murphy9s[500][20] = {"A", "KREDNQ", "C", "G", "H", "ILVM", "FYW", "P", "S/T" };
//group                     0    1        2    3    4    5       6      7    8
char murphy9r[] = "AKCGHIFPS";
char murphy9[500] = {0, 1, 1, 1, 2, 1, 1, 3, 4, 5, 5, 1, 5, 6, 7, 8, 8, 6, 6, 5};
//                   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V


//hsdm.4
char hsdm4s[500][20] = {"LIVFMYW", "C", "DNTSKEQRAGP", "H"};
//group               0          1    2             3
char hsdm4r[] = "LCDH";
char hsdm4[500] = {2, 2, 2, 2, 1, 2, 2, 2, 3, 0, 0, 2, 0, 0, 2, 2, 2, 0, 0, 0};
//              A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V

char aabet20s[500][20] = {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
char aabet20r[] = "ARNDCQEGHILKMFPSTWYV";
char aabet20[500] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

char aabet[500] = "ARNDCQEGHILKMFPSTWYV@";
//                 0123456789...
char aaindex[500] = {0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20};
//                   A   B  C  D  E   F  G  H  I  J   K   L   M   N  O   P   Q  R  S   T   U   V   W   X   Y   Z

char aa3bet[500][5] = {"ALA", "ARG", "ASN", "ASP", "CYS",
                       "GLN", "GLU", "GLY", "HIS", "ILE",
                       "LEU", "LYS", "MET", "PHE", "PRO",
                       "SER", "THR", "TRP", "TYR", "VAL",
                       "TER"   /* Terminal codon */
                      };

char aa3codon[500][5] = {"GGG", "GGA", "GGC", "GGT",
                         "GAG", "GAA", "GAC", "GAT",
                         "GCG", "GCA", "GCC", "GCT",
                         "GTG", "GTA", "GTC", "GTT",
                         "AGG", "AGA", "AGC", "AGT",
                         "AAG", "AAA", "AAC", "AAT",
                         "ACG", "ACA", "ACC", "ACT",
                         "ATG", "ATA", "ATC", "ATT",
                         "CGG", "CGA", "CGC", "CGT",
                         "CAG", "CAA", "CAC", "CAT",
                         "CCG", "CCA", "CCC", "CCT",
                         "CTG", "CTA", "CTC", "CTT",
                         "TGG", "TGA", "TGC", "TGT",
                         "TAG", "TAA", "TAC", "TAT",
                         "TCG", "TCA", "TCC", "TCT",
                         "TTG", "TTA", "TTC", "TTT"
                        };

char aa1codon[500] =
{
    7,7,7,7,
    6,6,3,3,
    0,0,0,0,
    19,19,19,19,
    1,1,15,15,
    11,11,2,2,
    16,16,16,16,
    12,9,9,9,
    1,1,1,1,
    5,5,8,8,
    14,14,14,14,
    10,10,10,10,
    17,20,4,4, /*TGA*/
    20,20,18,18, /*TAG, TAA */
    15,15,15,15,
    10,10,13,13
};

/*
 20 char gl_Ter = [ "TAA", "TAG", "TGA" ]
13	gl_Phe = [ "TTT", "TTC" ]
10	gl_Leu = [ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "CTN" ]
15	gl_Ser = [ "TCT", "TCC", "TCA", "TCG", "TCN", "AGT", "AGC" ]
18	gl_Tyr = [ "TAT", "TAC" ]
4	gl_Cys = [ "TGT", "TGC" ]

17	gl_Trp = [ "TGG" ]
14	gl_Pro = [ "CCT", "CCC", "CCA", "CCG", "CCN" ]
8	gl_His = [ "CAT", "CAC" ]
5	gl_Gln = [ "CAA", "CAG" ]
1	gl_Arg = [ "CGT", "CGC", "CGA", "CGG", "CGN", "AGA", "AGG" ]

9	gl_Ile = [ "ATT", "ATC", "ATA" ]
12	gl_Met = [ "ATG" ]
16	gl_Thr = [ "ACT", "ACC", "ACA", "ACG", "ACN" ]
2	gl_Asn = [ "AAT", "AAC" ]
11	gl_Lys = [ "AAA", "AAG" ]

19	gl_Val = [ "GTT", "GTC", "GTA", "GTG", "GTN" ]
0	gl_Ala = [ "GCT", "GCC", "GCA", "GCG", "GCN" ]
3	gl_Asp = [ "GAT", "GAC" ]
6	gl_Glu = [ "GAA", "GAG" ]
7	gl_Gly = [ "GGT", "GGC", "GGA", "GGG", "GGN" ]
*/

