# RAPSearch

##### Version 2.22 (with multi-threading support), released on Oct 27, 2014

##### Affiliation: School of Informatics and Computing, Indiana University, Bloomington

##### Developers:
 * Yongan Zhao <yongzhao@indiana.edu>
 * Yuzhen Ye <yye@indiana.edu>
 * Haixu Tang <hatang@indiana.edu>
 
The development of rapsearch was supported by NIH grant 1R01HG004908 to YY

rapsearch is free software under the terms of the GNU General Public License as
published by the Free Software Foundation.

#### Upgrading from RAPSearch1.0 to RAPSearch2.0

RAPSearch2 runs faster than RAPSearch1, uses less memory, and best of all,
RAPSearch2 supports multiple threads!!  You can define the number of threads that
you want to use by setting ``-z`` parameter (see below). **You need to re-run
prerapsearch to prepare new database files!!**

#### Before you start
RAPSearch means Reduced Alphabet based Protein similarity Search 

A linux/unix computer that has 4G memory should be sufficient for similarity
search against database of any size

RAPSearch2 on the web:

* [http://rapsearch2.sourceforge.net/](http://rapsearch2.sourceforge.net/)
* [http://omics.informatics.indiana.edu/mg/RAPSearch2](http://omics.informatics.indiana.edu/mg/RAPSearch2)


Please check the project home page for updates and newer version of the
rapsearch.

RAPSearch1 on the web:

* [http://omics.informatics.indiana.edu/mg/RAPSearch](http://omics.informatics.indiana.edu/mg/RAPSearch)

#### Installation

simply call: 

    ./install

The executable files:

* ``rapsearch``: the similarity search tool
* ``prerapsearch``: the program that generates similarity search database files
  to be used by rapsearch

Install RAPSearch2 on Mac OS X:
[https://github.com/zhaoyanswill/RAPSearch2/issues/3](https://github.com/zhaoyanswill/RAPSearch2/issues/3)
   
#### Prepare similarity search database files

    prerapsearch -d database(FASTA file) -n processed-db-file

###### Example 1:
   
    prerapsearch -d nogCOGdomN95.seq -n nogCOGdomN95_db
         
Input: ``nogCOGdomN95.seq``, a fasta file
	   
Outputs: ``nogCOGdomN95_db``, a binary file & a couple of other files

###### Example 2:
   
    prerapsearch -d nr_fasta -n nr_db
         
Input: ``nr_fasta`` (NCBI nr file in FASTA format)

Outputs: ``nr_db``

#### Similarity search

###### Example 1:
   
    rapsearch -q 4440037.3.dna.fa -d nogCOGdomN95_db -o 4440037.3.dna-vs-nogCOGdomN95 -z 4
       
Here ``-z`` is set to 4, so that four threads be used.
        
Input: 
	   
* ``4440037.3.dna.fa``: the query file, note if it is a file of short nucleotide sequences,      
* ``nogCOGdomN95_db``: the base name of the similarity search database 

Output:

* ``4440037.3.dna-vs-nogCOGdomN95.m8``: the similarty search result, one hit in one line, like ``-m 8`` output from blast note the only difference is that log10(E-value) is listed in the file, maximum 500 hits per query
* ``4440037.3.dna-vs-nogCOGdomN95.aln``: detailed alignments, maximum 100 alignments per query

###### Example 2:
    
    cat 4440037.4440037.3.dna | rapsearch -q stdin -d nogCOGdomN95_db -o 4440037.3.dna-vs-nogCOGdomN95 -z 4

Input: stdin of system, it's easy and convinient to incorporate RAPSearch into the pipeline

###### Example 3:
   
    rapsearch -q 4440037.3.dna.fa -d nogCOGdomN95_db -o 4440037.3.dna-vs-nogCOGdomN95 -z 4 -b 0

Output: Only output m8 file.

###### Example 4:
   
    rapsearch -q 4440037.3.dna.fa -d nogCOGdomN95_db -z 4 -u 1

Output: Only generate m8 file and output it to stdout of system. It's easy and
convinient to incorporate RAPSearch into the pipeline

**Notes**:  

* By default, rapsearch program uses only one thread; you may use `-z` to
change the threads.

* RAPSearch outputs hits with log\_10(E-value) < 1 (i.e., E-value of 10); you
may change this by setting log\_10(E-value) using `-e` option.

* You may also try `-v` and `-b` options to change the number of database
sequences that you want to output.

#### More notes 

   1. The sample files (e.g. ``4440037.3.dna.fa`` & ``nogCOGdomN95.seq``) are available
   at the rapsearch website; nr can be downloaded from NCBI ftp site.

   2. For big similarity search database (like nr), prerapsearch automatically
   splits the database into several files to reduce the memory requirement. If
   you need further reduce the memory usage, you can choose to use ``-s`` option to
   run prerapsearch
