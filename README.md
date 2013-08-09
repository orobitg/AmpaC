AmpaC
=====

A C version of the Antimicrobial Sequence Scanning System (AmpaC). The original one is implemented in perl and it can be found in http://tcoffee.crg.cat/apps/ampa/do

A theoretical approach to spot active regions in antimicrobial proteins

Quick Start
-----------

Clone the git repository on your computer with the following command:

	$ git clone git@github.com:orobitg/AmpaC.git

Compile AmpaC.c running the makefile:

	$ make

Usage
-----

Run AmpaC with the following parameters.

 	$ ampaC -in <input fasta file> [other options (optional)]

Available options:

 **-in      Fasta sequence input file (required)**
 **-t       Threshold value (default: 7)**
 **-w       Window size value (default: 0.225)**
 **-rf      File where store the program result**
 **-df      File where store the produced plot data (Only Nseq == 1)**
 **-gf      File where store the generated plot (Only Nseq == 1)**
 **-noplot  Skip plot creation step**
 **-help    This help information**

Dependencies 
------------

 * G++
 * Perl MATH Modules
    * MATH/CMD.pm
    * MATH/Round.p,
 * Bioperl


