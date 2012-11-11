Prerequisites
-------------

MERmaid only has three main dependencies:
* Boost C++: Base, MPI and filesystem
* An MPI framework: We only test it with openmpi
* Google [sparsehash](http://code.google.com/p/sparsehash)

Currently, we are only running it on 64-bit Linux. Mac OS X is known to have
problems building.

Building
--------

Just run `make` in the top-level directory.

Running
-------

    ./mermaid output-file-prefix input-files ..

`input-files` must be a list of FASTQ files using the Illumina 1.3 quality
score format

The outputs will be 
* A set of files containing lists of initial k-mers after the initial counting
  phase. These will be named output-file-prefix.ufx.<num>
* A file containing all the contigs. This will be named 
  output-file-prefix.contigs.0

Mermaid will override default configurations if they are specified in
"mermaid.conf". Copy mermaid.conf.example to mermaid.conf and adjust parameter
values to override defaults.

TODO
----

Verify that this README is correct
