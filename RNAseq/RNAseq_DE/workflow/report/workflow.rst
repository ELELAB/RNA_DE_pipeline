This workflow performs differential expression analysis on single- or paired-end RNA-seq data.
After adapter removal with `Cutadapt <http://cutadapt.readthedocs.io>`_, reads were mapped and gene counts were generated with `STAR <https://github.com/alexdobin/STAR>`_.
Gene counts of replicates were summed up using custom script. These counts are by default used as an input for DESeq2. Alternatively, `FeatureCounts <http://subread.sourceforge.net/>`_ output might be used as an input for further analysis. 
Integrated normalization and differential expression analysis was conducted with `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ following standard procedure as outlined in the manual.
