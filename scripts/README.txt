README for HiCvis.py

HiCvis is a domain visualization tool for Hi-C data, which will plot a heat map
of the Hi-C matrix along with an overlay of one or two sets of TADs/domains.

You will need to have the following Python packages installed:

    numpy
    matplotlib
    scipy
    seaborn

Required input: Hi-C matrix
===========================

HiCvis accepts two types of Hi-C formats.

Format 1. Tab-delimited text file (such as Dixon et al. format from
http://chromosome.sdsc.edu/mouse/hi-c/download.html)

Format 2. Sparse matrix format (Rao et al. data). When using this format you
must also enter the resolution of the data with the -r flag


Optional arguments:
-------------------

  -r Resolution         Hi-C Resolution (only needed if using Rao data format)
     	Data resolution needed to parse Rao data format

  -b startBound endBound
                        Bounds for viewing window (optional)
    For viewing only a selected window of the chromosome rather than the entire
    chromosome at once

  -d1 domainFile1       TAD file
  -d2 domainFile2       second TAD file (optional)
      TAD files should contain two columns representing the start and end of
      each domain. Any text in the file will be ignored. Currently HiCvis
      supports the formats output by Armatus, TADtree, and Dixon (if
      chromosomes are separated into individual files), as well as any other
      domain finder with output files as described.


  -dr1 domainResolution1
                        Resolution of domains in domainFile1
  -dr2 domainResolution2
                        Resolution of domains in domainFile2
    Some domain finders (ie Dixon, some versions of Armatus) output domain
    locations referring to the genomic locus rather than the Hi-C bin, so the
    resolution is needed to show these with the Hi-C data


  -l1 legendName1       Legend name for first set of domains
  -l2 legendName2       Legend name for second set of domains
	If given, these strings will be the legend labels for the respective domains


-o outputFile         Filename for saved image file
    If given, instead of displaying the resulting image it will be saved under
    this filename. If no extension is given, .png is the default. Supported
    formats: bmp, eps, gif, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz,
    tif, tiff.



Example runs:
=============

The Hi-C data file for the first example can be found by downloading the Human
IMR90 Fibroblast Normalized Matrices data set from
http://chromosome.sdsc.edu/mouse/hi-c/download.html, and using the normalized
matrix from chromosome 20 (filename: nij.chr20)

Running the following line will show chromosome 20 from Dixon et al.
(resolution = 40kb), comparing Armatus consensus domains with Dixon domains

    python HiCvis.py -i example/nij.chr20 -d1 example/chr20armatus.consensus.txt -d2 example/chr20dixondomains.txt -dr2 40000 -l1 ArmatusConsensus -l2 Dixon


The Hi-C data file for this second example can be found by downloading the
GSE63525_GM12878_primary_intrachromosomal_contact_matrices.tar.gz file from
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525 . In this example we
used 10kb resolution, so the file will be chr17_10kb.RAWobserved.  Using Rao
data for chromosome 17, comparing consensus Armatus domains vs gamma = 0 for a
window from locus 2500 to 3000:

    python HiCvis.py -i example/chr17_10kb.RAWobserved -r 10000 -d1 example/RaoChr17_10kb.consensus.txt -d2 example/RaoChr17_10kb.gamma.0.0.txt -l1 ArmatusConsensus -l2 ArmatusGamma=0 -b 2500 3000
