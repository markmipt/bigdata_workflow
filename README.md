# Massive proteogenomic reanalysis of publicly available proteomic datasets of human tissues in search for protein recoding via adenosine-to-inosine RNA editing

This repository hosts the source code for the comutational workflow described in:

> Massive proteogenomic reanalysis of publicly available proteomic datasets of human tissues in search for protein
> recoding via adenosine-to-inosine RNA editing. biorXiv, 2022. DOI: https://doi.org/10.1101/2022.11.10.515815

Installation
------------

```
pip install .
```

Usage
-----

```
bdrunner --help
```

As input, the script needs a directory with raw files and preprocessed FASTA databases.
As output, it produces tables with FDR-controlled variant peptides, as well as intermediate files.
The script contains nine consecutive steps as follows:

1. **Conversion.** Raw data are converted to mzML using MSconvert from ProteoWizard.
   Conversion is done using multiple filters: “peakPicking true 1-”, “MS2Deisotope”, “zeroSamples removeExtra” and “threshold absolute 1 most-intense”.
2. **Canonical search.** All mzML files are searched with IdentiPy search engine against the canonical protein database.
   IdentiPy search engine is run along with Dinosaur feature detection algorithm to extract MS1 peptide isotopic cluster intensities for identified MS2 spectra.
   Scavager postsearch algorithm is used to process the output of IdentiPy.
3. **Variant search.** All the data are processed again using IdentiPy/Scavager combination to search against the peptide-based variant database.
4. **Brute-force search.** One hundred of top-scored canonical proteins are chosen for every processed dataset file.
   These proteins and their decoy counterparts are saved as separate FASTA files (one for each run).
   Using these protein databases, the brute-force searches are performed using IdentiPy search engine.
   Briefly, brute-force approach checks every possible single amino acid change for all canonical tryptic peptides.
   The idea is that LC-MS/MS data is enriched with modified peptides from top scored proteins.
   We assume that the most identified "single amino acid substitutions" for these proteins are just artifacts produced by abundant modifications.
5. **Variant tables generation.** A table with variant peptides is generated for each file in the dataset.
   Variant peptides are sequences which belong only to variant proteins (target or decoy) and not shared with any canonical proteins.
6. **Group-specific FDR evaluation.** All tables of variants are combinedand the list of identifications is filtered to 5% group-specific FDR
   using posterior error probabilities calculated by Scavager.
7. **Prosit MS/MS prediction.** MS/MS peak intensities are predicted using PROSIT for 5% FDR filtered variant peptides
   and for 1000 random canonical peptides identified. Next, correlations between experimental and predicted MS/MS spectra are calculated
   and the threshold of correlation is calculated as the 5-th percentile for canonical peptides.
   Variant peptides are filtered where correlations are below the threshold.
   The idea here is that modifications mimicking amino acid substitutions should result in lower correlations between predicted and experimental spectra.
8. **DeepLC RT prediction.** Theoretical retention times are predicted using DeepLC for 5% FDR filtered variant peptides
   and for all identified canonical peptides.
   The RT prediction errors are calculated as the difference between predicted and experimental peptide retention times.
   The distribution of RT errors for canonical peptides is fitted as a sum of normal and uniform distributions.
   The mean shift and standard deviation are estimated and used to filter out variant peptides with RT errors outside of the 3 standard deviation range.
   The idea here is that modifications mimicking amino acid substitutions should result in a higher difference between predicted and experimental RT.
9. **Final table composition.** The final variant peptide table is composed, containing all information about the identifications,
   including the number of identified PSMs, best MS/MS correlation, lowest RT error, amino acid substitution rate in brute-force search, etc.
