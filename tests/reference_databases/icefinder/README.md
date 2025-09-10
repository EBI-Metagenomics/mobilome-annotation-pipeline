# Test Database Preparation for ICEFinder 2 Lite

ICEFinder 2 lite requires three reference databases for proper functionality:

## ICEFinder 2 specific HMM Models
These models are sourced directly from [ICEFinder 2](https://github.com/EBI-Metagenomics/icefinder2) and have been converted from HMMER 2 format to HMMER 3. The complete model set is used without trimming down as it is small enough not to be problem.

## MacSyFinder Models
The [MacSyFinder](https://macsyfinder.readthedocs.io/en/latest/) models were compiled by the ICEFinder 2 team. For the test database, a subset of these models has been selected.

## UniProt Database (Prokka format)
A custom UniProt database formatted following Prokka custom headers has been created. This format is required because ICEFinder 2 extracts specific annotations from the database during analysis.