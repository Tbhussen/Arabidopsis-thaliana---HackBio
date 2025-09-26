# UV Stress Analysis of Arabidopsis Thaliana

### Tamim Hussein
### HackBio Internship: Stage 2 Project

# Overview of the Task

This project investigates how different tissues of the Arabidopsis thaliana leaf—epidermis, mesophyll, and vasculature—respond to multiple abiotic stresses. Specifically, our focus is on determining the genes that respond to UV stress in the vasculature. The core analysis involves performing differential gene expression analysis between UV-C treated vasculature tissues and Water (control) samples.

# About Arabidopsis thaliana

Arabidopsis thaliana is a small flowering plant widely used as a model organism in plant biology. It is a member of the mustard (Brassicaceae) family, which includes cultivated species such as cabbage and radish. While Arabidopsis itself is not of major agronomic significance, its small genome size (~135 Mb) and ease of cultivation offer significant advantages for fundamental research in genetics and molecular biology. The haploid chromosome number is five.

# Pipeline Description

The analysis pipeline consisted of the following steps:

1. **Data Acquisition:** Raw sequencing data were obtained from the European Nucleotide Archive, consisting of six samples:
    - **Control (Water):** SRR12808527, SRR12808528, SRR12808529
    - **UV-C Treated:** SRR12808497, SRR12808498, SRR12808499

2. **Quality Control & Processing:** Quality checks and trimming were performed on the raw data files.

3. **Mapping to Reference Genome:** Processed reads were mapped to the Arabidopsis thaliana reference genome.

4. **Indexing & BAM File Generation:** The mapped reads were indexed and converted into BAM files.

5. **Feature Counting:** Using the STAR tool and reference annotation from Ensembl Plants, feature counts were generated from the BAM files.

6. **Differential Expression Analysis:** R code was used to perform differential expression analysis using DESeq2, yielding lists of upregulated and downregulated genes.

7. **Pathway Analysis:** Upregulated genes were analyzed using ShinyGO to identify pathways stimulated by UV-C treatment.

# Results

## Top 100 Upregulated Genes

```
AT3G19184  AT1G77655  AT4G22710
AT3G49620  AT1G80820  AT2G04070
AT1G56250  AT1G19250  AT1G43160
AT3G57460  AT1G51890  AT3G01830
AT1G32350  AT5G22540  AT2G37430
AT5G39670  AT5G67080  AT5G40000
AT2G35980  AT2G29470  AT5G25920
AT5G36925  AT5G40990  AT1G05990
AT1G64290  AT5G59680  AT1G47890
AT1G68765  AT1G35230  AT1G57560
AT1G13480  AT5G24640  AT1G05575
AT1G61120  AT3G48640  AT4G15975
AT2G23270  AT5G47740  AT1G21240
AT1G14540  AT3G14700  AT5G42380
AT1G20310  AT5G13080  AT2G18193
AT1G57630  AT5G38310  AT3G26830
AT4G15150  AT3G23220
```

## Top 5 Upregulated Pathways

1. **Photosynthesis-antenna proteins**
2. **Phenylalanine, tyrosine and tryptophan biosynthesis**
3. **Plant-pathogen interaction**
4. **Zeatin biosynthesis**
5. **Alpha-linolenic acid metabolism**

### Scientific Understanding of Pathway Stimulation

- **Photosynthesis-antenna proteins:** These proteins are responsible for capturing light energy and transferring it to the photosynthetic reaction centers. Upregulation under UV stress suggests an adaptive response to maximize light harvesting and protect against UV-induced damage.
- **Phenylalanine, tyrosine and tryptophan biosynthesis:** These aromatic amino acids serve as precursors for a wide array of secondary metabolites, many of which act in stress defense and signaling. Their increased biosynthesis indicates a heightened defense state under UV stress.
- **Plant-pathogen interaction:** Genes associated with this pathway may be upregulated due to cross-talk between abiotic (UV) and biotic stress response pathways, enhancing the plant’s overall resilience.
- **Zeatin biosynthesis:** Zeatin is a type of cytokinin, a plant hormone involved in cell division and growth. Its biosynthesis may be increased to promote recovery and regeneration following UV-induced cellular damage.
- **Alpha-linolenic acid metabolism:** This pathway leads to the production of signaling molecules like jasmonic acid, which are involved in plant stress and defense responses, further strengthening the plant's adaptive mechanisms.

---

_End of document._
