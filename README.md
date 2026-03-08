# mosquito-virome-jock
This is a repository for my work in mosquito virome analysis. This is an ongoing project and the pipeline is alway in development. Transcriptome analysis and smallRNA analysis.

## 🦟 Overview
This project focuses on the identification and characterization of viral sequences in mosquito samples through dual-analysis approaches:
- **Transcriptome analysis**: Detection of viral transcripts and viral diversity
- **Small RNA analysis**: Identification of viral-derived small RNAs and antiviral RNAi responses
- 
## 🚧 Project Status

**⚠️ Active Development**: This pipeline is under continuous development. Methods, parameters, and recommendations may change as we optimized for different mosquito species and viral detection strategies.

## 📋 Pipeline Workflow "(some tools need to be updated)"

### Transcriptome Analysis Pipeline

This is a standard **Metatranscriptomics** or **Viral Discovery** pipeline. The goal is to take raw RNA sequencing data, clean it, remove the host (mosquito) genetic material, and then assemble and identify the remaining microbial (especially viral) sequences.
### The Metatranscriptomics/Viral Discovery Pipeline
---
#### Step 1: Input Data
- **File:** `1- Archivo Fastq.`
- **Action:** You start with your raw sequencing data. This is typically a FASTQ file (or paired files: `sample_R1.fastq` and `sample_R2.fastq`).
---
#### Step 2: Initial Quality Control (QC)
- **Tools:** `fastqc` / `multiqc`
- **Action:** Run FastQC on the raw FASTQ files to get a report on base quality, GC content, adapter contamination, and duplicate levels. Then run MultiQC to aggregate all reports into one single HTML overview.
---
#### Step 3: Trimming
- **Tools:** `trimmomatic` (or `trim_galore`)
- **Action:** Based on the QC report, trim low-quality bases (quality scores) and remove adapter sequences. This ensures that downstream assemblies are accurate.
- **Prompt:** "Run Trimmomatic on my raw FASTQ files to trim low-quality bases and remove adapters. Use default parameters for now."
---
#### Step 4: Quality Control (Again)
- **Tools:** `fastqc` / `multiqc`
- **Action:** It is essential to re-run FastQC/MultiQC on the *trimmed* data to verify that the trimming was successful and that the data is now clean.
- **Prompt:** "Run FastQC and then MultiQC on the trimmed FASTQ files generated in Step 3."
---
#### Step 5: Host Depletion (Remove Mosquito Reads)
- **Goal:** Remove reads that belong to the host organism (*Aedes* mosquitoes) to enrich for viral or other microbial sequences.

    - **5.1: Obtain Host References**
        - **Action:** Search NCBI for all *Aedes* genomes (drafts and complete).
        - **Prompt for AI:** "Help me construct an NCBI E-utilities command (esearch/efetch) to download all nucleotide sequences for the genus 'Aedes' (taxid:7157), including complete genomes and scaffolds. Save them to `aedes_all_genomes.fasta`."

    - **5.2: Build a "Super-reference" and Index**
        - **Action:** Concatenate all downloaded *Aedes* sequences into a single FASTA file. Then, build a STAR index for this file. STAR requires a genome directory and a FASTA file.
        - **Prompt:** "I have a file `aedes_all_genomes.fasta`. Generate a STAR genome index for it. Use --runThreadN [number_of_cores] and --runMode genomeGenerate. Name the output directory `STAR_Aedes_index`."
        - 
    - **5.3: Align Reads to Host**
        - **Action:** Align the trimmed reads from Step 3 to the *Aedes* super-reference using STAR.
        - **Prompt:** "Align my trimmed reads (R1 and R2) to the Aedes index (`STAR_Aedes_index`) using STAR. Use the `--outReadsUnmapped Fastx` parameter. This is crucial as it will output the reads that did *not* map to the mosquito genome."

    - **5.4: Obtain Unmapped Reads (The "Good" Reads)**
        - **Action:** STAR will produce a file containing the unmapped reads. These represent the putative non-host (viral, bacterial, etc.) sequences that we want to assemble. The file is usually named `Unmapped.out.mate1` and `Unmapped.out.mate2`.
        - **Prompt:** "Where does STAR save the unmapped reads? I used the `--outReadsUnmapped Fastx` option."

---
#### **This are not the final tools tehe pipeline is still to be stablished from this step**
#### Step 6: *De Novo* Assembly
- **Goal:** Assemble the unmapped reads into longer contiguous sequences (contigs) to reconstruct viral genomes.
- **Tools:** `rnaspades.py` (RNA version of SPAdes), `metaspades.py` (Metagenomic version), `Trinity`, `megahit`
- **Action:** Run multiple assemblers. Different assemblers have different strengths (genome complexity, transcript expression levels). Using several increases the chance of recovering diverse viruses.
- **Prompt:** "Assemble the unmapped reads (from Step 5.4) using multiple tools. First, run rnaspades.py. Then, run megahit. Finally, run Trinity."

---

#### Step 6a: Merge Assemblies
- **Tool:** `CAP3`
- **Action:** Concatenate the contigs from all assemblers into one large file. Then, run CAP3 to cluster and merge overlapping contigs, reducing redundancy and creating a final, non-redundant set of longer contigs.
- **Prompt for AI:** "I have contig files from SPAdes, MEGAHIT, and Trinity. Concatenate them into `all_assemblies.fasta` and then run CAP3 (`cap3 all_assemblies.fasta`) to merge redundant sequences."

---

#### Step 7: Initial Homology Search (Nucleotide)
- **Tool:** `blastn`
- **Database:** NCBI `nt` (non-redundant nucleotide database)
- **Action:** Compare the merged contigs against the entire `nt` database. This will quickly identify known viruses, as well as other organisms (bacteria, fungi) that might have been present in the sample.
- **Prompt:** "Run blastn with the CAP3 output file (`all_assemblies.fasta.cap.contigs`) against the NCBI `nt` database. Output the results in tabular format (`-outfmt 6`)."

---

#### Step 8: Sensitive Homology Search (Protein)
- **Tool:** `DIAMOND` (BLASTX alternative)
- **Database:** NCBI `nr` (non-redundant protein database)
- **Action:** Perform a translated nucleotide-to-protein search. This is much more sensitive for detecting divergent or novel viruses, as protein sequences are more conserved than nucleotide sequences. DIAMOND is used because it is extremely fast.
- **Prompt:** "Run DIAMOND BLASTX with my CAP3 contigs against the NCBI `nr` protein database. Use the `--sensitive` mode. Save the output in BLAST tabular format (`-f 6`)."

By the end of Step 8, you will have two lists of hits (nucleotide and protein) that you can analyze to identify which contigs are likely viral and what they are similar to.
 A diagram of a timeline
<img width="912" height="1536" alt="PIPELINE - visual selection" src="https://github.com/user-attachments/assets/6817c7e0-ac3e-4af2-94ce-fea772319c14" />


