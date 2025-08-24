# VARIANT Web Interface - Upload Guide

## Overview
This guide will walk you through uploading your virus data files to the VARIANT web interface. You can now upload your own virus data without needing pre-configuration!

## Step-by-Step Upload Process

### Option 1: Using Pre-configured Viruses
1. **Open the web interface** at `http://localhost:8000`
2. **Click on "Data Upload"** tab (it should be the first tab)
3. **Select a pre-configured virus** from the dropdown (SARS-CoV-2, HIV-1, etc.)
4. **Upload your files** following the steps below

### Option 2: Creating Your Own Custom Virus (Recommended for New Users)
1. **Open the web interface** at `http://localhost:8000`
2. **Click on "Data Upload"** tab
3. **Create a Custom Virus:**
   - Enter a unique virus name (e.g., "MyCustomVirus")
   - Choose virus type:
     - **Single-segment virus**: Most viruses (SARS-CoV-2, HIV-1, etc.)
     - **Multi-segment virus**: Viruses like Influenza A (H3N2) with multiple segments
   - For multi-segment viruses, enter segment names (e.g., "segment_1,segment_2,segment_3")
   - Click "Create Custom Virus"
4. **Your custom virus will appear** in the virus selection dropdown
5. **Select your custom virus** and proceed with file uploads

## Required Files for Analysis

For each virus, you need to upload **3 types of files**:

### 1. Reference Genome
- **Format**: FASTA (.fasta, .fa, .fna)
- **Content**: Complete genome sequence of the reference strain
- **Example**: `NC_045512.fasta` (SARS-CoV-2 reference genome)

### 2. Proteome
- **Format**: FASTA (.fasta, .fa, .faa)
- **Content**: Protein sequences translated from the genome
- **Example**: `SARS-CoV-2_proteome.fasta`

### 3. MSA File (Multiple Sequence Alignment)
- **Format**: Text or FASTA (.txt, .fasta, .fa, .aln, .clustal)
- **Content**: Multiple sequence alignment of virus strains
- **Example**: `gisaid_hcov-19_20221201_20221230_China_0_msa.txt`

## File Upload Process

### For Single-Segment Viruses:
1. **Select your virus** from the dropdown
2. **Choose file type** (Reference Genome, Proteome, or MSA)
3. **Select your file** from your computer
4. **Click "Upload File"**
5. **Repeat for all 3 file types**

### For Multi-Segment Viruses:
1. **Select your virus** from the dropdown
2. **Choose the segment** you want to upload files for
3. **Upload files for each segment** (Reference Genome, Proteome, MSA)
4. **Repeat for all segments**

## Upload Checklist

The interface includes an interactive checklist that shows:
- ✅ **Reference Genome**: Uploaded
- ✅ **Proteome**: Uploaded  
- ✅ **MSA File**: Uploaded

**Remember**: You must upload all 3 file types (reference genome, proteome, MSA) before you can start analysis!

## File Organization

Your uploaded files are automatically organized with user-friendly names:
- **Single-segment viruses**: `data/VirusName/refs/` and `data/VirusName/clustalW/`
- **Multi-segment viruses**: `data/VirusName/SegmentName/refs/` and `data/VirusName/SegmentName/clustalW/`

### File Naming Convention
Files are automatically renamed for clarity:
- **Reference Genome**: `VirusName_reference_genome.fasta`
- **Proteome**: `VirusName_proteome.fasta` 
- **MSA File**: `VirusName_msa.txt`

If you upload multiple files of the same type, they'll be numbered (e.g., `VirusName_proteome_1.fasta`, `VirusName_proteome_2.fasta`).

## After Upload

Once all files are uploaded:
1. **Go to the "Analysis" tab**
2. **Select your virus** (or custom virus) - the same one you uploaded files for
3. **Configure analysis parameters**:
   - **Process all genomes**: Recommended (processes all sequences in your MSA file)
   - **Detect frameshifts**: Optional (identifies potential frameshifting sites)
4. **Start the analysis** - the system will automatically use your uploaded files
5. **Monitor progress** in the "History" tab
6. **Download results** from the "Files" tab

**Note**: You don't need to upload MSA files again in the Analysis tab - the system uses the files you already uploaded in the Data Upload tab.

## Tips for Success

- **Use descriptive virus names** for custom viruses
- **Ensure file formats are correct** (check file extensions)
- **Upload files in any order** - the system will organize them
- **Check the upload checklist** to ensure all files are uploaded
- **Use the "Files" tab** to verify your uploaded files

## Troubleshooting

- **File upload fails**: Check file format and size
- **Virus not appearing**: Refresh the page or check virus name spelling
- **Analysis won't start**: Ensure all 3 file types are uploaded
- **Results not showing**: Check the "History" tab for job status

## Example Workflow

1. **Create custom virus** "MySARSVariant"
2. **Upload reference genome** `my_reference.fasta`
3. **Upload proteome** `my_proteome.fasta`
4. **Upload MSA file** `my_alignment.txt`
5. **Start analysis** in the Analysis tab
6. **Download results** from the Files tab

**Happy analyzing! 🧬🔬**
