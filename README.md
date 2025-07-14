# countsMQ4Targets
Compute Average MQ0 Ratios and Generate Summary Report
Introduction
This script calculates the average MQ0 ratios across multiple samples for specified genomic regions and generates a summary report. It utilizes bedtools and samtools to process BAM files, extract coverage information, and compute the average total depth, average MQ0 depth, and average ratio for each region. Finally, the script outputs a file containing aggregated information from all samples.
Features
- Parse Sample ID File: Reads a list of sample IDs from a specified file.
- Process BAM Files:
  - Uses bedtools coverage to calculate coverage for each sample in specified BED regions.
  - Filters out MQ0 (mapping quality 0) reads and computes their coverage.
- Calculate Average Ratios: For each region, computes the average total depth, average MQ0 depth, and average ratio across all samples.
- Generate Summary Report: Outputs a file containing the average depth, average ratio, and the number of samples exceeding ratios of 0.1, 0.2, and 0.3 for each region.
Dependencies
- Python 3.x
- argparse module (usually installed with Python)
- os and sys modules (built-in modules)
- bedtools
- samtools
Usage
Command Line Arguments
bash
usage: scriptname.py [-h] -s SAMPLEIDS -b BAMDIR -e BEDFILE -o OUTPUTDIR
Compute summary of average MQ0 ratios across samples.
optional arguments:
  -h, --help         show this help message and exit
  -s, --sampleids   Path to sample IDs file
  -b, --bamdir      Base directory containing BAM files
  -e, --bedfile     Path to BED file
  -o, --outputdir   Output directory for results
Example
Assuming your sample ID file is sampleids.txt, BAM files are stored in /data/bams, the BED file is regions.bed, and you want to output the results to the results directory, you can run the following command:
bash
python scriptname.py -s sampleids.txt -b /data/bams -e regions.bed -o results
Script Execution Flow
1. Parse Sample ID File: Reads all sample IDs from sampleids.txt.
2. Iterate Over Each Sample:
   - Constructs the path to the BAM file.
   - Uses bedtools coverage to calculate full sample coverage and saves it to allcoverage{sampleid}.txt.
   - Filters out MQ0 reads, generates mq0{sampleid}.bam, and converts it to mq0{sampleid}.bed.
   - Uses bedtools coverage to calculate MQ0 coverage and saves it to mq0coverage{sampleid}.txt.
   - Calls the calculateaveragemqratio.py script to compute the average MQ0 ratio for each sample and saves it to averagemq0ratio{sampleid}.txt.
   - Parses this file and aggregates the data into the regiondata dictionary.
   - Cleans up intermediate files (e.g., allcoverage{sampleid}.txt, mq0{sampleid}.bam, etc.).
3. Generate Summary Report: Aggregates data from all samples into summaryaveragemq0ratio.txt and outputs it to the specified outputdir.
Output Files
- averagemq0ratio{sampleid}.txt: Average MQ0 ratio file for each sample, stored in the outputdir.
- summaryaveragemq0ratio.txt: Summary report for all samples, containing the average total depth, average MQ0 depth, average ratio, and the number of samples exceeding ratios of 0.1, 0.2, and 0.3 for each region.
Notes
- Ensure that bedtools and samtools are correctly installed and added to your system's PATH.
- Ensure that all sample BAM files exist in the specified bamdir and that the paths are correct.
- The script will automatically create the output directory (if it doesn't exist), but make sure you have the appropriate write permissions.
- Intermediate temporary files generated during computation will be deleted afterward to save disk space.
Contact
For any questions or suggestions, please contact huzhuocheng2024@sinh.ac.cn or submit an issue on the GitHub repository.
