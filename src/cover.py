import os
import sys
import argparse

def parse_sample_ids(sample_id_file):
    with open(sample_id_file, 'r') as f:
        sample_ids = [line.strip() for line in f if line.strip()]
    return sample_ids

def parse_average_ratio_file(file_path):
    ratios = {}
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3]
            avg_total_depth = float(parts[4])
            avg_mq0_depth = float(parts[5])
            avg_ratio = float(parts[6])
            key = (chrom, start, end, gene)
            if key not in ratios:
                ratios[key] = {
                    'total_depths': [],
                    'mq0_depths': [],
                    'ratios': []
                }
            ratios[key]['total_depths'].append(avg_total_depth)
            ratios[key]['mq0_depths'].append(avg_mq0_depth)
            ratios[key]['ratios'].append(avg_ratio)
    return ratios

def main():
    parser = argparse.ArgumentParser(description="Compute summary of average MQ0 ratios across samples.")
    parser.add_argument("-s", "--sample_ids", required=True, help="Path to sample IDs file")
    parser.add_argument("-b", "--bam_dir", required=True, help="Base directory containing BAM files")
    parser.add_argument("-e", "--bed_file", required=True, help="Path to BED file")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for results")

    args = parser.parse_args()

    sample_id_file = args.sample_ids
    bam_dir = args.bam_dir
    bed_file = args.bed_file
    output_dir = args.output_dir

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sample_ids = parse_sample_ids(sample_id_file)

    # Dictionary to store ratios for each region across all samples
    region_data = {}

    for sample_id in sample_ids:
        print(f"Processing sample: {sample_id}")
        bam_subdir = os.path.join(bam_dir, f"DLBCL_{sample_id}/03MarkBQSR/")
        bam_file = os.path.join(bam_subdir, f"BQSR_DLBCL_{sample_id}.bam")
        if not os.path.exists(bam_file):
            print(f"Warning: BAM file not found for sample {sample_id}: {bam_file}")
            continue

        # Compute all_coverage.txt
        all_coverage_file = f"all_coverage_{sample_id}.txt"
        cmd_all = f"bedtools coverage -d -a {bed_file} -b {bam_file} > {all_coverage_file}"
        os.system(cmd_all)

        # Filter MQ0 reads and compute mq0_coverage.txt
        mq0_bam = f"mq0_{sample_id}.bam"
        cmd_filter = f"samtools view -b -q 0 -e 'mapq == 0' {bam_file} > {mq0_bam}"
        os.system(cmd_filter)
        mq0_bed = f"mq0_{sample_id}.bed"
        cmd_bamtobed = f"bedtools bamtobed -i {mq0_bam} > {mq0_bed}"
        os.system(cmd_bamtobed)
        mq0_coverage_file = f"mq0_coverage_{sample_id}.txt"
        cmd_mq0 = f"bedtools coverage -d -a {bed_file} -b {mq0_bed} > {mq0_coverage_file}"
        os.system(cmd_mq0)

        # Compute average_mq0_ratio.txt for the sample
        average_ratio_file = os.path.join(output_dir, f"average_mq0_ratio_{sample_id}.txt")
        cmd_python = f"python calculate_average_mq_ratio.py -a {all_coverage_file} -m {mq0_coverage_file} -o {average_ratio_file}"
        os.system(cmd_python)

        # Parse the average_mq0_ratio file
        sample_ratios = parse_average_ratio_file(average_ratio_file)
        for region, data in sample_ratios.items():
            if region not in region_data:
                region_data[region] = {
                    'total_depths': [],
                    'mq0_depths': [],
                    'ratios': []
                }
            region_data[region]['total_depths'].extend(data['total_depths'])
            region_data[region]['mq0_depths'].extend(data['mq0_depths'])
            region_data[region]['ratios'].extend(data['ratios'])

        # Clean up intermediate files if needed
        os.remove(all_coverage_file)
        os.remove(mq0_bam)
        os.remove(mq0_bed)
        os.remove(mq0_coverage_file)

    # Write the summary to output_file
    output_file = os.path.join(output_dir, 'summary_average_mq0_ratio.txt')
    with open(output_file, 'w') as out:
        out.write("Chrom\tStart\tEnd\tGene\tAverage_Total_Depth\tAverage_MQ0_Depth\tAverage_Ratio\tSamples_above_0.1\tSamples_above_0.2\tSamples_above_0.3\n")
        for region, data in region_data.items():
            chrom, start, end, gene = region
            avg_total_depth = sum(data['total_depths']) / len(data['total_depths']) if data['total_depths'] else 0
            avg_mq0_depth = sum(data['mq0_depths']) / len(data['mq0_depths']) if data['mq0_depths'] else 0
            avg_ratio = sum(data['ratios']) / len(data['ratios']) if data['ratios'] else 0
            count_01 = sum(1 for r in data['ratios'] if r > 0.1)
            count_02 = sum(1 for r in data['ratios'] if r > 0.2)
            count_03 = sum(1 for r in data['ratios'] if r > 0.3)
            out.write(f"{chrom}\t{start}\t{end}\t{gene}\t{avg_total_depth:.2f}\t{avg_mq0_depth:.2f}\t{avg_ratio:.2f}\t{count_01}\t{count_02}\t{count_03}\n")

    print(f"Summary written to {output_file}")

if __name__ == "__main__":
    main()
