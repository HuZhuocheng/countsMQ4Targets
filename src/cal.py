import argparse

def parse_coverage_file(file_path):
    coverage = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3] if len(parts) > 3 else "Unknown"
            pos = int(parts[4])  # Position is in column 5 (1-based)
            depth = int(parts[5])
            key = (chrom, start, end, gene)
            if key not in coverage:
                coverage[key] = {}
            coverage[key][pos] = depth
    return coverage

def calculate_average_ratios(all_cov, mq0_cov, output_file):
    with open(output_file, 'w') as out:
        out.write("Chrom\tStart\tEnd\tGene\tAverage_Total_Depth\tAverage_MQ0_Depth\tAverage_Ratio\n")
        for key in all_cov:
            chrom, start, end, gene = key
            total_positions = 0
            total_all_depth = 0.0
            total_mq0_depth = 0.0
            total_ratio = 0.0
            for pos in all_cov[key]:
                all_depth = all_cov[key].get(pos, 0)
                mq0_depth = mq0_cov[key].get(pos, 0)
                total_all_depth += all_depth
                total_mq0_depth += mq0_depth
                if all_depth > 0:
                    ratio = mq0_depth / all_depth
                else:
                    ratio = 0 if mq0_depth == 0 else float('inf')
                total_ratio += ratio
                total_positions += 1
            if total_positions > 0:
                average_ratio = total_ratio / total_positions
                average_all_depth = total_all_depth / total_positions
                average_mq0_depth = total_mq0_depth / total_positions
            else:
                average_ratio = 0
                average_all_depth = 0
                average_mq0_depth = 0
            # Format average values to two decimal places
            out.write(f"{chrom}\t{start}\t{end}\t{gene}\t{average_all_depth:.2f}\t{average_mq0_depth:.2f}\t{average_ratio:.2f}\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate average MQ0 ratio for each region.")
    parser.add_argument("-a", "--all_coverage", required=True, help="Path to all_coverage.txt file")
    parser.add_argument("-m", "--mq0_coverage", required=True, help="Path to mq0_coverage.txt file")
    parser.add_argument("-o", "--output", required=True, help="Path to output average_mq0_ratio.txt file")

    args = parser.parse_args()

    all_cov = parse_coverage_file(args.all_coverage)
    mq0_cov = parse_coverage_file(args.mq0_coverage)
    calculate_average_ratios(all_cov, mq0_cov, args.output)

    print(f"Average MQ0 ratios written to {args.output}")

if __name__ == "__main__":
    main()
