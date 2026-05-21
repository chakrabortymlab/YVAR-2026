#!/bin/bash
#SBATCH --job-name=array_coverage         
#SBATCH --ntasks-per-node=1                      	
#SBATCH --cpus-per-task=16                
#SBATCH --mem=64G                         
#SBATCH --time=04:00:00                   
#SBATCH --output=array_coverage_%j.out    
#SBATCH --error=array_coverage_%j.err     

set -euo pipefail

module purge
module load GCC/13.2.0 SAMtools/1.21

# Create the header for our summary report
SUMMARY="array_coverage_summary.tsv"
echo -e "Strain\tContig\tArray_Size(bp)\tMean_Depth\tMedian_Depth\tMin_Depth\tMax_Depth\tBreadth(>0X)" > "$SUMMARY"

for STRAIN in ISO1 A3 A4; do
    BAM="${STRAIN}_array_strict.bam"
    
    if [ ! -f "$BAM" ]; then
        echo ">>> Warning: Cannot find $BAM. Skipping $STRAIN..."
        continue
    fi

    echo "========================================================="
    echo " Processing strain: $STRAIN"
    echo "========================================================="

    DEPTH_OUT="${STRAIN}_raw_depth.tsv"
    TRACK_OUT="${STRAIN}_10kb_coverage.bedGraph"

    # 1. Generate and save the raw base-by-base depth file
    echo "[1/3] Generating raw depth file: $DEPTH_OUT..."
    samtools depth -a "$BAM" > "$DEPTH_OUT"

    # 2. Calculate Summary Statistics using the raw depth file
    echo "[2/3] Calculating summary statistics..."
    awk -v strain="$STRAIN" '
    BEGIN {
        count = 0; sum = 0;
        min = 9999999; max = 0;
        cov_gt_0 = 0;
    }
    {
        contig = $1;
        d = $3;
        depths[count] = d;
        count++;
        sum += d;
        if (d < min) min = d;
        if (d > max) max = d;
        if (d > 0) cov_gt_0++;
    }
    END {
        if (count > 0) {
            mean = sum / count;
            breadth = (cov_gt_0 / count) * 100;
            asort(depths);
            if (count % 2 == 1) {
                median = depths[int(count/2) + 1];
            } else {
                median = (depths[count/2] + depths[count/2 + 1]) / 2;
            }
            printf "%s\t%s\t%d\t%.2f\t%.2f\t%d\t%d\t%.2f%%\n", strain, contig, count, mean, median, min, max, breadth;
        } else {
            printf "%s\tN/A\t0\t0\t0\t0\t0\t0%%\n", strain;
        }
    }' "$DEPTH_OUT" >> "$SUMMARY"

    # 3. Generate a 10kb Rolling Window Track from the raw depth file
    echo "[3/3] Generating 10kb rolling window BedGraph for IGV..."
    WINDOW=10000
    
    awk -v w="$WINDOW" '
    BEGIN { sum=0; count=0; start=1; }
    {
        sum += $3;
        count++;
        if (count == w) {
            print $1"\t"start-1"\t"$2"\t"sum/w;
            sum=0; count=0; start=$2+1;
        }
    }
    END {
        if (count > 0) {
            print $1"\t"start-1"\t"$2"\t"sum/count;
        }
    }' "$DEPTH_OUT" > "$TRACK_OUT"

    echo ">>> $STRAIN complete."
done

echo "---------------------------------------------------------"
