# scMOCHA Data Merging - Usage Guide

## Script: EXPRCELLADMERGE_IMPROVED.py

**Author:** Chunjie Liu
**Date:** 2025-08-15
**Description:** Efficiently merge 577 10x h5 files into one AnnData object with memory management

## Features

✅ **Memory-efficient batch processing**
✅ **Rich logging and progress tracking**
✅ **Typer CLI interface**
✅ **Error handling for failed samples**
✅ **Automatic QC metrics calculation**
✅ **Both raw and normalized outputs**
✅ **Memory usage monitoring**

## Test Results

Successfully tested with 3 samples:
- **Input:** 3 samples, 17,221 cells, 61,806 genes
- **Output Raw:** 17,221 cells × 61,806 genes (340 MB)
- **Output Normalized:** 17,083 cells × 36,381 genes (669 MB)
- **Peak Memory:** 1.74 GB
- **Processing Time:** ~20 seconds

## Usage Examples

### Test with small subset
```bash
conda activate renv
cd /home/liuc9/github/scMOCHA-data/analysis

# Test with 10 samples, batch size 5
python EXPRCELLADMERGE_IMPROVED.py --max-samples 10 --batch-size 5

# Test with 50 samples, batch size 25
python EXPRCELLADMERGE_IMPROVED.py --max-samples 50 --batch-size 25
```

### Full processing (577 samples)
```bash
# Recommended for high-memory machines (64GB+ RAM)
python EXPRCELLADMERGE_IMPROVED.py --batch-size 25

# Conservative approach for lower-memory machines
python EXPRCELLADMERGE_IMPROVED.py --batch-size 15

# Very conservative for limited memory
python EXPRCELLADMERGE_IMPROVED.py --batch-size 10
```

### Custom output directory
```bash
python EXPRCELLADMERGE_IMPROVED.py \
    --output-dir /path/to/custom/output \
    --batch-size 25
```

## Memory Recommendations

| RAM Available | Recommended Batch Size | Expected Processing Time |
| ------------- | ---------------------- | ------------------------ |
| 16 GB         | 10-15 samples          | 2-3 hours                |
| 32 GB         | 20-25 samples          | 1-2 hours                |
| 64 GB+        | 30-40 samples          | 45-60 minutes            |

## Output Files

The script creates two files in the output directory:

1. **`merged_raw_data.h5ad`** - Raw count data with QC metrics
2. **`merged_normalized_data.h5ad`** - Filtered and normalized data

## Memory Monitoring

The script automatically reports peak memory usage. Monitor your system during processing:

```bash
# In another terminal, monitor memory usage
watch -n 10 'free -h && ps aux | grep python | grep EXPRCELLADMERGE'
```

## Error Handling

- Failed samples are logged but don't stop processing
- Batch failures are handled gracefully
- Progress is tracked with Rich progress bars
- Detailed logging for debugging

## Next Steps

1. **Test first:** Run with `--max-samples 50` to verify everything works
2. **Monitor memory:** Check available RAM before full run
3. **Adjust batch size:** Based on your machine's capabilities
4. **Full run:** Execute with all 577 samples

## Technical Details

- Uses conda environment: `/home/liuc9/tools/anaconda3/envs/renv`
- Input data: 577 samples from `gse_srrid_srrdir.csv`
- Processing: Batch-wise loading, QC metrics, filtering, normalization
- Output format: AnnData h5ad files
- Memory management: Garbage collection between batches
