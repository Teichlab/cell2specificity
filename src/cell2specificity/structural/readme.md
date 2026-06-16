# TCR-pMHC Structure Prediction Container

Standalone Singularity/Apptainer container for predicting TCR-pMHC binding using AlphaFold2.

### 1. Download the container

```bash
apptainer pull oras://ghcr.io/izu0421/c2s_structure:latest
```

This downloads `c2s_structure_latest.sif` (~14GB).

### 2. Prepare your input

Create a TSV file with TCR-pMHC data. See [test.tsv](test.tsv) for format:

```tsv
pdbid	organism	mhc_class	mhc	peptide	va	ja	cdr3a	vb	jb	cdr3b	sample_id
1oga	human	1	A*02:01	GILGFVFTL	TRAV27*01	TRAJ42*01	CAGAGSQGNLIF	TRBV19*01	TRBJ2-7*01	CASSSRSSYEQYF	1oga
```

### 3. Run predictions

```bash
apptainer run --nv c2s_structure_latest.sif \
  --input-tsv test.tsv \
  --out-root output \
  --tcrdock-dir /opt/tcrdock \
  --data-dir /opt/tcrdock/params \
  --model-params /opt/tcrdock/params/tcrpmhc_run4_af_mhc_params_891.pkl \
  --model-prefix /opt/models/res_af2_v2_lowFNR \
  --output-csv output/predictions.csv
```

**Output:** Results in `output/` directory with structure predictions (PDB), metrics (pLDDT, PAE), and binding predictions.
