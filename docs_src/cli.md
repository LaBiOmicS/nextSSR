# Command Line Interface (CLI)

`nextSSR` includes a Command Line Interface built with `click` and `rich`.

```bash
nextssr [COMMAND] [OPTIONS]
```

---

## Commands

### 1. `run`
Execute SSR mining, primer design, and FAIR output package creation on a FASTA file.

```bash
nextssr run input.fasta -o output_dir --threads 8
```

Options:
- `-i, --input`: Target FASTA file (can be `.fasta`, `.fa`, `.fna` or `.gz`).
- `-o, --output`: Base output directory.
- `-c, --config`: Path to YAML/JSON configuration file.
- `-t, --threads`: Number of CPU threads (default: CPU core count).
- `--no-primers`: Disable PCR primer design.

---

### 2. `epcr`
Run *in silico* e-PCR simulation using designed primers against a genome FASTA file.

```bash
nextssr epcr --fasta genome.fasta --primers output_dir/primers/nextssr_primers.tsv
```

---

### 3. `init-config`
Generate a default `nextssr.yaml` configuration template file.

```bash
nextssr init-config -o my_config.yaml
```
