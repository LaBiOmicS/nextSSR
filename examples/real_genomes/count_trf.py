import subprocess
import os
import shutil
import glob

def parse_trf_output(trf_dat_file: str) -> int:
    if not os.path.exists(trf_dat_file):
        return 0
    count = 0
    with open(trf_dat_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 12 and parts[0].isdigit() and parts[1].isdigit():
                period = int(parts[2])
                if 1 <= period <= 6:
                    count += 1
    return count

def run_trf_count_benchmark():
    fasta_dir = os.path.abspath("examples/real_genomes/fasta")
    out_base_dir = os.path.abspath("examples/real_genomes/results")
    trf_bin = shutil.which("trf") or "trf"

    genomes = [
        {"id": "NC_000932.1", "name": "Arabidopsis thaliana Chloroplast", "file": "NC_000932.1_Arabidopsis_cp.fasta"},
        {"id": "NC_011033.1", "name": "Oryza sativa (Rice) Mitochondrion", "file": "NC_011033.1_Oryza_mt.fasta"},
        {"id": "NC_000913.3", "name": "Escherichia coli K-12 MG1655 Complete Genome", "file": "NC_000913.3_Ecoli_MG1655.fasta"},
        {"id": "S_cerevisiae_16Chr", "name": "Saccharomyces cerevisiae (Yeast 16 Chromosomes)", "file": "yeast_16_chromosomes.fasta"}
    ]

    trf_counts = {}

    for g in genomes:
        fasta_path = os.path.join(fasta_dir, g["file"])
        trf_work = os.path.join(out_base_dir, f"trf_count_{g['id']}")
        os.makedirs(trf_work, exist_ok=True)
        shutil.copy(fasta_path, os.path.join(trf_work, "seq.fasta"))

        # Run TRF with minscore 20 for short microsatellites
        subprocess.run([trf_bin, "seq.fasta", "2", "7", "7", "80", "10", "20", "500", "-h", "-ngs"], cwd=trf_work, capture_output=True, text=True)
        
        dat_files = glob.glob(os.path.join(trf_work, "*.dat"))
        total_trf_ssrs = 0
        for df in dat_files:
            total_trf_ssrs += parse_trf_output(df)
            
        trf_counts[g["id"]] = total_trf_ssrs
        shutil.rmtree(trf_work, ignore_errors=True)

    print("TRF SSR Counts (MinScore 20, Period 1-6bp):", trf_counts)
    return trf_counts

if __name__ == '__main__':
    run_trf_count_benchmark()
