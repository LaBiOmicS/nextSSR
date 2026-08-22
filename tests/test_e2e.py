import os
import shutil
import tempfile
import pytest
from nextssr.config import SSRConfig
from nextssr.finder import SSRFinder
from nextssr.compound import CompoundSSRProcessor
from nextssr.primer import PrimerDesigner
from nextssr.artifacts import ArtifactManager
from nextssr.models import ExecutionProvenance
from nextssr.provenance import FAIRProvenanceManager

def test_full_pipeline_end_to_end():
    # 1. Prepare sample FASTA with multiple SSR types
    fasta_content = (
        ">seq1 Complete E2E Test Sequence\n"
        "GATTACAAGCTAGCTACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAACGTACGTACGTACGT"
        "AAAAAAAAAAAAAAAAAAAACACACACACACAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        "TGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACT\n"
    )

    tmp_dir = tempfile.mkdtemp()
    try:
        fasta_path = os.path.join(tmp_dir, "test.fasta")
        out_dir = os.path.join(tmp_dir, "output")

        with open(fasta_path, 'w') as f:
            f.write(fasta_content)

        # 2. Config generation & loading
        config_path = os.path.join(tmp_dir, "nextssr.yaml")
        SSRConfig.generate_default_config(config_path)
        assert os.path.exists(config_path)

        cfg = SSRConfig.from_file(config_path, threads=2)
        assert cfg.threads == 2

        # 3. Finder & Compound processing
        finder = SSRFinder(cfg, flank_len=100)
        compound_proc = CompoundSSRProcessor(cfg)
        primer_designer = PrimerDesigner(opt_tm=58.0, min_product_size=50, max_product_size=250)

        res = finder.analyze_sequence("seq1", fasta_content.split("\n")[1])
        assert len(res.ssrs) > 0

        res = compound_proc.process(res)
        assert len(res.compounds) > 0

        for ssr in res.ssrs:
            ssr.primer_pair = primer_designer.design_primers(ssr.flank_5p, ssr.sequence, ssr.flank_3p)
            assert ssr.motif_class != ""
            assert ssr.weber_class in ["Perfect", "Imperfect"]

        # 4. Artifact Manager export
        prov = ExecutionProvenance(
            threads_used=2,
            device_used="CPU Multi-core",
            input_file_hash="fakehash",
            total_sequences=1,
            total_ssrs=len(res.ssrs),
            total_compounds=len(res.compounds),
            total_primers_designed=sum(1 for s in res.ssrs if s.primer_pair and s.primer_pair.status == "OK"),
            execution_time="0.05s"
        )

        art_mgr = ArtifactManager(out_dir)
        artifacts = art_mgr.save_artifacts([res], cfg, prov, fasta_path)

        # 5. Verify all generated artifacts exist and are non-empty
        for name, filepath in artifacts.items():
            assert os.path.exists(filepath), f"Missing artifact: {name} at {filepath}"
            assert os.path.getsize(filepath) > 0, f"Empty artifact file: {name}"

        # 6. Verify GFF3 content
        with open(artifacts["gff"], 'r') as f:
            gff_text = f.read()
            assert "##gff-version 3" in gff_text
            assert "SO:0000289" in gff_text
            assert "weber_class=" in gff_text

        # 7. Verify TSV content
        with open(artifacts["tsv"], 'r') as f:
            tsv_text = f.read()
            assert "Seq_ID" in tsv_text
            assert "Motif_Class" in tsv_text
            assert "Weber_Classification" in tsv_text
            assert "Forward_Primer" in tsv_text

        # 8. Verify FAIR RO-Crate content
        with open(artifacts["ro_crate"], 'r') as f:
            crate_text = f.read()
            assert "@context" in crate_text
            assert "ro-crate-metadata.json" in crate_text

    finally:
        shutil.rmtree(tmp_dir)
