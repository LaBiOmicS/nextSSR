import click
import os
import time
from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TimeRemainingColumn,
)
from rich.table import Table
from rich import box

from nextssr.config import SSRConfig
from nextssr.finder import SSRFinder
from nextssr.compound import CompoundSSRProcessor
from nextssr.primer import PrimerDesigner
from nextssr.utils import parse_fasta
from nextssr.models import ExecutionProvenance
from nextssr.provenance import FAIRProvenanceManager
from nextssr.artifacts import ArtifactManager

console = Console()

BANNER = r"""[bold cyan]
                     _  ____ ____  ____
 _ __   _____  _| |_/ ___/ ___||  _ \
| '_ \ / _ \ \/ / __\___ \___ \| |_) |
| | | |  __/>  <| |_ ___) |__) |  _ <
|_| |_|\___/_/\_\\__|____/____/|_| \_\
[/bold cyan]
[dim]Next-Generation High-Performance & FAIR-Compliant SSR Platform[/dim]
"""


class DefaultGroup(click.Group):
    """Click group allowing default subcommand fallback for direct FASTA file arguments."""

    def parse_args(self, ctx, args):
        if args and not args[0].startswith("-") and args[0] not in self.commands:
            args.insert(0, "run")
        return super().parse_args(ctx, args)


@click.group(cls=DefaultGroup)
def main():
    """nextSSR: Next-generation High-Performance & FAIR-compliant SSR identification and Primer Design platform."""
    pass


@main.command(name="run")
@click.argument("fasta_file", type=click.Path(exists=True))
@click.option(
    "--config",
    "-c",
    type=click.Path(exists=True),
    help="Path to nextSSR configuration file (nextssr.yaml / nextssr.json)",
)
@click.option(
    "--output-dir", "-o", default="results", help="Output directory for artifacts"
)
@click.option(
    "--threads",
    "-t",
    default=None,
    type=int,
    help="Number of parallel CPU worker threads/processes",
)
@click.option(
    "--design-primers/--no-primers",
    default=True,
    help="Automatically design PCR primers for identified SSRs",
)
@click.option(
    "--flank-len",
    default=150,
    help="Flanking sequence window size in bp (default: 150bp)",
)
@click.option(
    "--opt-tm", default=58.0, help="Optimal primer melting temperature Tm in °C"
)
@click.option(
    "--min-product-size", default=100, help="Minimum amplicon product size in bp"
)
@click.option(
    "--max-product-size", default=300, help="Maximum amplicon product size in bp"
)
@click.option(
    "--gpu/--no-gpu",
    default=False,
    help="Enable GPU hardware acceleration if CUDA is available",
)
@click.option(
    "--ro-crate/--no-ro-crate",
    default=True,
    help="Generate FAIR RO-Crate JSON-LD metadata",
)
def run_analysis(
    fasta_file: str,
    config: str,
    output_dir: str,
    threads: int,
    design_primers: bool,
    flank_len: int,
    opt_tm: float,
    min_product_size: int,
    max_product_size: int,
    gpu: bool,
    ro_crate: bool,
):
    """Run SSR identification and primer design on a FASTA sequence file."""
    start_time = time.time()
    console.print(Panel(BANNER, border_style="cyan"))

    # Load configuration
    if config:
        cfg = SSRConfig.from_file(config, threads=threads, use_gpu=gpu)
    else:
        cfg = SSRConfig(
            threads=threads or os.cpu_count() or 4,
            use_gpu=gpu,
            design_primers=design_primers,
            flank_len=flank_len,
            opt_tm=opt_tm,
            min_product_size=min_product_size,
            max_product_size=max_product_size,
            generate_ro_crate=ro_crate,
        )

    # Output parameters panel
    param_table = Table(show_header=False, box=box.SIMPLE)
    param_table.add_row("[bold]Input FASTA File:[/bold]", fasta_file)
    param_table.add_row("[bold]Output Directory:[/bold]", os.path.abspath(output_dir))
    param_table.add_row("[bold]Parallel Workers:[/bold]", f"{cfg.threads} CPU Cores")
    param_table.add_row(
        "[bold]Hardware Acceleration:[/bold]",
        "[green]GPU (CUDA Enabled)[/green]" if cfg.use_gpu else "CPU Vectorized",
    )
    param_table.add_row(
        "[bold]Primer Design Engine:[/bold]",
        (
            f"[green]Enabled[/green] (Tm: {cfg.opt_tm}°C | Amplicons: {cfg.min_product_size}-{cfg.max_product_size}bp)"
            if cfg.design_primers
            else "[yellow]Disabled[/yellow]"
        ),
    )
    param_table.add_row(
        "[bold]FAIR RO-Crate Provenance:[/bold]",
        (
            "[green]Enabled[/green]"
            if cfg.generate_ro_crate
            else "[yellow]Disabled[/yellow]"
        ),
    )

    console.print(
        Panel(
            param_table,
            title="[bold cyan]Run Parameters[/bold cyan]",
            border_style="dim",
        )
    )

    finder = SSRFinder(cfg, flank_len=cfg.flank_len)
    compound_proc = CompoundSSRProcessor(cfg)
    primer_designer = (
        PrimerDesigner(
            opt_tm=cfg.opt_tm,
            min_tm=cfg.min_tm,
            max_tm=cfg.max_tm,
            min_product_size=cfg.min_product_size,
            max_product_size=cfg.max_product_size,
        )
        if cfg.design_primers
        else None
    )

    # Multi-parallel streaming batch execution with Rich progress bar
    results = []
    fasta_stream = parse_fasta(fasta_file)

    total_ssrs = 0
    total_compounds = 0
    total_primers_designed = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("[cyan]Processing FASTA sequences...", total=None)

        for res in finder.analyze_batch_parallel(fasta_stream):
            res = compound_proc.process(res)

            if primer_designer:
                for ssr in res.ssrs:
                    pair = primer_designer.design_primers(
                        ssr.flank_5p, ssr.sequence, ssr.flank_3p
                    )
                    ssr.primer_pair = pair
                    if pair.status == "OK":
                        total_primers_designed += 1

            results.append(res)
            total_ssrs += len(res.ssrs)
            total_compounds += len(res.compounds)
            progress.update(
                task,
                description=f"[cyan]Analyzed {len(results)} sequences ({total_ssrs} SSRs found)...",
            )

    exec_time = f"{time.time() - start_time:.2f}s"

    # Save artifacts using ArtifactManager
    file_hash = FAIRProvenanceManager.compute_file_sha256(fasta_file)
    prov = ExecutionProvenance(
        threads_used=cfg.threads,
        device_used="GPU" if cfg.use_gpu else "CPU Multi-core",
        input_file_hash=file_hash,
        total_sequences=len(results),
        total_ssrs=total_ssrs,
        total_compounds=total_compounds,
        total_primers_designed=total_primers_designed,
        execution_time=exec_time,
        parameters_hash=cfg.get_hash(),
    )

    art_mgr = ArtifactManager(output_dir)
    artifacts = art_mgr.save_artifacts(results, cfg, prov, fasta_file)

    # Summary Results Table
    summary_table = Table(title="Execution Summary & Metrics", box=box.ROUNDED)
    summary_table.add_column("Metric", style="cyan", no_wrap=True)
    summary_table.add_column("Value", style="bold green")

    summary_table.add_row("Total Sequences Analyzed", str(len(results)))
    summary_table.add_row("Total SSRs Identified", str(total_ssrs))
    summary_table.add_row("Total Compound SSRs", str(total_compounds))
    summary_table.add_row(
        "Valid PCR Primer Pairs Designed", f"{total_primers_designed} / {total_ssrs}"
    )
    summary_table.add_row("Total Execution Time", exec_time)

    console.print(summary_table)

    # Artifacts panel
    art_table = Table(show_header=False, box=box.SIMPLE)
    art_table.add_row(
        "[bold]GFF3 Sequence Ontology Annotation:[/bold]", artifacts["gff"]
    )
    art_table.add_row("[bold]TSV Tabular & Primer Report:[/bold]", artifacts["tsv"])
    art_table.add_row("[bold]FAIR RO-Crate Provenance:[/bold]", artifacts["ro_crate"])
    art_table.add_row("[bold]Summary Text Statistics:[/bold]", artifacts["summary"])
    art_table.add_row("[bold]Execution Manifest:[/bold]", artifacts["manifest"])

    console.print(
        Panel(
            art_table,
            title="[bold green]Generated Run Artifacts[/bold green]",
            border_style="green",
        )
    )


@main.command(name="init-config")
@click.option(
    "--output", "-o", default="nextssr.yaml", help="Output configuration file path"
)
def init_config(output: str):
    """Generate a documented default nextSSR YAML configuration file."""
    config_path = SSRConfig.generate_default_config(output)
    console.print(
        f"[green]✓ Created default nextSSR configuration file at:[/green] [bold]{config_path}[/bold]"
    )


if __name__ == "__main__":
    main()
