"""
Entrypoint for bacterial-genome

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from .util import snake_base, print_version, copy_config, run_snakemake, OrderedCommands, print_citation


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option('--output', help='Output directory', type=click.Path(),
                     default='bacterial_genome.out', show_default=True),
        click.option('--configfile', default='config.yaml', help='Custom config file', show_default=True),
        click.option('--threads', help='Number of threads to use', default=1, show_default=True),
        click.option('--use-conda/--no-use-conda', default=True, help='Use conda for Snakemake rules',
                     show_default=True),
        click.option('--conda-prefix', default=snake_base(os.path.join('workflow', 'conda')),
                     help='Custom conda env directory', type=click.Path(), show_default=False),
        click.option('--snake-default', multiple=True,
                     default=['--rerun-incomplete', '--printshellcmds', '--nolock', '--show-failed-logs'],
                     help="Customise Snakemake runtime args", show_default=True),
        click.argument('snake_args', nargs=-1)
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(cls=OrderedCommands)
def cli():
    """For more options, run:
    bacterial_genome command --help"""
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
bacterial_genome run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           bacterial_genome run --input [file]
Specify threads:    bacterial_genome run ... --threads [threads]
Disable conda:      bacterial_genome run ... --no-use-conda 
Change defaults:    bacterial_genome run ... --snake-default="-k --nolock"
Add Snakemake args: bacterial_genome run ... --dry-run --keep-going --touch
Specify targets:    bacterial_genome run ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--input', '_input', help='Input file/directory', type=str, required=True)
@click.option('--preprocess', help="sequencing method", default='paired', show_default=True,
                     type=click.Choice(['paired', 'longread']))

@common_options
def run(_input, preprocess, configfile, output, threads, use_conda, conda_prefix, snake_default,
        snake_args, **kwargs):
    """Run bacterial-genome"""

    # copy default config file if missing
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))

    # Config to add or update in configfile
    merge_config = {
        'input': _input,
        'output': output,
        'sequencing': preprocess,
        }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'run.smk')),   # Full path to Snakefile
        configfile=configfile,
        merge_config=merge_config,
        threads=threads,
        use_conda=use_conda,
        conda_prefix=conda_prefix,
        snake_default_args=snake_default,
        snake_extra=snake_args,
    )


@click.command()
@click.option('--configfile', default='config.yaml', help='Copy template config to file', show_default=True)
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(config)
cli.add_command(citation)


def main():
    print_version()
    cli()


if __name__ == '__main__':
    main()
