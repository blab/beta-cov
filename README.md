# nextstrain.org/coronavirus

This is the [Nextstrain](https://nextstrain.org) build for coronaviruses, visible at
[nextstrain.org/coronavirus](https://nextstrain.org/coronavirus).

The build encompasses fetching data, preparing it for analysis, doing quality
control, performing analyses, and saving the results in a format suitable for
visualization (with [auspice][]).  This involves running components of
Nextstrain such as [fauna][] and [augur][].

All Zika-specific steps and functionality for the Nextstrain pipeline should be
housed in this repository.

## Usage

If you're unfamiliar with Nextstrain builds, you may want to follow our
[quickstart guide][] first and then come back here.

There are two main ways to run & visualise the output from this build:

The first, and easiest, way to run this pathogen build is using the [Nextstrain
command-line tool][nextstrain-cli]:
```
nextstrain build .
nextstrain view auspice/
```

See the [nextstrain-cli README][] for how to install the `nextstrain` command.

The second is to install augur & auspice using conda, following [these instructions](https://nextstrain.org/docs/getting-started/local-installation#install-augur--auspice-with-conda-recommended).
The build may then be run via:
```
snakemake
auspice --datasetDir auspice/
```

Build output goes into the directories `data/`, `results/` and `auspice/`.


## Configuration

Configuration takes place entirely with the `Snakefile`. This can be read top-to-bottom, each rule
specifies its file inputs and output and also its parameters. There is little redirection and each
rule should be able to be reasoned with on its own.


### fauna / RethinkDB credentials

This build starts by pulling sequences from our live [fauna][] database (a RethinkDB instance). This
requires environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY` to be set.

If you don't have access to our database, you can run the build using the
example data provided in this repository.  Before running the build, copy the
example sequences into the `data/` directory like so:

    mkdir -p data/
    cp example_data/zika.fasta data/


[Nextstrain]: https://nextstrain.org
[fauna]: https://github.com/nextstrain/fauna
[augur]: https://github.com/nextstrain/augur
[auspice]: https://github.com/nextstrain/auspice
[snakemake cli]: https://snakemake.readthedocs.io/en/stable/executable.html#all-options
[nextstrain-cli]: https://github.com/nextstrain/cli
[nextstrain-cli README]: https://github.com/nextstrain/cli/blob/master/README.md
[quickstart guide]: https://nextstrain.org/docs/getting-started/quickstart
