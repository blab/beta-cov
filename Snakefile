rule all:
    input:
        auspice_json = "auspice/beta-cov.json",

rule files:
    params:
        input_fasta = "data/beta-cov.fasta",
        include = "config/include.txt",
        dropped_strains = "config/dropped_strains.txt",
        reference = "config/betacov_reference.gb",
        auspice_config = "config/auspice_config.json",
        colors = "config/colors.tsv",
        description = "config/description.md"

files = rules.files.params

rule download:
    message: "Downloading sequences from fauna"
    output:
        sequences = "data/beta-cov.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location source locus authors url title journal puburl host virus_species"
    shell:
        """
        python3 ../fauna/vdb/download.py \
            --database vdb \
            --virus coronavirus \
            --fasta_fields {params.fasta_fields} \
            --resolve_method choose_genbank \
            --path $(dirname {output.sequences}) \
            --fstem $(basename {output.sequences} .fasta)
        """

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = rules.download.output.sequences
    output:
        sequences = "results/sequences.fasta",
        metadata = "results/metadata.tsv"
    params:
        fasta_fields = "strain virus accession date region country division city db segment authors url title journal paper_url host virus_species",
        prettify_fields = "region country division city"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields} \
            --prettify-fields {params.prettify_fields}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - excluding strains in {input.exclude}
          - minimum genome length of {params.min_length}
        """
    input:
        sequences = rules.parse.output.sequences,
        metadata = rules.parse.output.metadata,
        include = files.include,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered.fasta"
    params:
        group_by = "country",
        sequences_per_group = 50,
        min_length = 5000
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --include {input.include} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        mafft --retree 1 --nofft --thread 2 results/filtered.fasta > results/aligned.fasta
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --root ZBCoV \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        auspice_config = files.auspice_config,
        colors = files.colors,
        description = files.description
    output:
        auspice_json = rules.all.input.auspice_json
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --description {input.description} \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
