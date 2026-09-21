rule jbrowse_create:
    output:
        touch("results/jbrowse/create"),
    conda:
        "../envs/jbrowse.yml"
    log: 
        "results/jbrowse/create.log",
    message:
        "create jbrowse folder"
    params:
        output_dir="jbrowse",
    shell:
        """
        jbrowse create {params.output_dir} --force
        """


# ─────────────────────────────────────────────────────────────────────────────
# 1. Add Assembly
# ─────────────────────────────────────────────────────────────────────────────


rule faToTwoBit_fa:
    input:
        "results/genome/genome.fasta",
    output:
        "results/genome/genome.2bit",
    log:
        "results/genome/genome.fa_to_2bit.log",
    wrapper:
        "v7.1.0/bio/ucsc/faToTwoBit"


rule jbrowse_add_assembly:
    input:
        fa="results/genome/genome.2bit",
        jbrowse_created="results/jbrowse/create",
    output:
        config=temp("jbrowse/config_assembly.json"),
    conda:
        "../envs/jbrowse.yml"
    log:
        "results/jbrowse/add_assembly.log",
    resources:
        file_lock=1,
    params:
        s3_url=lambda wc, input: "{url}/{fa}".format(
            url=config["jbrowse"]["s3_url"], 
            fa="results/genome/genome.2bit"
        ),
        extra=config["jbrowse"]["add_assembly"]["extra"],
    message:
        "add genome assembly to jbrowse"
    shell:
        """
        jbrowse add-assembly {params.s3_url} --type twoBit --target {output.config} {params.extra} --force
        """

# ─────────────────────────────────────────────────────────────────────────────
# 2. Add Annotation
# ─────────────────────────────────────────────────────────────────────────────


rule sort_gff:
    input:
        gff="results/genome/genome.gff",
    output:
        gff="results/genome/genome.sorted.gff.gz",
    conda:
        "../envs/jbrowse.yml"
    log:
        "results/genome/genome_sort_gff.log",
    message:
        "sort gff3"
    shell:
        """
        jbrowse sort-gff {input.gff} | bgzip >{output.gff}
        """


rule index_gff:
    input:
        "results/genome/genome.sorted.gff.gz",
    output:
        "results/genome/genome.sorted.gff.gz.tbi",
    conda:
        "../envs/jbrowse.yml"
    log:
        "results/genome/genome_index_gff.log",
    message:
        "index gff3"
    shell:
        """
        tabix {input}
        """


rule jbrowse_add_anno:
    input:
        gff="results/genome/genome.sorted.gff.gz",
        gff_tbi="results/genome/genome.sorted.gff.gz.tbi",
        config="jbrowse/config_assembly.json",
    output:
        config=temp("jbrowse/config_anno.json"),
    conda:
        "../envs/jbrowse.yml"
    log:
        "results/jbrowse/add_anno.log",
    resources:
        file_lock=1,
    params:
        s3_url=lambda wc, input: "{url}/{gff}".format(
            url=config["jbrowse"]["s3_url"], 
            gff="results/genome/genome.sorted.gff.gz"
        ),
        extra=config["jbrowse"]["add_anno"]["extra"],
    message:
        "add genome annotation to jbrowse"
    shell:
        """
        cp {input.config} {output.config}
        jbrowse add-track {params.s3_url} --target {output.config} {params.extra} --force
        """


# ─────────────────────────────────────────────────────────────────────────────
# 3. Add BigWigs
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_add_bw:
    input:
        config="jbrowse/config_anno.json",
    output:
        config=temp("jbrowse/config_bw.json"),
    conda:
        "../envs/jbrowse.yml"
    log:
        "results/jbrowse/add_bw.log",
    resources:
        file_lock=1,
    params:
        s3_url_plus=expand(
            os.path.join(config["jbrowse"]["s3_url"], "results/deeptools/coverage/{sample}.plus.bw"),
            sample=samples.index,
        ),
        s3_url_minus=expand(
            os.path.join(config["jbrowse"]["s3_url"], "results/deeptools/coverage/{sample}.minus.bw"),
            sample=samples.index,
        ),
        extra=config["jbrowse"]["add_bw"]["extra"],
    message:
        "add plus bw tracks to jbrowse"
    shell:
        """
        cp {input.config} {output.config}
        for i in {params.s3_url_plus}; do
            jbrowse add-track $i \
                --target {output.config} \
                --name "${{i##*/}}" \
                --force \
                {params.extra}
        done

        for i in {params.s3_url_minus}; do
            jbrowse add-track $i \
                --target {output.config} \
                --name "${{i##*/}}" \
                --force \
                --config '{{"displays":[{{"type":"LinearWiggleDisplay","displayId":"my_bw-LinearWiggleDisplay","inverted":true}}]}}' \
                {params.extra}
        done
        """

# ─────────────────────────────────────────────────────────────────────────────
# 3. Add BigWigs
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_add_cram:
    input:
        config="jbrowse/config_bw.json",
    output:
        config="jbrowse/config.json",
    conda:
        "../envs/jbrowse.yml"
    log:
        "results/jbrowse/add_cram.log"
    resources:
        file_lock=1,
    params:
        s3_url=expand(
            os.path.join(config["jbrowse"]["s3_url"], "results/processed_alignment/cram/{sample}.cram"),
            sample=samples.index
        ),
        extra=config["jbrowse"]["add_cram"]["extra"],
    message:
        "add plus cram tracks to jbrowse"
    shell:
        """
        cp {input.config} {output.config}
        for i in {params.s3_url}; do
            jbrowse add-track $i \
                --indexFile $i.crai \
                --target {output.config} \
                --name "${{i##*/}}" \
                --force \
                --config '{{"displays":[{{"type":"LinearPileupDisplay", "colorBySetting": {{"type": "strand"}}}}]}}' \
                {params.extra}
        done
        """
