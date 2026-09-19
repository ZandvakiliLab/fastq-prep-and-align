rule jbrowse_create:
    output:
        touch("results/jbrowse/create"),
        html=os.path.join(config["jbrowse"]["dir"], "index.html"),
    conda:
        "../envs/jbrowse.yml"
    log: 
        "results/jbrowse/create.log",
    params:
        jbrowse_dir=lambda w, input: os.path.dirname("{output.html}")
    message:
        "create jbrowse folder"
    shell:
        """
        jbrowse create {params.jbrowse_dir} --force
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
        os.path.join(config["jbrowse"]["dir"], "index.html"),
    output:
        touch("results/jbrowse/add_assembly"),
        config=os.path.join(config["jbrowse"]["dir"], "config.json"),
    conda:
        "../envs/jbrowse.yml"
    log:
        "results/jbrowse/add_assembly.log",
    resources:
        file_lock=1,
    params:
        s3_url=lambda wc: "{url}/results/genome/genome.2bit".format(
            url=config["jbrowse"]["s3_url"]
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
        config=os.path.join(config["jbrowse"]["dir"], "config.json"),
    output:
        touch("results/jbrowse/add_anno"),
    conda:
        "../envs/jbrowse.yml"
    log:
        "results/jbrowse/add_anno.log",
    resources:
        file_lock=1,
    params:
        s3_url=lambda wc: "{url}/results/genome/genome.sorted.gff.gz".format(
            url=config["jbrowse"]["s3_url"]
        ),
        extra=config["jbrowse"]["add_anno"]["extra"],
    message:
        "add genome annotation to jbrowse"
    shell:
        """
        jbrowse add-track {params.s3_url} --target {input.config} {params.extra}
        """


# ─────────────────────────────────────────────────────────────────────────────
# 3. Add BigWigs
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_add_bw:
    input:
        config=os.path.join(config["jbrowse"]["dir"], "config.json"),
    output:
        touch(
            expand(
                "results/jbrowse/{sample}_{strand}_bw",
                sample=samples.index,
                strand=["plus", "minus"]
            )
        ),
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
        for i in {params.s3_url_plus}; do
            jbrowse add-track $i \
                --target {input.config} \
                --name "${{i##*/}}" \
                {params.extra}
        done

        for i in {params.s3_url_minus}; do
            jbrowse add-track $i \
                --target {input.config} \
                --name "${{i##*/}}" \
                --config '{{"displays":[{{"type":"LinearWiggleDisplay","displayId":"my_bw-LinearWiggleDisplay","inverted":true}}]}}' \
                {params.extra}
        done
        """

# ─────────────────────────────────────────────────────────────────────────────
# 3. Add BigWigs
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_add_cram:
    input:
        config=os.path.join(config["jbrowse"]["dir"], "config.json"),
    output:
        touch(
            expand(
                "results/jbrowse/{sample}_cram",
                sample=samples.index,
            )
        ),
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
        for i in {params.s3_url}; do
            jbrowse add-track $i \
                --indexFile $i.crai \
                --target {input.config} \
                --name "${{i##*/}}" \
                --config '{{"displays":[{{"type":"LinearPileupDisplay", "colorBySetting": {{"type": "strand"}}}}]}}' \
                {params.extra}
        done
        """
