skip_trimming = config["trimming"].get("skip", False)

if not skip_trimming:

    ruleorder: Trim_Adapters_PE > Trim_Adapters_SE


    rule Trim_Adapters_PE:
        input:
            get_fastq
        output:
            fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
            fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
            qc="results/trimmed/{sample}-{unit}.qc.txt"
        threads: 8
        params:
            adapters="-a {0} -A {0} {1}".format(config["trimming"]["adapter"], config["params"]["cutadapt-pe"]),
            extra="--minimum-length=30 -a A\{100\} -A A\{100\}"
        resources:
            mem_mb=64000
        # This is so that it only keeps a set of reads if both pairs are at least 20bp long
        # And, it performs 3' poly-A trimming on the reads
        log: "results/logs/cutadapt/{sample}-{unit}.log"
        wrapper:
            "v1.5.0/bio/cutadapt/pe"


    rule Trim_Adapters_SE:
        input:
            get_fastq
        output:
            fastq="results/trimmed/{sample}-{unit}.fastq.gz",
            qc="results/trimmed/{sample}-{unit}.qc.txt"
        threads: 8
        params:
            adapters="-a {0} {1}".format(config["trimming"]["adapter"], config["params"]["cutadapt-se"]),
            extra="--minimum-length=35 -a A\{100\}"
        resources:
            mem_mb=64000
        log: "results/logs/cutadapt/{sample}-{unit}.log"
        wrapper:
            "v1.5.0/bio/cutadapt/se"
