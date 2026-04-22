process COVERAGE {

    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    path bam
    path ref

    output:
    path "coverage.txt"
    path "coverage.png"

    script:
    """
    bedtools genomecov -ibam $bam -d > coverage.txt

    python - <<EOF

import matplotlib.pyplot as plt

x = []
y = []

with open("coverage.txt") as f:
    for line in f:
        chrom, pos, cov = line.strip().split()
        x.append(int(pos))
        y.append(int(cov))

plt.plot(x, y)
plt.xlabel("Position")
plt.ylabel("Coverage")
plt.savefig("coverage.png")
EOF
    """
}
