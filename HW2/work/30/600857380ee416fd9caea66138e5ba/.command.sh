#!/bin/bash -ue
bedtools genomecov -ibam aligned.bam -d > coverage.txt

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
