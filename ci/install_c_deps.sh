#!/bin/bash

python -m pip install numpy
curl -fsSL https://git.io/JerYI | GITHUB_USER=horta GITHUB_PROJECT=logaddexp DO_CMD=sudo bash
curl -fsSL https://git.io/JerYI | GITHUB_USER=EBI-Metagenomics GITHUB_PROJECT=imm DO_CMD=sudo bash
curl -fsSL https://git.io/JerYI | GITHUB_USER=EBI-Metagenomics GITHUB_PROJECT=nmm DO_CMD=sudo bash
