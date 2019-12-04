#!/bin/bash

if [ "$TRAVIS_OS_NAME" != "windows" ]; then
    export DO_CMD=sudo
fi

python -m pip install numpy
curl -fsSL https://git.io/JerYI | GITHUB_USER=horta GITHUB_PROJECT=logaddexp bash
curl -fsSL https://git.io/JerYI | GITHUB_USER=EBI-Metagenomics GITHUB_PROJECT=imm bash
curl -fsSL https://git.io/JerYI | GITHUB_USER=EBI-Metagenomics GITHUB_PROJECT=nmm bash

sudo ldconfig
