#!/usr/bin/env bash
sed 's/λ/$\\lambda$/g' README.md|sed 's/log(φ)/$log(\\varphi)$/g'|pandoc -s -V geometry:margin=1.5cm -o README.pdf
