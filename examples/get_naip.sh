#!/usr/bin/env bash
# Download ten NAIP acquisitions (2009 to 2020) of one quarter quad in North Dakota, 4 bands
# (red, green, blue, near infrared) at 5 m, from the opengeos/data sample collection on GitHub.
# NAIP imagery: USDA Farm Service Agency, public domain.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data
for d in 20090818 20100629 20120731 20140720 20150926 20160816 20170827 20180716 20190904 20200831; do
  curl -sSL -o "data/m_4709963_se_14_1_${d}.tif" \
    "https://raw.githubusercontent.com/opengeos/data/main/naip/m_4709963_se_14_1_${d}.tif"
done
echo "Saved $(ls data/*.tif | wc -l) images to examples/data/"
