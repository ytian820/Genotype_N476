# CHIP-CAD 1.0.5: add current R packages without upgrading existing tools.
FROM ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963

# Use CRAN's R binaries and r2u's CRAN package binaries for Ubuntu 22.04.
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        -o /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && echo 'deb [arch=amd64 signed-by=/etc/apt/trusted.gpg.d/cran_ubuntu_key.asc] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/' \
        > /etc/apt/sources.list.d/cran.list \
    && curl -fsSL https://r2u.stat.illinois.edu/ubuntu/dirk_eddelbuettel_pubkey.asc \
        -o /etc/apt/trusted.gpg.d/r2u_key.asc \
    && echo 'deb [arch=amd64 signed-by=/etc/apt/trusted.gpg.d/r2u_key.asc] https://r2u.stat.illinois.edu/ubuntu jammy main' \
        > /etc/apt/sources.list.d/r2u.list \
    && printf 'Package: *\nPin: release o=CRAN-Apt Project\nPin-Priority: 700\n' \
        > /etc/apt/preferences.d/99r2u

# Temporarily hold existing packages so adding R cannot upgrade existing tools.
RUN dpkg --get-selections > /tmp/base-selections \
    && dpkg-query -W -f '${binary:Package}\n' | xargs apt-mark hold > /dev/null \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        'r-base-core=4.6.2*' r-cran-tidyverse r-cran-data.table \
    && dpkg --set-selections < /tmp/base-selections \
    && rm /tmp/base-selections \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'library(tidyverse); library(data.table); stopifnot(getRversion() == "4.6.2", packageVersion("tidyverse") >= "2.0.0", packageVersion("data.table") >= "1.18.6.1"); print(R.version.string); print(packageVersion("tidyverse")); print(packageVersion("data.table"))'

LABEL org.opencontainers.image.version="1.0.5"
