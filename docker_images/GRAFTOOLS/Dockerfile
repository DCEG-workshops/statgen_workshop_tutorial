# syntax=docker/dockerfile:1
###############################################################################
# GRAF + GrafAnc — genetic relationship & ancestry inference
#
#   GRAF     https://github.com/ncbi/graf            prebuilt binary + Perl + data
#   GrafAnc  https://github.com/jimmy-penn/grafanc   C++ source (needs htslib)
#
# Both tools are by Yumi (Jimmy) Jin, share reference-panel data, and are often
# used together, so this image ships both on the PATH.
#
# Two upstream issues were found and handled while building/testing this image:
#
#   1. GrafAnc's Makefile lists libraries BEFORE the object files on the link
#      line. Modern GNU ld (Ubuntu defaults to --as-needed) then drops libhts,
#      so `make` fails with "undefined reference to hts_open/bcf_*". This file
#      compiles GrafAnc explicitly with libraries AFTER the objects instead.
#
#   2. The prebuilt GRAF 2.3 binary (built for a ~2006 GNU/Linux 2.6.18 target)
#      reproducibly stack-smashes inside SavePopulationResults() when run on a
#      current glibc using the shipped `-pop` example, matching 0 SNPs. GRAF
#      ships NO source, so it cannot be rebuilt here. GrafAnc is the actively
#      maintained successor for ancestry; prefer it. See notes in the runtime
#      stage. If you must run GRAF-pop, an older-glibc base (e.g. ubuntu:20.04)
#      is more likely to tolerate the overflow — test with your own data.
#
# Reproducible by default: upstream revisions are pinned below. Override with
#   --build-arg GRAF_REF=... --build-arg GRAFANC_REF=...
###############################################################################

###############################################################################
# BUILD (on Mac) -> PUBLISH (ghcr.io) -> RUN (Biowulf / Singularity)
#
# !! ARCHITECTURE !! Biowulf compute nodes are x86-64 (amd64), and the GRAF
# binary shipped by NCBI is x86-64-ONLY. If you build on an Apple-Silicon Mac,
# a plain `docker build` produces an arm64 image that will NOT run on Biowulf
# (and can't even run GRAF locally). Always cross-build for linux/amd64:
#
# Target registry path MUST be lowercase, even though the org is DCEG-workshops:
#     ghcr.io/dceg-workshops/graftools
#
#   # 1. Log in once (PAT classic with 'write:packages'; export as $GHCR_PAT):
#   echo "$GHCR_PAT" | docker login ghcr.io -u <your-github-username> --password-stdin
#
#   # 2. Cross-build FOR amd64 and push in one step (buildx ships with Docker
#   #    Desktop; on Apple Silicon this runs under QEMU emulation, so it's slow):
#   docker buildx build --platform linux/amd64 \
#       -t ghcr.io/dceg-workshops/graftools:1.0.0 \
#       -t ghcr.io/dceg-workshops/graftools:latest \
#       --push .
#   # (swap --push for --load to keep the image locally instead of pushing)
#
# New packages are PRIVATE by default — set the package visibility to Public in
# its GHCR settings if workshop attendees should pull without authenticating.
#
#   # 3. On Biowulf, in an interactive/batch job (NOT the login node), convert
#   #    the Docker image to a Singularity image (.sif) by pulling from GHCR:
#   sinteractive --gres=lscratch:20
#   module load singularity                      # or apptainer, if needed
#   export SINGULARITY_CACHEDIR=/lscratch/$SLURM_JOB_ID
#   export SINGULARITY_TMPDIR=/lscratch/$SLURM_JOB_ID
#   # if the package is still PRIVATE, authenticate the pull:
#   #   export SINGULARITY_DOCKER_USERNAME=<your-github-username>
#   #   export SINGULARITY_DOCKER_PASSWORD=<PAT with read:packages>
#   singularity pull docker://ghcr.io/dceg-workshops/graftools:1.0.0
#
#   # 4. Run tools with `exec` (Biowulf auto-binds /data, /lscratch, etc.):
#   singularity exec --cleanenv graftools_1.0.0.sif \
#       grafanc /data/$USER/input.vcf.gz /data/$USER/out.txt
#   singularity exec --cleanenv graftools_1.0.0.sif graf
#
# Singularity notes: it runs as YOUR user (the image's non-root `graf` user and
# /work are ignored), so read/write real Biowulf paths (/data, /lscratch), not
# in-image paths. It DOES honor the image's env, so GRAFPATH/PATH still work;
# --cleanenv stops your host env (e.g. a stray PYTHONPATH) from leaking in.
#
# To link the published package to the DCEG-workshops repo hosting this file
# (so it lists on the repo and inherits access), set the source label in the
# runtime stage to a SINGLE repo URL, e.g.
#     org.opencontainers.image.source=https://github.com/DCEG-workshops/<repo>
# (it currently records upstream provenance, which does not auto-link.)
#
# TIP: building via a GitHub Actions workflow avoids BOTH pain points here —
# GitHub's runners are native amd64 (no slow emulation) and use the repo's
# built-in GITHUB_TOKEN (no PAT). You'd then just `singularity pull` on Biowulf.
# Ask if you'd like the workflow file.
###############################################################################

ARG GRAF_REF=e5e86b11c826fb3edfa779e0876cc154b18a7448
ARG GRAFANC_REF=b5e51c78b5004a8388418ce56755bb2f71f24eb5
ARG HTSLIB_VERSION=1.20

###############################################################################
# Stage 1 — builder: build htslib, compile GrafAnc from source
###############################################################################
FROM ubuntu:24.04 AS builder

ARG GRAF_REF
ARG GRAFANC_REF
ARG HTSLIB_VERSION
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        git \
        wget \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- htslib (GrafAnc reads VCF/BCF through it) ------------------------------
# Local-file use only, so remote-access backends are disabled to stay lean.
RUN wget -q "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
    && tar xjf "htslib-${HTSLIB_VERSION}.tar.bz2" \
    && cd "htslib-${HTSLIB_VERSION}" \
    && ./configure --prefix="/opt/htslib/${HTSLIB_VERSION}" \
                   --disable-libcurl --disable-gcs --disable-s3 \
    && make -j"$(nproc)" \
    && make install \
    && cd / && rm -rf "htslib-${HTSLIB_VERSION}" "htslib-${HTSLIB_VERSION}.tar.bz2"

# ---- GrafAnc: fetch, decompress reference panel, compile --------------------
RUN git clone https://github.com/jimmy-penn/grafanc.git /opt/grafanc \
    && git -C /opt/grafanc checkout "${GRAFANC_REF}"

WORKDIR /opt/grafanc/cpp
# The allele-frequency panel must be decompressed before use (per README).
RUN gunzip -f data/AncSnpPopAFs.txt.gz

# Explicit compile: libraries AFTER objects (fixes the Makefile link-order bug),
# and -O2 (upstream ships only -g / no optimization for this compute-heavy tool).
RUN g++ -O2 -std=c++11 -pthread \
        -I. -I"/opt/htslib/${HTSLIB_VERSION}/include" \
        -o grafanc *.cpp \
        -L"/opt/htslib/${HTSLIB_VERSION}/lib" -lhts -lz -lm -lstdc++ \
    && LD_LIBRARY_PATH="/opt/htslib/${HTSLIB_VERSION}/lib" ./grafanc 2>&1 | grep -q "GrafAnc"

# Drop what the runtime image doesn't need (test data ~40 MB, git history, objects).
RUN rm -rf /opt/grafanc/.git /opt/grafanc/medium_testing /opt/grafanc/cpp/*.o

# ---- GRAF: prebuilt binary + data (no source to build) ----------------------
RUN git clone https://github.com/ncbi/graf.git /opt/graf \
    && git -C /opt/graf checkout "${GRAF_REF}" \
    && rm -rf /opt/graf/.git

###############################################################################
# Stage 2 — runtime
###############################################################################
FROM ubuntu:24.04 AS runtime

ARG HTSLIB_VERSION
# Build a lean image without the plotting toolchains:
#   docker build --build-arg WITH_PLOTTING=false -t graftools:lean .
ARG WITH_PLOTTING=true
ENV DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="GRAF + GrafAnc" \
      org.opencontainers.image.description="NCBI GRAF (relationship/ancestry) + GrafAnc (ancestry inference)" \
      org.opencontainers.image.source="https://github.com/ncbi/graf https://github.com/jimmy-penn/grafanc"

# Runtime libraries: GRAF binary needs libstdc++; GrafAnc needs htslib's deps
# (zlib/bz2/lzma). Perl runs both projects' Perl helpers.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 \
        zlib1g \
        libbz2-1.0 \
        liblzma5 \
        perl \
    && rm -rf /var/lib/apt/lists/*

# Optional plotting stacks:
#   GRAF     -> Perl GD   (PlotGraf.pl, PlotPopulations.pl)
#   GrafAnc  -> matplotlib + pandas (python/PlotGrafAnc.py)
RUN if [ "$WITH_PLOTTING" = "true" ]; then \
        apt-get update && apt-get install -y --no-install-recommends \
            libgd-graph-perl \
            libcgi-pm-perl \
            python3 \
            python3-matplotlib \
            python3-pandas \
        && rm -rf /var/lib/apt/lists/* ; \
    fi

# htslib shared library for the GrafAnc binary, registered with the loader.
COPY --from=builder /opt/htslib/${HTSLIB_VERSION}/lib/ /opt/htslib/${HTSLIB_VERSION}/lib/
RUN echo "/opt/htslib/${HTSLIB_VERSION}/lib" > /etc/ld.so.conf.d/htslib.conf && ldconfig

# Tools + reference data.
COPY --from=builder /opt/grafanc/ /opt/grafanc/
COPY --from=builder /opt/graf/    /opt/graf/

# Put executables on PATH.
#   GrafAnc: symlinked; /proc/self/exe resolves the link back to /opt/grafanc/cpp,
#   so it locates data/AncSnpPopAFs.txt on its own. GRAFPATH below is a fallback.
#   GRAF Perl plotters are symlinked (they have no sibling-module imports).
RUN ln -s /opt/grafanc/cpp/grafanc         /usr/local/bin/grafanc \
    && ln -s /opt/graf/bin/graf            /usr/local/bin/graf \
    && ln -s /opt/graf/bin/graf_dups       /usr/local/bin/graf_dups \
    && ln -s /opt/graf/bin/PlotGraf.pl     /usr/local/bin/PlotGraf.pl \
    && ln -s /opt/graf/bin/PlotPopulations.pl /usr/local/bin/PlotPopulations.pl \
    && chmod +x /opt/graf/bin/*.pl /opt/grafanc/perl/*.pl

# GrafAnc reference-data fallback lookup.
ENV GRAFPATH=/opt/grafanc/cpp
# GrafAnc's Python plotters import sibling modules, so expose their directory.
# Run them as:  python3 /opt/grafanc/python/PlotGrafAnc.py ...
ENV PYTHONPATH=/opt/grafanc/python

# NOTE on GRAF reference data: FP_SNPs.txt, the G1000 panel, etc. live in
# /opt/graf/data. GRAF's `-pop`/`-rel` commands may need to be run from that
# directory (or given explicit paths) per the GRAF docs.

# Run as non-root; mount your data under /work.
# Ubuntu 24.04 ships a stock `ubuntu` user already at UID 1000, so free it first
# (keeping UID 1000 makes bind-mounted host files owned by the typical first
# Linux user writable inside the container).
RUN userdel -r ubuntu 2>/dev/null || true \
    && useradd -m -u 1000 graf \
    && mkdir -p /work && chown graf:graf /work
USER graf
WORKDIR /work

CMD ["bash"]
