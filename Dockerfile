FROM mambaorg/micromamba:1.5.10-noble

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    git git-lfs sudo curl wget ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

RUN micromamba create -y -p /opt/conda/envs/bioenv \
    -c conda-forge -c bioconda \
    nextflow \
    ucsc-genepredtobed==482 \
    ucsc-gtftogenepred==482 \
    samtools=1.22.1 \
    duckdb-cli=1.3.2 \
    bedtools=2.31.1 \
    seqkit=2.10.1 \
    python=3.12.12 \
    && micromamba clean --all -y

COPY --from=ghcr.io/astral-sh/uv:0.10.9 /uv /usr/local/bin/uv

# Tell uv not to download or manage Python — use the conda-managed interpreter
ENV UV_PYTHON_DOWNLOADS=never \
    UV_PYTHON=/opt/conda/envs/bioenv/bin/python \
    UV_LINK_MODE=copy \
    UV_NO_CACHE=1

# Create venv in user-writable location
RUN uv venv /home/$MAMBA_USER/.venv --python /opt/conda/envs/bioenv/bin/python

# venv binaries take precedence; conda env provides nextflow, samtools, etc.
ENV VIRTUAL_ENV=/home/$MAMBA_USER/.venv \
    PATH=/home/$MAMBA_USER/.venv/bin:/opt/conda/envs/bioenv/bin:$PATH

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt .
RUN uv pip install -r requirements.txt

COPY --chown=$MAMBA_USER:$MAMBA_USER . .

RUN which python && python --version && which nextflow && nextflow -version

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python"]
CMD ["-c", "import time; time.sleep(float('inf'))"]