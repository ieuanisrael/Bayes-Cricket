# Bayes Cricket — R + rstan (no host setup beyond Docker)
#
# Build:
#   docker build -t bayes-cricket .
#
# Run example (fits model, prints rankings). Use -t for a TTY so logs stream:
#   docker run --rm -t bayes-cricket
#
# Interactive R (override entrypoint):
#   docker run --rm -it --entrypoint R bayes-cricket
#
# Mount your own innings CSV (expects data/innings.csv inside the mount):
#   docker run --rm -v /path/on/host/data:/app/data bayes-cricket
#
# More CPU for parallel chains:
#   docker run --rm -e MC_CORES=4 bayes-cricket

FROM rocker/r-ver:4.4.2

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    util-linux \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages('rstan', repos='https://cloud.r-project.org', dependencies=TRUE, Ncpus=2)"

WORKDIR /app

ENV BAYES_CRICKET_ROOT=/app
ENV MC_CORES=2

COPY stan/ ./stan/
COPY R/ ./R/
COPY data/ ./data/
COPY docker/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Run R under a pseudo-tty so logs stream in Docker Desktop / `docker logs` without `-t`
ENTRYPOINT ["/entrypoint.sh"]
