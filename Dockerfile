FROM --platform=linux/amd64 mambaorg/micromamba:2.0-bookworm-slim

COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y

# curl for the downloads below.
USER root
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# isotwas and the dependencies that have no conda package, at fixed versions.
COPY tools/install_r_source_packages.sh /opt/isoqtl_sim/tools/install_r_source_packages.sh
RUN export PATH=/opt/conda/bin:$PATH \
    && bash /opt/isoqtl_sim/tools/install_r_source_packages.sh fetch /tmp/r_src \
    && bash /opt/isoqtl_sim/tools/install_r_source_packages.sh install /tmp/r_src \
    && rm -rf /tmp/r_src

COPY R /opt/isoqtl_sim/R
COPY bin /opt/isoqtl_sim/bin
COPY tests /opt/isoqtl_sim/tests
COPY benchmark /opt/isoqtl_sim/benchmark
COPY config /opt/isoqtl_sim/config

# twas_sim at a fixed commit, with two documented fixes (see benchmark/fetch_twas_sim.sh).
RUN bash /opt/isoqtl_sim/benchmark/fetch_twas_sim.sh /opt/twas_sim \
    && chmod +x /opt/isoqtl_sim/bin/isoqtl-sim

ENV PATH=/opt/conda/bin:/opt/isoqtl_sim/bin:$PATH
ENV TWAS_SIM_DIR=/opt/twas_sim
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["isoqtl-sim", "--help"]
