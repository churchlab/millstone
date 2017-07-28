FROM ubuntu:16.04

RUN apt-get update && apt-get install -y --no-install-recommends \
          build-essential \
          ca-certificates \
          curl \
          git \
          libpq-dev \
          postgresql \
          python \
          python-dev \
          ssh \
          && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN curl -O https://bootstrap.pypa.io/get-pip.py && \
    python get-pip.py && \
    rm get-pip.py

RUN git clone https://github.com/churchlab/millstone.git && \

# Run tests. Used while setting up Docker.
RUN cd millstone && \
    pip install -r requirements/deploy.txt && \
    cd genome_designer && \
    ./scripts/run_unit_tests.sh
