FROM ubuntu:latest

WORKDIR /app

COPY mybash.sh .
COPY environment.txt .

# Install conda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

#install dependencies using conda
RUN conda create --name cris_env --file environment.txt

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "cris_env", "/bin/bash", "-c"]

CMD bash CRIS.sh