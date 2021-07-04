FROM ubuntu:20.04

WORKDIR /app

COPY environment.txt .
COPY CRIS_docker.sh .
COPY data ./data 

# Install miniconda
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
RUN conda install --file environment.txt

##backup
#sudo docker build --tag app:v1 . 
#sudo docker run -v $PWD/file.txt:/app/file.txt -t app:v1 bash main.sh file.txt

#Error: "ubuntu focal-backports InRelease: At least one invalid signature was encountered"
#docker system prune --force

#sudo docker run -v $PWD/SRR1814049_test.bam:/app/SRR1814049_test.bam -t app:v1 bash main.sh SRR1814049_test.bam

