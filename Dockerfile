FROM kbase/sdkbase2:python
MAINTAINER Sean Jungbluth <sjungbluth@lbl.gov>
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.


# To install all the dependencies
RUN apt-get update && apt-get install -y libgsl0-dev samtools git zip unzip bedtools bowtie2 wget python-pip && \
    apt-get install -y r-base r-cran-gplots

# Need a recent version of numpy to work with CONCOCT
RUN pip install --upgrade --no-deps --force-reinstall numpy

# To download the CONCOCT software from Github and install it and its requirements
RUN git clone https://github.com/BinPro/CONCOCT && cd /CONCOCT && \
    pip install -r requirements.txt && \
    python setup.py install && \
    cd .. && \
    cp -R CONCOCT /kb/deployment/bin/CONCOCT

WORKDIR /kb/module/lib/kb_concoct/bin/

RUN wget https://sourceforge.net/projects/bbmap/files/latest/download && tar -xvf download

RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2 && tar -xvf minimap2-* && cd minimap2* && make && cd ../ && rm minimap2-2.17.tar.bz2

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && unzip hisat2-* && rm hisat2-2.1.0-Linux_x86_64.zip

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

ENV PATH=/kb/module/lib/kb_concoct/bin:$PATH
ENV PATH=/kb/module/lib/kb_concoct/bin/bbmap:$PATH
ENV PATH=/kb/module/lib/kb_concoct/bin/minimap2-2.17/:$PATH
ENV PATH=/kb/module/lib/kb_concoct/bin/hisat2-2.1.0/:$PATH
ENV PATH=/kb/deployment/bin/CONCOCT/bin:/kb/deployment/bin/CONCOCT/scripts:$PATH

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
