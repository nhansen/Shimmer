FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Copy content from host to container
RUN mkdir /shimmer
COPY . /shimmer

RUN apt-get update && apt-get install -y perl r-base samtools

RUN perl -MCPAN -e 'CPAN::Shell->install("Module::Build")'

# Install R statmod
RUN R -e "install.packages(c('statmod'),dependencies=TRUE, repos = 'http://cran.rstudio.com/')"

# Install Shimmer
RUN cd /shimmer && perl /shimmer/Build.PL --install_base /shimmer
RUN /shimmer/Build
RUN /shimmer/Build install

ENV PATH="${PATH}:/shimmer/bin"
ENV PERL5LIB="/shimmer/bin/lib/perl5:${PERL5LIB}"
