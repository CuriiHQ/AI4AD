FROM debian:bullseye

USER root

#RUN apt-get update && apt-get install -y --no-install-recommends gnupg
#RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
#RUN echo "deb http://cloud.r-project.org/bin/linux/debian buster-cran35/" > /etc/apt/sources.list.d/cran.list

RUN apt-get update && \
    apt-get -y install r-base r-cran-rcpp r-cran-data.table libapt-pkg-dev less \
      wget \
      python3-pip \
      python3-numpy \
      python3-scipy \
      python3-matplotlib \
      python3-pandas && \
    apt-get clean

RUN R -e "install.packages('Matrix', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('glmnet', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('reticulate', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('foreach', repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('selectiveInference', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('ncvreg', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('ncpen', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('ggplot2', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('caret', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('mpath', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('stats', repos = 'http://cran.us.r-project.org')"

RUN pip3 install sklearn
RUN pip3 install --upgrade scipy
RUN pip3 install seaborn
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 2

RUN pip3 install notebook
RUN R -e "install.packages('IRkernel', repos = 'http://cran.us.r-project.org')"
RUN R -e "IRkernel::installspec(user = FALSE)"
