FROM ubuntu:18.04 as packages
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y \
        python-dev \
        python-pip \
        python-tk 

RUN pip install --upgrade distribute && pip install ipython==3.2.3 matplotlib==1.5.3 numpy==1.6.2 Pillow==3.4.0 pyparsing==2.4.0 python-dateutil==2.8.0 pytz==2019.1 scipy==0.17.1 six==1.12.0

