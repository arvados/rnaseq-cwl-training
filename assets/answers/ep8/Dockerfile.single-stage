FROM debian:10-slim
MAINTAINER Peter Amstutz <peter.amstutz@curii.com>

RUN apt-get update -qy
RUN apt-get install -qy build-essential wget unzip zlib1g-dev

# Install BWA 07.7.17
RUN wget https://github.com/lh3/bwa/archive/v0.7.17.zip && \
	unzip v0.7.17 && \
	cd bwa-0.7.17 && \
	make && \
	cp bwa /usr/bin && \
	cd .. && \
	rm -rf bwa-0.7.17
