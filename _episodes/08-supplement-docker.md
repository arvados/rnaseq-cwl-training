---
title: "Supplement: Creating Docker Images"
teaching: 10
exercises: 1
questions:
- "How do I create Docker images from scratch?"
- "What some best practices for Docker images?"
objectives:
- ""
keypoints:
- ""
---

Common Workflow Language supports running tasks inside software
containers.  Software container systems (such as Docker) create an
execution environment that is isolated from the host system, so that
software installed on the host system does not conflict with the
software installed inside the container.

Programs running inside a software container get a different (and
generally restricted) view of the system than processes running
outside the container.  One of the most important and useful features
is that the containerized program has a different view of the file
system.  A program running inside a container, searching for
libraries, modules, configuration files, data files, etc, only sees
the files defined inside the container.

This means that, usually, a given file _path_ refers to _different
actual files_ depending from the persective of being inside or outside
the container.  It is also possible to have a file from the host
system appear at some location inside the container, meaning that the
_same file_ appears at _different paths_ depending from the persective
of being inside or outside the container.

The complexity of translating between the container and its host
environment is handled by the Common Workflow Language runner.  As a
workflow author, you only need to worry about the environment _inside_
the container.

# What are Docker images?

The Docker image describes the starting conditions for the container.
Most importantly, this includes starting layout and contents of the
container's file system.  This file system is typically a lightweight
POSIX environment, providing a standard set of POSIX utilities like a
`sh`, `ls`, `cat`, etc and organized into standard POSIX directories
like `/bin` and `/lib`.

The image is is made up of multiple "layers".  Each layer modifies the
layer below it by adding, removing or modifying files to produce a new
layer.  This allows lower layers to be re-used.

# Writing a Dockerfile

In this example, we will build a Docker image containing the
Burrows-Wheeler Aligner (BWA) by Heng Li.  This is just for
demonstration, in practice you should prefer to use existing
containers from [BioContainers](https://biocontainers.pro/), which
includes `bwa`.

Each line of the Docker file consists of a COMMAND in all caps,
following by the parameters of that command.

The first line of the file will specify the base image that we are
going to build from.  As mentioned, images are divided up into
"layers", so this tells Docker what to use for the first layer.


```
FROM debian:10-slim
```

This starts from the lightweight ("slim") Debian 10 Docker image.

Docker images have a special naming scheme.

A bare name like "debian" or "ubuntu" means it is an official Docker
image.  It has an implied prefix of "library", so you may see the
image referred to as "library/debian".  Official images are published
on [Docker Hub](https://hub.docker.com/search?type=image&image_filter=official).

A name with two parts separated by a slash is published on Docker Hub
by someone else.  For example, `amazon/aws-cli` is published by
Amazon.  These can also be found on [Docker Hub](https://hub.docker.com/search?type=image).

A name with three parts separated by slashes means it is published on
a different container register.  For example,
`quay.io/biocontainers/subread` is published by `quay.io`.

Following image name, separated by a colon is the "tag".  This is
typically the version of the image.  If not provided, the default tag
is "latest".  In this example, the tag is "10-slim" indicating Debian
release 10.

> ## Best practice
>
> You should always include the tag to refer to a specific image
> version, or you might run into problems when "latest" changes.

The Docker file should also include a MAINTAINER (this is purely
metadata, it is stored in the image but not used for execution).

```
MAINTAINER Peter Amstutz <peter.amstutz@curii.com>
```

Next is the default user inside the image.  By making choosing root,
we can change anything inside the image (but not outside).

The body of the Dockerfile is a series of `RUN` commands.

Each command is run with `/bin/sh` inside the Docker container.

Each `RUN` command creates a new layer.

The `RUN` command can span multiple lines by using a trailing
backslash.

For the first command, we use `apt-get` to install some packages that
will be needed to compile `bwa`.  The `build-essential` package
installs `gcc`, `make`, etc.

```
RUN apt-get update -qy && \
	apt-get install -qy build-essential wget unzip
```

Now we do everything else: download the source code of bwa, unzip it,
make it, copy the resulting binary to `/usr/bin`, and clean up.

```
# Install BWA 07.7.17
RUN wget https://github.com/lh3/bwa/archive/v0.7.17.zip && \
	unzip v0.7.17 && \
	cd bwa-0.7.17 && \
	make && \
	cp bwa /usr/bin && \
	cd .. && \
	rm -rf bwa-0.7.17
```

Because each `RUN` command creates a new layer, having the build and
clean up in separate `RUN` commands would mean creating a layer that
includes the intermediate object files from the build.  These would
then be carried around as part of the container image forever, despite
being useless.  By doing the entire build and clean up in one `RUN`
command, only the final state of the file system, with the binary
copied to `/usr/bin`, is committed to a layer.

To build a Docker image from a Dockerfile, use `docker build`.

This command takes the name to use for the image with `-t`, and the
directory that it should find the `Dockerfile`:

```
docker build -t training/bwa .
```

> ## Exercise
>
> Create a `Dockerfile` based on this lesson and build it for yourself.
>
{: .challenge}

# Adding files to the image during the build

Using the `COPY` command, you can copy files from the source directory
(this is the directory your Dockerfile was located) into the image
during the build.  For example, you have a `requirements.txt` next to
Dockerfile:

```
COPY requirements.txt /tmp/
RUN pip install --requirement /tmp/requirements.txt
```

# Best practices for Docker images

Docker has published guidelines on building efficient images:

https://docs.docker.com/develop/develop-images/dockerfile_best-practices/

Some additional considerations when building images for use with Workflows:

## Store Dockerfiles in git, alongside workflow definitions

Dockerfiles are scripts and should be managed with version control
just like other kinds of code.

## Be specific about software versions

Instead of blindly installing the latest version of a package, or
checking out the `master` branch of a git repository and building from
that, be specific in your Dockerfile about what version of the
software you are installing.  This will greatly aid the
reproducibility of your Docker image builds.

## Tag your builds

Use meaningful tags on the Docker image so you can tell versions of
your Docker image apart as it is updated over time.  These can reflect
the version of the underlying software, or the version of the
Dockerfile itself.  These can be manually assigned version numbers
(e.g. 1.0, 1.1, 1.2, 2.0), timestamps (e.g. YYYYMMDD like 20220126) or
the hash of a git commit.

## Avoid putting reference data to Docker images

Bioinformatics tools often require large reference data sets to run.
These should be supplied externally (as workflow inputs) rather than
added to the container image.  This makes it easy to update reference
data instead of having to rebuild and re-upload a new Docker image
every time, which is much more time consuming.

## Small scripts can be inputs, too

If you have a small script, e.g. a self-contained Python script which
relies on modules installed inside the container, but is itself
contained in a single file, you can supply the script as a workflow
input.  This makes it easy to update the script instead of having to
rebuild and re-upload a new Docker image every time, which is much
more time consuming.

## Don't use ENTRYPOINT

The `ENTRYPOINT` Dockerfile command modifies the command line that is executed
inside the container.  This can result in confusion when the command
line that was supplied to the container and the command that actually
runs are different.
