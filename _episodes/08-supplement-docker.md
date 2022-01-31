---
title: "Supplement: Creating Docker Images for Workflows"
teaching: 10
exercises: 1
questions:
- "How do I create Docker images from scratch?"
- "What some best practices for Docker images?"
objectives:
- "Understand how to get started writing Dockerfiles"
keypoints:
- "Docker images contain the initial state of the filesystem for a container"
- "Docker images are made up of layers"
- "Dockerfiles consist of a series of commands to install software into the container."
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

Use the `-t` option to specify the name of the image.  Use `-f` if the
file isn't named exactly `Dockerfile`.  The last part is the directory
where it will find the `Dockerfile` and any files that are referenced
by `COPY` (described below).

```
docker build -t training/bwa -f Dockerfile.single-stage .
```

> ## Exercise
>
> Create a `Dockerfile` based on this lesson and build it for yourself.
>
> > ## Solution
> >
> > FROM debian:10-slim
> > MAINTAINER Peter Amstutz <peter.amstutz@curii.com>
> >
> > RUN apt-get update -qy
> > RUN apt-get install -qy build-essential wget unzip zlib1g-dev
> >
> > # Install BWA 07.7.17
> > RUN wget https://github.com/lh3/bwa/archive/v0.7.17.zip && \
> > 	unzip v0.7.17 && \
> > 	cd bwa-0.7.17 && \
> > 	make && \
> > 	cp bwa /usr/bin && \
> > 	cd .. && \
> > 	rm -rf bwa-0.7.17
> >
> {: .solution}
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

# Multi-stage builds

As noted, it is good practice to avoiding leaving files in the Docker
image that were required to build the program, but not to run it, as
those files are simply useless bloat.  Docker offers a more
sophisticated way to create clean builds by separating the build steps
from the creation of the final container.  These are called
"multi-stage" builds.

A multi stage build has multiple `FROM` lines.  Each `FROM` line is a
separate container build.  The last `FROM` in the file describes the
final container image that will be created.

The key benefit is that the different stages are independent, but you
can copy files from one stage to another.

Here is an example of the bwa build as a multi-stage build.  It is a
little bit more complicated, but the outcome is a smaller image,
because the "build-essential" tools are not included in the final
image.

```
# Build the base image.  This is the starting point for both the build
# stage and the final stage.
# the "AS base" names the image within the Dockerfile
FROM debian:10-slim AS base
MAINTAINER Peter Amstutz <peter.amstutz@curii.com>

# Install libz, because the bwa binary will depend on it.
# As it happens, this already included in the base Debian distribution
# because lots of things use libz specifically, but it is good practice
# to explicitly declare that we need it.
RUN apt-get update -qy
RUN apt-get install -qy zlib1g


# This is the builder image.  It has the commands to install the
# prerequisites and then build the bwa binary.
FROM base as builder
RUN apt-get install -qy build-essential wget unzip zlib1g-dev

# Install BWA 07.7.17
RUN wget https://github.com/lh3/bwa/archive/v0.7.17.zip
RUN unzip v0.7.17
RUN cd bwa-0.7.17 && \
    make && \
    cp bwa /usr/bin


# Build the final image.  It starts from base (where we ensured that
# libz was installed) and then copies the bwa binary from the builder
# image.  The result is the final image only has the compiled bwa
# binary, but not the clutter from build-essentials or from compiling
# the program.
FROM base AS final

# This is the key command, we use the COPY command described earlier,
# but instead of copying from the host, the --from option copies from
# the builder image.
COPY --from=builder /usr/bin/bwa /usr/bin/bwa
```

# Best practices for Docker images

[Docker has published guidelines on building efficient images.](https://docs.docker.com/develop/develop-images/dockerfile_best-practices/)

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

Similarly, be as specific as possible about the version of the base
image you want to use in your `FROM` command.  If you don't specify a
tag, the default tag is called "latest", which can change at any time.

## Tag your builds

Use meaningful tags on your own Docker image so you can tell versions
of your Docker image apart as it is updated over time.  These can
reflect the version of the underlying software, or a version you
assign to the Dockerfile itself.  These can be manually assigned
version numbers (e.g. 1.0, 1.1, 1.2, 2.0), timestamps (e.g. YYYYMMDD
like 20220126) or the hash of a git commit.

## Avoid putting reference data to Docker images

Bioinformatics tools often require large reference data sets to run.
These should be supplied externally (as workflow inputs) rather than
added to the container image.  This makes it easy to update reference
data instead of having to rebuild a new Docker image every time, which
is much more time consuming.

## Small scripts can be inputs, too

If you have a small script, e.g. a self-contained single-file Python
script which imports Python modules installed inside the container,
you can supply the script as a workflow input.  This makes it easy to
update the script instead of having to rebuild a new Docker image
every time, which is much more time consuming.

## Don't use ENTRYPOINT

The `ENTRYPOINT` Dockerfile command modifies the command line that is
executed inside the container.  This can produce confusion when the
command line that supplied to the container and the command that
actually runs are different.

## Be careful about the build cache

Docker build has a useful feature where if it has a record of the
exact `RUN` command against the exact base layer, it can re-use the
layer from cache instead of re-running it every time.  This is a great
time-saver during development, but can also be a source of
frustration: build steps often download files from the Internet.  If
the file being downloaded changes without the command being used to
download it changing, it will reuse the cached step with the old copy
of the file, instead of re-downloading it.  If this happens, use
`--no-cache` to force it to re-run the steps.

> ## Episode solution
> * <a href="../assets/answers/ep8/Dockerfile.single-stage">Dockerfile.single-stage</a>
> * <a href="../assets/answers/ep8/Dockerfile.multi-stage">Dockerfile.multi-stage</a>
{: .solution}
