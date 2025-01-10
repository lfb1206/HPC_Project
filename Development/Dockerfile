FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install development tools, Open MPI, Valgrind, and additional dependencies
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
        build-essential \
        openmpi-bin \
        libopenmpi-dev \
        valgrind \
        python3 \
        libpng-dev \
        libjpeg-dev \
        libgtk-3-dev \
        make && \
    rm -rf /var/lib/apt/lists/*

# Create a non-root user (e.g., "mpiuser")
RUN useradd -ms /bin/bash mpiuser

# Switch to this user
USER mpiuser

# Set working directory
WORKDIR /home/mpiuser

# Default command: Start a bash shell
CMD ["/bin/bash"]