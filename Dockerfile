# Use the Ubuntu rolling release as base image
FROM ubuntu:rolling

# Install necessary packages for C++ development
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    cmake \
    gdb \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /app

# Copy your C++ source code into the container (replace 'src' with your source code directory)
COPY . /app/

# Run any necessary build commands here (e.g., cmake and make)
# Example:
# RUN cmake /app/src
# RUN make

# Set the default command to run your application (if applicable)
# CMD ["/app/src/your_executable"]

# You can also leave the CMD line commented out and run your executable manually after starting the container

# Example usage:
# 1. Build the Docker image:
#    docker build -t my_cpp_app .
# 2. Run the Docker container:
#    docker run -it my_cpp_app
