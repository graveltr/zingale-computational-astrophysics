# Use the official GCC image from the Docker Hub
FROM gcc:latest

# Set the working directory inside the container
WORKDIR /usr/src/app

# Install any additional packages if necessary
RUN apt-get update && \
    apt-get install -y vim gdb cmake make

# Copy the current directory contents into the container at /usr/src/app
COPY . .
