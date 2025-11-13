
# Use official Python image as base
FROM python:3.9.19-slim-bookworm

# Set working directory
WORKDIR /app

# Fix 'Hash Sum Mismatch' bug for mac device
RUN echo "Acquire::http::Pipeline-Depth 0;" > /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::http::No-Cache true;" >> /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::BrokenProxy    true;" >> /etc/apt/apt.conf.d/99custom

# Install system dependencies needed for scientific packages
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    libhdf5-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements file
COPY requirements.txt .

RUN --mount=type=cache,target=/root/.cache/pip \
python -m pip install --upgrade pip setuptools wheel \
&& pip install -r requirements.txt