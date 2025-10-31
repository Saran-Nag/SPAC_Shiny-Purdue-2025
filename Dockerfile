# Use official Python image as base
# FROM python:3.9.19-slim-bookworm
# Issue #15: Fix python version conflicts
# Error message when building docker:
# Package 'datashader' requires a different Python: 3.9.19 not in '>=3.10'
FROM python:3.10-slim-bookworm

# Set working directory
WORKDIR /app

# Issue #12: Fix 'Hash Sum Mismatch' bug for mac device
RUN echo "Acquire::http::Pipeline-Depth 0;" > /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::http::No-Cache true;" >> /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::BrokenProxy    true;" >> /etc/apt/apt.conf.d/99custom

# This is a previous method to fix this bug. It no longer works for 3.9.19-slim-bookworm
# RUN sed -i 's|http://deb.debian.org|https://deb.debian.org|g' /etc/apt/sources.list

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

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . .

# Expose port 8000 (default for Shiny)
EXPOSE 8000

# Run the Shiny app
CMD ["python", "-m", "shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "8000"]