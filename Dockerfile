# Use official Python image as base
# Fix: python version conflicts (Issue #15 in Purdue-2025 fork)
# Error message when building docker:
# Package 'datashader' requires a different Python: 3.9.19 not in '>=3.10'
# FROM python:3.9.19-slim-bookworm
FROM python:3.10-slim-bookworm

# Set working directory
WORKDIR /app

# Fix: 'Hash Sum Mismatch' bug for mac device (Issue #12 in Purdue-2025 fork)
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

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . .

# Expose port 8000 (default for Shiny)
EXPOSE 8000

# Run the Shiny app
CMD ["python", "-m", "shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "8000"]