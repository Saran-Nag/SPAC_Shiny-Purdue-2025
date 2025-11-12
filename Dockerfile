# Use official Python image as base
FROM python:3.9.13-slim

# Set working directory
WORKDIR /app

# Install system dependencies needed for scientific packages
# Use BuildKit cache mount to speed up apt operations
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    sed -i 's|http://deb.debian.org|https://deb.debian.org|g' /etc/apt/sources.list \
    && apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    libhdf5-dev \
    pkg-config

# Copy requirements file ONLY (not the entire app yet)
COPY requirements.txt .

# Upgrade pip and install dependencies with pip cache mount
# This dramatically speeds up reinstalls even when requirements change slightly
RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install --upgrade pip setuptools wheel \
    && pip install -r requirements.txt

# Copy the application code (do this LAST so code changes don't rebuild dependencies)
COPY . .

# Expose port 8000 (default for Shiny)
EXPOSE 8000

# Run the Shiny app
CMD ["python", "-m", "shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "8000"]