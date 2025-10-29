# Use official Python image as base
FROM python:3.9.19-slim-bookworm

# Set working directory
WORKDIR /app

# Install system dependencies needed for scientific packages
RUN sed -i 's|http://deb.debian.org|https://deb.debian.org|g' /etc/apt/sources.list \
    && apt-get update && apt-get install -y \
# RUN apt-get update && apt-get install -y \
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