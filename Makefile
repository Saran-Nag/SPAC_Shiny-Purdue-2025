# SPAC Shiny Docker Commands
# Usage: make <command>

.PHONY: help build run stop clean logs shell restart

# Default target
help:
	@echo "SPAC Shiny Docker Management"
	@echo "============================"
	@echo ""
	@echo "Available commands:"
	@echo "  help      - Show this help message"
	@echo "  build     - Build the Docker image"
	@echo "  run       - Run the container (builds if needed)"
	@echo "  run-dev   - Run with volume mount for development"
	@echo "  stop      - Stop the running container"
	@echo "  restart   - Restart the container"
	@echo "  logs      - Show container logs"
	@echo "  shell     - Open a shell inside the container"
	@echo "  clean     - Remove container and image"
	@echo "  compose   - Use docker-compose (recommended)"
	@echo ""
	@echo "Quick start:"
	@echo "  make compose    # Start with docker-compose"
	@echo "  make logs       # View logs"
	@echo "  make stop       # Stop when done"
	@echo ""
	@echo "App will be available at: http://localhost:8001"

# Build the Docker image
build:
	@echo "Building SPAC Shiny Docker image..."
	docker build -t spac-shiny .

# Run the container (background)
run: build
	@echo "Starting SPAC Shiny container..."
	@echo "App will be available at: http://localhost:8001"
	docker run -d --name spac-shiny-app -p 8001:8000 spac-shiny

# Run with development volume mount
run-dev: build
	@echo "Starting SPAC Shiny container in development mode..."
	@echo "App will be available at: http://localhost:8001"
	docker run -d --name spac-shiny-app -p 8001:8000 -v $(PWD):/app spac-shiny

# Use docker-compose (recommended)
compose:
	@echo "Starting SPAC Shiny with docker-compose..."
	@echo "App will be available at: http://localhost:8001"
	docker-compose up -d

# Stop the container
stop:
	@echo "Stopping SPAC Shiny container..."
	-docker stop spac-shiny-app
	-docker rm spac-shiny-app
	-docker-compose down

# Restart the container
restart: stop run

# Show container logs
logs:
	@echo "Showing SPAC Shiny logs..."
	@echo "Press Ctrl+C to exit log view"
	-docker logs -f spac-shiny-app || docker-compose logs -f

# Open a shell inside the running container
shell:
	@echo "Opening shell in SPAC Shiny container..."
	docker exec -it spac-shiny-app /bin/bash

# Clean up everything
clean: stop
	@echo "Cleaning up SPAC Shiny Docker resources..."
	-docker rmi spac-shiny
	-docker system prune -f

# Check if container is running
status:
	@echo "SPAC Shiny Container Status:"
	@echo "=========================="
	@docker ps | grep spac-shiny || echo "No SPAC Shiny containers running"
	@echo ""
	@docker images | grep spac-shiny || echo "No SPAC Shiny images found"