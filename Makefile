# SPAC Shiny Docker Commands
# Usage: make <command>

.PHONY: help build run run-dev stop clean logs shell restart deploy status

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
	@echo "  deploy    - Deploy to Posit Connect (appshare-dev)"
	@echo "  status    - Show container status"
	@echo ""
	@echo "Quick start:"
	@echo "  make run        # Build and run the container"
	@echo "  make logs       # View logs"
	@echo "  make stop       # Stop when done"
	@echo "  make deploy     # Deploy to Posit Connect"
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



# Stop the container
stop:
	@echo "Stopping SPAC Shiny container..."
	-docker stop spac-shiny-app
	-docker rm spac-shiny-app

# Restart the container
restart: stop run

# Show container logs
logs:
	@echo "Showing SPAC Shiny logs..."
	@echo "Press Ctrl+C to exit log view"
	docker logs -f spac-shiny-app

# Open a shell inside the running container
shell:
	@echo "Opening shell in SPAC Shiny container..."
	docker exec -it spac-shiny-app /bin/bash

# Clean up everything
clean: stop
	@echo "Cleaning up SPAC Shiny Docker resources..."
	-docker rmi spac-shiny
	-docker system prune -f

# Deploy to Posit Connect
deploy:
	@echo "Deploying SPAC Shiny to Posit Connect..."
	@echo "Getting current commit hash for version tracking..."
	$(eval COMMIT_HASH := $(shell git rev-parse --short HEAD))
	@echo "Deploying version: $(COMMIT_HASH)"
	@echo "Sourcing bash profile and activating conda environment..."
	bash -c "source ~/.bash_profile && cd $(PWD) && conda activate spac-shiny-3-9-16 && rsconnect deploy shiny -n appshare-dev -t 'SPAC - Spatial Analysis Dashboard ($(COMMIT_HASH))' -v -a 4d6cc93f-3829-4935-8987-8c169549dbff ."
	@echo ""
	@echo "Deployment completed! Access your app at:"
	@echo "Dashboard: https://appshare-dev.cancer.gov/connect/#/apps/4d6cc93f-3829-4935-8987-8c169549dbff/access"
	@echo "Direct URL: https://appshare-dev.cancer.gov/content/4d6cc93f-3829-4935-8987-8c169549dbff/"

# Check if container is running
status:
	@echo "SPAC Shiny Container Status:"
	@echo "=========================="
	@docker ps | grep spac-shiny || echo "No SPAC Shiny containers running"
	@echo ""
	@docker images | grep spac-shiny || echo "No SPAC Shiny images found"