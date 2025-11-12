# Configuration
IMAGE_NAME = spac-shiny
BASE_IMAGE_NAME = spac-shiny-base
VERSION = latest
PORT = 8000

# Build the base image (run once or when requirements.txt changes)
build-base:
	@echo "Building base image with all dependencies..."
	@echo "This will take some time"
	DOCKER_BUILDKIT=1 docker build -f base.Dockerfile -t $(BASE_IMAGE_NAME):$(VERSION) .
	@echo "Base image built successfully!"

# Build the application image (fast - only copies your code)
build:
	@echo "Building application image..."
	DOCKER_BUILDKIT=1 docker build -t $(IMAGE_NAME):$(VERSION) .
	@echo "Application image built success"

# Run the application
run: build
	@echo "Starting application on port $(PORT)..."
	docker run -p $(PORT):$(PORT) $(IMAGE_NAME):$(VERSION)

# Complete rebuild (when requirements.txt changes)
rebuild: build-base build
	@echo "Complete rebuild finished!"

# First-time setup
setup: build-base build
	@echo "Setup complete!"
	@echo "Run 'make run' to start the application"

# Clean up application image only (keeps base)
clean:
	@echo "Cleaning up application image..."
	docker rmi $(IMAGE_NAME):$(VERSION) 2>/dev/null || true
	@echo "Done! Base image preserved."

# Deep clean (removes everything including base image)
clean-all:
	@echo "Cleaning up all images..."
	docker rmi $(IMAGE_NAME):$(VERSION) 2>/dev/null || true
	docker rmi $(BASE_IMAGE_NAME):$(VERSION) 2>/dev/null || true
	@echo "All images removed."

# Show image information
info:
	@echo "Docker images:"
	@docker images | grep -E "$(BASE_IMAGE_NAME)|$(IMAGE_NAME)|REPOSITORY" || echo "No images found"

# Help command
help:
	@echo "Available commands:"
	@echo "  make setup       - First-time setup (builds base + app)"
	@echo "  make build-base  - Build base image with dependencies"
	@echo "  make build       - Build application image"
	@echo "  make run         - Build and run the application"
	@echo "  make rebuild     - Rebuild everything (when requirements change)"
	@echo "  make clean       - Remove app image (keep base)"
	@echo "  make clean-all   - Remove all images"
	@echo "  make info        - Show image information"

.PHONY: build-base build run rebuild setup clean clean-all info help