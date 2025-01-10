SHELL := /bin/bash

BOLD ?= $(shell tput bold)
NORMAL ?= $(shell tput sgr0)

help:
	@echo "Available targets:"
	@echo "  make start   - Build and start Docker Compose, attach a shell"
	@echo "  make stop    - Stop containers"
	@echo "  make down    - Remove containers"

start:
	docker-compose up -d
	docker-compose exec mycontainer /bin/bash

stop:
	docker-compose stop

down:
	docker-compose down
