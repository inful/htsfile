SHELL := /usr/bin/env bash
GO    ?= go

PACKAGES := ./...

.PHONY: all
all: fmt vet test lint ## fmt, vet, test, lint

.PHONY: fmt
fmt: ## run gofmt
	$(GO) fmt $(PACKAGES)

.PHONY: vet
vet: ## run go vet
	$(GO) vet $(PACKAGES)

.PHONY: test
test: ## run tests (with race detector)
	$(GO) test -race -count=1 $(PACKAGES)

.PHONY: test-short
test-short: ## run tests in short mode
	$(GO) test -race -count=1 -short $(PACKAGES)

.PHONY: cover
cover: ## run tests with coverage
	$(GO) test -race -count=1 -coverprofile=coverage.txt -covermode=atomic $(PACKAGES)
	@$(GO) tool cover -func=coverage.txt | tail -n 1

.PHONY: cover-html
cover-html: cover ## produce coverage.html
	@$(GO) tool cover -html=coverage.txt -o coverage.html

.PHONY: bench
bench: ## run benchmarks
	$(GO) test -bench=. -benchmem -run=^$$ $(PACKAGES)

.PHONY: lint
lint: ## run golangci-lint
	@command -v golangci-lint >/dev/null 2>&1 || { echo "golangci-lint not installed"; exit 1; }
	golangci-lint run $(PACKAGES)

.PHONY: tidy
tidy: ## go mod tidy
	$(GO) mod tidy

.PHONY: help
help: ## show this help
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)