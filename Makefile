SHELL := /usr/bin/env bash
GO    ?= go

PACKAGES := ./...

.PHONY: all
all: fmt-check vet test lint ## fmt-check, vet, test, lint

.PHONY: fmt
fmt: ## rewrite Go files in place with gofmt
	gofmt -w .

.PHONY: fmt-check
fmt-check: ## fail if any Go file needs gofmt (mirrors CI)
	@test -z "$$(gofmt -l .)" || { echo "Files need gofmt:"; gofmt -l .; exit 1; }

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

.PHONY: ci
ci: fmt-check vet test lint ## full local pre-push check (mirrors CI)

.PHONY: tidy
tidy: ## go mod tidy
	$(GO) mod tidy

.PHONY: install-hooks
install-hooks: ## install the repo's git hooks (one-time per clone)
	git config core.hooksPath .githooks

.PHONY: help
help: ## show this help
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)