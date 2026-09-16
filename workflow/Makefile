# Phase 1 pipeline. Run from the repository root:  make -f workflow/Makefile all
PY      ?= python
CONFIGS := configs/baseline_breach1_lgp.yaml configs/baseline_breach2_alsek.yaml
DAMS    := data/processed/dam_breach1_lgp.npz data/processed/dam_breach2_alsek.npz

.PHONY: all preprocess run test sensitivity fcurve clean

all: run test

preprocess: $(DAMS)

data/processed/dam_%.npz: configs/baseline_%.yaml
	$(PY) -m gpbreach.cli preprocess $<

run: preprocess
	@for c in $(CONFIGS); do $(PY) -m gpbreach.cli run $$c; done

sensitivity: preprocess
	@for c in $(CONFIGS); do $(PY) -m gpbreach.cli sensitivity $$c; done

fcurve: preprocess
	@for c in $(CONFIGS); do $(PY) -m gpbreach.cli fcurve $$c; done

# Prefers pytest; falls back to the no-pytest runner (see environment.yml).
test: preprocess
	@if $(PY) -c "import pytest" 2>/dev/null; then $(PY) -m pytest; \
	else echo "pytest not installed; using workflow/run_tests_without_pytest.py"; \
	     $(PY) workflow/run_tests_without_pytest.py; fi

clean:
	rm -rf results/scratch/* data/interim/*
