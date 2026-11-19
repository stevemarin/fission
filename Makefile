# Common commands for the fission project. Run `make` or `make help` for the list.
#
# Variables you can override on the command line:
#   make show NUC=U235
#   make scan NUC=U238 E=20.87 WIDTHS=4
#   make test ARGS="-k slbw -x"

UV      ?= uv
PY      := $(UV) run python
NUC     ?= U238
E       ?= 6.674
WIDTHS  ?= 6
ARGS    ?=
SRC     := src tests docs/figures

# nuclides and thermal scattering laws needed for the PWR pin cell
PIN_NUCLIDES := H1 O16 U234 U235 U238 Zr90 Zr91 Zr92 Zr94 Zr96 B10 B11
PIN_THERMAL  := HinH2O OinUO2 UinUO2
TEMPS        ?= 293.6 600 900 1200
T            ?= 293.6
WATER_TEMPS  ?= 293.6 600 900
CASES        ?= 1A 1B 1C 1D
SHIELD       ?= nr

# Independent checking tools (NJOY, OpenMC): prebuilt conda-forge packages in
# an isolated micromamba environment under tools/ (no sudo; separate from uv).
TOOLS        := tools
MAMBA_VER    := 2.10.0-0
MAMBA        := $(TOOLS)/bin/micromamba
MAMBA_ROOT   := $(TOOLS)/mamba
TOOLS_ENV    := $(TOOLS)/env
NJOY_VER     := 2016.78
OPENMC_VER   := 0.16.0
IN_TOOLS     := $(MAMBA) run -r $(MAMBA_ROOT) -p $(TOOLS_ENV)
MAMBA_URL    := https://github.com/mamba-org/micromamba-releases/releases/download/$(MAMBA_VER)/micromamba-linux-64

.DEFAULT_GOAL := help
.PHONY: help sync test test-slow test-rust test-fast lint format format-check typecheck check fix \
        tools njoy openmc tools-check clean-tools data tables doppler-tables groups pincell clean-tables show scan rm doppler secondary library docs figures clean

help:  ## list targets
	@grep -E '^[a-z-]+:.*## ' $(MAKEFILE_LIST) | \
		awk -F':.*## ' '{printf "  \033[36m%-13s\033[0m %s\n", $$1, $$2}'

# ---------------------------------------------------------------- environment
sync:  ## install/update the project environment from uv.lock
	$(UV) sync

# ---------------------------------------------------------------- checking tools
$(MAMBA):
	@mkdir -p $(TOOLS)/bin
	curl -fsSL -o $@.tmp $(MAMBA_URL)
	curl -fsSL -o $@.sha256 $(MAMBA_URL).sha256
	echo "$$(cat $@.sha256)  $@.tmp" | sha256sum -c -
	chmod +x $@.tmp
	mv $@.tmp $@

$(TOOLS_ENV): | $(MAMBA)
	$(MAMBA) create -y -q -r $(MAMBA_ROOT) -p $(TOOLS_ENV)

tools: njoy openmc tools-check  ## install NJOY and OpenMC into tools/ (independent checks only)

njoy: | $(TOOLS_ENV)  ## install NJOY2016 (conda-forge njoy2016) into tools/env
	$(MAMBA) install -y -q -r $(MAMBA_ROOT) -p $(TOOLS_ENV) -c conda-forge --override-channels njoy2016=$(NJOY_VER)

openmc: | $(TOOLS_ENV)  ## install OpenMC (executable + Python API, conda-forge) into tools/env
	$(MAMBA) install -y -q -r $(MAMBA_ROOT) -p $(TOOLS_ENV) -c conda-forge --override-channels openmc=$(OPENMC_VER) njoy2016=$(NJOY_VER)

tools-check:  ## show the NJOY and OpenMC versions installed in tools/env
	@printf 'stop\n' | $(IN_TOOLS) njoy | grep -i -m1 'njoy'
	@$(IN_TOOLS) openmc --version | head -2
	@$(IN_TOOLS) python -c "import openmc; print('OpenMC Python API', openmc.__version__)"

clean-tools:  ## remove tools/ (micromamba, NJOY, OpenMC)
	rm -rf tools

# ---------------------------------------------------------------- quality
test:  ## run all tests (ARGS="..." passes extra pytest options)
	$(UV) run pytest -q $(ARGS)

test-slow:  ## run only the slow physics checks (full resonance reconstructions)
	$(UV) run pytest -q -m slow $(ARGS)

test-rust:  ## Rust unit tests (cargo)
	cargo test --quiet

test-fast:  ## run tests, stop at first failure
	$(UV) run pytest -q -x $(ARGS)

lint:  ## ruff lint
	$(UV) run ruff check $(SRC)

format:  ## ruff format (rewrites files)
	$(UV) run ruff format $(SRC)

format-check:  ## check formatting without changing files
	$(UV) run ruff format --check $(SRC)

typecheck:  ## mypy
	$(UV) run mypy $(SRC)

fix:  ## auto-fix lint issues, then format
	$(UV) run ruff check --fix $(SRC)
	$(UV) run ruff format $(SRC)

check: lint format-check typecheck test test-rust  ## everything CI would run (not slow)

# ---------------------------------------------------------------- data
data:  ## fetch the ENDF/B-VIII.1 evaluations needed for the pin cell
	$(PY) -m fission.data $(PIN_NUCLIDES) --thermal $(PIN_THERMAL)

show:  ## summarize an evaluation's header and MF3 (NUC=U238)
	$(PY) -m fission.endf $(NUC)

scan:  ## SLBW scan across a resonance (NUC=U238 E=6.674 WIDTHS=6)
	$(PY) -m fission.endf.slbw $(NUC) $(E) --widths $(WIDTHS)

tables:  ## build (or load) cached pointwise tables for the pin-cell nuclides
	$(PY) -m fission.endf.cache $(PIN_NUCLIDES) --list

doppler-tables:  ## build (or load) Doppler-broadened tables (TEMPS="293.6 600 900 1200")
	$(PY) -m fission.processing.cache $(PIN_NUCLIDES) --temps $(TEMPS)

groups:  ## build (or load) multigroup libraries, XMAS-172, free gas + H in H2O
	$(PY) -m fission.processing.library $(PIN_NUCLIDES) --temps $(TEMPS) --quiet
	$(PY) -m fission.processing.library H1 --thermal HinH2O --temps $(WATER_TEMPS) --quiet

pincell:  ## k-infinity of the homogenized VERA problem 1 cell (CASES="1A" SHIELD=nr|none)
	$(PY) -m fission.pincell.run $(CASES) --shielding $(SHIELD)

clean-tables:  ## delete all cached pointwise, broadened and multigroup tables
	rm -rf endf-data/pointwise endf-data/doppler endf-data/groups

rm:  ## Reich-Moore vs SLBW across a resonance (NUC=U238 E=6.674 WIDTHS=6)
	$(PY) -m fission.endf.reich_moore $(NUC) $(E) --widths $(WIDTHS)

library:  ## multigroup summary of one nuclide (NUC=U235 T=293.6)
	$(PY) -m fission.processing.library $(NUC) --temps $(T)

secondary:  ## nu, delayed groups, fission spectrum, angular data (NUC=U235)
	$(PY) -m fission.endf.distributions $(NUC)

doppler:  ## Doppler-broadened peak and area of one resonance vs T (NUC=U238 E=6.674)
	$(PY) -m fission.processing.doppler $(NUC) $(E)

# ---------------------------------------------------------------- documentation
DOCS_MATH ?= --mathml
DOCS_SRC  := docs/theory/metadata.yaml $(sort $(wildcard docs/theory/[0-9]*.md))

docs:  ## build docs/build/theory.html (single file, offline; DOCS_MATH=--mathjax for MathJax)
	@mkdir -p docs/build
	pandoc $(DOCS_SRC) --standalone --embed-resources $(DOCS_MATH) --toc --toc-depth=2 \
		--number-sections --resource-path=docs/theory --css=docs/theory/style.css \
		-o docs/build/theory.html
	@echo "wrote docs/build/theory.html"

figures:  ## regenerate documentation figures from ENDF data
	$(PY) docs/figures/make_figures.py

# ---------------------------------------------------------------- housekeeping
clean:  ## remove caches and Rust build output (keeps endf-data/)
	rm -rf .pytest_cache .mypy_cache .ruff_cache target
	find . -name __pycache__ -type d -prune -not -path './.venv/*' -exec rm -rf {} +
