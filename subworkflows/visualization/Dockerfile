# Interactive visualization image (sp_visualizer.py).
# Docker equivalent of visu.def -> published to ghcr.io/cbib/deconvolista-visu by
# .github/workflows/build-images.yml. The script itself is mounted at runtime
# (Singularity --bind), so the image only needs to ship the Python dependencies.
FROM python:3.11-slim

# Auto-link the published package to the repo on GHCR (provenance + visibility).
LABEL org.opencontainers.image.source=https://github.com/cbib/DeconvoliSTa

RUN pip install --no-cache-dir pandas bokeh pillow numpy

# Sanity check at build time (fails the build if a dep is missing).
RUN python -c "import pandas, bokeh, PIL, numpy; print('visu deps OK | bokeh', bokeh.__version__)"
