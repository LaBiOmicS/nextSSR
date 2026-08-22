# Reproducible Docker Container for nextSSR
FROM python:3.11-slim

LABEL maintainer="Fabiano Menegidio"
LABEL description="nextSSR: High-performance, FAIR-compliant SSR identification"
LABEL org.opencontainers.image.source="https://github.com/fabianomenegidio/nextSSR"

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install nextSSR package
COPY pyproject.toml README.md ./
COPY nextssr ./nextssr

RUN pip install --no-cache-dir .

ENTRYPOINT ["nextssr"]
CMD ["--help"]
