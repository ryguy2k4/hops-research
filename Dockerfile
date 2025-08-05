FROM python:3.12.8-slim

# Install system dependencies (add others if needed)
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt /tmp/
RUN pip install --upgrade pip && pip install -r /tmp/requirements.txt

# Create working directory
WORKDIR /project

# Copy everything into the container
COPY . /project

# Default command (overridable at runtime)
CMD ["snakemake", "--cores", "1"]
