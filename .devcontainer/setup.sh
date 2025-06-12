#!/bin/bash

set -e

echo "🚀 Setting up PowerLawsEt development environment..."

# Update package lists
sudo apt-get update

# Install essential build tools and libraries
sudo apt-get install -y \
    build-essential \
    git \
    curl \
    wget \
    ca-certificates \
    gnupg \
    lsb-release \
    pkg-config \
    cmake \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    libblas-dev \
    libhdf5-dev \
    hdf5-tools \
    libgl1-mesa-glx \
    libxrandr2 \
    libxinerama1 \
    libxcursor1 \
    libxi6 \
    libglu1-mesa \
    freeglut3-dev \
    mesa-common-dev \
    libglfw3-dev \
    libglew-dev

# Install Julia
echo "📦 Installing Julia..."
JULIA_VERSION="1.10.0"
JULIA_ARCH="linux-x86_64"
JULIA_URL="https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-${JULIA_VERSION}-${JULIA_ARCH}.tar.gz"

cd /tmp
wget -q $JULIA_URL
tar xzf julia-${JULIA_VERSION}-${JULIA_ARCH}.tar.gz
sudo mv julia-${JULIA_VERSION} /usr/local/julia
sudo ln -sf /usr/local/julia/bin/julia /usr/local/bin/julia

# Verify Julia installation
echo "✅ Julia version:"
julia --version

# Set up Julia environment
echo "🔧 Setting up Julia environment..."
cd /workspaces/PowerLawsEt

# Activate project and install dependencies
julia --project=. -e '
    using Pkg
    Pkg.instantiate()
    Pkg.precompile()
    
    # Install additional development packages
    Pkg.add([
        "Revise",
        "Pluto", 
        "PlutoUI",
        "BenchmarkTools",
        "ProfileView",
        "JuliaFormatter"
    ])
    
    # Precompile everything
    Pkg.precompile()
'

# Set up git configuration for the container
git config --global init.defaultBranch main
git config --global pull.rebase false
git config --global --add safe.directory /workspaces/PowerLawsEt

# Create useful aliases
cat >> ~/.bashrc << 'EOF'

# Julia aliases
alias jl='julia'
alias jlp='julia --project=.'
alias pluto='julia -e "using Pluto; Pluto.run(host=\"0.0.0.0\", port=1234, launch_browser=false)"'

# Git aliases  
alias gs='git status'
alias ga='git add'
alias gc='git commit'
alias gp='git push'
alias gl='git log --oneline'

echo "🎯 PowerLawsEt development environment ready!"
echo "💡 Use 'pluto' to start Pluto server on port 1234"
echo "💡 Use 'jlp' to start Julia with project environment"
EOF

echo "🎉 Setup complete! PowerLawsEt development environment is ready."
echo ""
echo "Quick start commands:"
echo "  pluto           - Start Pluto notebook server"
echo "  jlp             - Start Julia with project environment"
echo "  julia --project=. -e 'using PowerLawsEt; f, mainseq, fit = replicate_pub_whk_sim()'"