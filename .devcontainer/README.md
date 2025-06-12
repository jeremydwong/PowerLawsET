# Development Container Setup

This directory contains the configuration for GitHub Codespaces and VS Code Development Containers.

## Quick Start with GitHub Codespaces

1. Go to your GitHub repository
2. Click the green "Code" button
3. Select "Codespaces" tab
4. Click "Create codespace on main"

The environment will automatically set up with:
- Julia 1.10.0
- All project dependencies
- Pluto notebook server
- VS Code extensions for Julia development

## Available Ports

- **1234**: Pluto notebook server
- **8888**: Jupyter server (if needed)
- **9999**: GLMakie server

## Quick Commands

Once your codespace is running:

```bash
# Start Pluto notebook server
pluto

# Start Julia with project environment
jlp

# Run the main analysis
julia --project=. -e "using PowerLawsEt; f, mainseq, fit = replicate_pub_whk_sim()"
```

## Development Workflow

1. **Notebooks**: Use Pluto notebooks in `src/notebooks/`
2. **Code**: Edit Julia source files in `src/`
3. **Data**: Access data files in `data/reaching/`
4. **Testing**: Run tests with `julia --project=. -e "using Pkg; Pkg.test()"`

## VS Code Extensions Included

- Julia Language Support
- Jupyter Notebooks
- GitHub Copilot
- Python (for any mixed workflows)

## Local Development Container

If you prefer to use VS Code locally with Docker:

1. Install Docker Desktop
2. Install "Dev Containers" extension in VS Code
3. Open this repository in VS Code
4. Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on Mac)
5. Select "Dev Containers: Reopen in Container"