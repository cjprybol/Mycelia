# Mycelia CLI

A command line interface for the Mycelia genome graph construction and manipulation package.

## Installation

### Prerequisites

- Current Julia LTS or higher
- Mycelia package installed (`]add Mycelia`)

### Setup

1. Clone this repository or download the files:
   ```
   git clone https://github.com/cjprybol/Mycelia.git
   cd Mycelia
   ```

2. Initialize the project:
   ```
   julia --project -e "import Pkg; Pkg.instantiate()"
   ```

Setup the project (run from project root):
```
julia setup.jl install
```

Build both system image and app (run from project root):
```bash
julia setup.jl build
```

Or to build just one version:
```bash
julia bin/build_mycelia.jl sysimage
# or
julia bin/build_mycelia.jl app
```

Run the CLI tool (from project root):
```bash
# Using the system image
bin/mycelia command --args

# Using the standalone app
bin/mycelia-app command --args
```