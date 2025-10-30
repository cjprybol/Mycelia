# Contributing to Mycelia

Thank you for your interest in contributing to Mycelia! This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally
3. Create a new branch for your feature or fix
4. Make your changes following our coding standards

## Development Setup

```bash
## Clone your fork
git clone https://github.com/your-username/Mycelia.jl.git
cd Mycelia.jl

## Install dependencies
julia --project=. -e "using Pkg; Pkg.instantiate()"

## Run tests
julia --project=. -e "using Pkg; Pkg.test()"
```

## Code Standards

### Style Guidelines

- Follow Julia style conventions
- Use full module namespacing (e.g., `Test.@test`, not `@test`)
- Import packages at module level, not in individual files
- Add docstrings to all exported functions

### Testing

- Write tests for all new functionality
- Ensure all tests pass before submitting PR
- Include edge cases and error conditions
- Follow existing test organization structure

### Documentation

- Add docstrings using DocStringExtensions format
- Update relevant tutorial files
- Add examples for new features
- Keep documentation measured and accurate

## Submitting Changes

1. Ensure all tests pass
2. Update documentation as needed
3. Commit with clear, descriptive messages
4. Push to your fork
5. Submit a pull request

## Pull Request Guidelines

- Reference related issues
- Describe the changes made
- Include test results
- Be responsive to review feedback

## Questions?

- Open an issue for bugs or feature requests
- Start a discussion for general questions
- Check existing issues and documentation first

## Code of Conduct

Be respectful, inclusive, and professional in all interactions.

Thank you for contributing to Mycelia!
