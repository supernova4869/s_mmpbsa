## Documentation Overview

This is the official documentation for the s_mmpbsa project, providing detailed project introduction, installation guides, usage methods, API references, and more. The documentation is built using Sphinx and can generate multiple formats (HTML, PDF, EPUB, etc.).

## Documentation Structure

This documentation contains the following main sections:

- **Homepage**: Project overview and main feature introduction
- **Introduction**: Basic principles of the MM-PBSA method and advantages of s_mmpbsa
- **Installation**: System requirements and installation steps
- **Quick Start**: Basic usage流程 and examples
- **Usage Guide**: Detailed feature descriptions and parameter settings
- **API Documentation**: Developer reference documentation
- **FAQ**: Frequently asked questions and answers

## Building the Documentation

To build the documentation locally, you first need to install Python and Sphinx. Here are the steps to build the documentation:

### 1. Install Dependencies

First, install the Python dependencies required for building the documentation:

```bash
# Execute in the project root directory
pip install -r docs/requirements.txt
```

### 2. Build HTML Documentation

#### On Linux/Mac Systems:

```bash
cd docs
make html
```

Or use the provided simplified script:

```bash
cd docs
chmod +x build_docs.sh
./build_docs.sh
```

#### On Windows Systems:

```batch
cd docs
make.bat html
```

### 3. View Documentation

After building is complete, you can find the generated HTML documentation in the `docs/_build/html` directory. Open the `index.html` file with a browser to view the complete documentation.

## Building Documentation in Other Formats

In addition to HTML format, you can also build documentation in other formats:

### PDF Format

```bash
# Linux/Mac
cd docs
make latexpdf

# Windows
cd docs
make.bat latexpdf
```

Note: Building PDF documentation requires installing a LaTeX distribution (such as TeX Live, MiKTeX, etc.).

### EPUB Format

```bash
# Linux/Mac
cd docs
make epub

# Windows
cd docs
make.bat epub
```

## Viewing Documentation on ReadTheDocs

Project documentation is also hosted on the ReadTheDocs platform, and you can access it directly in your browser:

- Stable version: https://s-mmpbsa.readthedocs.io/en/stable/
- Development version: https://s-mmpbsa.readthedocs.io/en/latest/

## Contributing to Documentation

If you find errors in the documentation or have suggestions for improvement, we welcome your contributions:

1. Raise an Issue describing the problem or suggestion
2. Submit a Pull Request to fix the problem or add new content

The documentation is written in reStructuredText format, and the main content files are located in the `docs` directory.

## Contact Information

If you have any questions, please contact the project maintainer:

- Email: email@example.com
- QQ Group: 123456789