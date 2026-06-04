# StreamObs Documentation

This directory contains the source files for the StreamObs documentation, built using [Sphinx](https://www.sphinx-doc.org/).

## Documentation Structure

```
docs/
├── source/           # Documentation source files
│   ├── conf.py      # Sphinx configuration
│   ├── index.rst    # Main documentation page
│   ├── *.md         # Content pages (Markdown)
│   └── *.rst        # API reference (reStructuredText)
├── build/           # Built HTML files (generated, not tracked)
├── Makefile         # Build commands for Linux/macOS
├── make.bat         # Build commands for Windows
└── deploy_to_gh_pages.sh  # Deployment script
```

## Dependencies

To build the documentation, you need the following Python packages:

```bash
pip install sphinx sphinx-book-theme numpydoc myst-parser
```

Or install all documentation dependencies:

```bash
# From the repository root
pip install -r docs/requirements.txt
```

### Required Packages:
- **Sphinx** (>=8.0): Documentation generator
- **sphinx-book-theme**: Clean, modern documentation theme
- **numpydoc**: NumPy-style docstring support
- **myst-parser**: Markdown support in Sphinx

## Building the Documentation

### Local Build

To build the documentation locally for preview:

```bash
cd docs
make html
```

The built documentation will be in `docs/build/html/`. Open `docs/build/html/index.html` in your browser to view it.

### Clean Build

To remove all built files and rebuild from scratch:

```bash
cd docs
make clean
make html
```

### Common Build Issues

- **Import errors**: Make sure the `streamobs` package is in your `PYTHONPATH`
- **Missing dependencies**: Install all required packages listed above
- **Warnings about missing files**: Check that all files referenced in `index.rst` exist

## Editing Documentation

### Content Pages (Markdown)

Content pages are written in Markdown with MyST syntax:

- `about.md` - Package overview and information
- `quickstart.md` - Quick start guide
- `dependencies.md` - Installation instructions
- `citation.md` - Citation information

To edit, simply modify the `.md` files in `docs/source/`.

### API Reference (Auto-generated)

The API reference is automatically generated from docstrings in the Python code:

- Docstrings should follow the [NumPy documentation style](https://numpydoc.readthedocs.io/)
- Edit docstrings in the `.py` files in `streamobs/`
- The API pages (e.g., `streamobs.surveys.rst`) are auto-generated

### Main Index Page

The main landing page is `docs/source/index.rst` written in reStructuredText.

## Deploying Documentation

The documentation is hosted on GitHub Pages at: https://lsstdesc.github.io/streamobs/

### Manual Deployment

If you need to manually deploy documentation:

```bash
# Make sure you're on the main branch with all changes committed
cd docs
./deploy_to_gh_pages.sh
```

This script will:
1. Check you're on the correct branch
2. Build the documentation
3. Switch to the `gh-pages` branch
4. Copy built files to the root
5. Commit and push to GitHub

**Note**: The deployment script currently requires you to be on the `main` branch.

## Workflow Summary

1. Edit the relevant `.md` or `.rst` files in `docs/source/`
2. Build locally to preview: `make html`
3. Commit and push changes to your branch
4. Create a pull request to `main`
5. Once merged, documentation deploy manually the documentation by running
```bash
# Make sure you're on the main branch with all changes committed
cd docs
./deploy_to_gh_pages.sh
```

## Additional Resources

- [Sphinx Documentation](https://www.sphinx-doc.org/)
- [MyST Parser Guide](https://myst-parser.readthedocs.io/)
- [NumPy Docstring Guide](https://numpydoc.readthedocs.io/en/latest/format.html)
- [reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
