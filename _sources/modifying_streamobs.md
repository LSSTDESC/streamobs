# Contributing to StreamObs

Basic guidelines for contributors.

## GitHub workflow

### 1. Update your local `main` branch

```bash
git checkout main
git pull origin main
```

### 2. Create a new branch

```bash
git checkout -b MYBRANCH
```

Use a descriptive branch name (e.g. `feature/new-model` or `fix/orbit-bug`).

### 3. Make your changes

Commit your changes regularly:

```bash
git add .
git commit -m "Short description of the changes"
```

You may create multiple commits on the same branch.

### 4. Run the tests

Before pushing your changes, make sure all tests pass:

```bash
pytest
```

If you add new features or fix bugs, add corresponding tests.

### 5. Push your branch

```bash
git push origin MYBRANCH
```

### 6. Open a Pull Request
Once you did all the modifications you wanted, open a Pull Request on GitHub and request a review from one of the main contributors.

### 7. Update your local repository after merging

```bash
git checkout main
git pull origin main
```

## Testing

All changes merged into `main` must pass the test suite.

### Adding tests

* Add tests in the `tests/` directory.
* Place tests in the file corresponding to the modified module.
* Reuse fixtures defined in `conftest.py` whenever possible.
* Add regression tests for bug fixes.

### Test markers

Custom markers can be added using:

```python
@pytest.mark.mytag
```

Any new marker must also be declared in `pytest.ini`.

### Running tests

Run the full test suite:

```bash
pytest
```

Run a subset of tests:

```bash
pytest -k "mytag"
```