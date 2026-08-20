# Development setup

The following steps can be used to set up a clean virtual environment that is configured to automatically rebuild the library before the test suite is run, allowing for fast development iteration.

## 1. Install HDF4

Install HDF4, which is needed to run the PACE test suite.

  ```bash
  sudo apt-get install -y libhdf4-dev 
  ```

## 2. Clone and pull large files from Git LFS

```bash
git clone https://github.com/Light-and-Life-Lab/Atmospheric_Gas_Correction_Library.git
cd Atmospheric_Gas_Correction_Library
git lfs pull   # Certain large HDF files are used by the library and are too large to be stored in a normal GitHub repo, so Git LFS is needed to store them.
```

## 3. Create a fresh virtual environment

Create a fresh virtual environment, give it a name, and then activate it.

```bash
python3 -m venv {venv_name}
source {venv_name}/bin/activate
```

## 4. Install build dependencies manually

This step is required because the next step uses `--no-build-isolation`,
which skips pip's normal automatic build-dependency installation:

```bash
pip install scikit-build-core pybind11 setuptools_scm ninja cmake
```

## 5. Run pip install

The settings `editable.rebuild = true` and `build-dir = "build/{wheel_tag}"` are set in the pyproject.toml file, which sets up the library so that it will be rebuilt on import.
The built .so file will be placed in the directory build/{wheel_tag}, where {wheel_tag} is a directory whose name matches the wheel's naming convention (e.g. cp313-cp313-linux_x86_64).

```bash
pip install -e ".[test]" --no-build-isolation
```

## 6. Run pytest

The environment is now set up, and the following command can be run from the project root (Atmospheric_Gas_Correction_Library) to run the test suite:

```bash
pytest gas_corrections_lib/test
```

The rebuild-on-import settings that were configured as part of the previous steps will ensure that the library gets rebuilt right before the test suite runs.
This way, if you are making changes to the .cpp files and would like to quickly see the results without needing to manually rebuild, all you will need to do is run this pytest command, which will automatically rebuild and then run the tests.