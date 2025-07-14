import pytest
import sys
from pathlib import Path

# Add the bin directory to Python path
project_root = Path(__file__).parent.parent
bin_dir = project_root / "bin"
sys.path.insert(0, str(bin_dir))

@pytest.fixture(scope="session")
def project_root():
    """Return the project root directory"""
    return Path(__file__).parent.parent

@pytest.fixture(scope="session")
def bin_dir(project_root):
    """Return the bin directory"""
    return project_root / "bin"
