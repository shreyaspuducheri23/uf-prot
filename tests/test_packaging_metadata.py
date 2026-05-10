"""Sanity checks for packaging metadata."""
import tomllib
from pathlib import Path


def test_project_dependencies_are_scoped_under_project():
    data = tomllib.loads(Path("pyproject.toml").read_text())
    assert "dependencies" in data["project"]
    deps = data["project"]["dependencies"]
    assert isinstance(deps, list)
    assert any(dep.startswith("pandas") for dep in deps)


def test_wheel_targets_real_package():
    data = tomllib.loads(Path("pyproject.toml").read_text())
    packages = data["tool"]["hatch"]["build"]["targets"]["wheel"]["packages"]
    assert "scripts" in packages
    assert "code" not in packages
