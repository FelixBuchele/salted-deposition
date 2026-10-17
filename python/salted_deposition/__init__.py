"""Optional subprocess interface to SALTED; numerical work stays in C++."""

from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
import os
import shutil
import subprocess
from tempfile import TemporaryDirectory

__all__ = ["Simulator", "SimulationResult"]


@dataclass(frozen=True)
class SimulationResult:
    """Completed output directory and existing exports; absent STL poses are None."""

    output_dir: Path
    vtk_path: Path | None
    stl_poses_path: Path | None


def _quoted(value: object) -> str:
    """Encode one C++ configuration value; reject embedded line breaks or NULs."""
    text = str(value)
    if any(c in text for c in "\r\n\0"):
        raise ValueError("Configuration overrides cannot contain line breaks or NULs")
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _required_file(path: Path) -> Path:
    """Return an existing output file or raise FileNotFoundError after execution."""
    if not path.is_file():
        raise FileNotFoundError(f"SALTED did not produce the required output: {path}")
    return path


class Simulator:
    """Invoke an already compiled SALTED executable; never build or install it."""

    def __init__(self, executable: str | Path = "deposit") -> None:
        """Find an existing executable path or PATH entry; raise FileNotFoundError.

        Relative executable paths are anchored to the caller's current directory
        at construction. Python and SALTED must run in the same environment.
        """
        name = os.fspath(executable)
        if isinstance(executable, Path) or os.path.dirname(name):
            path = Path(name).absolute()
            if not path.is_file():
                raise FileNotFoundError(f"SALTED executable not found: {path}")
            self._executable = str(path)
        else:
            found = shutil.which(name)
            if found is None:
                raise FileNotFoundError(f"SALTED executable not found on PATH: {name}")
            self._executable = str(Path(found).absolute())

    def _invoke(
        self, operation: str, config: Path | None, output: Path,
        overrides: dict[str, object],
    ) -> None:
        """Run once with a small temporary overlay, parsed exclusively by C++."""
        command = [self._executable, operation]
        if config is not None:
            command.append(str(config))
        context = TemporaryDirectory(prefix="salted-") if overrides else nullcontext()
        with context as directory:
            if directory is not None:
                overlay = Path(directory) / "overrides.conf"
                overlay.write_text(
                    "".join(f"{key} = {_quoted(value)}\n"
                            for key, value in overrides.items()),
                    encoding="utf-8",
                )
                if config is not None:
                    command.append("--overrides")
                command.append(str(overlay))
            command.extend(["--output", str(output)])
            subprocess.run(command, check=True, shell=False)

    def run(
        self, config: str | Path, *, output: str | Path,
        count: int | None = None, seed: int | None = None,
        shape: str | None = None, particle: str | Path | None = None,
        export_stl_poses: bool | None = None,
    ) -> SimulationResult:
        """Run a configuration with selected overrides and return verified exports.

        None preserves config values. Shape and particle are mutually exclusive;
        shapes are sphere/spiky/cube/stick. Configuration-relative paths stay in
        C++; method paths use the caller's current directory. VTK is always written.
        STL poses follow config unless explicitly enabled/disabled. C++ diagnostics
        stream to stdout/stderr; nonzero exits raise subprocess.CalledProcessError.
        Missing required exports raise FileNotFoundError; invalid API choices raise
        ValueError/TypeError. Existing outputs are refused by C++, not overwritten.
        """
        if shape is not None and particle is not None:
            raise ValueError("Specify shape or particle, not both")
        if shape is not None and shape not in ("sphere", "spiky", "cube", "stick"):
            raise ValueError("shape must be sphere, spiky, cube or stick")
        if export_stl_poses is not None and not isinstance(export_stl_poses, bool):
            raise TypeError("export_stl_poses must be bool or None")
        overrides: dict[str, object] = {}
        if count is not None:
            overrides["count"] = count
        if seed is not None:
            overrides["seed"] = seed
        if shape is not None:
            overrides["particle"] = shape
        if particle is not None:
            overrides["particle"] = Path(particle).absolute()
        if export_stl_poses is not None:
            overrides["stl_poses"] = str(export_stl_poses).lower()
        destination = Path(output).absolute()
        self._invoke("run", Path(config).absolute(), destination, overrides)
        vtk = _required_file(destination / "deposition.vtk")
        poses = destination / "stl_poses.vtk"
        if export_stl_poses is True:
            _required_file(poses)
        return SimulationResult(destination, vtk, poses if poses.is_file() else None)

    def generate_particle(
        self, mesh: str | Path, *, output: str | Path,
        max_spheres: int | None = None, units: str | None = None,
        config: str | Path | None = None,
    ) -> Path:
        """Generate through GEMSS and return the verified reusable particle path.

        Supply explicit units (mm/cm/m/legacy), or a generation config declaring
        them; advanced settings stay in that C++ config. None preserves its values
        or C++ defaults. All method paths use the caller's current directory.
        C++ preserves target-derived physical properties. Disabled GEMSS, invalid
        inputs and existing outputs raise subprocess.CalledProcessError, with live
        C++ diagnostics. Missing output raises FileNotFoundError; omitted units
        without config raises ValueError. No mesh data is read by Python.
        """
        if config is None and units is None:
            raise ValueError("Declare units or supply a generation config with units")
        overrides = {"stl": Path(mesh).absolute()}
        if max_spheres is not None:
            overrides["maximum_spheres"] = max_spheres
        if units is not None:
            overrides["units"] = units
        destination = Path(output).absolute()
        self._invoke(
            "generate", Path(config).absolute() if config is not None else None,
            destination, overrides,
        )
        return _required_file(destination)
