"""Run the compiled SALTED program from Python; all numerical work stays in C++."""

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
    """Paths from a completed run: its output directory, sphere VTK and optional STL poses.

    Paths are absolute. stl_poses_path is None when no pose file was produced.
    """

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
    """Run an existing SALTED executable without building or installing it.

    C++ handles configuration and scientific calculations. Failed commands raise
    subprocess.CalledProcessError, with diagnostics shown in the terminal.
    Existing outputs are refused by SALTED; missing output files raise
    FileNotFoundError.
    """

    def __init__(self, executable: str | Path = "deposit") -> None:
        """Find executable by path or on PATH; raise FileNotFoundError if absent.

        Relative executable paths use the current directory at construction.
        Python and SALTED must run in the same environment.
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
        """Run one command, writing only overrides to a temporary C++ configuration.

        C++ combines these values with config; Python does not parse that file.
        """
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
        """Run config in output and return the paths of the generated exports.

        count and seed override the file only when supplied. Choose either a
        built-in shape (sphere/spiky/cube/stick) or a particle file; None keeps
        the configuration value. Conflicting choices raise ValueError.

        VTK is always written. export_stl_poses follows config when None; True
        requires a pose file or raises FileNotFoundError. The flag must be bool
        or None. Method paths use the current directory; paths inside config
        remain relative to that file. Execution errors propagate as in Simulator.
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
        """Generate a reusable particle file from STL through GEMSS; return its path.

        Declare mesh units (mm/cm/m/legacy), or supply a generation config that
        declares them; omitting both raises ValueError. max_spheres limits the
        sphere count. Advanced settings stay in config; None keeps its values
        or C++ defaults. Method paths use the current directory.

        C++ calculates physical properties from the target mesh and handles all
        geometry processing. Disabled GEMSS, invalid inputs and file errors
        propagate as described in Simulator.
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
