import copy
import math
from pathlib import Path


def normalize_r_checks(r_checks):
    """Validate, deduplicate, and sort one or more neighborhood radii."""
    if r_checks is None:
        raise ValueError("Detection method 'radius' requires --r_check to be set.")

    if isinstance(r_checks, (int, float)):
        radii = [float(r_checks)]
    else:
        radii = [float(radius) for radius in r_checks]

    if any(not math.isfinite(radius) or radius <= 0 for radius in radii):
        raise ValueError("All values supplied to --r_check must be finite positive numbers.")

    if len(set(radii)) != len(radii):
        raise ValueError("Duplicate values are not allowed in --r_check.")

    return sorted(radii)


def radius_token(radius):
    """Return a filesystem-safe, deterministic label for a radius."""
    value = repr(float(radius))
    if value.endswith('.0'):
        value = value[:-2]
    return value.replace('+', 'p').replace('-', 'm').replace('.', 'p')


def add_radius_suffix(file_path, radius):
    """Add a radius suffix immediately before a result file's extension."""
    path = Path(file_path)
    return str(path.with_name(f"{path.stem}_r_{radius_token(radius)}{path.suffix}"))


def build_radius_options(options):
    """Build isolated option objects for one or more complete SCRIN runs."""
    if options.detection_method != 'radius':
        if isinstance(options.r_check, (list, tuple)):
            if len(options.r_check) > 1:
                raise ValueError("Multiple --r_check values are only supported with detection method 'radius'.")
            if len(options.r_check) == 1:
                options.r_check = options.r_check[0]
        return [options]

    radii = normalize_r_checks(options.r_check)
    if len(radii) == 1:
        options.r_check = radii[0]
        return [options]

    run_options = []
    for radius in radii:
        radius_options = copy.copy(options)
        radius_options.r_check = radius
        radius_options.save_path = add_radius_suffix(options.save_path, radius)
        if options.intermediate_dir is not None:
            radius_options.intermediate_dir = str(
                Path(options.intermediate_dir) / f"r_{radius_token(radius)}"
            )
        run_options.append(radius_options)

    return run_options
