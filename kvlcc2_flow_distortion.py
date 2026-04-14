"""
kvlcc2_flow_distortion.py

Interpolation of time-averaged CFD velocity maps around a KVLCC2 tanker
at different drift angles, for ship maneuvering simulator applications.

Reference:
    Schiaveto Neto, L.A., Rosman, P.C.C., Tannuri, E.A. (2026).
    Associated manuscript under review (2026).

Ship-fixed coordinate system:
    Origin at ship center. x positive toward bow, y positive toward port.

All headings use the MATHEMATICAL convention:
    0 deg = positive x-axis, counterclockwise positive.
    If your simulator uses nautical convention (0=North, CW),
    convert first: heading_math = 90 - heading_nautical (in degrees).
"""

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
from typing import Tuple, Optional

__version__ = "1.0.0"

# KVLCC2 reference dimensions (full scale, meters)
KVLCC2_LPP = 320.0
KVLCC2_B = 58.0

# Available drift angles in the database (degrees).
# 0 to -90: direct CFD. Others: mirrored (see manuscript Sec. 2.1.4).
DRIFT_ANGLES = np.array([
    -180.0, -172.5, -165.0, -157.5, -150.0, -135.0, -120.0, -105.0, -90.0,
     -75.0,  -60.0,  -45.0,  -30.0,  -22.5,  -15.0,   -7.5,    0.0,
       7.5,   15.0,   22.5,   30.0,   45.0,   60.0,   75.0,   90.0,
     105.0,  120.0,  135.0,  150.0,  165.0,  180.0
], dtype=float)


def _format_angle(angle):
    s = "%05.1f" % abs(angle)
    return ("-" if angle < 0 else "") + s


def load_map(maps_dir, angle, lpp=KVLCC2_LPP, beam=KVLCC2_B):
    """Load a single CFD velocity map for a given drift angle.

    Coordinates are scaled from KVLCC2 to target vessel using lpp and beam.

    Returns: (x_axis, y_axis, U, V)
    """
    maps_dir = Path(maps_dir)
    a = _format_angle(angle)
    kw = dict(sep="\t", engine="python", header=None)

    X = pd.read_csv(maps_dir / f"_SURFACE_MAP_SMS_X_{a}.txt", **kw).to_numpy().T
    Y = pd.read_csv(maps_dir / f"_SURFACE_MAP_SMS_Y_{a}.txt", **kw).to_numpy().T
    U = pd.read_csv(maps_dir / f"_SURFACE_MAP_SMS_U_{a}.txt", **kw).to_numpy().T
    V = pd.read_csv(maps_dir / f"_SURFACE_MAP_SMS_V_{a}.txt", **kw).to_numpy().T

    x_axis = X[0, :] * lpp / KVLCC2_LPP
    y_axis = Y[:, 0] * beam / KVLCC2_B
    return x_axis, y_axis, U, V


def build_interpolators(maps_dir, lpp=KVLCC2_LPP, beam=KVLCC2_B,
                        drift_angles=None):
    """Build 3D interpolators for u and v velocity components.

    Returns:
        (u_interp, v_interp) - query with interp((angle, y, x))
        (angles, y_axis, x_axis) - grid axes
    """
    if drift_angles is None:
        drift_angles = DRIFT_ANGLES

    x_axis, y_axis = None, None
    U_list, V_list = [], []

    for angle in drift_angles:
        x_a, y_a, U, V = load_map(maps_dir, angle, lpp, beam)
        if x_axis is None:
            x_axis, y_axis = x_a, y_a
        U_list.append(U)
        V_list.append(V)

    U_all = np.stack(U_list, axis=0)
    V_all = np.stack(V_list, axis=0)

    u_interp = RegularGridInterpolator(
        (drift_angles, y_axis, x_axis), U_all,
        bounds_error=False, fill_value=np.nan, method="linear")
    v_interp = RegularGridInterpolator(
        (drift_angles, y_axis, x_axis), V_all,
        bounds_error=False, fill_value=np.nan, method="linear")

    return (u_interp, v_interp), (drift_angles, y_axis, x_axis)


def interp_velocity(u_interp, v_interp, drift_angle_deg, x, y):
    """Interpolate CFD velocity at a point in the ship-fixed frame.

    Returns (u, v) normalized by 1 m/s reference current. NaN if outside domain.
    """
    pt = np.array([[float(drift_angle_deg), float(y), float(x)]])
    return u_interp(pt).item(), v_interp(pt).item()


def global_to_ship_frame(x_global, y_global, x_ship, y_ship, heading_rad):
    """Transform global coordinates to ship-fixed frame (Eq. 12).

    WARNING: heading_rad must be in mathematical convention
    (0 = positive x-axis, CCW positive). To convert from nautical:
        heading_math_rad = np.deg2rad(90 - heading_nautical_deg)

    Returns (x_local, y_local).
    """
    c, s = np.cos(heading_rad), np.sin(heading_rad)
    dx = np.asarray(x_global) - x_ship
    dy = np.asarray(y_global) - y_ship
    return c * dx + s * dy, -s * dx + c * dy


def compute_distorted_current(u_interp, v_interp,
                              ship_heading_rad, ship_ul, ship_ut,
                              ambient_current_u, ambient_current_v,
                              tug_x_local, tug_y_local):
    """Compute distorted current at a tug position (Eqs. 10-13).

    WARNING: ship_heading_rad must be in mathematical convention.
    To convert from nautical:
        heading_math_rad = np.deg2rad(90 - heading_nautical_deg)

    Args:
        u_interp, v_interp: from build_interpolators()
        ship_heading_rad: ship heading, math convention (rad)
        ship_ul, ship_ut: ship velocity in ship-fixed frame (m/s)
        ambient_current_u, ambient_current_v: undisturbed current, global (m/s)
        tug_x_local, tug_y_local: tug position in ship-fixed frame (m)

    Returns:
        (curr_u, curr_v) in global frame (m/s). NaN if tug outside map domain.
    """
    # Eq. 10: ship velocity relative to water
    c_h = np.cos(-ship_heading_rad)
    s_h = np.sin(-ship_heading_rad)
    ul_w = ship_ul - (c_h * ambient_current_u - s_h * ambient_current_v)
    ut_w = ship_ut - (s_h * ambient_current_u + c_h * ambient_current_v)

    beta_deg = np.degrees(np.arctan2(ut_w, ul_w))
    speed_w = np.sqrt(ul_w**2 + ut_w**2)

    # Interpolate CFD map
    u_cfd, v_cfd = interp_velocity(u_interp, v_interp,
                                   beta_deg, tug_x_local, tug_y_local)
    if np.isnan(u_cfd) or np.isnan(v_cfd):
        return np.nan, np.nan

    # Eq. 13: scale, add ship velocity, rotate to global
    c_h2 = np.cos(ship_heading_rad)
    s_h2 = np.sin(ship_heading_rad)
    u_s = speed_w * u_cfd + ship_ul
    v_s = speed_w * v_cfd + ship_ut
    return c_h2 * u_s - s_h2 * v_s, s_h2 * u_s + c_h2 * v_s


# ===========================================================================
#  EXAMPLE
# ===========================================================================

if __name__ == "__main__":
    import sys

    maps_dir = "./maps"
    if not Path(maps_dir).is_dir():
        print(f"Error: '{maps_dir}' not found. Set maps_dir to your maps path.")
        sys.exit(1)

    print("Loading CFD maps...")
    (u_interp, v_interp), (angles, y_axis, x_axis) = build_interpolators(
        maps_dir=maps_dir)
    print(f"  {len(angles)} drift angles loaded")
    print(f"  Domain: x=[{x_axis.min():.0f}, {x_axis.max():.0f}] m, "
          f"y=[{y_axis.min():.0f}, {y_axis.max():.0f}] m")

    # Query: midship, 0.5B to port, beta = -45 deg
    u, v = interp_velocity(u_interp, v_interp, -45.0, 0.0, 0.5 * KVLCC2_B)
    print(f"\nAt beta=-45 deg, midship, 0.5B port:")
    print(f"  u={u:.4f}, v={v:.4f}, |V|={np.sqrt(u**2+v**2):.4f} (normalized)")
