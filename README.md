# KVLCC2 Flow Distortion Maps

Time-averaged CFD velocity maps of the distorted flow around a full-scale KVLCC2 tanker at different drift angles in shallow water, and a Python interpolation module for applying these maps in ship maneuvering simulations.

## Reference

If you use this dataset or code, please cite:

> Schiaveto Neto, L.A., Rosman, P.C.C., Tannuri, E.A. (2026). "CFD investigation of KVLCC2 flow distortion with IDDES at different drift angles and its impact on tugboat current loads in ship maneuvering simulations." 
Associated manuscript under review (2026).



## Description

Ship maneuvering simulators commonly represent currents as undisturbed fields, neglecting the flow distortion caused by nearby vessels. This dataset provides time-averaged velocity maps around a KVLCC2 tanker hull that can be used to correct the ambient current experienced by tugboats operating near a large vessel.

The CFD simulations were performed using the IDDES turbulence model (Star-CCM+) for a full-scale KVLCC2 (Lpp = 320 m, B = 58 m, T = 20.8 m) in shallow water (h/T ≈ 1.2), with an incoming current of 1 m/s. The maps are time-averaged over 750 s of simulation (from t = 250 s to t = 1000 s).

## Coordinate System

The maps are defined in a **ship-fixed reference frame**:

```
                    Bow (+x)
                      ↑
                      |
         Port (+y) ←--●--→ Starboard (-y)
                      |
                    Stern (-x)
```

- **Origin**: ship center (midship, centerline)
- **x-axis**: positive toward the bow
- **y-axis**: positive toward port
- **Drift angle (β)**: angle between the ship heading and the current direction. β = 0° for current flowing from bow to stern; β < 0° for current from the starboard side.

## Drift Angles

Direct CFD simulations were performed for drift angles: 0°, −7.5°, −15°, −22.5°, −30°, −45°, −60°, −75°, −90°.

Additional angles were obtained by mirroring:
- **0° to +90°**: mirrored about the vessel x-axis (symmetric).
- **|β| > 90°**: mirrored about the x-axis for angles beyond 90° (less accurate due to bow-stern asymmetry; see manuscript for discussion).

## File Structure

```
├── README.md
├── LICENSE
├── kvlcc2_flow_distortion.py    # Python interpolation module
├── requirements.txt
└── maps/                        # CFD velocity map files
    ├── _SURFACE_MAP_SMS_X_-090.0.txt
    ├── _SURFACE_MAP_SMS_Y_-090.0.txt
    ├── _SURFACE_MAP_SMS_U_-090.0.txt
    ├── _SURFACE_MAP_SMS_V_-090.0.txt
    ├── _SURFACE_MAP_SMS_X_-075.0.txt
    ├── ...
    └── (4 files per drift angle: X, Y, U, V)
```

### Map file format

Each drift angle has four tab-separated text files (no header):
- `_SURFACE_MAP_SMS_X_{angle}.txt` — x-coordinates grid (m, KVLCC2 scale)
- `_SURFACE_MAP_SMS_Y_{angle}.txt` — y-coordinates grid (m, KVLCC2 scale)
- `_SURFACE_MAP_SMS_U_{angle}.txt` — time-averaged x-velocity component (m/s, normalized by 1 m/s reference current)
- `_SURFACE_MAP_SMS_V_{angle}.txt` — time-averaged y-velocity component (m/s, normalized by 1 m/s reference current)

The grid is uniform with 5 m spacing, covering a domain of 320 m × 320 m (one ship length) centered at the ship.

## KVLCC2 Dimensions

| Parameter | Value |
|-----------|-------|
| Lpp       | 320 m |
| B (beam)  | 58 m  |
| T (draft) | 20.8 m |
| Cb        | 0.81  |
| h/T       | ≈ 1.2 |

## Installation

```bash
pip install numpy scipy pandas
```

## Quick Start

```python
from kvlcc2_flow_distortion import build_interpolators, interp_velocity

# Load maps and build interpolators (KVLCC2 dimensions)
(u_interp, v_interp), (angles, y_axis, x_axis) = build_interpolators(
    maps_dir="./maps"
)

# Query velocity at midship, 0.5B to port, drift angle = -45°
u, v = interp_velocity(u_interp, v_interp,
                       drift_angle_deg=-45.0,
                       x=0.0,          # midship
                       y=29.0)         # 0.5 * B = 29 m (port side)

print(f"u = {u:.3f} m/s, v = {v:.3f} m/s")
```

### Scaling to a different vessel

To apply the maps to a vessel with different dimensions, pass the target vessel's length and beam when building the interpolators. The coordinates are scaled linearly.

```python
# Bulk carrier: Lpp = 350 m, B = 65 m
(u_interp, v_interp), _ = build_interpolators(
    maps_dir="./maps", lpp=350.0, beam=65.0
)
```

### Full transformation (global frame)

To compute the distorted current at a tug position during a maneuvering simulation, use the `compute_distorted_current` function, which implements Equations 10–13 of the manuscript:

```python
from kvlcc2_flow_distortion import (
    build_interpolators, global_to_ship_frame, compute_distorted_current
)

(u_interp, v_interp), _ = build_interpolators(maps_dir="./maps")

# Ship state (from simulator)
ship_heading_rad = 1.5    # radians
ship_ul = 1.0             # longitudinal speed, m/s
ship_ut = 0.2             # transversal speed, m/s
ambient_u = 0.5           # ambient current, global x, m/s
ambient_v = 0.3           # ambient current, global y, m/s

# Tug position (global frame)
tug_x_global = 1050.0
tug_y_global = 520.0
ship_x = 1000.0
ship_y = 500.0

# Step 1: Transform tug position to ship-fixed frame
tug_x_local, tug_y_local = global_to_ship_frame(
    tug_x_global, tug_y_global, ship_x, ship_y, ship_heading_rad
)

# Step 2: Compute distorted current at tug (global frame, m/s)
curr_u, curr_v = compute_distorted_current(
    u_interp, v_interp,
    ship_heading_rad=ship_heading_rad,
    ship_ul=ship_ul, ship_ut=ship_ut,
    ambient_current_u=ambient_u, ambient_current_v=ambient_v,
    tug_x_local=tug_x_local, tug_y_local=tug_y_local
)

print(f"Distorted current at tug: u = {curr_u:.3f} m/s, v = {curr_v:.3f} m/s")
```

## Limitations

- The CFD simulations were performed for a single depth condition (h/T ≈ 1.2).
- Propeller wake and free-surface deformation are not included.
- Velocity maps for |β| > 90° are obtained by mirroring and are less accurate due to bow-stern asymmetry.
- The maps represent time-averaged quantities; instantaneous turbulent fluctuations are not captured.
- Scaling to different vessel dimensions is approximate and most appropriate for full-form hulls with similar block coefficients.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for the code. The CFD data files are released under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).
