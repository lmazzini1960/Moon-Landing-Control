# Moon-Landing-Control

## Main_feedback

`Main_feedback` reads the control (state-feedback) map from a `.mat` file and simulates a controlled lunar descent trajectory.

The feedback map is expressed in terms of the dimensionless variables `vxr` and `vyr`, defined as:

- `vxr = vx / sqrt(A0 * r)`
- `vyr = vy / sqrt(A0 * r)`

where:

- `vx` is the horizontal velocity,
- `vy` is the vertical velocity,
- `A0` is the reference thrust acceleration,
- `r` is the current altitude.

The script loads the precomputed feedback table from `Landing_K_2D_C_MEM6g.mat`, which contains:

- the state grid (`vxr`, `vyr`),
- the feedback control law (`alfa0`, i.e. thrust angle),
- the estimated remaining time-to-go (`V`).

As output, the code provides the time histories of:

- altitude,
- horizontal velocity,
- vertical velocity,
- normalized velocities (`vxr`, `vyr`),
- thrust angle,
- estimated time-to-go.

The script also compares the feedback-guided trajectory with the corresponding extremal trajectory.
