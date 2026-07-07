# Moon-Landing-Control

**Main_feedback** reads the control (state-feedback) map from a `.mat` file and simulates a controlled lunar descent trajectory. The feedback map is defined in terms of the dimensionless variables \(v_{xr}\) and \(v_{yr}\), where the horizontal and vertical velocities are normalized by \(\sqrt{A_0 r}\), with \(A_0\) the reference thrust acceleration and \(r\) the current altitude.

The script first loads the precomputed feedback table (`Landing_K_2D_C_MEM6g.mat`), containing the state grid \((v_{xr}, v_{yr})\), the corresponding control law (thrust angle \(\alpha\)), and the estimated remaining time-to-go \(V\). 

As output, the code provides the time histories of altitude, horizontal and vertical velocities, normalized velocities, thrust angle, and estimated time-to-go, together with a comparison between the feedback-guided solution and the corresponding extremal trajectory.
