def mroll(start, stop, current_view, axis=(0,1,0), degrees=360.0):
    # snapshot current view
    v0 = current_view
    frames = stop - start + 1

    for i in range(frames):
        frac = float(i) / frames
        angle = math.radians(frac * degrees)

        # rotate axis-angle -> rotation matrix
        ax, ay, az = axis
        c = math.cos(angle)
        s = math.sin(angle)
        C = 1 - c
        R = [
            ax*ax*C + c,     ax*ay*C - az*s, ax*az*C + ay*s,
            ay*ax*C + az*s,  ay*ay*C + c,    ay*az*C - ax*s,
            az*ax*C - ay*s,  az*ay*C + ax*s, az*az*C + c,
        ]

        # apply rotation to the camera orientation part of the view
        # view[0:9] is the 3×3 rotation matrix, rest unchanged
        new_view = list(v0)
        #old = v0[0:9]
        old=[0,0,0,0,0,0,0,0,0]
        new = [
            R[0]*old[0] + R[1]*old[3] + R[2]*old[6],
            R[0]*old[1] + R[1]*old[4] + R[2]*old[7],
            R[0]*old[2] + R[1]*old[5] + R[2]*old[8],

            R[3]*old[0] + R[4]*old[3] + R[5]*old[6],
            R[3]*old[1] + R[4]*old[4] + R[5]*old[7],
            R[3]*old[2] + R[4]*old[5] + R[5]*old[8],

            R[6]*old[0] + R[7]*old[3] + R[8]*old[6],
            R[6]*old[1] + R[7]*old[4] + R[8]*old[7],
            R[6]*old[2] + R[7]*old[5] + R[8]*old[8],
        ]
        new_view[0:9] = new
        return new_view
        #cmd.mview("store", start + i, view=new_view)