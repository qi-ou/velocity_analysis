from merge_tif import *

def calc_N_from_EU(path):
    frame = os.path.basename(path)
    orientation = frame[3]
    track = frame[:4]

    n_file = os.path.join(path, "{}.geo.N.ml10.tif".format(frame))
    e_file = os.path.join(path, "{}.geo.E.ml10.tif".format(frame))
    u_file = os.path.join(path, "{}.geo.U.ml10.tif".format(frame))
    # n = OpenTif(n_file)
    e = OpenTif(e_file)
    u = OpenTif(u_file)

    if orientation == "A":
        inc_rad = np.arccos(u.data)
        head_rad = np.arcsin(n.data / np.sin(inc_rad))
    elif orientation == "D":
        inc_rad = np.arccos(u.data)
        head_rad = np.arcsin(- n.data / np.sin(inc_rad)) - np.pi
    else:
        raise ValueError("The 4th character of frameID is neither A nor D, please check your frame name.")

    azi_N = np.cos(head_rad)
    azi_E = np.sin(head_rad)
    azi_U = np.zeros(head_rad.shape)
    export_tif(N, e, "../NEU/{}_azi_N.tif".format(track))
    export_tif(azi_E, n, "../NEU/{}_azi_E.tif".format(track))
    export_tif(azi_U, n, "../NEU/{}_azi_U.tif".format(track))
    export_tif(n.data, n, "../NEU/{}_rng_N.tif".format(track))
    export_tif(e.data, n, "../NEU/{}_rng_E.tif".format(track))
    export_tif(u.data, n, "../NEU/{}_rng_U.tif".format(track))
