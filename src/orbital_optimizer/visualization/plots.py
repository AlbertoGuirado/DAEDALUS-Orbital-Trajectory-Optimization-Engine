"""Export scientific validation figures with Matplotlib's noninteractive backend."""
from pathlib import Path
import numpy as np
from ..simulation.propagator import PropagationResult
from ..orbital.mechanics import specific_energy, specific_angular_momentum, eccentricity_vector


def save_validation_plots(result: PropagationResult, reference: np.ndarray, output: Path) -> None:
    """Save SI trajectory/reference and invariant-error PNGs; no interactive UI.

    Input states are (N,6), position m and velocity m/s. Positions are displayed
    in km and time in minutes only at this presentation boundary. Invariants
    follow epsilon=v²/2-mu/r, h=r cross v, e=v cross h/mu-r/r.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8,7))
    ax = fig.add_subplot(111,projection="3d")
    xyz, analytic = result.states[:,:3]/1000, reference[:,:3]/1000
    ax.plot(*xyz.T, label=f"Numerical {result.method}", linewidth=2)
    ax.plot(*analytic.T, linestyle="--", label="Analytic circular solution", linewidth=1)
    ax.scatter(0,0,0,color="black",s=30,label="Central point mass")
    ax.scatter(*xyz[0],color="tab:green",label="Initial state")
    bound = np.max(np.abs(xyz))*1.1
    ax.set(xlim=(-bound,bound),ylim=(-bound,bound),zlim=(-bound,bound),
           xlabel="x (km)",ylabel="y (km)",zlabel="z (km)",title="Two-body orbital propagation")
    ax.set_box_aspect((1,1,1))
    ax.legend(loc="upper left",fontsize=8)
    fig.tight_layout()
    fig.savefig(output/"trajectory.png",dpi=160)
    plt.close(fig)

    energy = np.array([specific_energy(s,result.mu) for s in result.states])
    h = np.array([specific_angular_momentum(s) for s in result.states])
    e = np.array([eccentricity_vector(s,result.mu) for s in result.states])
    values = [np.abs(energy-energy[0])/abs(energy[0]),
              np.linalg.norm(h-h[0],axis=1)/np.linalg.norm(h[0]),
              np.linalg.norm(e-e[0],axis=1),
              np.linalg.norm(result.states[:,:3]-reference[:,:3],axis=1)]
    labels = ["Relative energy drift", "Relative angular momentum drift",
              "Absolute eccentricity-vector drift", "Position error (m)"]
    fig, axes = plt.subplots(2,2,figsize=(10,7))
    for ax, value, label in zip(axes.flat,values,labels):
        ax.plot(result.times/60,value)
        ax.set(xlabel="Elapsed time (min)",ylabel=label)
        ax.ticklabel_format(axis="y",style="sci",scilimits=(-3,3))
        ax.grid(alpha=.3)
    fig.suptitle("Physical consistency and analytic error (sampled)")
    fig.tight_layout()
    fig.savefig(output/"consistency.png",dpi=160)
    plt.close(fig)
