"""
Plot f(x, y) as a single blue color whose transparency (alpha)
encodes the value of f at each point.

The key idea:
  - Every pixel gets the SAME RGB color (blue).
  - f(x, y) is normalized to [0, 1] and used as the alpha channel.
  - This is done by building an RGBA image manually and showing it
    with imshow (imshow understands a 4-channel array as RGBA).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def f(x, y):
    """Define your function here. Must return values in a bounded range."""
    return np.exp(-(x**2 + y**2))          # example: a 2D Gaussian bump


def draw_sphere(ax, center, radius, color=(0.5, 0.5, 0.5),
                 light_dir=(-1, 1, 1.5), resolution=100):
    """
    Draw a small plain shaded sphere centered at `center` with the given
    `radius`, using simple Lambertian shading so it reads as 3D.
    `color` is the base RGB of the sphere.
    """
    cx, cy = center
    lin = np.linspace(-1, 1, resolution)
    U, V = np.meshgrid(lin, lin)
    R2 = U**2 + V**2

    with np.errstate(invalid="ignore"):
        W = np.sqrt(np.clip(1 - R2, 0, 1))

    inside = R2 <= 1.0

    light = np.array(light_dir, dtype=float)
    light /= np.linalg.norm(light)
    shade = U * light[0] + V * light[1] + W * light[2]
    shade = np.clip(shade, 0, 1)
    shade = 0.25 + 0.75 * shade  # keep some ambient light, avoid pure black

    rgba = np.zeros((resolution, resolution, 4))
    rgba[..., 0] = color[0] * shade
    rgba[..., 1] = color[1] * shade
    rgba[..., 2] = color[2] * shade
    rgba[..., 3] = np.where(inside, 1.0, 0.0)

    extent = (cx - radius, cx + radius, cy - radius, cy + radius)
    im = ax.imshow(rgba, extent=extent, origin="lower", zorder=5)
    clip_circle = Circle((cx, cy), radius, transform=ax.transData)
    im.set_clip_path(clip_circle)


def draw_arrow(ax, center, radius, angle=None, color="black", rng=None, zorder=6):
    """
    Draw an arrow starting at `center`, pointing in a random direction
    (or a given `angle` in radians), with:
      - length = 3 * radius
      - width  = diameter / 3 = (2 * radius) / 3
    """
    if angle is None:
        rng = rng or np.random.default_rng()
        angle = rng.uniform(0, 2 * np.pi)

    cx, cy = center
    length = 3 * radius
    width = (2 * radius) / 3

    dx = length * np.cos(angle)
    dy = length * np.sin(angle)

    ax.arrow(
        cx, cy, dx, dy,
        width=width,
        head_width=width * 2.5,
        head_length=width * 2.5,
        length_includes_head=True,
        color=color,
        zorder=zorder,
    )


def generate_random_spheres(
    n,
    center=(0.0, 0.0),
    std=1.0,
    sphere_radius=0.2,
    color=(0.5, 0.5, 0.5),
    seed=None,
):
    """
    Randomly generate `n` sphere positions drawn from an isotropic 2D
    Gaussian distribution.

    center : (x, y) mean of the distribution
    std    : standard deviation (same in x and y) — controls how spread
             out the spheres are around `center`
    sphere_radius : radius to assign to every sphere (float), or a
             (low, high) tuple to draw each sphere's radius uniformly
             from that range
    color  : RGB to assign to every sphere, or a list of RGB tuples to
             pick from at random for each sphere
    seed   : optional int for reproducible randomness

    Returns a list of dicts compatible with `plot_alpha_map(..., spheres=...)`.
    """
    rng = np.random.default_rng(seed)

    positions = rng.normal(loc=center, scale=std, size=(n, 2))
    angles = rng.uniform(0, 2 * np.pi, size=n)

    if isinstance(sphere_radius, (tuple, list)):
        radii = rng.uniform(sphere_radius[0], sphere_radius[1], size=n)
    else:
        radii = np.full(n, sphere_radius)

    if isinstance(color, (list, tuple)) and len(color) > 0 and isinstance(color[0], (list, tuple)):
        colors = [color[rng.integers(0, len(color))] for _ in range(n)]
    else:
        colors = [color] * n

    spheres = [
        {
            "center": tuple(positions[i]),
            "radius": float(radii[i]),
            "color": colors[i],
            "angle": float(angles[i]),
        }
        for i in range(n)
    ]
    return spheres


def plot_alpha_map(
    f,
    x_range=(-3, 3),
    y_range=(-3, 3),
    resolution=500,
    blue=(0.0, 0.0, 1.0),   # RGB of the color used everywhere (pure blue)
    background="white",
    normalize=True,
    spheres=None,           # list of dicts: {"center": (x, y), "radius": r, "color": (r,g,b)}
):
    # Build the grid
    x = np.linspace(x_range[0], x_range[1], resolution)
    y = np.linspace(y_range[0], y_range[1], resolution)
    X, Y = np.meshgrid(x, y)

    # Evaluate f on the grid
    Z = f(X, Y)

    # Normalize f to [0, 1] to use as alpha
    if normalize:
        z_min, z_max = np.nanmin(Z), np.nanmax(Z)
        if z_max > z_min:
            alpha = (Z - z_min) / (z_max - z_min)
        else:
            alpha = np.zeros_like(Z)
    else:
        # Assume Z is already in [0, 1]; just clip for safety
        alpha = np.clip(Z, 0, 1)

    # Build an RGBA image: same blue everywhere, alpha varies per pixel
    rgba = np.zeros((resolution, resolution, 4))
    rgba[..., 0] = blue[0]
    rgba[..., 1] = blue[1]
    rgba[..., 2] = blue[2]
    rgba[..., 3] = alpha

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_facecolor(background)
    ax.imshow(
        rgba,
        extent=(x_range[0], x_range[1], y_range[0], y_range[1]),
        origin="lower",
        aspect="equal",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("f(x, y) shown as alpha of a fixed blue")

    # Draw the requested spheres (and their arrows) on top of the map
    if spheres:
        for s in spheres:
            draw_sphere(
                ax,
                center=s["center"],
                radius=s.get("radius", 0.2),
                color=s.get("color", (0.5, 0.5, 0.5)),
            )
            draw_arrow(
                ax,
                center=s["center"],
                radius=s.get("radius", 0.2),
                angle=s.get("angle", None),
            )

    plt.tight_layout()
    plt.savefig("alpha_plot.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    x_range = (-3, 3)
    y_range = (-3, 3)
    plot_side = x_range[1] - x_range[0]
    r = plot_side / 50  # sphere radius ~ 1/50 of the plot side

    spheres = generate_random_spheres(
        n=15,
        center=(0.0, 0.0),
        std=1.0,
        sphere_radius=r,
        color=[(0.8, 0.1, 0.1), (0.1, 0.6, 0.1), (0.9, 0.9, 0.1)],
        seed=42,
    )
    plot_alpha_map(f, x_range=x_range, y_range=y_range, spheres=spheres)
