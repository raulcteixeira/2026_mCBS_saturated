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
from matplotlib.patches import FancyArrow, Rectangle
from matplotlib.transforms import Affine2D

def beam(w, theta, lmbda, x, y):
    """
    Compute the intensity of a Gaussian beam with waist w and tilt angle theta
    at the points (x, y).

    Parameters:
    - w: waist of the Gaussian beam
    - theta: tilt angle in radians
    - x, y: coordinates where the intensity is computed

    Returns:
    - Intensity of the Gaussian beam at (x, y)
    """
    # Rotate coordinates by theta
    x_rot = x * np.cos(theta) + y * np.sin(theta)
    y_rot = -x * np.sin(theta) + y * np.cos(theta)

    # Compute the Gaussian intensity profile
    field = np.exp(- x_rot**2 / w**2) * np.exp(1j * 2*np.pi * y_rot/lmbda)
    
    return field

def superposition (x, y, w, theta, lmbda):
    # Example: a 2D Gaussian bump
    return np.abs(beam(w, theta, lmbda, x, y) + beam(w, - theta, lmbda, x, y))**2  # Squared to represent intensity

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

    positions_polar = rng.uniform(low = 0, high = 1, size=(n, 2))
    positions = positions_polar.copy()
    positions[:, 0] = np.sqrt(positions_polar[:, 0]) * std *np.cos(positions[:, 1]*2*np.pi) + center[0]
    positions[:, 1] = np.sqrt(positions_polar[:, 0]) * std *np.sin(positions[:, 1]*2*np.pi) + center[1]

    if isinstance(sphere_radius, (tuple, list)):
        radii = rng.uniform(sphere_radius[0], sphere_radius[1], size=n)
    else:
        radii = np.full(n, sphere_radius)

    if isinstance(color, (list, tuple)) and len(color) > 0 and isinstance(color[0], (list, tuple)):
        colors = [color[rng.integers(0, len(color))] for _ in range(n)]
    else:
        colors = [color] * n

    spheres = [
        {"center": tuple(positions[i]), "radius": float(radii[i]), "color": colors[i]}
        for i in range(n)
    ]
    return spheres

class Arrow2D:
    def __init__(self, position, direction, length, width):
        """
        position: (x, y) tuple - tail of the arrow
        direction: (dx, dy) tuple - direction vector (will be normalized)
        length: float - total length of the arrow
        width: float - shaft width (head scales proportionally)
        """
        self.position = np.array(position, dtype=float)
        direction = np.array(direction, dtype=float)
        norm = np.linalg.norm(direction)
        if norm == 0:
            raise ValueError("direction vector cannot be zero")
        self.direction = direction / norm
        self.length = length
        self.width = width

    def draw(self, ax=None, color="black", **kwargs):
        if ax is None:
            ax = plt.gca()

        dx, dy = self.direction * self.length
        x, y = self.position

        head_length = np.min([self.length * 0.2, 0.3])  # Limit head length to avoid overly large heads

        arrow = FancyArrow(
            x, y, dx, dy,
            width=self.width,
            head_width=self.width * 3,
            head_length=head_length,
            length_includes_head=True,
            color=color,
            **kwargs
        )
        ax.add_patch(arrow)
        return arrow

class Rectangle2D:
    def __init__(self, position, width, height, angle=0.0, anchor="center"):
        """
        position: (x, y) tuple
        width: float - rectangle width
        height: float - rectangle height
        angle: float - rotation in degrees (counterclockwise)
        anchor: "center" or "corner" - whether `position` refers to the
                rectangle's center or its bottom-left corner (before rotation)
        """
        self.position = np.array(position, dtype=float)
        self.width = width
        self.height = height
        self.angle = angle
        self.anchor = anchor

    def draw(self, ax=None, color="blue", fill=False, **kwargs):
        if ax is None:
            ax = plt.gca()

        x, y = self.position

        if self.anchor == "center":
            corner = (x - self.width / 2, y - self.height / 2)
        elif self.anchor == "corner":
            corner = (x, y)
        else:
            raise ValueError("anchor must be 'center' or 'corner'")

        rect = Rectangle(
            corner, self.width, self.height,
            edgecolor='black', facecolor = color, fill=fill, **kwargs
        )

        # Rotate about the given position point
        t = Affine2D().rotate_deg_around(x, y, self.angle) + ax.transData
        rect.set_transform(t)

        ax.add_patch(rect)
        return rect

def plot_superposition(
    w,
    theta,
    lmbda,
    x_range=(-10, 0),
    y_range=(-3, 3),
    resolution=500,
    blue=(0.0, 0.0, 1.0),   # RGB of the color used everywhere (pure blue)
    background="white",
    normalize=True,
):
    # Build the grid
    x = np.linspace(x_range[0], x_range[1], resolution)
    y = np.linspace(y_range[0], y_range[1], resolution)
    X, Y = np.meshgrid(x, y)

    # Evaluate f on the grid
    Z = superposition(X, Y, w, theta, lmbda)

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

    atom_cloud_center = (-w/3*np.tan(theta), 0.0)
    atom_cloud_radius = w*np.tan(theta)/10

    spheres = generate_random_spheres(
        n=60,
        center=atom_cloud_center,
        std=atom_cloud_radius,
        sphere_radius=0.06,
        color=(0.0, 0.5, 0.0)
    )
    for sphere in spheres:
        draw_sphere(ax, sphere["center"], sphere["radius"], color=sphere["color"])

    arrow_pos = 8
    arrow_length = 1.0
    a1 = Arrow2D(position=(-arrow_pos*np.sin(theta), -arrow_pos*np.cos(theta)), direction=(np.sin(theta), np.cos(theta)), length=arrow_length, width=0.1)
    a1.draw(ax, color="blue")

    arrow_pos2 = arrow_pos - arrow_length
    a2 = Arrow2D(position=(-arrow_pos2*np.sin(theta), arrow_pos2*np.cos(theta)), direction=(-np.sin(theta), np.cos(theta)), length=arrow_length, width=0.1)
    a2.draw(ax, color="blue")

    # N_arrow_scat = 10
    # arrow_pos3 = arrow_pos*0.5
    # for jj in range(N_arrow_scat):
    #     theta_em = 0.85*np.pi*(1 - 0.7*jj/(N_arrow_scat - 1))*(N_arrow_scat - 1)/(N_arrow_scat)
    #     a3 = Arrow2D(position=(-arrow_pos3*np.sin(theta_em), -arrow_pos3*np.cos(theta_em)), direction=(-np.sin(theta_em), -np.cos(theta_em)), length=arrow_length*0.7, width=0.05)
    #     a3.draw(ax, color="green")

    theta_em = np.pi/2-2*(np.pi/2-theta)
    a3 = Arrow2D(position=(atom_cloud_center[0], atom_cloud_center[1]), direction=(-np.sin(theta_em), np.cos(theta_em)), length=4*arrow_length, width=0.05)
    a3.draw(ax, color="green")

    a4 = Arrow2D(position=(-atom_cloud_center[0], atom_cloud_center[1]), direction=(-np.sin(theta_em), np.cos(theta_em)), length=7.8*arrow_length, width=0.05)
    a4.draw(ax, color="green")

    a5 = Arrow2D(position=(atom_cloud_center[0], atom_cloud_center[1]), direction=(np.sin(theta_em), np.cos(theta_em)), length=7.8*arrow_length, width=0.05)
    a5.draw(ax, color="green")

    mirror_thickness = 0.2
    mirror_height = y_range[1] - y_range[0]
    mirror_center_y = y_range[0] + mirror_height / 2
    mirror = Rectangle2D(position=(mirror_thickness/2, mirror_center_y), width=mirror_thickness, height=mirror_height)
    mirror.draw(ax, color=(0.7,0.7,1.0), fill=True)

    #ax.set_xlabel("x")
    #ax.set_ylabel("y")
    #ax.set_title("f(x, y) shown as alpha of a fixed blue")
    ax.set_axis_off()
    ax.set_xlim((x_range[0], x_range[1] + mirror_thickness + 0.01))
    ax.set_ylim((y_range[0] - 0.01, y_range[1] + 0.01))
    plt.tight_layout()
    plt.savefig("alpha_plot.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    plot_superposition(w=1.0, theta=0.9*np.pi/2, lmbda = 0.4)
    