# SolidificationAlCu benchmark

Diffusion-controlled directional solidification of an Al-Cu alloy (SolidificationAlCu).

A single dendrite growth is simulated under directional solidification. Material
parameters for an Al - 4.0 at.% Cu alloy:

| Parameter                                   | Value                         |
|---------------------------------------------|-------------------------------|
| Melting point of pure Al                    | T_melt  = 933.6 K             |
| Entropy change during solidification        | ΔS = 10⁶ J/(m³·K)             |
| Diffusion coefficient of Cu in liquid Al    | D_liquid = 3×10⁻⁹ m²/s        |
| Partition coefficient                       | k = C_solid / C_liquid = 7.14 |
| Liquidus slope                              | 2.6 K/at%                     |
| Temperature gradient in Z direction         | Gz = 1000 K/m                 |
| Cooling rate                                | dT/dt = -0.4 K/s              |
| Pulling velocity                            | V = 0.0004 m/s                |
| Grid size                                   | 71×71×301 points              |
| Grid spacing                                | dx = 1×10⁻⁶ m                 |
| Time step                                   | dt = 5×10⁻⁵ s                 |

The simulation box starts to move with the growing dendrite once the tip reaches
Nz/3. After the initial transient, a steady dendrite shape grows at constant
velocity equal to the pulling velocity.

## Usage
```bash
./SolidificationAlCu 
```
The `ProjectInput.opi` file contains the necessary input for the simulation. 
Results will be placed in the `VTK` directory and can be visualised 
using ParaView (www.paraview.org).

## Postprocessing

The video post processing can be started with:
```bash
./post.sh
```
### Workflow:
1. **Rendering VTK files** to images (`post-render.py`)  
2. **Cropping images** for clean presentation (`post-crop.sh`)  
3. **Compiling images into a video** (`post-video.sh`)

### Output directories
* `Images` for png images
* `Videos` for webM video

### Requirements
* ParaView 5.11+ (for post-render.py)
* Bash and ImageMagick (for post-crop.sh)
* ffmpeg (for post-video.sh)
