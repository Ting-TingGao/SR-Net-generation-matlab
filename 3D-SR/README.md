# 3D Scale-Rich Network Generation (MATLAB)
This code generates a **3D scale-rich network** inside the unit cube $([0,1]^3$) by repeatedly slicing polyhedra with random planes. At step *t*, a void polyhedron is chosen with probability proportional to its volume and split by a **mid plane** plus two **offset planes** whose thickness follows
$\lambda_t = \lambda_0 \, t^{-\alpha}$, ($D_t$ in code).
The result is a set of plane patches (polyhedra), their centers, and a contact network adjacency matrix **A**.

## Requirements
- MATLAB R2020a+ (R2021a or newer recommended)

The script loops over **alpha_list × D0_list × run_time** and stops once it adds $N$ planes. For each run, it saves:
	
  •	*_upper_planes.mat — upper_plane_all{t} (upper offset plane vertices)
  
  •	*_down_planes.mat — down_plane_all{t} (lower offset plane vertices)
  
  •	*_planes.mat — planes_coords{t} (mid plane polygon vertices)
  
  •	*_A.mat — adjacency matrix **A** from detectPolyhedronContacts
  
  •	*_final_volume.mat — total volume **V** of all plane polyhedra
  
  •	*_area.mat — per-plane “area” proxy ( $\text{volume}_i / \lambda_i$ ), (in code is $\text{volume}_i / D_i$, $D_i = D_0i^{-\alpha}$ ).

It also creates a PNG visualization (unit cube wireframe + plane patches + selected centroid links).
