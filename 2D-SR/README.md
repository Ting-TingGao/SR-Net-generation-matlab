# 2D Scale-Rich Network Generation (MATLAB)

This project builds a **2D scale-rich (SR) network** inside the unit square $([0,1]^2)$. At each step *t*, it picks a void polygon (probability ∝ area), inserts a **line** of thickness $\lambda_t = \lambda_0\,(t+1)^{-\alpha}$, splits the polygon, and records the new line. Then analyze the resulting network (lengths, degrees, thickness statistics) and produce log–log distributions with power-law fits.


## Requirements
- MATLAB R2020a+ (R2021a or newer recommended)

## Run
To generate 2D SR network structure, run **Main.m**.

## Outputs
  •	Line_list — (t × 12) numeric array;
  
  •	Line_polygon — {t×1} cell of polyshape;
  
  •	meta — struct with options and history (e.g., meta.opts.alpha);
  
  •	A — (t × t) adjacency;
  
  •	Thickness, length, and degree distributions (log-binned) + power-law fits.
  

Tuning Tips
	•	Bins & fit ranges: use bins_* and start_*_fit to control log-binning and the fitted tail.
