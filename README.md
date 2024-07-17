# kidney_jamming
Code and raw data for Prahl et al. (2024) "Jamming of nephron-forming niches in the developing mouse kidney creates cyclical mechanical stresses". Nature Materials.

Raw data are uploaded as excel files labeled by the relevant figure number.

Code notes:

Figure 1:
1) Open Nov3_E16_reconstructions_example.3dm and curvature_mapping.gh in Rhino and grasshopper within Rhino, respectively. curvature_mapping is a grasshopper workflow that takes a rhino surface as input and computes curvature map.
2) tip_voronoi_example.m generates Voronoi, delaunay, and shape index data and overlays from a csv input of UB tip coordinates (x,y as columns).
3) tip_pos_overlay_example.m generates images of UB tip positions and surface height map from a csv input of UB tip coordinates (x,y as columns) and a .stl file describing the kidney surface.
