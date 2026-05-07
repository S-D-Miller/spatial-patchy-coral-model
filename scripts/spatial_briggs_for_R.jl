using DifferentialEquations, RecursiveArrayTools, LinearAlgebra

C_growth = zeros(size(init_C)[1])
Mi_growth = zeros(size(init_C)[1])
Mi_disp = zeros(size(init_C)[1])
turf = zeros(size(init_C)[1])

p = VectorOfArray([Rc, Rm, pm, Gtc, Gtv, Gti, y, Dc, Dv, Di, w, neighbor_mat_m, neighbor_mat_c, dispersal_mat, C_growth, Mi_growth, Mi_disp, turf])
u0 = vcat(init_C', init_Mv', init_Mi')