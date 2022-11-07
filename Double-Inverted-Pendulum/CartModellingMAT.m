%% Cart Model
syms M g x dx ddx F B
syms n_g K_g n_m k_t V_m r_mp k_m R_m
F = (n_g*K_g*n_m*k_t*(V_m*r_mp-K_g*k_m*dx))/(R_m*(r_mp^2));

ddq = subs((F-B*dx)/M);