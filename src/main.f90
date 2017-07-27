program main
USE Declarations 
implicit none

!============================ Read mesh ======================
call mesh

!=============  Allocate the dimensions to the variables ======================
allocate (Fgrav(N-1,M-1),P(N-1,M-1),alpha(N-1,M-1),nn(N-1,M-1),mm(N-1,M-1))
allocate (K(N-1,M-1),K_south(N-1,M-1),Fpres(N-1,M-1))
allocate (K_east(N-1,M-1),K_west(N-1,M-1),K_north(N-1,M-1),Ks(N-1,M-1),Fcond(N-1,M-1))
allocate (k_eff(N-1,M-1),k_eff_west(N-1,M-1),k_eff_east(N-1,M-1),k_eff_south(N-1,M-1),k_eff_north(N-1,M-1))
allocate (Qin(N-1),outflow(N-1),wrc(N-1,M-1),sphericity(N-1,M-1),dendricity(N-1,M-1))
allocate (Fgravw(N-1, M-1),Fgravn(N-1, M-1),Fgravs(N-1, M-1),grain_opt(N-1,M-1))
allocate (Fprese(N-1, M-1),Fpresw(N-1, M-1),Fpresn(N-1, M-1),Fpress(N-1, M-1),Fgrave(N-1, M-1))
allocate (Fconde(N-1, M-1),Fcondw(N-1, M-1),Fcondn(N-1, M-1),Fconds(N-1, M-1))
allocate (Tss(N-1),Fconv(N-1,M-1), alpha_wetting(N-1,M-1), theta_r(N-1,M-1,100))
allocate (water_content(N-1,M-1,2),grain_classic(N-1,M-1),porosity(N-1,M-1),T(N-1,M-1,2))
allocate (S(N-1,M-1),density(N-1,M-1),dry_density(N-1,M-1),Hn(N-1),rain(N-1), theta_s(N-1,M-1,100))
allocate (Pid(N-1,M-1), Pdi(N-1,M-1), theta_di(N-1,M-1,100), theta_id(N-1,M-1,100), order_scanning(N-1,M-1))

!================= Read model inputs ====================
call read_inputs

Dt = Dt_ini
rain = rain / 3600._dp ! from [m/hr] to [m/s]

!================= Run model ========================
call SMPP



end
