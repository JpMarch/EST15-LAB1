using DataFrames, LaTeXStrings, Plots, Unitful, XLSX

Resultados = Dict{Any,Any}()

function nuplot(v::Vector, name::LaTeXString=L"")
    return plot(ustrip.(dados[:,1]),ustrip.(v),label = name)
end

function nuplot!(v::Vector, name::LaTeXString=L"")
    return plot!(ustrip.(dados[:,1]),ustrip.(v),label = name)
end

## Calibração

    # Extensômetro

    R =  120 * u"Ω" # Resistência
    Γ = 2.08        # Gage factor

    # Preencher durante o experimento

    R_c = 9.78e3 * u"Ω"

    ε = -5747e-6

    # Cálculo experimental

    ε_cal = - R / (Γ * (R + R_c))

    σ_rel = (ε_cal - ε) / ε_cal

    Resultados["Calibração"] = Dict{Any,Any}("ε" => ε, "ε_cal" => ε_cal, "σ_rel" => σ_rel)

## Experimento 1

    # Barra

    a =    4.0 * u"inch" # Distância do apoio ao ponto de aplicação da força
    t =  0.255 * u"inch" # Espessura
    b =    1.0 * u"inch" # Largura
    E = 10.5e6 * u"psi"  # Módulo de elasticidade

    g = 9.81 * u"m" * u"s^-2" # Aceleração da gravidade

    # Preencher durante o experimento

    M = 1.5 * u"kg"
    P = M * g

    ε_a =   110e-6
    ε_b = - 115e-6
    ε_c =   225e-6

    # Valores teóricos

    I = b * t ^ 3 / 12 # Momento de inércia da seção transversal

    ε_teo = upreferred(P * a * t / (2 * E * I))

    ε_teo_a = ε_teo
    ε_teo_b = - ε_teo
    ε_teo_c = 2 * ε_teo

    # Cálculo experimental

    σ_rel_a = (ε_a - ε_teo_a) / ε_teo
    σ_rel_b = (ε_b + ε_teo_b) / ε_teo
    σ_rel_c = (ε_c - ε_teo_c) / ε_teo

    Resultados["Experimento 1"] = Dict{Any,Any}(     "ε_a" => ε_a,         "ε_b" => ε_b,         "ε_c" => ε_c,
                                                 "ε_teo_a" => ε_teo_a, "ε_teo_b" => ε_teo_b, "ε_teo_c" => ε_teo_c,
                                                 "σ_rel_a" => σ_rel_a, "σ_rel_b" => σ_rel_b, "σ_rel_c" => σ_rel_c)

## Experimento 2

    # Barra

    d =   0.75 * u"inch" # Diâmetro da barra
    E = 10.5e6 * u"psi"  # Módulo de elasticidade
    L =     10 * u"inch" # Distância entre os extensômetros

    # Preenchercher durante o experimento

    M_real = 2.154 * u"kg"
    ε = 217e-6

    # Valores teóricos

    I = π * d ^ 4 / 64 # Momento de inércia da seção transversal da barra

    # Cálculo exsperimental

    P = uconvert(u"N", ε * E * I / (L * d))
    M_medida = upreferred(P/g)
    σ_rel = (M_real - M_medida)/M_real

    Resultados["Experimento 2"] = Dict{Any,Any}(  "M_real" => M_real, "ε" => ε, "P" => P,
                                                "M_medida" => M_medida, "σ_rel" => σ_rel)

## Experimento 3

    # Cilindro

    d =   3.5 * u"inch" # Diâmetro interno
    e =  0.12 * u"inch" # Espessura da parede
    E =  29e6 * u"psi"  # Módulo de elasticidade
    ν = 0.283           # Coeficiente de Poisson
    a = d / 2           # Raio interno
    b = a + e           # Raio externo

    # Preencher durante o experimento

    θ_a = 25
    θ_b = θ_a + 45
    θ_c = θ_b + 45

    M = [
            (cosd(θ_a))^2 (sind(θ_a))^2 sind(θ_a)*cosd(θ_a)
            (cosd(θ_b))^2 (sind(θ_b))^2 sind(θ_b)*cosd(θ_b)
            (cosd(θ_c))^2 (sind(θ_c))^2 sind(θ_c)*cosd(θ_c)
    ]

    dados = DataFrame(XLSX.readdata((@__DIR__)*"/dados.xls","Plan1","B50:E349"),[:P,:ε_a,:ε_b,:ε_c])

    plt1 = nuplot(dados[:,2],L"\varepsilon_a")
           nuplot!(dados[:,3],L"\varepsilon_b")
           nuplot!(dados[:,4],L"\varepsilon_c")
           xlabel!("P (psi)")
           ylabel!(L"\mathbf{\varepsilon}")

    ε_1_vec = Vector{Any}(undef,length(dados[:,1]))
    ε_2_vec = Vector{Any}(undef,length(dados[:,1]))

    σ_1_vec = Vector{Any}(undef,length(dados[:,1]))
    σ_2_vec = Vector{Any}(undef,length(dados[:,1]))
    σ_t_vec = Vector{Any}(undef,length(dados[:,1]))
    σ_l_vec = Vector{Any}(undef,length(dados[:,1]))

    P_vec = Vector{Any}(undef,length(dados[:,1]))

    i_370 = 1

    for i in 1:length(dados[:,1])

        P   = dados[i,1] * u"psi"
        if abs(ustrip(P) - 370) < abs(dados[i_370,1] - 370)
            i_370 = i
        end
        P_vec[i] = P

        ε_a = dados[i,2] * 1e-6
        ε_b = dados[i,3] * 1e-6
        ε_c = dados[i,4] * 1e-6

        # Valores teóricos

        σ_t = 2 * a ^ 2 * P / (b ^ 2 - a ^ 2)
        σ_l = σ_t / 2

        σ_t_vec[i] = σ_t; σ_l_vec[i] = σ_l

        # Cálculo experimental

        ε_x, ε_y, γ_xy = M\[ε_a,ε_b,ε_c]

        ε_1 = (ε_x + ε_y) / 2 + sqrt(((ε_x - ε_y) / 2) ^ 2 + (γ_xy / 2) ^ 2)
        ε_2 = (ε_x + ε_y) / 2 - sqrt(((ε_x - ε_y) / 2) ^ 2 + (γ_xy / 2) ^ 2)

        ε_1_vec[i] = ε_1; ε_2_vec[i] = ε_2

        σ_1 = (E / (1 - ν ^ 2)) * (ε_1 + ν * ε_2)
        σ_2 = (E / (1 - ν ^ 2)) * (ε_2 + ν * ε_1)

        σ_1_vec[i] = σ_1; σ_2_vec[i] = σ_2

    end

    σ_rel_1 = (σ_1_vec - σ_t_vec) ./ σ_t_vec
    σ_rel_2 = (σ_2_vec - σ_l_vec) ./ σ_l_vec

    plt2 = nuplot(σ_t_vec,L"\sigma_1\; \textrm{teórico}")
           nuplot!(σ_1_vec,L"\sigma_1\; \textrm{experimental}")
           nuplot!(σ_l_vec,L"\sigma_2\; \textrm{experimental}")
           nuplot!(σ_2_vec,L"\sigma_2\; \textrm{experimental}")
           xlabel!("P (psi)")
           ylabel!("σ (psi)")

    P_370_m = P_vec[i_370]
    ε_a_370_m = dados.ε_a[i_370]
    ε_b_370_m = dados.ε_b[i_370]
    ε_c_370_m = dados.ε_c[i_370]
    ε_1_370_m = ε_1_vec[i_370]
    ε_2_370_m = ε_2_vec[i_370]
    σ_1_370_m = σ_1_vec[i_370]
    σ_t_370_m = σ_t_vec[i_370]
    σ_2_370_m = σ_2_vec[i_370]
    σ_l_370_m = σ_l_vec[i_370]

    P_370_p = P_vec[i_370 + 1]
    ε_a_370_p = dados.ε_a[i_370 + 1]
    ε_b_370_p = dados.ε_b[i_370 + 1]
    ε_c_370_p = dados.ε_c[i_370 + 1]
    ε_1_370_p = ε_1_vec[i_370 + 1]
    ε_2_370_p = ε_2_vec[i_370 + 1]
    σ_1_370_p = σ_1_vec[i_370 + 1]
    σ_t_370_p = σ_t_vec[i_370 + 1]
    σ_2_370_p = σ_2_vec[i_370 + 1]
    σ_l_370_p = σ_l_vec[i_370 + 1]

    ε_a_370 = ε_a_370_m + (370*u"psi" - P_370_m)/(P_370_p - P_370_m) * (ε_a_370_p - ε_a_370_m)
    ε_b_370 = ε_b_370_m + (370*u"psi" - P_370_m)/(P_370_p - P_370_m) * (ε_b_370_p - ε_b_370_m)
    ε_c_370 = ε_c_370_m + (370*u"psi" - P_370_m)/(P_370_p - P_370_m) * (ε_c_370_p - ε_c_370_m)
    ε_1_370 = ε_1_370_m + (370*u"psi" - P_370_m)/(P_370_p - P_370_m) * (ε_1_370_p - ε_1_370_m)
    ε_2_370 = ε_2_370_m + (370*u"psi" - P_370_m)/(P_370_p - P_370_m) * (ε_2_370_p - ε_2_370_m)
    σ_1_370 = σ_1_370_m + (370*u"psi" - P_370_m)/(P_370_p - P_370_m) * (σ_1_370_p - σ_1_370_m)
    σ_t_370 = σ_t_370_m + (370*u"psi" - P_370_m)/(P_370_p - P_370_m) * (σ_t_370_p - σ_t_370_m)
    σ_2_370 = σ_2_370_m + (370*u"psi" - P_370_m)/(P_370_p - P_370_m) * (σ_2_370_p - σ_2_370_m)
    σ_l_370 = σ_l_370_m + (370*u"psi" - P_370_m)/(P_370_p - P_370_m) * (σ_l_370_p - σ_l_370_m)
    
    σ_rel_1 = (σ_t_370 - σ_1_370) / σ_t_370
    σ_rel_2 = (σ_l_370 - σ_2_370) / σ_l_370

    Resultados["Experimento 3"] = Dict{Any,Any}("Deformações" => plt1, "Tensões" => plt2,
                                                # "P_370_m" => P_370_m, "P_370_p" => P_370_p,
                                                # "σ_1_370_m" => σ_1_370_m, "σ_1_370_p" => σ_1_370_p,
                                                # "σ_t_370_m" => σ_t_370_m, "σ_t_370_p" => σ_t_370_p,
                                                # "σ_2_370_m" => σ_2_370_m, "σ_2_370_p" => σ_2_370_p,
                                                # "σ_l_370_m" => σ_l_370_m, "σ_l_370_p" => σ_l_370_p,
                                                  "σ_1_370" => σ_1_370,     "σ_t_370" => σ_t_370,
                                                  "σ_2_370" => σ_2_370,     "σ_l_370" => σ_l_370,
                                                  "σ_rel_1" => σ_rel_1,     "σ_rel_2" => σ_rel_2)

    for key in ["Calibração","Experimento 1", "Experimento 2","Experimento 3"]
        exp = Resultados[key]
        println(key,"\n")
        for key1 in keys(exp)
            if (key1 != "Tensões") && (key1 != "Deformações")
                println(key1," = ",exp[key1])
            else
                display(exp[key1])
            end
        end
        println("\n")
    end
