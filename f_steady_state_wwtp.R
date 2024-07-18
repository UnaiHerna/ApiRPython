
# Start of code ------------------------------------------------------------------------------------

f_steady_state_wwtp <- function(mltss_sp, so_aer_sp, q_int, tss_eff_sp, temp) {
    
    # Check input argument mltss_sp
    is_ok <- checkmate::check_double(x = mltss_sp, # set-point of MLTSS in aerobic reactor
                                     lower = 2500,
                                     upper = 4000,
                                     finite = TRUE,
                                     any.missing = FALSE,
                                     len = 1L,
                                     null.ok = FALSE)
    if (!is.logical(x = is_ok)) {
        return(paste0("Error in input argument mltss_sp. ", is_ok, "."))
    }
    remove(is_ok)
    
    # Check input argument so_aer_sp
    is_ok <- checkmate::check_double(x = so_aer_sp, # set-point of dossilved O2 in aerobic reactor
                                     lower = 0.5,
                                     upper = 3.5,
                                     finite = TRUE,
                                     any.missing = FALSE,
                                     len = 1L,
                                     null.ok = FALSE)
    if (!is.logical(x = is_ok)) {
        return(paste0("Error in input argument so_aer_sp. ", is_ok, "."))
    }
    remove(is_ok)
    
    # Check input argument q_int
    is_ok <- checkmate::check_double(x = q_int, # internal flowrate
                                     lower = 20648,
                                     upper = 20648 * 5,
                                     finite = TRUE,
                                     any.missing = FALSE,
                                     len = 1L,
                                     null.ok = FALSE)
    if (!is.logical(x = is_ok)) {
        return(paste0("Error in input argument q_int ", q_int, "."))
    }
    remove(is_ok)
    
    # Check input argument tss_eff_sp
    is_ok <- checkmate::check_double(x = tss_eff_sp, # TSS in effluent
                                     lower = 0,
                                     upper = 20,
                                     finite = TRUE,
                                     any.missing = FALSE,
                                     len = 1L,
                                     null.ok = FALSE)
    if (!is.logical(x = is_ok)) {
        return(paste0("Error in input argument tss_eff_sp ", tss_eff_sp, "."))
    }
    remove(is_ok)
    
    # Check input argument temp
    is_ok <- checkmate::check_double(x = temp, # water temperature
                                     lower = 13,
                                     upper = 25,
                                     finite = TRUE,
                                     any.missing = FALSE,
                                     len = 1L,
                                     null.ok = FALSE)
    if (!is.logical(x = is_ok)) {
        return(paste0("Error in input argument temp ", temp, "."))
    }
    remove(is_ok)
    
    # Load functions
    f_ode_wwtp <- function(t, y, pars) {
        
        f_activated_sludge_unit <- function(v_conc, v_conc_inflow, v_param) {
            
            # State variables: g/m3
            salk <- v_conc[1L]
            si <- v_conc[2L]
            snd <- v_conc[3L]
            snh <- v_conc[4L]
            sno <- v_conc[5L]
            so <- v_conc[6L]
            ss <- v_conc[7L]
            xba <- v_conc[8L]
            xbh <- v_conc[9L]
            xi <- v_conc[10L]
            xnd <- v_conc[11L]
            xp <- v_conc[12L]
            xs <- v_conc[13L]
            
            # Inflow parameters
            flowrate <- unname(obj = v_param["flowrate"])
            
            # Hydraulic parameters
            vol <- unname(obj = v_param["vol"])
            
            # Aeration parameters
            so_sat <- unname(obj = v_param["so_sat"])
            kla <- unname(obj = v_param["kla"])
            
            # Stoichiometric parameters
            y_a <- unname(obj = v_param["y_a"])
            y_h <- unname(obj = v_param["y_h"])
            f_p <- unname(obj = v_param["f_p"])
            
            # Composition parameters
            i_x_b <- unname(obj = v_param["i_x_b"])
            i_x_p <- unname(obj = v_param["i_x_p"])
            
            # Kinetic parameters
            k_nh_a <- unname(obj = v_param["k_nh_a"])
            k_nh_h <- unname(obj = v_param["k_nh_h"])
            k_no <- unname(obj = v_param["k_no"])
            k_oa <- unname(obj = v_param["k_oa"])
            k_oh <- unname(obj = v_param["k_oh"])
            k_s <- unname(obj = v_param["k_s"])
            k_x_temp <- unname(obj = v_param["k_x_temp"])
            b_a_temp <- unname(obj = v_param["b_a_temp"])
            b_h_temp <- unname(obj = v_param["b_h_temp"])
            k_a_temp <- unname(obj = v_param["k_a_temp"])
            k_h_temp <- unname(obj = v_param["k_h_temp"])
            mu_a_temp <- unname(obj = v_param["mu_a_temp"])
            mu_h_temp <- unname(obj = v_param["mu_h_temp"])
            n_g <- unname(obj = v_param["n_g"])
            n_h <- unname(obj = v_param["n_h"])
            remove(v_param)
            
            # Monod equations: 1
            monod_snh_a_act <- snh / (k_nh_a + snh)
            monod_snh_h_act <- snh / (k_nh_h + snh)
            monod_sno_act <- sno / (k_no + sno)
            monod_so_h_act <- so / (k_oh + so)
            monod_so_a_act <- so / (k_oa + so)
            monod_so_des <- k_oh / (k_oh + so)
            monod_ss_act <- ss / (k_s + ss)
            monod_xs_act <- xs / (k_x_temp * xbh + xs)
            remove(k_nh_a, k_nh_h, k_no, k_oa, k_oh, k_s, k_x_temp)
            
            # Process rates: g/(m3 d)
            # 1. Aerobic growth of Xbh
            # 2. Anoxic growth of Xbh
            # 3. Aerobic growth of Xba
            # 4. Decay of Xbh
            # 5. Decay of Xba
            # 6. Ammonification of soluble organic nitrogen
            # 7. Hydrolisis of entraped organics
            # 8. Hydrolisis of entraped organic nitrogen
            # 9. Aeration
            rho_1 <- mu_h_temp * monod_ss_act * monod_so_h_act * monod_snh_h_act * xbh
            rho_2 <- mu_h_temp * monod_ss_act * monod_so_des * monod_sno_act * monod_snh_h_act * n_g *
                xbh
            rho_3 <- mu_a_temp * monod_snh_a_act * monod_so_a_act * xba
            rho_4 <- b_h_temp * xbh
            rho_5 <- b_a_temp * xba
            rho_6 <- k_a_temp * snd * xbh
            rho_7 <- k_h_temp * monod_xs_act *
                (monod_so_h_act + n_h * monod_so_des * monod_sno_act) * xbh
            rho_8 <- rho_7 * xnd / xs
            rho_9 <- kla * (so_sat - so)
            remove(monod_snh_a_act, monod_snh_h_act, monod_sno_act, monod_so_h_act, monod_so_a_act,
                   monod_so_des, monod_ss_act, monod_xs_act)
            remove(b_a_temp, b_h_temp, k_a_temp, k_h_temp, mu_h_temp, n_g, n_h)
            remove(kla, so_sat)
            remove(salk, si, snd, snh, sno, so, ss, xba, xbh, xi, xnd, xp, xs)
            
            # ODEs for transformations: g/(m3 d)
            v_ode_transf <- vector(mode = "double", length = length(x = v_conc)) * NA_real_
            # Salk
            v_ode_transf[1L] <- - i_x_b/14 * rho_1 +
                ((1 - y_h)/(14 * 2.86 * y_h) - i_x_b/14) * rho_2 +
                (-i_x_b/14 - 1/(7*y_a)) * rho_3 +
                1/14 * rho_6
            # Si
            v_ode_transf[2L] <- 0
            # Snd
            v_ode_transf[3L] <- - rho_6 +
                rho_8
            # Snh
            v_ode_transf[4L] <- - i_x_b * rho_1 +
                (- i_x_b) * rho_2 +
                (- i_x_b - 1/y_a) * rho_3 + 
                rho_6
            # Sno
            v_ode_transf[5L] <- (-(1 - y_h)/(2.86 * y_h)) * rho_2 +
                1/y_a * rho_3
            # So
            v_ode_transf[6L] <- (-(1 - y_h)/y_h) * rho_1 +
                (-(4.57 - y_a)/y_a) * rho_3 +
                rho_9
            # Ss
            v_ode_transf[7L] <- -1/y_h * rho_1 +
                (-1/y_h) * rho_2 +
                rho_7
            # Xba
            v_ode_transf[8L] <- rho_3 +
                (-rho_5)
            # Xbh
            v_ode_transf[9L] <- rho_1 +
                rho_2 +
                (-rho_4)
            # Xi
            v_ode_transf[10L] <- 0
            # Xnd
            v_ode_transf[11L] <- (i_x_b - f_p * i_x_p) * rho_4 +
                (i_x_b - f_p * i_x_p) * rho_5 +
                (-rho_8)
            # Xp
            v_ode_transf[12L] <- f_p  * rho_4 +
                f_p  * rho_5
            # Xs
            v_ode_transf[13L] <- (1 - f_p) * rho_4 +
                (1 - f_p) * rho_5 + 
                (-rho_7)
            remove(rho_1, rho_2, rho_3, rho_4, rho_5, rho_6, rho_7, rho_8, rho_9)
            remove(f_p, i_x_b, i_x_p, y_a, y_h)
            
            # ODEs for transportation: g/(m3 d)
            v_ode_transp <- flowrate/vol * (v_conc_inflow - v_conc)
            remove(flowrate, vol, v_conc_inflow, v_conc)
            
            # ODEs total: g/(m3 d)
            v_ode <- v_ode_transf + v_ode_transp
            v_ode <- unname(obj = v_ode)
            remove(v_ode_transf, v_ode_transp)
            
            # End function
            v_ode
            
        }
        
        # Time not used
        remove(t)
        
        # State variables. Without names works faster
        # 1L:13L: anoxic reactor: salk, si, snd, snh, sno, so, ss, xba, xbh, xi, xnd, xp, xs
        # 14L:26L: aerobic reactor: salk, si, snd, snh, sno, so, ss, xba, xbh, xi, xnd, xp, xs
        # 27: integral of the error of the MLTSS control
        v_conc_anx <- unname(obj = y[1L:13L])
        v_conc_aer <- unname(obj = y[14L:26L])
        ie_mltss <- unname(obj = y[27L])
        ie_kla_aer <- unname(obj = y[28L])
        remove(y)
        
        # Control of MLTSS in aerobic reactor
        mltss_sp <- unname(obj = pars["mltss_sp"]) # set-point
        tss_pcod <- unname(obj = pars["tss_pcod"])
        mltss <- sum(v_conc_aer[8L:13L]) * tss_pcod
        error <- mltss_sp - mltss
        q_was_0 <- unname(obj = pars["q_was_0"]) # m3/d
        kp_q_was <- unname(obj = pars["kp_q_was"]) # m6/(gSS d)
        ti_q_was <- unname(obj = pars["ti_q_was"]) # d
        q_was <- q_was_0 + kp_q_was * error + kp_q_was/ti_q_was * ie_mltss
        if (q_was > 5000) {
            q_was <- 5000
        } else if (q_was < 0) {
            q_was <- 0
        }
        ode_ie_mltss <- error
        remove(mltss_sp, mltss, error, q_was_0, kp_q_was, ti_q_was, ie_mltss)
        
        # Flowrates (sum inputs == sum outputs)
        # Combiner: in (q_inf, q_int, q_ext); out (q_combiner)
        # Anoxic reactor and aerobic reactor: in and out (q_combiner)
        # Splitter 1: in (q_combiner); out (q_int, q_splitter)
        # Secondary settler: in (q_splitter); out (q_eff, q_under)
        # Splitter 2: in (q_under); out (q_was, q_ext) 
        
        # Calculate flowrates
        q_inf <- unname(obj = pars["q_inf"]) # influent flowrate
        q_int <- unname(obj = pars["q_int"]) # internal flowrate
        q_under <- unname(obj = pars["q_under"]) # flowrate extracted from secondary settler
        q_combiner <- q_inf + q_int + q_under - q_was
        q_splitter_1 <- q_inf + q_under - q_was # internal recirculation
        q_splitter_2 <- q_under - q_was # wastage
        q_eff <- q_inf - q_was # effluent
        
        # Secondary settler
        tss_eff_sp <- unname(obj = pars["tss_eff_sp"])
        fns_aux <- 0.005 # non-sedimentable PCOD
        tss_aux <- sum(v_conc_aer[8L:13L]) * fns_aux * tss_pcod
        f_ns <- fns_aux * tss_eff_sp / tss_aux
        v_conc_eff <- v_conc_aer # solubles (from 1L:7L)
        v_conc_eff[8L:13L] <- v_conc_aer[8L:13L] * f_ns # particulates (from 8L:13L)
        tss_eff <- sum(v_conc_eff[8L:13L]) * tss_pcod
        v_conc_under <- (v_conc_aer * q_splitter_1 - v_conc_eff * q_eff)/q_under
        remove(tss_eff_sp, fns_aux, tss_aux, f_ns, q_splitter_1, q_under)
        
        # Effluent
        # v_eff <- c(q_eff, v_conc_eff)
        remove(q_eff, v_conc_eff)
        
        # Wastage
        sludge_prod <- sum(v_conc_under[8L:13L]) * tss_pcod * q_was / 1000 # kgSS/d
        remove(tss_pcod, q_was)
        
        # Influent
        v_conc_inf <- c(salk = unname(obj = pars["salk_inf"]),
                        si = unname(obj = pars["si_inf"]),
                        snd = unname(obj = pars["snd_inf"]),
                        snh = unname(obj = pars["snh_inf"]),
                        sno = unname(obj = pars["sno_inf"]),
                        so = unname(obj = pars["so_inf"]),
                        ss = unname(obj = pars["ss_inf"]),
                        xba = unname(obj = pars["xba_inf"]),
                        xbh = unname(obj = pars["xbh_inf"]),
                        xi = unname(obj = pars["xi_inf"]),
                        xnd = unname(obj = pars["xnd_inf"]),
                        xp = unname(obj = pars["xp_inf"]),
                        xs = unname(obj = pars["xs_inf"]))
        
        # Combiner 1: influent + internal recirculation + external recirculation
        v_conc_combiner <- (v_conc_inf * q_inf + v_conc_aer * q_int + v_conc_under * q_splitter_2) / 
            q_combiner
        remove(q_inf, q_int, q_splitter_2, v_conc_inf)
        
        # Anoxic reactor
        v_param_anx <- c(
            # Inflow
            flowrate = q_combiner,
            # Hydraulic parameters
            vol = unname(obj = pars["vol_anx"]),
            # Aeration
            so_sat = unname(obj = pars["so_sat"]),
            kla = unname(obj = pars["kla_anx"]),
            # Stoichiometric parameters
            y_a = unname(obj = pars["y_a"]),
            y_h = unname(obj = pars["y_h"]),
            f_p = unname(obj = pars["f_p"]),
            # Composition parameters
            i_x_b = unname(obj = pars["i_x_b"]),
            i_x_p = unname(obj = pars["i_x_p"]),
            # Kinetic parameters
            k_nh_a = unname(obj = pars["k_nh_a"]),
            k_nh_h = unname(obj = pars["k_nh_h"]),
            k_no = unname(obj = pars["k_no"]),
            k_oa = unname(obj = pars["k_oa"]),
            k_oh = unname(obj = pars["k_oh"]),
            k_s = unname(obj = pars["k_s"]),
            k_x_temp = unname(obj = pars["k_x_temp"]),
            b_a_temp = unname(obj = pars["b_a_temp"]),
            b_h_temp = unname(obj = pars["b_h_temp"]),
            k_a_temp = unname(obj = pars["k_a_temp"]),
            k_h_temp = unname(obj = pars["k_h_temp"]),
            mu_a_temp = unname(obj = pars["mu_a_temp"]),
            mu_h_temp = unname(obj = pars["mu_h_temp"]),
            n_g = unname(obj = pars["n_g"]),
            n_h = unname(obj = pars["n_h"]))
        v_ode_anx <- f_activated_sludge_unit(v_conc = v_conc_anx,
                                             v_conc_inflow = v_conc_combiner,
                                             v_param = v_param_anx)
        remove(v_conc_combiner, v_param_anx)
        
        # Control of kla in aerobic reactor
        so_aer_sp <- unname(obj = pars["so_aer_sp"]) # set-point
        error <- so_aer_sp - v_conc_aer[6L]
        kla_aer_0 <- unname(obj = pars["kla_aer_0"]) # 1/d
        kp_kla_aer <- unname(obj = pars["kp_kla_aer"]) # m3/(gO2 d)
        ti_kla_aer <- unname(obj = pars["ti_kla_aer"]) # d
        kla_aer <- kla_aer_0 + kp_kla_aer * error + kp_kla_aer/ti_kla_aer * ie_kla_aer
        if (kla_aer > 5000) {
            kla_aer <- 5000
        } else if (kla_aer < 0) {
            kla_aer <- 0
        }
        ode_ie_kla_aer <- error
        remove(so_aer_sp, error, kla_aer_0, kp_kla_aer, ti_kla_aer, ie_kla_aer)
        
        # Aerobic reactor
        v_param_aer <- c(
            # Inflow
            flowrate = q_combiner,
            # Hydraulic parameters
            vol = unname(obj = pars["vol_aer"]),
            # Aeration
            so_sat = unname(obj = pars["so_sat"]),
            kla = kla_aer,
            # Stoichiometric parameters
            y_a = unname(obj = pars["y_a"]),
            y_h = unname(obj = pars["y_h"]),
            f_p = unname(obj = pars["f_p"]),
            # Composition parameters
            i_x_b = unname(obj = pars["i_x_b"]),
            i_x_p = unname(obj = pars["i_x_p"]),
            # Kinetic parameters
            k_nh_a = unname(obj = pars["k_nh_a"]),
            k_nh_h = unname(obj = pars["k_nh_h"]),
            k_no = unname(obj = pars["k_no"]),
            k_oa = unname(obj = pars["k_oa"]),
            k_oh = unname(obj = pars["k_oh"]),
            k_s = unname(obj = pars["k_s"]),
            k_x_temp = unname(obj = pars["k_x_temp"]),
            b_a_temp = unname(obj = pars["b_a_temp"]),
            b_h_temp = unname(obj = pars["b_h_temp"]),
            k_a_temp = unname(obj = pars["k_a_temp"]),
            k_h_temp = unname(obj = pars["k_h_temp"]),
            mu_a_temp = unname(obj = pars["mu_a_temp"]),
            mu_h_temp = unname(obj = pars["mu_h_temp"]),
            n_g = unname(obj = pars["n_g"]),
            n_h = unname(obj = pars["n_h"]))
        v_ode_aer <- f_activated_sludge_unit(v_conc = v_conc_aer,
                                             v_conc_inflow = v_conc_anx,
                                             v_param = v_param_aer)
        remove(v_conc_anx, v_conc_aer, v_param_aer, pars, q_combiner, f_activated_sludge_unit)
        
        # System of ODEs. The order should be the same as "y"
        v_ode <- c(v_ode_anx, v_ode_aer, ode_ie_mltss, ode_ie_kla_aer)
        remove(v_ode_anx, v_ode_aer, ode_ie_mltss, ode_ie_kla_aer)
        
        # mean absolute rate of change to evaluate if steady state has been reached (rootSolve::stode)
        precis <- mean(x = abs(x = v_ode))
        
        # End function
        list(v_ode,
             precis = precis,
             kla_aer = kla_aer,
             tss_eff = tss_eff,
             sludge_prod = sludge_prod)
        
    }
    
    # Influent
    q_inf <- 20648L # m3/d
    salk_inf <- 200 # molar units: alkalinity
    si_inf <- 35.2 # gCOD/m3: soluble inert organic matter
    snd_inf <- 5.28 # gN/m3: soluble biodegradable organic nitrogen
    snh_inf <- 25.6 # gN/m3: NH4 + NH3 nitrogen
    sno_inf <- 0.01 # gN/m3: nitrate and nitrite nitrogen
    so_inf <- 0.01 # gO2/m3: oxygen
    ss_inf <- 88 # gCOD/m3: readily biodegradable substrate
    xba_inf <- 0.01 # gCOD/m3: active autotrophic biomass
    xbh_inf <- 0.01 # gCOD/m3: active heterotrophic biomass
    xi_inf <- 77.04 # gCOD/m3: particulate inert organic matter
    xnd_inf <- 6.12 # gN/m3: particulate biodegradable organic nitrogen
    xp_inf <- 0.01 # gCOD/m3: particulate products arising from biomass decay
    xs_inf <- 179.76 # gCOD/m3: slowly biodegradable substrate
    
    # Hydraulic parameters
    vol_anx <- 2600L # m3
    vol_aer <- 5250L # m3
    
    # Stoichiometric parameters
    y_a <- 0.24 # gCOD/gN, yield of autotrophic biomass
    y_h <- 0.67 # gCOD/gCOD, yield of heterotrophic biomass
    f_p <- 0.08 # 1, fraction of biomass converted to inert matter
    
    # Composition parameters
    i_x_b <- 0.086 # gN/gCOD, mass of nitrogen per mass of COD in biomass
    i_x_p <- 0.06 # gN/gCOD, mass of nitrogen per mass of COD in products formed
    
    # Kinetic parameters at 20 ºC
    k_nh_a <- 1 # gNH3-N/m3, ammonia half-saturation coefficient For autotrophic biomass
    k_nh_h <- 0.05 # g/m3, saturation coefficient of heterotrophs for ammonium
    k_no <- 0.5 # gNO3-N/m3, nitrate half-saturation coefficient for denitrifying Xbh
    k_oa <- 0.4 # gO2/m3, oxygen half-saturation coefficient for autotrophic biomass
    k_oh <- 0.2 # gO2/m3, oxygen half-saturation coefficient for heterotrophic biomass
    k_s <- 20 # gCOD/m3, half-saturation coefficient for heterotrophic biomass
    k_x <- 0.03 # gCOD/gCOD, half saturation coeff. for hydrolysis of slowly biodegradable substrate
    b_a <- 0.15 # 1/d, decay coefficient for autotrophic biomass
    b_h <- 0.62 # 1/d, decay coefficient for heterotrophic biomass
    k_a <- 0.08 # m3/(gCOD*d), maximum specific ammonification rate
    k_h <- 3 # gCOD/(gCOD*d), maximum specific hydrolysis rate
    mu_a <- 0.8 # 1/d, maximum specific growth rate for autotrophic biomass
    mu_h <- 6 # 1/d, maximum specific growth rate for heterotrophic biomass
    n_g <- 0.8 # 1, correction factor for anoxic growth of heteritrophs
    n_h <- 0.4 # 1, correction factor for anoxic hydrolysis
    
    # Temperature correction factors. Unit: 1. theta = k_T/k_Tref^(1/(T-Tref))
    theta_k_x <- 1.116
    theta_b_a <- 1.116
    theta_b_h <- 1.12
    theta_k_a <- 1.072
    theta_k_h <- 1.116
    theta_mu_a <- 1.103
    theta_mu_h <- 1.072
    
    # Kinetic parameters at T
    k_x_temp <- k_x * theta_k_x^(temp - 20) # gCOD/gCOD
    b_a_temp <- b_a * theta_b_a^(temp - 20) # 1/d
    b_h_temp <- b_h * theta_b_h^(temp - 20) # 1/d
    k_a_temp <- k_a * theta_k_a^(temp - 20) # m3/(gCOD*d)
    k_h_temp <- k_h * theta_k_h^(temp - 20) # gCOD/(gCOD*d)
    mu_a_temp <- mu_a * theta_mu_a^(temp - 20) # 1/d
    mu_h_temp <- mu_h * theta_mu_h^(temp - 20) # 1/d
    remove(k_x, b_a, b_h, k_a, k_h, mu_a, mu_h,
           theta_k_x, theta_b_a, theta_b_h, theta_k_a, theta_k_h, theta_mu_a, theta_mu_h)
    
    # Aeration
    temp_0 <- 273.15 # ºC
    temp_k <- temp + temp_0
    so_sat <- 290326 * exp(x = -66.7354 +
                               87.4755/(temp_k / 100) +
                               24.4526 * log(x = temp_k / 100, base = exp(x = 1L))) # WEST formula
    remove(temp_k, temp_0, temp)
    
    # Vector of parameters
    v_pars <- c(
        # Influent
        q_inf = q_inf,
        salk_inf = salk_inf,
        si_inf = si_inf,
        snd_inf = snd_inf,
        snh_inf = snh_inf,
        sno_inf = sno_inf,
        so_inf = so_inf,
        ss_inf = ss_inf,
        xba_inf = xba_inf,
        xbh_inf = xbh_inf,
        xi_inf = xi_inf,
        xnd_inf = xnd_inf,
        xp_inf = xp_inf,
        xs_inf = xs_inf,
        # Hydraulic parameters
        vol_anx = vol_anx,
        vol_aer = vol_aer,
        # Stoichiometric parameters
        y_a = y_a,
        y_h = y_h,
        f_p = f_p,
        # Composition parameters
        i_x_b = i_x_b,
        i_x_p = i_x_p,
        # Kinetic parameters
        k_nh_a = k_nh_a,
        k_nh_h = k_nh_h,
        k_no = k_no,
        k_oa = k_oa,
        k_oh = k_oh,
        k_s = k_s,
        k_x_temp = k_x_temp,
        b_a_temp = b_a_temp,
        b_h_temp = b_h_temp,
        k_a_temp = k_a_temp,
        k_h_temp = k_h_temp,
        mu_a_temp = mu_a_temp,
        mu_h_temp = mu_h_temp,
        n_g = n_g,
        n_h = n_h,
        # Aeration
        so_sat = so_sat,
        kla_anx = 3,
        # Internal recirculation
        q_int = q_int, # m3
        # Secondary settler
        q_under = q_inf, # m3
        # Effluent
        tss_eff_sp = tss_eff_sp, # g SS/m3
        tss_pcod = 0.75, # g SS/g COD
        # Control MLTSS
        mltss_sp = mltss_sp, # g SS/m3
        q_was_0 = 358.0192, # m3/d
        kp_q_was = -0.015, # m6/(gSS d)
        ti_q_was = 1, # d
        # Control OD in aerobic reactor
        so_aer_sp = so_aer_sp, # gO2/m3
        kla_aer_0 = 149.7310, # m3/d
        kp_kla_aer = 50, # m6/(gSS d)
        ti_kla_aer = 0.1) # d
    remove(q_inf, salk_inf, si_inf, snd_inf, snh_inf, sno_inf, so_inf, ss_inf, xba_inf, xbh_inf,
           xi_inf, xnd_inf, xp_inf, xs_inf,
           vol_anx, vol_aer,
           y_a, y_h, f_p, i_x_b, i_x_p, k_nh_a, k_nh_h, k_no, k_oa, k_oh, k_s, k_x_temp, b_a_temp,
           b_h_temp, k_a_temp, k_h_temp, mu_a_temp, mu_h_temp, n_g, n_h,
           so_sat)
    remove(mltss_sp, so_aer_sp, tss_eff_sp, q_int)
    
    # Initial values of state variables
    v_y0 <- c(
        # Anoxic tank
        "salk_anx" = 198.836587557366,
        "si_anx" = 35.2,
        "snd_anx" = 0.629074186863072,
        "snh_anx" = 9.96600864422233,
        "sno_anx" = 0.442517574623713,
        "so_anx" = 0.00891858439206792,
        "ss_anx" = 5.99271552562978,
        "xba_anx" = 85.827419398825,
        "xbh_anx" = 1971.00150750895,
        "xi_anx" = 2000.20358007567,
        "xnd_anx" = 4.92033971386844,
        "xp_anx" = 553.09231231093,
        "xs_anx" = 89.974348762338,
        # Aerobic tank
        "salk_aer" = 197.749977589894,
        "si_aer" = 35.2,
        "snd_aer" = 0.820141779069356,
        "snh_aer" = 2.06620444843421,
        "sno_aer" = 7.97651818992103,
        "so_aer" = 2,
        "ss_aer" = 2.71848921856332,
        "xba_aer" = 87.2292631528311,
        "xbh_aer" = 1993.47111696308,
        "xi_aer" = 2000.20358007578,
        "xnd_aer" = 1.76033897135142,
        "xp_aer" = 557.926622904095,
        "xs_aer" = 26.0757446129945,
        # Control (integral of error)
        "ie_mltss" = 0,
        "ie_kla_aer" = 0)
    
    # Run simulation
    m_sim_vode <- deSolve::vode(y = v_y0,
                                times = c(0, 1000),
                                func = f_ode_wwtp,
                                parms = v_pars,
                                rtol = 1e-6, # def: 1e-6
                                atol = 1e-6, # def: 1e-6
                                jactype = "fullint", # def: "fullint"
                                verbose = FALSE,
                                ynames = FALSE, # FALSE to increase speed
                                maxsteps = 5000L) # def: 5000L
    remove(v_y0, f_ode_wwtp)
    
    # Take results
    v_conc_anx <- m_sim_vode[nrow(x = m_sim_vode), 2L:14L] # last instant of time
    v_conc_aer <- m_sim_vode[nrow(x = m_sim_vode), 15L:27L] # last instant of time
    mltss <- sum(v_conc_aer[8L:13L]) * v_pars["tss_pcod"][[1L]] # MLTSS in aerobic reactor
    kla_aer <- m_sim_vode[nrow(x = m_sim_vode), "kla_aer"][[1L]]
    tss_eff <- m_sim_vode[nrow(x = m_sim_vode), "tss_eff"][[1L]] # TSS effluent
    sludge_prod <- m_sim_vode[nrow(x = m_sim_vode), "sludge_prod"][[1L]] # sludge production, kgSS/d
    remove(v_pars)
    
    # Evaluate if steady state has been reached
    precis <- m_sim_vode[nrow(x = m_sim_vode), "precis"][[1L]]
    remove(m_sim_vode)
    
    # End function
    list(v_conc_anx = v_conc_anx,
         v_conc_aer = v_conc_aer,
         precis = precis,
         mltss = mltss,
         kla_aer = kla_aer,
         tss_eff = tss_eff,
         sludge_prod = sludge_prod)
    
}

# End of code --------------------------------------------------------------------------------------
