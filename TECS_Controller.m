function [d_thr, d_elv, error_T_dot, error_L_dot] = TECS_Controller(c)
    
    m = c(1);
    v_sp = c(2);
    h_sp = c(3);
    gamma_sp = c(4);
    v = c(5);
    h = c(6);
    gamma = c(7);
    v_dot = c(8);
    error_T = c(9);
    error_L = c(10);


    % TECS PID Constants
    KTP = -10;
    KTI = -0.01;
    KTD = 0;
    
    KLP = 0.0001;
    KLI = 5e-8;
    KLD = 0;
    
    KV = 0.5;
    
    
    % TECS Setpoints
    sp_energy = totalEnergy(m, v_sp, h_sp);
    sp_energyRate = totalEnergyRate(KV*v_sp, gamma_sp);
    sp_balance = energyBal(m, v_sp, h_sp);
    sp_balanceRate = energyBalRate(KV*v_sp, gamma_sp); 

    
    % TECS Current state
    energy = totalEnergy(m, v, h);
    energyRate = totalEnergyRate(KV*v, gamma);
    balance = energyBal(m, v, h);
    balanceRate = energyBalRate(KV*v, gamma);

    
    % TECS Errors
    err_energy = (sp_energy - energy)/1e6;
    err_energyRate = sp_energyRate - energyRate;
    err_balance = (sp_balance - balance)/1e6;
    err_balanceRate = sp_balanceRate - balanceRate;
    
    error_T = error_T + err_energy;
    error_L = error_L + err_balance;
    
        
    % TECS PID
    d_thr = err_energy *KTP + KTI*error_T + KTD*err_energyRate;
    d_elv = err_balance*KLP + KLI*error_L + KLD*err_balanceRate;
    
    
    
    
    error_T_dot = err_energy;
    error_L_dot = err_balance;


end



function energy = totalEnergy(m, v, h) 
    energy = 0.5 * m*v^2  + m*h*9.81;
end

function energyRate = totalEnergyRate(v_dot, gamma) 
    energyRate = v_dot/9.81  + gamma;
end

function balance = energyBal(m, v, h)
    balance = 0.5 * m*v^2  - m*h*9.81;
    % Negative means more potential than kinetic - Likely case
end

function balanceRate = energyBalRate(v_dot, gamma) 
    balanceRate = v_dot/9.81  - gamma;
end