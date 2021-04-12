function [d_thr, d_elv, error_T_dot, error_L_dot] = TECS_Controller(c)
    
    g = 9.81;
    m = c(1);
    v_sp = c(2);
    h_sp = c(3);
    gamma_sp = c(4);
    v = c(5);
    h = c(6);
    gamma = c(7);
    error_T = c(8);
    error_L = c(9);

    % TECS PID Constants
    KTP = -10;
    KTI = -0.001;
    KTD = 0;
    
    KLP = 0.00015;
    KLI = 5e-8;
    KLD = 0.94;
    
    KV = 0.22;
    
    
    % TECS Errors
    err_energy = (0.5 * m*(v_sp^2 - v^2)  + m*g*(h_sp-h))/1e6;
    err_energyRate = KV*(v_sp - v)/g  + (gamma_sp-gamma);
    err_balance = (0.5 * m*(v_sp^2 - v^2)  - m*g*(h_sp-h))/1e6;
    err_balanceRate = KV*(v_sp - v)/g  - (gamma_sp-gamma);
    
    
    error_T = error_T + err_energy;
    error_L = error_L + err_balance;
    
        
    % TECS PID
    d_thr = err_energy *KTP + KTI*error_T + KTD*err_energyRate;
    d_elv = err_balance*KLP + KLI*error_L + KLD*err_balanceRate;
    
   
    error_T_dot = err_energy;
    error_L_dot = err_balance;
end
