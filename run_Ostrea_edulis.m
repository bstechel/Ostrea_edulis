for i = 1:1
    
close all; 
global pets

pets = {'Ostrea_edulis'};
check_my_pet(pets); 

estim_options('default'); 
estim_options('max_step_number',5000); 
estim_options('max_fun_evals',5000);  

estim_options('pars_init_method',2);
estim_options('results_output', 2);
estim_options('method', 'no');

estim_pars; 


end
