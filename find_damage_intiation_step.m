function find_damage_intiation_step

ubar = 0.0107;
damage_error = 1;
damage_tolerance = 0.00001;
ubar_vs_damage = [];
while  damage_error > damage_tolerance 
    damage_data = main_code(ubar);
    min_damage_value = 10;
    for guass_point = 1 : 25600 
        damage = damage_data(guass_point,1);
        if(damage > 0)
            min_damage_value = min(min_damage_value,damage);
        end 
    end 
    disp(['displacement ' num2str(ubar)]);
    disp(['minimum damage value ' num2str(min_damage_value)]);
    damage_error = (min_damage_value - 0.001)/0.001;
    disp(['damage error ' num2str(damage_error)]);
    ubar_vs_damage = [ubar_vs_damage; ubar min_damage_value];
    if(damage_error > damage_tolerance)
        ubar = ubar - 0.00002;
    end 
end 