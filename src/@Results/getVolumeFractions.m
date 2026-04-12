function [phi,z,t] = getVolumeFractions(obj)
    z = obj.SandFilter.GridPoints.Centers;
    t = obj.Frames.Time;

    subphase_type = {["Particles" ,"Matrix"]; 
                     ["Particles", "Enclosed"]; 
                     ["Liquids", "Enclosed"]; 
                     ["Particles", "Flowing"];
                     ["Liquids", "Flowing"]};
    subphase_index = {'M'; 'Pe'; 'Le'; 'Pf'; 'Lf'};
    sub_phases = [subphase_type, subphase_index];

    phi_sum = cell(length(sub_phases(:,2)),1);
    [phi_sum{1:end}] = deal(zeros(length(z),length(t)));
    phi_sum = table(phi_sum{:}, 'VariableNames' ,sub_phases(:,2));

    for s = 1:size(sub_phases, 1)
        component_type = sub_phases{s, 1}(1);
        volume_type = sub_phases{s, 1}(2);
        phase_name = sub_phases{s, 2};

        rho = [obj.Model.(component_type).Density];
        names = [obj.Model.(component_type).Name];
        for i = 1:length(rho)
            phi_sum.(phase_name) = phi_sum.(phase_name) + ...
                obj.Frames.Concentrations{names(i), volume_type}{:} / rho(i);
        end
    end
    
    phiWe = obj.Frames.Concentrations{"Water", "Enclosed"}{:} / obj.Model.WaterDensity;
    phiM = phi_sum.("M");
    phie = phi_sum.("Pe") + phi_sum.("Le") + phiWe;
    phib = phiM + phie;

    phi = struct('Biofilm',phib, 'Enclosed', phie, 'Matrix',phiM); 
end