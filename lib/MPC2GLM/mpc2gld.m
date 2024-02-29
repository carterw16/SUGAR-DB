function str = mpc2gld(mpc,omega,varargin)
%%% convert a matpower case structure to a gridlabd string.
%%% Since gridlabd does not support PV buses, all PV buses are converted to
%%% PQ. Phase shifts are not supported currently but off-nominal tap ratios
%%% are.
%%% 
%%% KNOWN ISSUES:
%%%     Base KV MUST be specified in the mpc.bus matrics, otherwise all
%%%     voltages will be set to 0 in the glm model, which wll result in
%%%     errors.
%%%
%%% INPUTS:
%%%        mpc: matpower case
%%%        omega: power radian frequency (2*pi*f)
%%% OPTIONAL NAME VALUE PAIRS:
%%%        exclude_buses: default empty.
%%%                       vector of bus numbers to exclude from the model.
%%%                       example: mpc2gld(mpc,omega,'exclude_buses',[5,7])
%%%                       nodes 5 and 7 will not be created. However, edges
%%%                       connecting to nodes 5 and 7 will still be
%%%                       written!!! The assumption is that these nodes
%%%                       will be added elsewhere as part of a feeder
%%%                       model.
%%%         no_shunt: default true.
%%%                   If false, no susceptances will be added to the glm.
%%% OUTPUTS:
%%%        str: formated string in glm format (note: preamble such as
%%%        module loading/clock etc. is provided.)
%%%
%%% written July, 2017 by Eran Schweitzer (eranschweitzer@gmail.com) at PNNL

idx = find(strcmp(varargin,'exclude_buses'));
if ~isempty(idx)
    exclude_buses = varargin{idx + 1};
else
    exclude_buses = [];
end

idx = find(strcmp(varargin,'with_shunt'));
if ~isempty(idx)
    shunt = varargin{idx+1};
else
    shunt = true;
end


define_constants;
str = '';
nmap = sparse(mpc.bus(:,BUS_I),1,1:size(mpc.bus,1));

for n = 1:size(mpc.bus,1)
    if ~any(mpc.bus(n,BUS_I) == exclude_buses)
        if (abs(mpc.bus(n,PD)) > 0) || (abs(mpc.bus(n,QD)) > 0) || (abs(mpc.bus(n,GS)) > 0) || (abs(mpc.bus(n,BS)) > 0)
            str = strcat(str,load_str(mpc.bus(n,:)));
        else 
            str = strcat(str,node_str(mpc.bus(n,:)));
        end
    end
    
%     if (abs(mpc.bus(n,BS)) > 0) %&& (abs(mpc.bus(n,GS)) == 0)
%         str = strcat(str,cap_str(mpc.bus(n,:)));
%     end
end

branch_mask = false(size(mpc.branch,1),1);
sorted_branches = [min(mpc.branch(:,1:2),[],2), max(mpc.branch(:,1:2),[],2) ];
for b = 1:size(mpc.branch,1)
    if branch_mask(b)
        continue
    end
    if mpc.branch(b,BR_STATUS)
        par_mask = parallel_mask(sorted_branches,b);
        
        par_mask(b:end) = false;
        branch_mask(b) = true;
        dummy_num = sum(par_mask);
        if dummy_num > 0
            %parallel branch, add dummy node    
%             str = strcat(str,dummy_node(mpc.bus(mpc.branch(b,F_BUS),:),mpc.branch(b,T_BUS),dummy_num));
            str = strcat(str,dummy_node(mpc.bus(nmap(mpc.branch(b,F_BUS)),:),mpc.branch(b,T_BUS),dummy_num));
        end
        
%         branch_mask(par_mask) = true;
%         par_mask = par_mask & logical(mpc.branch(b,BR_STATUS)); %remove any out of service branches
        
        % a transformer has any of the following properties:
        % 1) different voltages at the ends of the branch
        % 2) a tap setting not equal to 0 (even tap of 1 will be used as a
        % transfomer)
        xfmr_check = mpc.bus(nmap(mpc.branch(b,F_BUS)),BASE_KV) ~= mpc.bus(nmap(mpc.branch(b,T_BUS)),BASE_KV);
        xfmr_check = xfmr_check || (mpc.branch(b,TAP) ~= 0);
        if ~xfmr_check
            Zbase = mpc.bus(nmap(mpc.branch(b,F_BUS)),BASE_KV)^2/mpc.baseMVA;
%             z = mpc.branch(par_mask,BR_R) + 1i*mpc.branch(par_mask,BR_X);
%             z = sum(z.^(-1))^(-1)*Zbase;
            z = Zbase*(mpc.branch(b,BR_R) + 1i*mpc.branch(b,BR_X));
%             c = 1e9*sum(mpc.branch(par_mask,BR_B))/(omega*Zbase);
            if shunt
                c = 1e9*mpc.branch(b,BR_B)/(omega*Zbase);
            else
                c = 0;
            end
            str = strcat(str,line_str(mpc.branch(b,:),z,c,dummy_num));
        else
            % treat as a transformer
%             z = mpc.branch(par_mask,BR_R) + 1i*mpc.branch(par_mask,BR_X);
%             z = sum(z.^(-1))^(-1);
            z = mpc.branch(b,BR_R) + 1i*mpc.branch(b,BR_X);
%             bshunt = sum(mpc.branch(par_mask,BR_B));
            if shunt
                bshunt = mpc.branch(b,BR_B);
            else
                bshunt = 0;
            end
            if bshunt ~= 0
                zshunt = 1/(1i*bshunt);
                str = strcat(str,xfmr_shunt(mpc.branch(b,F_BUS),mpc.branch(b,T_BUS),...
                mpc.bus(nmap(mpc.branch(b,F_BUS)),BASE_KV),mpc.bus(nmap(mpc.branch(b,T_BUS)),BASE_KV),...
                mpc.baseMVA,zshunt,dummy_num));
            else
                zshunt = 0;
            end
            str = strcat(str,xfmr_str(mpc.branch(b,:), z, zshunt, mpc.baseMVA,...
                mpc.bus(nmap(mpc.branch(b,F_BUS)),BASE_KV),mpc.bus(nmap(mpc.branch(b,T_BUS)),BASE_KV),dummy_num));
        end
    end
end

for g = unique(mpc.gen(:,GEN_BUS)).'
%     fprintf('adding generator(s) on bus %d\n',g)
    rows = find(mpc.gen(:,GEN_BUS) == g);
    for gnum = 1:length(rows)
        if mpc.gen(rows(gnum),GEN_STATUS) > 0
           str = strcat(str,gen_str(mpc.gen(rows(gnum),:),gnum,mpc.bus(nmap(g),BASE_KV))); 
        end
    end
end

end

function str = dummy_node(bus_row,tobus,dummy_num)
    str = node_str(bus_row,sprintf('_to%d_dummy%d',tobus,dummy_num));
    str = strcat(str,sprintf(['\nobject switch {\n',...
        'name dummy_switch%d_bus%d_to%d;\n',...
        'phases ABC;\n',...
        'from bus%d;\n',...
        'to bus%d_to%d_dummy%d;\n',...
        'status CLOSED;\n',...
        '}\n'],...
        dummy_num,bus_row(1),tobus,bus_row(1),bus_row(1),tobus,dummy_num));
end
function mask = parallel_mask(branches,bid)
    mask = (branches(:,1) == branches(bid,1)) & (branches(:,2) == branches(bid,2)); 
end

function str = gen_str(gen_row,gnum,baseKV)
    %generators will simply be treated like negative loads 
    
    define_constants;
    S = -( gen_row(PG) + 1i*gen_row(QG) ); %negative since treated as load in gld.
    str = sprintf(['\nobject load {\n',...
        '\tname gen%d_bus%d;\n',...
        '\tparent bus%d;\n',...
        '\tphases ABC;\n',...
        '\tnominal_voltage %0.3f kV;\n',...
        '\tconstant_power_A %s MVA;\n',...
        '\tconstant_power_B %s MVA;\n',...
        '\tconstant_power_C %s MVA;\n',...
        '}\n'],gnum,gen_row(GEN_BUS),gen_row(GEN_BUS),baseKV/sqrt(3),...
        complex_string(S/3),complex_string(S/3),complex_string(S/3));
end

function str = xfmr_shunt(from,to,baseKV1,baseKV2,baseMVA,zpu,dummy_num)
    
    if dummy_num > 0
        from = sprintf('%d_to%d_dummy%d',from,to,dummy_num);
    else
        from = sprintf('%d',from);
    end
    
    Zbase1 = baseKV1^2/baseMVA;
    Zbase2 = baseKV2^2/baseMVA;
    str = sprintf(['\nobject load {\n',...
        '\tname xfmr_%s_%d_shunt_from;\n',...
        '\tparent bus%s;\n',...
        '\tphases ABC;\n',...
        '\tnominal_voltage %0.3f kV;\n',...
        '\tconstant_impedance_A %s Ohm;\n',...
        '\tconstant_impedance_B %s Ohm;\n',...
        '\tconstant_impedance_C %s Ohm;\n',...
        '}\n'],from,to,from,baseKV1/sqrt(3),...
        complex_string(2*zpu*Zbase1),complex_string(2*zpu*Zbase1),complex_string(2*zpu*Zbase1));
    
    str = strcat(str,sprintf(['\nobject load {\n',...
        '\tname xfmr_%s_%d_shunt_to;\n',...
        '\tparent bus%d;\n',...
        '\tphases ABC;\n',...
        '\tnominal_voltage %0.3f kV;\n',...
        '\tconstant_impedance_A %s Ohm;\n',...
        '\tconstant_impedance_B %s Ohm;\n',...
        '\tconstant_impedance_C %s Ohm;\n',...
        '}\n'],from,to,to,baseKV2/sqrt(3),...
        complex_string(2*zpu*Zbase2),complex_string(2*zpu*Zbase2),complex_string(2*zpu*Zbase2)));
end

function str = xfmr_str(branch_row,z, zshunt, baseMVA,baseKV1,baseKV2,dummy_num)
    

    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch; %#ok<ASGLU>
    
    if dummy_num > 0
        dummy_name = sprintf('_dummy%d',dummy_num);
        dummy_from = sprintf('_to%d_dummy%d',branch_row(T_BUS),dummy_num);
    else
        dummy_name = '';
        dummy_from = '';
    end

    
    tap = branch_row(TAP);
    if tap == 0
        tap = 1;
    end
    
%     if baseKV1 > baseKV2
        from = branch_row(F_BUS);
        to   = branch_row(T_BUS);
%     else
%         to   = branch_row(F_BUS);
%         from = branch_row(T_BUS);   
%     end
    str = sprintf(['\n',...
    'object transformer_configuration {\n',...
    '\tname trans_branch_config_%d_%d%s;\n',...
    '\tconnect_type WYE_WYE;\n',...
    '\tpower_rating %0.1f MVA;\n',...
    '\tprimary_voltage %0.2f kV;\n',...
    '\tsecondary_voltage %0.2f kV;\n',...
    '\tresistance %0.6f;\n',...
    '\treactance %0.6f;\n'],...
        branch_row(F_BUS),branch_row(T_BUS),dummy_name,baseMVA,...
        baseKV1*tap,baseKV2,max(real(z),1e-6),imag(z));
    
%     if zshunt ~= 0
%        str = strcat(str,sprintf('\n\tshunt_impedance %s;\n',complex_string(zshunt))); 
%     end
    
    str = strcat(str,sprintf('\n}\n'));
    
    
    str = strcat(str,sprintf(['\nobject transformer {\n',...
    '\tname trans_branch_%d_%d%s;\n',...
    '\tphases ABC;\n',...
    '\tfrom bus%d%s;\n',...
    '\tto bus%d;\n',...
    '\tconfiguration trans_branch_config_%d_%d%s;\n',...
    '}\n'],...
    branch_row(F_BUS),branch_row(T_BUS),dummy_name,...
    from,dummy_from,to,...
    branch_row(F_BUS),branch_row(T_BUS),dummy_name));

end

function str = line_str(branch_row,z,c,dummy_num)

    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch; %#ok<ASGLU>

    if dummy_num > 0
        dummy_name = sprintf('_dummy%d',dummy_num);
        dummy_from = sprintf('_to%d_dummy%d',branch_row(T_BUS),dummy_num);
    else
        dummy_name = '';
        dummy_from = '';
    end

%     z = (branch_row(BR_R) + 1i*branch_row(BR_X))*baseKV^2/baseMVA;
%     c = 1e9*branch_row(BR_B)*baseMVA/(baseKV^2*omega);
    str = sprintf(['\nobject line_configuration {\n',...
            '\tname line_config_%d_%d%s;\n',...
            '\tz11 %s Ohm/mile;\n',...
            '\tz22 %s Ohm/mile;\n',...
            '\tz33 %s Ohm/mile;\n',...
            '\tc11 %0.4f nF/mile;\n',...
            '\tc22 %0.4f nF/mile;\n',...
            '\tc33 %0.4f nF/mile;\n',...
            '}\n'],branch_row(F_BUS),branch_row(T_BUS),dummy_name,...
            complex_string(z),complex_string(z),complex_string(z),...
            c,c,c);
        
     str = strcat(str,sprintf(['\nobject overhead_line {\n',...
             '\tname line_%d_%d%s;\n',...
             '\tphases ABC;\n',...
             '\tfrom bus%d%s;\n',...
             '\tto bus%d;\n',...
             '\tlength 1 mile;\n',...
             '\tconfiguration line_config_%d_%d%s;\n}\n'],...
             branch_row(F_BUS),branch_row(T_BUS),dummy_name,...
             branch_row(F_BUS),dummy_from,branch_row(T_BUS),...
             branch_row(F_BUS),branch_row(T_BUS),dummy_name));
end

function str = load_str(bus_row)
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus; %#ok<ASGLU>
    str = sprintf(['\nobject load {\n',...
        '\tname bus%d;\n',...
        '\tphases ABC;\n',...
        '\tnominal_voltage %0.3f kV;\n',...
        '\tbustype %s;\n'],bus_row(BUS_I),bus_row(BASE_KV)/sqrt(3),bustype(bus_row(BUS_TYPE)));
    
    v = bus_row(VM)*exp(1i*bus_row(VA)*pi/180)*bus_row(BASE_KV)/sqrt(3);
    a = exp(1i*120*pi/180);
    str = strcat(str,sprintf(['\n\tvoltage_A %s kV;\n',...
        '\tvoltage_B %s kV;\n',...
        '\tvoltage_C %s kV;\n'],...
        complex_string(v,'format','polar'),...
        complex_string(v*a^2,'format','polar'),...
        complex_string(v*a,'format','polar')));
    
    if (abs(bus_row(PD)) > 0) || (abs(bus_row(QD)) > 0) 
        S = bus_row(PD) + 1i*bus_row(QD);
        str = strcat(str,sprintf(['\n\tconstant_power_A %s MVA;\n',...
                '\tconstant_power_B %s MVA;\n',...
                '\tconstant_power_C %s MVA;\n'],...
                complex_string(S/3),complex_string(S/3),complex_string(S/3)));
    end
    
    if (abs(bus_row(GS)) > 0) || (abs(bus_row(BS)) > 0)
        S = bus_row(GS) - 1i*bus_row(BS); %power in MVA
        Z = bus_row(BASE_KV)^2/(conj(S)); %impedance in Ohm.
%         if real(Z) == 0
%             Z = Z + 3e6;
%         end
        str = strcat(str,sprintf(['\n\tconstant_impedance_A %s Ohm;\n',...
                '\tconstant_impedance_B %s Ohm;\n',...
                '\tconstant_impedance_C %s Ohm;\n'],...
                complex_string(Z),complex_string(Z),complex_string(Z)));    
    end
    
    str = strcat(str,sprintf('\n}\n'));
end

function str = cap_str(bus_row)
    
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus; %#ok<ASGLU>
    
    str = sprintf(['\nobject capacitor {\n',...
        '\tname cap%d;\n',...
        '\tparent bus%d;\n',...
        '\tphases ABC;\n',...
        '\tphases_connected ABC;\n',...
        '\tnominal_voltage %0.3f kV;\n',...
        '\tcapacitor_A %0.3f MVAr;\n',...
        '\tcapacitor_B %0.3f MVAr;\n',...
        '\tcapacitor_C %0.3f MVAr;\n',...
        '}\n'],...
        bus_row(BUS_I),bus_row(BUS_I),bus_row(BASE_KV)/sqrt(3),...
        bus_row(BS)/3,bus_row(BS)/3,bus_row(BS)/3);
        
end

function str = node_str(bus_row,dummy_name)
    if nargin == 1
        dummy_name = '';
    end
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus; %#ok<ASGLU>
    str = sprintf(['\nobject node {\n',...
        '\tname bus%d%s;\n',...
        '\tphases ABC;\n',...
        '\tnominal_voltage %0.3f kV;\n',...
        '\tbustype %s;\n'],bus_row(BUS_I),dummy_name,bus_row(BASE_KV)/sqrt(3),bustype(bus_row(BUS_TYPE)));
    
    v = bus_row(VM)*exp(1i*bus_row(VA)*pi/180)*bus_row(BASE_KV)/sqrt(3);
    a = exp(1i*120*pi/180);
    str = strcat(str,sprintf(['\n\tvoltage_A %s kV;\n',...
        '\tvoltage_B %s kV;\n',...
        '\tvoltage_C %s kV;\n'],...
        complex_string(v,'format','polar'),...
        complex_string(v*a^2,'format','polar'),...
        complex_string(v*a,'format','polar')));
    
    str = strcat(str,sprintf('\n}\n'));
end

function str = bustype(n)
    switch n
        case 1
            str = 'PQ';
        case 2
            str = 'PQ'; % PV buses are not implemented in gridlabd!
        case 3
            str = 'SWING';
    end
end

function str = complex_string(v,varargin)
    idx = find(strcmp(varargin,'format'),1);
    if ~isempty(idx)
        fmt = varargin{idx + 1};
    else
        fmt = 'rect';
    end
    if strcmp(fmt,'rect')
        if imag(v) < 0
            str = sprintf('%0.4f%0.4fj',real(v),imag(v));
        else
            str = sprintf('%0.4f+%0.4fj',real(v),imag(v));
        end
    elseif strcmp(fmt,'polar')
        if angle(v) < 0
            str = sprintf('%0.4f%0.4fd',abs(v),angle(v)*180/pi);
        else
            str = sprintf('%0.4f+%0.4fd',abs(v),angle(v)*180/pi);
        end
    end
end