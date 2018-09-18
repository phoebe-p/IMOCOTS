function[layer_indices] = needed_layers(SC)
    top_ind = [];
    bot_ind = [];
    calc_layer = [];
    SC_layer_index = fliplr(2:(length(SC.layer) - 1)); % exclude incidence & transmission medium
    for i = SC_layer_index % exclude incidence & transmission medium
        if ~isempty(SC.layer(i).type)
            if SC.layer(i).type == 2
                calc_layer(end+1) = i;
                top_ind(end+1) = SC.layer(i).grating_layer_indices(2);
                bot_ind(end+1) = SC.layer(i).grating_layer_indices(1);
                
                % if the layers above and below another layer are
                % full-field, then can easily calculate power flow in the
                % intermediate layer as well.
            elseif SC.layer(i).type ~= 2 & SC.layer(i-1).type == 2 & SC.layer(i+1).type == 2
                calc_layer(end+1) = i;
                top_ind(end+1) = SC.layer(i+1).grating_layer_indices(1);
                bot_ind(end+1) = top_ind(end) - 1 - SC.layer(i).n_strata;
            end
        end
    end
    layer_indices.calc_layer = calc_layer;
    layer_indices.t = top_ind;
    layer_indices.b = bot_ind;
end
