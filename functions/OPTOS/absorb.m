function [absorb_bulk, absorb_front, absorb_back, R, T] = absorb(iter_or_geom,B,C,D_si,v0,points,power_threshold )	

if iter_or_geom == 'iter'	%**********************************************************************************************
							%Calculation of absorptance via iteration
							%**********************************************************************************************
	disp('Iterative calculation of sheet absorptance...')
	v_rear = D_si*v0;
	iter_index = 0;
	abs_iter_bulk = sum(v0(1:end-1))-sum(v_rear(1:end-1));
    %abs_iter_front = v0(end-1);
    %abs_iter_back = v0(end); % should always be zero
    R = 1-sum(v0);
    T = 0;
    
	while(sum(v_rear(1:end-2))>= power_threshold)
        %size(C)
        %size(v_rear)
		vstrich_rear = C*v_rear; % after interaction with rear
        %a = sum(vstrich_rear(1:end-1))
        
        T = T + sum(v_rear) - sum(vstrich_rear);
        
		v_front = D_si*vstrich_rear; % at front, after passing through bulk,
        % before interacting with front
         %b = sum(v_front(1:end-1))

		vstrich_front = B*v_front; % at front, after interacting with front
         %c = sum(vstrich_front(1:end-1))
         
		v_rear = D_si*vstrich_front;
        %d = sum(v_rear)% at back, after passing through bulk, 
        % before interacting with back
		abs_iter_bulk = abs_iter_bulk + sum(vstrich_rear(1:end-1)) - sum(v_front(1:end-1)) + sum(vstrich_front(1:end-1)) - sum(v_rear(1:end-1));
		%abs_iter_front = abs_iter_front + v_rear(end-1);
        %abs_iter_back = abs_iter_back + v_rear(end);
        
        R = R + sum(v_front) - sum(vstrich_front);
        
        iter_index = iter_index+1;
        %        v_rear(end-1)
        %v_rear(end)
	end
	absorb_bulk = abs_iter_bulk;
    %absorb_front = abs_iter_front;
    %absorb_back = abs_iter_back;
    absorb_front = v_rear(end-1);
    absorb_back = v_rear(end);
	%iter_index;
elseif iter_or_geom == 'geom' 	%**********************************************************************************************
								%Calculation of absorptance via geometrical series
								%**********************************************************************************************
	disp('Calculation of sheet absorptance based on geometric series...')	
	% Gesmatabsorption mittels geometrischer Reihe berechnet. Es ergibt sich die Gesamtabsorption nach unendlich vielen Strahldurchgängen. Keine Ortsinformation!
	to_sum = inv(sparse(eye(size(points,2)))-B*D_si*C*D_si)*v0-inv(sparse(eye(size(points,2)))...
	-D_si*B*D_si*C)*D_si*v0 + inv(sparse(eye(size(points,2)))-C*D_si*B*D_si)*C*D_si*v0...
	-inv(sparse(eye(size(points,2)))-D_si*C*D_si*B)*D_si*C*D_si*v0;
save('to_sum', 'to_sum');
    
    absorb_bulk = sum(to_sum(1:end-2));
    
    % DOES NOT WORK:
    absorb_front = to_sum(end-1);
    absorb_back = to_sum(end);
end
disp('Calculation of sheet absorptance done.')	