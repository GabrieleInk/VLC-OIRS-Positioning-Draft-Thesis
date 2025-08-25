function p=WILS(Coord_Ancore,distanze,initial_guess,zk,w) % Iterative Weighted Least Squares
W=diag(w);
M=size(Coord_Ancore,1);

sp_ = Coord_Ancore'; 
pr_ = distanze'; 

gu_(1)=initial_guess(1); gu_(2)=initial_guess(2); % initial guesses for node position 
% Calculating rao: the pseudo-range
for j_=1:M
    rao_(j_)=((gu_(1)-sp_(1,j_))^2+(gu_(2)-sp_(2,j_))^2 +(zk-sp_(3,j_))^2)^.5;
end

erro_ = 1;
exit_ = 0;

while erro_ > .000001 && exit_<200
%while erro_ > .01 && exit_<10 %When the threshold is unde
exit_=exit_+1;
    alpha_ = [];
    for j_ = 1:M
        for k_ = 1:2
            alpha_(j_,k_) = (gu_(k_)-sp_(k_,j_))/(rao_(j_)); % find first 2 columns of alpha matrix
        end
    end
    %alpha_ = [alpha_ 1*ones(M,1)];
    
    drao_ = pr_ - rao_'; %** find delta rao includes clock bias
    %dl_ = pinv(alpha_)*drao_; % Equation (2.16)
    dl_ = (alpha_' * W * alpha_) \ (alpha_' * W * drao_);
    gu_ = gu_ + dl_(1:2)'; %**find new position
    erro_ = dl_(1)^2+dl_(2)^2; % find error
    for j_ = 1:M
        rao_(j_) = ((gu_(1)-sp_(1,j_))^2+(gu_(2)-sp_(2,j_))^2 +(zk-sp_(3,j_))^2)^.5; % find new rao without clock bias
    end
end

if(gu_(1) < 0 || gu_(1) > 5 || gu_(2) < 0 || gu_(2) > 5 )
p = [initial_guess(1:2) 0]; % The gu is localized outside the room, return initial guess
else
p=[gu_ 0]; % The gu is localized inside the room, return ILS result
end

