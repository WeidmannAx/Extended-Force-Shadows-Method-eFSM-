function [x_Bias,Sigma_t_1, w_t]=Bias_Estimation(Insole, i, t0, Estim, Area, left, param, Sigma_t_1, w_t, polynom_degree)
Affiliation.Heel = [1, 1, 1, 1, 0.3465, 0.1636]; Affiliation.Arch = [0.6535, 0.8364, 0.7716, 0.6438];
Affiliation.Meta = [0.2284, 0.3562, 0.6508, 0.7621, 0.726, 0.5974, 0.4823]; Affiliation.Toes = [0.3492, 0.2379, 0.274, 0.4026, 0.5177, 1, 1, 1];

Affiliation.Heel = [Affiliation.Heel, zeros(1,10)];
Affiliation.Arch = [0, 0, 0, 0, Affiliation.Arch, zeros(1,8)];
Affiliation.Meta = [zeros(1,6), Affiliation.Meta, 0, 0, 0];
Affiliation.Toes = [zeros(1,8), Affiliation.Toes];

PadNr = 16;
%% initialization / building eFSM force vector with corresponding Moticon sensor pad affiliation
if left
    y_m_t = (Insole.pressuresLeft(i+t0, :) .* (Area / 9.81));
    L_Omega_k = [Estim(1), Estim(2), sum(Estim(3:4)), sum(Estim(5:6))];
    y_FSM_tilde = zeros(1,16);
    if sum(y_m_t .* Affiliation.Heel) ~= 0
        y_FSM_tilde = y_FSM_tilde + (y_m_t .* Affiliation.Heel)*L_Omega_k(1)/sum(y_m_t .* Affiliation.Heel);
    end
    if sum(y_m_t .* Affiliation.Arch) ~= 0
        y_FSM_tilde = y_FSM_tilde + (y_m_t .* Affiliation.Arch)*L_Omega_k(2)/sum(y_m_t .* Affiliation.Arch);
    end
    if sum(y_m_t .* Affiliation.Meta) ~= 0
        y_FSM_tilde = y_FSM_tilde + (y_m_t .* Affiliation.Meta)*L_Omega_k(3)/sum(y_m_t .* Affiliation.Meta);
    end
    if sum(y_m_t .* Affiliation.Toes) ~= 0
        y_FSM_tilde = y_FSM_tilde + (y_m_t .* Affiliation.Toes)*L_Omega_k(4)/sum(y_m_t .* Affiliation.Toes);
    end
    y_FSM = y_FSM_tilde;
    y_k = (Insole.pressuresLeft(i+t0, :) .* (Area / 9.81))' - y_FSM';
    X = (Insole.pressuresLeft(i+t0, :) .* (Area / 9.81))';
    Y = y_k;
else
    y_m_t = (Insole.pressuresRight(i+t0, :) .* (Area / 9.81));
    L_Omega_k = [Estim(1), Estim(2), sum(Estim(3:4)), sum(Estim(5:6))];
    
    y_FSM_tilde = zeros(1,16);
    if sum(y_m_t .* Affiliation.Heel) ~= 0
        y_FSM_tilde = y_FSM_tilde + (y_m_t .* Affiliation.Heel)*L_Omega_k(1)/sum(y_m_t .* Affiliation.Heel);
    end
    if sum(y_m_t .* Affiliation.Arch) ~= 0
        y_FSM_tilde = y_FSM_tilde + (y_m_t .* Affiliation.Arch)*L_Omega_k(2)/sum(y_m_t .* Affiliation.Arch);
    end
    if sum(y_m_t .* Affiliation.Meta) ~= 0
        y_FSM_tilde = y_FSM_tilde + (y_m_t .* Affiliation.Meta)*L_Omega_k(3)/sum(y_m_t .* Affiliation.Meta);
    end
    if sum(y_m_t .* Affiliation.Toes) ~= 0
        y_FSM_tilde = y_FSM_tilde + (y_m_t .* Affiliation.Toes)*L_Omega_k(4)/sum(y_m_t .* Affiliation.Toes);
    end
    y_FSM = y_FSM_tilde;
    y_k = (Insole.pressuresRight(i+t0, :) .* (Area / 9.81))' - y_FSM';
    X = (Insole.pressuresRight(i+t0, :) .* (Area / 9.81))';
    Y = y_k;
end
%% Regression Algorithm
Phi = ones(PadNr,polynom_degree+1);
for k = 1:polynom_degree
    Phi(:,k+1) = Phi(:,k).*X;
end
for k = 1:16
    Sigma_t_p = Sigma_t_1((k-1)*(polynom_degree+1)+1:k*(polynom_degree+1),:);
    Sigma_t = pinv(Phi(k,:)' * param.filter_inv_meas_variance(k,k) * Phi(k,:) + pinv(Sigma_t_p));
    w_t(:,k) = Sigma_t * (Phi(k,:)' * param.filter_inv_meas_variance(k,k) * Y(k) + pinv(Sigma_t_p) * w_t(:,k));
    Sigma_t_1((k-1)*(polynom_degree+1)+1:k*(polynom_degree+1),:) = Sigma_t;
    x_Bias(k,1) = Phi(k,:)*w_t(:,k);
end

end