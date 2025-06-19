function loss = calc_loss_E(my_FIM)
% c-optimality loss: c' * inv(FIM) * c
    loss = -lambda_min(my_FIM);
end