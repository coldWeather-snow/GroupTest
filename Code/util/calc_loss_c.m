function loss = calc_loss_c(my_FIM, c)
% c-optimality loss: c' * inv(FIM) * c
    loss = c' * inv(my_FIM) * c;
end