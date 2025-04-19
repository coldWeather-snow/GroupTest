function loss = calc_loss_A(my_FIM)
% A-optimality loss: trace(inv(FIM))
    loss = trace(inv(my_FIM));
end