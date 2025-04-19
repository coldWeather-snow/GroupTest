function loss = calc_loss_D(my_FIM)
% D-optimality loss: 1 / det(FIM)
    loss = 1 / det(my_FIM);
end