function dX = get_ionic_fun_dX(ionic_fun, t, x, p, S)
X = cell(1,4);
[X{:}] = ionic_fun(t, x, p, S);
dX = X{end};