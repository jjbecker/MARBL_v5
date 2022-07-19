function [my_cell] = global_tracer_units(lciso_on)
%tracer_names MARBL "units" of all tracers

% if MARBL is shutdown then mex_marbl_driver('tracer_sname', i) crashes...
%    just do a lot of typing...

my_cell{ 1} = '(mol)';
my_cell{ 2} = '(mol)';
my_cell{ 3} = '(mol)';
my_cell{ 4} = '(mol)';
my_cell{ 5} = '(mol)';
my_cell{ 6} = '(mol)';
my_cell{ 7} = '(mol)';
my_cell{ 8} = '(mol)';
my_cell{ 9} = '(mol)';
my_cell{10} = '(eq)';
my_cell{11} = '(eq)';
my_cell{12} = '(mol)';
my_cell{13} = '(mol)';
my_cell{14} = '(mol)';
my_cell{15} = '(mol)';
my_cell{16} = '(mol)';
my_cell{17} = '(mol)';
my_cell{18} = '(mol)';
my_cell{19} = '(g)';
my_cell{20} = '(mol)';
my_cell{21} = '(mol)';
my_cell{22} = '(mol)';
my_cell{23} = '(mol)';
my_cell{24} = '(g)';
my_cell{25} = '(mol)';
my_cell{26} = '(mol)';
my_cell{27} = '(mol)';
my_cell{28} = '(mol)';
my_cell{29} = '(g)';
my_cell{30} = '(mol)';
my_cell{31} = '(mol)';
my_cell{32} = '(mol)';
if lciso_on == 1
    my_cell{33} = '(mol)';
    my_cell{34} = '(mol)';
    my_cell{35} = '(mol)';
    my_cell{36} = '(mol)';
    my_cell{37} = '(mol)';
    my_cell{38} = '(mol)';
    my_cell{39} = '(mol)';
    my_cell{40} = '(mol)';
    my_cell{41} = '(mol)';
    my_cell{42} = '(mol)';
    my_cell{43} = '(mol)';
    my_cell{44} = '(mol)';
    my_cell{45} = '(mol)';
    my_cell{46} = '(mol)';
end
end

