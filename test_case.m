function [stat,mesh] = test_case(mesh,case_type)

if case_type == 2
    [stat,mesh] = tc2(mesh);
elseif case_type ==5
    [stat,mesh] = tc5(mesh);
elseif case_type ==6
    [stat,mesh] = tc6(mesh);
end