function elem_buff = select_cell(elem,cells_to_simulate,n_cells)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if iscell(elem)
    elem_buff = cell(n_cells,1);
    for i = 1:n_cells
        elem_buff{i} = elem{cells_to_simulate(i)};
    end

else
    elem_buff = zeros(n_cells,1);
    for i = 1:n_cells
        elem_buff(i) = elem(cells_to_simulate(i));
    end
end

end

