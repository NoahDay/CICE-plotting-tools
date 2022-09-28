function perimeter = find_perimeter(mat,R,C)
    % Returns the sum of perimeter of shapes formed with 1s
    perimeter = 0;
    for i = 1:R
        for j = 1:C
            if mat(i,j)
                perimeter = perimeter + (4 - num_of_neighbour(mat, i, j));
            end

        end
    end
end

