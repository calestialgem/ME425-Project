function x = ask(msg, x_min, x_max)
    while true
        x = input(msg);
        if isempty(x)
            fprintf(2, "No input!\n");
        elseif ischar(x) || isstring(x)
            fprintf(2, "The input is a string!\n")
        elseif ~isscalar(x)
            fprintf(2, "The input is not a scalar!\n")
        elseif x < x_min
            fprintf(2, "The input is too small!\n");
        elseif x > x_max
            fprintf(2, "The input is too big!\n");
        else
            break;
        end
    end
end
