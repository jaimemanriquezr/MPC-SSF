function isEqual = isequal(A, B)
% ISEQUAL  Field-by-field comparison; function handles compare by expression (sym).
% Scalar struct or object fields (e.g. CohesionSubModel) are compared one level down.
    isEqual = true;
    for field = string(fieldnames(A).')
        fieldA = A.(field);
        fieldB = B.(field);

        if isstruct(fieldA) || (isobject(fieldA) && isscalar(fieldA) && isscalar(fieldB))
            if ~strcmp(class(fieldA), class(fieldB))
                isEqual = false;
                return;
            end
            for subField = string(fieldnames(fieldA).')
                if ~sameValue(fieldA.(subField), fieldB.(subField))
                    isEqual = false;
                    return;
                end
            end
        elseif ~sameValue(fieldA, fieldB)
            isEqual = false;
            return;
        end
    end
end

function same = sameValue(a, b)
    if isa(a, "function_handle") && isa(b, "function_handle")
        same = isequal(sym(a), sym(b));
    else
        same = isequaln(a, b);   % NaN defaults (e.g. CardinalTemperatures) count as equal
    end
end
