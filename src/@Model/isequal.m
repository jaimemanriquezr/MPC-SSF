function isEqual = isequal(A, B)
    isEqual = true;
    for field = string(fieldnames(A).')
        fieldA = A.(field);
        fieldB = B.(field);

        if isequal(class(fieldA), "function_handle")
            fieldA = sym(fieldA);
            fieldB = sym(fieldB);
        end

        if isa(field, "struct")
            for subField = string(fieldnames(fieldA).')
                subFieldA = fieldA.(subField);
                subFieldB = fieldB.(subField);
                if isequal(class(subFieldA), "function_handle")
                    subFieldA = sym(subFieldA);
                    subFieldB = sym(subFieldB);
                end

                if ~isequal(subFieldA,subFieldB)
                    isEqual = false;
                    return;
                end
            end
        end

        if ~isequal(fieldA, fieldB)
            isEqual = false;
            return;
        end
    end
end