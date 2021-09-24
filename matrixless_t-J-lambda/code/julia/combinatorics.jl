function combinations(n::Int64, k::Int64)::Int64
    if n <= 0 || k < 0 || n < k
        return 0
    end
    if k > n - k
        k = n - k
    end
    result::Int64 = 1
    for i::Int64 in 1:k
        result *= n - i + 1
        result = div(result, i)
    end
    return result
end

function binary_permutation_index(value::Int64)::Int64
    bit_length::Int64 = 0
    digit_count::Int64 = 0

    target_value::Int64 = value
    while target_value != 0
        bit_length += 1;
        digit_count += target_value & 1
        target_value = target_value >> 1
    end

    if digit_count == bit_length
        return 1
    end

    target_value = 1 << (bit_length -= 1)
    result::Int64 = 1 + combinations(bit_length, digit_count)

    while ((digit_count -= 1) != 0) && (digit_count != bit_length)
        lower_bound::Int64 = 0;
        upper_bound::Int64 = bit_length;

        while upper_bound != lower_bound
            position::Int64 = div(upper_bound - lower_bound, 2) + lower_bound
            if (target_value + (1 << position)) <= value
                lower_bound = position + 1
            else
                upper_bound = position
            end
        end

        target_value += 1 << (upper_bound -= 1)
        bit_length = upper_bound

        result += combinations(bit_length, digit_count)
    end

    return result
end

function get_binary_permutation(index::Int64, n::Int64, k::Int64)::Int64
    if k > n || index > combinations(n, k) || index < 1
        throw("Requested permutation does not exist!")
    end

    result::Int64 = 0
    for exponent::Int64 in 0:(k-1)
        result += 1 << exponent
    end

    upper_bound::Int64 = n
    for digit::Int64 in 1:k
        exponent::Int64 = k - digit
        result -= 1 << exponent
        while (exponent <= upper_bound) && (binary_permutation_index(result + (1 << exponent)) <= index)
            exponent += 1
        end
        result += 1 << (exponent -= 1)
        upper_bound = exponent - 1
    end

    return result
end
