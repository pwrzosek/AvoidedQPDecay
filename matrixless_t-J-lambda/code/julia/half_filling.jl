function hf_run(input::OrderedDict; howmany = 1, krylov_dim = 30, which = :SR)
    global system = make_system(input)
    @time global hf_basis = hf_make_basis(system)
    global ik = 2.0 * pi * im * system.momentum / system.size

    return eigsolve(hf_model, length(hf_basis), howmany, which, ComplexF64, krylovdim = krylov_dim, ishermitian = true)
end

function hf_make_basis(system::System)::Basis
    if isodd(system.size)
        error("Requested 'system.size' is odd! Only even sizes are supported.")
    end
    if system.magnetization < 0 || 2 * system.magnetization > system.size
        error("Wrong magnetization sector in the input file! Possible subspaces are denoted by integers from 0 to N/2.")
    end

    n_spins_up::Int64 = div(system.size, 2) - system.magnetization
    subspace_size::Int64 = binomial(system.size, n_spins_up)
    state::Int64 = n_spins_up == 0 ? 0 : sum(n -> 1 << n, 0 : (n_spins_up - 1))
    mask = sum(1 << k for k in 0:2:(system.size - 1))

    basis::Basis = Basis()

    index = 0
    for _ in 1:subspace_size
        magnon_state = sublattice_rotation(state, mask)
        if has_momentum(magnon_state, system)
            rep_state = hf_get_representative(magnon_state, system)
            if get(basis, rep_state, nothing) === nothing
                push!(basis, rep_state => (index += 1))
            end
        end
        state = get_next_state(state, system)
    end

    return basis
end

function sublattice_rotation(state::Int64, mask::Int64)::Int64
    return xor(state, mask)
end

function has_momentum(state::Int64, system::System)::Bool
    ### system.size must divide system.momentum times periodicity
    return rem(system.momentum * get_periodicity(state, system), system.size) == 0
end

function get_periodicity(state::Int64, system::System)::Int64
    ### initialize some constants for faster evaluation
    l::Int = system.size
    highest_bit::Int = 1 << (l - 1)
    highest_value::Int = (1 << l) - 1

    ### initialize state translation
    state_translation = state

    ### periodicity is smallest positive number of translations
    ### that transform state onto itself
    periodicity::Int64 = 1
    while state != (state_translation = bitmov(state_translation, l, false, hb = highest_bit, hv = highest_value))
        periodicity += 1
    end

    return periodicity
end

@inline bitmov(s::Int, l::Int, f::Bool = false; hb::Int = 1 << (l - 1), hv::Int = (1 << l) - 1) = f ? 2s - div(s, hb) * hv : div(s, 2) + rem(s, 2) * hb

function hf_get_representative(state::Int64, system::System)::Int64
    l::Int = system.size
    highest_bit::Int = 1 << (l - 1)
    highest_value::Int = (1 << l) - 1

    result = state
    new_state = state
    while state != (new_state = bitmov(new_state, l, false, hb = highest_bit, hv = highest_value))
        if result > new_state
            result = new_state
        end
    end

    return result
end

function get_state_info(state::Int64, system::System)::Tuple{Bool, Int64, Int64, Int64}
    l::Int = system.size
    highest_bit::Int = 1 << (l - 1)
    highest_value::Int = (1 << l) - 1

    representative = state

    ### loop over translations of state
    new_state = state
    distance::Int64 = 0
    periodicity::Int64 = 1
    while state != (new_state = bitmov(new_state, l, false, hb = highest_bit, hv = highest_value))
        if representative > new_state
            representative = new_state
            distance = periodicity
        end
        periodicity += 1
    end
    has_k = rem(system.momentum * periodicity, system.size) == 0
    return (has_k, representative, periodicity, distance)
end

function get_next_state(state::Int64, system::System)::Int64
    count = 0
    ### loop over <i,j> site pairs (bonds)
    for i in 0:(system.size - 2)
        j = i + 1

        ### get numeric value at i and j bit positions
        i_value, j_value = 1 << i, 1 << j

        ### check if there is bit 1 at site i
        if state & i_value > 0
            ### check if there is bit 1 at site j
            if state & j_value > 0
                ### clear bit at site i
                state &= ~i_value

                ### raise counter of cleared bits
                count += 1
            else
                ### swap bits at i and j
                state = xor(state, i_value + j_value)

                ### add cleared bits to lowest positions
                state += count == 0 ? 0 : sum(n -> 1 << n, 0:(count-1))

                ### return new state
                return state
            end
        end
    end
    ### if loop ends without returning new state we return -1
    ### as there is no next state in the requestes subspace
    return -1
end

function hf_model(x::Vector{ComplexF64})::Vector{ComplexF64}
    y = zeros(ComplexF64, length(x))
    for (state, index) in hf_basis
        periodicity = get_periodicity(state, system)
        diagonal::ComplexF64 = 0;
        off_diagonal::ComplexF64 = 0;
        for i in 1:system.size
            j = mod1(i + 1, system.size)
            i_value, j_value = (1 << (i - 1)), (1 << (j - 1))
            i_bit, j_bit = div(state & i_value, i_value), div(state & j_value, j_value)
            if i_bit == j_bit
                new_state = xor(state, i_value + j_value)
                has_k, rep_state, rep_periodicity, distance = get_state_info(new_state, system)
                if has_k
                    coefficient = 0.5 * exp(ik * distance) * sqrt(periodicity / rep_periodicity)
                    off_diagonal += coefficient * x[hf_basis[rep_state]]
                end
            end
            i_bit, j_bit = div(state & i_value, i_value), div(state & j_value, j_value)
            diagonal -= 0.25 - 0.5 * (i_bit + j_bit) + system.interaction * i_bit * j_bit
        end
        y[index] = off_diagonal + diagonal * x[index]
    end
    return y
end
